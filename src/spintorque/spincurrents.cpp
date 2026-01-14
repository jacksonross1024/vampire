//------------------------------------------------------------------------------
// spincurrents.cpp
//
// Minimal 1D (along st::internal::stz) transient spin-accumulation
// solver, with:
//  - implicit (backward-Euler) diffusion step (only dzz) on a fine z-grid
//  - explicit Heun step for local precession/relaxation + drift-divergence +
//    optional demagnetisation-driven source
//  - Neumann outer boundaries (no Dirichlet spin sinks), internal Robin-like
//    interface coupling via st::internal::r_int_edge
//
// Notes:
//  - This implementation is designed for stack-by-stack (x,y) independence:
//    no lateral gradients/flows are taken.
//  - Laser-driven charge transient solver (coarse grid) for Eqs. (1)-(3),
//    producing a 1D charge current Jc(z) used in the spin drift term (Eq. 4).
//------------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <algorithm>

#include "internal.hpp"
#include "material.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

namespace st {
namespace internal {

// Physical constants (SI)
static constexpr double kHbar = 1.05457162e-34;      // J s
static constexpr double kMuB  = 9.27400968e-24;      // J/T
static constexpr double kE    = 1.60217662e-19;      // C


// Planck constant and speed of light (SI)
static constexpr double kPlanck = 6.62607015e-34;     // J s
static constexpr double kC0     = 2.99792458e8;       // m/s

static inline double gaussian_pulse(const double t, const double t0, const double fwhm){
   if(fwhm <= 0.0) return 0.0;
   // sigma from FWHM: FWHM = 2*sqrt(2*ln2)*sigma
   const double sigma = fwhm / (2.0*std::sqrt(2.0*std::log(2.0)));
   const double x = (t - t0) / sigma;
   return std::exp(-0.5*x*x);
}

static inline double laser_power_density(const double z_m, const double t_s){
   if(!sc1d_laser_enable || sc1d_laser_Q0 <= 0.0) return 0.0;
   const double d = std::max(sc1d_optical_absorption_length, 1e-18);
   const double env_t = gaussian_pulse(t_s, sc1d_laser_t0, sc1d_laser_fwhm);
   return sc1d_laser_Q0 * std::exp(-z_m / d) * env_t;
}

// Matches the constant used in spinaccumulation.cpp
static constexpr double kAtomcellVolume = 7.666015625e-30; // m^3
static constexpr double kInvMuB = 1.0 / kMuB;
static constexpr double kInvE   = 1.0 / kE;

// Stable-axis freezing (compile-time flag, no user interface by request)
static constexpr bool   kFreezeAxisWhenMsmall = true;
static constexpr double kFreezeFracOfM0       = 0.05;  // update axis only when |m| > frac * |m0|
// Interpolate material/axis parameters between neighbouring coarse microcells.
// Keep OFF by default to avoid unintended smoothing across sharp interfaces.
static constexpr bool   kInterpolateBetweenMicrocells = false;
static constexpr double kEps                  = 1e-30;

//------------------------------------------------------------------------------
// Local helpers
//------------------------------------------------------------------------------

static inline double norm3(const double x, const double y, const double z) {
   return std::sqrt(x*x + y*y + z*z);
}

static inline void normalize3(double &x, double &y, double &z) {
   const double n = norm3(x,y,z);
   if(n > kEps) { x/=n; y/=n; z/=n; }
}

static inline void cross3(const double ax, const double ay, const double az,
                          const double bx, const double by, const double bz,
                          double &cx, double &cy, double &cz) {
   cx = ay*bz - az*by;
   cy = az*bx - ax*bz;
   cz = ax*by - ay*bx;
}

// Thomas tridiagonal solver for a single RHS.
// a: lower diag (a[0] unused), b: diag, c: upper diag (c[n-1] unused)
static inline void thomas_solve(const std::vector<double>& a,
                                const std::vector<double>& b,
                                const std::vector<double>& c,
                                std::vector<double>& d) {
   const int n = static_cast<int>(d.size());
   std::vector<double> cp(n,0.0);
   std::vector<double> dp(n,0.0);

   double denom = b[0];
   if(std::fabs(denom) < kEps) denom = (denom >= 0 ? kEps : -kEps);
   cp[0] = c[0] / denom;
   dp[0] = d[0] / denom;

   for(int i=1;i<n;i++){
      denom = b[i] - a[i]*cp[i-1];
      if(std::fabs(denom) < kEps) denom = (denom >= 0 ? kEps : -kEps);
      cp[i] = (i < n-1) ? (c[i] / denom) : 0.0;
      dp[i] = (d[i] - a[i]*dp[i-1]) / denom;
   }

   d[n-1] = dp[n-1];
   for(int i=n-2;i>=0;i--){
      d[i] = dp[i] - cp[i]*d[i+1];
   }
}

//------------------------------------------------------------------------------
// Charge transient solver (coarse grid): Eqs. (1)-(3) in 1D
//
// We evolve non-equilibrium charge density ns and excess charge density ne on the
// coarse microcell grid (z-only), then build Jc(z-edges) for use in the spin
// accumulation drift term.
//
// Simplifications (per latest implementation plan):
//  - Seebeck term and explicit Te-gradient terms are not included yet.
//  - Temperature-gradient driven superdiffusive term is folded into ve = De/d.
//  - All quantities are treated as charge densities (C/m^3) and charge currents (A/m^2).
//  - Boundaries enforce zero normal charge current (Jc=0) and zero normal ns flux (Jns=0).
//------------------------------------------------------------------------------
static inline void update_charge_transients_1d_stack(const int start_cell,
                                                    const double dt,
                                                    const double t_s,
                                                    double* ns,          // size ncz
                                                    double* ne,          // size ncz
                                                    double* Jc_edge,     // size ncz+1
                                                    const int ncz){
   const double dz = micro_cell_thickness * 1.0e-10; // m

   // Edge-averaged diffusion constants and ve = De/d (all in SI)
   std::vector<double> De_cell(ncz,0.0), ve_cell(ncz,0.0);
   for(int k=0;k<ncz;k++){
      const int cell = start_cell + k;
      const double De = std::max(0.0, diffusion[cell]); // m^2/s
      De_cell[k] = De;
      const double dopt = std::max(sc1d_optical_absorption_length, 1e-18);
      ve_cell[k] = De / dopt; // m/s
   }

   std::vector<double> De_edge(ncz+1,0.0), ve_edge(ncz+1,0.0);
   // boundary edges set to 0 to enforce no-flux for ns and no normal Jc at surfaces
   for(int e=1;e<ncz;e++){
      const double DL = De_cell[e-1];
      const double DR = De_cell[e];
      // harmonic mean for diffusion at material steps
      if(DL > 0.0 && DR > 0.0) De_edge[e] = 2.0*DL*DR/(DL+DR);
      else De_edge[e] = 0.0;
      ve_edge[e] = 0.5*(ve_cell[e-1] + ve_cell[e]);
      if(ve_edge[e] < 0.0) ve_edge[e] = 0.0;
   }

   // Solve ns implicitly (backward Euler) with advection-diffusion-reaction + laser source
   std::vector<double> a(ncz,0.0), b(ncz,0.0), c(ncz,0.0), rhs(ncz,0.0);

   const double tau_s = std::max(sc1d_tau_s, 0.0);
   const double inv_tau = (tau_s > 0.0) ? (1.0/tau_s) : 0.0;

   for(int k=0;k<ncz;k++){
      const double DeL = De_edge[k];
      const double DeR = De_edge[k+1];
      const double veL = ve_edge[k];
      const double veR = ve_edge[k+1];

      a[k] = -(dt/dz)*veL - dt*DeL/(dz*dz);
      b[k] = 1.0 + dt*inv_tau + (dt/dz)*veR + dt*(DeL + DeR)/(dz*dz);
      c[k] = -dt*DeR/(dz*dz);

      // Laser source in Eq. (2): + e*eta*Q*lambda/(h*c)
      const double zc = (k + 0.5) * dz; // m
      const double Q = laser_power_density(zc, t_s); // W/m^3
      const double src = (kE * sc1d_laser_eta * Q * sc1d_laser_wavelength) / (kPlanck * kC0); // C/(m^3 s)

      rhs[k] = ns[k] + dt*src;
   }
   // Boundary fixes (De_edge[0]=De_edge[ncz]=0, ve_edge[0]=ve_edge[ncz]=0) already yield correct.
   thomas_solve(a,b,c,rhs);
   for(int k=0;k<ncz;k++) ns[k] = rhs[k];

   // Build non-equilibrium current edges Jns = ve*ns - De*d(ns)/dz (upwind for ve>=0), with Jns=0 at boundaries
   std::vector<double> Jns_edge(ncz+1,0.0);
   for(int e=1;e<ncz;e++){
      const int kL = e-1;
      const int kR = e;
      const double ve = ve_edge[e];
      const double De = De_edge[e];
      const double adv = ve * ns[kL]; // upwind (ve assumed >=0)
      const double diff = -De * (ns[kR] - ns[kL]) / dz;
      Jns_edge[e] = adv + diff; // A/m^2
   }

   // Solve ne implicitly from continuity Eq. (3): d ne/dt = -d Jc/dz, where Jc = -De*d(ne)/dz + Jns
   for(int k=0;k<ncz;k++){
      const double DeL = De_edge[k];
      const double DeR = De_edge[k+1];
      a[k] = -dt*DeL/(dz*dz);
      b[k] = 1.0 + dt*(DeL + DeR)/(dz*dz);
      c[k] = -dt*DeR/(dz*dz);

      const double divJns = (Jns_edge[k+1] - Jns_edge[k]) / dz;
      rhs[k] = ne[k] - dt*divJns;
   }
   thomas_solve(a,b,c,rhs);
   for(int k=0;k<ncz;k++) ne[k] = rhs[k];

   // Build total charge current edges Jc, enforcing Jc=0 at boundaries
   Jc_edge[0]   = 0.0;
   Jc_edge[ncz] = 0.0;
   for(int e=1;e<ncz;e++){
      const int kL = e-1;
      const int kR = e;
      const double De = De_edge[e];
      const double diff = -De * (ne[kR] - ne[kL]) / dz;
      Jc_edge[e] = diff + Jns_edge[e];
   }
}


//------------------------------------------------------------------------------
// Initialisation
//------------------------------------------------------------------------------

void initialise_spincurrents_1d(){
   // Determine fine subdivisions per coarse microcell thickness (Angstrom -> Angstrom)
   const double dzA = std::max(0.01, sc1d_fine_dz);
   const double DzA = std::max(0.01, micro_cell_thickness);
   int nsub = static_cast<int>(std::ceil(DzA / dzA));
   nsub = std::max(1, std::min(128, nsub));

   sc1d_nsub = nsub;
   sc1d_nf   = nsub * num_microcells_per_stack;

   // Local stacks owned by this rank
   sc1d_local_stacks.clear();
   #ifdef MPICF
      sc1d_local_stacks = mpi_stack_list_y;
   #else
      sc1d_local_stacks.resize(num_stacks_y);
      for(int s=0; s<num_stacks_y; ++s) sc1d_local_stacks[s]=s;
   #endif

   // Map global stack -> local index
   sc1d_stack_local_index.assign(num_stacks_y, -1);
   for(int i=0; i<(int)sc1d_local_stacks.size(); ++i){
      sc1d_stack_local_index[sc1d_local_stacks[i]] = i;
   }

   // Fine-grid state per local stack
   sc1d_Sfine.assign( (size_t)sc1d_local_stacks.size() * (size_t)sc1d_nf * 3u, 0.0 );

   // Coarse-grid charge transient state per local stack
   const int ncz = num_microcells_per_stack;
   sc1d_ns_coarse.assign( (size_t)sc1d_local_stacks.size() * (size_t)ncz, 0.0 );
   sc1d_ne_coarse.assign( (size_t)sc1d_local_stacks.size() * (size_t)ncz, 0.0 );
   sc1d_Jc_edge_coarse.assign( (size_t)sc1d_local_stacks.size() * (size_t)(ncz+1), 0.0 );

   // Demag tracking arrays per coarse microcell
   const int total_microcells = num_x_stacks * num_y_stacks * num_microcells_per_stack;
   sc1d_m_prev_mag.assign(total_microcells, 0.0);
   sc1d_m0_mag.assign(total_microcells, 0.0);
   sc1d_mhat_ref.assign(total_microcells*3u, 0.0);
   for(int c=0;c<total_microcells;c++){
      sc1d_mhat_ref[3*c + 2] = 1.0; // default axis +z
   }
}

//------------------------------------------------------------------------------
// Main update
//------------------------------------------------------------------------------

void calculate_spin_accumulation_1d(){
   // If solver was not initialised (e.g., enable flag set after initialise), do it now.
   if(sc1d_nf <= 0 || sc1d_Sfine.empty()){
      initialise_spincurrents_1d();
   }

   // Clear outputs
   std::fill(spin_torque.begin(), spin_torque.end(), 0.0);
   std::fill(sa_final.begin(),    sa_final.end(),    0.0);
   std::fill(ns_final.begin(),    ns_final.end(),    0.0);
   std::fill(ne_final.begin(),    ne_final.end(),    0.0);
   std::fill(jc_final.begin(),    jc_final.end(),    0.0);
   std::fill(js_final.begin(),    js_final.end(),    0.0);

   const double dt = mp::dt;
   if(dt <= 0.0) return;
   const double dt_half = 0.5*dt;
   // Increment internal step counter (used for laser timing and charge/ spin stride logic)
   ++sc1d_step_counter;
   const double t_s = (double)sc1d_step_counter * dt;

   const int ncz  = num_microcells_per_stack;
   const int nsub = sc1d_nsub;
   const int nf   = sc1d_nf;

   // Fine spacing in metres
   const double dz = (micro_cell_thickness * 1.0e-10) / (double)nsub;

   // Temporary per-stack buffers
   std::vector<double> D(nf,0.0), Bc(nf,0.0), Bd(nf,0.0), lsf(nf,0.0), lphi(nf,0.0), Jsd(nf,0.0), chi(nf,0.0), dmdt(nf,0.0);
   std::vector<double> mhat(3*nf,0.0), msrc(3*nf,0.0);
   std::vector<double> alpha_edge(std::max(0,nf-1),0.0);
   std::vector<double> a_tr(nf,0.0), b_tr(nf,0.0), c_tr(nf,0.0);
   std::vector<double> rhs(nf,0.0);

   // Drift flux at edges (vector, size nf+1)
   std::vector<double> Jd_edge(3*(nf+1),0.0);

   // Charge current at fine edges (A/m^2), built from coarse solver; size nf+1
   std::vector<double> Jc_edge_fine(nf+1,0.0);

   auto Sref = [&](double* Sbase, int i, int comp) -> double& { return Sbase[3*i + comp]; };

   // Loop over local stacks
   for(int ls=0; ls<(int)sc1d_local_stacks.size(); ++ls){
      const int stack = sc1d_local_stacks[ls];
      const int start_cell = stack_index_y[stack];
      double* Sbase = sc1d_Sfine.data() + (size_t)ls * (size_t)nf * 3u;

// --- Charge transient update on coarse grid (Eqs. 1–3) ---
double* ns_coarse = sc1d_ns_coarse.data() + (size_t)ls * (size_t)ncz;
double* ne_coarse = sc1d_ne_coarse.data() + (size_t)ls * (size_t)ncz;
double* Jc_coarse = sc1d_Jc_edge_coarse.data() + (size_t)ls * (size_t)(ncz + 1);

const bool do_charge_update = (sc1d_charge_stride <= 1) || ((sc1d_step_counter % (unsigned long)sc1d_charge_stride) == 0UL);
if(do_charge_update){
   update_charge_transients_1d_stack(start_cell, dt, t_s, ns_coarse, ne_coarse, Jc_coarse, ncz);
}

// Build fine-edge Jc by linear interpolation within each coarse microcell (no lateral flow)
// Boundaries use Jc=0 as enforced by the coarse solver.
for(int e=0; e<=nf; ++e) Jc_edge_fine[e] = 0.0;
for(int k=0;k<ncz;k++){
   const double Jl = Jc_coarse[k];
   const double Jr = Jc_coarse[k+1];
   for(int j=0;j<=nsub;j++){
      const int e = k*nsub + j;
      if(e > nf) continue;
      const double f = (double)j / (double)nsub;
      Jc_edge_fine[e] = (1.0 - f)*Jl + f*Jr;
   }
}

      // Coarse dm/dt per microcell in this stack
      std::vector<double> dm_dt_coarse(ncz,0.0);

      for(int k=0;k<ncz;k++){
         const int cell = start_cell + k;
         const double mx = m[3*cell+0];
         const double my = m[3*cell+1];
         const double mz = m[3*cell+2];
         const double mm = norm3(mx,my,mz);

         if(sc1d_m0_mag[cell] <= 0.0 && mm > 0.0) sc1d_m0_mag[cell] = mm;

         const double thresh = kFreezeFracOfM0 * std::max(sc1d_m0_mag[cell], 1e-12);

         // Update stable axis if permitted
         if(!kFreezeAxisWhenMsmall || mm > thresh){
            if(mm > kEps){
               sc1d_mhat_ref[3*cell+0] = mx/mm;
               sc1d_mhat_ref[3*cell+1] = my/mm;
               sc1d_mhat_ref[3*cell+2] = mz/mm;
            }
         }

         // Avoid an artificial "initial transient" when the solver starts:
         // on the first call, initialise the history and set dm/dt = 0.
         if(sc1d_step_counter <= 1UL){
            dm_dt_coarse[k] = 0.0;
            sc1d_m_prev_mag[cell] = mm;
         } else {
            dm_dt_coarse[k] = (mm - sc1d_m_prev_mag[cell]) / dt;
            sc1d_m_prev_mag[cell] = mm;
         }
      }

      // Build fine-grid interpolated parameters and axes (no interpolation across interfaces)
      for(int k=0;k<ncz;k++){
         const int cell1 = start_cell + k;
         const int cell2 = (k < ncz-1) ? (cell1 + 1) : cell1;
         const bool is_interface = (k < ncz-1) ? (r_int_edge[cell1] > 0.0) : true;

         const double mx1 = m[3*cell1+0], my1=m[3*cell1+1], mz1=m[3*cell1+2];
         const double mx2 = m[3*cell2+0], my2=m[3*cell2+1], mz2=m[3*cell2+2];

         const double bc1 = beta_cond[cell1], bd1 = beta_diff[cell1], d1 = diffusion[cell1];
         const double bc2 = beta_cond[cell2], bd2 = beta_diff[cell2], d2 = diffusion[cell2];

         const double lsf1 = lambda_sdl[cell1], lsf2 = lambda_sdl[cell2];
         const double lph1 = lambda_phi[cell1], lph2 = lambda_phi[cell2];

         const double jsd1 = sd_exchange[cell1], jsd2 = sd_exchange[cell2];
         const double chi1 = chi_demag[cell1], chi2 = chi_demag[cell2];

         const double dm1 = dm_dt_coarse[k];
         const double dm2 = (k < ncz-1) ? dm_dt_coarse[k+1] : dm1;

         const double sx = sc1d_mhat_ref[3*cell1+0];
         const double sy = sc1d_mhat_ref[3*cell1+1];
         const double sz = sc1d_mhat_ref[3*cell1+2];

         for(int j=0;j<nsub;j++){
            const int i = k*nsub + j;
            double w = 0.0;
            if(kInterpolateBetweenMicrocells && k < ncz-1 && !is_interface){
               w = (j + 0.5) / (double)nsub;
            }

            const double mx = (1.0-w)*mx1 + w*mx2;
            const double my = (1.0-w)*my1 + w*my2;
            const double mz = (1.0-w)*mz1 + w*mz2;

            double hx = mx, hy = my, hz = mz;
            const double hm = norm3(hx,hy,hz);
            if(hm > kEps){ hx/=hm; hy/=hm; hz/=hm; }
            else { hx = sx; hy = sy; hz = sz; }

            mhat[3*i+0]=hx; mhat[3*i+1]=hy; mhat[3*i+2]=hz;
            msrc[3*i+0]=sx; msrc[3*i+1]=sy; msrc[3*i+2]=sz;

            Bc[i]   = (1.0-w)*bc1  + w*bc2;
            Bd[i]   = (1.0-w)*bd1  + w*bd2;
            D[i]    = (1.0-w)*d1   + w*d2;
            lsf[i]  = (1.0-w)*lsf1 + w*lsf2;
            lphi[i] = (1.0-w)*lph1 + w*lph2;
            Jsd[i]  = (1.0-w)*jsd1 + w*jsd2;
            chi[i]  = (1.0-w)*chi1 + w*chi2;
            dmdt[i] = (1.0-w)*dm1  + w*dm2;
         }
      }

      // Build alpha_edge for diffusion operator, including Robin-like interface resistance
      if(nf >= 2){
         for(int i=0;i<nf-1;i++){
            double Rint = 0.0;
            if( (i % nsub) == (nsub-1) ){
               const int k = i / nsub;
               if(k < ncz-1){
                  const int cell = start_cell + k;
                  Rint = r_int_edge[cell];
               }
            }
            const double Di = D[i];
            const double Dj = D[i+1];
            if(Di < kEps || Dj < kEps){
               alpha_edge[i] = 0.0;
               continue;
            }
            const double Rseries = dz/(2.0*Di) + Rint + dz/(2.0*Dj);
            if(Rseries <= 0.0){
               alpha_edge[i] = 0.0;
               continue;
            }
            const double C = 1.0 / Rseries;     // m/s
            alpha_edge[i] = C / dz;             // 1/s
         }
      }

      // Assemble tridiagonal matrix for backward-Euler half-step diffusion
      if(nf == 1){
         a_tr[0]=0.0; b_tr[0]=1.0; c_tr[0]=0.0;
      } else {
         a_tr[0] = 0.0;
         b_tr[0] = 1.0 + dt_half*alpha_edge[0];
         c_tr[0] = -dt_half*alpha_edge[0];
         for(int i=1;i<nf-1;i++){
            a_tr[i] = -dt_half*alpha_edge[i-1];
            b_tr[i] = 1.0 + dt_half*(alpha_edge[i-1] + alpha_edge[i]);
            c_tr[i] = -dt_half*alpha_edge[i];
         }
         a_tr[nf-1] = -dt_half*alpha_edge[nf-2];
         b_tr[nf-1] = 1.0 + dt_half*alpha_edge[nf-2];
         c_tr[nf-1] = 0.0;
      }

      // Helper: diffusion half-step with boundary flux chosen so that total spin current is Neumann (J_total=0)
      auto diffusion_half_step = [&](double* S){
         // Boundary drift flux
         const double jdl0x = -Bc[0]    * Jc_edge_fine[0]  * mhat[0];
         const double jdl0y = -Bc[0]    * Jc_edge_fine[0]  * mhat[1];
         const double jdl0z = -Bc[0]    * Jc_edge_fine[0]  * mhat[2];
         const double jdrNx = -Bc[nf-1] * Jc_edge_fine[nf] * mhat[3*(nf-1)+0];
         const double jdrNy = -Bc[nf-1] * Jc_edge_fine[nf] * mhat[3*(nf-1)+1];
         const double jdrNz = -Bc[nf-1] * Jc_edge_fine[nf] * mhat[3*(nf-1)+2];

         // Diffusive boundary flux cancels drift: J_diff = -J_drift
         const double Jleft[3]  = { -jdl0x, -jdl0y, -jdl0z };
         const double Jright[3] = { -jdrNx, -jdrNy, -jdrNz };

         for(int comp=0; comp<3; ++comp){
            for(int i=0;i<nf;i++) rhs[i] = Sref(S,i,comp);
            rhs[0]     += dt_half * (Jleft[comp] / dz);
            rhs[nf-1]  -= dt_half * (Jright[comp] / dz);
            thomas_solve(a_tr,b_tr,c_tr,rhs);
            for(int i=0;i<nf;i++) Sref(S,i,comp) = rhs[i];
         }
      };

      // First diffusion half-step
      diffusion_half_step(Sbase);

      // Build drift fluxes at edges for explicit step
      // Jd_edge[0] and Jd_edge[nf] use node values; interior edges use averaged (normalized) mhat and averaged beta
      Jd_edge[0] = -Bc[0]*Jc_edge_fine[0]*mhat[0];
      Jd_edge[1] = -Bc[0]*Jc_edge_fine[0]*mhat[1];
      Jd_edge[2] = -Bc[0]*Jc_edge_fine[0]*mhat[2];
      for(int i=0;i<nf-1;i++){
         const double bc = 0.5*(Bc[i] + Bc[i+1]);
         double hx = mhat[3*i+0] + mhat[3*(i+1)+0];
         double hy = mhat[3*i+1] + mhat[3*(i+1)+1];
         double hz = mhat[3*i+2] + mhat[3*(i+1)+2];
         normalize3(hx,hy,hz);
         Jd_edge[3*(i+1)+0] = -bc*Jc_edge_fine[i+1]*hx;
         Jd_edge[3*(i+1)+1] = -bc*Jc_edge_fine[i+1]*hy;
         Jd_edge[3*(i+1)+2] = -bc*Jc_edge_fine[i+1]*hz;
      }
      Jd_edge[3*nf+0] = -Bc[nf-1]*Jc_edge_fine[nf]*mhat[3*(nf-1)+0];
      Jd_edge[3*nf+1] = -Bc[nf-1]*Jc_edge_fine[nf]*mhat[3*(nf-1)+1];
      Jd_edge[3*nf+2] = -Bc[nf-1]*Jc_edge_fine[nf]*mhat[3*(nf-1)+2];

      // Explicit Heun step for local terms + drift divergence + demag source
      std::vector<double> k1(3*nf,0.0), k2(3*nf,0.0), Spred(3*nf,0.0);

      auto compute_rhs = [&](const double* S, std::vector<double>& kout){
         for(int i=0;i<nf;i++){
            // drift divergence term (no lateral flow)
            const double divx = -(Jd_edge[3*(i+1)+0] - Jd_edge[3*i+0]) / dz;
            const double divy = -(Jd_edge[3*(i+1)+1] - Jd_edge[3*i+1]) / dz;
            const double divz = -(Jd_edge[3*(i+1)+2] - Jd_edge[3*i+2]) / dz;

            // demag-driven source: (-chi * d|m|/dt) * m_ref
            const double src_amp = -chi[i] * dmdt[i];
            const double srcx = src_amp * msrc[3*i+0];
            const double srcy = src_amp * msrc[3*i+1];
            const double srcz = src_amp * msrc[3*i+2];

            // relaxation times
            const double denom = std::max(1e-12, 1.0 - Bc[i]*Bd[i]);
            const double BBp = 1.0 / std::sqrt(denom);
            const double lambda_sf = lsf[i] * BBp;
            const double Di = std::max(D[i], 1e-30);
            const double tau_sf = (lambda_sf > 0.0) ? (lambda_sf*lambda_sf / Di) : 1e30;
            const double tau_phi = (lphi[i]  > 0.0) ? (lphi[i]*lphi[i]   / Di) : 1e30;

            // exchange precession frequency
            const double omega = (Jsd[i] > 0.0) ? (Jsd[i] / (2.0*kHbar)) : 0.0;

            const double sx = S[3*i+0];
            const double sy = S[3*i+1];
            const double sz = S[3*i+2];

            const double hx = mhat[3*i+0];
            const double hy = mhat[3*i+1];
            const double hz = mhat[3*i+2];

            // S_perp = S - (S·h)h
            const double sdot = sx*hx + sy*hy + sz*hz;
            const double spx = sx - sdot*hx;
            const double spy = sy - sdot*hy;
            const double spz = sz - sdot*hz;

            // precession term: -omega (S x h)
            double cx,cy,cz;
            cross3(sx,sy,sz, hx,hy,hz, cx,cy,cz);
            const double prex = -omega * cx;
            const double prey = -omega * cy;
            const double prez = -omega * cz;

            const double relx = -sx / tau_sf;
            const double rely = -sy / tau_sf;
            const double relz = -sz / tau_sf;

            const double dephx = -spx / tau_phi;
            const double dephy = -spy / tau_phi;
            const double dephz = -spz / tau_phi;

            kout[3*i+0] = divx + srcx + prex + relx + dephx;
            kout[3*i+1] = divy + srcy + prey + rely + dephy;
            kout[3*i+2] = divz + srcz + prez + relz + dephz;
         }
      };

      // k1
      compute_rhs(Sbase, k1);
      for(int i=0;i<3*nf;i++) Spred[i] = Sbase[i] + dt * k1[i];

      // k2
      compute_rhs(Spred.data(), k2);

      for(int i=0;i<3*nf;i++) Sbase[i] = Sbase[i] + 0.5*dt*(k1[i] + k2[i]);

      // Second diffusion half-step
      diffusion_half_step(Sbase);

      // Coarse average and torque mapping
      for(int k=0;k<ncz;k++){
         const int cell = start_cell + k;
         
         // Store coarse charge quantities
         ns_final[cell] = ns_coarse[k];
         ne_final[cell] = ne_coarse[k];
         jc_final[cell] = 0.5 * (Jc_coarse[k] + Jc_coarse[k+1]);

         double ax=0.0, ay=0.0, az=0.0;
         double jsx_sum=0.0, jsy_sum=0.0, jsz_sum=0.0;

         for(int j=0;j<nsub;j++){
            const int i = k*nsub + j;
            ax += Sbase[3*i+0];
            ay += Sbase[3*i+1];
            az += Sbase[3*i+2];

            // Calculate local spin current Js = J_drift + J_diff
            // 1. Drift: - beta * Jc * mhat
            const double Jc_loc = 0.5*(Jc_edge_fine[i] + Jc_edge_fine[i+1]);
            const double jdr_x = -Bc[i] * Jc_loc * mhat[3*i+0];
            const double jdr_y = -Bc[i] * Jc_loc * mhat[3*i+1];
            const double jdr_z = -Bc[i] * Jc_loc * mhat[3*i+2];

            // 2. Diffusion: - D * dS/dz
            double dSdz_x, dSdz_y, dSdz_z;
            if(i==0){ // Left boundary
               dSdz_x = (Sbase[3*(i+1)+0] - Sbase[3*i+0])/dz;
               dSdz_y = (Sbase[3*(i+1)+1] - Sbase[3*i+1])/dz;
               dSdz_z = (Sbase[3*(i+1)+2] - Sbase[3*i+2])/dz;
            } else if(i==nf-1){ // Right boundary
               dSdz_x = (Sbase[3*i+0] - Sbase[3*(i-1)+0])/dz;
               dSdz_y = (Sbase[3*i+1] - Sbase[3*(i-1)+1])/dz;
               dSdz_z = (Sbase[3*i+2] - Sbase[3*(i-1)+2])/dz;
            } else { // Central difference
               dSdz_x = (Sbase[3*(i+1)+0] - Sbase[3*(i-1)+0])/(2.0*dz);
               dSdz_y = (Sbase[3*(i+1)+1] - Sbase[3*(i-1)+1])/(2.0*dz);
               dSdz_z = (Sbase[3*(i+1)+2] - Sbase[3*(i-1)+2])/(2.0*dz);
            }
            
            const double jdiff_x = -D[i] * dSdz_x;
            const double jdiff_y = -D[i] * dSdz_y;
            const double jdiff_z = -D[i] * dSdz_z;

            jsx_sum += (jdr_x + jdiff_x);
            jsy_sum += (jdr_y + jdiff_y);
            jsz_sum += (jdr_z + jdiff_z);
         }
         const double inv = 1.0 / (double)nsub;
         ax *= inv; ay *= inv; az *= inv;

         sa_final[3*cell+0] = ax;
         sa_final[3*cell+1] = ay;
         sa_final[3*cell+2] = az;

         js_final[3*cell+0] = jsx_sum * inv;
         js_final[3*cell+1] = jsy_sum * inv;
         js_final[3*cell+2] = jsz_sum * inv;

         // Convert accumulation to effective field contribution (consistent with existing code path)
         if(cell_natom[cell] > 0.5){
            spin_torque[3*cell+0] = kAtomcellVolume * sd_exchange[cell] * ax * kInvE * kInvMuB;
            spin_torque[3*cell+1] = kAtomcellVolume * sd_exchange[cell] * ay * kInvE * kInvMuB;
            spin_torque[3*cell+2] = kAtomcellVolume * sd_exchange[cell] * az * kInvE * kInvMuB;
         }
      }
   }

#ifdef MPICF
   // Gather per-cell fields across ranks (each stack computed on exactly one rank)
   const int total_cells = num_x_stacks * num_y_stacks * num_microcells_per_stack;
   MPI_Allreduce(MPI_IN_PLACE, spin_torque.data(), 3*total_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, sa_final.data(),    3*total_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, ns_final.data(),    total_cells,   MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, ne_final.data(),    total_cells,   MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, jc_final.data(),    total_cells,   MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, js_final.data(),    3*total_cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}

} // namespace internal
} // namespace st
