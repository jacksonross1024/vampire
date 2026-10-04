//------------------------------------------------------------------------------
// spincurrents.cpp
//
// Minimal 1D (along st::internal::stz) transient spin-accumulation
// solver, with:
//  - implicit (backward-Euler) diffusion step (only dzz) on a fine z-grid
//  - explicit Heun step for local precession/relaxation + drift-divergence.
//    Equilibrium amplitude is sa_inf (direction follows m-hat). Demagnetization
//    dumps -chi*(d|m|/dt)*m-hat into S; spin-flip returns the excess to sa_inf.
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
#include "ltmp.hpp"
#include "errors.hpp"
#include "atoms.hpp"
#if ST_SC1D_TIMINGS
#include "stopwatch.hpp"
#else
struct stopwatch_t {
   void start() {}
   double elapsed_seconds() { return 0.0; }
};
#endif
#include "vio.hpp"

namespace st {
namespace internal {

// True while calculate_spin_accumulation_1d's total stopwatch is running.
// Init-time charge work is then already inside that total.
static bool sc1d_step_timer_active = false;

static void account_charge_time(const double elapsed, const double interp_before){
#if ST_SC1D_TIMINGS
   const double interp_dt = sc1d_time_interp - interp_before;
   sc1d_time_charge += elapsed - interp_dt;
   if(!sc1d_step_timer_active) sc1d_time_total += elapsed;
#else
   (void)elapsed;
   (void)interp_before;
#endif
}

static void close_sc1d_step_timer(stopwatch_t& sw_total){
#if ST_SC1D_TIMINGS
   sc1d_time_total += sw_total.elapsed_seconds();
   sc1d_time_other = sc1d_time_total
      - (sc1d_time_io + sc1d_time_charge + sc1d_time_spin_current
         + sc1d_time_spin_acc + sc1d_time_interp + sc1d_time_bcast);
#else
   (void)sw_total;
#endif
   sc1d_step_timer_active = false;
}

// Physical constants (SI)
static constexpr double kHbar = 1.05457162e-34;      // J s
static constexpr double kMuB  = 9.27400968e-24;      // J/T
static constexpr double kE    = 1.60217662e-19;      // C


// Planck constant and speed of light (SI)
static constexpr double kPlanck = 6.62607015e-34;     // J s
static constexpr double kC0     = 2.99792458e8;       // m/s
static constexpr double kEps0   = 8.854187817e-12;    // F/m
static constexpr double kBoltzmann = 1.380649e-23;    // J/K
static constexpr double kMetalN = 1.0e28;             // m^-3, Einstein n for σ fallback

static inline double gaussian_pulse(const double t, const double t0, const double fwhm, const double equilibration_offset_s){
   if(fwhm <= 0.0) return 0.0;
   
   // sigma from FWHM: FWHM = 2*sqrt(2*ln2)*sigma
   const double sigma = fwhm / (2.0*std::sqrt(2.0*std::log(2.0)));
   
   // Offset the pulse by the equilibration period (in seconds)
   const double x = (t - t0 - equilibration_offset_s) / sigma;
   
   return std::exp(-0.5*x*x);
}

//------------------------------------------------------------------------------
// Laser power density function
// Can use either built-in laser or ltmp attenuated laser depending on flag
//------------------------------------------------------------------------------
// Drive phase is the ASD integration after sim:equilibration-time-steps.
// Spin-current equilibration and ASD equilibration both leave the laser off.
static bool sc1d_drive_active(){
   return sim::time > sim::equilibration_time;
}

// Seconds from the start of the drive phase, matching the time ltmp is given.
static double drive_time_s(){
   if(!sc1d_drive_active()) return 0.0;
   return mp::dt_SI * static_cast<double>(sim::time - sim::equilibration_time);
}

static inline double laser_power_density(const double z_m, const double t_s, const unsigned long step_counter, const double dt_si){
   (void)t_s;
   (void)step_counter;
   if(!sc1d_laser_enable || sc1d_laser_Q0 <= 0.0) return 0.0;
   if(!sc1d_drive_active()) return 0.0;
   
   // Built-in laser function (always used now that thermal gradients are in spin-currents module)
   // Laser enters from top of stack (z_max), attenuates as it propagates downward
   // z=0 is at bottom, z=z_max is at top
   // Form: Q = Q0 * exp(-distance_from_top / d) * gaussian(t - t0)
   // t is measured from the start of the drive phase.
   const double d = std::max(sc1d_optical_absorption_length, 1e-18);
   const double z_max = (double)st::internal::num_microcells_per_stack * st::internal::micro_cell_thickness * 1.0e-10; // total stack height in meters
   const double distance_from_top = z_max - z_m; // distance from top surface
   const double t_drive = (double)(sim::time - sim::equilibration_time) * dt_si;
   const double env_t = gaussian_pulse(t_drive, sc1d_laser_t0, sc1d_laser_fwhm, 0.0);
   const double attenuation = std::exp(-distance_from_top / d);
   
   return sc1d_laser_Q0 * attenuation * env_t;
}

// Superdiffusion reads the ltmp pump. That value is already multiplied by the
// cell attenuation (the absorption profile, or Beer-Lambert if no profile).
// The spin-current laser stays the source only when superdiffusion is off.
static double ns_source_power_density(const double z_m){
   if(!sc1d_drive_active()) return 0.0;
   if(sc1d_superdiffusive_enable){
      return ltmp::get_attenuated_laser_power(z_m * 1.0e10, drive_time_s());
   }
   return laser_power_density(z_m, 0.0, sc1d_step_counter, mp::dt_SI);
}

//------------------------------------------------------------------------------
// Electron and phonon temperatures for the spin solver.
//
// ltmp owns the cell temperatures and the Langevin field (fields.cpp).
// A laser pulse without ltmp owns a single Te and Tp on the global TTM
// (sim::TTTe, sim::TTTp). This solver only reads those values.
// program 6 is laser-pulse. program 19 is laser-electrical-pulse.
//------------------------------------------------------------------------------

static bool laser_pulse_program(){
   return program::program == 6 || program::program == 19;
}

static bool spatial_ltmp_temperatures(){
   return ltmp::is_enabled();
}

static bool global_ttm_temperatures(){
   return !ltmp::is_enabled() && laser_pulse_program();
}

static double finite_temperature(const double T){
   if(!(T > 0.0) || !std::isfinite(T)) return sc1d_reference_temperature;
   return T;
}

// z_m is the fine-cell centre in metres, origin at the bottom of the stack.
static double get_local_electron_temperature(const double z_m){
   if(spatial_ltmp_temperatures()){
      return finite_temperature(ltmp::get_electron_temperature_at_z(z_m * 1.0e10));
   }
   if(global_ttm_temperatures()) return finite_temperature(sim::TTTe);
   return sc1d_reference_temperature;
}

static double get_local_phonon_temperature(const double z_m){
   if(spatial_ltmp_temperatures()){
      return finite_temperature(ltmp::get_phonon_temperature_at_z(z_m * 1.0e10));
   }
   if(global_ttm_temperatures()) return finite_temperature(sim::TTTp);
   return sc1d_reference_temperature;
}

// A gradient exists only on the ltmp grid. The global TTM is one Te and one Tp.
static bool sc1d_has_electron_temperature(){
   return spatial_ltmp_temperatures();
}

// spin-currents-1d-thermal-effects. D, lambda, tau_s, S and sigma would be
// replaced here by functions of Te and Tp. Those laws are not specified yet,
// so the constants are returned unchanged.
static inline double apply_thermal_effects(const double value, const double Te, const double Tp){
   if(!sc1d_thermal_effects) return value;
   (void)Te;
   (void)Tp;
   return value;
}

static inline double scale_diffusion(const double D0, const double Te, const double Tp){
   return apply_thermal_effects(D0, Te, Tp);
}

static inline double scale_lambda_sdl(const double lambda0, const double Te, const double Tp){
   return apply_thermal_effects(lambda0, Te, Tp);
}

static inline double scale_tau_s(const double tau0, const double Te, const double Tp){
   return apply_thermal_effects(tau0, Te, Tp);
}

static inline double scale_seebeck(const double S0, const double Te, const double Tp){
   return apply_thermal_effects(S0, Te, Tp);
}

//------------------------------------------------------------------------------
// Superdiffusive transport: Enhanced hot electron velocity calculation
// Based on Option 1 from SUPERDIFFUSION_MODEL_ANALYSIS.md
// Uses real temperature gradients from ltmp (not exponential approximation)
//
// Physical distinction from tau_s in Eq. 2:
//  - tau_s (Eq. 2): Population relaxation time - controls decay of ns (hot electron density)
//    This is the rate at which hot electrons thermalize and are removed from the ns population
//  - tau_e (here): Energy relaxation time - controls decay of velocity (energy loss)
//    This is the rate at which hot electrons lose energy and slow down, even while still "hot"
//
// These are related but distinct:
//  - tau_e ~ 50-200 fs: Fast energy loss (electrons cool down quickly)
//  - tau_s ~ 200 fs - ps: Slower population decay (electrons eventually thermalize/absorbed)
//  - Typically tau_e < tau_s: Electrons lose energy faster than they disappear
//
// In the superdiffusive current Jns = ve*ns:
//  - ns decays with tau_s (population relaxation in Eq. 2)
//  - ve decays with tau_e (energy relaxation, velocity decrease)
//  - Both effects contribute to the time evolution of Jns
//------------------------------------------------------------------------------
static double superdiffusive_speed_factor(const double Te){
   if(!sc1d_superdiffusive_enable) return 1.0;

   const double Te_ref_safe = std::max(sc1d_reference_temperature, 1.0);
   const double Te_ratio = std::max(0.01, Te / Te_ref_safe);
   const double temp_fac = std::sqrt(Te_ratio);

   // tau_e is sc1d_tau_e: the internal default, or spin-currents-1d-tau-e.
   // It is not the laser-electrical hot-electron lifetime.
   // Decay starts at the ltmp envelope peak, on the drive-phase clock.
   double time_fac = 1.0;
   if(sc1d_drive_active()){
      const double t_from_peak = drive_time_s() - ltmp::laser_peak_time();
      if(t_from_peak > 0.0){
         const double tau_e = std::max(sc1d_tau_e, 1e-18);
         time_fac = std::exp(-t_from_peak / tau_e);
      }
   }
   const double fac = temp_fac * time_fac;
   return std::max(fac, 0.1);
}

static double signed_hot_electron_velocity(const double De,
                                           const double Te,
                                           const double dTe_dz,
                                           const bool have_Te){
   double ve = 0.0;
   // Superdiffusion requires ltmp, so the scale is the electron-temperature
   // gradient on that grid: ve = -De (dTe/dz) / Te.
   if(have_Te && Te > 1.0){
      ve = -De * dTe_dz / Te;
   } else if(!sc1d_superdiffusive_enable){
      const double dopt = std::max(sc1d_optical_absorption_length, 1e-18);
      const double vmag = (sc1d_v_e0 > 0.0) ? std::fabs(sc1d_v_e0) : (De / dopt);
      ve = -vmag;
   }
   ve *= superdiffusive_speed_factor(Te);
   return ve;
}

static double einstein_conductivity(const double D, const double Te){
   const double T = std::max(Te, 1.0);
   return (kMetalN * kE * kE * std::max(D, 0.0)) / (kBoltzmann * T);
}

// Matches the constant used in spinaccumulation.cpp
static constexpr double kAtomcellVolume = 2.89e-30; // m^3
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

static double lerp(const double a, const double b, const double w){
   return (1.0 - w)*a + w*b;
}

static double microcell_interpolation_weight(const int j,
                                             const int nsub,
                                             const int k,
                                             const int ncz,
                                             const bool is_interface){
   if(kInterpolateBetweenMicrocells && k < ncz-1 && !is_interface){
      return (j + 0.5) / (double)nsub;
   }
   return 0.0;
}

static double edge_harmonic_mean(const double left, const double right){
   if(left > 0.0 && right > 0.0){
      return 2.0*left*right / (left + right);
   }
   return 0.0;
}

static double edge_seebeck_coefficient(const double SL, const double SR){
   if(SL != 0.0 && SR != 0.0 && SL*SR < 0.0){
      return 0.5*(SL + SR);
   }
   if(SL != 0.0 && SR != 0.0){
      const double sum = SL + SR;
      if(std::fabs(sum) > 1e-12){
         return 2.0*SL*SR / sum;
      }
      return 0.5*(SL + SR);
   }
   if(SL != 0.0) return SL;
   if(SR != 0.0) return SR;
   return 0.0;
}

//------------------------------------------------------------------------------
// Fine-grid 1D charge solver: ns, ne, V, Jc
// Screened piece: σ E + spin-diffusion charge, Poisson ∂z(ε E) = ne.
// Serban Eq. (1): Jc also carries −σ S ∇Te + J_C^(S) − De ∇ne, with
// J_C^(S) = ve ns − De ∇ns. That part uses the wall condition Jc = 0
// and is added after the screened solve, so Eq. (4) still sees it.
//------------------------------------------------------------------------------

struct Block2 {
   double a00;
   double a01;
   double a10;
   double a11;
};

struct Vec2 {
   double x;
   double y;
};

static Block2 block2_zero(){
   Block2 A;
   A.a00 = 0.0; A.a01 = 0.0; A.a10 = 0.0; A.a11 = 0.0;
   return A;
}

static Block2 block2_mul(const Block2& A, const Block2& B){
   Block2 C;
   C.a00 = A.a00*B.a00 + A.a01*B.a10;
   C.a01 = A.a00*B.a01 + A.a01*B.a11;
   C.a10 = A.a10*B.a00 + A.a11*B.a10;
   C.a11 = A.a10*B.a01 + A.a11*B.a11;
   return C;
}

static Vec2 block2_mul_vec(const Block2& A, const Vec2& v){
   Vec2 r;
   r.x = A.a00*v.x + A.a01*v.y;
   r.y = A.a10*v.x + A.a11*v.y;
   return r;
}

static Block2 block2_inv(const Block2& A){
   const double det = A.a00*A.a11 - A.a01*A.a10;
   double d = det;
   if(std::fabs(d) < kEps) d = (d >= 0.0 ? kEps : -kEps);
   const double idet = 1.0 / d;
   Block2 I;
   I.a00 =  A.a11 * idet;
   I.a01 = -A.a01 * idet;
   I.a10 = -A.a10 * idet;
   I.a11 =  A.a00 * idet;
   return I;
}

static Block2 block2_sub(const Block2& A, const Block2& B){
   Block2 C;
   C.a00 = A.a00 - B.a00;
   C.a01 = A.a01 - B.a01;
   C.a10 = A.a10 - B.a10;
   C.a11 = A.a11 - B.a11;
   return C;
}

static Vec2 vec2_sub(const Vec2& a, const Vec2& b){
   Vec2 r;
   r.x = a.x - b.x;
   r.y = a.y - b.y;
   return r;
}

static void block_thomas(const std::vector<Block2>& L,
                         const std::vector<Block2>& D,
                         const std::vector<Block2>& U,
                         std::vector<Vec2>& rhs){
   const int n = static_cast<int>(rhs.size());
   std::vector<Block2> Cp(n, block2_zero());
   std::vector<Vec2> dp(n);

   Block2 Dinv = block2_inv(D[0]);
   Cp[0] = block2_mul(Dinv, U[0]);
   dp[0] = block2_mul_vec(Dinv, rhs[0]);

   for(int i=1;i<n;i++){
      const Block2 LCp = block2_mul(L[i], Cp[i-1]);
      const Block2 Di = block2_sub(D[i], LCp);
      Dinv = block2_inv(Di);
      if(i < n-1){
         Cp[i] = block2_mul(Dinv, U[i]);
      } else {
         Cp[i] = block2_zero();
      }
      const Vec2 Ldp = block2_mul_vec(L[i], dp[i-1]);
      const Vec2 rhs_i = vec2_sub(rhs[i], Ldp);
      dp[i] = block2_mul_vec(Dinv, rhs_i);
   }

   rhs[n-1] = dp[n-1];
   for(int i=n-2;i>=0;i--){
      const Vec2 Cpx = block2_mul_vec(Cp[i], rhs[i+1]);
      rhs[i].x = dp[i].x - Cpx.x;
      rhs[i].y = dp[i].y - Cpx.y;
   }
}

static void fill_fine_charge_transport(const int nf,
                                       const double dz,
                                       const double t_s,
                                       const double* D_cell,
                                       const double* sigma0_cell,
                                       double* Te_cell,
                                       double* De_cell,
                                       double* ve_cell,
                                       double* sigma_cell){
   (void)t_s;
   const bool have_Te = sc1d_has_electron_temperature();
   stopwatch_t sw_te;
   sw_te.start();
   for(int i=0;i<nf;i++){
      const double zc = (i + 0.5) * dz;
      const double Te = get_local_electron_temperature(zc);
      const double Tp = get_local_phonon_temperature(zc);
      Te_cell[i] = Te;
      const double D0 = std::max(0.0, D_cell[i]);
      De_cell[i] = scale_diffusion(D0, Te, Tp);
      // A negative mat value is unset: Einstein σ from De and Te.
      // An explicit 0 is an insulator and stays 0. A positive value is used as written.
      double sigma = (sigma0_cell[i] < 0.0) ? einstein_conductivity(De_cell[i], Te) : sigma0_cell[i];
      sigma_cell[i] = apply_thermal_effects(sigma, Te, Tp);
   }
   sc1d_time_interp += sw_te.elapsed_seconds();

   for(int i=0;i<nf;i++){
      double dTe_dz = 0.0;
      if(nf == 1){
         dTe_dz = 0.0;
      } else if(i == 0){
         dTe_dz = (Te_cell[1] - Te_cell[0]) / dz;
      } else if(i == nf-1){
         dTe_dz = (Te_cell[nf-1] - Te_cell[nf-2]) / dz;
      } else {
         dTe_dz = (Te_cell[i+1] - Te_cell[i-1]) / (2.0 * dz);
      }
      ve_cell[i] = signed_hot_electron_velocity(De_cell[i], Te_cell[i], dTe_dz, have_Te);
   }
}

static void fill_fine_charge_faces(const int nf,
                                   const double* De_cell,
                                   const double* ve_cell,
                                   const double* sigma_cell,
                                   double* De_edge,
                                   double* ve_edge,
                                   double* sigma_edge,
                                   double* eps_edge){
   for(int e=0;e<=nf;e++){
      De_edge[e] = 0.0;
      ve_edge[e] = 0.0;
      sigma_edge[e] = 0.0;
      eps_edge[e] = 0.0;
   }
   for(int e=1;e<nf;e++){
      De_edge[e] = edge_harmonic_mean(De_cell[e-1], De_cell[e]);
      sigma_edge[e] = edge_harmonic_mean(sigma_cell[e-1], sigma_cell[e]);
      ve_edge[e] = 0.5*(ve_cell[e-1] + ve_cell[e]);
      eps_edge[e] = kEps0;
   }
}

static void assemble_ns_fine_tridiagonal(const int nf,
                                         const double dt,
                                         const double dz,
                                         const double t_s,
                                         const double* De_edge,
                                         const double* ve_edge,
                                         const double* ns,
                                         const double* Te_cell,
                                         std::vector<double>& a,
                                         std::vector<double>& b,
                                         std::vector<double>& c,
                                         std::vector<double>& rhs){
   (void)t_s;
   double Te_avg = 0.0;
   double Tp_avg = 0.0;
   for(int i=0;i<nf;i++){
      Te_avg += Te_cell[i];
      Tp_avg += get_local_phonon_temperature((i + 0.5) * dz);
   }
   Te_avg /= std::max(1, nf);
   Tp_avg /= std::max(1, nf);
   const double tau_s = scale_tau_s(std::max(sc1d_tau_s, 0.0), Te_avg, Tp_avg);
   const double inv_tau = (tau_s > 0.0) ? (1.0/tau_s) : 0.0;
   const double idz = dt / dz;
   const double idz2 = dt / (dz*dz);

   for(int i=0;i<nf;i++){
      a[i] = 0.0;
      b[i] = 1.0 + dt*inv_tau;
      c[i] = 0.0;

      const double DeL = De_edge[i];
      const double DeR = De_edge[i+1];
      a[i] -= idz2 * DeL;
      b[i] += idz2 * (DeL + DeR);
      c[i] -= idz2 * DeR;

      const double veL = ve_edge[i];
      const double veR = ve_edge[i+1];
      if(i > 0){
         if(veL >= 0.0){
            a[i] -= idz * veL;
         } else {
            b[i] -= idz * veL;
         }
      }
      if(i < nf-1){
         if(veR >= 0.0){
            b[i] += idz * veR;
         } else {
            c[i] += idz * veR;
         }
      }

      const double zc = (i + 0.5) * dz;
      const double Q = ns_source_power_density(zc);
      const double src = (kE * sc1d_laser_eta * Q * sc1d_laser_wavelength) / (kPlanck * kC0);
      rhs[i] = ns[i] + dt*src;
   }
}

static void build_jns_fine(const int nf,
                           const double dz,
                           const double* De_edge,
                           const double* ve_edge,
                           const double* ns,
                           double* Jns_edge){
   Jns_edge[0] = 0.0;
   Jns_edge[nf] = 0.0;
   for(int e=1;e<nf;e++){
      const double ve = ve_edge[e];
      const double ns_up = (ve >= 0.0) ? ns[e-1] : ns[e];
      const double adv = ve * ns_up;
      const double diff = -De_edge[e] * (ns[e] - ns[e-1]) / dz;
      Jns_edge[e] = adv + diff;
   }
}

static void build_jknown_fine(const int nf,
                              const double dz,
                              const double* De_edge,
                              const double* Bd_cell,
                              const double* S,
                              const double* mhat,
                              const double* Jns_edge,
                              double* Jknown){
   // Serban Eq. (1) puts J_C^(S) = −De[(∇Te/Te) ns + ∇ns] inside Jc, and
   // Eq. (4) multiplies that full Jc by P. The screened solve would cancel
   // it with σE in one step, so Jns is integrated with the open-circuit ne
   // instead. Only the spin-diffusion charge stays here.
   (void)Jns_edge;
   Jknown[0] = 0.0;
   Jknown[nf] = 0.0;
   for(int e=1;e<nf;e++){
      const double Bd_f = 0.5*(Bd_cell[e-1] + Bd_cell[e]);
      const double sL = S[3*(e-1)+0]*mhat[3*(e-1)+0] + S[3*(e-1)+1]*mhat[3*(e-1)+1] + S[3*(e-1)+2]*mhat[3*(e-1)+2];
      const double sR = S[3*e+0]*mhat[3*e+0] + S[3*e+1]*mhat[3*e+1] + S[3*e+2]*mhat[3*e+2];
      const double Jspin = (De_edge[e] * Bd_f * kE / kMuB) * (sR - sL) / dz;
      Jknown[e] = Jspin;
   }
}

static void solve_ne_V_implicit(const int nf,
                                const double dt,
                                const double dz,
                                const double* De_edge,
                                const double* sigma_edge,
                                const double* eps_edge,
                                const double* Jknown,
                                const double* ne,
                                double* ne_out,
                                double* V_out){
   std::vector<Block2> L(nf, block2_zero());
   std::vector<Block2> D(nf, block2_zero());
   std::vector<Block2> U(nf, block2_zero());
   std::vector<Vec2> rhs(nf);
   const double idz = dt / dz;

   for(int i=0;i<nf;i++){
      const double hL = De_edge[i] / dz;
      const double hR = De_edge[i+1] / dz;
      const double gL = sigma_edge[i] / dz;
      const double gR = sigma_edge[i+1] / dz;
      const double aeL = eps_edge[i] / (dz*dz);
      const double aeR = eps_edge[i+1] / (dz*dz);

      D[i].a00 = 1.0;
      if(i > 0){
         D[i].a00 += idz * hL;
         L[i].a00 = -idz * hL;
         D[i].a01 += idz * gL;
         L[i].a01 = -idz * gL;
      }
      if(i < nf-1){
         D[i].a00 += idz * hR;
         U[i].a00 = -idz * hR;
         D[i].a01 += idz * gR;
         U[i].a01 = -idz * gR;
      }

      const double JkR = Jknown[i+1];
      const double JkL = Jknown[i];
      rhs[i].x = ne[i] - idz * (JkR - JkL);

      if(i == 0){
         D[i].a10 = 0.0;
         D[i].a11 = 1.0;
         L[i].a10 = 0.0;
         L[i].a11 = 0.0;
         U[i].a10 = 0.0;
         U[i].a11 = 0.0;
         rhs[i].y = 0.0;
      } else {
         D[i].a10 = -1.0;
         D[i].a11 = aeL + aeR;
         L[i].a10 = 0.0;
         L[i].a11 = -aeL;
         U[i].a10 = 0.0;
         U[i].a11 = -aeR;
         rhs[i].y = 0.0;
      }
   }

   block_thomas(L, D, U, rhs);

   for(int i=0;i<nf;i++){
      ne_out[i] = rhs[i].x;
      V_out[i] = rhs[i].y;
   }
   V_out[0] = 0.0;
}

static void build_jc_fine(const int nf,
                          const double dz,
                          const double* De_edge,
                          const double* sigma_edge,
                          const double* Jknown,
                          const double* ne,
                          const double* V,
                          double* Jc_edge){
   Jc_edge[0] = 0.0;
   Jc_edge[nf] = 0.0;
   for(int e=1;e<nf;e++){
      const double E = -(V[e] - V[e-1]) / dz;
      const double diff_ne = -De_edge[e] * (ne[e] - ne[e-1]) / dz;
      Jc_edge[e] = sigma_edge[e] * E + diff_ne + Jknown[e];
   }
}

// Steady open circuit on the fine grid. Both walls have Jc = 0, so in one
// dimension ∂z Jc = 0 forces Jc = 0 on every face:
//   σ E - De ∂z ne + J_known = 0,   ∂z(ε E) = ne,   V_0 = 0.
// J_known (spin diffusion, and hot electrons when that drive is on) is
// cancelled by E inside this solve. Seebeck is integrated separately.
// The matrix is the face balance for
// V, with ne eliminated by the discrete Poisson law, and is pentadiagonal.
// The previous block-Thomas step kept σE and J_known as separate 1e14 terms
// and leaked ~1e-4 of J_known into Jc.
static void solve_screened_open_circuit(const int nf,
                                        const double dz,
                                        const double* De_edge,
                                        const double* sigma_edge,
                                        const double* eps_edge,
                                        const double* Jknown,
                                        double* ne,
                                        double* V,
                                        double* Jc_edge){
   if(nf <= 0) return;
   V[0] = 0.0;
   if(nf == 1){
      ne[0] = 0.0;
      Jc_edge[0] = 0.0;
      Jc_edge[1] = 0.0;
      return;
   }

   const int m = nf - 1;
   const double idz2 = 1.0 / (dz * dz);
   std::vector<double> al(m, 0.0), bl(m, 0.0), diag(m, 0.0), au(m, 0.0), au2(m, 0.0), rhs(m, 0.0);

   auto add_coeff = [&](const int row, const int vindex, const double coef){
      if(row < 0 || row >= m || vindex <= 0 || vindex >= nf) return;
      const int delta = (vindex - 1) - row;
      if(delta == -2) al[row] += coef;
      else if(delta == -1) bl[row] += coef;
      else if(delta == 0) diag[row] += coef;
      else if(delta == 1) au[row] += coef;
      else if(delta == 2) au2[row] += coef;
   };
   auto add_ne = [&](const int row, const int cell, const double scale){
      const double aL = eps_edge[cell] * idz2;
      const double aR = eps_edge[cell + 1] * idz2;
      add_coeff(row, cell - 1, scale * (-aL));
      add_coeff(row, cell,     scale * (aL + aR));
      add_coeff(row, cell + 1, scale * (-aR));
   };

   for(int face = 1; face < nf; ++face){
      const int row = face - 1;
      const double sig = sigma_edge[face];
      if(!(sig > 0.0)){
         add_coeff(row, face, 1.0);
         add_coeff(row, face - 1, -1.0);
         rhs[row] = 0.0;
         continue;
      }
      const double r = De_edge[face] / sig;
      add_coeff(row, face - 1, 1.0);
      add_coeff(row, face, -1.0);
      add_ne(row, face, -r);
      add_ne(row, face - 1, r);
      rhs[row] = -Jknown[face] * dz / sig;
   }

   for(int i = 0; i < m; ++i){
      if(i >= 2){
         const double piv = diag[i - 2];
         const double w = (std::fabs(piv) > 0.0) ? al[i] / piv : 0.0;
         bl[i] -= w * au[i - 2];
         diag[i] -= w * au2[i - 2];
         rhs[i] -= w * rhs[i - 2];
      }
      if(i >= 1){
         const double piv = diag[i - 1];
         const double w = (std::fabs(piv) > 0.0) ? bl[i] / piv : 0.0;
         diag[i] -= w * au[i - 1];
         au[i] -= w * au2[i - 1];
         rhs[i] -= w * rhs[i - 1];
      }
   }

   std::vector<double> x(m, 0.0);
   for(int i = m - 1; i >= 0; --i){
      double s = rhs[i];
      if(i + 1 < m) s -= au[i] * x[i + 1];
      if(i + 2 < m) s -= au2[i] * x[i + 2];
      double piv = diag[i];
      if(std::fabs(piv) < 1e-30) piv = (piv >= 0.0 ? 1e-30 : -1e-30);
      x[i] = s / piv;
   }

   for(int i = 1; i < nf; ++i) V[i] = x[i - 1];
   for(int i = 0; i < nf; ++i){
      const double aL = eps_edge[i] * idz2;
      const double aR = eps_edge[i + 1] * idz2;
      const double Vm = (i > 0) ? V[i - 1] : 0.0;
      if(i == nf - 1) ne[i] = aL * (V[i] - Vm);
      else ne[i] = aL * (V[i] - Vm) + aR * (V[i] - V[i + 1]);
   }

   Jc_edge[0] = 0.0;
   Jc_edge[nf] = 0.0;
   for(int e = 1; e < nf; ++e){
      const double E = -(V[e] - V[e - 1]) / dz;
      const double diff_ne = -De_edge[e] * (ne[e] - ne[e - 1]) / dz;
      const double Jc = sigma_edge[e] * E + diff_ne + Jknown[e];
      Jc_edge[e] = std::isfinite(Jc) ? Jc : 0.0;
   }
}

static void solve_open_circuit_charge_steady(double* ns,
                                             double* ne,
                                             double* V,
                                             double* Jc_edge,
                                             const int nf,
                                             const double dz,
                                             const double t_s,
                                             const double* D_cell,
                                             const double* sigma0_cell,
                                             const double* Bd_cell,
                                             const double* seebeck_cell,
                                             const double* S,
                                             const double* mhat){
   for(int i=0;i<nf;i++) ns[i] = 0.0;

   std::vector<double> Te_cell(nf, sc1d_reference_temperature);
   std::vector<double> De_cell(nf, 0.0);
   std::vector<double> ve_cell(nf, 0.0);
   std::vector<double> sigma_cell(nf, 0.0);
   fill_fine_charge_transport(nf, dz, t_s, D_cell, sigma0_cell,
                              Te_cell.data(), De_cell.data(), ve_cell.data(), sigma_cell.data());

   std::vector<double> De_edge(nf+1, 0.0);
   std::vector<double> ve_edge(nf+1, 0.0);
   std::vector<double> sigma_edge(nf+1, 0.0);
   std::vector<double> eps_edge(nf+1, 0.0);
   fill_fine_charge_faces(nf, De_cell.data(), ve_cell.data(), sigma_cell.data(),
                          De_edge.data(), ve_edge.data(), sigma_edge.data(), eps_edge.data());

   std::vector<double> Jns_edge(nf+1, 0.0);
   std::vector<double> Jknown(nf+1, 0.0);
   (void)seebeck_cell;
   build_jknown_fine(nf, dz, De_edge.data(), Bd_cell,
                     S, mhat, Jns_edge.data(), Jknown.data());
   solve_screened_open_circuit(nf, dz, De_edge.data(), sigma_edge.data(), eps_edge.data(),
                               Jknown.data(), ne, V, Jc_edge);
}

// Serban Eq. (1) and (3). One excess density answers both drives:
//   Jc = −σ S ∇Te − Jns − De ∇ne,   ∂ne/∂t = −∂z Jc,   Jc = 0 on the walls.
// Jns = ve ns − De ∇ns is the hot-electron flux in charge-magnitude units.
// The source fills ns with +e times the number density, so the carriers'
// charge current is −Jns. Positive Jc is positive charge toward z = max.
// Hot electrons flow toward the cold substrate (ve < 0), which is Jc > 0.
// Eq. (4) then takes this Jc. S is the absolute thermopower. Co and Pt are
// negative, so −S ∇Te points the same way. z = 0 is the bottom.
// Backward Euler. The screened solve does not see this.
static void step_seebeck_open_circuit(const int nf,
                                      const double dt,
                                      const double dz,
                                      const double* De_edge,
                                      const double* sigma_edge,
                                      const double* Te_cell,
                                      const double* seebeck_cell,
                                      const double* Jns_edge,
                                      double* ne_seebeck,
                                      double* J_seebeck){
   std::vector<double> Jdrive(nf + 1, 0.0);
   if(sc1d_seebeck_enable && sc1d_drive_active() && nf > 1 && dz > 0.0){
      for(int e = 1; e < nf; ++e){
         const double TeL = Te_cell[e - 1];
         const double TeR = Te_cell[e];
         if(!(std::isfinite(TeL) && std::isfinite(TeR) && TeL > 0.0 && TeR > 0.0)) continue;
         const double zL = (e - 0.5) * dz;
         const double zR = (e + 0.5) * dz;
         const double TpL = get_local_phonon_temperature(zL);
         const double TpR = get_local_phonon_temperature(zR);
         const double SL = scale_seebeck(seebeck_cell[e - 1], TeL, TpL);
         const double SR = scale_seebeck(seebeck_cell[e], TeR, TpR);
         const double S_avg = edge_seebeck_coefficient(SL, SR);
         const double sig = sigma_edge[e];
         if(!(sig > 0.0) || !std::isfinite(S_avg)) continue;
         Jdrive[e] = -sig * S_avg * (TeR - TeL) / dz;
      }
   }
   if(sc1d_superdiffusive_enable && sc1d_drive_active() && Jns_edge != nullptr && nf > 1){
      for(int e = 1; e < nf; ++e){
         if(std::isfinite(Jns_edge[e])) Jdrive[e] -= Jns_edge[e];
      }
   }

   J_seebeck[0] = 0.0;
   if(nf <= 1 || dz <= 0.0 || dt <= 0.0){
      if(nf >= 0) J_seebeck[nf] = 0.0;
      return;
   }

   std::vector<double> a(nf, 0.0);
   std::vector<double> b(nf, 0.0);
   std::vector<double> c(nf, 0.0);
   std::vector<double> rhs(nf, 0.0);
   const double idz = dt / dz;
   const double idz2 = dt / (dz * dz);
   auto alpha = [&](const int e){
      if(e <= 0 || e >= nf) return 0.0;
      return idz2 * std::max(De_edge[e], 0.0);
   };

   {
      const double aR = alpha(1);
      b[0] = 1.0 + aR;
      c[0] = -aR;
      rhs[0] = ne_seebeck[0] - idz * Jdrive[1];
   }
   for(int i = 1; i < nf - 1; ++i){
      const double aL = alpha(i);
      const double aR = alpha(i + 1);
      a[i] = -aL;
      b[i] = 1.0 + aL + aR;
      c[i] = -aR;
      rhs[i] = ne_seebeck[i] - idz * (Jdrive[i + 1] - Jdrive[i]);
   }
   {
      const int i = nf - 1;
      const double aL = alpha(nf - 1);
      a[i] = -aL;
      b[i] = 1.0 + aL;
      rhs[i] = ne_seebeck[i] + idz * Jdrive[nf - 1];
   }

   thomas_solve(a, b, c, rhs);
   for(int i = 0; i < nf; ++i){
      ne_seebeck[i] = std::isfinite(rhs[i]) ? rhs[i] : 0.0;
   }

   J_seebeck[nf] = 0.0;
   for(int e = 1; e < nf; ++e){
      const double back = -std::max(De_edge[e], 0.0) * (ne_seebeck[e] - ne_seebeck[e - 1]) / dz;
      const double J = Jdrive[e] + back;
      J_seebeck[e] = std::isfinite(J) ? J : 0.0;
   }
}

static void update_charge_transients_1d_stack(double* ns,
                                              double* ne,
                                              double* V,
                                              double* Jc_edge,
                                              double* ne_seebeck,
                                              const int nf,
                                              const double dz,
                                              const double dt,
                                              const double t_s,
                                              const double* D_cell,
                                              const double* sigma0_cell,
                                              const double* Bd_cell,
                                              const double* seebeck_cell,
                                              const double* S,
                                              const double* mhat){
   std::vector<double> Te_cell(nf, sc1d_reference_temperature);
   std::vector<double> De_cell(nf, 0.0);
   std::vector<double> ve_cell(nf, 0.0);
   std::vector<double> sigma_cell(nf, 0.0);
   fill_fine_charge_transport(nf, dz, t_s, D_cell, sigma0_cell,
                              Te_cell.data(), De_cell.data(), ve_cell.data(), sigma_cell.data());

   std::vector<double> De_edge(nf+1, 0.0);
   std::vector<double> ve_edge(nf+1, 0.0);
   std::vector<double> sigma_edge(nf+1, 0.0);
   std::vector<double> eps_edge(nf+1, 0.0);
   fill_fine_charge_faces(nf, De_cell.data(), ve_cell.data(), sigma_cell.data(),
                          De_edge.data(), ve_edge.data(), sigma_edge.data(), eps_edge.data());

   std::vector<double> a(nf, 0.0);
   std::vector<double> b(nf, 0.0);
   std::vector<double> c(nf, 0.0);
   std::vector<double> rhs(nf, 0.0);
   assemble_ns_fine_tridiagonal(nf, dt, dz, t_s, De_edge.data(), ve_edge.data(),
                                ns, Te_cell.data(), a, b, c, rhs);
   thomas_solve(a, b, c, rhs);
   for(int i=0;i<nf;i++) ns[i] = rhs[i];

   std::vector<double> Jns_edge(nf+1, 0.0);
   build_jns_fine(nf, dz, De_edge.data(), ve_edge.data(), ns, Jns_edge.data());

   std::vector<double> Jknown(nf+1, 0.0);
   build_jknown_fine(nf, dz, De_edge.data(), Bd_cell,
                     S, mhat, Jns_edge.data(), Jknown.data());

   solve_screened_open_circuit(nf, dz, De_edge.data(), sigma_edge.data(), eps_edge.data(),
                               Jknown.data(), ne, V, Jc_edge);

   if(ne_seebeck != nullptr && (sc1d_seebeck_enable || sc1d_superdiffusive_enable)){
      std::vector<double> Jopen(nf + 1, 0.0);
      step_seebeck_open_circuit(nf, dt, dz, De_edge.data(), sigma_edge.data(),
                                Te_cell.data(), seebeck_cell, Jns_edge.data(),
                                ne_seebeck, Jopen.data());
      for(int e = 0; e <= nf; ++e) Jc_edge[e] += Jopen[e];
      for(int i = 0; i < nf; ++i) ne[i] += ne_seebeck[i];
   }
}

//------------------------------------------------------------------------------
// Initialisation helpers
//------------------------------------------------------------------------------

static int compute_fine_subdivisions(){
   const double dzA = std::max(0.01, sc1d_fine_dz);
   const double DzA = std::max(0.01, micro_cell_thickness);
   int nsub = static_cast<int>(std::ceil(DzA / dzA));
   nsub = std::max(1, std::min(128, nsub));
   return nsub;
}

static int column_owner_rank(const int stack){
#ifdef MPICF
   const int nrank = std::max(1, vmpi::num_processors);
   return stack % nrank;
#else
   (void)stack;
   return 0;
#endif
}

static bool rank_owns_column(const int stack){
#ifdef MPICF
   return column_owner_rank(stack) == vmpi::my_rank;
#else
   (void)stack;
   return true;
#endif
}

static void collect_local_stacks(){
   // Column s belongs to rank (s % nprocs): the first rank in rank order.
   // Spare ranks own nothing. Extra columns wrap, so one rank can own several.
   // The atom domain is not this split. m is already summed, and the ltmp
   // profile is the same on every rank, so the owner can step the column
   // without holding its atoms. The cell fields are broadcast afterwards.
   sc1d_local_stacks.clear();
   sc1d_stack_local_index.assign(num_stacks_y, -1);
   for(int s=0; s<num_stacks_y; ++s){
      if(!rank_owns_column(s)) continue;
      sc1d_stack_local_index[s] = (int)sc1d_local_stacks.size();
      sc1d_local_stacks.push_back(s);
   }
#ifdef MPICF
   if(vmpi::my_rank == 0)
#endif
   {
      zlog << zTs() << "1D spin-current columns " << num_stacks_y
           << ". Rank r owns columns r, r+nranks, ...; spare ranks own none."
           << std::endl;
   }
}

static void allocate_sc1d_arrays(){
   // Every rank stores every column. Only the owner steps a column; the
   // unused slices stay zero until the field broadcast fills them.
   const std::size_t nstack = (std::size_t)std::max(0, num_stacks_y);
   const std::size_t expected_size = nstack * (std::size_t)sc1d_nf * 3u;
   if(sc1d_Sfine.size() != expected_size){
      sc1d_Sfine.assign(expected_size, 0.0);
   }
   if(sc1d_k_prev.size() != expected_size){
      sc1d_k_prev.assign(expected_size, 0.0);
   }

   const int nf = sc1d_nf;
   const std::size_t nfn = nstack * (std::size_t)nf;
   sc1d_ns_fine.assign(nfn, 0.0);
   sc1d_ne_fine.assign(nfn, 0.0);
   sc1d_ne_seebeck_fine.assign(nfn, 0.0);
   sc1d_V_fine.assign(nfn, 0.0);
   sc1d_Jc_edge_fine.assign(nstack * (std::size_t)(nf+1), 0.0);
   sc1d_sigma_fine.assign(nfn, 0.0);
   sc1d_seebeck_fine.assign(nfn, 0.0);
   sc1d_mat_fine.assign(nfn, -1);
   sc1d_sa_demag_fine.assign(nfn, 0.0);

   sc1d_Bc_fine.assign(nfn, 0.0);
   sc1d_Bd_fine.assign(nfn, 0.0);
   sc1d_D_fine.assign(nfn, 0.0);
   sc1d_lsf_fine.assign(nfn, 0.0);
   sc1d_lphi_fine.assign(nfn, 0.0);
   sc1d_Jsd_fine.assign(nfn, 0.0);
   sc1d_chi_fine.assign(nfn, 0.0);
   sc1d_sa_inf_fine.assign(nfn, 0.0);
   sc1d_alpha_edge.assign(nstack * (std::size_t)std::max(0,nf-1), 0.0);
}

static void initialise_magnetization_history(){
   const int total_microcells = num_x_stacks * num_y_stacks * num_microcells_per_stack;
   sc1d_m_prev_mag.assign(total_microcells, 0.0);
   sc1d_m0_mag.assign(total_microcells, 0.0);
   sc1d_mhat_ref.assign(total_microcells*3u, 0.0);

   for(int c=0;c<total_microcells;c++){
      const double mx = m[3*c+0];
      const double my = m[3*c+1];
      const double mz = m[3*c+2];
      const double mm = std::sqrt(mx*mx + my*my + mz*mz);
      if(mm > 1e-12){
         sc1d_mhat_ref[3*c+0] = mx/mm;
         sc1d_mhat_ref[3*c+1] = my/mm;
         sc1d_mhat_ref[3*c+2] = mz/mm;
      } else {
         sc1d_mhat_ref[3*c+0] = 0.0;
         sc1d_mhat_ref[3*c+1] = 0.0;
         sc1d_mhat_ref[3*c+2] = 1.0;
      }
   }
}

// Linear interpolation of a coarse-cell field onto a fine-cell centre.
// Coarse values sit at cell centres. Below the first centre and above the
// last centre the end cell is used. Robin resistance is separate and sits
// on the coarse face.
static double lerp_coarse_centre(const double* field,
                                 const int start_cell,
                                 const int ncz,
                                 const double s){
   if(ncz <= 1 || s <= 0.5) return field[start_cell];
   const double s_top = (double)ncz - 0.5;
   if(s >= s_top) return field[start_cell + ncz - 1];
   const int k = static_cast<int>(std::floor(s - 0.5));
   const double w = s - ((double)k + 0.5);
   return (1.0 - w)*field[start_cell + k] + w*field[start_cell + k + 1];
}

// Fine-cell constants are the linear interpolation of the coarse averages
// between neighbouring coarse-cell centres.
static void prolong_coarse_constants(const int start_cell,
                                     const int nsub,
                                     const int nf,
                                     double* Bc_fine,
                                     double* Bd_fine,
                                     double* D_fine,
                                     double* lsf_fine,
                                     double* lphi_fine,
                                     double* Jsd_fine,
                                     double* chi_fine,
                                     double* sa_inf_fine,
                                     double* sigma_fine,
                                     double* seebeck_fine){
   const int nsub_safe = std::max(1, nsub);
   const int ncz = std::max(1, num_microcells_per_stack);
   for(int i=0;i<nf;i++){
      const double s = ((double)i + 0.5) / (double)nsub_safe;
      Bc_fine[i] = lerp_coarse_centre(beta_cond.data(), start_cell, ncz, s);
      Bd_fine[i] = lerp_coarse_centre(beta_diff.data(), start_cell, ncz, s);
      D_fine[i] = lerp_coarse_centre(diffusion.data(), start_cell, ncz, s);
      lsf_fine[i] = lerp_coarse_centre(lambda_sdl.data(), start_cell, ncz, s);
      lphi_fine[i] = lerp_coarse_centre(lambda_phi.data(), start_cell, ncz, s);
      Jsd_fine[i] = lerp_coarse_centre(sd_exchange.data(), start_cell, ncz, s);
      chi_fine[i] = lerp_coarse_centre(chi_demag.data(), start_cell, ncz, s);
      sa_inf_fine[i] = lerp_coarse_centre(sa_infinity.data(), start_cell, ncz, s);
      // Conductivity is not blended. A negative sentinel means Einstein and a
      // 0 is an insulator; interpolating the two would conduct part of the silica.
      const int k_own = std::max(0, std::min(ncz - 1, static_cast<int>(std::floor(s))));
      sigma_fine[i] = conductivity[start_cell + k_own];
      seebeck_fine[i] = lerp_coarse_centre(seebeck_coefficient.data(), start_cell, ncz, s);
   }
}

static void assemble_alpha_edges_from_coarse(const int start_cell,
                                             const int nsub,
                                             const int nf,
                                             const double dz,
                                             const double* D_fine,
                                             double* alpha_edge_stack){
   if(nf < 2) return;
   const int nsub_safe = std::max(1, nsub);
   for(int i=0;i<nf-1;i++){
      double Rint = 0.0;
      // The face between the last subcell of coarse cell k and the first of k+1.
      if(((i + 1) % nsub_safe) == 0){
         const int cell = start_cell + (i / nsub_safe);
         if(cell >= 0 && cell < (int)r_int_edge.size()) Rint = r_int_edge[cell];
      }
      const double Di = D_fine[i];
      const double Dj = D_fine[i+1];
      if(Di < kEps || Dj < kEps){
         alpha_edge_stack[i] = 0.0;
         continue;
      }
      const double Rseries = dz/(2.0*Di) + Rint + dz/(2.0*Dj);
      if(Rseries <= 0.0){
         alpha_edge_stack[i] = 0.0;
         continue;
      }
      alpha_edge_stack[i] = (1.0 / Rseries) / dz;
   }
}

static void initialise_fine_spin_if_empty(const int start_cell,
                                         const int nsub,
                                         const int nf,
                                         const double dz,
                                         const double* D_fine,
                                         const double* lsf_fine,
                                         const double* alpha_edge,
                                         const double* sa_inf_fine,
                                         double* Sbase){
   bool needs_init = true;
   for(int i=0; i<nf*3; ++i){
      if(std::abs(Sbase[i]) > 1e-20){
         needs_init = false;
         break;
      }
   }
   if(!needs_init) return;

   std::vector<double> sa_eq(3*(std::size_t)nf, 0.0);
   std::vector<double> relax(nf, 0.0);
   for(int i=0;i<nf;i++){
      const int k = i / std::max(1, nsub);
      const int cell = start_cell + k;
      const double hx = sc1d_mhat_ref[3*cell+0];
      const double hy = sc1d_mhat_ref[3*cell+1];
      const double hz = sc1d_mhat_ref[3*cell+2];
      const double sa = sa_inf_fine[i];
      sa_eq[3*i+0] = sa * hx;
      sa_eq[3*i+1] = sa * hy;
      sa_eq[3*i+2] = sa * hz;
      const double lambda = lsf_fine[i];
      if(lambda > 0.0 && D_fine[i] > 0.0){
         relax[i] = D_fine[i] / (lambda * lambda);
      }
   }

   // Steady diffusion-relaxation, Jc = 0: d/dz (De dS/dz) = De (S - sa_eq) / λsf².
   // sa_eq is sa_inf along the local m-hat. sa_inf is the linear blend of the
   // coarse averages. Robin mixing sits on the coarse face.
   // The solve spreads those steps into a continuous tail, so Js = -De dS/dz
   // is already the gradient of that profile and is not zero everywhere.
   // Neumann ends.
   std::vector<double> a(nf, 0.0);
   std::vector<double> b(nf, 0.0);
   std::vector<double> c(nf, 0.0);
   std::vector<double> rhs(nf, 0.0);
   for(int comp=0; comp<3; ++comp){
      for(int i=0;i<nf;i++){
         const double aL = (i > 0) ? alpha_edge[i-1] : 0.0;
         const double aR = (i < nf-1) ? alpha_edge[i] : 0.0;
         a[i] = -aL;
         b[i] = aL + aR + relax[i];
         c[i] = -aR;
         rhs[i] = relax[i] * sa_eq[3*i+comp];
         if(std::fabs(b[i]) < kEps){
            b[i] = 1.0;
            rhs[i] = sa_eq[3*i+comp];
         }
      }
      thomas_solve(a, b, c, rhs);
      for(int i=0;i<nf;i++){
         Sbase[3*i+comp] = rhs[i];
      }
   }
}

static void build_atom_fine_map(const int nf, const double dzA){
   sc1d_atom_ls.assign(num_local_atoms, -1);
   sc1d_atom_fine_lo.assign(num_local_atoms, -1);
   sc1d_atom_fine_hi.assign(num_local_atoms, -1);

   std::vector<double> zlist;
   zlist.reserve((std::size_t)num_local_atoms);
   for(int atom=0; atom<num_local_atoms; ++atom){
      zlist.push_back(atoms::z_coord_array[atom]);
   }
   std::sort(zlist.begin(), zlist.end());
   std::vector<double> zplane;
   const double ztol = 1.0e-4;
   for(int p=0; p<(int)zlist.size(); ++p){
      if(zplane.empty() || std::fabs(zlist[p] - zplane.back()) > ztol){
         zplane.push_back(zlist[p]);
      }
   }
#ifdef MPICF
   {
      const int nrank = vmpi::num_processors;
      const int nloc = (int)zplane.size();
      std::vector<int> counts(nrank, 0);
      std::vector<int> displs(nrank, 0);
      MPI_Allgather(&nloc, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
      int nall = 0;
      for(int r=0; r<nrank; ++r){
         displs[r] = nall;
         nall += counts[r];
      }
      if(nall > 0){
         std::vector<double> allz((std::size_t)nall, 0.0);
         const double send_dummy = 0.0;
         const double* send = (nloc > 0) ? zplane.data() : &send_dummy;
         MPI_Allgatherv(send, nloc, MPI_DOUBLE, allz.data(), counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
         std::sort(allz.begin(), allz.end());
         zplane.clear();
         for(int p=0; p<nall; ++p){
            if(zplane.empty() || std::fabs(allz[p] - zplane.back()) > ztol){
               zplane.push_back(allz[p]);
            }
         }
      }
   }
#endif
   const int nplane = (int)zplane.size();

   for(int atom=0; atom<num_local_atoms; ++atom){
      const int cell = atom_st_index[atom];
      if(cell < 0 || cell >= (int)cell_stack_index.size()) continue;
      const int stack = cell_stack_index[cell] - 1;
      if(stack < 0 || stack >= num_stacks_y) continue;

      const int mat = atoms::type_array[atom];
      if(mat >= 0 && mat < (int)mp.size()){
         if(mp[mat].sd_exchange <= 0.0) continue;
      }

      const double za = atoms::z_coord_array[atom];
      int plane = 0;
      double best = 1.0e99;
      for(int p=0;p<nplane;p++){
         const double d = std::fabs(za - zplane[p]);
         if(d < best){
            best = d;
            plane = p;
         }
      }
      double z_lo = za - 0.5*dzA;
      double z_hi = za + 0.5*dzA;
      if(nplane > 1){
         if(plane > 0) z_lo = 0.5*(zplane[plane-1] + zplane[plane]);
         if(plane < nplane-1) z_hi = 0.5*(zplane[plane] + zplane[plane+1]);
      }

      int i_lo = static_cast<int>(std::floor(z_lo / dzA));
      int i_hi = static_cast<int>(std::floor((z_hi - 1.0e-12) / dzA));
      if(i_lo < 0) i_lo = 0;
      if(i_hi < 0) i_hi = 0;
      if(i_lo >= nf) i_lo = nf-1;
      if(i_hi >= nf) i_hi = nf-1;
      if(i_hi < i_lo) i_hi = i_lo;

      sc1d_atom_ls[atom] = stack;
      sc1d_atom_fine_lo[atom] = i_lo;
      sc1d_atom_fine_hi[atom] = i_hi;
   }
}

void initialise_spincurrents_1d(){
   const int nsub = compute_fine_subdivisions();
   sc1d_nsub = nsub;
   sc1d_nf   = nsub * num_microcells_per_stack;

   collect_local_stacks();
   allocate_sc1d_arrays();
   initialise_magnetization_history();

   const int nf = sc1d_nf;
   const double dz = (micro_cell_thickness * 1.0e-10) / (double)nsub;
   const double dzA = micro_cell_thickness / (double)nsub;

   // Constants are the linear interpolation of the coarse-cell averages
   // between neighbouring coarse centres. Only the owner integrates that column.
   for(int stack=0; stack<num_stacks_y; ++stack){
      const int start_cell = stack_index_y[stack];

      double* Bc_fine = sc1d_Bc_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* Bd_fine = sc1d_Bd_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* D_fine = sc1d_D_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* lsf_fine = sc1d_lsf_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* lphi_fine = sc1d_lphi_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* Jsd_fine = sc1d_Jsd_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* chi_fine = sc1d_chi_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* sa_inf_fine = sc1d_sa_inf_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* sigma_fine = sc1d_sigma_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* seebeck_fine = sc1d_seebeck_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* alpha_edge_stack = sc1d_alpha_edge.data() + (std::size_t)stack * (std::size_t)std::max(0,nf-1);

      prolong_coarse_constants(start_cell, nsub, nf, Bc_fine, Bd_fine, D_fine, lsf_fine,
                               lphi_fine, Jsd_fine, chi_fine, sa_inf_fine, sigma_fine, seebeck_fine);
      assemble_alpha_edges_from_coarse(start_cell, nsub, nf, dz, D_fine, alpha_edge_stack);
      if(!rank_owns_column(stack)) continue;

      double* Sbase = sc1d_Sfine.data() + (std::size_t)stack * (std::size_t)nf * 3u;
      initialise_fine_spin_if_empty(start_cell, nsub, nf, dz, D_fine, lsf_fine,
                                    alpha_edge_stack, sa_inf_fine, Sbase);

      std::vector<double> mhat_init(3*(std::size_t)nf, 0.0);
      for(int i=0;i<nf;i++){
         const int k = i / std::max(1, nsub);
         const int cell = start_cell + k;
         mhat_init[3*(std::size_t)i+0] = sc1d_mhat_ref[3*cell+0];
         mhat_init[3*(std::size_t)i+1] = sc1d_mhat_ref[3*cell+1];
         mhat_init[3*(std::size_t)i+2] = sc1d_mhat_ref[3*cell+2];
      }
      double* ns_fine = sc1d_ns_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* ne_fine = sc1d_ne_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* V_fine = sc1d_V_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* Jc_edge_fine = sc1d_Jc_edge_fine.data() + (std::size_t)stack * (std::size_t)(nf + 1);
      stopwatch_t sw_charge;
      const double interp_before = sc1d_time_interp;
      sw_charge.start();
      solve_open_circuit_charge_steady(ns_fine, ne_fine, V_fine, Jc_edge_fine,
                                       nf, dz, 0.0, D_fine, sigma_fine, Bd_fine, seebeck_fine,
                                       Sbase, mhat_init.data());
      account_charge_time(sw_charge.elapsed_seconds(), interp_before);
   }

   build_atom_fine_map(nf, dzA);

   sc1d_atom_use_phonon.assign(num_local_atoms, false);
   for(int atom=0; atom<num_local_atoms; ++atom){
      const int mat = atoms::type_array[atom];
      if(mat >= 0 && mat < (int)mp::material.size()){
         sc1d_atom_use_phonon[atom] = mp::material[mat].couple_to_phonon_temperature;
      }
   }

}
//------------------------------------------------------------------------------
// Per-step 1D spin-accumulation helpers
//------------------------------------------------------------------------------

static void update_magnetization_rates(const int start_cell,
                                       const int ncz,
                                       const double dt,
                                       double* dm_dt_coarse){
   for(int k=0;k<ncz;k++){
      const int cell = start_cell + k;
      const double mx = m[3*cell+0];
      const double my = m[3*cell+1];
      const double mz = m[3*cell+2];
      const double mm = norm3(mx,my,mz);

      if(sc1d_m0_mag[cell] <= 0.0 && mm > 0.0) sc1d_m0_mag[cell] = mm;

      const double thresh = kFreezeFracOfM0 * std::max(sc1d_m0_mag[cell], 1e-12);

      if(!kFreezeAxisWhenMsmall || mm > thresh){
         if(mm > kEps){
            sc1d_mhat_ref[3*cell+0] = mx/mm;
            sc1d_mhat_ref[3*cell+1] = my/mm;
            sc1d_mhat_ref[3*cell+2] = mz/mm;
         }
      }

      if(sc1d_step_counter <= 1UL){
         dm_dt_coarse[k] = 0.0;
         sc1d_m_prev_mag[cell] = mm;
      } else {
         dm_dt_coarse[k] = (mm - sc1d_m_prev_mag[cell]) / dt;
         sc1d_m_prev_mag[cell] = mm;
      }
   }
}

static void fill_fine_time_dependent_fields(const int start_cell,
                                            const int ncz,
                                            const int nsub,
                                            const double* dm_dt_coarse,
                                            double* mhat,
                                            double* msrc,
                                            double* dmdt){
   for(int k=0;k<ncz;k++){
      const int cell1 = start_cell + k;

      const double mx1 = m[3*cell1+0], my1=m[3*cell1+1], mz1=m[3*cell1+2];

      const double dm1 = dm_dt_coarse[k];

      const double sx = sc1d_mhat_ref[3*cell1+0];
      const double sy = sc1d_mhat_ref[3*cell1+1];
      const double sz = sc1d_mhat_ref[3*cell1+2];

      for(int j=0;j<nsub;j++){
         const int i = k*nsub + j;

         double hx = mx1, hy = my1, hz = mz1;
         const double hm = norm3(hx,hy,hz);
         if(hm > kEps){ hx/=hm; hy/=hm; hz/=hm; }
         else { hx = sx; hy = sy; hz = sz; }

         mhat[3*i+0]=hx; mhat[3*i+1]=hy; mhat[3*i+2]=hz;
         msrc[3*i+0]=sx; msrc[3*i+1]=sy; msrc[3*i+2]=sz;
         dmdt[i] = dm1;
      }
   }
}

static void assemble_diffusion_tridiagonal(const int nf,
                                           const double dt_half,
                                           const double* alpha_edge,
                                           std::vector<double>& a_tr,
                                           std::vector<double>& b_tr,
                                           std::vector<double>& c_tr){
   if(nf == 1){
      a_tr[0]=0.0; b_tr[0]=1.0; c_tr[0]=0.0;
      return;
   }
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

static void apply_diffusion_half_step(double* S,
                                      const int nf,
                                      const double dz,
                                      const double dt_half,
                                      const double* Bc,
                                      const double* mhat,
                                      const double* Jc_edge_fine,
                                      const std::vector<double>& a_tr,
                                      const std::vector<double>& b_tr,
                                      const std::vector<double>& c_tr){
   const double jdl0x = -Bc[0]    * Jc_edge_fine[0]  * mhat[0];
   const double jdl0y = -Bc[0]    * Jc_edge_fine[0]  * mhat[1];
   const double jdl0z = -Bc[0]    * Jc_edge_fine[0]  * mhat[2];
   const double jdrNx = -Bc[nf-1] * Jc_edge_fine[nf] * mhat[3*(nf-1)+0];
   const double jdrNy = -Bc[nf-1] * Jc_edge_fine[nf] * mhat[3*(nf-1)+1];
   const double jdrNz = -Bc[nf-1] * Jc_edge_fine[nf] * mhat[3*(nf-1)+2];

   const double Jleft[3]  = { -jdl0x, -jdl0y, -jdl0z };
   const double Jright[3] = { -jdrNx, -jdrNy, -jdrNz };

   std::vector<double> rhs(nf, 0.0);
   for(int comp=0; comp<3; ++comp){
      for(int i=0;i<nf;i++){
         rhs[i] = S[3*i + comp];
      }
      rhs[0]     += dt_half * (Jleft[comp] / dz);
      rhs[nf-1]  -= dt_half * (Jright[comp] / dz);

      thomas_solve(a_tr,b_tr,c_tr,rhs);

      for(int i=0;i<nf;i++){
         S[3*i + comp] = rhs[i];
      }
   }
}

static void build_drift_edges(const int nf,
                              const double* Bc,
                              const double* mhat,
                              const double* Jc_edge_fine,
                              double* Jd_edge){
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
}

static void compute_local_rhs(const double* S,
                              std::vector<double>& kout,
                              const int nf,
                              const double dz,
                              const double* Bc,
                              const double* Bd,
                              const double* D,
                              const double* lsf,
                              const double* lphi,
                              const double* Jsd,
                              const double* chi,
                              const double* sa_inf,
                              const double* mhat,
                              const double* msrc,
                              const double* dmdt,
                              const double* Jd_edge){
   (void)msrc;
   for(int i=0;i<nf;i++){
      const double divx = -(Jd_edge[3*(i+1)+0] - Jd_edge[3*i+0]) / dz;
      const double divy = -(Jd_edge[3*(i+1)+1] - Jd_edge[3*i+1]) / dz;
      const double divz = -(Jd_edge[3*(i+1)+2] - Jd_edge[3*i+2]) / dz;

      const double denom = std::max(1e-12, 1.0 - Bc[i]*Bd[i]);
      const double BBp = 1.0 / std::sqrt(denom);
      double lambda0 = lsf[i];
      double D0 = D[i];
      if(sc1d_thermal_effects){
         const double zc = (i + 0.5) * dz;
         const double Te = get_local_electron_temperature(zc);
         const double Tp = get_local_phonon_temperature(zc);
         lambda0 = scale_lambda_sdl(lambda0, Te, Tp);
         D0 = scale_diffusion(D0, Te, Tp);
      }
      const double lambda_sf = lambda0 * BBp;
      const double Di = std::max(D0, 1e-30);
      const double inv_lambda_sf_sq = (lambda_sf > 0.0) ? (1.0 / (lambda_sf*lambda_sf)) : 0.0;
      const double inv_lambda_phi_sq = (lphi[i] > 0.0) ? (1.0 / (lphi[i]*lphi[i])) : 0.0;

      const double omega = (Jsd[i] > 0.0) ? (Jsd[i] / (2.0*kHbar)) : 0.0;

      const double sx = S[3*i+0];
      const double sy = S[3*i+1];
      const double sz = S[3*i+2];

      const double hx = mhat[3*i+0];
      const double hy = mhat[3*i+1];
      const double hz = mhat[3*i+2];

      const double sdot = sx*hx + sy*hy + sz*hz;
      const double spx = sx - sdot*hx;
      const double spy = sy - sdot*hy;
      const double spz = sz - sdot*hz;

      double cx,cy,cz;
      cross3(sx,sy,sz, hx,hy,hz, cx,cy,cz);
      const double prex = -omega * cx;
      const double prey = -omega * cy;
      const double prez = -omega * cz;

      // Attractor length is the material sa_inf. Direction follows m-hat.
      // Demagnetization dumps spin into S along m-hat; spin-flip returns the
      // excess to sa_inf. chi is demag-spin-coupling, (C/m^3) per muB.
      const double sa_eq_amp = sa_inf[i];
      const double sa_eq_x = sa_eq_amp * hx;
      const double sa_eq_y = sa_eq_amp * hy;
      const double sa_eq_z = sa_eq_amp * hz;
      const double demag = -chi[i] * dmdt[i];
      const double relx = Di * (-(sx - sa_eq_x) * inv_lambda_sf_sq);
      const double rely = Di * (-(sy - sa_eq_y) * inv_lambda_sf_sq);
      const double relz = Di * (-(sz - sa_eq_z) * inv_lambda_sf_sq);

      const double dephx = Di * (-spx * inv_lambda_phi_sq);
      const double dephy = Di * (-spy * inv_lambda_phi_sq);
      const double dephz = Di * (-spz * inv_lambda_phi_sq);

      kout[3*i+0] = divx + prex + relx + dephx + demag*hx;
      kout[3*i+1] = divy + prey + rely + dephy + demag*hy;
      kout[3*i+2] = divz + prez + relz + dephz + demag*hz;
   }
}

static void integrate_local_terms(double* S,
                                  double* k_prev,
                                  const int nf,
                                  const double dt,
                                  const bool first_step,
                                  const double dz,
                                  const double* Bc,
                                  const double* Bd,
                                  const double* D,
                                  const double* lsf,
                                  const double* lphi,
                                  const double* Jsd,
                                  const double* chi,
                                  const double* sa_inf,
                                  const double* mhat,
                                  const double* msrc,
                                  const double* dmdt,
                                  const double* Jd_edge){
   const int n3 = 3*nf;
   std::vector<double> k1(n3,0.0);

   if(first_step){
      std::vector<double> k2(n3,0.0);
      std::vector<double> Spred(n3,0.0);
      compute_local_rhs(S, k1, nf, dz,
                        Bc, Bd, D, lsf, lphi, Jsd, chi, sa_inf,
                        mhat, msrc, dmdt, Jd_edge);
      for(int i=0;i<n3;i++){
         Spred[i] = S[i] + dt * k1[i];
      }
      compute_local_rhs(Spred.data(), k2, nf, dz,
                        Bc, Bd, D, lsf, lphi, Jsd, chi, sa_inf,
                        mhat, msrc, dmdt, Jd_edge);
      for(int i=0;i<n3;i++){
         S[i] = S[i] + 0.5*dt*(k1[i] + k2[i]);
         k_prev[i] = k1[i];
      }
   } else {
      compute_local_rhs(S, k1, nf, dz,
                        Bc, Bd, D, lsf, lphi, Jsd, chi, sa_inf,
                        mhat, msrc, dmdt, Jd_edge);
      for(int i=0;i<n3;i++){
         S[i] = S[i] + dt * (1.5*k1[i] - 0.5*k_prev[i]);
         k_prev[i] = k1[i];
      }
   }
}

static void fine_cell_spin_current(const double* Sbase,
                                   const int i,
                                   const int nf,
                                   const double dz,
                                   const double Jc_loc,
                                   const double Bc_i,
                                   const double D_i,
                                   const double* mhat,
                                   double& jsx,
                                   double& jsy,
                                   double& jsz){
   const double jdr_x = -Bc_i * Jc_loc * mhat[3*i+0];
   const double jdr_y = -Bc_i * Jc_loc * mhat[3*i+1];
   const double jdr_z = -Bc_i * Jc_loc * mhat[3*i+2];

   double dSdz_x = 0.0;
   double dSdz_y = 0.0;
   double dSdz_z = 0.0;
   if(i==0){
      dSdz_x = (Sbase[3*(i+1)+0] - Sbase[3*i+0])/dz;
      dSdz_y = (Sbase[3*(i+1)+1] - Sbase[3*i+1])/dz;
      dSdz_z = (Sbase[3*(i+1)+2] - Sbase[3*i+2])/dz;
   } else if(i==nf-1){
      dSdz_x = (Sbase[3*i+0] - Sbase[3*(i-1)+0])/dz;
      dSdz_y = (Sbase[3*i+1] - Sbase[3*(i-1)+1])/dz;
      dSdz_z = (Sbase[3*i+2] - Sbase[3*(i-1)+2])/dz;
   } else {
      dSdz_x = (Sbase[3*(i+1)+0] - Sbase[3*(i-1)+0])/(2.0*dz);
      dSdz_y = (Sbase[3*(i+1)+1] - Sbase[3*(i-1)+1])/(2.0*dz);
      dSdz_z = (Sbase[3*(i+1)+2] - Sbase[3*(i-1)+2])/(2.0*dz);
   }

   if(!std::isfinite(dSdz_x)) dSdz_x = 0.0;
   if(!std::isfinite(dSdz_y)) dSdz_y = 0.0;
   if(!std::isfinite(dSdz_z)) dSdz_z = 0.0;

   const double max_grad = 1e20;
   const double dSdz_mag = std::sqrt(dSdz_x*dSdz_x + dSdz_y*dSdz_y + dSdz_z*dSdz_z);
   if(dSdz_mag > max_grad && dSdz_mag > kEps){
      const double scale = max_grad / dSdz_mag;
      dSdz_x *= scale;
      dSdz_y *= scale;
      dSdz_z *= scale;
   }

   double jdiff_x = -D_i * dSdz_x;
   double jdiff_y = -D_i * dSdz_y;
   double jdiff_z = -D_i * dSdz_z;
   if(!std::isfinite(jdiff_x)) jdiff_x = 0.0;
   if(!std::isfinite(jdiff_y)) jdiff_y = 0.0;
   if(!std::isfinite(jdiff_z)) jdiff_z = 0.0;

   jsx = jdr_x + jdiff_x;
   jsy = jdr_y + jdiff_y;
   jsz = jdr_z + jdiff_z;
}

static void coarsen_spin_and_map_torque(const int start_cell,
                                        const int ncz,
                                        const int nsub,
                                        const int nf,
                                        const double dz,
                                        const double* Sbase,
                                        const double* ns_fine,
                                        const double* ne_fine,
                                        const double* Jc_edge_fine,
                                        const double* Bc,
                                        const double* D,
                                        const double* mhat){
   // Two passes so the Js construction and the fine-to-coarse averages are
   // timed separately. Each sum is unchanged.
   stopwatch_t sw;
   sw.start();
   for(int k=0;k<ncz;k++){
      const int cell = start_cell + k;

      double ns_sum = 0.0;
      double ne_sum = 0.0;
      double jc_sum = 0.0;
      double ax=0.0, ay=0.0, az=0.0;

      for(int j=0;j<nsub;j++){
         const int i = k*nsub + j;
         ns_sum += ns_fine[i];
         ne_sum += ne_fine[i];
         jc_sum += 0.5*(Jc_edge_fine[i] + Jc_edge_fine[i+1]);
         ax += Sbase[3*i+0];
         ay += Sbase[3*i+1];
         az += Sbase[3*i+2];
      }
      const double inv = 1.0 / (double)nsub;
      ns_final[cell] = ns_sum * inv;
      ne_final[cell] = ne_sum * inv;
      jc_final[cell] = jc_sum * inv;
      ax *= inv;
      ay *= inv;
      az *= inv;

      sa_final[3*cell+0] = ax;
      sa_final[3*cell+1] = ay;
      sa_final[3*cell+2] = az;

      if(cell_natom[cell] > 0.5){
         if(sc1d_step_counter <= sc1d_relax_steps){
            ax = 0.0;
            ay = 0.0;
            az = 0.0;
         }

         spin_torque[3*cell+0] = kAtomcellVolume * sd_exchange[cell] * ax * kInvE * kInvMuB;
         spin_torque[3*cell+1] = kAtomcellVolume * sd_exchange[cell] * ay * kInvE * kInvMuB;
         spin_torque[3*cell+2] = kAtomcellVolume * sd_exchange[cell] * az * kInvE * kInvMuB;
      }
   }
   sc1d_time_interp += sw.elapsed_seconds();

   sw.start();
   for(int k=0;k<ncz;k++){
      const int cell = start_cell + k;
      double jsx_sum=0.0, jsy_sum=0.0, jsz_sum=0.0;

      for(int j=0;j<nsub;j++){
         const int i = k*nsub + j;
         const double Jc_loc = 0.5*(Jc_edge_fine[i] + Jc_edge_fine[i+1]);
         double jsx = 0.0;
         double jsy = 0.0;
         double jsz = 0.0;
         fine_cell_spin_current(Sbase, i, nf, dz, Jc_loc, Bc[i], D[i], mhat, jsx, jsy, jsz);
         jsx_sum += jsx;
         jsy_sum += jsy;
         jsz_sum += jsz;
      }
      const double inv = 1.0 / (double)nsub;
      js_final[3*cell+0] = jsx_sum * inv;
      js_final[3*cell+1] = jsy_sum * inv;
      js_final[3*cell+2] = jsz_sum * inv;
   }
   sc1d_time_spin_current += sw.elapsed_seconds();
}

// Every rank receives every column's cell fields in one sum. Owners contribute
// the updated slice; everyone else contributes zero.
static void broadcast_cell_fields(){
   const int nf = sc1d_nf;
   const int n3 = nf * 3;
   if(nf > 0 && n3 > 0){
      for(int s=0; s<num_stacks_y; ++s){
         if(rank_owns_column(s)) continue;
         const std::size_t off = (std::size_t)s * (std::size_t)n3;
         if(off + (std::size_t)n3 > sc1d_Sfine.size()) break;
         std::fill(sc1d_Sfine.begin() + static_cast<std::ptrdiff_t>(off),
                   sc1d_Sfine.begin() + static_cast<std::ptrdiff_t>(off + (std::size_t)n3),
                   0.0);
      }
   }
#ifdef MPICF
   if(vmpi::num_processors <= 1) return;
   const std::size_t npack = sc1d_Sfine.size() + sa_final.size() + js_final.size()
      + ns_final.size() + ne_final.size() + jc_final.size() + spin_torque.size();
   if(npack == 0) return;
   std::vector<double> packed(npack, 0.0);
   std::size_t off = 0;
   const std::vector<double>* parts[7] = {
      &sc1d_Sfine, &sa_final, &js_final, &ns_final, &ne_final, &jc_final, &spin_torque
   };
   for(int p=0; p<7; ++p){
      const std::vector<double>& src = *parts[p];
      std::copy(src.begin(), src.end(), packed.begin() + static_cast<std::ptrdiff_t>(off));
      off += src.size();
   }
   MPI_Allreduce(MPI_IN_PLACE, packed.data(), static_cast<int>(packed.size()),
                 MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   off = 0;
   std::vector<double>* dests[7] = {
      &sc1d_Sfine, &sa_final, &js_final, &ns_final, &ne_final, &jc_final, &spin_torque
   };
   for(int p=0; p<7; ++p){
      std::vector<double>& dst = *dests[p];
      std::copy(packed.begin() + static_cast<std::ptrdiff_t>(off),
                packed.begin() + static_cast<std::ptrdiff_t>(off + dst.size()),
                dst.begin());
      off += dst.size();
   }
#endif
}

//------------------------------------------------------------------------------
// Main update
//------------------------------------------------------------------------------

void calculate_spin_accumulation_1d(){
   // st::initialise runs before ltmp::initialise, so this check waits until
   // the first step, when a requested local-temperature pulse is initialised.
   if(sc1d_superdiffusive_enable && !ltmp::is_enabled()){
      terminaltextcolor(RED);
      std::cerr << "Error: spin-currents-1d-superdiffusive-enable requires the local temperature pulse module." << std::endl;
      zlog << zTs() << "Error: spin-currents-1d-superdiffusive-enable requires the local temperature pulse module." << std::endl;
      terminaltextcolor(WHITE);
      err::vexit();
   }

   stopwatch_t sw_total;
   sw_total.start();
   sc1d_step_timer_active = true;
   stopwatch_t sw;

   if(sc1d_nf <= 0 || sc1d_Sfine.empty()){
      initialise_spincurrents_1d();
   }

   std::fill(spin_torque.begin(), spin_torque.end(), 0.0);
   std::fill(sa_final.begin(),    sa_final.end(),    0.0);
   std::fill(ns_final.begin(),    ns_final.end(),    0.0);
   std::fill(ne_final.begin(),    ne_final.end(),    0.0);
   std::fill(jc_final.begin(),    jc_final.end(),    0.0);
   std::fill(js_final.begin(),    js_final.end(),    0.0);

   const double dt = mp::dt_SI;
   if(dt <= 0.0){
      close_sc1d_step_timer(sw_total);
      return;
   }
   const double dt_half = 0.5*dt;

   ++sc1d_step_counter;
   const double t_s = (double)sc1d_step_counter * dt;

   const int ncz  = num_microcells_per_stack;
   const int nsub = sc1d_nsub;
   const int nf   = sc1d_nf;
   const double dz = (micro_cell_thickness * 1.0e-10) / (double)nsub;

   std::vector<double> dmdt(nf,0.0);
   std::vector<double> mhat(3*nf,0.0);
   std::vector<double> msrc(3*nf,0.0);
   std::vector<double> a_tr(nf,0.0);
   std::vector<double> b_tr(nf,0.0);
   std::vector<double> c_tr(nf,0.0);
   std::vector<double> Jd_edge(3*(nf+1),0.0);
   std::vector<double> dm_dt_coarse(ncz,0.0);

   const bool first_step = (sc1d_step_counter == 1UL) ||
                           (sc1d_step_counter == sc1d_relax_steps + 1UL);

   for(int ls=0; ls<(int)sc1d_local_stacks.size(); ++ls){
      const int stack = sc1d_local_stacks[ls];
      const int start_cell = stack_index_y[stack];
      double* Sbase = sc1d_Sfine.data() + (size_t)stack * (size_t)nf * 3u;
      double* ns_fine = sc1d_ns_fine.data() + (size_t)stack * (size_t)nf;
      double* ne_fine = sc1d_ne_fine.data() + (size_t)stack * (size_t)nf;
      double* V_fine = sc1d_V_fine.data() + (size_t)stack * (size_t)nf;
      double* Jc_edge_fine = sc1d_Jc_edge_fine.data() + (size_t)stack * (size_t)(nf + 1);
      const double* seebeck_fine = sc1d_seebeck_fine.data() + (size_t)stack * (size_t)nf;

      const double* Bc = sc1d_Bc_fine.data() + (size_t)stack * (size_t)nf;
      const double* Bd = sc1d_Bd_fine.data() + (size_t)stack * (size_t)nf;
      const double* D = sc1d_D_fine.data() + (size_t)stack * (size_t)nf;
      const double* lsf = sc1d_lsf_fine.data() + (size_t)stack * (size_t)nf;
      const double* lphi = sc1d_lphi_fine.data() + (size_t)stack * (size_t)nf;
      const double* Jsd = sc1d_Jsd_fine.data() + (size_t)stack * (size_t)nf;
      const double* chi = sc1d_chi_fine.data() + (size_t)stack * (size_t)nf;
      const double* sa_inf = sc1d_sa_inf_fine.data() + (size_t)stack * (size_t)nf;
      const double* sigma0 = sc1d_sigma_fine.data() + (size_t)stack * (size_t)nf;
      const double* alpha_edge = sc1d_alpha_edge.data() + (size_t)stack * (size_t)std::max(0,nf-1);
      double* k_prev = sc1d_k_prev.data() + (size_t)stack * (size_t)nf * 3u;

      update_magnetization_rates(start_cell, ncz, dt, dm_dt_coarse.data());
      sw.start();
      fill_fine_time_dependent_fields(start_cell, ncz, nsub,
                                      dm_dt_coarse.data(),
                                      mhat.data(), msrc.data(), dmdt.data());
      sc1d_time_interp += sw.elapsed_seconds();

      const bool do_charge_update = (sc1d_charge_stride <= 1) || ((sc1d_step_counter % (unsigned long)sc1d_charge_stride) == 0UL);
      if(do_charge_update){
         const double interp_before = sc1d_time_interp;
         sw.start();
         double* ne_seebeck = sc1d_ne_seebeck_fine.data() + (size_t)stack * (size_t)nf;
         update_charge_transients_1d_stack(ns_fine, ne_fine, V_fine, Jc_edge_fine, ne_seebeck,
                                           nf, dz, dt, t_s,
                                           D, sigma0, Bd, seebeck_fine, Sbase, mhat.data());
         account_charge_time(sw.elapsed_seconds(), interp_before);
      }

      sw.start();
      assemble_diffusion_tridiagonal(nf, dt_half, alpha_edge, a_tr, b_tr, c_tr);
      apply_diffusion_half_step(Sbase, nf, dz, dt_half, Bc, mhat.data(), Jc_edge_fine, a_tr, b_tr, c_tr);
      sc1d_time_spin_acc += sw.elapsed_seconds();
      sw.start();
      build_drift_edges(nf, Bc, mhat.data(), Jc_edge_fine, Jd_edge.data());
      sc1d_time_spin_current += sw.elapsed_seconds();
      sw.start();
      integrate_local_terms(Sbase, k_prev, nf, dt, first_step, dz,
                            Bc, Bd, D, lsf, lphi, Jsd, chi, sa_inf,
                            mhat.data(), msrc.data(), dmdt.data(), Jd_edge.data());
      apply_diffusion_half_step(Sbase, nf, dz, dt_half, Bc, mhat.data(), Jc_edge_fine, a_tr, b_tr, c_tr);
      sc1d_time_spin_acc += sw.elapsed_seconds();
      coarsen_spin_and_map_torque(start_cell, ncz, nsub, nf, dz, Sbase,
                                  ns_fine, ne_fine, Jc_edge_fine,
                                  Bc, D, mhat.data());
   }

   // Owners hold the updated slices. Unowned spin slices still contain the
   // previous broadcast, so clear them before the sum. Coarse cell fields
   // were zeroed at the start of the step and written only by their owner.
   sw.start();
   broadcast_cell_fields();
   sc1d_time_bcast += sw.elapsed_seconds();

   sw.start();
   output_microcell_data();
   sc1d_time_io += sw.elapsed_seconds();

   close_sc1d_step_timer(sw_total);
#if ST_SC1D_TIMINGS
   report_sc1d_timing();
#endif
}

} // namespace internal
} // namespace st
