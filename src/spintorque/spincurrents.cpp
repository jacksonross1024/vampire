//------------------------------------------------------------------------------
// spincurrents.cpp
//
// Minimal 1D (along st::internal::stz) transient spin-accumulation
// solver, with:
//  - implicit (backward-Euler) anisotropic diffusion on a fine z-grid
//  - explicit Heun step for local precession/relaxation + drift-divergence.
//    Equilibrium amplitude is sa_inf (direction follows m-hat). This is the
//    full equilibrium spin accumulation, the spin-split Fermi density, not the
//    micromagnetic deviation from equilibrium. Demagnetization
//    dumps -chi*d(|m|/msat)/dt*m-hat into S; spin-flip returns the excess to sa_inf.
//  - Neumann outer boundaries (no Dirichlet spin sinks), internal Robin-like
//    interface coupling via st::internal::r_int_edge
//
// Notes:
//  - This implementation is designed for stack-by-stack (x,y) independence:
//    no lateral gradients/flows are taken.
//  - Laser-driven charge transient on the fine grid for Eqs. (1)-(3),
//    producing a 1D charge current Jc(z) used in the spin drift term (Eq. 4).
//  - Material constants are occupied on the fine grid and prolonged from
//    those cells. Only the magnetisation is prolonged again each step.
//------------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "internal.hpp"
#include "material.hpp"
#include "vio.hpp"

#include "program.hpp"
#include "sim.hpp"
#include "vmpi.hpp"
#include "ltmp.hpp"
#include "errors.hpp"
#include "atoms.hpp"
#include "stopwatch.hpp"
#include "vio.hpp"

namespace st {
namespace internal {

// True while calculate_spin_accumulation_1d's total stopwatch is running.
// Init-time charge work is then already inside that total.
static bool sc1d_step_timer_active = false;

static void account_charge_time(const double elapsed, const double interp_before){
   if(!sc1d_step_timer_active) return;
   const double interp_dt = sc1d_time_interp - interp_before;
   sc1d_time_charge += elapsed - interp_dt;
}

static void close_sc1d_step_timer(stopwatch_t& sw_total){
   const double step_dt = sw_total.elapsed_seconds();
   sc1d_time_total += step_dt;
   sc1d_time_other = sc1d_time_total
      - (sc1d_time_io + sc1d_time_charge + sc1d_time_spin_current
         + sc1d_time_diffusion + sc1d_time_eom + sc1d_time_interp
         + sc1d_time_bcast + sc1d_time_coarsen);
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
// Thermal-window density n = g(Ef) k_B T from the Co DoS at 300 K.
// With material-file electrical-conductivity set, this is only the Einstein fallback.
static constexpr double kMetalN = 2.856e27;           // m^-3, Einstein n for σ fallback

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

// 0 K is a valid temperature. A non-finite value is left as it is and the
// initialisation smoke test fails the run.
static double finite_temperature(const double T){
   return T;
}

// Chosen once in configure_temperature_source, after ltmp::initialise.
// 0: reference temperature. 1: one global TTM temperature. 2: ltmp profile.
static int sc1d_temperature_source = 0;

// ltmp stores one temperature per vertical cell. It has no fine-grid gradient.
// delta_temperature_array is the heat-flow divergence, not dT/dz.
// Sample the temperature once per spin step and keep dT/dz on the fine faces.
// Array lengths and the temperature source are fixed for the run.
static void refresh_fine_temperature_profile(){
   const int nf = sc1d_nf;
   const double dz = (micro_cell_thickness * 1.0e-10) / (double)sc1d_nsub;
   if(sc1d_temperature_source == 2){
      const double z_scale = dz * 1.0e10;
      for(int i = 0; i < nf; ++i){
         const double zA = ((double)i + 0.5) * z_scale;
         sc1d_Te_fine[(std::size_t)i] = finite_temperature(ltmp::get_electron_temperature_at_z(zA));
         sc1d_Tp_fine[(std::size_t)i] = finite_temperature(ltmp::get_phonon_temperature_at_z(zA));
      }
   } else if(sc1d_temperature_source == 1){
      const double Te = finite_temperature(sim::TTTe);
      const double Tp = finite_temperature(sim::TTTp);
      std::fill(sc1d_Te_fine.begin(), sc1d_Te_fine.end(), Te);
      std::fill(sc1d_Tp_fine.begin(), sc1d_Tp_fine.end(), Tp);
   } else {
      const double Tref = sc1d_reference_temperature;
      std::fill(sc1d_Te_fine.begin(), sc1d_Te_fine.end(), Tref);
      std::fill(sc1d_Tp_fine.begin(), sc1d_Tp_fine.end(), Tref);
   }
   sc1d_dTe_dz_edge[0] = 0.0;
   sc1d_dTe_dz_edge[(std::size_t)nf] = 0.0;
   sc1d_dTp_dz_edge[0] = 0.0;
   sc1d_dTp_dz_edge[(std::size_t)nf] = 0.0;
   const double idz = 1.0 / dz;
   for(int e = 1; e < nf; ++e){
      sc1d_dTe_dz_edge[(std::size_t)e] = (sc1d_Te_fine[(std::size_t)e] - sc1d_Te_fine[(std::size_t)(e - 1)]) * idz;
      sc1d_dTp_dz_edge[(std::size_t)e] = (sc1d_Tp_fine[(std::size_t)e] - sc1d_Tp_fine[(std::size_t)(e - 1)]) * idz;
   }
}

// Bounds and the temperature source. Called once from initialise_spincurrents_1d,
// which runs after ltmp::initialise. The step path does not repeat this.
static void configure_temperature_source(){
   if(sc1d_superdiffusive_enable && !ltmp::is_enabled()){
      terminaltextcolor(RED);
      std::cerr << "Error: spin-currents-1d-superdiffusive-enable requires the local temperature pulse module." << std::endl;
      zlog << zTs() << "Error: spin-currents-1d-superdiffusive-enable requires the local temperature pulse module." << std::endl;
      terminaltextcolor(WHITE);
      err::vexit();
   }
   // program 6 is laser-pulse. program 19 is laser-electrical-pulse.
   if(ltmp::is_enabled()) sc1d_temperature_source = 2;
   else if(program::program == 6 || program::program == 19) sc1d_temperature_source = 1;
   else sc1d_temperature_source = 0;

   const int nf = sc1d_nf;
   const bool sized = nf > 0
      && (int)sc1d_Te_fine.size() == nf
      && (int)sc1d_Tp_fine.size() == nf
      && (int)sc1d_dTe_dz_edge.size() == nf + 1
      && (int)sc1d_dTp_dz_edge.size() == nf + 1
      && micro_cell_thickness > 0.0
      && sc1d_nsub >= 1;
   if(!sized){
      terminaltextcolor(RED);
      std::cerr << "Error: spin-current temperature arrays are not sized for the fine grid." << std::endl;
      zlog << zTs() << "Error: spin-current temperature arrays are not sized for the fine grid." << std::endl;
      terminaltextcolor(WHITE);
      err::vexit();
   }
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

   // tau_pw is sim:laser-pulse-time. Decay starts at the ltmp envelope peak.
   // Before the peak the exponential factor is 1.
   double time_fac = 1.0;
   if(sc1d_drive_active()){
      const double t_from_peak = drive_time_s() - ltmp::laser_peak_time();
      if(t_from_peak > 0.0){
         const double tau_pw = (sim::pump_time > 0.0) ? sim::pump_time : sc1d_tau_e;
         time_fac = std::exp(-t_from_peak / std::max(tau_pw, 1e-30));
      }
   }
   const double fac = temp_fac * time_fac;
   return std::max(fac, 0.1);
}

static double signed_hot_electron_velocity(const double De,
                                           const double Te,
                                           const double dTe_dz,
                                           const bool have_Te){
   // Constant speed. The input number is a speed; the direction is toward
   // the substrate, away from the cap. Same direction as the decay and the
   // temperature-gradient drift. No De, no post-peak factor.
   if(sc1d_superdiffusive_type == 1){
      return -std::fabs(sc1d_superdiffusive_speed);
   }
   // Lifetime decay. |ve| is 1e6 m/s at the start of the drive and 1e5 m/s
   // one lifetime later: v(t) = -1e6 * 10^(-t/tau). The sign is toward the
   // substrate, the same direction as the temperature-gradient drift.
   // De and the mixing sheet do not enter.
   if(sc1d_superdiffusive_type == 2){
      const double tau = std::max(sc1d_superdiffusive_lifetime, 1e-30);
      const double t = drive_time_s();
      return -1.0e6 * std::pow(10.0, -t / tau);
   }

   double ve = 0.0;
   // Gradient, with a local Te:
   //   ve = -De * sqrt(Te/Tref) * (dTe/dz)/Te * exp(-(t-t_peak)/tau_pw)
   // The sqrt and the exponential are superdiffusive_speed_factor. That
   // product is floored at 0.1. De is the material-file value.
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

static constexpr double kEps = 1e-30;

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
   static std::vector<double> cp;
   static std::vector<double> dp;
   if((int)cp.size() != n){
      cp.assign((std::size_t)n, 0.0);
      dp.assign((std::size_t)n, 0.0);
   }

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


// gamma = beta_c * beta_d. The smoke test rejects a polarisation outside
// [-1, 1] and a product that is not below 1. sigma <= 0 is not a gamma:
// 0 is the insulator and a negative stored value is the Einstein sentinel.
double spin_gamma(const double Bc, const double Bd){
   return Bc * Bd;
}

static double edge_harmonic_mean(const double left, const double right){
   if(left > 0.0 && right > 0.0){
      return 2.0*left*right / (left + right);
   }
   return 0.0;
}

// Face thermopower. Each half-cell sees half of a linear temperature drop,
// so the face value is the arithmetic mean. A harmonic mean would dip wherever
// the two materials differ.
static double edge_seebeck_coefficient(const double SL, const double SR){
   return 0.5 * (SL + SR);
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

// 3x3 block Thomas for the anisotropic spin-diffusion half-step.
// Same elimination as block_thomas. The inverse is the cofactor formula,
// with the same pivot guard as block2_inv.
struct Block3 {
   double a00;
   double a01;
   double a02;
   double a10;
   double a11;
   double a12;
   double a20;
   double a21;
   double a22;
};

struct Vec3 {
   double x;
   double y;
   double z;
};

static Block3 block3_zero(){
   Block3 A;
   A.a00 = 0.0; A.a01 = 0.0; A.a02 = 0.0;
   A.a10 = 0.0; A.a11 = 0.0; A.a12 = 0.0;
   A.a20 = 0.0; A.a21 = 0.0; A.a22 = 0.0;
   return A;
}

static Block3 block3_identity(){
   Block3 A = block3_zero();
   A.a00 = 1.0;
   A.a11 = 1.0;
   A.a22 = 1.0;
   return A;
}

static Block3 block3_add(const Block3& A, const Block3& B){
   Block3 C;
   C.a00 = A.a00 + B.a00; C.a01 = A.a01 + B.a01; C.a02 = A.a02 + B.a02;
   C.a10 = A.a10 + B.a10; C.a11 = A.a11 + B.a11; C.a12 = A.a12 + B.a12;
   C.a20 = A.a20 + B.a20; C.a21 = A.a21 + B.a21; C.a22 = A.a22 + B.a22;
   return C;
}

static Block3 block3_sub(const Block3& A, const Block3& B){
   Block3 C;
   C.a00 = A.a00 - B.a00; C.a01 = A.a01 - B.a01; C.a02 = A.a02 - B.a02;
   C.a10 = A.a10 - B.a10; C.a11 = A.a11 - B.a11; C.a12 = A.a12 - B.a12;
   C.a20 = A.a20 - B.a20; C.a21 = A.a21 - B.a21; C.a22 = A.a22 - B.a22;
   return C;
}

static Block3 block3_scale(const Block3& A, const double s){
   Block3 C;
   C.a00 = A.a00*s; C.a01 = A.a01*s; C.a02 = A.a02*s;
   C.a10 = A.a10*s; C.a11 = A.a11*s; C.a12 = A.a12*s;
   C.a20 = A.a20*s; C.a21 = A.a21*s; C.a22 = A.a22*s;
   return C;
}

static Block3 block3_mul(const Block3& A, const Block3& B){
   Block3 C;
   C.a00 = A.a00*B.a00 + A.a01*B.a10 + A.a02*B.a20;
   C.a01 = A.a00*B.a01 + A.a01*B.a11 + A.a02*B.a21;
   C.a02 = A.a00*B.a02 + A.a01*B.a12 + A.a02*B.a22;
   C.a10 = A.a10*B.a00 + A.a11*B.a10 + A.a12*B.a20;
   C.a11 = A.a10*B.a01 + A.a11*B.a11 + A.a12*B.a21;
   C.a12 = A.a10*B.a02 + A.a11*B.a12 + A.a12*B.a22;
   C.a20 = A.a20*B.a00 + A.a21*B.a10 + A.a22*B.a20;
   C.a21 = A.a20*B.a01 + A.a21*B.a11 + A.a22*B.a21;
   C.a22 = A.a20*B.a02 + A.a21*B.a12 + A.a22*B.a22;
   return C;
}

static Vec3 block3_mul_vec(const Block3& A, const Vec3& v){
   Vec3 r;
   r.x = A.a00*v.x + A.a01*v.y + A.a02*v.z;
   r.y = A.a10*v.x + A.a11*v.y + A.a12*v.z;
   r.z = A.a20*v.x + A.a21*v.y + A.a22*v.z;
   return r;
}

static Vec3 vec3_sub(const Vec3& a, const Vec3& b){
   Vec3 r;
   r.x = a.x - b.x;
   r.y = a.y - b.y;
   r.z = a.z - b.z;
   return r;
}

static Block3 block3_inv(const Block3& A){
   const double c00 = A.a11*A.a22 - A.a12*A.a21;
   const double c01 = A.a12*A.a20 - A.a10*A.a22;
   const double c02 = A.a10*A.a21 - A.a11*A.a20;
   const double c10 = A.a02*A.a21 - A.a01*A.a22;
   const double c11 = A.a00*A.a22 - A.a02*A.a20;
   const double c12 = A.a01*A.a20 - A.a00*A.a21;
   const double c20 = A.a01*A.a12 - A.a02*A.a11;
   const double c21 = A.a02*A.a10 - A.a00*A.a12;
   const double c22 = A.a00*A.a11 - A.a01*A.a10;
   double det = A.a00*c00 + A.a01*c01 + A.a02*c02;
   if(std::fabs(det) < kEps) det = (det >= 0.0 ? kEps : -kEps);
   const double idet = 1.0 / det;
   Block3 I;
   I.a00 = c00 * idet; I.a01 = c10 * idet; I.a02 = c20 * idet;
   I.a10 = c01 * idet; I.a11 = c11 * idet; I.a12 = c21 * idet;
   I.a20 = c02 * idet; I.a21 = c12 * idet; I.a22 = c22 * idet;
   return I;
}

static void block3_thomas(const std::vector<Block3>& L,
                          const std::vector<Block3>& D,
                          const std::vector<Block3>& U,
                          std::vector<Vec3>& rhs){
   const int n = static_cast<int>(rhs.size());
   if(n <= 0) return;
   static std::vector<Block3> Cp;
   static std::vector<Vec3> dp;
   if((int)Cp.size() != n){
      Cp.assign((std::size_t)n, block3_zero());
      dp.resize((std::size_t)n);
   }

   Block3 Dinv = block3_inv(D[0]);
   Cp[0] = block3_mul(Dinv, U[0]);
   dp[0] = block3_mul_vec(Dinv, rhs[0]);

   for(int i=1;i<n;i++){
      const Block3 LCp = block3_mul(L[i], Cp[i-1]);
      const Block3 Di = block3_sub(D[i], LCp);
      Dinv = block3_inv(Di);
      if(i < n-1){
         Cp[i] = block3_mul(Dinv, U[i]);
      } else {
         Cp[i] = block3_zero();
      }
      const Vec3 Ldp = block3_mul_vec(L[i], dp[i-1]);
      const Vec3 rhs_i = vec3_sub(rhs[i], Ldp);
      dp[i] = block3_mul_vec(Dinv, rhs_i);
   }

   rhs[n-1] = dp[n-1];
   for(int i=n-2;i>=0;i--){
      const Vec3 Cpx = block3_mul_vec(Cp[i], rhs[i+1]);
      rhs[i].x = dp[i].x - Cpx.x;
      rhs[i].y = dp[i].y - Cpx.y;
      rhs[i].z = dp[i].z - Cpx.z;
   }
}

// D_charge is the sheet-scaled diffusivity, D_bulk * (σ↑+σ↓) / σ_bulk.
// D_hot is the material diffusivity. The mixing sheet does not enter ve.
static void fill_fine_charge_transport(const int nf,
                                       const double dz,
                                       const double t_s,
                                       const double* D_charge,
                                       const double* D_hot,
                                       const double* sigma0_cell,
                                       double* Te_cell,
                                       double* De_cell,
                                       double* De_hot_cell,
                                       double* ve_cell,
                                       double* sigma_cell){
   (void)t_s;
   (void)dz;
   const bool have_Te = (sc1d_temperature_source == 2);
   for(int i=0;i<nf;i++){
      const double Te = sc1d_Te_fine[(std::size_t)i];
      const double Tp = sc1d_Tp_fine[(std::size_t)i];
      Te_cell[i] = Te;
      De_cell[i] = scale_diffusion(std::max(0.0, D_charge[i]), Te, Tp);
      De_hot_cell[i] = scale_diffusion(std::max(0.0, D_hot[i]), Te, Tp);
      // sigma0_cell is σ↑+σ↓ after the sheet. Einstein values on occupied
      // cells were resolved before prolongation, then the sheet rescaled them.
      sigma_cell[i] = apply_thermal_effects(sigma0_cell[i], Te, Tp);
   }

   for(int i=0;i<nf;i++){
      double dTe_dz = 0.0;
      if(i == 0) dTe_dz = sc1d_dTe_dz_edge[1];
      else if(i == nf - 1) dTe_dz = sc1d_dTe_dz_edge[(std::size_t)(nf - 1)];
      else dTe_dz = 0.5 * (sc1d_dTe_dz_edge[(std::size_t)i] + sc1d_dTe_dz_edge[(std::size_t)(i + 1)]);
      ve_cell[i] = signed_hot_electron_velocity(De_hot_cell[i], Te_cell[i], dTe_dz, have_Te);
   }
}

static void fill_fine_charge_faces(const int nf,
                                   const double* De_cell,
                                   const double* De_hot_cell,
                                   const double* ve_cell,
                                   const double* sigma_cell,
                                   double* De_edge,
                                   double* De_hot_edge,
                                   double* ve_edge,
                                   double* sigma_edge,
                                   double* eps_edge){
   for(int e=0;e<=nf;e++){
      De_edge[e] = 0.0;
      De_hot_edge[e] = 0.0;
      ve_edge[e] = 0.0;
      sigma_edge[e] = 0.0;
      eps_edge[e] = 0.0;
   }
   for(int e=1;e<nf;e++){
      De_edge[e] = edge_harmonic_mean(De_cell[e-1], De_cell[e]);
      De_hot_edge[e] = edge_harmonic_mean(De_hot_cell[e-1], De_hot_cell[e]);
      sigma_edge[e] = edge_harmonic_mean(sigma_cell[e-1], sigma_cell[e]);
      // The harmonic mean is 0 when either cell has no diffusivity. That face
      // cannot carry hot electrons, so it does not inherit the neighbour speed.
      if(De_hot_edge[e] > 0.0) ve_edge[e] = 0.5*(ve_cell[e-1] + ve_cell[e]);
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
      Tp_avg += sc1d_Tp_fine[(std::size_t)i];
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
   // Spin-diffusion charge only. Hot electrons are a separate piece, J_hot = -Jns.
   (void)Jns_edge;
   Jknown[0] = 0.0;
   Jknown[nf] = 0.0;
   for(int e=1;e<nf;e++){
      const double Bd_f = 0.5*(Bd_cell[e-1] + Bd_cell[e]);
      const double sL = S[3*(e-1)+0]*mhat[3*(e-1)+0] + S[3*(e-1)+1]*mhat[3*(e-1)+1] + S[3*(e-1)+2]*mhat[3*(e-1)+2];
      const double sR = S[3*e+0]*mhat[3*e+0] + S[3*e+1]*mhat[3*e+1] + S[3*e+2]*mhat[3*e+2];
      // S is already C/m^3, so De * d(S·m)/dz is A/m^2. No e/mu_B.
      const double Jspin = (De_edge[e] * Bd_f) * (sR - sL) / dz;
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

// One column of the steady charge problem. J0 is the series current.
// Open circuit passes 0. A bias run calls this twice, at 0 and at 1.
// Conducting faces keep sigma*E + J_cd + J_imp = J0. A barrier face is either
// the jump V_L - V_R = J0/G or a discrete Laplace row that makes V linear.
struct ChargeBarrier {
   int vL;
   int vR;
   int cell_lo;
   int cell_hi;
   int face_lo;
   int face_hi;
   int jump_face;
   int n_cells;
   double Gup;
   double Gdn;
   double lup;
   double ldn;
   double ge;
   double gp;
   double li;
   double pe;
   double pp;
   double G;
   double P;
   double sfrac;
   double mx;
   double my;
   double mz;
   bool thick;
};

static void abort_charge(const char* msg){
   terminaltextcolor(RED);
   std::cerr << "Error - " << msg << std::endl;
   terminaltextcolor(WHITE);
   zlog << zTs() << "Error - " << msg << std::endl;
   err::vexit();
}

// The electrical endpoint can be an empty fine cell. Walk out to the nearest
// cell that actually carries a material, which is the electrode.
static int electrode_material(const int* mat, const int start, const int dir, const int nf){
   if(mat == nullptr || dir == 0) return -1;
   for(int i = start; i >= 0 && i < nf; i += dir){
      if(mat[i] >= 0) return mat[i];
   }
   return -1;
}

static bool lookup_tunnel(const int a, const int b, double& gup, double& gdn, double& lup, double& ldn){
   gup = 0.0;
   gdn = 0.0;
   lup = 0.0;
   ldn = 0.0;
   if(a < 0 || b < 0) return false;
   // Materials declared after the tunnel line grow mp. Rebuild the pair
   // index to that size or the conductance is read from the old stride.
   ensure_tunnel_matrix_size(mp.size());
   const std::size_t nmat = mp.size();
   if(nmat == 0) return false;
   const std::size_t idx = static_cast<std::size_t>(a) * nmat + static_cast<std::size_t>(b);
   if(idx >= tunnel_pair_set.size() || !tunnel_pair_set[idx]) return false;
   gup = tunnel_g_up_pair[idx];
   gdn = tunnel_g_down_pair[idx];
   if(idx < tunnel_l_up_pair.size()) lup = tunnel_l_up_pair[idx];
   if(idx < tunnel_l_down_pair.size()) ldn = tunnel_l_down_pair[idx];
   return true;
}

// Magnetisation the oxide reads from one electrode. The jump itself stays on
// the insulating segment. Step reads the contact cell, which already holds the
// prolonged magnetisation when that cell is empty. Linear and smoothstep
// average from that cell back into the electrode, over the prolongation decay.
static void electrode_magnetisation(const double* mhat, const int contact, const int away,
                                    const int nf, const double dz,
                                    double& mx, double& my, double& mz){
   mx = mhat[3 * contact];
   my = mhat[3 * contact + 1];
   mz = mhat[3 * contact + 2];
   if(sc1d_prolongation == sc1d_prolong_step) return;
   const double decay = sc1d_prolongation_decay * 1.0e-10;
   if(!(decay > 0.0) || !(dz > 0.0) || away == 0) return;
   double sx = 0.0;
   double sy = 0.0;
   double sz = 0.0;
   double wsum = 0.0;
   for(int i = contact; i >= 0 && i < nf; i += away){
      const double dist = std::fabs(static_cast<double>(i - contact)) * dz;
      const double x = dist / decay;
      if(!(x < 1.0)) break;
      double w = 1.0 - x;
      if(sc1d_prolongation == sc1d_prolong_smooth){
         const double t = 1.0 - x;
         w = t * t * (3.0 - 2.0 * t);
      }
      sx += w * mhat[3 * i];
      sy += w * mhat[3 * i + 1];
      sz += w * mhat[3 * i + 2];
      wsum += w;
   }
   if(!(wsum > 0.0)) return;
   mx = sx / wsum;
   my = sy / wsum;
   mz = sz / wsum;
}

static double unit_or_zero(double& x, double& y, double& z){
   const double n2 = x * x + y * y + z * z;
   // A vanishing electrode has no direction. A non-finite norm is left in place.
   if(!(n2 > 1.0e-24)){
      if(n2 >= 0.0){
         x = 0.0;
         y = 0.0;
         z = 0.0;
         return 0.0;
      }
      return std::sqrt(n2);
   }
   const double n = std::sqrt(n2);
   x /= n;
   y /= n;
   z /= n;
   return n;
}

static void load_inelastic(const int a, const int b, ChargeBarrier& seg){
   seg.ge = 0.0;
   seg.gp = 0.0;
   seg.li = 0.0;
   seg.pe = 0.0;
   seg.pp = 0.0;
   if(a < 0 || b < 0) return;
   ensure_tunnel_matrix_size(mp.size());
   const std::size_t nmat = mp.size();
   if(nmat == 0) return;
   const std::size_t idx = static_cast<std::size_t>(a) * nmat + static_cast<std::size_t>(b);
   if(idx < tunnel_gi_e_pair.size()) seg.ge = tunnel_gi_e_pair[idx];
   if(idx < tunnel_gi_p_pair.size()) seg.gp = tunnel_gi_p_pair[idx];
   if(idx < tunnel_li_pair.size()) seg.li = tunnel_li_pair[idx];
   if(idx < tunnel_pe_pair.size()) seg.pe = tunnel_pe_pair[idx];
   if(idx < tunnel_pp_pair.size()) seg.pp = tunnel_pp_pair[idx];
}

static double thermal_power(const double T, const double power){
   if(!(power > 0.0)) return 1.0;
   if(!(T > 0.0)) return 0.0;
   return std::pow(T / 300.0, power);
}

// Mean electron and phonon temperature on the oxide. The charge solver's
// reference temperature is the fallback, so a run without a temperature
// profile still evaluates the room-temperature inelastic channel.
static void barrier_temperatures(const int stack, const ChargeBarrier& seg, const int nf,
                                 double& Te, double& Tp){
   (void)stack;
   const double fallback = sc1d_reference_temperature;
   Te = fallback;
   Tp = fallback;
   if(nf <= 0) return;
   double se = 0.0;
   double sp = 0.0;
   int n = 0;
   auto add = [&](const int i){
      if(i < 0 || i >= nf) return;
      se += sc1d_Te_fine[static_cast<std::size_t>(i)];
      sp += sc1d_Tp_fine[static_cast<std::size_t>(i)];
      ++n;
   };
   if(seg.thick && seg.cell_lo >= 0){
      for(int i = seg.cell_lo; i <= seg.cell_hi; ++i) add(i);
   } else {
      add(seg.vL);
      add(seg.vR);
   }
   if(n > 0){
      Te = se / static_cast<double>(n);
      Tp = sp / static_cast<double>(n);
   }
}

static double channel_conductance(const double g0, const double lambda, const double thickness){
   if(!(g0 > 0.0)) return 0.0;
   if(!(lambda > 0.0) || !(thickness > 0.0)) return g0;
   return g0 * std::exp(-thickness / lambda);
}

// Two spin channels. Each decays with its own length across the oxide thickness.
// Parallel is g_up + g_down. Antiparallel is the series combination of the two
// channels. The angle weight is the magnetisation sampled from each electrode.
static void set_barrier_angle(ChargeBarrier& seg, const double* mhat, const int* mat,
                              const int stack, const int nf, const double dz){
   seg.G = 0.0;
   seg.P = 0.0;
   seg.sfrac = 0.0;
   seg.mx = 0.0;
   seg.my = 0.0;
   seg.mz = 1.0;
   const double thickness = (seg.thick && seg.n_cells > 0) ? (static_cast<double>(seg.n_cells) * dz) : 0.0;
   const double gup = channel_conductance(seg.Gup, seg.lup, thickness);
   const double gdn = channel_conductance(seg.Gdn, seg.ldn, thickness);
   const double Gp = gup + gdn;
   const double Gap = (Gp > 0.0) ? (2.0 * gup * gdn / Gp) : 0.0;
   seg.P = (Gp > 0.0) ? ((gup - gdn) / Gp) : 0.0;
   double Te = 0.0;
   double Tp = 0.0;
   barrier_temperatures(stack, seg, nf, Te, Tp);
   const double g_inel =
      channel_conductance(seg.ge, seg.li, thickness) * thermal_power(Te, seg.pe)
      + channel_conductance(seg.gp, seg.li, thickness) * thermal_power(Tp, seg.pp);
   if(!(Gp > 0.0) && !(g_inel > 0.0)) return;

   double lx = 0.0;
   double ly = 0.0;
   double lz = 1.0;
   double rx = 0.0;
   double ry = 0.0;
   double rz = 1.0;
   electrode_magnetisation(mhat, seg.vL, -1, nf, dz, lx, ly, lz);
   electrode_magnetisation(mhat, seg.vR, +1, nf, dz, rx, ry, rz);
   const double nL = unit_or_zero(lx, ly, lz);
   const double nR = unit_or_zero(rx, ry, rz);
   double cos_th = 1.0;
   if(nL > 0.0 && nR > 0.0 && seg.vL != seg.vR){
      cos_th = lx * rx + ly * ry + lz * rz;
      if(cos_th > 1.0) cos_th = 1.0;
      if(cos_th < -1.0) cos_th = -1.0;
   }
   const double wp = 0.5 * (1.0 + cos_th);
   const double wa = 0.5 * (1.0 - cos_th);
   const double g_el = Gp * wp + Gap * wa;
   seg.G = g_el + g_inel;
   // The inelastic channel is unpolarised, so only the elastic imbalance
   // contributes to the spin flux. Heating lowers the polarisation by
   // increasing the charge conductance in the denominator.
   const double Gs = (gup - gdn) * wp;
   seg.sfrac = (seg.G > 0.0) ? (Gs / seg.G) : 0.0;

   double mx = lx;
   double my = ly;
   double mz = lz;
   if(mat != nullptr && mat[seg.vR] == reference_layer && nR > 0.0){
      mx = rx;
      my = ry;
      mz = rz;
   }
   if(mat != nullptr && mat[seg.vL] == reference_layer && nL > 0.0){
      mx = lx;
      my = ly;
      mz = lz;
   }
   if(!(nL > 0.0) && !(nR > 0.0)){
      mx = 0.0;
      my = 0.0;
      mz = 1.0;
   }
   seg.mx = mx;
   seg.my = my;
   seg.mz = mz;
}

static std::vector<ChargeBarrier> find_charge_barriers(const int nf, const double* sigma, const int* mat){
   std::vector<ChargeBarrier> barriers;
   if(nf <= 1 || sigma == nullptr) return barriers;
   std::vector<char> owned(static_cast<std::size_t>(nf), 0);
   int i = 0;
   while(i < nf){
      if(sigma[i] == 0.0){
         const int b = i;
         while(i < nf && sigma[i] == 0.0){
            ++i;
         }
         const int c = i - 1;
         ChargeBarrier seg;
         seg.vL = (b > 0) ? b - 1 : 0;
         seg.vR = (c + 1 < nf) ? c + 1 : nf - 1;
         seg.cell_lo = b;
         seg.cell_hi = c;
         seg.n_cells = c - b + 1;
         seg.thick = true;
         seg.face_lo = std::max(b, 1);
         seg.face_hi = std::min(c + 1, nf - 1);
         seg.jump_face = seg.face_hi;
         seg.Gup = 0.0;
         seg.Gdn = 0.0;
         seg.lup = 0.0;
         seg.ldn = 0.0;
         seg.ge = 0.0;
         seg.gp = 0.0;
         seg.li = 0.0;
         seg.pe = 0.0;
         seg.pp = 0.0;
         seg.G = 0.0;
         seg.P = 0.0;
         seg.sfrac = 0.0;
         seg.mx = 0.0;
         seg.my = 0.0;
         seg.mz = 1.0;
         const int matL = electrode_material(mat, seg.vL, -1, nf);
         const int matR = electrode_material(mat, seg.vR, +1, nf);
         lookup_tunnel(matL, matR, seg.Gup, seg.Gdn, seg.lup, seg.ldn);
         load_inelastic(matL, matR, seg);
         if(seg.face_lo <= seg.face_hi){
            for(int e = seg.face_lo; e <= seg.face_hi; ++e) owned[static_cast<std::size_t>(e)] = 1;
            barriers.push_back(seg);
         }
      } else {
         ++i;
      }
   }
   if(!tunnel_enabled || mat == nullptr) return barriers;
   for(int e = 1; e < nf; ++e){
      if(owned[static_cast<std::size_t>(e)]) continue;
      if(!(sigma[e - 1] > 0.0) || !(sigma[e] > 0.0)) continue;
      double gup = 0.0;
      double gdn = 0.0;
      double lup = 0.0;
      double ldn = 0.0;
      if(!lookup_tunnel(electrode_material(mat, e - 1, -1, nf),
                        electrode_material(mat, e, +1, nf), gup, gdn, lup, ldn)) continue;
      ChargeBarrier seg;
      seg.vL = e - 1;
      seg.vR = e;
      seg.cell_lo = -1;
      seg.cell_hi = -1;
      seg.n_cells = 0;
      seg.thick = false;
      seg.face_lo = e;
      seg.face_hi = e;
      seg.jump_face = e;
      seg.Gup = gup;
      seg.Gdn = gdn;
      seg.lup = lup;
      seg.ldn = ldn;
      load_inelastic(electrode_material(mat, e - 1, -1, nf),
                     electrode_material(mat, e, +1, nf), seg);
      seg.G = 0.0;
      seg.P = 0.0;
      seg.sfrac = 0.0;
      seg.mx = 0.0;
      seg.my = 0.0;
      seg.mz = 1.0;
      barriers.push_back(seg);
   }
   return barriers;
}

// One metal interval V[lo+1..hi], with V[i] known for i <= lo. A wide barrier
// that begins at hi is not an unknown: its first oxide node is the linear
// profile V[hi] - (J0/G)/nspan, which folds back onto V[hi]. The block stays
// pentadiagonal. Oxide nodes of an earlier barrier are already in V.
static void solve_metal_charge_block(const int lo,
                                     const int hi,
                                     const int nf,
                                     const double dz,
                                     const double J0,
                                     const double* De_edge,
                                     const double* sigma_edge,
                                     const double* eps_edge,
                                     const double* J_diffusion,
                                     const double* J_seebeck,
                                     const double* J_hot,
                                     const std::vector<int>& face_kind,
                                     const std::vector<int>& face_barrier,
                                     const std::vector<char>& wide_face,
                                     const std::vector<ChargeBarrier>& barriers,
                                     const int nxt_b,
                                     double* V){
   const int n = hi - lo;
   if(n <= 0) return;
   const double idz2 = 1.0 / (dz * dz);
   std::vector<double> al((std::size_t)n, 0.0);
   std::vector<double> bl((std::size_t)n, 0.0);
   std::vector<double> diag((std::size_t)n, 0.0);
   std::vector<double> au((std::size_t)n, 0.0);
   std::vector<double> au2((std::size_t)n, 0.0);
   std::vector<double> rhs((std::size_t)n, 0.0);
   const ChargeBarrier* nxt = (nxt_b >= 0) ? &barriers[static_cast<std::size_t>(nxt_b)] : nullptr;
   const double nspan = nxt ? static_cast<double>(nxt->vR - nxt->vL) : 1.0;
   const double delta = (nxt && nxt->G > 0.0) ? (J0 / nxt->G) : 0.0;

   struct Placer {
      int lo, hi, n, nf;
      const ChargeBarrier* nxt;
      double delta, nspan;
      double* V;
      double* al;
      double* bl;
      double* diag;
      double* au;
      double* au2;
      double* rhs;
      void place(const int row, const int vindex, const double coef){
         if(row < 0 || row >= n || vindex <= 0 || vindex >= nf || coef == 0.0) return;
         if(nxt && vindex == hi + 1){
            place(row, hi, coef);
            rhs[row] += coef * (delta / nspan);
            return;
         }
         if(vindex <= lo){
            rhs[row] -= coef * V[vindex];
            return;
         }
         const int col = vindex - (lo + 1);
         const int span = col - row;
         if(span == -2) al[row] += coef;
         else if(span == -1) bl[row] += coef;
         else if(span == 0) diag[row] += coef;
         else if(span == 1) au[row] += coef;
         else if(span == 2) au2[row] += coef;
         else abort_charge("metal charge block is outside the pentadiagonal stencil");
      }
   };
   Placer place_into{lo, hi, n, nf, nxt, delta, nspan, V,
                     al.data(), bl.data(), diag.data(), au.data(), au2.data(), rhs.data()};

   auto add_ne = [&](const int row, const int cell, const double scale){
      const double aL = eps_edge[cell] * idz2;
      const double aR = eps_edge[cell + 1] * idz2;
      place_into.place(row, cell - 1, scale * (-aL));
      place_into.place(row, cell, scale * (aL + aR));
      place_into.place(row, cell + 1, scale * (-aR));
   };

   for(int row = 0; row < n; ++row){
      const int face = lo + 1 + row;
      const int kind = face_kind[static_cast<std::size_t>(face)];
      if(kind != 0){
         if(wide_face[static_cast<std::size_t>(face)]) abort_charge("wide barrier face entered a metal charge block");
         const ChargeBarrier& seg = barriers[static_cast<std::size_t>(face_barrier[static_cast<std::size_t>(face)])];
         if(kind == 1){
            place_into.place(row, seg.vL, 1.0);
            place_into.place(row, seg.vR, -1.0);
            rhs[static_cast<std::size_t>(row)] += (seg.G > 0.0) ? (J0 / seg.G) : 0.0;
         } else {
            place_into.place(row, face - 1, -1.0);
            place_into.place(row, face, 2.0);
            place_into.place(row, face + 1, -1.0);
         }
         continue;
      }
      const double sig = sigma_edge[face];
      if(!(sig > 0.0)){
         place_into.place(row, face, 1.0);
         place_into.place(row, face - 1, -1.0);
         continue;
      }
      const double r = De_edge[face] / sig;
      place_into.place(row, face - 1, 1.0);
      place_into.place(row, face, -1.0);
      add_ne(row, face, -r);
      add_ne(row, face - 1, r);
      double Jimp = 0.0;
      if(J_diffusion != nullptr) Jimp += J_diffusion[face];
      if(J_seebeck != nullptr) Jimp += J_seebeck[face];
      if(J_hot != nullptr) Jimp += J_hot[face];
      rhs[static_cast<std::size_t>(row)] += (J0 - Jimp) * dz / sig;
   }

   for(int i = 0; i < n; ++i){
      if(i >= 2){
         const double piv = diag[static_cast<std::size_t>(i - 2)];
         const double w = (std::fabs(piv) > 0.0) ? al[static_cast<std::size_t>(i)] / piv : 0.0;
         bl[static_cast<std::size_t>(i)] -= w * au[static_cast<std::size_t>(i - 2)];
         diag[static_cast<std::size_t>(i)] -= w * au2[static_cast<std::size_t>(i - 2)];
         rhs[static_cast<std::size_t>(i)] -= w * rhs[static_cast<std::size_t>(i - 2)];
      }
      if(i >= 1){
         const double piv = diag[static_cast<std::size_t>(i - 1)];
         const double w = (std::fabs(piv) > 0.0) ? bl[static_cast<std::size_t>(i)] / piv : 0.0;
         diag[static_cast<std::size_t>(i)] -= w * au[static_cast<std::size_t>(i - 1)];
         au[static_cast<std::size_t>(i)] -= w * au2[static_cast<std::size_t>(i - 1)];
         rhs[static_cast<std::size_t>(i)] -= w * rhs[static_cast<std::size_t>(i - 1)];
      }
   }
   std::vector<double> sol((std::size_t)n, 0.0);
   for(int i = n - 1; i >= 0; --i){
      double s = rhs[static_cast<std::size_t>(i)];
      if(i + 1 < n) s -= au[static_cast<std::size_t>(i)] * sol[static_cast<std::size_t>(i + 1)];
      if(i + 2 < n) s -= au2[static_cast<std::size_t>(i)] * sol[static_cast<std::size_t>(i + 2)];
      double piv = diag[static_cast<std::size_t>(i)];
      if(std::fabs(piv) < 1e-30) piv = (piv >= 0.0 ? 1e-30 : -1e-30);
      sol[static_cast<std::size_t>(i)] = s / piv;
   }
   for(int i = 0; i < n; ++i) V[lo + 1 + i] = sol[static_cast<std::size_t>(i)];
}

// Thick oxide: drop the interior nodes. Laplace on those faces is the linear
// profile between the electrodes, and the jump is V_L - V_R = J0/G. Each metal
// segment is then a pentadiagonal Thomas solve. Same discrete system as the
// dense factorisation, including the unpolarised inelastic part of G.
static void solve_charge_wide_schur(const int nf,
                                    const double dz,
                                    const double J0,
                                    const double* De_edge,
                                    const double* sigma_edge,
                                    const double* eps_edge,
                                    const double* J_diffusion,
                                    const double* J_seebeck,
                                    const double* J_hot,
                                    const std::vector<ChargeBarrier>& barriers,
                                    double* V){
   for(int i = 0; i < nf; ++i) V[i] = 0.0;
   std::vector<int> face_kind(static_cast<std::size_t>(nf), 0);
   std::vector<int> face_barrier(static_cast<std::size_t>(nf), -1);
   std::vector<char> wide_face(static_cast<std::size_t>(nf), 0);
   std::vector<int> wide_ids;
   for(int b = 0; b < (int)barriers.size(); ++b){
      const ChargeBarrier& seg = barriers[static_cast<std::size_t>(b)];
      const bool is_wide = (seg.vR - seg.vL > 2);
      if(is_wide) wide_ids.push_back(b);
      if(seg.face_lo > seg.face_hi) continue;
      for(int e = seg.face_lo; e <= seg.face_hi; ++e){
         face_barrier[static_cast<std::size_t>(e)] = b;
         face_kind[static_cast<std::size_t>(e)] = (e == seg.jump_face) ? 1 : 2;
         if(is_wide) wide_face[static_cast<std::size_t>(e)] = 1;
      }
   }
   std::sort(wide_ids.begin(), wide_ids.end(), [&](const int a, const int b){
      return barriers[static_cast<std::size_t>(a)].vL < barriers[static_cast<std::size_t>(b)].vL;
   });

   int lo = 0;
   for(std::size_t w = 0; w < wide_ids.size(); ++w){
      const int b = wide_ids[w];
      const ChargeBarrier& seg = barriers[static_cast<std::size_t>(b)];
      solve_metal_charge_block(lo, seg.vL, nf, dz, J0, De_edge, sigma_edge, eps_edge,
                               J_diffusion, J_seebeck, J_hot, face_kind, face_barrier,
                               wide_face, barriers, b, V);
      const double delta = (seg.G > 0.0) ? (J0 / seg.G) : 0.0;
      const int vL = seg.vL;
      const int vR = seg.vR;
      V[vR] = V[vL] - delta;
      const double nspan = static_cast<double>(vR - vL);
      for(int j = 1; j < vR - vL; ++j){
         V[vL + j] = V[vL] + (static_cast<double>(j) / nspan) * (V[vR] - V[vL]);
      }
      lo = vR;
   }
   solve_metal_charge_block(lo, nf - 1, nf, dz, J0, De_edge, sigma_edge, eps_edge,
                            J_diffusion, J_seebeck, J_hot, face_kind, face_barrier,
                            wide_face, barriers, -1, V);
}

// Potentials only. V[0] stays 0. ne and the face currents are rebuilt by the caller
// after the series current is known. A thick barrier is reduced by the Schur
// split above. A thin interface still fits the pentadiagonal stencil.
static void solve_charge_potentials(const int nf,
                                    const double dz,
                                    const double J0,
                                    const double* De_edge,
                                    const double* sigma_edge,
                                    const double* eps_edge,
                                    const double* J_diffusion,
                                    const double* J_seebeck,
                                    const double* J_hot,
                                    const std::vector<ChargeBarrier>& barriers,
                                    double* V){
   if(nf <= 0) return;
   V[0] = 0.0;
   if(nf == 1) return;

   const int m = nf - 1;
   const double idz2 = 1.0 / (dz * dz);
   bool wide = false;
   std::vector<int> face_kind(static_cast<std::size_t>(nf), 0);
   std::vector<int> face_barrier(static_cast<std::size_t>(nf), -1);
   for(int b = 0; b < (int)barriers.size(); ++b){
      const ChargeBarrier& seg = barriers[static_cast<std::size_t>(b)];
      if(seg.vR - seg.vL > 2) wide = true;
      if(seg.face_lo > seg.face_hi) continue;
      for(int e = seg.face_lo; e <= seg.face_hi; ++e){
         face_barrier[static_cast<std::size_t>(e)] = b;
         face_kind[static_cast<std::size_t>(e)] = (e == seg.jump_face) ? 1 : 2;
      }
   }
   if(wide){
      solve_charge_wide_schur(nf, dz, J0, De_edge, sigma_edge, eps_edge,
                              J_diffusion, J_seebeck, J_hot, barriers, V);
      return;
   }

   static std::vector<double> al, bl, diag, au, au2, rhs, x;
   if((int)al.size() != m){
      al.assign((std::size_t)m, 0.0);
      bl.assign((std::size_t)m, 0.0);
      diag.assign((std::size_t)m, 0.0);
      au.assign((std::size_t)m, 0.0);
      au2.assign((std::size_t)m, 0.0);
      rhs.assign((std::size_t)m, 0.0);
      x.assign((std::size_t)m, 0.0);
   } else {
      std::fill(al.begin(), al.end(), 0.0);
      std::fill(bl.begin(), bl.end(), 0.0);
      std::fill(diag.begin(), diag.end(), 0.0);
      std::fill(au.begin(), au.end(), 0.0);
      std::fill(au2.begin(), au2.end(), 0.0);
      std::fill(rhs.begin(), rhs.end(), 0.0);
   }
   auto add_coeff = [&](const int row, const int vindex, const double coef){
      if(row < 0 || row >= m || vindex <= 0 || vindex >= nf) return;
      const int delta = (vindex - 1) - row;
      if(delta == -2) al[row] += coef;
      else if(delta == -1) bl[row] += coef;
      else if(delta == 0) diag[row] += coef;
      else if(delta == 1) au[row] += coef;
      else if(delta == 2) au2[row] += coef;
      else abort_charge("barrier coupling does not fit the pentadiagonal charge matrix");
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
      const int kind = face_kind[static_cast<std::size_t>(face)];
      if(kind != 0){
         const ChargeBarrier& seg = barriers[static_cast<std::size_t>(face_barrier[static_cast<std::size_t>(face)])];
         if(kind == 1){
            add_coeff(row, seg.vL, 1.0);
            add_coeff(row, seg.vR, -1.0);
            rhs[row] = (seg.G > 0.0) ? (J0 / seg.G) : 0.0;
         } else {
            add_coeff(row, face - 1, -1.0);
            add_coeff(row, face, 2.0);
            add_coeff(row, face + 1, -1.0);
            rhs[row] = 0.0;
         }
         continue;
      }
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
      double Jimp = 0.0;
      if(J_diffusion != nullptr) Jimp += J_diffusion[face];
      if(J_seebeck != nullptr) Jimp += J_seebeck[face];
      if(J_hot != nullptr) Jimp += J_hot[face];
      rhs[row] = (J0 - Jimp) * dz / sig;
   }

   {
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
      for(int i = m - 1; i >= 0; --i){
         double s = rhs[i];
         if(i + 1 < m) s -= au[i] * x[i + 1];
         if(i + 2 < m) s -= au2[i] * x[i + 2];
         double piv = diag[i];
         if(std::fabs(piv) < 1e-30) piv = (piv >= 0.0 ? 1e-30 : -1e-30);
         x[i] = s / piv;
      }
      for(int i = 1; i < nf; ++i) V[i] = x[i - 1];
   }
}


// Harmonic mean of the prolonged bulk conductivity. The analytic sheet is not included.
static void bulk_sigma_edge(const double* sigma_bulk, const int nf, double* sigma_edge){
   for(int e = 0; e <= nf; ++e) sigma_edge[e] = 0.0;
   if(sigma_bulk == nullptr) return;
   for(int e = 1; e < nf; ++e){
      sigma_edge[e] = edge_harmonic_mean(sigma_bulk[e - 1], sigma_bulk[e]);
   }
}

// J_seebeck = -sigma * S * dTe/dz. The face value is the arithmetic mean.
// sigma is the bulk face conductivity. Walls stay 0.
static void build_j_seebeck_edge(const int nf,
                                 const double dz,
                                 const double* sigma_edge,
                                 const double* Te_cell,
                                 const double* seebeck_cell,
                                 const double* dTe_dz,
                                 double* J_seebeck){
   (void)dz;
   for(int e = 0; e <= nf; ++e) J_seebeck[e] = 0.0;
   if(!(sc1d_seebeck_enable && sc1d_drive_active())) return;
   for(int e = 1; e < nf; ++e){
      const double sig = sigma_edge[e];
      if(!(sig > 0.0)) continue;
      const double TeL = Te_cell[e - 1];
      const double TeR = Te_cell[e];
      const double TpL = sc1d_Tp_fine[(std::size_t)(e - 1)];
      const double TpR = sc1d_Tp_fine[(std::size_t)e];
      const double SL = scale_seebeck(seebeck_cell[e - 1], TeL, TpL);
      const double SR = scale_seebeck(seebeck_cell[e], TeR, TpR);
      const double S_avg = edge_seebeck_coefficient(SL, SR);
      J_seebeck[e] = -sig * S_avg * dTe_dz[e];
   }
}

// J_hot = -Jns. build_jns_fine stores the electron flux. The minus is positive
// charge toward +z. Walls stay 0.
static void charge_from_electron_flux(const int nf, const double* Jns, double* J_hot){
   if(J_hot == nullptr || nf < 0) return;
   if(Jns == nullptr) return;
   for(int e = 0; e <= nf; ++e) J_hot[e] = -Jns[e];
}

static std::vector<std::vector<ChargeBarrier>> saved_barriers;
static std::vector<double> saved_J0;

static void screen_open_circuit_faces(const int stack,
                                      const int nf,
                                      const double dz,
                                      const double* sigma_cell,
                                      const int* mat_cell,
                                      const double* De_edge,
                                      const double* sigma_edge,
                                      const double* eps_edge,
                                      const double* Bd_cell,
                                      const double* S,
                                      const double* mhat,
                                      const double* Jns,
                                      double* J_hot,
                                      double* J_seebeck,
                                      double* J_spin,
                                      double* J_field,
                                      double* J_cd,
                                      double* ne,
                                      double* V,
                                      double* Jc){
   charge_from_electron_flux(nf, Jns, J_hot);
   build_jknown_fine(nf, dz, De_edge, Bd_cell, S, mhat, Jns, J_spin);

   std::vector<ChargeBarrier> barriers = find_charge_barriers(nf, sigma_cell, mat_cell);
   for(std::size_t b = 0; b < barriers.size(); ++b){
      set_barrier_angle(barriers[b], mhat, mat_cell, stack, nf, dz);
   }

   std::vector<char> barrier_face(static_cast<std::size_t>(nf + 1), 0);
   std::vector<char> barrier_cell(static_cast<std::size_t>(std::max(nf, 0)), 0);
   for(std::size_t b = 0; b < barriers.size(); ++b){
      const ChargeBarrier& seg = barriers[b];
      if(seg.face_lo <= seg.face_hi){
         for(int e = seg.face_lo; e <= seg.face_hi; ++e) barrier_face[static_cast<std::size_t>(e)] = 1;
      }
      if(seg.thick && seg.cell_lo >= 0){
         for(int i = seg.cell_lo; i <= seg.cell_hi; ++i) barrier_cell[static_cast<std::size_t>(i)] = 1;
      }
   }

   std::vector<double> V_part(static_cast<std::size_t>(std::max(nf, 0)), 0.0);
   solve_charge_potentials(nf, dz, 0.0, De_edge, sigma_edge, eps_edge,
                           J_spin, J_seebeck, J_hot, barriers, V_part.data());

   double J0 = 0.0;
   if(!sc1d_bias_enable){
      for(int i = 0; i < nf; ++i) V[i] = V_part[static_cast<std::size_t>(i)];
   } else {
      std::vector<double> V_unit(static_cast<std::size_t>(nf), 0.0);
      solve_charge_potentials(nf, dz, 1.0, De_edge, sigma_edge, eps_edge,
                              nullptr, nullptr, nullptr, barriers, V_unit.data());
      const double denom = V_unit[static_cast<std::size_t>(nf - 1)];
      const double vpart = V_part[static_cast<std::size_t>(nf - 1)];
      const double dV = sc1d_voltage_top - sc1d_voltage_bottom;
      J0 = (dV - vpart) / denom;
      for(int i = 0; i < nf; ++i){
         V[i] = V_part[static_cast<std::size_t>(i)] + J0 * V_unit[static_cast<std::size_t>(i)];
      }
   }

   const double idz2 = (dz > 0.0) ? (1.0 / (dz * dz)) : 0.0;
   for(int i = 0; i < nf; ++i){
      const double aL = eps_edge[i] * idz2;
      const double aR = eps_edge[i + 1] * idz2;
      const double Vm = (i > 0) ? V[i - 1] : 0.0;
      if(i == nf - 1) ne[i] = aL * (V[i] - Vm);
      else ne[i] = aL * (V[i] - Vm) + aR * (V[i] - V[i + 1]);
   }

   if(J_field != nullptr){
      J_field[0] = 0.0;
      J_field[nf] = 0.0;
   }
   if(J_cd != nullptr){
      J_cd[0] = 0.0;
      J_cd[nf] = 0.0;
   }
   if(Jc != nullptr){
      Jc[0] = J0;
      Jc[nf] = J0;
   }
   for(int e = 1; e < nf; ++e){
      if(barrier_face[static_cast<std::size_t>(e)]){
         if(J_field != nullptr) J_field[e] = 0.0;
         if(J_cd != nullptr) J_cd[e] = 0.0;
         if(J_spin != nullptr) J_spin[e] = 0.0;
         if(J_seebeck != nullptr) J_seebeck[e] = 0.0;
         if(J_hot != nullptr) J_hot[e] = 0.0;
         if(Jc != nullptr) Jc[e] = J0;
         continue;
      }
      double cd = 0.0;
      if(De_edge != nullptr && dz > 0.0){
         cd = -De_edge[e] * (ne[e] - ne[e - 1]) / dz;
      }
      if(J_cd != nullptr) J_cd[e] = cd;
      double Jr = 0.0;
      const double sig = sigma_edge[e];
      if(sig > 0.0 && dz > 0.0){
         const double E = -(V[e] - V[e - 1]) / dz;
         Jr = sig * E;
      }
      if(J_field != nullptr) J_field[e] = Jr;
      if(Jc != nullptr){
         const double js = (J_spin != nullptr) ? J_spin[e] : 0.0;
         const double jsee = (J_seebeck != nullptr) ? J_seebeck[e] : 0.0;
         const double jhot = (J_hot != nullptr) ? J_hot[e] : 0.0;
         Jc[e] = Jr + cd + js + jsee + jhot;
      }
   }

   double* joule = sc1d_joule_fine.data() + static_cast<std::size_t>(stack) * static_cast<std::size_t>(nf);
   for(int i = 0; i < nf; ++i) joule[i] = 0.0;
   for(std::size_t b = 0; b < barriers.size(); ++b){
      const ChargeBarrier& seg = barriers[b];
      const double drop = V[seg.vL] - V[seg.vR];
      if(seg.thick && seg.n_cells > 0 && dz > 0.0){
         const double p = std::fabs(J0 * drop) / (static_cast<double>(seg.n_cells) * dz);
         for(int i = seg.cell_lo; i <= seg.cell_hi; ++i) joule[i] = p;
      } else if(!seg.thick && dz > 0.0){
         const double half = 0.5 * std::fabs(J0 * drop) / dz;
         joule[seg.vL] += half;
         joule[seg.vR] += half;
      }
   }
      for(int i = 0; i < nf; ++i){
         if(barrier_cell[static_cast<std::size_t>(i)]) continue;
         double acc = 0.0;
         int nface_p = 0;
         const int faces[2] = {i, i + 1};
         for(int f = 0; f < 2; ++f){
            const int e = faces[f];
            if(e < 1 || e >= nf) continue;
            if(barrier_face[static_cast<std::size_t>(e)]) continue;
            if(!(sigma_edge[e] > 0.0) || !(dz > 0.0)) continue;
            const double E = -(V[e] - V[e - 1]) / dz;
            const double pe = sigma_cell[i] * E * E;
            acc += pe;
            ++nface_p;
         }
         if(nface_p > 0) joule[i] += acc / static_cast<double>(nface_p);
      }

   if(sc1d_bias_enable){
      for(int i = 0; i < nf; ++i) V[i] += sc1d_voltage_bottom;
   }

   saved_J0[static_cast<std::size_t>(stack)] = J0;
   saved_barriers[static_cast<std::size_t>(stack)] = barriers;
}

// One column of charge workspace. The grid length is fixed after initialisation,
// so the same storage is reused by the steady solve and by every charge step.
struct FineChargeScratch {
   std::vector<double> Te, De, De_hot, ve, sigma;
   std::vector<double> De_edge, De_hot_edge, ve_edge, sigma_edge, eps_edge;
   std::vector<double> a, b, c, rhs, Jns, sigma_see;
   int nf = -1;
   void fit(const int n){
      if(nf == n) return;
      nf = n;
      Te.resize((std::size_t)n);
      De.resize((std::size_t)n);
      De_hot.resize((std::size_t)n);
      ve.resize((std::size_t)n);
      sigma.resize((std::size_t)n);
      De_edge.resize((std::size_t)n + 1u);
      De_hot_edge.resize((std::size_t)n + 1u);
      ve_edge.resize((std::size_t)n + 1u);
      sigma_edge.resize((std::size_t)n + 1u);
      eps_edge.resize((std::size_t)n + 1u);
      a.resize((std::size_t)n);
      b.resize((std::size_t)n);
      c.resize((std::size_t)n);
      rhs.resize((std::size_t)n);
      Jns.resize((std::size_t)n + 1u);
      sigma_see.resize((std::size_t)n + 1u);
   }
};

static FineChargeScratch charge_scratch;

static void solve_open_circuit_charge_steady(const int stack,
                                             double* ns,
                                             double* ne,
                                             double* V,
                                             double* J_seebeck,
                                             double* J_field,
                                             double* J_hot,
                                             double* J_cd,
                                             double* J_spin,
                                             double* Jc_edge,
                                             const int nf,
                                             const double dz,
                                             const double t_s,
                                             const double* sigma_bulk,
                                             const double* D_charge,
                                             const double* D_hot,
                                             const double* sigma0_cell,
                                             const double* Bd_cell,
                                             const double* seebeck_cell,
                                             const double* S,
                                             const double* mhat){
   for(int i=0;i<nf;i++) ns[i] = 0.0;

   fill_fine_charge_transport(nf, dz, t_s, D_charge, D_hot, sigma0_cell,
                              charge_scratch.Te.data(), charge_scratch.De.data(),
                              charge_scratch.De_hot.data(),
                              charge_scratch.ve.data(), charge_scratch.sigma.data());
   fill_fine_charge_faces(nf, charge_scratch.De.data(), charge_scratch.De_hot.data(),
                          charge_scratch.ve.data(), charge_scratch.sigma.data(),
                          charge_scratch.De_edge.data(), charge_scratch.De_hot_edge.data(),
                          charge_scratch.ve_edge.data(),
                          charge_scratch.sigma_edge.data(), charge_scratch.eps_edge.data());
   bulk_sigma_edge(sigma_bulk, nf, charge_scratch.sigma_see.data());
   build_j_seebeck_edge(nf, dz, charge_scratch.sigma_see.data(), charge_scratch.Te.data(), seebeck_cell,
                        sc1d_dTe_dz_edge.data(), J_seebeck);
   build_jns_fine(nf, dz, charge_scratch.De_hot_edge.data(), charge_scratch.ve_edge.data(), ns, charge_scratch.Jns.data());
   const int* mat_cell = sc1d_mat_fine.data() + static_cast<std::size_t>(stack) * static_cast<std::size_t>(nf);
   screen_open_circuit_faces(stack, nf, dz, charge_scratch.sigma.data(), mat_cell,
                             charge_scratch.De_edge.data(), charge_scratch.sigma_edge.data(),
                             charge_scratch.eps_edge.data(),
                             Bd_cell, S, mhat, charge_scratch.Jns.data(),
                             J_hot, J_seebeck, J_spin, J_field, J_cd, ne, V, Jc_edge);
}

static void update_charge_transients_1d_stack(const int stack,
                                              double* ns,
                                              double* ne,
                                              double* V,
                                              double* Jc_edge,
                                              double* J_seebeck,
                                              double* J_field,
                                              double* J_hot,
                                              double* J_cd,
                                              double* J_spin,
                                              const int nf,
                                              const double dz,
                                              const double dt,
                                              const double t_s,
                                              const double* sigma_bulk,
                                              const double* D_charge,
                                              const double* D_hot,
                                              const double* sigma0_cell,
                                              const double* Bd_cell,
                                              const double* seebeck_cell,
                                              const double* S,
                                              const double* mhat){
   fill_fine_charge_transport(nf, dz, t_s, D_charge, D_hot, sigma0_cell,
                              charge_scratch.Te.data(), charge_scratch.De.data(),
                              charge_scratch.De_hot.data(),
                              charge_scratch.ve.data(), charge_scratch.sigma.data());
   fill_fine_charge_faces(nf, charge_scratch.De.data(), charge_scratch.De_hot.data(),
                          charge_scratch.ve.data(), charge_scratch.sigma.data(),
                          charge_scratch.De_edge.data(), charge_scratch.De_hot_edge.data(),
                          charge_scratch.ve_edge.data(),
                          charge_scratch.sigma_edge.data(), charge_scratch.eps_edge.data());
   assemble_ns_fine_tridiagonal(nf, dt, dz, t_s, charge_scratch.De_hot_edge.data(), charge_scratch.ve_edge.data(),
                                ns, charge_scratch.Te.data(), charge_scratch.a, charge_scratch.b,
                                charge_scratch.c, charge_scratch.rhs);
   thomas_solve(charge_scratch.a, charge_scratch.b, charge_scratch.c, charge_scratch.rhs);
   for(int i=0;i<nf;i++) ns[i] = charge_scratch.rhs[i];

   build_jns_fine(nf, dz, charge_scratch.De_hot_edge.data(), charge_scratch.ve_edge.data(), ns, charge_scratch.Jns.data());
   bulk_sigma_edge(sigma_bulk, nf, charge_scratch.sigma_see.data());
   build_j_seebeck_edge(nf, dz, charge_scratch.sigma_see.data(), charge_scratch.Te.data(), seebeck_cell,
                        sc1d_dTe_dz_edge.data(), J_seebeck);
   const int* mat_cell = sc1d_mat_fine.data() + static_cast<std::size_t>(stack) * static_cast<std::size_t>(nf);
   screen_open_circuit_faces(stack, nf, dz, charge_scratch.sigma.data(), mat_cell,
                             charge_scratch.De_edge.data(), charge_scratch.sigma_edge.data(),
                             charge_scratch.eps_edge.data(),
                             Bd_cell, S, mhat, charge_scratch.Jns.data(),
                             J_hot, J_seebeck, J_spin, J_field, J_cd, ne, V, Jc_edge);
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
   // Column s belongs to rank (s % nprocs). The atom domain is a different
   // split. Magnetisation is prolonged from the atoms, and the ltmp profile
   // is the same on every rank, so the owner can step the column without
   // holding its atoms. The cell fields are broadcast afterwards.
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
      zlog << zTs() << "1D spin diffusion "
           << (sc1d_diffusion == sc1d_diff_bdf2 ? "bdf2" : "backward-euler")
           << std::endl;
   }
}

static std::vector<Block3> step_A_face;
static std::vector<Block3> step_L_blk;
static std::vector<Block3> step_D_blk;
static std::vector<Block3> step_U_blk;
static std::vector<double> step_Jd_edge;

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
   if(sc1d_diffusion == sc1d_diff_bdf2){
      if(sc1d_S_bdf2_h1.size() != expected_size) sc1d_S_bdf2_h1.assign(expected_size, 0.0);
      if(sc1d_S_bdf2_h2.size() != expected_size) sc1d_S_bdf2_h2.assign(expected_size, 0.0);
      sc1d_bdf2_primed = false;
   } else {
      sc1d_S_bdf2_h1.clear();
      sc1d_S_bdf2_h2.clear();
   }

   const int nf = sc1d_nf;
   const std::size_t nfn = nstack * (std::size_t)nf;
   sc1d_ns_fine.assign(nfn, 0.0);
   sc1d_ne_fine.assign(nfn, 0.0);
   sc1d_V_fine.assign(nfn, 0.0);
   sc1d_Jc_edge_fine.assign(nstack * (std::size_t)(nf+1), 0.0);
   sc1d_joule_fine.assign(nfn, 0.0);
   sc1d_mat_fine.assign(nfn, -1);
   sc1d_sigma_fine.assign(nfn, 0.0);
   sc1d_sigma_bulk_fine.assign(nfn, 0.0);
   sc1d_Rint_face.assign(nfn, 0.0);
   sc1d_dmdt_fine.assign(nfn, 0.0);
   sc1d_seebeck_fine.assign(nfn, 0.0);
   sc1d_beta_seebeck_fine.assign(nfn, 0.0);
   sc1d_phot_fine.assign(nfn, 0.0);
   sc1d_J_seebeck_edge.assign(nstack * (std::size_t)(nf+1), 0.0);
   sc1d_J_restore_edge.assign(nstack * (std::size_t)(nf+1), 0.0);
   sc1d_J_hot_edge.assign(nstack * (std::size_t)(nf+1), 0.0);
   sc1d_J_cd_edge.assign(nstack * (std::size_t)(nf+1), 0.0);
   sc1d_J_sd_edge.assign(nstack * (std::size_t)(nf+1), 0.0);
   sc1d_Js_fine.assign(nstack * (std::size_t)nf * 3u, 0.0);
   sc1d_mhat_fine.assign(nstack * (std::size_t)nf * 3u, 0.0);
   sc1d_Te_fine.assign((std::size_t)nf, sc1d_reference_temperature);
   sc1d_Tp_fine.assign((std::size_t)nf, sc1d_reference_temperature);
   sc1d_dTe_dz_edge.assign((std::size_t)(nf + 1), 0.0);
   sc1d_dTp_dz_edge.assign((std::size_t)(nf + 1), 0.0);

   sc1d_Bc_fine.assign(nfn, 0.0);
   sc1d_Bd_fine.assign(nfn, 0.0);
   sc1d_D_fine.assign(nfn, 0.0);
   sc1d_Dmat_fine.assign(nfn, 0.0);
   sc1d_lsf_fine.assign(nfn, 0.0);
   sc1d_lphi_fine.assign(nfn, 0.0);
   sc1d_Jsd_fine.assign(nfn, 0.0);
   sc1d_chi_fine.assign(nfn, 0.0);
   sc1d_sa_inf_fine.assign(nfn, 0.0);
   saved_J0.assign(nstack, 0.0);
   saved_barriers.assign(nstack, std::vector<ChargeBarrier>());
   charge_scratch.fit(nf);
   const int nface = (nf > 1) ? (nf - 1) : 0;
   step_A_face.assign(static_cast<std::size_t>(nface), block3_zero());
   step_L_blk.assign(static_cast<std::size_t>(std::max(0, nf)), block3_zero());
   step_D_blk.assign(static_cast<std::size_t>(std::max(0, nf)), block3_zero());
   step_U_blk.assign(static_cast<std::size_t>(std::max(0, nf)), block3_zero());
   step_Jd_edge.assign(3u * static_cast<std::size_t>(nf + 1), 0.0);
}


// R_half = (dz/(2D)) [ I + (gamma/(1-gamma)) m m^T ].
// Perpendicular eigenvalues are dz/(2D). The parallel one is that over (1-gamma).
static Block3 half_cell_resistance(const double D,
                                  const double gamma,
                                  const double dz,
                                  const double mx,
                                  const double my,
                                  const double mz){
   const double s = dz / (2.0 * D);
   const double gfac = gamma / (1.0 - gamma);
   const double sm = s * gfac;
   Block3 R;
   R.a00 = s + sm*mx*mx;
   R.a01 = sm*mx*my;
   R.a02 = sm*mx*mz;
   R.a10 = R.a01;
   R.a11 = s + sm*my*my;
   R.a12 = sm*my*mz;
   R.a20 = R.a02;
   R.a21 = R.a12;
   R.a22 = s + sm*mz*mz;
   return R;
}

// A_face[i] is G/dz between fine cells i and i+1, with G = inv(R_face) in m/s.
// Same series split as the old scalar alpha. A closed face (D below kEps) is
// zero and is not inverted. Rint_face is the baked analytic Robin resistance.
static void assemble_diffusion_faces(const int nf,
                                    const double dz,
                                    const double* D_fine,
                                    const double* Bc_fine,
                                    const double* Bd_fine,
                                    const double* mhat,
                                    const double* Rint_face,
                                    std::vector<Block3>& A_face){
   if(nf < 2) return;
   const double idz = 1.0 / dz;
   for(int i=0;i<nf-1;i++){
      const double Di = D_fine[i];
      const double Dj = D_fine[i+1];
      if(Di < kEps || Dj < kEps){
         A_face[(std::size_t)i] = block3_zero();
         continue;
      }
      double Rint = 0.0;
      if(Rint_face != nullptr) Rint = Rint_face[i];
      // The smallest eigenvalue is the old scalar series resistance. A non-positive
      // value is a closed face, same as the previous alpha = 0 branch.
      const double Rseries = dz/(2.0*Di) + Rint + dz/(2.0*Dj);
      if(Rseries <= 0.0){
         A_face[(std::size_t)i] = block3_zero();
         continue;
      }
      const double gi = spin_gamma(Bc_fine[i], Bd_fine[i]);
      const double gj = spin_gamma(Bc_fine[i+1], Bd_fine[i+1]);
      Block3 R = block3_add(
         half_cell_resistance(Di, gi, dz, mhat[3*i+0], mhat[3*i+1], mhat[3*i+2]),
         half_cell_resistance(Dj, gj, dz, mhat[3*(i+1)+0], mhat[3*(i+1)+1], mhat[3*(i+1)+2]));
      R.a00 += Rint;
      R.a11 += Rint;
      R.a22 += Rint;
      A_face[(std::size_t)i] = block3_scale(block3_inv(R), idz);
   }
}

static void add_relax(Block3& A, const double relax){
   A.a00 += relax;
   A.a11 += relax;
   A.a22 += relax;
}

// Steady operator: D_i S_i - A_{i-1} S_{i-1} - A_i S_{i+1}, plus relax on D.
static void assemble_steady_blocks(const int nf,
                                  const double* relax,
                                  const std::vector<Block3>& A_face,
                                  std::vector<Block3>& L,
                                  std::vector<Block3>& D,
                                  std::vector<Block3>& U){
   const Block3 Z = block3_zero();
   if(nf <= 0) return;
   if(nf == 1){
      L[0] = Z;
      D[0] = block3_zero();
      add_relax(D[0], relax[0]);
      U[0] = Z;
      return;
   }
   L[0] = Z;
   D[0] = A_face[0];
   add_relax(D[0], relax[0]);
   U[0] = block3_scale(A_face[0], -1.0);
   for(int i=1;i<nf-1;i++){
      L[(std::size_t)i] = block3_scale(A_face[(std::size_t)(i-1)], -1.0);
      D[(std::size_t)i] = block3_add(A_face[(std::size_t)(i-1)], A_face[(std::size_t)i]);
      add_relax(D[(std::size_t)i], relax[i]);
      U[(std::size_t)i] = block3_scale(A_face[(std::size_t)i], -1.0);
   }
   const int last = nf - 1;
   L[(std::size_t)last] = block3_scale(A_face[(std::size_t)(last-1)], -1.0);
   D[(std::size_t)last] = A_face[(std::size_t)(last-1)];
   add_relax(D[(std::size_t)last], relax[last]);
   U[(std::size_t)last] = Z;
}

// Backward Euler: I + dt_half * A, Neumann ends. Built once per stack and
// reused by both Strang half-steps.
static void assemble_diffusion_half_blocks(const int nf,
                                          const double dt_half,
                                          const std::vector<Block3>& A_face,
                                          std::vector<Block3>& L,
                                          std::vector<Block3>& D,
                                          std::vector<Block3>& U){
   const Block3 Z = block3_zero();
   const Block3 I = block3_identity();
   if(nf <= 0) return;
   if(nf == 1){
      L[0] = Z;
      D[0] = I;
      U[0] = Z;
      return;
   }
   const Block3 dt0 = block3_scale(A_face[0], dt_half);
   L[0] = Z;
   D[0] = block3_add(I, dt0);
   U[0] = block3_scale(dt0, -1.0);
   for(int i=1;i<nf-1;i++){
      const Block3 dtL = block3_scale(A_face[(std::size_t)(i-1)], dt_half);
      const Block3 dtR = block3_scale(A_face[(std::size_t)i], dt_half);
      L[(std::size_t)i] = block3_scale(dtL, -1.0);
      D[(std::size_t)i] = block3_add(I, block3_add(dtL, dtR));
      U[(std::size_t)i] = block3_scale(dtR, -1.0);
   }
   const int last = nf - 1;
   const Block3 dtL = block3_scale(A_face[(std::size_t)(last-1)], dt_half);
   L[(std::size_t)last] = block3_scale(dtL, -1.0);
   D[(std::size_t)last] = block3_add(I, dtL);
   U[(std::size_t)last] = Z;
}

static void initialise_fine_spin_if_empty(const int nf,
                                         const double dz,
                                         const double* Bc_fine,
                                         const double* Bd_fine,
                                         const double* D_fine,
                                         const double* lsf_fine,
                                         const double* sa_inf_fine,
                                         const double* mhat,
                                         const double* Rint_face,
                                         double* Sbase){
   if(nf <= 0) return;

   std::vector<Vec3> sa_eq((std::size_t)nf);
   std::vector<double> relax((std::size_t)nf, 0.0);
   for(int i=0;i<nf;i++){
      const double sa = sa_inf_fine[i];
      sa_eq[(std::size_t)i].x = sa * mhat[3*i+0];
      sa_eq[(std::size_t)i].y = sa * mhat[3*i+1];
      sa_eq[(std::size_t)i].z = sa * mhat[3*i+2];
      const double lambda = lsf_fine[i];
      if(lambda > 0.0 && D_fine[i] > 0.0){
         const double gamma = spin_gamma(Bc_fine[i], Bd_fine[i]);
         const double lambda_sf = lambda / std::sqrt(1.0 - gamma);
         relax[(std::size_t)i] = D_fine[i] / (lambda_sf * lambda_sf);
      }
   }

   // Steady diffusion-relaxation at Jc = 0, on the same faces as the half-step.
   // lambda_sf = lambda_sdl / sqrt(1-gamma), so the decay length is lambda_sdl.
   // sa_eq follows the same unit axis as the time stepper. Neumann ends.
   // A zero diagonal (no diffusion and no relaxation) falls back to sa_eq.
   const int nface = (nf > 1) ? (nf - 1) : 0;
   std::vector<Block3> A_face((std::size_t)nface, block3_zero());
   assemble_diffusion_faces(nf, dz, D_fine, Bc_fine, Bd_fine, mhat, Rint_face, A_face);

   std::vector<Block3> L((std::size_t)nf, block3_zero());
   std::vector<Block3> Dblk((std::size_t)nf, block3_zero());
   std::vector<Block3> U((std::size_t)nf, block3_zero());
   assemble_steady_blocks(nf, relax.data(), A_face, L, Dblk, U);

   std::vector<Vec3> rhs((std::size_t)nf);
   const Block3 Z = block3_zero();
   const Block3 I = block3_identity();
   for(int i=0;i<nf;i++){
      const double rel = relax[(std::size_t)i];
      rhs[(std::size_t)i].x = rel * sa_eq[(std::size_t)i].x;
      rhs[(std::size_t)i].y = rel * sa_eq[(std::size_t)i].y;
      rhs[(std::size_t)i].z = rel * sa_eq[(std::size_t)i].z;
      const double tr = Dblk[(std::size_t)i].a00 + Dblk[(std::size_t)i].a11 + Dblk[(std::size_t)i].a22;
      if(std::fabs(tr) < kEps){
         L[(std::size_t)i] = Z;
         Dblk[(std::size_t)i] = I;
         U[(std::size_t)i] = Z;
         rhs[(std::size_t)i] = sa_eq[(std::size_t)i];
      }
   }
   block3_thomas(L, Dblk, U, rhs);
   for(int i=0;i<nf;i++){
      Sbase[3*i+0] = rhs[(std::size_t)i].x;
      Sbase[3*i+1] = rhs[(std::size_t)i].y;
      Sbase[3*i+2] = rhs[(std::size_t)i].z;
   }
}

static void build_atom_fine_map(const int nf, const double dzA){
   sc1d_atom_ls.assign(num_local_atoms, -1);
   sc1d_atom_fine.assign(num_local_atoms, -1);
   sc1d_atom_fine_lo.assign(num_local_atoms, -1);
   sc1d_atom_fine_hi.assign(num_local_atoms, -1);
   sc1d_plane_z.clear();
   sc1d_plane_lo.clear();
   sc1d_plane_hi.clear();

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
   sc1d_plane_z = zplane;
   sc1d_plane_lo.assign(static_cast<std::size_t>(nplane), -1);
   sc1d_plane_hi.assign(static_cast<std::size_t>(nplane), -1);

   for(int atom=0; atom<num_local_atoms; ++atom){
      const int cell = atom_st_index[atom];
      if(cell < 0 || cell >= (int)cell_stack_index.size()) continue;
      const int stack = cell_stack_index[cell] - 1;
      if(stack < 0 || stack >= num_stacks_y) continue;

      const double za = atoms::z_coord_array[atom];
      int fine = static_cast<int>(std::floor(za / dzA));
      if(fine < 0) fine = 0;
      if(fine >= nf) fine = nf - 1;
      sc1d_atom_ls[atom] = stack;
      sc1d_atom_fine[atom] = fine;

      const int mat = atoms::type_array[atom];
      if(mat >= 0 && mat < (int)mp.size()){
         // Zero sd exchange marks a material with no s-d coupling. A negative
         // value is a reversed exchange and still maps onto the fine grid.
         if(mp[mat].sd_exchange == 0.0) continue;
      }

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

      sc1d_atom_fine_lo[atom] = i_lo;
      sc1d_atom_fine_hi[atom] = i_hi;
      if(plane >= 0 && plane < nplane){
         sc1d_plane_lo[static_cast<std::size_t>(plane)] = i_lo;
         sc1d_plane_hi[static_cast<std::size_t>(plane)] = i_hi;
      }
   }
}

static void fail_nonfinite(const char* name, const std::size_t index, const double value){
   terminaltextcolor(RED);
   std::cerr << "Error - non-finite " << name << " at index " << index
             << " = " << value
             << " after spin-current initialisation. The running solver does not replace this."
             << std::endl;
   terminaltextcolor(WHITE);
   zlog << zTs() << "Error - non-finite " << name << " at index " << index
        << " = " << value << " after spin-current initialisation." << std::endl;
   err::vexit();
}

static void require_finite_values(const char* name, const std::vector<double>& field){
   for(std::size_t i = 0; i < field.size(); ++i){
      if(std::isfinite(field[i])) continue;
      fail_nonfinite(name, i, field[i]);
   }
}

static void require_resolved_conductivity(){
   for(std::size_t i = 0; i < sc1d_sigma_fine.size(); ++i){
      const double sigma = sc1d_sigma_fine[i];
      if(!std::isfinite(sigma)) fail_nonfinite("sigma", i, sigma);
      if(sigma < 0.0){
         terminaltextcolor(RED);
         std::cerr << "Error - unresolved conductivity sentinel " << sigma
                   << " on fine cell " << i << "." << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Error - unresolved conductivity sentinel " << sigma
              << " on fine cell " << i << "." << std::endl;
         err::vexit();
      }
   }
}

static void require_spin_channels(){
   const std::size_t n = sc1d_Bc_fine.size();
   for(std::size_t i = 0; i < n; ++i){
      const double Bc = sc1d_Bc_fine[i];
      const double Bd = sc1d_Bd_fine[i];
      const double g = Bc * Bd;
      const bool polarisation = !(std::fabs(Bc) <= 1.0) || !(std::fabs(Bd) <= 1.0);
      const bool channel = !(g < 1.0) || !std::isfinite(Bc) || !std::isfinite(Bd) || !std::isfinite(g);
      if(!polarisation && !channel) continue;
      terminaltextcolor(RED);
      std::cerr << "Error - spin polarisation beta_c = " << Bc
                << ", beta_d = " << Bd
                << ", 1 - beta_c*beta_d = " << (1.0 - g)
                << " on fine cell " << i << ". ";
      if(polarisation) std::cerr << "Each polarisation must lie in [-1, 1]." << std::endl;
      else std::cerr << "The parallel diffusivity is not positive." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error - beta_c = " << Bc << ", beta_d = " << Bd
           << " on fine cell " << i << " is not a physical spin-diffusion channel." << std::endl;
      err::vexit();
   }
}

// Two passes, once, after the grid, the magnetisation prolongation, and the
// steady charge solve. The time step does not repeat these tests.
static void smoke_test_sc1d(){
   if(sc1d_bias_enable && sc1d_nf < 2){
      terminaltextcolor(RED);
      std::cerr << "Error - applied voltage needs at least two fine cells." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error - applied voltage needs at least two fine cells." << std::endl;
      err::vexit();
   }

   // Pass 1. Prolonged constants and the magnetisation profile.
   require_finite_values("beta_c", sc1d_Bc_fine);
   require_finite_values("beta_d", sc1d_Bd_fine);
   require_finite_values("D", sc1d_D_fine);
   require_finite_values("De", sc1d_Dmat_fine);
   require_finite_values("lambda_sdl", sc1d_lsf_fine);
   require_finite_values("lambda_phi", sc1d_lphi_fine);
   require_finite_values("Jsd", sc1d_Jsd_fine);
   require_finite_values("chi", sc1d_chi_fine);
   require_finite_values("sa_inf", sc1d_sa_inf_fine);
   require_finite_values("sigma_bulk", sc1d_sigma_bulk_fine);
   require_finite_values("seebeck", sc1d_seebeck_fine);
   require_finite_values("beta_seebeck", sc1d_beta_seebeck_fine);
   require_finite_values("P_hot", sc1d_phot_fine);
   require_finite_values("msat", sc1d_msat_fine);
   require_finite_values("mhat", sc1d_mhat_fine);
   require_finite_values("dmdt", sc1d_dmdt_fine);
   require_finite_values("anchor weight", sc1d_anchor_t);
   require_finite_values("Te", sc1d_Te_fine);
   require_finite_values("Tp", sc1d_Tp_fine);
   require_resolved_conductivity();
   require_spin_channels();

   // Pass 2. Faces and the steady charge state written during initialisation.
   require_finite_values("Rint", sc1d_Rint_face);
   require_finite_values("dTe/dz", sc1d_dTe_dz_edge);
   require_finite_values("dTp/dz", sc1d_dTp_dz_edge);
   require_finite_values("Jc", sc1d_Jc_edge_fine);
   require_finite_values("J_seebeck", sc1d_J_seebeck_edge);
   require_finite_values("J_field", sc1d_J_restore_edge);
   require_finite_values("J_hot", sc1d_J_hot_edge);
   require_finite_values("J_charge", sc1d_J_cd_edge);
   require_finite_values("J_spin", sc1d_J_sd_edge);
   require_finite_values("V", sc1d_V_fine);
   require_finite_values("ns", sc1d_ns_fine);
   require_finite_values("ne", sc1d_ne_fine);
   require_finite_values("S", sc1d_Sfine);
   require_finite_values("Joule", sc1d_joule_fine);
}

void initialise_spincurrents_1d(){
   const int nsub = compute_fine_subdivisions();
   sc1d_nsub = nsub;
   sc1d_nf   = nsub * num_microcells_per_stack;

   collect_local_stacks();
   allocate_sc1d_arrays();

   const int nf = sc1d_nf;
   const double dz = (micro_cell_thickness * 1.0e-10) / (double)nsub;
   const double dzA = micro_cell_thickness / (double)nsub;
   configure_temperature_source();
   refresh_fine_temperature_profile();
   build_atom_fine_map(nf, dzA);
   build_fine_material_grid();
   cache_prolong_maps();
   prolong_fine_magnetisation(mp::dt_SI);

   // Constants are baked. Only the owner integrates that column.
   for(int stack=0; stack<num_stacks_y; ++stack){
      if(!rank_owns_column(stack)) continue;

      double* Bc_fine = sc1d_Bc_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* Bd_fine = sc1d_Bd_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* D_fine = sc1d_D_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* lsf_fine = sc1d_lsf_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* sa_inf_fine = sc1d_sa_inf_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* sigma_fine = sc1d_sigma_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* seebeck_fine = sc1d_seebeck_fine.data() + (std::size_t)stack * (std::size_t)nf;
      const double* mhat = sc1d_mhat_fine.data() + (std::size_t)stack * (std::size_t)nf * 3u;
      const double* Rint = sc1d_Rint_face.data() + (std::size_t)stack * (std::size_t)nf;

      double* Sbase = sc1d_Sfine.data() + (std::size_t)stack * (std::size_t)nf * 3u;
      initialise_fine_spin_if_empty(nf, dz, Bc_fine, Bd_fine, D_fine, lsf_fine,
                                    sa_inf_fine, mhat, Rint, Sbase);
      double* ns_fine = sc1d_ns_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* ne_fine = sc1d_ne_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* V_fine = sc1d_V_fine.data() + (std::size_t)stack * (std::size_t)nf;
      double* Jc_edge_fine = sc1d_Jc_edge_fine.data() + (std::size_t)stack * (std::size_t)(nf + 1);
      double* J_seebeck = sc1d_J_seebeck_edge.data() + (std::size_t)stack * (std::size_t)(nf + 1);
      double* J_restore = sc1d_J_restore_edge.data() + (std::size_t)stack * (std::size_t)(nf + 1);
      double* J_hot = sc1d_J_hot_edge.data() + (std::size_t)stack * (std::size_t)(nf + 1);
      double* J_cd = sc1d_J_cd_edge.data() + (std::size_t)stack * (std::size_t)(nf + 1);
      double* J_spin = sc1d_J_sd_edge.data() + (std::size_t)stack * (std::size_t)(nf + 1);
      stopwatch_t sw_charge;
      const double interp_before = sc1d_time_interp;
      sw_charge.start();
      const double* Dmat = sc1d_Dmat_fine.data() + (std::size_t)stack * (std::size_t)nf;
      const double* sigma_bulk = sc1d_sigma_bulk_fine.data() + (std::size_t)stack * (std::size_t)nf;
      solve_open_circuit_charge_steady(stack, ns_fine, ne_fine, V_fine, J_seebeck, J_restore,
                                       J_hot, J_cd, J_spin, Jc_edge_fine,
                                       nf, dz, 0.0, sigma_bulk, D_fine, Dmat, sigma_fine, Bd_fine, seebeck_fine,
                                       Sbase, mhat);
      account_charge_time(sw_charge.elapsed_seconds(), interp_before);
   }


   sc1d_atom_use_phonon.assign(num_local_atoms, false);
   for(int atom=0; atom<num_local_atoms; ++atom){
      const int mat = atoms::type_array[atom];
      if(mat >= 0 && mat < (int)mp::material.size()){
         sc1d_atom_use_phonon[atom] = mp::material[mat].couple_to_phonon_temperature;
      }
   }

   cache_sc1d_torque_atoms();
   smoke_test_sc1d();
   output_sc1d_fine_constants();
}
//------------------------------------------------------------------------------
// Per-step 1D spin-accumulation helpers
//------------------------------------------------------------------------------

static void apply_diffusion_half_step(double* S,
                                      const int nf,
                                      const double dz,
                                      const double dt_half,
                                      const double* Bc,
                                      const double* mhat,
                                      const double* Jc_edge_fine,
                                      const std::vector<Block3>& L,
                                      const std::vector<Block3>& D,
                                      const std::vector<Block3>& U){
   if(nf <= 0) return;
   const double jdl0x = -Bc[0]    * Jc_edge_fine[0]  * mhat[0];
   const double jdl0y = -Bc[0]    * Jc_edge_fine[0]  * mhat[1];
   const double jdl0z = -Bc[0]    * Jc_edge_fine[0]  * mhat[2];
   const double jdrNx = -Bc[nf-1] * Jc_edge_fine[nf] * mhat[3*(nf-1)+0];
   const double jdrNy = -Bc[nf-1] * Jc_edge_fine[nf] * mhat[3*(nf-1)+1];
   const double jdrNz = -Bc[nf-1] * Jc_edge_fine[nf] * mhat[3*(nf-1)+2];

   const double Jleft[3]  = { -jdl0x, -jdl0y, -jdl0z };
   const double Jright[3] = { -jdrNx, -jdrNy, -jdrNz };

   std::vector<Vec3> rhs((std::size_t)nf);
   for(int i=0;i<nf;i++){
      rhs[(std::size_t)i].x = S[3*i+0];
      rhs[(std::size_t)i].y = S[3*i+1];
      rhs[(std::size_t)i].z = S[3*i+2];
   }
   rhs[0].x += dt_half * (Jleft[0] / dz);
   rhs[0].y += dt_half * (Jleft[1] / dz);
   rhs[0].z += dt_half * (Jleft[2] / dz);
   rhs[(std::size_t)(nf-1)].x -= dt_half * (Jright[0] / dz);
   rhs[(std::size_t)(nf-1)].y -= dt_half * (Jright[1] / dz);
   rhs[(std::size_t)(nf-1)].z -= dt_half * (Jright[2] / dz);

   block3_thomas(L, D, U, rhs);

   for(int i=0;i<nf;i++){
      S[3*i+0] = rhs[(std::size_t)i].x;
      S[3*i+1] = rhs[(std::size_t)i].y;
      S[3*i+2] = rhs[(std::size_t)i].z;
   }
}

// A_op from the same stencil as assemble_diffusion_half_blocks. A closed face
// (insulator, D = 0, including MgO after prolongation) is a zero block, so it
// drops out here exactly as it does for backward Euler.
static void apply_face_operator(const int nf,
                                const std::vector<Block3>& A_face,
                                const std::vector<Vec3>& S,
                                std::vector<Vec3>& out){
   if(nf <= 0) return;
   out.resize((std::size_t)nf);
   if(nf == 1){
      out[0].x = 0.0;
      out[0].y = 0.0;
      out[0].z = 0.0;
      return;
   }
   const Vec3 left0 = block3_mul_vec(A_face[0], S[0]);
   const Vec3 right0 = block3_mul_vec(A_face[0], S[1]);
   out[0].x = left0.x - right0.x;
   out[0].y = left0.y - right0.y;
   out[0].z = left0.z - right0.z;
   for(int i = 1; i < nf - 1; ++i){
      const Block3 both = block3_add(A_face[(std::size_t)(i - 1)], A_face[(std::size_t)i]);
      const Vec3 from_left = block3_mul_vec(A_face[(std::size_t)(i - 1)], S[(std::size_t)(i - 1)]);
      const Vec3 centre = block3_mul_vec(both, S[(std::size_t)i]);
      const Vec3 to_right = block3_mul_vec(A_face[(std::size_t)i], S[(std::size_t)(i + 1)]);
      out[(std::size_t)i].x = -from_left.x + centre.x - to_right.x;
      out[(std::size_t)i].y = -from_left.y + centre.y - to_right.y;
      out[(std::size_t)i].z = -from_left.z + centre.z - to_right.z;
   }
   const int last = nf - 1;
   const Vec3 from_left = block3_mul_vec(A_face[(std::size_t)(last - 1)], S[(std::size_t)(last - 1)]);
   const Vec3 centre = block3_mul_vec(A_face[(std::size_t)(last - 1)], S[(std::size_t)last]);
   out[(std::size_t)last].x = -from_left.x + centre.x;
   out[(std::size_t)last].y = -from_left.y + centre.y;
   out[(std::size_t)last].z = -from_left.z + centre.z;
}

// Neumann data at the outer ends. Same signs as the backward-Euler half-step.
// Interior closed faces, including the oxide, are not ends of this column.
static void diffusion_boundary_source(const int nf,
                                      const double dz,
                                      const double* Bc,
                                      const double* mhat,
                                      const double* Jc_edge_fine,
                                      std::vector<Vec3>& g){
   g.assign((std::size_t)nf, Vec3());
   if(nf <= 0) return;
   for(int i = 0; i < nf; ++i){
      g[(std::size_t)i].x = 0.0;
      g[(std::size_t)i].y = 0.0;
      g[(std::size_t)i].z = 0.0;
   }
   const double jdl0x = -Bc[0]    * Jc_edge_fine[0]  * mhat[0];
   const double jdl0y = -Bc[0]    * Jc_edge_fine[0]  * mhat[1];
   const double jdl0z = -Bc[0]    * Jc_edge_fine[0]  * mhat[2];
   const double jdrNx = -Bc[nf-1] * Jc_edge_fine[nf] * mhat[3*(nf-1)+0];
   const double jdrNy = -Bc[nf-1] * Jc_edge_fine[nf] * mhat[3*(nf-1)+1];
   const double jdrNz = -Bc[nf-1] * Jc_edge_fine[nf] * mhat[3*(nf-1)+2];
   g[0].x = -jdl0x / dz;
   g[0].y = -jdl0y / dz;
   g[0].z = -jdl0z / dz;
   g[(std::size_t)(nf - 1)].x = jdrNx / dz;
   g[(std::size_t)(nf - 1)].y = jdrNy / dz;
   g[(std::size_t)(nf - 1)].z = jdrNz / dz;
}

// One-step TR-BDF2 on a Strang half-step. The trapezoidal stage saves the
// intermediate state; the BDF2 stage uses only that state and this half-step's
// own input. A previous full step is not history: the explicit precession
// between half-steps is not part of this operator, and extrapolating it
// (4/3 S - 1/3 S_prev) makes a tilted moment run away. With A = 0 the pair of
// stages returns S, including on insulator cells whose faces were closed.
static void apply_diffusion_tr_bdf2_half_step(double* S,
                                              const int nf,
                                              const double dz,
                                              const double dt_half,
                                              const std::vector<Block3>& A_face,
                                              const double* Bc,
                                              const double* mhat,
                                              const double* Jc_edge_fine,
                                              std::vector<Block3>& L,
                                              std::vector<Block3>& D,
                                              std::vector<Block3>& U){
   if(nf <= 0) return;
   const double gamma = 2.0 - std::sqrt(2.0);
   const double c_tr = 0.5 * gamma;
   const double c_bdf = (1.0 - gamma) / (2.0 - gamma);
   const double w_star = 1.0 / (gamma * (2.0 - gamma));
   const double w_n = -((1.0 - gamma) * (1.0 - gamma)) / (gamma * (2.0 - gamma));

   std::vector<Vec3> S0((std::size_t)nf);
   std::vector<Vec3> AS((std::size_t)nf);
   std::vector<Vec3> g((std::size_t)nf);
   std::vector<Vec3> rhs((std::size_t)nf);
   for(int i = 0; i < nf; ++i){
      S0[(std::size_t)i].x = S[3 * i + 0];
      S0[(std::size_t)i].y = S[3 * i + 1];
      S0[(std::size_t)i].z = S[3 * i + 2];
   }
   diffusion_boundary_source(nf, dz, Bc, mhat, Jc_edge_fine, g);
   apply_face_operator(nf, A_face, S0, AS);
   for(int i = 0; i < nf; ++i){
      rhs[(std::size_t)i].x = S0[(std::size_t)i].x - c_tr * dt_half * AS[(std::size_t)i].x + gamma * dt_half * g[(std::size_t)i].x;
      rhs[(std::size_t)i].y = S0[(std::size_t)i].y - c_tr * dt_half * AS[(std::size_t)i].y + gamma * dt_half * g[(std::size_t)i].y;
      rhs[(std::size_t)i].z = S0[(std::size_t)i].z - c_tr * dt_half * AS[(std::size_t)i].z + gamma * dt_half * g[(std::size_t)i].z;
   }
   assemble_diffusion_half_blocks(nf, c_tr * dt_half, A_face, L, D, U);
   block3_thomas(L, D, U, rhs);

   std::vector<Vec3> Sstar(rhs);
   for(int i = 0; i < nf; ++i){
      rhs[(std::size_t)i].x = w_star * Sstar[(std::size_t)i].x + w_n * S0[(std::size_t)i].x + c_bdf * dt_half * g[(std::size_t)i].x;
      rhs[(std::size_t)i].y = w_star * Sstar[(std::size_t)i].y + w_n * S0[(std::size_t)i].y + c_bdf * dt_half * g[(std::size_t)i].y;
      rhs[(std::size_t)i].z = w_star * Sstar[(std::size_t)i].z + w_n * S0[(std::size_t)i].z + c_bdf * dt_half * g[(std::size_t)i].z;
   }
   assemble_diffusion_half_blocks(nf, c_bdf * dt_half, A_face, L, D, U);
   block3_thomas(L, D, U, rhs);

   for(int i = 0; i < nf; ++i){
      S[3 * i + 0] = rhs[(std::size_t)i].x;
      S[3 * i + 1] = rhs[(std::size_t)i].y;
      S[3 * i + 2] = rhs[(std::size_t)i].z;
   }
}

// A null edge array means that term is absent. A stored value is used as it is.
static double face_current(const double* a, const int e){
   if(a == nullptr) return 0.0;
   return a[e];
}

static double cell_factor(const double* a, const int i){
   if(a == nullptr) return 0.0;
   return a[i];
}

// P_charge matches P_field. P_spin is zero, and that term is still written.
// J_spin is inside the P_field product so its screening return is not polarized again.
// An open-circuit reduction
//   J_drift · m = (P_field - P_S) J_seebeck + (P_field - P_hot) J_hot
// is an analytic example only. It needs a vanishing five-term residual and one
// polarisation on the face. The residual is stored, and the face factors are
// averages of the two cells, so the drift keeps every term below.
static double drift_along_m(const double p_field,
                            const double p_charge,
                            const double p_seebeck,
                            const double p_hot,
                            const double j_field,
                            const double j_spin,
                            const double j_charge,
                            const double j_seebeck,
                            const double j_hot){
   const double p_spin = 0.0;
   return -p_field * (j_field + j_spin)
        - p_charge * j_charge
        - p_seebeck * j_seebeck
        - p_hot * j_hot
        - p_spin * j_spin;
}

static void build_drift_edges(const int nf,
                              const double* Bc,
                              const double* Bs,
                              const double* Ph,
                              const double* mhat,
                              const double* J_field,
                              const double* J_spin,
                              const double* J_charge,
                              const double* J_seebeck,
                              const double* J_hot,
                              double* Jd_edge){
   auto along = [&](const int e, const double p_field, const double p_see, const double p_hot){
      return drift_along_m(p_field, p_field, p_see, p_hot,
                           face_current(J_field, e), face_current(J_spin, e),
                           face_current(J_charge, e), face_current(J_seebeck, e),
                           face_current(J_hot, e));
   };
   const double j0 = along(0, Bc[0], cell_factor(Bs, 0), cell_factor(Ph, 0));
   Jd_edge[0] = j0*mhat[0];
   Jd_edge[1] = j0*mhat[1];
   Jd_edge[2] = j0*mhat[2];
   for(int i=0;i<nf-1;i++){
      const double pf = 0.5*(Bc[i] + Bc[i+1]);
      const double ps = 0.5*(cell_factor(Bs, i) + cell_factor(Bs, i+1));
      const double ph = 0.5*(cell_factor(Ph, i) + cell_factor(Ph, i+1));
      double hx = mhat[3*i+0] + mhat[3*(i+1)+0];
      double hy = mhat[3*i+1] + mhat[3*(i+1)+1];
      double hz = mhat[3*i+2] + mhat[3*(i+1)+2];
      normalize3(hx,hy,hz);
      const double j = along(i+1, pf, ps, ph);
      Jd_edge[3*(i+1)+0] = j*hx;
      Jd_edge[3*(i+1)+1] = j*hy;
      Jd_edge[3*(i+1)+2] = j*hz;
   }
   const int iN = nf - 1;
   const double jN = along(nf, Bc[iN], cell_factor(Bs, iN), cell_factor(Ph, iN));
   Jd_edge[3*nf+0] = jN*mhat[3*iN+0];
   Jd_edge[3*nf+1] = jN*mhat[3*iN+1];
   Jd_edge[3*nf+2] = jN*mhat[3*iN+2];
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

      const double denom = 1.0 - spin_gamma(Bc[i], Bd[i]);
      const double BBp = 1.0 / std::sqrt(denom);
      double lambda0 = lsf[i];
      double D0 = D[i];
      if(sc1d_thermal_effects){
         const double Te = sc1d_Te_fine[(std::size_t)i];
         const double Tp = sc1d_Tp_fine[(std::size_t)i];
         lambda0 = scale_lambda_sdl(lambda0, Te, Tp);
         D0 = scale_diffusion(D0, Te, Tp);
      }
      const double lambda_sf = lambda0 * BBp;
      const double Di = std::max(D0, 1e-30);
      const double inv_lambda_sf_sq = (lambda_sf > 0.0) ? (1.0 / (lambda_sf*lambda_sf)) : 0.0;
      const double inv_lambda_phi_sq = (lphi[i] > 0.0) ? (1.0 / (lphi[i]*lphi[i])) : 0.0;

      const double omega = Jsd[i] / (2.0*kHbar);

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
      // excess to sa_inf. chi is demag-spin-coupling in C/m^3. dmdt is
      // d(|m|/msat)/dt, with msat the cell sum of atomic moments in mu_B.
      const double sa_eq_amp = sa_inf[i];
      const double sa_eq_x = sa_eq_amp * hx;
      const double sa_eq_y = sa_eq_amp * hy;
      const double sa_eq_z = sa_eq_amp * hz;
      const double demag = (chi[i] == 0.0) ? 0.0 : -chi[i] * dmdt[i];
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
                                   const double* Jd_edge,
                                   const double Bc_i,
                                   const double Bd_i,
                                   const double D_i,
                                   const double* mhat,
                                   double& jsx,
                                   double& jsy,
                                   double& jsz){
   const double jdr_x = 0.5*(Jd_edge[3*i+0] + Jd_edge[3*(i+1)+0]);
   const double jdr_y = 0.5*(Jd_edge[3*i+1] + Jd_edge[3*(i+1)+1]);
   const double jdr_z = 0.5*(Jd_edge[3*i+2] + Jd_edge[3*(i+1)+2]);

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

   const double mx = mhat[3*i+0];
   const double my = mhat[3*i+1];
   const double mz = mhat[3*i+2];
   const double gamma = spin_gamma(Bc_i, Bd_i);
   const double parallel = mx*dSdz_x + my*dSdz_y + mz*dSdz_z;
   const double along = D_i * gamma * parallel;
   const double jdiff_x = -D_i * dSdz_x + along * mx;
   const double jdiff_y = -D_i * dSdz_y + along * my;
   const double jdiff_z = -D_i * dSdz_z + along * mz;

   jsx = jdr_x + jdiff_x;
   jsy = jdr_y + jdiff_y;
   jsz = jdr_z + jdiff_z;
}

static void coarsen_spin_and_map_torque(const int stack,
                                        const int start_cell,
                                        const int ncz,
                                        const int nsub,
                                        const int nf,
                                        const double dz,
                                        const double* Sbase,
                                        const double* ns_fine,
                                        const double* ne_fine,
                                        const double* Jc_edge_fine,
                                        const double* Jd_edge,
                                        const double* Bc,
                                        const double* Bd,
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
   sc1d_time_coarsen += sw.elapsed_seconds();

   sw.start();
   for(int k=0;k<ncz;k++){
      const int cell = start_cell + k;
      double jsx_sum=0.0, jsy_sum=0.0, jsz_sum=0.0;

      for(int j=0;j<nsub;j++){
         const int i = k*nsub + j;
         double jsx = 0.0;
         double jsy = 0.0;
         double jsz = 0.0;
         fine_cell_spin_current(Sbase, i, nf, dz, Jd_edge, Bc[i], Bd[i], D[i], mhat, jsx, jsy, jsz);
         const std::size_t slot = ((std::size_t)stack * (std::size_t)nf + (std::size_t)i) * 3u;
         sc1d_Js_fine[slot + 0u] = jsx;
         sc1d_Js_fine[slot + 1u] = jsy;
         sc1d_Js_fine[slot + 2u] = jsz;
         jsx_sum += jsx;
         jsy_sum += jsy;
         jsz_sum += jsz;
      }
      const double inv = 1.0 / (double)nsub;
      js_final[3*cell+0] = jsx_sum * inv;
      js_final[3*cell+1] = jsy_sum * inv;
      js_final[3*cell+2] = jsz_sum * inv;
   }
   sc1d_time_coarsen += sw.elapsed_seconds();
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
   if(sc1d_step_counter == 0UL){
      sc1d_time_total = 0.0;
      sc1d_time_io = 0.0;
      sc1d_time_charge = 0.0;
      sc1d_time_spin_current = 0.0;
      sc1d_time_spin_acc = 0.0;
      sc1d_time_diffusion = 0.0;
      sc1d_time_eom = 0.0;
      sc1d_time_interp = 0.0;
      sc1d_time_bcast = 0.0;
      sc1d_time_coarsen = 0.0;
      sc1d_time_other = 0.0;
   }
   stopwatch_t sw_total;
   sw_total.start();
   sc1d_step_timer_active = true;
   stopwatch_t sw;

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
   refresh_fine_temperature_profile();
   // Relax holds the atomic spins fixed, so the profile from initialise stays
   // valid. Equilibration and drive move m, and those steps prolong again.
   if(sc1d_step_counter > sc1d_relax_steps){
      sw.start();
      prolong_fine_magnetisation(dt);
      sc1d_time_interp += sw.elapsed_seconds();
   }

   const bool first_step = (sc1d_step_counter == 1UL) ||
                           (sc1d_step_counter == sc1d_relax_steps + 1UL);
   const bool want_bdf2 = (sc1d_diffusion == sc1d_diff_bdf2);

   for(int ls=0; ls<(int)sc1d_local_stacks.size(); ++ls){
      const int stack = sc1d_local_stacks[ls];
      const int start_cell = stack_index_y[stack];
      double* Sbase = sc1d_Sfine.data() + (size_t)stack * (size_t)nf * 3u;
      double* ns_fine = sc1d_ns_fine.data() + (size_t)stack * (size_t)nf;
      double* ne_fine = sc1d_ne_fine.data() + (size_t)stack * (size_t)nf;
      double* V_fine = sc1d_V_fine.data() + (size_t)stack * (size_t)nf;
      double* Jc_edge_fine = sc1d_Jc_edge_fine.data() + (size_t)stack * (size_t)(nf + 1);
      double* J_seebeck = sc1d_J_seebeck_edge.data() + (size_t)stack * (size_t)(nf + 1);
      double* J_restore = sc1d_J_restore_edge.data() + (size_t)stack * (size_t)(nf + 1);
      double* J_hot = sc1d_J_hot_edge.data() + (size_t)stack * (size_t)(nf + 1);
      double* J_cd = sc1d_J_cd_edge.data() + (size_t)stack * (size_t)(nf + 1);
      double* J_spin = sc1d_J_sd_edge.data() + (size_t)stack * (size_t)(nf + 1);
      const double* seebeck_fine = sc1d_seebeck_fine.data() + (size_t)stack * (size_t)nf;
      const double* Bs = sc1d_beta_seebeck_fine.data() + (size_t)stack * (size_t)nf;
      const double* Ph = sc1d_phot_fine.data() + (size_t)stack * (size_t)nf;

      double* Bc = sc1d_Bc_fine.data() + (size_t)stack * (size_t)nf;
      double* Bd = sc1d_Bd_fine.data() + (size_t)stack * (size_t)nf;
      double* D = sc1d_D_fine.data() + (size_t)stack * (size_t)nf;
      double* lsf = sc1d_lsf_fine.data() + (size_t)stack * (size_t)nf;
      double* lphi = sc1d_lphi_fine.data() + (size_t)stack * (size_t)nf;
      double* Jsd = sc1d_Jsd_fine.data() + (size_t)stack * (size_t)nf;
      const double* chi = sc1d_chi_fine.data() + (size_t)stack * (size_t)nf;
      const double* sa_inf = sc1d_sa_inf_fine.data() + (size_t)stack * (size_t)nf;
      double* sigma_fine = sc1d_sigma_fine.data() + (size_t)stack * (size_t)nf;
      double* k_prev = sc1d_k_prev.data() + (size_t)stack * (size_t)nf * 3u;
      const double* mhat = sc1d_mhat_fine.data() + (size_t)stack * (size_t)nf * 3u;
      const double* dmdt = sc1d_dmdt_fine.data() + (size_t)stack * (size_t)nf;
      const double* Rint = sc1d_Rint_face.data() + (size_t)stack * (size_t)nf;
      const double* sigma_bulk = sc1d_sigma_bulk_fine.data() + (size_t)stack * (size_t)nf;

      const bool do_charge_update = (sc1d_charge_stride <= 1) || ((sc1d_step_counter % (unsigned long)sc1d_charge_stride) == 0UL);
      if(do_charge_update){
         const double interp_before = sc1d_time_interp;
         sw.start();
         const double* Dmat = sc1d_Dmat_fine.data() + (size_t)stack * (size_t)nf;
         update_charge_transients_1d_stack(stack, ns_fine, ne_fine, V_fine, Jc_edge_fine,
                                           J_seebeck, J_restore, J_hot, J_cd, J_spin,
                                           nf, dz, dt, t_s, sigma_bulk,
                                           D, Dmat, sigma_fine, Bd, seebeck_fine, Sbase, mhat);
         account_charge_time(sw.elapsed_seconds(), interp_before);
      }

      sw.start();
      assemble_diffusion_faces(nf, dz, D, Bc, Bd, mhat, Rint, step_A_face);
      if(want_bdf2){
         apply_diffusion_tr_bdf2_half_step(Sbase, nf, dz, dt_half, step_A_face, Bc, mhat, Jc_edge_fine, step_L_blk, step_D_blk, step_U_blk);
      } else {
         assemble_diffusion_half_blocks(nf, dt_half, step_A_face, step_L_blk, step_D_blk, step_U_blk);
         apply_diffusion_half_step(Sbase, nf, dz, dt_half, Bc, mhat, Jc_edge_fine, step_L_blk, step_D_blk, step_U_blk);
      }
      sc1d_time_diffusion += sw.elapsed_seconds();
      sw.start();
      build_drift_edges(nf, Bc, Bs, Ph, mhat, J_restore, J_spin, J_cd, J_seebeck, J_hot, step_Jd_edge.data());
      {
         const double J0_tunnel = saved_J0[static_cast<std::size_t>(stack)];
         const std::vector<ChargeBarrier>& segs = saved_barriers[static_cast<std::size_t>(stack)];
         for(std::size_t b = 0; b < segs.size(); ++b){
            const ChargeBarrier& seg = segs[b];
            if(seg.face_lo > seg.face_hi) continue;
            for(int e = seg.face_lo; e <= seg.face_hi; ++e){
               step_Jd_edge[3 * e + 0] = seg.sfrac * J0_tunnel * seg.mx;
               step_Jd_edge[3 * e + 1] = seg.sfrac * J0_tunnel * seg.my;
               step_Jd_edge[3 * e + 2] = seg.sfrac * J0_tunnel * seg.mz;
            }
         }
      }
      sc1d_time_spin_current += sw.elapsed_seconds();
      sw.start();
      integrate_local_terms(Sbase, k_prev, nf, dt, first_step, dz,
                            Bc, Bd, D, lsf, lphi, Jsd, chi, sa_inf,
                            mhat, mhat, dmdt, step_Jd_edge.data());
      sc1d_time_eom += sw.elapsed_seconds();
      sw.start();
      if(want_bdf2){
         apply_diffusion_tr_bdf2_half_step(Sbase, nf, dz, dt_half, step_A_face, Bc, mhat, Jc_edge_fine, step_L_blk, step_D_blk, step_U_blk);
      } else {
         apply_diffusion_half_step(Sbase, nf, dz, dt_half, Bc, mhat, Jc_edge_fine, step_L_blk, step_D_blk, step_U_blk);
      }
      sc1d_time_diffusion += sw.elapsed_seconds();
      coarsen_spin_and_map_torque(stack, start_cell, ncz, nsub, nf, dz, Sbase,
                                  ns_fine, ne_fine, Jc_edge_fine, step_Jd_edge.data(),
                                  Bc, Bd, D, mhat);
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
}

} // namespace internal

double mean_joule_power(const double z_lo_angstrom, const double z_hi_angstrom){
   if(!internal::sc1d_enable) return 0.0;
   const int nf = internal::sc1d_nf;
   const int nstack = internal::num_stacks_y;
   if(nf <= 0 || nstack <= 0 || internal::sc1d_joule_fine.empty()) return 0.0;
   if(!(z_hi_angstrom > z_lo_angstrom)) return 0.0;
   const int nsub = std::max(1, internal::sc1d_nsub);
   const double dzA = internal::micro_cell_thickness / static_cast<double>(nsub);
   double sum = 0.0;
   int n = 0;
   for(int stack = 0; stack < nstack; ++stack){
      const double* joule = internal::sc1d_joule_fine.data() + static_cast<std::size_t>(stack) * static_cast<std::size_t>(nf);
      for(int i = 0; i < nf; ++i){
         const double zc = (static_cast<double>(i) + 0.5) * dzA;
         if(zc < z_lo_angstrom || zc >= z_hi_angstrom) continue;
         sum += joule[i];
         ++n;
      }
   }
   if(n == 0) return 0.0;
   return sum / static_cast<double>(n);
}

} // namespace st
