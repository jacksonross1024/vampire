//------------------------------------------------------------------------------
// spincurrents.cpp
//
// Minimal 1D (along st::internal::stz) transient spin-accumulation
// solver, with:
//  - implicit (backward-Euler) anisotropic diffusion on a fine z-grid
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

// Hold the diffusion, drift and relaxation axis when |m| is inside the thermal
// floor 3 |m0| / sqrt(2 N). Compile-time flag, no user interface.
static constexpr bool   kFreezeAxisWhenMsmall = true;
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

// gamma = beta_c * beta_d, clamped so the parallel resistance stays finite.
// A non-magnet has beta ~ 0.001, so gamma is ~1e-6 and the tensor is isotropic.
static double spin_gamma(const double Bc, const double Bd){
   const double gmax = 1.0 - 1.0e-6;
   double g = Bc * Bd;
   if(g < 0.0) g = 0.0;
   if(g > gmax) g = gmax;
   return g;
}

// Fraction t in [0, 1] between coarse centres k and k+1. Outside the first
// and last centres the end cell is held and t is 0. Step jumps at the coarse
// face (t = 1/2). Smoothstep is flat at each centre.
static double prolongation_t(const double s, const int ncz, int& k_left){
   k_left = 0;
   if(ncz <= 1) return 0.0;
   if(s <= 0.5) return 0.0;
   if(s >= (double)ncz - 0.5){
      k_left = ncz - 1;
      return 0.0;
   }
   k_left = static_cast<int>(std::floor(s - 0.5));
   if(k_left < 0) k_left = 0;
   if(k_left > ncz - 2) k_left = ncz - 2;
   double t = s - ((double)k_left + 0.5);
   if(t < 0.0) t = 0.0;
   if(t > 1.0) t = 1.0;
   if(sc1d_prolongation == sc1d_prolong_step){
      if(t < 0.5) t = 0.0;
      else t = 1.0;
   } else if(sc1d_prolongation == sc1d_prolong_smooth){
      t = t*t*(3.0 - 2.0*t);
   }
   return t;
}

static double prolong_scalar(const double* field,
                             const int start_cell,
                             const int ncz,
                             const double s){
   int k = 0;
   const double t = prolongation_t(s, ncz, k);
   const double left = field[start_cell + k];
   if(t == 0.0 || k >= ncz - 1) return left;
   return lerp(left, field[start_cell + k + 1], t);
}

// True when a finite floor exists. With the freeze flag off, every finite
// moment is accepted and thresh is 0. N < 1 or |m0| == 0 keeps the reference.
static bool axis_floor(const int cell, double& thresh){
   thresh = 0.0;
   if(!kFreezeAxisWhenMsmall) return true;
   if(cell < 0 || cell >= (int)cell_natom.size()) return false;
   if(cell >= (int)sc1d_m0_mag.size()) return false;
   const double N = cell_natom[cell];
   const double m0 = sc1d_m0_mag[cell];
   if(N < 1.0 || !(m0 > 0.0)) return false;
   thresh = 3.0 * m0 / std::sqrt(2.0 * N);
   return true;
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
   std::vector<Block3> Cp(n, block3_zero());
   std::vector<Vec3> dp(n);

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
      // sigma0_cell is already a conductivity. The coarse cell resolved the
      // Einstein sentinel, and that number was prolonged with D and m.
      sigma_cell[i] = apply_thermal_effects(sigma0_cell[i], Te, Tp);
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
      sc1d_m0_mag[c] = mm;
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

// Resolve the coarse conductivity, then prolong that number with the same
// weight as D and m. A negative sentinel is the Einstein value at the coarse
// centre. Zero stays an insulator. The face harmonic mean is unchanged.
static void prolong_conductivity(const int start_cell,
                                 const int nsub,
                                 const int nf,
                                 double* sigma_fine){
   const int nsub_safe = std::max(1, nsub);
   const int ncz = std::max(1, num_microcells_per_stack);
   const double Dz = micro_cell_thickness * 1.0e-10;
   std::vector<double> sig((std::size_t)ncz, 0.0);
   for(int k=0;k<ncz;k++){
      const int cell = start_cell + k;
      double stored = 0.0;
      double Dc = 0.0;
      if(cell >= 0 && cell < (int)conductivity.size()) stored = conductivity[cell];
      if(cell >= 0 && cell < (int)diffusion.size()) Dc = diffusion[cell];
      if(stored < 0.0){
         const double zc = ((double)k + 0.5) * Dz;
         const double Te = get_local_electron_temperature(zc);
         sig[(std::size_t)k] = einstein_conductivity(std::max(0.0, Dc), Te);
      } else {
         sig[(std::size_t)k] = stored;
      }
   }
   for(int i=0;i<nf;i++){
      const double s = ((double)i + 0.5) / (double)nsub_safe;
      sigma_fine[i] = prolong_scalar(sig.data(), 0, ncz, s);
   }
}

// Fine-cell constants use the selected prolongation between coarse centres.
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
      Bc_fine[i] = prolong_scalar(beta_cond.data(), start_cell, ncz, s);
      Bd_fine[i] = prolong_scalar(beta_diff.data(), start_cell, ncz, s);
      D_fine[i] = prolong_scalar(diffusion.data(), start_cell, ncz, s);
      lsf_fine[i] = prolong_scalar(lambda_sdl.data(), start_cell, ncz, s);
      lphi_fine[i] = prolong_scalar(lambda_phi.data(), start_cell, ncz, s);
      Jsd_fine[i] = prolong_scalar(sd_exchange.data(), start_cell, ncz, s);
      chi_fine[i] = prolong_scalar(chi_demag.data(), start_cell, ncz, s);
      sa_inf_fine[i] = prolong_scalar(sa_infinity.data(), start_cell, ncz, s);
      seebeck_fine[i] = prolong_scalar(seebeck_coefficient.data(), start_cell, ncz, s);
   }
   prolong_conductivity(start_cell, nsub, nf, sigma_fine);
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
   const double gfac = (gamma > 0.0) ? (gamma / (1.0 - gamma)) : 0.0;
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
// zero and is not inverted. Interface resistance is R_int on the coarse face only.
static void assemble_diffusion_faces(const int start_cell,
                                    const int nsub,
                                    const int nf,
                                    const double dz,
                                    const double* D_fine,
                                    const double* Bc_fine,
                                    const double* Bd_fine,
                                    const double* mhat,
                                    std::vector<Block3>& A_face){
   if(nf < 2) return;
   const int nsub_safe = std::max(1, nsub);
   const double idz = 1.0 / dz;
   for(int i=0;i<nf-1;i++){
      const double Di = D_fine[i];
      const double Dj = D_fine[i+1];
      if(Di < kEps || Dj < kEps){
         A_face[(std::size_t)i] = block3_zero();
         continue;
      }
      double Rint = 0.0;
      if(((i + 1) % nsub_safe) == 0){
         const int cell = start_cell + (i / nsub_safe);
         if(cell >= 0 && cell < (int)r_int_edge.size()) Rint = r_int_edge[cell];
      }
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

static void initialise_fine_spin_if_empty(const int start_cell,
                                         const int nsub,
                                         const int nf,
                                         const double dz,
                                         const double* Bc_fine,
                                         const double* Bd_fine,
                                         const double* D_fine,
                                         const double* lsf_fine,
                                         const double* sa_inf_fine,
                                         const double* mhat,
                                         double* Sbase){
   bool needs_init = true;
   for(int i=0; i<nf*3; ++i){
      if(std::abs(Sbase[i]) > 1e-20){
         needs_init = false;
         break;
      }
   }
   if(!needs_init || nf <= 0) return;

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
   assemble_diffusion_faces(start_cell, nsub, nf, dz, D_fine, Bc_fine, Bd_fine, mhat, A_face);

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

static void fill_fine_time_dependent_fields(const int start_cell,
                                            const int ncz,
                                            const int nsub,
                                            const double* dm_dt_coarse,
                                            double* mhat,
                                            double* msrc,
                                            double* dmdt);

void initialise_spincurrents_1d(){
   const int nsub = compute_fine_subdivisions();
   sc1d_nsub = nsub;
   sc1d_nf   = nsub * num_microcells_per_stack;

   collect_local_stacks();
   allocate_sc1d_arrays();
   initialise_magnetization_history();

   const int nf = sc1d_nf;
   const int ncz = num_microcells_per_stack;
   const double dz = (micro_cell_thickness * 1.0e-10) / (double)nsub;
   const double dzA = micro_cell_thickness / (double)nsub;

   // Constants and m use the same prolongation between coarse centres.
   // Only the owner integrates that column.
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

      prolong_coarse_constants(start_cell, nsub, nf, Bc_fine, Bd_fine, D_fine, lsf_fine,
                               lphi_fine, Jsd_fine, chi_fine, sa_inf_fine, sigma_fine, seebeck_fine);
      if(!rank_owns_column(stack)) continue;

      std::vector<double> dm_coarse((std::size_t)ncz, 0.0);
      std::vector<double> dmdt_fine((std::size_t)nf, 0.0);
      std::vector<double> mhat_init(3*(std::size_t)nf, 0.0);
      std::vector<double> msrc_init(3*(std::size_t)nf, 0.0);
      fill_fine_time_dependent_fields(start_cell, ncz, nsub, dm_coarse.data(),
                                      mhat_init.data(), msrc_init.data(), dmdt_fine.data());

      double* Sbase = sc1d_Sfine.data() + (std::size_t)stack * (std::size_t)nf * 3u;
      initialise_fine_spin_if_empty(start_cell, nsub, nf, dz, Bc_fine, Bd_fine, D_fine, lsf_fine,
                                    sa_inf_fine, mhat_init.data(), Sbase);
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

      double thresh = 0.0;
      const bool have_floor = axis_floor(cell, thresh);
      if(have_floor && mm >= thresh && mm > kEps){
         sc1d_mhat_ref[3*cell+0] = mx/mm;
         sc1d_mhat_ref[3*cell+1] = my/mm;
         sc1d_mhat_ref[3*cell+2] = mz/mm;
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

// Above the floor the coarse vector is the live moment, length included, so an
// antiparallel join passes through zero. Below the floor it is the reference
// axis at length |m0|: m is a sum over atoms, and a collapsed length would
// let the neighbour set the geometric weight.
static void coarse_axis_vector(const int cell, double& vx, double& vy, double& vz){
   const double mx = m[3*cell+0];
   const double my = m[3*cell+1];
   const double mz = m[3*cell+2];
   const double mm = norm3(mx, my, mz);
   double thresh = 0.0;
   const bool ok = axis_floor(cell, thresh);
   if(ok && mm >= thresh && mm > kEps){
      vx = mx;
      vy = my;
      vz = mz;
      return;
   }
   double scale = 1.0;
   if(cell >= 0 && cell < (int)sc1d_m0_mag.size() && sc1d_m0_mag[cell] > 0.0){
      scale = sc1d_m0_mag[cell];
   }
   vx = sc1d_mhat_ref[3*cell+0] * scale;
   vy = sc1d_mhat_ref[3*cell+1] * scale;
   vz = sc1d_mhat_ref[3*cell+2] * scale;
}

static void fill_fine_time_dependent_fields(const int start_cell,
                                            const int ncz,
                                            const int nsub,
                                            const double* dm_dt_coarse,
                                            double* mhat,
                                            double* msrc,
                                            double* dmdt){
   const int nsub_safe = std::max(1, nsub);
   const int nf = ncz * nsub_safe;
   for(int i=0;i<nf;i++){
      const double s = ((double)i + 0.5) / (double)nsub_safe;
      int k = 0;
      const double t = prolongation_t(s, ncz, k);
      const int k2 = (t > 0.0 && k < ncz - 1) ? (k + 1) : k;
      const double w = (k2 == k) ? 0.0 : t;
      const int c0 = start_cell + k;
      const int c1 = start_cell + k2;

      double v0x = 0.0, v0y = 0.0, v0z = 0.0;
      double v1x = 0.0, v1y = 0.0, v1z = 0.0;
      coarse_axis_vector(c0, v0x, v0y, v0z);
      if(c1 == c0){
         v1x = v0x; v1y = v0y; v1z = v0z;
      } else {
         coarse_axis_vector(c1, v1x, v1y, v1z);
      }
      const double px = lerp(v0x, v1x, w);
      const double py = lerp(v0y, v1y, w);
      const double pz = lerp(v0z, v1z, w);
      const double hm = norm3(px, py, pz);

      const double hx = lerp(sc1d_mhat_ref[3*c0+0], sc1d_mhat_ref[3*c1+0], w);
      const double hy = lerp(sc1d_mhat_ref[3*c0+1], sc1d_mhat_ref[3*c1+1], w);
      const double hz = lerp(sc1d_mhat_ref[3*c0+2], sc1d_mhat_ref[3*c1+2], w);

      double f0 = 0.0;
      double f1 = 0.0;
      const bool ok0 = axis_floor(c0, f0);
      const bool ok1 = axis_floor(c1, f1);
      bool above = false;
      if(!kFreezeAxisWhenMsmall){
         above = true;
      } else if(ok0 && ok1){
         above = hm >= lerp(f0, f1, w);
      } else if(ok0){
         above = (w < 1.0) && (hm >= f0);
      } else if(ok1){
         above = (w > 0.0) && (hm >= f1);
      }

      double ux = sc1d_mhat_ref[3*c0+0];
      double uy = sc1d_mhat_ref[3*c0+1];
      double uz = sc1d_mhat_ref[3*c0+2];
      if(above && hm > kEps){
         ux = px / hm;
         uy = py / hm;
         uz = pz / hm;
      } else {
         const double hr = norm3(hx, hy, hz);
         if(hr > kEps){
            ux = hx / hr;
            uy = hy / hr;
            uz = hz / hr;
         }
      }

      mhat[3*i+0] = ux;
      mhat[3*i+1] = uy;
      mhat[3*i+2] = uz;
      msrc[3*i+0] = ux;
      msrc[3*i+1] = uy;
      msrc[3*i+2] = uz;
      dmdt[i] = lerp(dm_dt_coarse[k], dm_dt_coarse[k2], w);
   }
}

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
                                   const double Bd_i,
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

   const double mx = mhat[3*i+0];
   const double my = mhat[3*i+1];
   const double mz = mhat[3*i+2];
   const double gamma = spin_gamma(Bc_i, Bd_i);
   const double parallel = mx*dSdz_x + my*dSdz_y + mz*dSdz_z;
   const double along = D_i * gamma * parallel;
   double jdiff_x = -D_i * dSdz_x + along * mx;
   double jdiff_y = -D_i * dSdz_y + along * my;
   double jdiff_z = -D_i * dSdz_z + along * mz;
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
         fine_cell_spin_current(Sbase, i, nf, dz, Jc_loc, Bc[i], Bd[i], D[i], mhat, jsx, jsy, jsz);
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
   const int nface = (nf > 1) ? (nf - 1) : 0;
   std::vector<Block3> A_face((std::size_t)nface, block3_zero());
   std::vector<Block3> L_blk((std::size_t)std::max(0, nf), block3_zero());
   std::vector<Block3> D_blk((std::size_t)std::max(0, nf), block3_zero());
   std::vector<Block3> U_blk((std::size_t)std::max(0, nf), block3_zero());
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
      double* sigma_fine = sc1d_sigma_fine.data() + (size_t)stack * (size_t)nf;
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
         prolong_conductivity(start_cell, nsub, nf, sigma_fine);
         double* ne_seebeck = sc1d_ne_seebeck_fine.data() + (size_t)stack * (size_t)nf;
         update_charge_transients_1d_stack(ns_fine, ne_fine, V_fine, Jc_edge_fine, ne_seebeck,
                                           nf, dz, dt, t_s,
                                           D, sigma_fine, Bd, seebeck_fine, Sbase, mhat.data());
         account_charge_time(sw.elapsed_seconds(), interp_before);
      }

      sw.start();
      assemble_diffusion_faces(start_cell, nsub, nf, dz, D, Bc, Bd, mhat.data(), A_face);
      assemble_diffusion_half_blocks(nf, dt_half, A_face, L_blk, D_blk, U_blk);
      apply_diffusion_half_step(Sbase, nf, dz, dt_half, Bc, mhat.data(), Jc_edge_fine, L_blk, D_blk, U_blk);
      sc1d_time_spin_acc += sw.elapsed_seconds();
      sw.start();
      build_drift_edges(nf, Bc, mhat.data(), Jc_edge_fine, Jd_edge.data());
      sc1d_time_spin_current += sw.elapsed_seconds();
      sw.start();
      integrate_local_terms(Sbase, k_prev, nf, dt, first_step, dz,
                            Bc, Bd, D, lsf, lphi, Jsd, chi, sa_inf,
                            mhat.data(), msrc.data(), dmdt.data(), Jd_edge.data());
      apply_diffusion_half_step(Sbase, nf, dz, dt_half, Bc, mhat.data(), Jc_edge_fine, L_blk, D_blk, U_blk);
      sc1d_time_spin_acc += sw.elapsed_seconds();
      coarsen_spin_and_map_torque(start_cell, ncz, nsub, nf, dz, Sbase,
                                  ns_fine, ne_fine, Jc_edge_fine,
                                  Bc, Bd, D, mhat.data());
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
