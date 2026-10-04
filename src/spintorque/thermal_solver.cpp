//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

// Vampire headers
#include "internal.hpp"
#include "sim.hpp"
#include "errors.hpp"
#include "vio.hpp"

namespace st {
namespace internal {

static double gaussian_pulse_thermal(const double t, const double t0, const double fwhm, const double equilibration_offset_s) {
   if(fwhm <= 0.0) return 0.0;
   const double sigma = fwhm / (2.0 * std::sqrt(2.0 * std::log(2.0)));
   const double x = (t - t0 - equilibration_offset_s) / sigma;
   return std::exp(-0.5 * x * x);
}

static double harmonic_mean(const double a, const double b) {
   if(a + b <= 1e-20) return 0.0;
   return 2.0 * a * b / (a + b);
}

//-----------------------------------------------------------------------------
// Einstein model phonon heat capacity
//-----------------------------------------------------------------------------
inline double einstein_model_phonon_heat_capacity(double T_D_over_T, double phonon_heat_capacity) {
   if(T_D_over_T > 24.0) T_D_over_T = 24.0;
   return phonon_heat_capacity * 4.0 * M_PI * M_PI * M_PI * M_PI / (5.0 * T_D_over_T * T_D_over_T * T_D_over_T);
}

//-----------------------------------------------------------------------------
// Phonon temperature projector step (predictor-corrector for T-dependent Cp)
//-----------------------------------------------------------------------------
inline double phonon_temperature_projector_step(double Tp_init, double delta_Jp, int cell, double sink_deltaTp = 0.0) {
   const double einstein_temp = sc1d_T_Debye[cell];
   const double phonon_heat_capacity = sc1d_Cp[cell];
   if(einstein_temp == 0.0) return phonon_heat_capacity;
   
   double T_D_over_T = einstein_temp / Tp_init;
   
   double phonon_temp_dep;
   if(T_D_over_T >= 12.0) {
      phonon_temp_dep = einstein_model_phonon_heat_capacity(T_D_over_T, phonon_heat_capacity);
   } else {
      int idx = static_cast<int>(std::floor(T_D_over_T * 1000.0));
      idx = std::max(0, std::min(idx, 11999));
      phonon_temp_dep = phonon_heat_capacity * sc1d_debye_table[idx];
   }
   
   double projected_Tp = Tp_init + delta_Jp / phonon_temp_dep - sink_deltaTp;
   double projected_T_D_over_T = einstein_temp / projected_Tp;
   
   double corrected_phonon_temp_dep;
   if(projected_T_D_over_T >= 12.0) {
      corrected_phonon_temp_dep = einstein_model_phonon_heat_capacity(projected_T_D_over_T, phonon_heat_capacity);
   } else {
      int idx = static_cast<int>(std::floor(projected_T_D_over_T * 1000.0));
      idx = std::max(0, std::min(idx, 11999));
      corrected_phonon_temp_dep = phonon_heat_capacity * sc1d_debye_table[idx];
   }
   
   return (corrected_phonon_temp_dep + phonon_temp_dep) * 0.5;
}

static void fill_stack_laser_source(const int ncz,
                                    const double dz,
                                    const double time_s,
                                    const double dt_si,
                                    std::vector<double>& laser_source){
   for(int k = 0; k < ncz; ++k){
      laser_source[k] = 0.0;
   }
   if(!sc1d_laser_enable || sc1d_laser_Q0 <= 0.0) return;
   if(sim::time <= sim::equilibration_time) return;

   const double d = std::max(sc1d_optical_absorption_length, 1e-18);
   const double z_max = (double)ncz * micro_cell_thickness * 1e-10;
   const double equilibration_offset_s = 0.0;
   (void)dt_si;

   for(int k = 0; k < ncz; ++k){
      const double z_cell = (k + 0.5) * dz;
      const double distance_from_top = z_max - z_cell;
      const double env_t = gaussian_pulse_thermal(time_s, sc1d_laser_t0, sc1d_laser_fwhm, equilibration_offset_s);
      const double attenuation = std::exp(-distance_from_top / d);
      laser_source[k] = sc1d_laser_Q0 * attenuation * env_t;
   }
}

static void add_thermal_neighbor_flux(const int cell,
                                      const int ncell,
                                      const int base_idx,
                                      const double Te_old,
                                      const double Tp_old,
                                      const double dz_sq,
                                      double& dTe_diff,
                                      double& dTp_diff){
   const int nidx = base_idx + ncell;
   const double nTe = sc1d_sqrt_Te[nidx] * sc1d_sqrt_Te[nidx];
   const double nTp = sc1d_sqrt_Tp[nidx] * sc1d_sqrt_Tp[nidx];

   // κe = κ0 Te/Tp. Harmonic-mean interface κ is the conservative FD of ∇·(κ∇Te)
   // and already includes ∇κ·∇Te, so no extra (∇Te)² term.
   const double kappa_e_self = sc1d_kappa_e[cell] * (Te_old / Tp_old);
   const double kappa_e_neighbor = sc1d_kappa_e[ncell] * (nTe / nTp);
   const double kappa_e_harmonic = harmonic_mean(kappa_e_self, kappa_e_neighbor);
   const double kappa_p_harmonic = harmonic_mean(sc1d_kappa_p[cell], sc1d_kappa_p[ncell]);

   dTe_diff += (nTe - Te_old) * kappa_e_harmonic / dz_sq;
   dTp_diff += (nTp - Tp_old) * kappa_p_harmonic / dz_sq;
}

static void advance_two_temperature_cell(const int stack_idx,
                                         const int cell,
                                         const int idx,
                                         const double dt,
                                         const double laser_source,
                                         const double delta_Te_diff,
                                         const double delta_Tp_diff){
   const double Te_old = sc1d_sqrt_Te[idx] * sc1d_sqrt_Te[idx];
   const double Tp_old = sc1d_sqrt_Tp[idx] * sc1d_sqrt_Tp[idx];

   const double coupling_Te = sc1d_G[cell] * (Tp_old - Te_old);
   const double coupling_Tp = sc1d_G[cell] * (Te_old - Tp_old);

   const double dTe_coupling_source = (coupling_Te + laser_source) * dt / (sc1d_Ce[cell] * Te_old);
   const double Te_new = Te_old + delta_Te_diff + dTe_coupling_source;

   const double delta_Jp_coupling = coupling_Tp * dt;
   const double phonon_temp_dep_coupling = phonon_temperature_projector_step(Tp_old, delta_Jp_coupling, cell);
   const double dTp_coupling = delta_Jp_coupling / phonon_temp_dep_coupling;
   const double Tp_new = Tp_old + delta_Tp_diff + dTp_coupling;

   if(!std::isfinite(Te_new) || Te_new <= 0.0){
      terminaltextcolor(RED);
      std::cerr << "Error: Invalid electron temperature in thermal solver!" << std::endl;
      std::cerr << "  Stack: " << stack_idx << ", Cell: " << cell << std::endl;
      std::cerr << "  Te_new: " << Te_new << " K" << std::endl;
      std::cerr << "  Te_old: " << Te_old << " K" << std::endl;
      terminaltextcolor(WHITE);
      err::vexit();
   }
   if(!std::isfinite(Tp_new) || Tp_new <= 0.0){
      terminaltextcolor(RED);
      std::cerr << "Error: Invalid phonon temperature in thermal solver!" << std::endl;
      std::cerr << "  Stack: " << stack_idx << ", Cell: " << cell << std::endl;
      std::cerr << "  Tp_new: " << Tp_new << " K" << std::endl;
      std::cerr << "  Tp_old: " << Tp_old << " K" << std::endl;
      terminaltextcolor(WHITE);
      err::vexit();
   }

   sc1d_sqrt_Te[idx] = std::sqrt(Te_new);
   sc1d_sqrt_Tp[idx] = std::sqrt(Tp_new);
   sc1d_Te[idx] = Te_new;
   sc1d_Tp[idx] = Tp_new;
}

//-----------------------------------------------------------------------------
// Update thermal gradients for one stack
//-----------------------------------------------------------------------------
void update_thermal_gradients_stack(int stack_idx, double time_s, double dt_si) {
   if(!sc1d_thermal_gradients_initialised) return;

   const int ncz = num_microcells_per_stack;
   const double dz = micro_cell_thickness * 1e-10;
   const double dt = dt_si;
   const double equilibration_temperature = sim::Teq;
   const double Tcool = sim::HeatSinkCouplingConstant;
   const int base_idx = stack_idx * ncz;
   const double dz_sq = dz * dz;

   std::vector<double> laser_source(ncz, 0.0);
   fill_stack_laser_source(ncz, dz, time_s, dt_si, laser_source);

   std::vector<double> delta_Te_diff(ncz, 0.0);
   std::vector<double> delta_Tp_diff(ncz, 0.0);

   if(ncz > 1 && dz_sq > 1e-20){
      const int cell = 0;
      const int idx = base_idx + cell;
      const double Te_old = sc1d_sqrt_Te[idx] * sc1d_sqrt_Te[idx];
      const double Tp_old = sc1d_sqrt_Tp[idx] * sc1d_sqrt_Tp[idx];
      double dTe_diff = 0.0;
      double dTp_diff = 0.0;
      add_thermal_neighbor_flux(cell, 1, base_idx, Te_old, Tp_old, dz_sq,
                                dTe_diff, dTp_diff);
      delta_Te_diff[cell] = dTe_diff * dt / (sc1d_Ce[cell] * Te_old);
      const double delta_Jp = dTp_diff * dt;
      const double sink_deltaTp = (Tp_old - equilibration_temperature) * Tcool * dt;
      const double phonon_temp_dep = phonon_temperature_projector_step(Tp_old, delta_Jp, cell, sink_deltaTp);
      delta_Tp_diff[cell] = delta_Jp / phonon_temp_dep - sink_deltaTp;
   }

   for(int k = 1; k < ncz - 1; ++k){
      const int cell = k;
      const int idx = base_idx + cell;
      const double Te_old = sc1d_sqrt_Te[idx] * sc1d_sqrt_Te[idx];
      const double Tp_old = sc1d_sqrt_Tp[idx] * sc1d_sqrt_Tp[idx];
      double dTe_diff = 0.0;
      double dTp_diff = 0.0;
      add_thermal_neighbor_flux(cell, k - 1, base_idx, Te_old, Tp_old, dz_sq,
                                dTe_diff, dTp_diff);
      add_thermal_neighbor_flux(cell, k + 1, base_idx, Te_old, Tp_old, dz_sq,
                                dTe_diff, dTp_diff);
      delta_Te_diff[cell] = dTe_diff * dt / (sc1d_Ce[cell] * Te_old);
      const double delta_Jp = dTp_diff * dt;
      const double phonon_temp_dep = phonon_temperature_projector_step(Tp_old, delta_Jp, cell);
      delta_Tp_diff[cell] = delta_Jp / phonon_temp_dep;
   }

   if(ncz > 1 && dz_sq > 1e-20){
      const int cell = ncz - 1;
      const int idx = base_idx + cell;
      const double Te_old = sc1d_sqrt_Te[idx] * sc1d_sqrt_Te[idx];
      const double Tp_old = sc1d_sqrt_Tp[idx] * sc1d_sqrt_Tp[idx];
      double dTe_diff = 0.0;
      double dTp_diff = 0.0;
      add_thermal_neighbor_flux(cell, ncz - 2, base_idx, Te_old, Tp_old, dz_sq,
                                dTe_diff, dTp_diff);
      delta_Te_diff[cell] = dTe_diff * dt / (sc1d_Ce[cell] * Te_old);
      const double delta_Jp = dTp_diff * dt;
      const double phonon_temp_dep = phonon_temperature_projector_step(Tp_old, delta_Jp, cell);
      delta_Tp_diff[cell] = delta_Jp / phonon_temp_dep;
   }

   for(int k = 0; k < ncz; ++k){
      const int idx = base_idx + k;
      advance_two_temperature_cell(stack_idx, k, idx, dt, laser_source[k],
                                   delta_Te_diff[k], delta_Tp_diff[k]);
   }
}

} // end of namespace internal
} // end of namespace st
