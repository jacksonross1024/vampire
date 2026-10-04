//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>
#include <cmath>
#include <vector>

// Vampire headers
#include "ltmp.hpp"
#include "vmpi.hpp"
#include "sim.hpp"
#include "errors.hpp"
#include "vio.hpp"
// Local temperature pulse headers
#include "internal.hpp"

namespace ltmp{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Function to calculate the local temperature using the two temperature model
      //
      // Pump assumes uniform heating and penetration depth of 10 nm
      // (see main program in src/program/temperature_pulse.cpp for more info)
      //-----------------------------------------------------------------------------

   inline double einstein_model_phonon_heat_capacity(double T_D_over_T, double phonon_heat_capacity) {
      if(T_D_over_T > 24.0) T_D_over_T = 24.0;
      return phonon_heat_capacity*4.0*M_PI*M_PI*M_PI*M_PI/(5.0*T_D_over_T*T_D_over_T*T_D_over_T);
      }

   static void abort_invalid_ltmp_temperature(const char* quantity, int cell, double T_new, double T_old,
                                              double pump, double Ce, double Cp, double dt){
      terminaltextcolor(RED);
      std::cerr << "Error: " << quantity << " is non-physical in ltmp (not capped)." << std::endl;
      std::cerr << "  Cell: " << cell << std::endl;
      std::cerr << "  T_new: " << T_new << " K   T_old: " << T_old << " K" << std::endl;
      std::cerr << "  pump: " << pump << " W/m^3   dt: " << dt << " s" << std::endl;
      std::cerr << "  Ce: " << Ce << " J/K^2/m^3   Cp: " << Cp << " J/K/m^3" << std::endl;
      std::cerr << "  Simulation terminated." << std::endl;
      terminaltextcolor(WHITE);
      err::vexit();
   }

   static double phonon_heat_capacity_at_T(const double T, const int cell){
      const double Cp0 = ltmp::internal::phonon_heat_capacity[cell];
      const double TD = ltmp::internal::Einstein_temperature[cell];
      if(!(TD > 0.0)) return Cp0;
      if(err::check==true && (!std::isfinite(T) || !(T > 0.0))){
         abort_invalid_ltmp_temperature("phonon temperature (Debye Cp)", cell, T, T, 0.0, 0.0, Cp0, ltmp::internal::dt);
      }
      const double x = TD / T;
      if(x >= 12.0) return einstein_model_phonon_heat_capacity(x, Cp0);
      const int idx = static_cast<int>(x * 1000.0);
      return Cp0 * ltmp::internal::Debeye_phonon_constant[idx];
   }

   inline double phonon_temperature_projector_step(double Tp_init, double delta_Jp, int cell ){
      const double phonon_temp_dep = phonon_heat_capacity_at_T(Tp_init, cell);
      const double projected_Tp = Tp_init + delta_Jp/phonon_temp_dep;
      const double corrected = phonon_heat_capacity_at_T(projected_Tp, cell);
      return 0.5*(corrected + phonon_temp_dep);
   }

   inline double phonon_temperature_projector_step(double Tp_init, double delta_Jp, int cell, double sink_deltaTp ){
      const double phonon_temp_dep = phonon_heat_capacity_at_T(Tp_init, cell);
      const double projected_Tp = Tp_init + delta_Jp/phonon_temp_dep - sink_deltaTp;
      const double corrected = phonon_heat_capacity_at_T(projected_Tp, cell);
      return 0.5*(corrected + phonon_temp_dep);
   }

   // Harmonic mean: 2 k_i k_j / (k_i + k_j). Vanishes if either side is an insulator.
   inline double harmonic_mean_kappa(const double k1, const double k2){
      const double ksum = k1 + k2;
      return (ksum == 0.0) ? 0.0 : 2.0 * k1 * k2 / ksum;
   }

   static void accumulate_cell_diffusion(const unsigned int cell, const double Te_old, const double Tp_old,
                                         double& dTe_diff, double& dTp_diff){
      dTe_diff = 0.0;
      dTp_diff = 0.0;
      for(int id = ltmp::internal::cell_neighbour_start_index[cell]; id < ltmp::internal::cell_neighbour_end_index[cell]; ++id) {
         const int ncell = ltmp::internal::cell_neighbour_list[id];
         const double nTe = root_temperature_array[2*ncell+0]*root_temperature_array[2*ncell+0];
         const double nTp = root_temperature_array[2*ncell+1]*root_temperature_array[2*ncell+1];
         const double dx = ltmp::internal::cell_position_array[3*ncell+0]- ltmp::internal::cell_position_array[3*cell+0];
         const double dy = ltmp::internal::cell_position_array[3*ncell+1]- ltmp::internal::cell_position_array[3*cell+1];
         const double dz = ltmp::internal::cell_position_array[3*ncell+2]- ltmp::internal::cell_position_array[3*cell+2];
         const double dr = dx*dx + dy*dy + dz*dz;
         // κe = κ0 Te/Tp (Chen–Beraun). Conservative face flux
         // ∇·(κ∇Te) uses harmonic-mean interface κ, which already contains ∇κ·∇Te.
         const double kappa_e = harmonic_mean_kappa(
            ltmp::internal::electron_thermal_conductivity[ncell] * (nTe / nTp),
            ltmp::internal::electron_thermal_conductivity[cell] * (Te_old / Tp_old));
         const double kappa_p = harmonic_mean_kappa(ltmp::internal::phonon_thermal_conductivity[ncell],
                                                    ltmp::internal::phonon_thermal_conductivity[cell]);
         dTe_diff += (nTe - Te_old) * kappa_e / (dr * 1e-20);
         dTp_diff += (nTp - Tp_old) * kappa_p / (dr * 1e-20);
      }
   }

   // Te_old is T = (sqrt(T))^2 from root_temperature_array. Ce is gamma in Ce = gamma Te.
   inline double electron_temperature_step(const double Te_old, const double dU_e, const double Ce){
      if(Te_old>1.0) return Te_old + dU_e/(Ce*Te_old);
      else           return Te_old + dU_e/Ce;
   }

   static void apply_coupling_and_pump(const unsigned int cell, const double pump, const double dt){
      const double Te_old = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0];
      const double Tp_old = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1];
      const double Ce = ltmp::internal::electron_heat_capacity[cell];

      const double coupling_Te = ltmp::internal::electron_phonon_coupling_constant[cell]*(Tp_old - Te_old);
      const double coupling_Tp = ltmp::internal::electron_phonon_coupling_constant[cell]*(Te_old - Tp_old);
      const double source_Te = pump*attenuation_array[cell];
      const double Te_new = electron_temperature_step(Te_old, (coupling_Te + source_Te)*dt, Ce);

      const double delta_Jp = coupling_Tp*dt;
      const double phonon_temp_dep = phonon_temperature_projector_step(Tp_old, delta_Jp, static_cast<int>(cell));
      const double Tp_new = Tp_old + delta_Jp/phonon_temp_dep;

      if(err::check==true){
         if(!std::isfinite(Te_new) || Te_new < 0.0){
            abort_invalid_ltmp_temperature("electron temperature", static_cast<int>(cell), Te_new, Te_old, source_Te, Ce, phonon_temp_dep, dt);
         }
         if(!std::isfinite(Tp_new) || Tp_new < 0.0){
            abort_invalid_ltmp_temperature("phonon temperature", static_cast<int>(cell), Tp_new, Tp_old, pump, Ce, phonon_temp_dep, dt);
         }
      }

      root_temperature_array[2*cell+0] = std::sqrt(Te_new);
      root_temperature_array[2*cell+1] = std::sqrt(Tp_new);
   }

   // Face conductance g_ij = κ_ij / Δr²  (W/m³/K). Same harmonic-mean κ as explicit.
   static void compute_face_conductances(std::vector<double>& ge, std::vector<double>& gp){
      const unsigned int nface = ltmp::internal::cell_neighbour_list.size();
      ge.assign(nface, 0.0);
      gp.assign(nface, 0.0);
      const unsigned int ncells = ltmp::internal::attenuation_array.size();
      for(unsigned int cell=0; cell<ncells; ++cell){
         const double Te = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0];
         const double Tp = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1];
         for(int id = ltmp::internal::cell_neighbour_start_index[cell]; id < ltmp::internal::cell_neighbour_end_index[cell]; ++id){
            const int ncell = ltmp::internal::cell_neighbour_list[id];
            const double nTe = root_temperature_array[2*ncell+0]*root_temperature_array[2*ncell+0];
            const double nTp = root_temperature_array[2*ncell+1]*root_temperature_array[2*ncell+1];
            const double dx = ltmp::internal::cell_position_array[3*ncell+0]- ltmp::internal::cell_position_array[3*cell+0];
            const double dy = ltmp::internal::cell_position_array[3*ncell+1]- ltmp::internal::cell_position_array[3*cell+1];
            const double dz = ltmp::internal::cell_position_array[3*ncell+2]- ltmp::internal::cell_position_array[3*cell+2];
            const double dr = dx*dx + dy*dy + dz*dz;
            const double inv_dr2 = 1.0 / (dr * 1e-20);
            const double kappa_e = harmonic_mean_kappa(
               ltmp::internal::electron_thermal_conductivity[ncell] * (nTe / nTp),
               ltmp::internal::electron_thermal_conductivity[cell] * (Te / Tp));
            const double kappa_p = harmonic_mean_kappa(ltmp::internal::phonon_thermal_conductivity[ncell],
                                                       ltmp::internal::phonon_thermal_conductivity[cell]);
            ge[id] = kappa_e * inv_dr2;
            gp[id] = kappa_p * inv_dr2;
         }
      }
   }

   static bool cells_are_1d_chain(const unsigned int ncells){
      for(unsigned int cell=0; cell<ncells; ++cell){
         for(int id = ltmp::internal::cell_neighbour_start_index[cell]; id < ltmp::internal::cell_neighbour_end_index[cell]; ++id){
            const int ncell = ltmp::internal::cell_neighbour_list[id];
            if(ncell != static_cast<int>(cell)-1 && ncell != static_cast<int>(cell)+1) return false;
         }
      }
      return true;
   }

   // Thomas algorithm for tridiagonal a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i.
   // Overwrites d with the solution. a[0] and c[n-1] unused.
   static bool thomas_solve(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c,
                            std::vector<double>& d, const int n){
      if(n <= 0) return false;
      if(n == 1){
         if(b[0] == 0.0) return false;
         d[0] /= b[0];
         return true;
      }
      for(int i=1; i<n; ++i){
         if(b[i-1] == 0.0) return false;
         const double w = a[i] / b[i-1];
         b[i] -= w * c[i-1];
         d[i] -= w * d[i-1];
      }
      if(b[n-1] == 0.0) return false;
      d[n-1] /= b[n-1];
      for(int i=n-2; i>=0; --i){
         d[i] = (d[i] - c[i] * d[i+1]) / b[i];
      }
      return true;
   }

   static void implicit_matvec(const std::vector<double>& C, const std::vector<double>& gface,
                               const double dt, const std::vector<double>& x, std::vector<double>& Ax,
                               const unsigned int ncells){
      for(unsigned int cell=0; cell<ncells; ++cell){
         double diag = C[cell];
         double off = 0.0;
         for(int id = ltmp::internal::cell_neighbour_start_index[cell]; id < ltmp::internal::cell_neighbour_end_index[cell]; ++id){
            const double gdt = dt * gface[id];
            diag += gdt;
            off += gdt * x[ltmp::internal::cell_neighbour_list[id]];
         }
         Ax[cell] = diag * x[cell] - off;
      }
   }

   static void implicit_cg_solve(const std::vector<double>& C, const std::vector<double>& gface,
                                 const double dt, std::vector<double>& x, const std::vector<double>& rhs,
                                 const unsigned int ncells){
      std::vector<double> r(ncells), p(ncells), Ap(ncells);
      implicit_matvec(C, gface, dt, x, Ap, ncells);
      double rnorm = 0.0;
      double dnorm = 0.0;
      for(unsigned int i=0; i<ncells; ++i){
         r[i] = rhs[i] - Ap[i];
         p[i] = r[i];
         rnorm += r[i]*r[i];
         dnorm += rhs[i]*rhs[i];
      }
      const double tol2 = 1.0e-20 * (dnorm > 0.0 ? dnorm : 1.0);
      const int maxit = static_cast<int>(ncells) + 2;
      for(int it=0; it<maxit; ++it){
         if(rnorm <= tol2) break;
         implicit_matvec(C, gface, dt, p, Ap, ncells);
         double pAp = 0.0;
         for(unsigned int i=0; i<ncells; ++i) pAp += p[i]*Ap[i];
         if(pAp == 0.0) break;
         const double alpha = rnorm / pAp;
         double rnorm_new = 0.0;
         for(unsigned int i=0; i<ncells; ++i){
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
            rnorm_new += r[i]*r[i];
         }
         const double beta = rnorm_new / rnorm;
         for(unsigned int i=0; i<ncells; ++i) p[i] = r[i] + beta * p[i];
         rnorm = rnorm_new;
      }
   }

   // Backward Euler on ∇·(κ∇T) only. C and κ lagged at T* (after coupling/pump).
   // (C_i + dt Σ g_ij) T_i^{n+1} - dt Σ g_ij T_j^{n+1} = C_i T_i*
   static void implicit_diffuse_field(const std::vector<double>& C, const std::vector<double>& gface,
                                      const std::vector<double>& Tstar, std::vector<double>& Tnew,
                                      const double dt, const unsigned int ncells, const bool use_thomas){
      Tnew = Tstar;
      std::vector<double> rhs(ncells);
      for(unsigned int i=0; i<ncells; ++i) rhs[i] = C[i] * Tstar[i];
      if(use_thomas){
         std::vector<double> a(ncells, 0.0), b(ncells, 0.0), c(ncells, 0.0), d(ncells, 0.0);
         for(unsigned int i=0; i<ncells; ++i){
            b[i] = C[i];
            d[i] = rhs[i];
            for(int id = ltmp::internal::cell_neighbour_start_index[i]; id < ltmp::internal::cell_neighbour_end_index[i]; ++id){
               const int j = ltmp::internal::cell_neighbour_list[id];
               const double gdt = dt * gface[id];
               b[i] += gdt;
               if(j == static_cast<int>(i)-1) a[i] = -gdt;
               else if(j == static_cast<int>(i)+1) c[i] = -gdt;
            }
         }
         if(!thomas_solve(a, b, c, d, static_cast<int>(ncells))){
            implicit_cg_solve(C, gface, dt, Tnew, rhs, ncells);
            return;
         }
         Tnew.swap(d);
      }
      else{
         implicit_cg_solve(C, gface, dt, Tnew, rhs, ncells);
      }
   }

   static void apply_implicit_diffusion(const unsigned int ncells, const double dt, const double pump){
      std::vector<double> ge, gp;
      compute_face_conductances(ge, gp);

      std::vector<double> Te_star(ncells), Tp_star(ncells), Ce_vol(ncells), Cp_vol(ncells);
      for(unsigned int cell=0; cell<ncells; ++cell){
         Te_star[cell] = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0];
         Tp_star[cell] = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1];
         const double gamma = ltmp::internal::electron_heat_capacity[cell];
         Ce_vol[cell] = (Te_star[cell] > 1.0) ? gamma * Te_star[cell] : gamma;
         Cp_vol[cell] = phonon_heat_capacity_at_T(Tp_star[cell], static_cast<int>(cell));
      }

      const bool use_thomas = cells_are_1d_chain(ncells);
      std::vector<double> Te_new, Tp_new;
      implicit_diffuse_field(Ce_vol, ge, Te_star, Te_new, dt, ncells, use_thomas);
      implicit_diffuse_field(Cp_vol, gp, Tp_star, Tp_new, dt, ncells, use_thomas);

      const double heat_sink_bottom = substrate_cool_bottom ? 1.0 : 0.0;
      const double heat_sink_top = substrate_cool_bottom ? 0.0 : 1.0;

      for(unsigned int cell=0; cell<ncells; ++cell){
         double sink_deltaTp = 0.0;
         if(ncells == 1){
            const double sink = substrate_cool_bottom ? heat_sink_bottom : heat_sink_top;
            sink_deltaTp = sink*(Tp_star[cell]-equilibration_temperature)*Tcool*dt;
         } else if(cell == 0){
            sink_deltaTp = heat_sink_bottom*(Tp_star[cell]-equilibration_temperature)*Tcool*dt;
         } else if(cell == ncells-1){
            sink_deltaTp = heat_sink_top*(Tp_star[cell]-equilibration_temperature)*Tcool*dt;
         }
         const double Tp_final = Tp_new[cell] - sink_deltaTp;
         if(err::check==true){
            if(!std::isfinite(Te_new[cell]) || Te_new[cell] < 0.0){
               abort_invalid_ltmp_temperature("electron temperature (implicit diffusion)", static_cast<int>(cell),
                                              Te_new[cell], Te_star[cell], pump, Ce_vol[cell], Cp_vol[cell], dt);
            }
            if(!std::isfinite(Tp_final) || Tp_final < 0.0){
               abort_invalid_ltmp_temperature("phonon temperature (implicit diffusion)", static_cast<int>(cell),
                                              Tp_final, Tp_star[cell], pump, Ce_vol[cell], Cp_vol[cell], dt);
            }
         }
         root_temperature_array[2*cell+0] = std::sqrt(Te_new[cell]);
         root_temperature_array[2*cell+1] = std::sqrt(Tp_final);
      }
   }

      void calculate_local_temperature_pulse(const double time_from_start) {

         const double gaussian = ltmp::internal::temporal_laser_envelope(time_from_start);
         double pump = 0.0;

         if(!ltmp::internal::use_serban_laser_envelope && ltmp::internal::use_sc1d_laser_params) {
            pump = ltmp::internal::sc1d_laser_Q0 * gaussian;
         } else {
            const double i_pump_time = 1.0/ltmp::internal::pump_time;
            const double two_delta_sqrt_pi_ln_2 = 0.9394372787;
            pump = 1e10*ltmp::internal::pump_power*two_delta_sqrt_pi_ln_2*gaussian*i_pump_time/penetration_depth;
         }

         if(sim::enable_laser_torque_fields) sim::laser_torque_strength = gaussian;

         const double dt = ltmp::internal::dt;
         const unsigned int ncells = ltmp::internal::attenuation_array.size();

         for(unsigned int cell=0; cell < ncells; ++cell) {
            apply_coupling_and_pump(cell, pump, dt);
         }

         if(ltmp::internal::implicit_diffusion){
            apply_implicit_diffusion(ncells, dt, pump);
         }
         else{
         const double heat_sink_bottom = substrate_cool_bottom ? 1.0 : 0.0;
         const double heat_sink_top = substrate_cool_bottom ? 0.0 : 1.0;

         for(unsigned int cell=0; cell < ncells; ++cell) {
            const double Te_old = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0];
            const double Tp_old = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1];
            double dTe_diff = 0.0;
            double dTp_diff = 0.0;
            accumulate_cell_diffusion(cell, Te_old, Tp_old, dTe_diff, dTp_diff);
            delta_temperature_array[2*cell+0] = dTe_diff;
            delta_temperature_array[2*cell+1] = dTp_diff;
         }

         for(unsigned int cell=0; cell < ncells; ++cell) {
            const double Te_old = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0];
            const double Tp_old = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1];
            const double Ce = ltmp::internal::electron_heat_capacity[cell];
            const double dTe_diff = delta_temperature_array[2*cell+0];
            const double dTp_diff = delta_temperature_array[2*cell+1];

            double sink_deltaTp = 0.0;
            if(ncells == 1){
               const double sink = substrate_cool_bottom ? heat_sink_bottom : heat_sink_top;
               sink_deltaTp = sink*(Tp_old-equilibration_temperature)*Tcool*dt;
            } else if(cell == 0){
               sink_deltaTp = heat_sink_bottom*(Tp_old-equilibration_temperature)*Tcool*dt;
            } else if(cell == ncells-1){
               sink_deltaTp = heat_sink_top*(Tp_old-equilibration_temperature)*Tcool*dt;
            }

            const double Te_new = electron_temperature_step(Te_old, dTe_diff*dt, Ce);
            const double delta_Jp = dTp_diff*dt;
            const double phonon_temp_dep = phonon_temperature_projector_step(Tp_old, delta_Jp, static_cast<int>(cell), sink_deltaTp);
            const double Tp_new = Tp_old + delta_Jp/phonon_temp_dep - sink_deltaTp;

            if(err::check==true){
               if(!std::isfinite(Te_new) || Te_new < 0.0){
                  abort_invalid_ltmp_temperature("electron temperature (diffusion)", static_cast<int>(cell),
                                                 Te_new, Te_old, pump, Ce, phonon_temp_dep, dt);
               }
               if(!std::isfinite(Tp_new) || Tp_new < 0.0){
                  abort_invalid_ltmp_temperature("phonon temperature (diffusion)", static_cast<int>(cell),
                                                 Tp_new, Tp_old, pump, Ce, phonon_temp_dep, dt);
               }
            }

            root_temperature_array[2*cell+0] = std::sqrt(Te_new);
            root_temperature_array[2*cell+1] = std::sqrt(Tp_new);
         }
         }

         if(ltmp::internal::output_microcell_data && (sim::time % 100 == 0)){
            ltmp::internal::write_cell_temperature_data();
         }

         return;

      }

   } // end of namespace internal
} // end of namespace ltmp
