//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
// Laser-driven charge current from a charge-conserving two-filter ODE of the
// Serban laser envelope
//
//   G(t) = exp( -8 ((t - tau_pw) / tau_pw)^2 )
//   tau_s dS/dt + S = G
//   lambda dB/dt + B = G
//   Jc = (A0 / tau_pw) (S - B)
//
// A0 is not an input. The negative lobe carries a fraction eta_t of the
// excited sheet charge at a fixed optical wavelength of 800 nm:
//
//   Sigma = e * F * lambda / (h*c)
//   |integral_{Jc<0} Jc dt| = eta_t * Sigma
//
// F is sim:laser-pulse-power, the fluence in J/m^2. One electron per photon.
// Jc is in A/m^2. LET STT in fields.cpp: Tesla values × Jc/Jc_ref;
// torkance and Eq. 6 × Js = P (μB/e) Jc.
//
// Optional heating is independent of temperature-pulse (program 6). If LTMP is
// enabled, cell Te/Tp are advanced there. Otherwise this program integrates the
// global TTM (Ce ∝ Te) and writes sim::TTTe / sim::TTTp.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>

// Vampire headers
#include "errors.hpp"
#include "ltmp.hpp"
#include "material.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"

// program module headers
#include "internal.hpp"

namespace pg = program;
namespace pgi = program::internal;

namespace {

// Optical wavelength used to turn fluence into an excited sheet charge.
// Change this constant if the laser wavelength changes.
constexpr double laser_wavelength = 800.0e-9; // m
constexpr double planck = 6.62607015e-34;     // J s
constexpr double light_speed = 2.99792458e8;  // m/s
constexpr double electron_charge = 1.602176634e-19; // C

// ∫_{S<B} (B-S) dt for the unit filters. |Q_neg| = (A0/tau_pw) * this weight.
double negative_lobe_weight(const double tau_pw, const double tau_s, const double lambda, const double dt){

   if(!(tau_pw > 0.0) || !(dt > 0.0) || !(tau_s > 0.0) || !(lambda > 0.0)) return 0.0;

   const double t_end = std::max(40.0 * lambda, std::max(30.0 * tau_pw, 20.0e-12));
   const int64_t n = static_cast<int64_t>(std::ceil(t_end / dt));
   const double sd = std::exp(-dt / tau_s);
   const double bd = std::exp(-dt / lambda);
   double S = 0.0;
   double B = 0.0;
   double acc = 0.0;
   double prev = 0.0;

   for(int64_t i = 0; i < n; ++i){
      const double x = (static_cast<double>(i) * dt - tau_pw) / tau_pw;
      const double g = std::exp(-8.0 * x * x);
      S = g + (S - g) * sd;
      B = g + (B - g) * bd;
      const double sample = (B > S) ? (B - S) : 0.0;
      if(i > 0) acc += 0.5 * (prev + sample) * dt;
      prev = sample;
   }
   return acc;
}

} // end of anonymous namespace

namespace program{

void laser_electrical_pulse(){

   if(err::check==true){std::cout << "program::laser_electrical_pulse has been called" << std::endl;}

   const bool use_ltmp = (pgi::laser_electrical_enable_temperature && ltmp::is_enabled());
   const bool use_ttm  = (pgi::laser_electrical_enable_temperature && !use_ltmp);

   if(use_ltmp) ltmp::set_serban_laser_envelope(true);

   // Equilibrate at Teq with zero current
   const double temp = sim::temperature;
   sim::temperature = sim::Teq;
   sim::TTTe = sim::Teq;
   sim::TTTp = sim::Teq;
   pg::laser_electrical_current = 0.0;
   pg::laser_electrical_S = 0.0;
   pg::laser_electrical_B = 0.0;
   if(sim::local_temperature==true){
      for(unsigned int mat=0;mat<mp::material.size();mat++){
         if(mp::material[mat].couple_to_phonon_temperature==true) mp::material[mat].temperature=sim::TTTp;
         else mp::material[mat].temperature=sim::TTTe;
      }
   }

   sim::actual_H_field = sim::equilibrium_H_field;
   sim::actual_H_vector[0] = sim::equilibrium_H_vector[0];
   sim::actual_H_vector[1] = sim::equilibrium_H_vector[1];
   sim::actual_H_vector[2] = sim::equilibrium_H_vector[2];

   while(sim::time < sim::equilibration_time){

      for(uint64_t tt=0; tt < sim::partial_time; tt++){
         if(use_ltmp) ltmp::update_localised_temperature(0.0);
         sim::integrate(1);
      }

      stats::update();
      vout::data();
   }

   // Restore constant temperature if heating is disabled
   if(pgi::laser_electrical_enable_temperature != true){
      sim::temperature = temp;
      sim::TTTe = temp;
      sim::TTTp = temp;
   }
   else{
      sim::temperature = sim::Teq;
      sim::TTTe = sim::Teq;
      sim::TTTp = sim::Teq;
   }
   if(sim::local_temperature==true){
      for(unsigned int mat=0;mat<mp::material.size();mat++){
         if(mp::material[mat].couple_to_phonon_temperature==true) mp::material[mat].temperature=sim::TTTp;
         else mp::material[mat].temperature=sim::TTTe;
      }
   }

   const uint64_t start_time = sim::time;
   pg::laser_electrical_S = 0.0;
   pg::laser_electrical_B = 0.0;

   const double tau_pw = sim::pump_time;
   const double dt = mp::dt_SI;
   // eta_t * e * F * lambda / (h c), with F = laser-pulse-power in J/m^2
   const double sheet_charge = pgi::laser_electrical_eta_t * electron_charge
      * sim::pump_power * laser_wavelength / (planck * light_speed);
   const double lobe_weight = negative_lobe_weight(tau_pw, pgi::laser_electrical_tau_s,
                                                   pgi::laser_electrical_lambda, dt);
   const double A0_over_tau = (lobe_weight > 0.0) ? sheet_charge / lobe_weight : 0.0;
   const double A0 = A0_over_tau * tau_pw;
   zlog << zTs() << "Laser-electrical Jc: eta_t = " << pgi::laser_electrical_eta_t
        << ", fluence = " << sim::pump_power << " J/m^2, wavelength = 800 nm, A0 = "
        << A0 << " C/m^2" << std::endl;
   const double S_decay = std::exp(-dt / pgi::laser_electrical_tau_s);
   const double B_decay = std::exp(-dt / pgi::laser_electrical_lambda);
   const double Gep = sim::TTG;
   const double Ce = sim::TTCe;
   const double Cl = sim::TTCl;
   const double heat_sink_dt = sim::HeatSinkCouplingConstant * dt;
   const double Teq = sim::Teq;
   const double pump_prefactor = 93943727.87 * sim::pump_power / tau_pw;
   const int ncells = use_ltmp ? ltmp::get_num_cells() : 0;
   const int top = (ncells > 0) ? ncells - 1 : 0;

   sim::actual_H_field = sim::applied_H_field;
   sim::actual_H_vector[0] = sim::applied_H_vector[0];
   sim::actual_H_vector[1] = sim::applied_H_vector[1];
   sim::actual_H_vector[2] = sim::applied_H_vector[2];

   while(sim::time < sim::total_time+start_time){

      for(uint64_t tt=0; tt < sim::partial_time; tt++){

         const double time_from_start = dt * double(sim::time-start_time);

         // Serban envelope G(t), peak at t = tau_pw
         const double x = (time_from_start - tau_pw) / tau_pw;
         const double Genv = std::exp(-8.0 * x * x);

         // tau_s S' + S = G and lambda B' + B = G, exact step with G held over dt
         pg::laser_electrical_S = Genv + (pg::laser_electrical_S - Genv) * S_decay;
         pg::laser_electrical_B = Genv + (pg::laser_electrical_B - Genv) * B_decay;
         pg::laser_electrical_current = A0_over_tau * (pg::laser_electrical_S - pg::laser_electrical_B);

         if(sim::enable_laser_torque_fields) sim::laser_torque_strength = Genv;

         if(use_ltmp){
            ltmp::update_localised_temperature(time_from_start);
            sim::TTTe = ltmp::get_electron_temperature(top);
            sim::TTTp = ltmp::get_phonon_temperature(top);
            sim::temperature = sim::TTTe;
            if(sim::local_temperature==true){
               for(unsigned int mat=0;mat<mp::material.size();mat++){
                  if(mp::material[mat].couple_to_phonon_temperature==true) mp::material[mat].temperature=sim::TTTp;
                  else mp::material[mat].temperature=sim::TTTe;
               }
            }
         }
         else if(use_ttm){
            // Global TTM as in temperature_pulse.cpp (Ce ∝ Te for Te > 1 K)
            const double pump = pump_prefactor * Genv;
            const double Te = sim::TTTe;
            const double Tp = sim::TTTp;
            if(Te>1.0) sim::TTTe = (-Gep*(Te-Tp)+pump)*dt/(Ce*Te) + Te;
            else sim::TTTe =       (-Gep*(Te-Tp)+pump)*dt/Ce + Te;
            sim::TTTp =            ( Gep*(Te-Tp)     )*dt/Cl + Tp - (Tp-Teq)*heat_sink_dt;
            sim::temperature = sim::TTTe;
            if(sim::local_temperature==true){
               for(unsigned int mat=0;mat<mp::material.size();mat++){
                  if(mp::material[mat].couple_to_phonon_temperature==true) mp::material[mat].temperature=sim::TTTp;
                  else mp::material[mat].temperature=sim::TTTe;
               }
            }
         }

         sim::integrate(1);
      }

      stats::update();
      vout::data();
   }

   return;

}

} // end of namespace program
