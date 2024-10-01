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

// Vampire headers
#include "ltmp.hpp"
#include "vmpi.hpp"
#include "sim.hpp"
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
      void calculate_local_temperature_pulse(const double time_from_start){

         const double i_pump_time = 1.0/ltmp::internal::pump_time;
         const double reduced_time = (time_from_start - 2.0*ltmp::internal::pump_time)*i_pump_time;
         const double four_ln_2 = 2.77258872224; // 4 ln 2
         // 2/(delta sqrt(pi/ln 2))*0.1, delta = 10 nm, J/m^2 -> mJ/cm^2 (factor 0.1)
         const double two_delta_sqrt_pi_ln_2 = 0.9394372787;
         const double gaussian = exp(-four_ln_2*reduced_time*reduced_time);
         const double pump= 1e10*ltmp::internal::pump_power*two_delta_sqrt_pi_ln_2*gaussian*i_pump_time/penetration_depth;

         // if(sim::enable_laser_torque_fields) sim::laser_torque_strength = gaussian;

         const double G  = ltmp::internal::TTG;
         const double Ce = ltmp::internal::TTCe;
         const double Cl = ltmp::internal::TTCl;
         const double dt = ltmp::internal::dt;

         // Precalculate heat transfer constant k*L/V (J/K/m^3/s) (divide by Angstroms^2)
         const double dTe_diff_prefactor = 5*ltmp::internal::electron_thermal_conductivity/(2*ltmp::internal::micro_cell_size[0]*ltmp::internal::micro_cell_size[1]*1.e-20);
         const double dTp_diff_prefactor = ltmp::internal::phonon_thermal_conductivity/(2*ltmp::internal::micro_cell_size[0]*ltmp::internal::micro_cell_size[1]*1.e-20);
         // std::cout << dTe_diff_prefactor << ", " << dTp_diff_prefactor << std::endl;
         // Determine change in Te and Tp
         
         double Te = root_temperature_array[2*0+0]*root_temperature_array[2*0+0];
         double Tp = root_temperature_array[2*0+1]*root_temperature_array[2*0+1];

         // calculate heat transfer from neighbouring cells
         double dTe_diff = 0.0;
         double dTp_diff = 0.0;
         for(int id = ltmp::internal::cell_neighbour_start_index[0]; id < ltmp::internal::cell_neighbour_end_index[0]; ++id){
            const int ncell = ltmp::internal::cell_neighbour_list[id]; // neighbour cell id
            double nTe = root_temperature_array[2*ncell+0]*root_temperature_array[2*ncell+0];
            double nTp = root_temperature_array[2*ncell+1]*root_temperature_array[2*ncell+1];
            dTe_diff += nTe - Te;
            dTp_diff += nTp - Tp;
            // std::cout <<"first " << Te << ", " << nTe << ", " << dTe_diff << std::endl;
         }

         if(Te < 1.0) delta_temperature_array[2*0+0] = (G*(Tp-Te) + pump*attenuation_array[0] + 2*dTe_diff*dTe_diff_prefactor)*dt/(Ce);
         else delta_temperature_array[2*0+0] = (G*(Tp-Te) + pump*attenuation_array[0] + 2*dTe_diff*dTe_diff_prefactor)*dt/(Ce*Te);
         delta_temperature_array[2*0+1] = (G*(Te-Tp)                                + 2*dTp_diff*dTp_diff_prefactor)*dt/Cl;
         
         for(unsigned int cell=1; cell<ltmp::internal::attenuation_array.size()-1; ++cell) {

            Te = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0];
            Tp = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1];

            // calculate heat transfer from neighbouring cells
            dTe_diff = 0.0;
            dTp_diff = 0.0;
            // double ddTe1 = 0.0;
            
            int id = ltmp::internal::cell_neighbour_start_index[cell];
            int ncell = ltmp::internal::cell_neighbour_list[id];
            double ddTe1 = root_temperature_array[2*ncell+0]*root_temperature_array[2*ncell+0];
            double ddTp1 = root_temperature_array[2*ncell+1]*root_temperature_array[2*ncell+1];

            for( id; id < ltmp::internal::cell_neighbour_end_index[cell]; ++id){
               ncell = ltmp::internal::cell_neighbour_list[id]; // neighbour cell id
               double nTe = root_temperature_array[2*ncell+0]*root_temperature_array[2*ncell+0];
               double nTp = root_temperature_array[2*ncell+1]*root_temperature_array[2*ncell+1];
               dTe_diff += nTe - Te;
               dTp_diff += nTp - Tp;
               // std::cout << cell << ", " <<  Te << ", " << nTe << ", " << dTe_diff << std::endl;
            }
            double ddTe2 = root_temperature_array[2*ncell+0]*root_temperature_array[2*ncell+0];
            double ddTp2 = root_temperature_array[2*ncell+1]*root_temperature_array[2*ncell+1];

            if(Te < 1.0) delta_temperature_array[2*cell+0] = (G*(Tp-Te) + pump*attenuation_array[cell] + (dTe_diff+0.5*(ddTe2-ddTe1))*dTe_diff_prefactor)*dt/(Ce);
            else delta_temperature_array[2*cell+0] = (G*(Tp-Te) + pump*attenuation_array[cell] + (dTe_diff+0.5*(ddTe2-ddTe1))*dTe_diff_prefactor)*dt/(Ce*Te);
            delta_temperature_array[2*cell+1] = (G*(Te-Tp)                                + (dTp_diff+0.5*(ddTp2-ddTp1))*dTp_diff_prefactor)*dt/Cl;

         } // end of cell loop

         int cell = ltmp::internal::attenuation_array.size()-1;
         Te = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0];
         Tp = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1];

         // calculate heat transfer from neighbouring cells
         dTe_diff = 0.0;
         dTp_diff = 0.0;
         for(int id=ltmp::internal::cell_neighbour_start_index[cell]; id<ltmp::internal::cell_neighbour_end_index[cell]; ++id){
            const int ncell = ltmp::internal::cell_neighbour_list[id]; // neighbour cell id
            double nTe = root_temperature_array[2*ncell+0]*root_temperature_array[2*ncell+0];
            double nTp = root_temperature_array[2*ncell+1]*root_temperature_array[2*ncell+1];
            dTe_diff += nTe - Te;
            dTp_diff += nTp - Tp;
            
         }
         if(Te < 1.0) delta_temperature_array[2*cell+0] = (G*(Tp-Te) + pump*attenuation_array[cell] + 2*dTe_diff*dTe_diff_prefactor)*dt/(Ce);
         else delta_temperature_array[2*cell+0] = (G*(Tp-Te) + pump*attenuation_array[cell] + 2*dTe_diff*dTe_diff_prefactor)*dt/(Ce*Te);
         delta_temperature_array[2*cell+1] = (G*(Te-Tp)                                + 2*dTp_diff*dTp_diff_prefactor)*dt/Cl - (Tp-substrate_temperature)*tau_s*dt;

         // Calculate new electron and lattice temperatures
         for(unsigned int cell=0; cell<ltmp::internal::attenuation_array.size(); ++cell){

            const double Te = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0] + delta_temperature_array[2*cell+0];
            const double Tp = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1] + delta_temperature_array[2*cell+1];
            // std::cout << Te << ", " <<  Tp << std::endl;
            root_temperature_array[2*cell+0] = sqrt(Te);
            root_temperature_array[2*cell+1] = sqrt(Tp);
         }

         // optionally output cell data
         // if(ltmp::internal::output_microcell_data) ltmp::internal::write_cell_temperature_data();

         return;

      }

   } // end of namespace internal
} // end of namespace ltmp
