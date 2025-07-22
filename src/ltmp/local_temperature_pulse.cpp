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

   inline double einstein_model_phonon_heat_capacity(double T_D_over_T, double phonon_heat_capacity) {
      if(T_D_over_T > 24.0) T_D_over_T = 24.0;
      return phonon_heat_capacity*4.0*M_PI*M_PI*M_PI*M_PI/(5.0*T_D_over_T*T_D_over_T*T_D_over_T);
      }



   inline double phonon_temperature_projector_step(double Tp_init, double delta_Jp, int cell ){
      const double einstein_temp = ltmp::internal::Einstein_temperature[cell]; //for readability and speed
      const double phonon_heat_capacity = ltmp::internal::phonon_heat_capacity[cell];
      if(einstein_temp == 0.0) return phonon_heat_capacity;
      
      double T_D_over_T = einstein_temp / Tp_init;

      double phonon_temp_dep = (T_D_over_T >= 12.0) ? einstein_model_phonon_heat_capacity(T_D_over_T, phonon_heat_capacity) : phonon_heat_capacity*ltmp::internal::Debeye_phonon_constant.at(int(floor(T_D_over_T*1000.0)));
      
      double projected_Tp = Tp_init + delta_Jp/phonon_temp_dep;
      double projected_T_D_over_T = einstein_temp / projected_Tp;
      double corrected_phonon_temp_dep = (projected_T_D_over_T >= 12.0) ? einstein_model_phonon_heat_capacity(projected_T_D_over_T, phonon_heat_capacity) : phonon_heat_capacity*ltmp::internal::Debeye_phonon_constant.at(int(floor(projected_T_D_over_T*1000.0)));
      return (corrected_phonon_temp_dep+phonon_temp_dep)*0.5;

   }

   inline double phonon_temperature_projector_step(double Tp_init, double delta_Jp, int cell, double sink_deltaTp ){
      const double einstein_temp = ltmp::internal::Einstein_temperature[cell]; //for readability and speed
      const double phonon_heat_capacity = ltmp::internal::phonon_heat_capacity[cell];
      if(einstein_temp == 0.0) return phonon_heat_capacity;

      double T_D_over_T = einstein_temp / Tp_init;

      double phonon_temp_dep = (T_D_over_T >= 12.0) ? einstein_model_phonon_heat_capacity(T_D_over_T, phonon_heat_capacity) : phonon_heat_capacity*ltmp::internal::Debeye_phonon_constant.at(int(floor(T_D_over_T*1000.0)));
      
      double projected_Tp = Tp_init + delta_Jp/phonon_temp_dep - sink_deltaTp;
      double projected_T_D_over_T = einstein_temp / projected_Tp;
      double corrected_phonon_temp_dep = (projected_T_D_over_T >= 12.0) ? einstein_model_phonon_heat_capacity(projected_T_D_over_T, phonon_heat_capacity) : phonon_heat_capacity*ltmp::internal::Debeye_phonon_constant.at(int(floor(projected_T_D_over_T*1000.0)));
      return (corrected_phonon_temp_dep+phonon_temp_dep)*0.5;

   }

      void calculate_local_temperature_pulse(const double time_from_start) {

         const double i_pump_time = 1.0/ltmp::internal::pump_time;
         const double reduced_time = (time_from_start - 2.0*ltmp::internal::pump_time)*i_pump_time;
         const double four_ln_2 = 2.77258872224; // 4 ln 2
         // 2/(delta sqrt(pi/ln 2))*0.1, delta = 10 nm, J/m^2 -> mJ/cm^2 (factor 0.1)
         const double two_delta_sqrt_pi_ln_2 = 0.9394372787;
         const double gaussian = exp(-four_ln_2*reduced_time*reduced_time);
         const double pump= 1e10*ltmp::internal::pump_power*two_delta_sqrt_pi_ln_2*gaussian*i_pump_time/penetration_depth;

         if(sim::enable_laser_torque_fields) sim::laser_torque_strength = gaussian;

         const double dt = ltmp::internal::dt;

         for(unsigned int cell=0; cell < ltmp::internal::attenuation_array.size(); ++cell) {

            const double Te = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0];// + delta_temperature_array[2*cell+0];
            const double Tp = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1];// + delta_temperature_array[2*cell+1];
            // if(cell == ltmp::internal::attenuation_array.size()-1) std::cout << root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1] << ", " << delta_temperature_array[2*cell+1] << std::endl;
            // std::cout << Tp << std::endl;
            const double delta_Jp = (ltmp::internal::electron_phonon_coupling_constant[cell]*(Te-Tp))*dt;
            const double phonon_temp_dep = phonon_temperature_projector_step(Tp, delta_Jp, cell );
            root_temperature_array[2*cell+0] = sqrt(Te + (ltmp::internal::electron_phonon_coupling_constant[cell]*(Tp-Te) + pump*attenuation_array[cell])*dt/(ltmp::internal::electron_heat_capacity[cell]*Te));
            root_temperature_array[2*cell+1] = sqrt(Tp + delta_Jp/ phonon_temp_dep);

         }

         double Te = root_temperature_array[2*0+0]*root_temperature_array[2*0+0];
         double Tp = root_temperature_array[2*0+1]*root_temperature_array[2*0+1];

         // double dTe_diff_prefactor = 0.0;
         // double dTp_diff_prefactor = 0.0;
         // // calculate heat transfer from neighbouring cells
         double dTe_diff = 0.0;
         double dTp_diff = 0.0;
         for(int id = ltmp::internal::cell_neighbour_start_index[0]; id < ltmp::internal::cell_neighbour_end_index[0]; ++id) {
            const int ncell = ltmp::internal::cell_neighbour_list[id]; // neighbour cell id
            double nTe = root_temperature_array[2*ncell+0]*root_temperature_array[2*ncell+0];
            double nTp = root_temperature_array[2*ncell+1]*root_temperature_array[2*ncell+1];
            
            dTe_diff += (nTe - Te)*(ltmp::internal::electron_thermal_conductivity[ncell]+ltmp::internal::electron_thermal_conductivity[0])*0.5;
            dTp_diff += (nTp - Tp)*(ltmp::internal::phonon_thermal_conductivity[ncell]+ltmp::internal::phonon_thermal_conductivity[0])*0.5;
            // dTe_diff_prefactor += ltmp::internal::electron_thermal_conductivity[0] + ltmp::internal::electron_thermal_conductivity[ncell]);
            // dTp_diff_prefactor += 0.5*(Tp*ltmp::internal::phonon_thermal_conductivity[0] + ltmp::internal::phonon_thermal_conductivity[ncell]);
         }
         double delta_Jp = dTp_diff*dt;
         double phonon_temp_dep = phonon_temperature_projector_step(Tp, delta_Jp, 0 );

         delta_temperature_array[2*0+0] = ( 1*dTe_diff)*dt/(ltmp::internal::electron_heat_capacity[0]*Te);
         delta_temperature_array[2*0+1] = delta_Jp/phonon_temp_dep;
         // std::cout << delta_temperature_array[2*0+1] << ", " << dTp_diff*dt/ltmp::internal::phonon_heat_capacity[0] << ", " << dTp_diff << std::endl;
         for(unsigned int cell=1; cell < ltmp::internal::attenuation_array.size()-1; ++cell) {

            Te = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0];
            Tp = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1];

            // calculate heat transfer from neighbouring cells
            dTe_diff = 0.0;
            dTp_diff = 0.0;
            // double ddTe1 = 0.0;
            
            int id = ltmp::internal::cell_neighbour_start_index[cell];
            int ncell;
            double ddTe1 = root_temperature_array[2*ltmp::internal::cell_neighbour_list[id]+0]*root_temperature_array[2*ltmp::internal::cell_neighbour_list[id]+0];
            double ddTp1 = root_temperature_array[2*ltmp::internal::cell_neighbour_list[id]+1]*root_temperature_array[2*ltmp::internal::cell_neighbour_list[id]+1];
            for( id; id < ltmp::internal::cell_neighbour_end_index[cell]; ++id) {
               ncell = ltmp::internal::cell_neighbour_list[id]; // neighbour cell id
               double nTe = root_temperature_array[2*ncell+0]*root_temperature_array[2*ncell+0];
               double nTp = root_temperature_array[2*ncell+1]*root_temperature_array[2*ncell+1];
               // dTe_diff += (nTe - Te)*(ltmp::internal::electron_thermal_conductivity[ncell]+ltmp::internal::electron_thermal_conductivity[cell])*0.5;
               dTp_diff += (nTp - Tp)*(ltmp::internal::phonon_thermal_conductivity[ncell]+ltmp::internal::phonon_thermal_conductivity[cell])*0.5;
               
            }
            double ddTe2 = root_temperature_array[2*ncell+0]*root_temperature_array[2*ncell+0];
            double ddTp2 = root_temperature_array[2*ncell+1]*root_temperature_array[2*ncell+1];

            // delta_Jp = (0.5*dTp_diff+0.25*(ddTp2+ddTp1-2.0*Tp)*(ltmp::internal::phonon_thermal_conductivity[  ltmp::internal::cell_neighbour_list[ltmp::internal::cell_neighbour_start_index[cell]]]+ltmp::internal::phonon_thermal_conductivity[  ncell]+2.0*ltmp::internal::phonon_thermal_conductivity[cell]))*dt;
            delta_Jp = 0.5*dTp_diff*dt;
            phonon_temp_dep = phonon_temperature_projector_step(Tp, delta_Jp, cell );
            // std::cout << phonon_temp_dep << ", " << (ltmp::internal::electron_phonon_coupling_constant[cell]*(Te-Tp))*dt/ phonon_temp_dep <<  std::endl; 
            
            delta_temperature_array[2*cell+0] = ((0.5*dTe_diff+0.25*(ddTe2+ddTe1-2.0*Te)*(ltmp::internal::electron_thermal_conductivity[ltmp::internal::cell_neighbour_list[ltmp::internal::cell_neighbour_start_index[cell]]]+ltmp::internal::electron_thermal_conductivity[ncell]+2.0*ltmp::internal::electron_thermal_conductivity[cell])))*dt/(ltmp::internal::electron_heat_capacity[cell]*Te);
            delta_temperature_array[2*cell+1] = delta_Jp/phonon_temp_dep;
            
         } // end of cell loop

         int cell = ltmp::internal::attenuation_array.size()-1;
         Te = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0];
         Tp = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1];

         // calculate heat transfer from neighbouring cells
         dTe_diff = 0.0;
         dTp_diff = 0.0;
         for(int id=ltmp::internal::cell_neighbour_start_index[cell]; id<ltmp::internal::cell_neighbour_end_index[cell]; ++id) {
            const int ncell = ltmp::internal::cell_neighbour_list[id]; // neighbour cell id
            double nTe = root_temperature_array[2*ncell+0]*root_temperature_array[2*ncell+0];
            double nTp = root_temperature_array[2*ncell+1]*root_temperature_array[2*ncell+1];
            dTe_diff += (nTe - Te)*(ltmp::internal::electron_thermal_conductivity[ncell]+ltmp::internal::electron_thermal_conductivity[cell])*0.5;
            dTp_diff += (nTp - Tp)*(ltmp::internal::phonon_thermal_conductivity[ncell]+ltmp::internal::phonon_thermal_conductivity[cell])*0.5;
            // std::cout << "cell: " << cell << ", my Tp: " << Tp << ", ncell Tp: " << nTp << ", ncell: " <<  ncell << ", dTp: " << 2.0*dTp_diff*dt/ltmp::internal::phonon_heat_capacity[cell] << std::endl;
         }
         delta_Jp = (1*dTp_diff)*dt;
         phonon_temp_dep = phonon_temperature_projector_step(Tp, delta_Jp, cell, (Tp-equilibration_temperature)*Tcool*dt );
         //if(phonon_temp_dep <= 0.0) {std::cout << phonon_temp_dep << std::endl; exit(1);}
         delta_temperature_array[2*cell+0] = (1*dTe_diff)*dt/(ltmp::internal::electron_heat_capacity[cell]*Te);
         delta_temperature_array[2*cell+1] = delta_Jp/phonon_temp_dep - (Tp-equilibration_temperature)*Tcool*dt;
        // std::cout << " substrate cooling: " << - (Tp-equilibration_temperature)*Tcool*dt << ", deltaTp: " << - (Tp-equilibration_temperature) << ", " << Tcool << std::endl;
         // Calculate new electron and lattice temperatures
         for(unsigned int cell=0; cell < ltmp::internal::attenuation_array.size(); ++cell) {

            const double Te = root_temperature_array[2*cell+0]*root_temperature_array[2*cell+0] + delta_temperature_array[2*cell+0];
            const double Tp = root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1] + delta_temperature_array[2*cell+1];
            // if(cell == ltmp::internal::attenuation_array.size()-1) std::cout << root_temperature_array[2*cell+1]*root_temperature_array[2*cell+1] << ", " << delta_temperature_array[2*cell+1] << std::endl;
            root_temperature_array[2*cell+0] = sqrt(Te);
            root_temperature_array[2*cell+1] = sqrt(Tp);
         }

         // optionally output cell data
         // if(ltmp::internal::output_microcell_data) ltmp::internal::write_cell_temperature_data();

         return;

      }

   } // end of namespace internal
} // end of namespace ltmp
