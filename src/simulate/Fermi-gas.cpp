//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2021 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the Geta General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Geta
//  General Public License for more details.
//
//  You should have received a copy of the Geta General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
///
/// @file
/// @brief Contains functions to model the conduction electrons as a Fermi Gas
///
/// @details None
///
/// @section notes Implementation Notes
/// This is a list of other notes, not related to functionality but rather to implementation.
/// Also include references, formulae and other notes here.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2021. All Rights Reserved.
///
/// @section info File Information
/// @author  J L Ross, jackson.ross@york.ac.uk
/// @version 1.0
/// @date    21/09/21
/// @internal
///	Created:		21/09/21
///	Revision:	  ---
///=====================================================================================
///
// Standard Libraries
#include <algorithm>
#include <cmath>
#include <math.h>
#include <iostream>

#include "internal.hpp"

// Vampire Header files
#include "anisotropy.hpp"
#include "atoms.hpp"
#include "exchange.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "constants.hpp"
#include "cells.hpp"


#include "sim.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "create_atoms_class.hpp"

namespace sim{

    //Fermi function defined as 
    //
    // F(eps) = 1 / [ exp(beta(eps - mu))]
    //
    void initialize_fermi_gas() {

    fermi_size = atoms::num_atoms;
  /*  #ifdef MPICF
             MPI_Allreduce(MPI_IN_PLACE, &fermi_size , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif
    */    
        fermi_electrons = fermi_size * conduction_electrons * vmpi::num_processors;
      //  std::cout << "fermi_electrons: " << fermi_electrons << std::endl;
        //calculate volume
        //fermi_volume = cells::atomic_volume * fermi_size; need cell initialization constants extracted
      // fermi_volume = cs::system_dimensions[0] * cs::system_dimensions[1] * cs::system_dimensions[2] * 1e-9 * 1e-9 * 1e-9; //square and rectangle only
        fermi_volume = 27e-27;
      //  std::cout << "fermi_volume "  << fermi_volume << std::endl;

        mu_0 = ((constants::h * constants::h) / ( 2 * constants::m_e)) * pow(((3 * fermi_electrons) / (8 * M_PI * fermi_volume)), 0.6666666667);
        E_0 = 3 * fermi_electrons * mu_0 / 5;
        P_0 = 2 * E_0 / (3 * fermi_volume);
        fermi_density = fermi_electrons / fermi_volume;
  
    }

    void fermi_calculations( double sim_temp) {
  
        temperature_variable = sim_temp;
    
        fermi_distribution_array.resize(2 * fermi_size, 0.0);
        beta_variable = 1 / (temperature_variable * constants::kB);
 
        eta = 1 / (beta_variable * mu_0);
        mu = mu_0 * (1 - (((M_PI * M_PI) / 12) * eta * eta) + (((M_PI * M_PI * M_PI) / 24) * eta * eta * eta));
        const long double fermi_range = 2 * mu; //16 * constants::kB * temperature_variable;
        
      
        long double StepSize = fermi_range / (2 * fermi_size);
        unsigned int esp_index = 0;
        long double state_counter = 0;
        std::ofstream distribution;

      
        
            
        distribution.open("Fermi_Statistics/Distribution/" + std::to_string(int(temperature_variable)));
    

        for (long double eps_variable = 0; eps_variable < fermi_range; eps_variable += StepSize) {
            
            state_counter += fermi_distribution_array[esp_index] = fermi_distribution(eps_variable);

            distribution << eps_variable << "    " << fermi_distribution_array[esp_index] << std::endl;
            esp_index += 1;
            
        }
   
}

long double fermi_distribution(long double eps) {
        long double x_var = beta_variable * (eps - mu);
        long double inverse_function = exp(x_var) + 1;
        long double function = 1 / inverse_function;

        return function;
}



std::string output_fermi_energy() {
        std::ostringstream res;

        if(vout::custom_precision){
            res.precision(vout::precision);
          if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
        }
        vout::fixed_width_output result(res,vout::fw_size); 

        const double energy = E_0 * (1 + (((5 * M_PI * M_PI) / 12) * eta * eta) - (((5 * M_PI * M_PI * M_PI) / 24) * eta * eta * eta));  
            result << energy;

        return result.str();
}

std::string output_fermi_pressure() {
        std::ostringstream res;

        if(vout::custom_precision){
            res.precision(vout::precision);
          if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
        }
        vout::fixed_width_output result(res,vout::fw_size); 

        const double pressure = P_0 * (1 + (((5 * M_PI * M_PI) / 12) * eta * eta) - (((5 * M_PI * M_PI * M_PI) / 24) * eta * eta * eta));  
            result << pressure;

        return result.str();
}

std::string output_fermi_Cv() {
        std::ostringstream res;

        if(vout::custom_precision){
            res.precision(vout::precision);
          if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
        }
        vout::fixed_width_output result(res,vout::fw_size); 

        const double Cv = M_PI * M_PI * fermi_electrons * constants::kB * eta / 2; 
        result << Cv;

        return result.str();
}

std::string output_relativistic_fermi_energy() {
        std::ostringstream res;

        if(vout::custom_precision){
            res.precision(vout::precision);
          if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
        }
        vout::fixed_width_output result(res,vout::fw_size); 

        const double energy = pow((3 / M_PI), 0.66666667) * 3 * constants::h * constants::h * fermi_volume * pow((fermi_density), 1.666666666667) / (40 * constants::m_e);
        result << energy;

        return result.str();    
}
std::string output_relativistic_fermi_pressure() {
        std::ostringstream res;

        if(vout::custom_precision){
            res.precision(vout::precision);
          if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
        }
        vout::fixed_width_output result(res,vout::fw_size); 

        const double pressure = 2 * ((pow((3 / M_PI), 0.66666667) * 3 * constants::h * constants::h * fermi_volume * pow((fermi_density), 1.666666666667) / (40 * constants::m_e)) / (3 * fermi_volume));
        result << pressure;

        return result.str();    

}
} 


