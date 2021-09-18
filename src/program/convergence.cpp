//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//  {J L Ross 2021}
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <iostream>
#include <sstream>

// Vampire headers
#include "stats.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "errors.hpp"


//This module part of a method to measure and assist in the convergence of equilibrium and loop averaging style
//code methods

//It is a subpart of the stats.hpp class
namespace stats{

  //  convergence_statistic_t (std::string n, double criteria, unsigned int convergence_check_value) 
  
    // ------------------------
    // Function to tell convergence class to begin convergence analysis
    //-------------------------
    void convergence_statistic_t::do_converge() {
        if (!sim::calculate_program_convergence) return;
        convergence_method = true;
        counter = 0;
    }

    //-------------------------
    // End of equilibrium check; if flag sent without convergence, sends 0 counter
    //--------------------------
    void convergence_statistic_t::end_convergence() {
        if (!sim::calculate_program_convergence) return;
        if (!converged) {
         //   vout::screen_output_list.push_back(72);
       //     vout::file_output_list.push_back(72);
            vout::write_out(std::cout, convergence_output_list);
            std::cout << "Equilibrium may not have been reached. More steps may be necessary." << std::endl;
        }
        convergence_method = false;
        converged = false;
    }
    //-------------------------
    // Function to return the status past convergence status of a simulation
    //-----------------------------
    bool convergence_statistic_t::did_converge() {
                if (!sim::calculate_program_convergence) return false;

        return converged;
    }

    //Getter for convergence method
    bool convergence_statistic_t::get_method() {
                if (!sim::calculate_program_convergence) return false;

        return convergence_method;
    }

    unsigned int convergence_statistic_t::get_converged_counter() {
                if (!sim::calculate_program_convergence) return 0;

        return converged_counter;
    }
    //-----------------------------
    //Feed the exit criteria and check frequency to the class. Check they have both been given
    // ---------------------------
    void convergence_statistic_t::initialize(double criteria, unsigned int convergence_check_value) {
            convergence_criteria = criteria;
            convergence_check = convergence_check_value;

            convergence_energy.resize(2, 0.0);
            convergence_magnetisation.resize(2, 0.0);

            output_convergence_counter = false;

            converged = false;
            energy_converged = false;
            convergence_method = false;
            magnetisation_converged = false;
            
            output_convergence();
           // convergence_output_list.resize(1, 0);

            if( !convergence_criteria || !convergence_check) {
                terminaltextcolor(RED);
                std::cerr << "Programmer Error - poor convergence criteria passed to convergence module - please initialize first." << std::endl;
                terminaltextcolor(WHITE);
                zlog << zTs() << "Programmer Error - poor convergence criteria passed to convergence module - please initialize first." << std::endl;
                err::vexit();
            }
            
     }
            

    //-----------------------------
    // Function to test if a simulation has converged. 
    //-----------------------------
    void  convergence_statistic_t::is_converged(unsigned int n_steps) {
                if (!sim::calculate_program_convergence) return;

        if (err::check) std::cout << "convergence check has been called" << std::endl;
        //test the convergence of the magnetisation and energy of a simulation
        if ((abs(convergence_magnetisation[0] - convergence_magnetisation[1])) < (convergence_magnetisation[0] * convergence_criteria)) magnetisation_converged = true;
        if ((abs(convergence_energy[0] - convergence_energy[1])) < (abs(convergence_energy[0]) * convergence_criteria)) energy_converged = true;
        
        //if both are converged, the convergence is switched
        if (energy_converged && magnetisation_converged) {
            converged = true;                   //convergence set to true
            convergence_method = false;         // the method set to false so no further analysis is done
            converged_counter = counter;        // count is saved 
            counter = 0;                        // counter reset for next equilibrium check
            energy_converged = false;           //reset energy
            magnetisation_converged = false;    //reset magnetisation
           // if(output_convergence_counter)      vout::write_out(zmag, convergence_output_list);
            // std::cout << "Convergence may have occured. Exiting equilibration steps on " << converged_counter << " out of " << sim::equilibration_time << std::endl;             // adds convergence rate to output list if flag is enabled
     
         //   std::cout << "Convergence may have occured. Exiting equilibration steps on " << converged_counter << "out of " << sim::equilibration_time << std::endl;             // adds convergence rate to output list if flag is enabled
         //   std::cout << (abs(convergence_energy[0] - convergence_energy[1])) << " vs " << (abs(convergence_energy[0]) * convergence_criteria) << std::endl;
           //  std::cout << (abs(convergence_magnetisation[0] - convergence_magnetisation[1])) << " vs " << (convergence_magnetisation[0] * convergence_criteria) << std::endl;
        }
      //  return converged;
    }

    //-----------------------------
    // Function to initialize the data sets and assign to the convergence arrays. Not meant to be called outside of the class
    //-----------------------------
    void convergence_statistic_t::initialize_first_pass() {
   ////     std::cout << "fail"<< std::endl;
         update_data(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, atoms::m_spin_array, atoms::type_array, sim::temperature);
//std::cout << "fail"<< std::endl;
        
        convergence_magnetisation[0] = stats::system_magnetization.get_normalized_mean_system_magnetization_length();
     //  std::cout << "fail"<< std::endl;
        convergence_energy[0]        = stats::system_energy.get_mean_total_system_energy();
    //std::cout << "fail"<< std::endl;
    }


    //-----------------------------
    // Function to initialize the data sets and assign to the convergence arrays. Not meant to be called outside of the class
    //-----------------------------
    void convergence_statistic_t::initialize_second_pass() {

        update_data(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, atoms::m_spin_array, atoms::type_array, sim::temperature);

        convergence_magnetisation[1] = stats::system_magnetization.get_normalized_mean_system_magnetization_length();
        convergence_energy[1]        = stats::system_energy.get_mean_total_system_energy();
        
    }

    //-----------------------------
    // Function to update the data of the statistics class necessary for convergence
    //-----------------------------
    void convergence_statistic_t::update_data(const std::vector<double>& sx, // spin unit vector
               const std::vector<double>& sy,
               const std::vector<double>& sz,
               const std::vector<double>& mm,
               const std::vector<int>& mat,
               const double temperature) {

                   // set calcualtion flags in case not set by output or system
                   //no need to calculate more than necessary
        stats::calculate_system_energy = true;
        stats::calculate_system_magnetization= true;

        stats::update(sx, sy, sz, mm, mat, temperature );
    }

    //-----------------------------
    // Function to update the convergence arrays. Internal function not meant to be called outside of class
    //-----------------------------
    void convergence_statistic_t::update_convergence() {
                if (!sim::calculate_program_convergence) return;

        if (err::check) std::cout << "update convergence has been called" << std::endl;
        //First updates the data
        update_data(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, atoms::m_spin_array, atoms::type_array, sim::temperature);
  
        //adjusts the arrays once more
        //ought to be called after the first and second initialization
        convergence_magnetisation[0] = convergence_magnetisation[1];
        convergence_magnetisation[1] = stats::system_magnetization.get_normalized_mean_system_magnetization_length();
        convergence_energy[0]        = convergence_energy[1];
        convergence_energy[1]        = stats::system_energy.get_mean_total_system_energy();
    }
    
    //-----------------------------
    // Function to update the convergence arrays and counter. This is the one to call outside the class
    //-----------------------------
    void convergence_statistic_t::update_counter() {
                if (!sim::calculate_program_convergence) return;

        //if convergence setting has been turned off, this will avoid the convergence
        if (err::check == true)   std::cout << "update counter has been called" << std::endl;

        //checks convergence on criteria of convergence-check
        if ( (counter % convergence_check == 0) && (counter != 0)) {
         //   std::cout << "counter multiples of concergence-check " << counter << std::endl;
            is_converged(counter);
            update_convergence();
       //   std::cout << "magnetisation0: " << convergence_magnetisation[0] << " energy0: " << convergence_energy[0] << std::endl;
        //  std::cout << "magnetisation1: " << convergence_magnetisation[1] << " energy1: " << convergence_energy[1] << std::endl;
        }
        // second initialization
        else if (counter == (convergence_check - 1)) {
        //    std::cout << "counter 99" << counter << std::endl;
            initialize_second_pass();
         //   std::cout << "magnetisation1: " << convergence_magnetisation[1] << " energy1: " << convergence_energy[1] << std::endl;
        }
        // first initialization
        else if (counter == 0) {
         //   std::cout << "Counter 0: " << counter << std::endl;
            initialize_first_pass();
         //   std::cout << "magnetisation0: " << convergence_magnetisation[0] << " energy0: " << convergence_energy[0] << std::endl;
        }
        //move counter forward one
        counter += 1;
    }

    //-----------------------------
    // Function to output the counter value at which convergence was achieved
    //-----------------------------
    std::string convergence_statistic_t::output_convergence_rate(bool header) {
        
     //   std::cout << counter << std::endl;
           // result string stream
         std::ostringstream res;
        unsigned int count;
        if (converged) count = converged_counter;
        else count = counter;

        if(vout::custom_precision){
            res.precision(vout::precision);
          if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
        }
        vout::fixed_width_output result(res,vout::fw_size); 
        if (header) result << name + std::to_string(counter);
        else result << counter;

        return result.str();

    }

    void convergence_statistic_t::output_convergence() {
        output_convergence_counter = sim::output_convergence_counter;
        if (!convergence_output_list.size()) convergence_output_list.push_back(72);
        
    }

}

