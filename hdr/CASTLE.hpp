



#ifndef CASTLE_H_
#define CASTLE_H_


// ==============================================
//   Coupled Atomistic and Spintronic Thermal Lattice Ensemble
//
//         =========       ========     =========   ============   ||           =========
//        ||             ||        ||   ||               ||        ||          ||
//        ||             ||        ||   ||               ||        ||          ||
//        ||             ||        ||   ||               ||        ||          ||
//        ||             || ====== ||   =========        ||        ||           =========
//        ||             ||        ||           ||       ||        ||          ||
//        ||             ||        ||           ||       ||        ||          ||
//        ||             ||        ||           ||       ||        ||          ||
//         =========                    =========                   =========   =========
//

//
//   This file is part of the VAMPIRE open source package and its subset CASTLE under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans and J L Ross 2021. All rights reserved.
//
//   Email: richard.evans@york.ac.uk and jackson.ross@york.ac.uk
//
//
//==============================================

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <ctime>
#include <random>
//#include <filesystem>
#include <omp.h> 


#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "constants.hpp"
#include "errors.hpp"
#include "stopwatch.hpp"
#include "unitcell.hpp"
#include "program.hpp" 
#include "cells.hpp"


namespace CASTLE {

    // input and material parameters
    extern bool CASTLE_program;
    extern bool CASTLE_output_data; //output position and velocity data
    extern bool equilibrium_step;

    extern int lattice_atoms; //number of lattice atoms
    extern int conduction_electrons; //number of conduction electrons
    extern double temperature;

    extern int velocity_verlet_step(double time_step);
    extern double lattice_height; //Angstroms
    extern double lattice_width;  //Angsrtoms
    extern double lattice_depth;  //Angstroms
    extern double atomic_size;    //Angstroms
    extern double screening_depth;//Angstroms
    extern double x_unit_size;
    extern double y_unit_size;
    extern double z_unit_size;

    //simulation variables
    extern double total_time_steps;
    extern double loop_time;
    extern int    CASTLE_output_rate; //output velocity and position data at this multiple
    extern int chosen_electron;
    extern double dt;
    extern double v_f; //meters
    extern double E_f; //meters
    extern double E_f_A;
    extern double mu_f; //meters
    extern double n_f; //meters
    extern double atomic_mass;
    extern double mu_r; //inverse reduced mass in reduced units
    extern double combined_mass; //inverse with reduced units
    extern double Tr; // inverse seconds

    extern double e_a_neighbor_cutoff;
    extern double e_a_coulomb_cutoff;
    extern double e_e_neighbor_cutoff;
    extern double e_e_coulomb_cutoff;
    extern double a_a_neighbor_cutoff;
    // extern double a_a_coulomb_cutoff;

   // extern int num_cells;

    //integration variables
    extern long long int current_time_step;
    extern double CASTLE_real_time;
   
    extern std::vector<double> atom_anchor_position;
    extern std::vector<double> atom_position;
   // extern std::vector<double> new_atom_position;
    //extern std::vector<double> atom_velocity; //Angstroms
    // extern std::vector<double> new_atom_velocity;
    // extern std::vector<double> atom_force;   //Angstroms
    // extern std::vector<double> new_atom_force;
    extern std::vector<double> atom_potential;
    extern std::vector<double> new_atom_potential;
    extern std::vector<std::vector<int> > atomic_nearest_electron_list;
    extern std::vector<std::vector<int> > atomic_nearest_atom_list;

    extern std::vector<double> electron_position; //Angstroms
    extern std::vector<double> new_electron_position;
    extern std::vector<double> electron_velocity; //Angstroms
    extern std::vector<double> new_electron_velocity;
    extern std::vector<double> electron_force;   //Angstroms
    extern std::vector<double> new_electron_force;
    extern std::vector<double> electron_potential; //A-1
    extern std::vector<double> new_electron_potential;
    extern std::vector<std::vector<int> > electron_nearest_electron_list;
    extern std::vector<std::vector<int> > electron_nearest_atom_list;
    extern std::vector<double> mean_radius;
    //outputs
   
    extern double TEPE; //Angstroms
    extern double TEKE; //Angstroms
    extern double TLE; //Angstroms
    // extern double TLKE; //Angstroms
    
    extern double MEPE; //meters
    extern double MEKE; //meters
    // extern double MLKE; //meters
    extern double MLE; //meters

    extern int x_flux;
    extern int y_flux;
    extern int z_flux;
    extern std::string time_stamp;
    extern std::ofstream lattice_output;
  
    extern std::ofstream electron_position_output_down;
    extern std::ofstream electron_velocity_output;
    extern std::ofstream mean_data;
    
    //control functions
    extern void initialize();
    extern void initialize_positions();
    extern void initialize_lattice();
    extern void initialize_electrons();
    extern void initialize_atoms();
    extern void initialize_forces();
    extern void initialize_electron_interactions();
    extern void initialize_atomic_interactions();
    extern void initialize_electron_atom_interactions();
    extern void initialize_velocities();
    extern void create();
    extern void output_data();

    extern void setup_output();
    extern void update_position();
    extern void update_dynamics();

    extern void e_a_coulomb(const int& e, const int& array_index, double& e_x_force, double& e_y_force, double& e_z_force, double& EPE);
              //  double& a_x_force, double& a_y_force, double& a_z_force, double& EPE, double& LPE);
   
   extern void neighbor_e_a_coulomb(const int& e, const int& array_index, double& e_x_force, double& e_y_force, double& e_z_force, double& EPE);
                // double& a_x_force, double& a_y_force, double& a_z_force, double& EPE, double& LPE);
    
    extern void e_e_coulomb(const int e, const int array_index, double& e_x_force, double& e_y_force, double& e_z_force, \
                double& EPE);
    extern void neighbor_e_e_coulomb(const int e, const int array_index, double& e_x_force, double& e_y_force, double& e_z_force, \
                double& EPE);
    
    extern void a_a_coulomb(const int a, const int array_index, \
                double& a_x_force, double& a_y_force, double& a_z_force, double& LPE);
    extern void neighbor_a_a_coulomb(const int a, const int array_index, \
                double& a_x_force, double& a_y_force, double& a_z_force, double& LPE);

    extern void update_velocity(int& array_index, double& EKE);
 
 /*
    extern double e_a_scattering(int e, int a, const double& l_x, const double& l_y, const double& l_z);
    extern double e_p_scattering(int e, int a, const double& x_distance, const double& y_distance, const double& z_distance);
    extern double electron_applied_voltage(int array_index, double& x_force, double& y_force, double& z_force);
    extern double reinitialize_electron_conserve_momentum(std::vector<double>& captured_electron_list);
*/
}



#endif //CASTLE_H_
