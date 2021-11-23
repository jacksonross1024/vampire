



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
//#include <omp.h> 


#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "constants.hpp"
#include "errors.hpp"
#include "stopwatch.hpp"



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

    //simulation variables
    extern double total_time_steps;
    extern double loop_time;
    extern int    CASTLE_output_rate; //output velocity and position data at this multiple
    extern long double chosen_electron;
    extern long double dt;
    extern long double v_f; //meters
    extern long double E_f; //meters
    extern long double E_f_A;
    extern long double mu_f; //meters
    extern long double n_f; //meters
    extern long double atomic_mass;
    extern long double mu_r; //inverse reduced mass in reduced units
    extern long double combined_mass; //inverse with reduced units


    extern long double e_a_neighbor_cutoff;
    extern long double e_a_coulomb_cutoff;
    extern long double e_e_neighbor_cutoff;
    extern long double e_e_coulomb_cutoff;
    extern long double a_a_neighbor_cutoff;
    extern long double a_a_coulomb_cutoff;

   // extern int num_cells;

    //integration variables
    extern int current_time_step;
    extern double CASTLE_real_time;
   
    extern std::vector<double> atom_position;
    extern std::vector<double> new_atom_position;
    extern std::vector<long double> atom_velocity; //Angstroms
    extern std::vector<long double> new_atom_velocity;
    extern std::vector<long double> atom_force;   //Angstroms
    extern std::vector<long double> new_atom_force;
    extern std::vector<long double> atom_potential;
    extern std::vector<std::vector<int> > atomic_nearest_electron_list;
    extern std::vector<std::vector<int> > atomic_nearest_atom_list;

    extern std::vector<long double> electron_position; //Angstroms
    extern std::vector<long double> new_electron_position;
    extern std::vector<long double> electron_velocity; //Angstroms
    extern std::vector<long double> new_electron_velocity;
    extern std::vector<long double> electron_force;   //Angstroms
    extern std::vector<long double> new_electron_force;
    extern std::vector<long double> electron_potential; //A-1
    extern std::vector<long double> new_electron_potential;
    extern std::vector<std::vector<int> > electron_nearest_electron_list;
    extern std::vector<std::vector<int> > electron_nearest_atom_list;

    //outputs
   
    extern long double TEPE; //Angstroms
    extern long double TEKE; //Angstroms
    extern long double TLPE; //Angstroms
    extern long double TLKE; //Angstroms
    
    extern long double MEPE; //meters
    extern long double MEKE; //meters
    extern long double MLKE; //meters
    extern long double MLPE; //meters

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

    extern void e_a_coulomb(const int a, const int array_index, long double& e_x_force, long double& e_y_force, long double& e_z_force, \
                long double& a_x_force, long double& a_y_force, long double& a_z_force, long double& EPE, long double& LPE);
    extern void neighbor_e_a_coulomb(const int a, const int array_index, long double& e_x_force, long double& e_y_force, long double& e_z_force, \
                long double& a_x_force, long double& a_y_force, long double& a_z_force, long double& EPE, long double& LPE);
    
    extern void e_e_coulomb(const int e, const int array_index, long double& e_x_force, long double& e_y_force, long double& e_z_force, \
                long double& EPE);
    extern void neighbor_e_e_coulomb(const int e, const int array_index, long double& e_x_force, long double& e_y_force, long double& e_z_force, \
                long double& EPE);
    
    extern void a_a_coulomb(const int a, const int array_index, \
                long double& a_x_force, long double& a_y_force, long double& a_z_force, long double& LPE);
    extern void neighbor_a_a_coulomb(const int a, const int array_index, \
                long double& a_x_force, long double& a_y_force, long double& a_z_force, long double& LPE);

    extern void update_velocity(int array_index, long double& EKE, long double& LKE);
 
 /*
    extern long double e_a_scattering(int e, int a, const long double& l_x, const long double& l_y, const long double& l_z);
    extern long double e_p_scattering(int e, int a, const long double& x_distance, const long double& y_distance, const long double& z_distance);
    extern long double electron_applied_voltage(int array_index, long double& x_force, long double& y_force, long double& z_force);
    extern long double reinitialize_electron_conserve_momentum(std::vector<long double>& captured_electron_list);
*/
}



#endif //CASTLE_H_
