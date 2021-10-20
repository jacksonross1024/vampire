



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
  //  extern double total_electrons; //lattice + conduction electrons
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
    extern long double dt;
    extern long double v_f; //meters
    extern long double E_f; //meters
    extern long double mu_f; //meters
    extern long double n_f; //meters

   // extern int num_cells;

    //integration variables
    extern int current_time_step;
    extern double CASTLE_real_time;
   
    extern std::vector<double> atom_position;
    extern std::vector<long double> electron_position; //Angstroms
    extern std::vector<long double> new_electron_position;
    extern std::vector<long double> electron_velocity; //Angstroms
    extern std::vector<long double> new_electron_velocity;
    extern std::vector<long double> electron_force;   //Angstroms
    extern std::vector<long double> new_force_array;
    extern std::vector<double> mean_data_array;
    extern std::vector<std::vector<int> > nearest_neighbor_list;
   // extern std::vector<double> lattice_electrons;
   // extern std::vector<bool>   conduction_electron_spin;
   // extern std::vector<bool>   lattice_electron_spin;
   // extern std::vector<int> electrons_per_cell; //size of number of cells
   // extern std::vector<int>  electron_cell;       // cell address for each electron
    extern std::vector<std::vector<bool> > symmetry_list;

    //outputs
    extern long double TKE;
    extern long double TPE;
    extern long double MPE; //meters
    extern long double MKE; //meters
    extern int x_flux;
    extern std::vector <long double> velocity_length_hist;
  //  extern int total_spin_up;
  //  extern int total_spin_down;
    extern std::ofstream lattice_output;
  //  extern std::ofstream electron_position_output_up;
    extern std::ofstream electron_position_output_down;
    extern std::ofstream electron_velocity_output;
    extern std::ofstream mean_data;
  //  extern std::ofstream electron_spin_output;

    
    //control functions
    extern void initialize();
    extern void initialize_lattice();
    extern void initialize_electrons();
    extern void initialize_forces();
    extern void create();
    extern void output_data();

    extern void setup_output();
    extern void update_position();
    extern void update_dynamics();
    extern long double electron_e_a_coulomb(int array_index, long double& x_force, long double& y_force, long double& z_force, const long double& x, const long double& y, const long double& z);
    extern long double electron_e_e_coulomb(int e, int array_index, long double& x_force, long double& y_force, long double& z_force, const long double& x, const long double& y, const long double& z);
    extern long double neighbor_e_e_coulomb(int e, int array_index, long double& x_force, long double& y_force, long double& z_force, const long double& x, const long double& y, const long double& z);
    extern long double update_velocity(int array_index);
    extern long double electron_applied_voltage(int array_index, long double& x_force, long double& y_force, long double& z_force);

/*  one day your time will come
    struct atom {
        double x_position;
        double y_position;
        double z_positon;

        double x_velocity;
        double y_velocity;
        double z_velocity;

        double x_acc;
        double y_acc;
        double z_acc;

        std::vector<double> nearest_atom;
        std::vector<double> nearest_electron;

    }; */ 

}



#endif //CASTLE_H_
