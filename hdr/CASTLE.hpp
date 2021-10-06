



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
    extern double conduction_electrons; //number of conduction electrons


    extern int velocity_verlet_step(double dt);
    extern   double lattice_height;
    extern double lattice_width;
    extern    double lattice_depth;

    extern   double atomic_size;
    extern double screening_depth;
    extern double v_f;
    extern double E_f;

    extern double TKE;
    extern double TPE;
    extern int total_spin_up;
    extern int total_spin_down;

    extern std::vector<std::vector<bool> > symmetry_list;

    extern double temperature;
    extern double mu_f;

    extern  double total_time_steps;
    extern double loop_time;
    extern  double dt;
    extern int current_time_step;

    extern int CASTLE_output_rate; //output velocity and position data at this multiple

    extern std::vector<double> electron_position; //superarray of past locations for each step
    extern std::vector<double> new_electron_position;
    extern std::vector<double> electron_velocity; //superarray for past velocities for each step
    extern std::vector<double> new_electron_velocity;
    extern std::vector<double> electron_acc;  
    extern std::vector<double> new_acc_array;
    extern std::vector<double> mean_data_array;
    extern std::vector<double> lattice_electrons;
    extern std::vector<bool>   conduction_electron_spin;
    extern std::vector<bool> lattice_electron_spin;

    extern std::vector<double> atom_position;
    
    extern void initialize(const double num_electrons, const double num_atoms);
   
    extern void create();
    extern void update_data();
    extern void create_gif();

    extern std::ofstream lattice_output;
    extern std::ofstream electron_position_output_up;
    extern std::ofstream electron_position_output_down;
    extern std::ofstream electron_velocity_output;
    extern std::ofstream mean_data;
    extern std::ofstream electron_spin_output;
    

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
