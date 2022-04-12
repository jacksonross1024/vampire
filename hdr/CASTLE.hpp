



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
#include <forward_list>


#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "constants.hpp"
#include "errors.hpp"
#include "stopwatch.hpp"
#include "random.hpp"



namespace CASTLE {

    // input and material parameters
    extern bool CASTLE_program;
    extern bool CASTLE_output_data; //output position and velocity data
    
    extern bool equilibrium_step;
    extern bool applied_voltage_sim;
    extern bool heat_pulse_sim;

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
    extern int    CASTLE_output_rate; 
    extern int CASTLE_MD_rate; //output velocity and position data at this multiple
  
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
    extern double a_heat_capacity; // AJ/K/mol -> AJ/K
    extern double e_heat_capacity; // AJ/K/mol -> AJ/K
    extern double a_heat_capacity_i; // K/AJ
    extern double e_heat_capacity_i; // K/AJ
    extern double a_specific_heat; //AJ/K/particle
    extern double e_specific_heat;
    extern double a_specific_heat_i;
    extern double e_specific_heat_i;
    extern double zero_pt_lattice_e;
    extern double phonon_energy;
    extern double new_phonon_energy;

    extern bool ee_coupling;
    extern bool ea_coupling;
    extern double ea_rate;
    extern double ee_rate;
    extern double ee_coupling_strength;
    extern double ea_coupling_strength;


    extern double e_a_neighbor_cutoff;
    extern double e_a_coulomb_cutoff;
    extern double e_e_neighbor_cutoff;
    extern double e_e_coulomb_cutoff;
    extern double a_a_neighbor_cutoff;
    //extern double a_a_neighbor_cutoff;
    // extern double a_a_coulomb_cutoff;

    extern double applied_voltage;
    extern double power_density;
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
    extern std::vector<bool> external_interaction_list;
    extern int external_interaction_list_count;
    extern std::vector<double> mean_radius;
    extern std::vector<std::vector<int> > electron_ee_scattering_list;
    extern std::vector<std::vector<int> > electron_ea_scattering_list;

    extern std::vector<double> phonon_distribution;
    extern MTRand_closed uniform_random;
    extern MTRand_int32 int_random;
    extern std::vector<MTRand_closed> omp_uniform_random;
    extern std::vector<MTRand_int32> omp_int_random;
   // extern std::vector<MTRand> omp_gaussian_random;
    //outputs
   
    extern double TEPE; //Angstroms
    extern double TEKE; //Angstroms
    extern double TLE; //Angstroms
    extern double Tp;
    extern double Te;
    
    extern double MEPE; //meters
    extern double MEKE; //meters
    // extern double MLKE; //meters
    extern double MLE; //meters

    extern int x_flux;
    extern int y_flux;
    extern int z_flux;
    extern int e_a_scattering_count;
    extern int e_e_scattering_count;
    extern int a_a_scattering_count;

    extern double TTMe;
    extern double d_TTMe;
    extern double TTMp;
    extern double d_TTMp;
    extern double G;
    
    extern std::string time_stamp;
    extern std::ofstream lattice_output;
  
    extern std::ofstream electron_position_output_down;
    extern std::ofstream electron_velocity_output;
    extern std::ofstream mean_data;
    extern std::ofstream temp_data;
    
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

    extern void e_a_coulomb(const int& e, const int& array_index);
              //  double& a_x_force, double& a_y_force, double& a_z_force, double& EPE, double& LPE);
   
   extern void neighbor_e_a_coulomb(const int& e, const int& array_index);
                // double& a_x_force, double& a_y_force, double& a_z_force, double& EPE, double& LPE);
    
    extern void e_e_coulomb(const int& e, const int& array_index);
    extern void neighbor_e_e_coulomb(const int& e, const int& array_index);
    
    extern void a_a_coulomb(const int a, const int array_index, \
                double& a_x_force, double& a_y_force, double& a_z_force, double& LPE);
    extern void neighbor_a_a_coulomb(const int a, const int array_index, \
                double& a_x_force, double& a_y_force, double& a_z_force, double& LPE);

    extern void electron_thermal_field(const int& e, const int& array_index, const double& EKE);
 
    extern double electron_applied_voltage(const int& e, const int& array_index, double& external_potential);
 
    extern void aa_scattering();
    extern void ea_scattering(const int& e, const int& array_index);
    extern void ee_scattering();
 /*
    extern double e_a_scattering(int e, int a, const double& l_x, const double& l_y, const double& l_z);
    extern double e_p_scattering(int e, int a, const double& x_distance, const double& y_distance, const double& z_distance);
    
    extern double reinitialize_electron_conserve_momentum(std::vector<double>& captured_electron_list);
*/
    extern double M_B_distrib(const double& epsilon, const double& beta);
    extern double B_E_distrib(const double& epsilon);
    extern void create_phonon_distribution(std::vector<double>& distribution, const double& beta);
    extern void create_phonon_distribution(const std::string& name, std::vector<double>& distribution, const double& beta);
    extern double return_phonon_distribution(const std::vector<double>& distribution);
}



#endif //CASTLE_H_
