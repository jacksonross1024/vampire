



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
#include <algorithm>
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

    // input variables
    extern int  omp_threads;
    extern bool CASTLE_program;
    extern bool CASTLE_output_data; //output position and velocity data
    
    extern bool equilibrium_step;
    extern bool applied_voltage_sim;
    extern bool heat_pulse_sim;

    extern int lattice_atoms;
    extern int conduction_electrons; //number of conduction electrons
    extern int layer_1;
    extern int layer_2;
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
  
    extern int full_int_var;
    extern int half_int_var;
    extern double boundary_conditions_cutoff;
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
    extern int ee_density;
    extern int dos_size;
    extern double dos_occ_1;
    extern double dos_occ_2;
    extern double local_dos_occ;

    extern bool ee_coupling;
    extern bool ea_coupling;
    extern double ea_rate;
    extern double ee_rate;
    extern double ee_coupling_strength;
    extern double ea_coupling_strength;
    extern double ee_scattering_angle;

    extern double e_e_integration_cutoff;
    extern double e_e_neighbor_cutoff;
    extern double e_e_coulomb_cutoff;

    extern int x_omp_cells;// = int(floor(lattice_width / 15.0));
    extern int y_omp_cells;// = int(floor(lattice_depth / 15.0));
    extern int z_omp_cells;// = int(floor(lattice_height/ 15.0));

    extern int total_cells;// = x_omp_cells*y_omp_cells*z_omp_cells;

    extern double x_step_size;// = lattice_width / double(x_omp_cells);
    extern double y_step_size;// = lattice_depth / double(y_omp_cells);
    extern double z_step_size;// = lattice_height/ double(z_omp_cells);
 
    extern double applied_voltage;
    extern double power_density;

    //integration variables
    extern long long int current_time_step;
    extern double CASTLE_real_time;
    extern int cells_per_thread;

    extern std::vector<double> electron_position; //Angstroms
    extern std::vector<double> electron_velocity; //Angstroms 
    extern std::vector<double> electron_potential; //A-1
    extern std::vector<std::vector< int> > ee_dos_hist;
    extern std::vector<std::vector< int> > global_e_dos_1;
    extern std::vector<std::vector< int> > global_e_dos_2;

    // extern std::vector<bool> electron_transport_list;
    extern std::vector<std::vector< int> > electron_integration_list;
    extern std::vector<std::vector< int> > electron_nearest_electron_list;
    extern std::vector<std::vector<double> > electron_ee_scattering_list;
    extern std::vector<std::vector< int> > cell_lattice_coordinate;
    extern std::vector<std::vector< int> > cell_integration_lists;
    extern std::vector<std::vector< int> > old_cell_integration_lists;
    extern std::vector<std::vector< int> > cell_nearest_neighbor_list;
    extern std::vector<std::vector<std::vector< int> > > lattice_cell_coordinate;
    extern std::vector<std::vector< int> > lattice_cells_per_omp;
    extern std::vector< int> escaping_electrons;
    extern std::vector<std::vector< int> > relaxation_time_hist_ee_1;
    extern std::vector<std::vector< int> > relaxation_time_hist_ee_2;
    // extern std::vector<std::vector< int> > relaxation_time_hist_ea;
    extern std::vector< int> flux_index;
    extern std::vector<MTRand> omp_uniform_random;
    extern std::vector<MTRand_int32> omp_int_random;
    //outputs
   
    extern double TEKE_1; //Angstroms
    extern double TLE_1; //Angstroms
    extern double Tp_1;
    extern double Te_1;
    extern double d_Tp_1;
    extern double d_Te_1;

    extern double TEKE_2; //Angstroms
    extern double TLE_2; //Angstroms
    extern double Tp_2;
    extern double Te_2;
    extern double d_Tp_2;
    extern double d_Te_2;

    extern long int x_flux;
    extern long int y_flux;
    extern long int z_flux;
    extern double p_x;
    extern double p_y;
    extern double p_z;
    extern int ee_transport_scattering_count_1;
    extern int ee_core_scattering_count_1;
    extern int ea_transport_scattering_count_1;
    extern int ea_core_scattering_count_1;
    extern int ee_transport_scattering_count_2;
    extern int ee_core_scattering_count_2;
    extern int ea_transport_scattering_count_2;
    extern int ea_core_scattering_count_2;

    extern double TTMe_1;
    extern double d_TTMe_1;
    extern double TTMp_1;
    extern double d_TTMp_1;
     extern double TTMe_2;
    extern double d_TTMe_2;
    extern double TTMp_2;
    extern double d_TTMp_2;
    extern double G;
    
    extern double transport_cutoff_1;
    extern double transport_cutoff_2;
    extern double core_cutoff;
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
    extern void initialize_cell_omp();
    extern void create();
    extern void output_data();

    extern void setup_output();
    extern void update_position();
    extern void update_dynamics();

//     extern void e_a_coulomb(const int& e, const int& array_index);
//               //  double& a_x_force, double& a_y_force, double& a_z_force, double& EPE, double& LPE);   
//    extern void neighbor_e_a_coulomb(const int& e, const int& array_index);
//                 // double& a_x_force, double& a_y_force, double& a_z_force, double& EPE, double& LPE);
    
    extern void e_e_coulomb(const int e, const int array_index);
    extern void neighbor_e_e_coulomb(const int e, const int array_index);
    
    // extern void a_a_coulomb(const int a, const int array_index, \
    //             double& a_x_force, double& a_y_force, double& a_z_force, double& LPE);
    // extern void neighbor_a_a_coulomb(const int a, const int array_index, \
    //             double& a_x_force, double& a_y_force, double& a_z_force, double& LPE);

    extern void electron_thermal_field(const int e, const int array_index, const double EKE, const int thread);
    extern double electron_applied_voltage(const int e, const int array_index, double& external_potential);
 
    extern void ea_scattering(const int e, const int array_index, const int thread);
    extern void ee_scattering();
    extern double M_B_distrib(const double& epsilon, const double& beta);
    extern double B_E_distrib(const double& epsilon);
    extern void create_phonon_distribution(std::vector<double>& distribution, const double& beta);
    extern void create_phonon_distribution(const std::string& name, std::vector<double>& distribution, const double& beta);
    extern void create_fermi_distribution(const std::string& name, std::vector<double>& distribution, const double beta);
   //  double return_phonon_distribution(const double& epsilon, const double& beta);
    
    extern double return_fermi_distribution(const double epsilon, const double beta);
    extern double return_BE_integrand(const double phonon_e, const double temperature);
}



#endif //CASTLE_H_
