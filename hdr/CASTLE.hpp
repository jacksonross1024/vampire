



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

    // input variables
    extern unsigned int  omp_threads;
    extern bool CASTLE_program;
    extern bool CASTLE_output_data; //output position and velocity data
    
    extern bool equilibrium_step;
    extern bool applied_voltage_sim;
    extern bool heat_pulse_sim;

    extern unsigned int lattice_atoms; //number of lattice atoms
    extern unsigned int conduction_electrons; //number of conduction electrons
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
    extern unsigned int    CASTLE_output_rate; 
    extern unsigned int CASTLE_MD_rate; //output velocity and position data at this multiple
  
    extern unsigned int full_int_var;
    extern unsigned int half_int_var;
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

    extern unsigned int x_omp_cells;// = int(floor(lattice_width / 15.0));
    extern unsigned int y_omp_cells;// = int(floor(lattice_depth / 15.0));
    extern unsigned int z_omp_cells;// = int(floor(lattice_height/ 15.0));

    extern unsigned int total_cells;// = x_omp_cells*y_omp_cells*z_omp_cells;

    extern double x_step_size;// = lattice_width / double(x_omp_cells);
    extern double y_step_size;// = lattice_depth / double(y_omp_cells);
    extern double z_step_size;// = lattice_height/ double(z_omp_cells);
 
    extern double applied_voltage;
    extern double power_density;

    //integration variables
    extern long long int current_time_step;
    extern double CASTLE_real_time;
    extern unsigned int cells_per_thread;
    extern std::vector<double> atom_anchor_position;
    extern std::vector<double> atom_position;

    extern std::vector<double> electron_position; //Angstroms
    extern std::vector<double> electron_velocity; //Angstroms 
    extern std::vector<double> electron_potential; //A-1

    extern std::vector<bool> electron_transport_list;
    extern std::vector<std::vector<unsigned int> > electron_integration_list;
    extern std::vector<std::vector<unsigned int> > electron_nearest_electron_list;
    extern std::vector<std::vector<unsigned int> > electron_nearest_atom_list;
    extern std::vector<std::vector<uint32_t> > electron_ee_scattering_list;
    extern std::vector<std::vector<unsigned int> > electron_ea_scattering_list;
    extern std::vector<std::vector< int> > cell_lattice_coordinate;
    extern std::vector<std::vector<unsigned int> > cell_integration_lists;
    extern std::vector<std::vector<unsigned int> > old_cell_integration_lists;
    extern std::vector<std::vector< unsigned int > > cell_nearest_neighbor_list;
    extern std::vector<std::vector<std::vector<unsigned int> > > lattice_cell_coordinate;
    extern std::vector<std::vector<unsigned int> > temp_Map;
    extern std::vector<std::vector<unsigned int> > lattice_cells_per_omp;
    extern  std::vector<unsigned int> escaping_electrons;
    extern MTRand_closed uniform_random;
    extern MTRand_int32 int_random;
    extern std::vector<MTRand_closed> omp_uniform_random;
    extern std::vector<MTRand_int32> omp_int_random;
    //outputs
   
    extern double TEPE; //Angstroms
    extern double TEKE; //Angstroms
    extern double TLE; //Angstroms
    extern double Tp;
    extern double Te;
    
    extern double MEPE; //meters
    extern double MEKE; //meters
    extern double MLE; //meters

    extern long int x_flux;
    extern long int y_flux;
    extern long int z_flux;
    extern double p_x;
    extern double p_y;
    extern double p_z;
    extern unsigned int e_a_scattering_count;
    extern unsigned int e_e_scattering_count;
    extern unsigned int transport_scattering_count;
    extern unsigned int core_scattering_count;

    extern double TTMe;
    extern double d_TTMe;
    extern double TTMp;
    extern double d_TTMp;
    extern double G;
    
    extern double transport_cutoff;
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

    extern void electron_thermal_field(const int e, const int array_index, const double EKE);
 
    extern double electron_applied_voltage(const int e, const int array_index, double& external_potential);
 
    extern void ea_scattering(const int e, const int array_index);
    extern void ee_scattering();
    extern int ee_energy_conserved(const int electron, const int electron_collision, const double deltaE);
    extern int ee_final_momentum_conserved(const int electron, const int electron_collision, const double deltaE, const double e_energy, const double d_e_energy);
    extern int ee_elastic(const int electron, const int electron_collision,  const double e_energy, const double d_e_energy, const double probability);
    extern double M_B_distrib(const double& epsilon, const double& beta);
    extern double B_E_distrib(const double& epsilon);
    extern void create_phonon_distribution(std::vector<double>& distribution, const double& beta);
    extern void create_phonon_distribution(const std::string& name, std::vector<double>& distribution, const double& beta);
    extern void create_fermi_distribution(const std::string& name, std::vector<double>& distribution, const double& beta);
   //  double return_phonon_distribution(const double& epsilon, const double& beta);
    
    extern double return_phonon_distribution(const double epsilon, const double beta);
}



#endif //CASTLE_H_
