

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
    extern  int  omp_threads;
    extern bool CASTLE_program;
    extern bool CASTLE_output_data; //output position and velocity data
    
    extern bool equilibrium_step;
    extern bool applied_voltage_sim;
    extern bool heat_pulse_sim;

    extern double lattice_atoms; //number of lattice atoms
    extern int conduction_electrons; //number of conduction electrons
    extern double temperature;

    extern int velocity_verlet_step(double time_step);
    extern double lattice_height; //Angstroms
    extern double lattice_width;  //Angsrtoms
    extern double lattice_depth;  //Angstroms
    extern double atomic_size;    //Angstroms
    extern double screening_depth;//Angstroms

    extern double min_as; //AJ
    extern double max_as; //AJ

    extern double x_unit_size;
    extern double y_unit_size;
    extern double z_unit_size;

    //simulation variables
    extern double total_time_steps;
    extern double loop_time;
    extern int    CASTLE_output_rate; 
    extern int CASTLE_MD_rate; //output velocity and position data at this multiple
  
    extern int full_int_var;
    extern std::vector<int> half_int_var;
    // [2];
    extern double boundary_conditions_cutoff;
    extern double dt;
    extern double v_f; //meters
    extern double E_f; //meters
    extern double E_f_A;
    extern double mu_f; //meters
    extern double n_f; //atoms/nm^3
    extern double atomic_mass;
    extern double mu_r; //inverse reduced mass in reduced units
    extern double combined_mass; //inverse with reduced units
    extern double total_e_scaling;
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
    extern double i_phonon_energy;
    extern int ee_density;
    extern int dos_size;
    extern double dos_occ;
    extern double local_dos_occ;
    extern std::vector<double> dos_standard; //e-/energy_hist_step (for full lattice)
    extern std::vector<double> dWdE_standard;
    extern std::vector<double> dWdE_standard_i;
    extern double step_size;

    extern bool ee_coupling;
    extern bool ea_coupling;
    extern double ea_rate;
    extern double ee_rate;
    extern double ee_coupling_strength;
    extern double ea_coupling_strength;
    extern double ee_scattering_angle;
    extern double q_sq; 
    extern double q_offset;
    extern double m_eff_ratio;

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
    extern double photon_energy;
    //integration variables
    extern long long int current_time_step;
    extern double CASTLE_real_time;
    extern int cells_per_thread;
    // extern std::vector<double> atom_anchor_position;
    // extern std::vector<double> atom_position;

    extern std::vector<double> electron_position; //Angstroms
    extern std::vector<double> electron_velocity; //Angstroms 
    extern std::vector<double> electron_potential; //A-1
    // extern std::vector<std::vector< int> > ee_dos_hist;
    extern std::vector<std::vector< int> > global_e_dos;

    // extern std::vector<bool> electron_transport_list;
    extern std::vector<std::vector< int> > electron_integration_list;
    // extern std::vector<std::vector< int> > electron_nearest_electron_list;
    // extern std::vector<std::vector< int> > electron_nearest_atom_list;
    extern std::vector<std::vector<double> > electron_ee_scattering_list;
    extern std::vector<std::vector< int> > electron_ea_scattering_list;
    extern std::vector<std::vector< int> > cell_lattice_coordinate;
    extern std::vector<std::vector< int> > cell_integration_lists;
    extern std::vector<std::vector< int> > old_cell_integration_lists;
    extern std::vector<std::vector< int> > cell_nearest_neighbor_list;
    extern std::vector<std::vector< int> > cell_lr_neighbor_list;
    extern std::vector<std::vector<std::vector< int> > > lattice_cell_coordinate;
    extern std::vector<std::vector< int> > temp_Map;
    extern std::vector<std::vector< int> > lattice_cells_per_omp;
    extern std::vector< int> escaping_electrons;
    extern std::vector<std::vector< int> > relaxation_time_hist_ee;
    // extern std::vector<std::vector< int> > relaxation_time_hist_ea;
    extern std::vector< double> flux_index;
    extern std::vector<MTRand> omp_uniform_random;
    extern std::vector<MTRand_int32> omp_int_random;
    //outputs
   
    extern double total_TEKE;
    extern double TEPE; //Angstroms
    extern double TEKE; //Angstroms
    extern double TLE; //Angstroms
    extern double Tp;
    extern double Te;
    extern double d_Tp;
    extern double d_Te;
    
    extern double MEPE; //meters
    extern double MEKE; //meters
    extern double MLE; //meters

    extern int x_flux;
    extern int y_flux;
    extern int z_flux;
    extern double p_x;
    extern double p_y;
    extern double p_z;
    extern  int e_a_scattering_count;
    extern  int e_e_scattering_count;
    extern  int ee_transport_scattering_count;
    extern  int ee_core_scattering_count;
    extern  int ea_transport_scattering_count;
    extern  int ea_core_scattering_count;

    extern double TTMe;
    extern double d_TTMe;
    extern double TTMp;
    extern double d_TTMp;
    extern double G;
    extern std::vector<double> global_tau_ep;
    extern std::vector<double> global_tau_ee;
    
    extern double transport_cutoff;
    extern double core_cutoff;
    extern double DoS_cutoff;
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

    extern void update_position();
    extern void update_dynamics();
    extern void e_e_coulomb(const int e, const int array_index);
    
    extern void neighbor_e_e_coulomb(const int e, const int array_index);
    extern void electron_thermal_field(const int e, const int array_index, const double EKE, const int thread);
 
    extern double electron_applied_voltage(const int e, const int array_index, const double potential);
 
    extern void ea_scattering(const int e, const int array_index, const int thread);
    extern void ee_scattering();
    extern void elastic_scattering(int thread, int e, int array_index, int d_e, int array_index_i, double e_energy, double d_e_energy );
    extern void inelastic_scattering(int thread, int e, int array_index, int d_e, int array_index_i, double e_energy, double d_e_energy);
    extern double B_E_distrib(const double& epsilon);
    extern void create_fermi_distribution(const std::string& name, std::vector<double>& distribution, const double temp);
    extern void create_defined_fermi_distribution(const std::string& name, std::vector<double>& distribution, const double temp);
  
    extern double return_fermi_distribution(const double energy, const double temp); //energy = e_i - E_f_A
    extern double return_BE_distribution(const double phonon_e, const double temperature);

    extern double return_dWdE(const double e_energy) ; //energy -> momentum 
   
    extern double return_dWdE_i(const double e_mom); // momentum -> energy
  
    extern double return_vel(const double energy); //inverse slope of dWdE / hbar_r

    extern double return_m_e_r(const double energy);
    extern double k_sq();
}



#endif //CASTLE_H_
