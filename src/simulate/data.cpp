//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "sim.hpp"

// Internal sim header
#include "internal.hpp"
#include "CASTLE.hpp"

namespace sim{
   //----------------------------------------------------------------------------
   // Shared variables used with main vampire code
   //---------------------------------------------------------------------------
   integrator_t integrator = llg_heun; // variable to specify integrator


   std::vector < double > track_field_x;
   std::vector < double > track_field_y;
   std::vector < double > track_field_z;

   double track_Ms = 0.1;
   double track_bit_width = 800;
   double track_bit_size = 1000;
   double track_bit_depth = 600;
   double cross_track_velocity = 0.0;
   double down_track_velocity = 0.0;
   double track_pos_x;
   double track_pos_z;
   double Ms;


   double initial_down_track_position = 0.0;
   double initial_cross_track_position = 0.0;

   int track_num_bits_per_track = 1;
   int track_num_tracks = 1;

   double track_bit_gap = 10.0;
	 double track_track_gap = 10.0;

   // distance of tracks from read head
   double track_fly_height = 100.0; // Angstroms

   double LFA_scan_field_step = 0.01;
   bool LFA = false;

   bool track_ms_file = false;

   int num_monte_carlo_preconditioning_steps(0);

   uint64_t time         = 0; // time step counter
   uint64_t total_time   = 10000; // total time steps (non-loop code)
   uint64_t loop_time    = 10000; // loop time steps (hysteresis/temperature loops)
   uint64_t partial_time = 1000; // same as time-step-increment
   uint64_t equilibration_time = 0; // equilibration time steps

   bool calculate_program_convergence           = false;
   double convergence_criteria;
   unsigned int convergence_check;
   bool output_convergence_counter              = false;

   int domain_wall_axis = 0;
   double domain_wall_position = 0.25;
   double domain_wall_discretisation = 10;
   double domain_wall_centre = 0;
   double domain_wall_width = 10.0;
   std::vector <bool > anti_PBC(3,false);
   std::vector < double > domain_wall_second_vector_x(100,0);
   std::vector < double > domain_wall_second_vector_y(100,0);
   std::vector < double > domain_wall_second_vector_z(100,1.0);


   //Fermi Gas
   bool calculate_fermi_distribution;
   std::vector<long double> fermi_distribution_array;
   double temperature_variable;
   long double beta_variable;
   double fermi_size;
   long double eps_variable;
   long double mu;

   long double mu_0; //mu at T = 0
	long double E_0; //E at T = 0
   double P_0;
	double fermi_electrons; //number of conduction electrons; perfect gas = num_atoms
   double conduction_electrons; //conduction electrons per atom. perfect gas = 1
	double fermi_volume; //volume to calculate the electron density
   double eta;     // = 1 / Beta * mu_0
   double fermi_density;

   bool fermi_pressure = false;
   bool fermi_energy = false;
   bool fermi_Cv = false;
   bool fermi_function= false;

  
   namespace internal{

      //----------------------------------------------------------------------------
      // Shared variables used within sim module
      //----------------------------------------------------------------------------
      bool enable_spin_torque_fields = false; // flag to enable spin torque fields

      std::vector<sim::internal::mp_t> mp; // array of material properties

      std::vector<double> stt_asm; // array of spin transfer torque asymmetry
      std::vector<double> stt_rj; // array of adiabatic spin torques
      std::vector<double> stt_pj; // array of non-adiabatic spin torques
      std::vector<double> stt_polarization_unit_vector(3,0.0); // stt spin polarization direction

      std::vector<double> sot_asm; // array of spin orbit torque asymmetry
      std::vector<double> sot_rj; // array of adiabatic spin torques
      std::vector<double> sot_pj; // array of non-adiabatic spin torques
      std::vector<double> sot_polarization_unit_vector(3,0.0); // sot spin polarization direction

   } // end of internal namespace

   //------------------------------------------------------------------------
   // getter functions to give access to internal variables
   //------------------------------------------------------------------------
   std::vector<double> get_stt_polarization_unit_vector(){
      return sim::internal::stt_polarization_unit_vector;
   }

   std::vector<double> get_stt_rj(){
      return sim::internal::stt_rj;
   }

   std::vector<double> get_stt_pj(){
      return sim::internal::stt_pj;
   }

} // end of sim namespace

namespace CASTLE {

     // input and material parameters
   bool CASTLE_program;
   bool CASTLE_output_data; //output position and velocity data
   bool equilibrium_step;

   int lattice_atoms; //number of lattice atoms
   int conduction_electrons; //number of conduction electrons
   double temperature;

   int velocity_verlet_step(double time_step);
   double lattice_height; //Angstroms
   double lattice_width;  //Angsrtoms
   double lattice_depth;  //Angstroms
   double atomic_size;    //Angstroms
   double screening_depth;//Angstroms

    //simulation variables
   double total_time_steps;
   double loop_time;
   int    CASTLE_output_rate; //output velocity and position data at this multiple
   long double chosen_electron;
   long double dt;
   long double v_f; //meters
   long double E_f; //meters
   long double E_f_A;
   long double mu_f; //meters
   long double n_f; //meters
   long double atomic_mass;


   long double e_a_neighbor_cutoff;
   long double e_a_coulomb_cutoff;
   long double e_e_neighbor_cutoff;
   long double e_e_coulomb_cutoff;
   long double a_a_neighbor_cutoff;
   long double a_a_coulomb_cutoff;

    //integration variables
   int current_time_step;
   double CASTLE_real_time;
   
   std::vector<double> atom_position;
   std::vector<double> new_atom_position;
   std::vector<long double> atom_velocity; //Angstroms
   std::vector<long double> new_atom_velocity;
   std::vector<long double> atom_force;   //Angstroms
   std::vector<long double> new_atom_force;
   std::vector<std::vector<int> > atomic_nearest_electron_list;
   std::vector<std::vector<int> > atomic_nearest_atom_list;

   std::vector<long double> electron_position; //Angstroms
   std::vector<long double> new_electron_position;
   std::vector<long double> electron_velocity; //Angstroms
   std::vector<long double> new_electron_velocity;
   std::vector<long double> electron_force;   //Angstroms
   std::vector<long double> new_electron_force;
   std::vector<long double> electron_potential; //A-1
   std::vector<long double> new_electron_potential;
   std::vector<std::vector<int> > electron_nearest_electron_list;
   std::vector<std::vector<int> > electron_nearest_atom_list;

    //outputs
   long double TEPE; //Angstroms
   long double TEKE; //Angstroms
   long double TLPE; //Angstroms
   long double TLKE; //Angstroms
    
   long double MEPE; //meters
   long double MEKE; //meters
   long double MLKE; //meters
   long double MLPE; //meters

   int x_flux;
   int y_flux;
   int z_flux;
   std::string time_stamp;
   std::ofstream lattice_output;
  
   std::ofstream electron_position_output_down;
   std::ofstream electron_velocity_output;
   std::ofstream mean_data;
    
    //control functions
   void initialize();
   void initialize_positions();
   void initialize_lattice();
   void initialize_electrons();
   void initialize_atoms();
   void initialize_forces();
   void initialize_electron_interactions();
   void initialize_atomic_interactions();
   void initialize_electron_atom_interactions();
   void initialize_velocities();
   void create();
   void output_data();

   void setup_output();
   void update_position();
   void update_dynamics();

   void e_a_coulomb(const int e, const int array_index, long double& e_x_force, long double& e_y_force, long double& e_z_force, \
                    long double& a_x_force, long double& a_y_force, long double& a_z_force, long double& EPE, long double& LPE);
   void neighbor_e_a_coulomb(const int e, const int array_index, long double& e_x_force, long double& e_y_force, long double& e_z_force, \
                             long double& a_x_force, long double& a_y_force, long double& a_z_force, long double& EPE, long double& LPE);
    
   void e_e_coulomb(const int e, const int array_index, long double& e_x_force, long double& e_y_force, long double& e_z_force, \
                    long double& EPE);
   void neighbor_e_e_coulomb(const int e, const int array_index, long double& e_x_force, long double& e_y_force, long double& e_z_force, \
                             long double& EPE);
    
   void a_a_coulomb(const int a, const int array_index, \
                   long double& a_x_force, long double& a_y_force, long double& a_z_force, long double& LPE);
   void neighbor_a_a_coulomb(const int a, const int array_index, \
                            long double& a_x_force, long double& a_y_force, long double& a_z_force, long double& LPE);

   void update_velocity(int array_index, long double& EKE, long double& LKE);

/*
   long double e_a_scattering(int e, int a, const long double& l_x, const long double& l_y, const long double& l_z);
   long double e_p_scattering(int e, int a, const long double& x_distance, const long double& y_distance, const long double& z_distance);
   long double electron_applied_voltage(int array_index, long double& x_force, long double& y_force, long double& z_force);
   long double reinitialize_electron_conserve_momentum(std::vector<long double>& captured_electron_list);

*/

}
