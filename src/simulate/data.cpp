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

        bool CASTLE_program;
      bool CASTLE_output_data = false;
      bool equilibrium_step;
   int velocity_verlet_step(double dt);

   int CASTLE_output_rate;


     int lattice_atoms;
    double conduction_electrons;
    int core_electrons;

    double lattice_height;
    double lattice_width;
    double lattice_depth;

    double atomic_size;
    double screening_depth;
    double v_f;
    double mu_f;
    double E_f;

    double TKE;
    double TPE;
    int total_spin_up;
    int total_spin_down;
    std::vector<std::vector<bool> > symmetry_list;

   double temperature;

     double total_time_steps;
     double dt;
     int current_time_step;
     double loop_time;

   std::vector<double> electron_position; //superarray of past locations for each step
   std::vector<double> new_electron_position;
   std::vector<double> electron_velocity; //superarray for past velocities for each step
   std::vector<double> new_electron_velocity;
   std::vector<double> electron_acc;  
   std::vector<double> new_acc_array;
   std::vector<double> atom_position;
   std::vector<double> mean_data_array;
   std::vector<double> lattice_electrons;
   std::vector<bool> conduction_electron_spin;
   std::vector<bool> lattice_electron_spin;

   //std::vector<struct electron> electron_list; //This is probably not the best way to do this
    

   void initialize(const double num_electrons, const double num_atoms);

   void create();
   void create_gif();

    std::ofstream lattice_output;
    std::ofstream electron_position_output_up;
    std::ofstream electron_position_output_down;
    std::ofstream electron_velocity_output;
    std::ofstream mean_data;
    std::ofstream electron_spin_output;




}
