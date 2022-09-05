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

   //CASTLE Input
   bool applied_voltage_sim = false;
   bool heat_pulse_sim = false;
   double applied_voltage = 0.0;
   double heat_pulse = 0.0;

   double CASTLE_x_vector = 0.0;
	double CASTLE_y_vector = 0.0;
	double CASTLE_z_vector = 0.0;

   bool ee_coupling = false;
   bool ea_coupling = false;

   int CASTLE_MD_rate = 1;

   double ee_coupling_strength = 0.0;
   double ea_coupling_strength = 0.0;
   double ee_scattering_angle = 0.1;
   

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
   unsigned int  omp_threads = 1;
   bool CASTLE_program;
   bool CASTLE_output_data; //output position and velocity data
   bool equilibrium_step;
   bool applied_voltage_sim = false;
   bool heat_pulse_sim = false;

   unsigned int lattice_atoms; //number of lattice atoms
   int conduction_electrons; //number of conduction electrons
   double temperature;

   int velocity_verlet_step(double time_step);
   double lattice_height; //Angstroms
   double lattice_width;  //Angsrtoms
   double lattice_depth;  //Angstroms
   double atomic_size;    //Angstroms
   double screening_depth;//Angstroms

   double x_unit_size;
   double y_unit_size;
   double z_unit_size;

    //simulation variables
   double total_time_steps;
   double loop_time;
   unsigned int    CASTLE_output_rate; //output velocity and position data at this multiple
   unsigned int CASTLE_MD_rate;

   unsigned int full_int_var;
   unsigned int half_int_var;
   double boundary_conditions_cutoff;
   double dt;
   double v_f; //meters
   double E_f; //meters
   double E_f_A;
   double mu_f; //meters
   double n_f; //meters
   double atomic_mass;
   double mu_r; //inverse reduced mass in reduced units
   double combined_mass; //inverse with reduced units
   double Tr; // inverse seconds
   double a_heat_capacity;
   double e_heat_capacity;
   double a_heat_capacity_i;
   double e_heat_capacity_i;
   double a_specific_heat; //AJ/K/particle
   double e_specific_heat;
   double a_specific_heat_i;
   double e_specific_heat_i;
   double zero_pt_lattice_e;
   double phonon_energy;
   int ee_density;

   double e_e_integration_cutoff;
   double e_e_neighbor_cutoff;
   double e_e_coulomb_cutoff;

   double applied_voltage = 0.0;
   double power_density = 0.0;
   
   bool ee_coupling = false;
   bool ea_coupling = false;
   double ee_rate;
   double ea_rate;
   double ee_coupling_strength;
   double ea_coupling_strength;
   double ee_scattering_angle;

   unsigned int x_omp_cells;// = int(floor(lattice_width / 15.0));
   unsigned int y_omp_cells;// = int(floor(lattice_depth / 15.0));
   unsigned int z_omp_cells;// = int(floor(lattice_height/ 15.0));

  unsigned int total_cells;// = x_omp_cells*y_omp_cells*z_omp_cells;

  double x_step_size;// = lattice_width / double(x_omp_cells);
  double y_step_size;// = lattice_depth / double(y_omp_cells);
  double z_step_size;// = lattice_height/ double(z_omp_cells);

    //integration variables
    long long int current_time_step;
   double CASTLE_real_time;
   unsigned int cells_per_thread;
  // std::vector<double> atom_anchor_position;
  // std::vector<double> atom_position;

   std::vector<double> electron_position; //Angstroms
   std::vector<double> electron_velocity; //Angstroms
   std::vector<double> electron_potential; //A
   std::vector<bool> electron_transport_list;
   std::vector<std::vector<unsigned int> > electron_integration_list;
   std::vector<std::vector<unsigned int> > electron_nearest_electron_list;
   std::vector<std::vector<unsigned int> > electron_nearest_atom_list;
   std::vector<std::vector<uint32_t> > electron_ee_scattering_list;
   std::vector<std::vector<unsigned int> > electron_ea_scattering_list;
   std::vector<std::vector< int> > cell_lattice_coordinate;
   std::vector<std::vector<unsigned int > > cell_integration_lists;
   std::vector<std::vector<unsigned  int> > old_cell_integration_lists;
   std::vector<std::vector<unsigned int > > cell_nearest_neighbor_list;
   std::vector<std::vector < std::vector <unsigned int > > > lattice_cell_coordinate;
   std::vector<std::vector<unsigned int> > temp_Map;
   std::vector<std::vector<unsigned int> > lattice_cells_per_omp;
   std::vector<unsigned int> escaping_electrons;

//    MTRand_closed uniform_random;
//    MTRand_int32 int_random;
    std::vector<MTRand> omp_uniform_random(32);
    std::vector<MTRand_int32> omp_int_random(32);
//     std::vector<MTRand_int32> omp_int_random(32);
//   double = omp_int_random[omp_get_thread_num()]();
//    std::srand(std::time(nullptr));
//     std::random_device rd;  //Will be used to obtain a seed for the random number engine
//     std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
//     std::uniform_int_distribution<> test_int(0, conduction_electrons -1);
   //std::vector<MTRand> omp_gaussian_random(omp_get_num_threads());
    //outputs
   
   double TEPE; //Angstroms
   double TEKE; //Angstroms
   double TLE;
   double Tp;
   double Te;
    
   double MEPE; //meters
   double MEKE; //meters
   double MLE;

   long int x_flux;
   long int y_flux;
   long int z_flux;
   double p_x;
   double p_y;
   double p_z;
   unsigned int e_a_scattering_count;
   unsigned int e_e_scattering_count;
   unsigned int transport_scattering_count;
   unsigned int core_scattering_count;

   double TTMe;
   double d_TTMe;
   double TTMp;
   double d_TTMp;
   double G;
   double transport_cutoff;
   double core_cutoff;

   std::string time_stamp;
   std::ofstream lattice_output;
  
   std::ofstream electron_position_output_down;
   std::ofstream electron_velocity_output;
   std::ofstream mean_data;
   std::ofstream temp_data;
    
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
   void initialize_cell_omp();
   void create();
   void output_data();

   void setup_output();
   void update_position();
   void update_dynamics();

   // void e_a_coulomb(const int& e, const int& array_index);
   //            //  double& a_x_force, double& a_y_force, double& a_z_force, double& EPE, double& LPE);
   
   // void neighbor_e_a_coulomb(const int& e, const int& array_index);
   //              // double& a_x_force, double& a_y_force, double& a_z_force, double& EPE, double& LPE);
    
   void e_e_coulomb(const int e, const int array_index);
   void neighbor_e_e_coulomb(const int e, const int array_index);
   
   // void a_a_coulomb(const int a, const int array_index, \
   //              double& a_x_force, double& a_y_force, double& a_z_force, double& LPE);
   // void neighbor_a_a_coulomb(const int a, const int array_index, \
   //              double& a_x_force, double& a_y_force, double& a_z_force, double& LPE);

   void electron_thermal_field(const int e, const int array_index, const double EKE, const int thread);
 
   double electron_applied_voltage(const int e, const int array_index, double& external_potential);

   // void aa_scattering();
   void ea_scattering(const int e, const int array_index, const int thread);
   void ee_scattering();
    int ee_energy_conserved(const int electron, const int electron_collision, const double deltaE);
     int ee_final_momentum_conserved(const int electron, const int electron_collision, const double deltaE, const double e_energy, const double d_e_energy);
   int ee_elastic(const int electron, const int electron_collision,  const double e_energy, const double d_e_energy, const double probability);

   double B_E_distrib(const double& epsilon);
   double M_B_distrib(const double& epsilon, const double& beta);
   void create_phonon_distribution(std::vector<double>& distribution, const double& beta);
   void create_phonon_distribution(const std::string& name, std::vector<double>& distribution, const double& beta);
   void create_fermi_distribution(const std::string& name, std::vector<double>& distribution, const double& beta);
   double return_phonon_distribution(const double& epsilon, const double& beta);// {  return (1.0/(exp(epsilon/beta) + 1.0));}
}
