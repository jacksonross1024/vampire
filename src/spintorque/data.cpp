//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "spintorque.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace st{
   namespace internal{
      //-----------------------------------------------------------------------------
      // Shared variables used for the spin torque calculation
      //-----------------------------------------------------------------------------
      bool enabled=false;  // disable spin torque calculation
      bool TMRenable=false; // disable TMR calculation, GMR is set as dafault
      bool fbc = true;
      //double micro_cell_size= 5*3.00; /// lateral size of spin torque microcells
      std::vector<double> micro_cell_size(3);
      std::string microcell_decomp_type;
      double micro_cell_thickness = 3.00; /// thickness of spin torque microcells (atomistic)

      int num_local_atoms; /// number of local atoms (ignores halo atoms in parallel simulation)
      int remove_nm = -1;
      int current_direction =2; /// direction for current x->0, y->1, z->2
      //   std::vector< std::vector< micro_cell_t > > stack;
      std::vector<int> atom_st_index; // mc which atom belongs to
      std::vector<int> cell_stack_index; 
      std::vector<double> x_field_array; // arrays to store atomic spin torque field
      std::vector<double> y_field_array;
      std::vector<double> z_field_array;

      int stx=0; // indices for x,y,z in the spin torque coordinate system (default z)
      int sty=1;
      int stz=2;

      int num_stacks_x;  // total number of stacks
      int num_stacks_y;  // total number of stacks
      int num_x_stacks; // number of stacks in x
      int num_y_stacks; // number of stack in y
      int num_microcells_per_stack; // number of microcells per stack

      int config_file_counter = 0; // spin torque config file counter

      int free_layer = 0;       /// index of free layer in magnetic tunnel junction
      int reference_layer = 0;  /// index of reference layer in magnetic tunnel junction

      double je = 0.0; // = 1.0e11; // current (C/s/m^2)
      double initial_beta;
      double initial_theta = 0.2;
      
      double rel_angle;  // relative angle between 2 FMs for TMR calculation
      
      int ST_output_rate;
      std::string output_torque_data;
      std::vector<double> initial_m(3);
     
      std::vector<double> init_stack_mag;
      std::vector<double> stack_init_mag;
      std::vector<int> stack_index_y; // start of stack in microcell arrays
      std::vector<int> stack_index_x; // start of stack in microcell arrays
      std::vector<std::vector<int> > cell_index_x;

      //STT
      std::vector<double> beta_cond; /// spin polarisation (conductivity) Beta B
      std::vector<double> beta_diff; /// spin polarisation (diffusion) Beta' Bp
      std::vector<double> sa_infinity; /// intrinsic spin accumulation
      std::vector<double> lambda_sdl; /// spin diffusion length
      std::vector<double> diffusion; /// diffusion constant Do
      std::vector<double> sd_exchange; /// electron(s)-spin(d) exchange interaction
      std::vector<double> a; // a parameter for spin accumulation
      std::vector<double> b; // b parameter for spin accumulation

      //SOT
      std::vector<double> sot_beta_cond; /// spin polarisation (conductivity) Beta B
      std::vector<double> sot_beta_diff; /// spin polarisation (diffusion) Beta' Bp
      std::vector<double> sot_sa_infinity; /// intrinsic spin accumulation
      std::vector<double> sot_lambda_sdl; /// spin diffusion length
      std::vector<double> sot_diffusion; /// diffusion constant Do
      std::vector<double> sot_sd_exchange; /// electron(s)-spin(d) exchange interaction
      std::vector<double> sot_a; // a parameter for spin accumulation
      std::vector<double> sot_b; // b parameter for spin accumulation
      std::vector<double> spin_acc_sign;
      std::vector<bool> sot_sa_source;
      bool sot_sa = false;
      //sot sa parameters 

      std::vector<double> coeff_ast; // adiabatic spin torque
      std::vector<double> coeff_nast; // non-adiabatic spin torque
      std::vector<double> cell_natom;
      
      // three-vector arrays
      std::vector<double> pos; /// stack position
      std::vector<double> m; // magnetisation
      std::vector<double> j_final_up_y;
      std::vector<double> j_final_down_y;
      std::vector<double> j_int_up_y; // spin current
      std::vector<double> j_int_down_y;
      std::vector<double> j_init_up_y; // spin current
      std::vector<double> j_init_down_y;
      std::vector<double> j_final_up_x;
      std::vector<double> sa_final; // spin accumulation
      // std::vector<double> sa_sot_final; // spin accumulation
      std::vector<double> sa_int; // spin accumulation
      // std::vector<double> sa_sot_init;
      std::vector<double> spin_torque; // spin torque energy (J)
      std::vector<double> ast; // adiabatic spin torque
      std::vector<double> nast; // non-adiabatic spin torque
      std::vector<double> total_ST; // non-adiabatic spin torque
      std::vector<double> magx_mat; // magnetisation of material
      std::vector<double> magy_mat;
      std::vector<double> magz_mat;
      
      std::vector<double> sa_sum;
      std::vector<double> j_final_up_y_sum;
      std::vector<double> j_final_down_y_sum;
      std::vector<double> j_final_up_x_sum;
      bool sot_check = false;
      
      std::vector<double> coeff_ast_sum;
      std::vector<double> coeff_nast_sum;
      std::vector<double> ast_sum;
      std::vector<double> nast_sum;
      std::vector<double> total_ST_sum;
      std::vector<int> cell_natom_sum;
      std::vector<int> mpi_stack_list_y;
      std::vector<int> mpi_stack_list_x;

      // array of material properties
      std::vector<st::internal::mp_t> mp;

      // default material properties
      st::internal::mp_t default_properties;

         // stopwatch.start();
   } // end of internal namespace
   
      double spin_acc_time = 0.0;
} // end of st namespace
