//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2022. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "program.hpp"

// program module headers
#include "internal.hpp"

namespace program{

   //---------------------------------------------------------------------------
   // Externally visible variables
   //---------------------------------------------------------------------------
   int program = 18; // program type to be run in vampire
   double fractional_electric_field_strength = 0.0; // factor controlling strength of stt/sot and voltage
   double laser_electrical_current = 0.0; // A/m^2
   double laser_electrical_S = 0.0;
   double laser_electrical_B = 0.0;

   namespace internal{

          int num_mag_cat;
		int num_mag_types;
		
		int num_dw_cells_x;
		int num_dw_cells_y;
		int num_dw_cells_z;

      int num_dw_cells;

      std::vector  < double > mag;
      std::vector <double > atom_to_cell_array;
		std::vector <int > cell_to_lattice_array;
      std::vector < int > num_atoms_in_cell;
      //------------------------------------------------------------------------
      // Shared variables inside program module
      //------------------------------------------------------------------------

      bool enabled = true; // bool to enable module

      //------------------------------------------------------------------------
      // Field pulse program
      //------------------------------------------------------------------------
      double field_pulse_time = 1.0e-12; // time constant for field pulse program (default 1 ps)

      //------------------------------------------------------------------------
      // Electrial pulse program
      //------------------------------------------------------------------------
      double electrical_pulse_time      = 0.0;//1.0e-9; // length of electrical pulses (1 ns default)
      double electrical_pulse_rise_time = 0.0;    // linear rise time for electrical pulse (0.0 default)
      double electrical_pulse_fall_time = 0.0;    // linear fall time for electrical pulse (0.0 default)
      int num_electrical_pulses         = 1;

      //------------------------------------------------------------------------
      // Laser electrical pulse (Serban two-filter ODE)
      // eta_t is the fraction of the excited sheet charge, e*F*lambda/(h*c) at
      // 800 nm, that crosses in the negative Jc lobe. 0.363 with F = 30 J/m^2
      // recovers the Fig. 1(a) fit A0 = 15.65 C/m^2 at tau_s = 69 fs, lambda = 0.78 ps.
      //------------------------------------------------------------------------
      bool laser_electrical_enable_temperature = true;
      double laser_electrical_eta_t = 0.363;
      double laser_electrical_tau_s = 69.0e-15;
      double laser_electrical_lambda = 0.78e-12;
      double laser_electrical_current_ref = 1.0e12; // 1 TA/m^2

      //------------------------------------------------------------------------
      // Material specific program parameters
      //------------------------------------------------------------------------
      std::vector<internal::mp_t> mp; // array of material properties

      double exchange_stiffness_min_constraint_angle = 0.0;
      double exchange_stiffness_max_constraint_angle   = 180.01; // degrees
      double exchange_stiffness_delta_constraint_angle =  5; // 22.5 degrees


   } // end of internal namespace

} // end of program namespace
