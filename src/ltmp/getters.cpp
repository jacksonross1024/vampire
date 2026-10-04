//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------
//
// Getter functions for external modules (e.g., spintorque) to access
// local temperature pulse data.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <cmath>
#include <algorithm>

// Vampire headers
#include "ltmp.hpp"
#include "errors.hpp"
#include "vio.hpp"

// Local temperature pulse headers
#include "internal.hpp"

namespace ltmp {

//-----------------------------------------------------------------------------
// Temporal laser envelope G(t). Serban form if requested, else sc1d, else default.
//-----------------------------------------------------------------------------
namespace internal {
double temporal_laser_envelope(const double time_from_start){

   if(use_serban_laser_envelope){
      const double x = (time_from_start - pump_time) / pump_time;
      return std::exp(-8.0 * x * x);
   }

   if(use_sc1d_laser_params){
      const double sigma = sc1d_laser_fwhm / (2.0*std::sqrt(2.0*std::log(2.0)));
      const double x = (time_from_start - sc1d_laser_t0) / sigma;
      return std::exp(-0.5*x*x);
   }

   const double i_pump_time = 1.0 / pump_time;
   const double reduced_time = (time_from_start - 3.0 * pump_time) * i_pump_time;
   const double four_ln_2 = 2.77258872224;
   return std::exp(-four_ln_2 * reduced_time * reduced_time);
}
} // end of internal namespace

//-----------------------------------------------------------------------------
// Drive ltmp with the Serban laser envelope
//-----------------------------------------------------------------------------
void set_serban_laser_envelope(const bool enabled){
   ltmp::internal::use_serban_laser_envelope = enabled;
}

//-----------------------------------------------------------------------------
// Get electron temperature (K) for a given cell index
//-----------------------------------------------------------------------------
double get_electron_temperature(int cell) {
   if(!ltmp::internal::enabled || !ltmp::internal::initialised) {
      return 300.0; // default room temperature if not enabled
   }
   if(cell < 0 || cell >= ltmp::internal::num_cells) {
      return 300.0;
   }
   // root_temperature_array stores sqrt(Te), sqrt(Tp) pairs
   const double sqrt_Te = ltmp::internal::root_temperature_array[2*cell + 0];
   const double Te = sqrt_Te * sqrt_Te;

   if(err::check==true){
      if(!std::isfinite(Te) || Te < 0.0) {
         terminaltextcolor(RED);
         std::cerr << "Error: Invalid electron temperature detected in ltmp getter!" << std::endl;
         std::cerr << "  Cell: " << cell << std::endl;
         std::cerr << "  Te: " << Te << " K" << std::endl;
         std::cerr << "  sqrt_Te: " << sqrt_Te << std::endl;
         std::cerr << "  Simulation terminated." << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
   }

   return Te;
}

//-----------------------------------------------------------------------------
// Get phonon temperature (K) for a given cell index
//-----------------------------------------------------------------------------
double get_phonon_temperature(int cell) {
   if(!ltmp::internal::enabled || !ltmp::internal::initialised) {
      return 300.0;
   }
   if(cell < 0 || cell >= ltmp::internal::num_cells) {
      return 300.0;
   }
   const double sqrt_Tp = ltmp::internal::root_temperature_array[2*cell + 1];
   const double Tp = sqrt_Tp * sqrt_Tp;

   if(err::check==true){
      if(!std::isfinite(Tp) || Tp < 0.0) {
         terminaltextcolor(RED);
         std::cerr << "Error: Invalid phonon temperature detected in ltmp getter!" << std::endl;
         std::cerr << "  Cell: " << cell << std::endl;
         std::cerr << "  Tp: " << Tp << " K" << std::endl;
         std::cerr << "  sqrt_Tp: " << sqrt_Tp << std::endl;
         std::cerr << "  Simulation terminated." << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
   }

   return Tp;
}

//-----------------------------------------------------------------------------
// Get electron temperature (K) at position z (Angstroms)
// For vertical discretisation, find the cell containing this z position
//-----------------------------------------------------------------------------
double get_electron_temperature_at_z(double z_angstrom) {
   if(!ltmp::internal::enabled || !ltmp::internal::initialised) {
      return 300.0;
   }
   
   // Find the cell index for this z position
   if(ltmp::internal::vertical_discretisation) {
      const double dz = ltmp::internal::micro_cell_size[2];
      int cell = static_cast<int>(z_angstrom / dz);
      cell = std::max(0, std::min(cell, ltmp::internal::num_cells - 1));
      return get_electron_temperature(cell);
   }
   
   // If no vertical discretisation, return the single temperature
   return get_electron_temperature(0);
}

//-----------------------------------------------------------------------------
// Get phonon temperature (K) at position z (Angstroms)
//-----------------------------------------------------------------------------
double get_phonon_temperature_at_z(double z_angstrom) {
   if(!ltmp::internal::enabled || !ltmp::internal::initialised) {
      return 300.0;
   }
   
   if(ltmp::internal::vertical_discretisation) {
      const double dz = ltmp::internal::micro_cell_size[2];
      int cell = static_cast<int>(z_angstrom / dz);
      cell = std::max(0, std::min(cell, ltmp::internal::num_cells - 1));
      return get_phonon_temperature(cell);
   }
   
   return get_phonon_temperature(0);
}

//-----------------------------------------------------------------------------
// Get attenuated laser power density (W/m^3) at position z (Angstroms)
// Returns the instantaneous absorbed laser power at this depth
// Uses spin-currents laser parameters if set, otherwise uses ltmp parameters
//-----------------------------------------------------------------------------
double get_attenuated_laser_power(double z_angstrom, double time_from_start) {
   if(!ltmp::internal::enabled || !ltmp::internal::initialised) {
      return 0.0;
   }
   
   double gaussian = ltmp::internal::temporal_laser_envelope(time_from_start);
   double base_power_density = 0.0;
   
   if(!ltmp::internal::use_serban_laser_envelope && ltmp::internal::use_sc1d_laser_params) {
      base_power_density = ltmp::internal::sc1d_laser_Q0 * gaussian;
   } else {
      const double i_pump_time = 1.0 / ltmp::internal::pump_time;
      const double two_delta_sqrt_pi_ln_2 = 0.9394372787;
      base_power_density = 1e10 * ltmp::internal::pump_power * two_delta_sqrt_pi_ln_2 * gaussian * i_pump_time / ltmp::internal::penetration_depth;
   }
   
   // Find attenuation for this z position (always use ltmp's attenuation_array)
   if(ltmp::internal::vertical_discretisation) {
      const double dz = ltmp::internal::micro_cell_size[2];
      int cell = static_cast<int>(z_angstrom / dz);
      cell = std::max(0, std::min(cell, ltmp::internal::num_cells - 1));
      
      const double result = base_power_density * ltmp::internal::attenuation_array[cell];
      return result;
   }
   
   // If no vertical discretisation, use uniform attenuation
   return base_power_density;
}

//-----------------------------------------------------------------------------
// Peak of the envelope passed to the two-temperature step.
//-----------------------------------------------------------------------------
double laser_peak_time() {
   if(ltmp::internal::use_serban_laser_envelope) return ltmp::internal::pump_time;
   if(ltmp::internal::use_sc1d_laser_params) return ltmp::internal::sc1d_laser_t0;
   return 3.0 * ltmp::internal::pump_time;
}

//-----------------------------------------------------------------------------
// Set spin-currents laser parameters for use in ltmp heat flow
// Called when spin-currents laser is enabled but not using builtin laser
// optical_absorption_length (m) overrides ltmp's default penetration_depth
//-----------------------------------------------------------------------------
void set_spin_currents_laser_params(double Q0, double t0, double fwhm, double optical_absorption_length) {
   ltmp::internal::use_sc1d_laser_params = true;
   ltmp::internal::sc1d_laser_Q0 = Q0;
   ltmp::internal::sc1d_laser_t0 = t0;
   ltmp::internal::sc1d_laser_fwhm = fwhm;
   ltmp::internal::sc1d_optical_absorption_length = optical_absorption_length;
   
   // Override ltmp's penetration_depth with spin-currents optical absorption length
   // Convert from meters to Angstroms (ltmp uses Angstroms)
   if(optical_absorption_length > 0.0) {
      const double old_penetration_depth = ltmp::internal::penetration_depth;
      ltmp::internal::penetration_depth = optical_absorption_length * 1e10; // m -> Angstroms
      
      // Recalculate attenuation_array using the new penetration_depth
      // This ensures the spatial profile uses the spin-currents absorption length
      if(ltmp::internal::initialised && ltmp::internal::vertical_discretisation && 
         !ltmp::internal::cell_position_array.empty()) {
         
         // Find maximum z (top of stack)
         double system_dimensions_z = 0.0;
         for(unsigned int cell = 0; cell < ltmp::internal::num_cells; ++cell) {
            const double z = ltmp::internal::cell_position_array[3*cell + 2];
            if(z > system_dimensions_z) system_dimensions_z = z;
         }
         
         // Recalculate vertical attenuation for each cell
         for(unsigned int cell = 0; cell < ltmp::internal::attenuation_array.size(); ++cell) {
            const double z = ltmp::internal::cell_position_array[3*cell + 2];
            const double z_from_surface = system_dimensions_z - z;
            // Calculate new vertical attenuation with updated penetration_depth
            const double vattn_new = std::exp(-z_from_surface / ltmp::internal::penetration_depth);
            
            // Preserve lateral attenuation if it exists
            if(ltmp::internal::lateral_discretisation) {
               // Extract lateral component from current attenuation_array
               // Current value is vattn_old * lattn, so lattn = current / vattn_old
               const double z_from_surface_old = system_dimensions_z - z;
               const double vattn_old = std::exp(-z_from_surface_old / old_penetration_depth);
               if(vattn_old > 1e-10) {
                  const double lattn = ltmp::internal::attenuation_array[cell] / vattn_old;
                  ltmp::internal::attenuation_array[cell] = vattn_new * lattn;
               } else {
                  ltmp::internal::attenuation_array[cell] = vattn_new;
               }
            } else {
               ltmp::internal::attenuation_array[cell] = vattn_new;
            }
         }
      }
   }
}

//-----------------------------------------------------------------------------
// Get the current laser envelope value (0 to 1, gaussian temporal profile)
//-----------------------------------------------------------------------------
double get_laser_envelope(double time_from_start) {
   if(!ltmp::internal::enabled || !ltmp::internal::initialised) {
      return 0.0;
   }
   return ltmp::internal::temporal_laser_envelope(time_from_start);
}

//-----------------------------------------------------------------------------
// Get number of temperature cells
//-----------------------------------------------------------------------------
int get_num_cells() {
   if(!ltmp::internal::enabled || !ltmp::internal::initialised) {
      return 0;
   }
   return ltmp::internal::num_cells;
}

//-----------------------------------------------------------------------------
// Get cell z-position (Angstroms)
//-----------------------------------------------------------------------------
double get_cell_z_position(int cell) {
   if(!ltmp::internal::enabled || !ltmp::internal::initialised) {
      return 0.0;
   }
   if(cell < 0 || cell >= ltmp::internal::num_cells) {
      return 0.0;
   }
   return ltmp::internal::cell_position_array[3*cell + 2];
}

} // end of ltmp namespace
