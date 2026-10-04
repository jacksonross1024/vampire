//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <cmath>
#include <algorithm>

// Vampire headers
#include "spintorque.hpp"
#include "internal.hpp"
#include "atoms.hpp"
#include "material.hpp"
#include "random.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "ltmp.hpp"

namespace st {

//-----------------------------------------------------------------------------
// Get thermal fields from thermal gradients
//-----------------------------------------------------------------------------
void get_thermal_fields(std::vector<double>& thermal_x,
                        std::vector<double>& thermal_y,
                        std::vector<double>& thermal_z,
                        int start_idx, int end_idx) {
   if(!ltmp::is_enabled()) {
      for(int atom = start_idx; atom < end_idx; ++atom) {
         thermal_x[atom] = 0.0;
         thermal_y[atom] = 0.0;
         thermal_z[atom] = 0.0;
      }
      return;
   }

   for(int atom = start_idx; atom < end_idx; ++atom) {
      thermal_x[atom] = mtrandom::gaussian();
      thermal_y[atom] = mtrandom::gaussian();
      thermal_z[atom] = mtrandom::gaussian();
   }

   for(int atom = start_idx; atom < end_idx; ++atom) {
      const double zA = atoms::z_coord_array[atom];
      const int mat = atoms::type_array[atom];
      bool use_phonon = false;
      if(atom < (int)st::internal::sc1d_atom_use_phonon.size()){
         use_phonon = st::internal::sc1d_atom_use_phonon[atom];
      } else {
         use_phonon = mp::material[mat].couple_to_phonon_temperature;
      }
      const double T = use_phonon ? ltmp::get_phonon_temperature_at_z(zA)
                                  : ltmp::get_electron_temperature_at_z(zA);
      double rootT = (T > 0.0 && std::isfinite(T)) ? std::sqrt(T) : 0.0;

      const double H_th_sigma = mp::material[mat].H_th_sigma;
      const double alpha = mp::material[mat].temperature_rescaling_alpha;
      const double Tc = mp::material[mat].temperature_rescaling_Tc;
      if(Tc > 0.0) {
         const double Tloc = rootT * rootT;
         if(Tloc < Tc) {
            const double root_Tc = std::sqrt(Tc);
            rootT = root_Tc * std::pow(rootT / root_Tc, alpha);
         }
      }

      const double field_magnitude = H_th_sigma * rootT;
      thermal_x[atom] *= field_magnitude;
      thermal_y[atom] *= field_magnitude;
      thermal_z[atom] *= field_magnitude;
   }
}

//-----------------------------------------------------------------------------
// Check if thermal gradients are enabled and initialized
//-----------------------------------------------------------------------------
bool thermal_gradients_enabled() {
   return st::internal::sc1d_thermal_gradients_enable && ltmp::is_enabled();
}

} // end of namespace st
