//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard F L Evans 2020. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <fstream>
#include <cmath>

// Vampire headers
#include "anisotropy.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------

   namespace internal{

      //---------------------------------------------------------------------------------
      // Function to add fourth order uniaxial anisotropy along vector e
      //
      //  Higher order anisotropies generally need to be described using orthogonal
      //  functions. The usual form (a series in S leads to cross pollution of terms,
      //  giving strange temperature dependencies.
      //
      //  The anisotropies are described with a minimal orthogonal set expansion,
      //  preserving the orthogonality of different orders while being simple to
      //  implement and understand. Explicity the energies are described by normalising
      //  the inner summation of the 2,4,6 order spherical harmonics to the prefactor
      //  of the highest order term with an abritrary shift so that E(0) = 0.
      //
      //  k2(sz) = k2 (1 - sz^2)
      //  k4(sz) = k4 (sz^4 - 30sz^2 / 35 - 5/35)
      //  k6(sz) = k6 (sz^6 - 315*sz^4/231 + 105sz^2/231 - 5/231)
      //
      //  The field induced by the harmonics is given by the first derivative w.r.t. sz.
      //  This can be projected onto any arbritrary direction ex,ey,ez allowing higher
      //  order anisotropy terms along any direction. This direction is shared with the
      //  other uniaxial anisotropy coefficients since they should not be used
      //  simultaneously.
      //
      //--------------------------------------------------------------------------------------------------------------
      void biaxial_fourth_order_simple_fields(std::vector<double>& spin_array_x,
                                        std::vector<double>& spin_array_y,
                                        std::vector<double>& spin_array_z,
                                        std::vector<int>&    atom_material_array,
                                        std::vector<double>& field_array_x,
                                        std::vector<double>& field_array_y,
                                        std::vector<double>& field_array_z,
                                        const int start_index,
                                        const int end_index){

         // if not enabled then do nothing
         if(!internal::enable_biaxial_fourth_order_simple) return;

         // Loop over all atoms between start and end index
         for(int atom = start_index; atom < end_index; atom++) {

            // get atom material
            const int mat = atom_material_array[atom];

            const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];

            const double u1x = internal::mp[mat].u1_vector[0];
            const double u1y = internal::mp[mat].u1_vector[1];
            const double u1z = internal::mp[mat].u1_vector[2];

            const double u2x = internal::mp[mat].u2_vector[0];
            const double u2y = internal::mp[mat].u2_vector[1];
            const double u2z = internal::mp[mat].u2_vector[2];

            const double u3x = -1.0*internal::mp[mat].u1_vector[0];
            const double u3y = -1.0*internal::mp[mat].u1_vector[1];
            const double u3z = -1.0*internal::mp[mat].u1_vector[2];

            const double u4x = -1.0*internal::mp[mat].u2_vector[0];
            const double u4y = -1.0*internal::mp[mat].u2_vector[1];
            const double u4z = -1.0*internal::mp[mat].u2_vector[2];

            // get reduced anisotropy constant ku/mu_s
            const double ku4 = internal::ku4[mat];

            const double sdotu1 = (sx*u1x + sy*u1y + sz*u1z);
            const double sdotu13 = sdotu1*sdotu1*sdotu1;

            const double sdotu2 = (sx*u2x + sy*u2y + sz*u2z);
            const double sdotu23 = sdotu2*sdotu2*sdotu2;

            const double sdotu3 = (sx*u3x + sy*u3y + sz*u3z);
            const double sdotu33 = sdotu3*sdotu3*sdotu3;

            const double sdotu4 = (sx*u4x + sy*u4y + sz*u4z);
            const double sdotu43 = sdotu4*sdotu4*sdotu4;

            // calculate field (double negative from scale factor and negative derivative)

            field_array_x[atom] += (2*ku4)*(u1x*sdotu13+u2x*sdotu23+u3x*sdotu33+u4x*sdotu43);
            field_array_y[atom] += (2*ku4)*(u1y*sdotu13+u2y*sdotu23+u3y*sdotu33+u4y*sdotu43);
            field_array_z[atom] += (2*ku4)*(u1z*sdotu13+u2z*sdotu23+u3z*sdotu33+u4z*sdotu43);

            // field_array_x[atom] += (ku4)*(sx*sx*sx + 4.0*sy*sy*sx);
            // field_array_y[atom] += (2*ku4)*(u1y*sdotu13+u2y*sdotu23+u3y*sdotu33+u4y*sdotu43);
            // field_array_z[atom] += (2*ku4)*(u1z*sdotu13+u2z*sdotu23+u3z*sdotu33+u4z*sdotu43);



         }

         return;

      }

      //---------------------------------------------------------------------------------
      // Function to add fourth order uniaxial anisotropy
      //---------------------------------------------------------------------------------
      // file scope constant

      double biaxial_fourth_order_simple_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz){

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double ku4 = internal::ku4[mat];



            const double u1x = internal::mp[mat].u1_vector[0];
            const double u1y = internal::mp[mat].u1_vector[1];
            const double u1z = internal::mp[mat].u1_vector[2];

            const double u2x = internal::mp[mat].u2_vector[0];
            const double u2y = internal::mp[mat].u2_vector[1];
            const double u2z = internal::mp[mat].u2_vector[2];

         const double sdotu1 = (sx*u1x + sy*u1y + sz*u1z);
         const double sdotu14 = sdotu1*sdotu1*sdotu1*sdotu1;

         const double sdotu2 = (sx*u2x + sy*u2y + sz*u2z);
         const double sdotu24 = sdotu2*sdotu2*sdotu2*sdotu2;

         return -(ku4/2.0)*(sdotu14+sdotu24);

      }

   } // end of internal namespace

} // end of anisotropy namespace
