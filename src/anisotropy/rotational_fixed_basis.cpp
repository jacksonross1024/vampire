//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins and Richard Evans 2020. All rights reserved.
//
//   Email: sarah.jenkins@york.ac.uk richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

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
      // Function to add fourth order rotational anisotropy of the form
      //
      // E_4r = sin^3 theta cos (4 phi)
      //
      // In cartesian coordinates this expands to
      //
      // E_4r = 1 + Sz^4 - 8Sx^2 + 8Sx^2Sz^2 + 8Sx^4 - 2Sz^2
      //      = 1 + Sz^4 - 8Sy^2 + 8Sy^2Sz^2 + 8Sy^4 - 2Sz^2
      // E_4r = 1 - 8*Sx^2  + 8*Sx^4

      // The associated internal field (-dE/dS) is then
      //
      // Hx = 16 Sx (1 - Sz^2 - 2Sx^2)
      // Hy = 16 Sy (1 - Sz^2 - 2Sy^2)
      // Hz =  4 Sz (1 - Sz^2 - 4Sx^2)
      //--------------------------------------------------------------------------------------------------------------
      void rotational_fourth_order_fields_fixed_basis(std::vector<double>& spin_array_x,
                                                      std::vector<double>& spin_array_y,
                                                      std::vector<double>& spin_array_z,
                                                      std::vector<int>&    atom_material_array,
                                                      std::vector<double>& field_array_x,
                                                      std::vector<double>& field_array_y,
                                                      std::vector<double>& field_array_z,
                                                      const int start_index,
                                                      const int end_index){

         //if not enabled then do nothing
         if(!internal::enable_fourth_order_rotational) return;

         // Loop over all atoms between start and end index
         for(int atom = start_index; atom < end_index; atom++){

      //       // get atom material
      //       const int mat = atom_material_array[atom];
      //       const double sx = spin_array_x[atom]; // store spin direction in temporary variables
      //       const double sy = spin_array_y[atom];
      //       const double sz = spin_array_z[atom];
      //       // get reduced anisotropy constant ku/mu_s
      //       const double k4r = internal::k4r[mat];
      //       const double sx2 = sx*sx;
      //       const double sy2 = sy*sy;
      //       const double sz2 = sz*sz;

      //       field_array_x[atom] += k4r * 1.0 * sx * (1.0 - sz2 - 2.0 * sx2)*1;
      //       //field_array_x[atom] += k4r * 4 * sx*(sx2 - 3*sy2);
      //       field_array_y[atom] += k4r * 1.0 * sy * (1.0 - sz2 - 2.0 * sy2)*1;
      //       // field_array_y[atom] += k4r * 4 * sy*(sy2 - 3*sx2);
      //       field_array_z[atom] += k4r * 0.50 * sz * (1.0 - 2.0 * sz2 - 4.0 * sx2 - 4.0 * sy2)*1;
      //       // field_array_z[atom] += 0.0;
      //   }
      //    return;

         const double sx = spin_array_x[atom]; // store spin direction in temporary variables
            const double sy = spin_array_y[atom];
            const double sz = spin_array_z[atom];

            // get atom material
            const int mat = atom_material_array[atom];

            const double fx = sqrt(2.0)*0.5;//internal::kr_vector[mat].x;
            const double fy = sqrt(2.0)*0.5;//internal::kr_vector[mat].y;
            const double fz = 0.0;//internal::kr_vector[mat].z;

            const double gx = -sqrt(2.0)*0.5;//internal::kl_vector[mat].x;
            const double gy = sqrt(2.0)*0.5;// internal::kl_vector[mat].y;
            const double gz = 0.0;//internal::kl_vector[mat].z;

            // calculate S_x and S_x^3 parts
            const double Sx = sx * fx + sy * fy + sz * fz;
            const double Sx2 = Sx * Sx;

            // calculate S_y and S_y^3 parts
            const double Sy = sx * gx + sy * gy + sz * gz;
            const double Sy2 = Sy * Sy;

            // get reduced anisotropy constant ku/mu_s
            const double four_k4r4 = 4.0 * internal::k4r[mat];

            // calculate full form to add to field
            const double fullx = four_k4r4 * Sx * (Sx2 - 3 * Sy2);
            const double fully = four_k4r4 * Sy * (Sy2 - 3 * Sx2);

            field_array_x[atom] += fullx * fx;
            field_array_y[atom] += fullx * fy;
            field_array_z[atom] += fullx * fz;

            // sum y-component of field, where y-direction is represented by gx, gy, gz
            field_array_x[atom] += fully * gx;
            field_array_y[atom] += fully * gy;
            field_array_z[atom] += fully * gz;

         }

         return;

      }
      

      //---------------------------------------------------------------------------------
      // Function to add fourth order rotational anisotropy in spherical coord (cartesian has problems still)
      //   // E_4r = sin^3 theta cos (4 phi)
      //---------------------------------------------------------------------------------
      double rotational_fourth_order_energy_fixed_basis(
         const int atom,
         const int mat,
         const double sx,
         const double sy,
         const double sz){

         // get reduced anisotropy constant ku/mu_s
         // const double k4r = internal::k4r[mat];

         // //S already normalised to 1
         // double theta = atan2(sy,sx);
         //    if(theta != theta) theta = 0.0;
         // double phi = acos(sz);

         // const double energy = k4r*sin(theta)*sin(theta)*sin(theta)*sin(theta)*cos(4*phi);
         // //dE/dS_x = k4r*(-16 s_x +16 s_x*s_z^2 + 32 s_x^3)
         // //        =-k4r*8.0*s_x(1-s_z^2 - 2s_x^2)
         // //
         // //
         // return energy;
          const double k4r4 = internal::k4r[mat];

            const double fx = sqrt(2.0)*0.5;//internal::kr_vector[mat].x;
            const double fy = sqrt(2.0)*0.5;//internal::kr_vector[mat].y;
            const double fz = 0.0;//internal::kr_vector[mat].z;

            const double gx = -sqrt(2.0)*0.5;//internal::kl_vector[mat].x;
            const double gy = sqrt(2.0)*0.5;// internal::kl_vector[mat].y;
            const double gz = 0.0;//internal::kl_vector[mat].z;

         // calculate sin^4{theta}cos{4phi} = sin^4{theta} * ( 8 * cos^4{phi} - 8 * cos^2{phi} + 1 )
         //                                 = 8 * Sx^4 - 8 * sin^2{theta} * Sx^2 + sin^4{theta}
         //                                 = Sx^4 - 6 Sx^2 * Sy^2 + Sy^4
         const double Sx = sx * fx + sy * fy + sz * fz;
         const double Sx2 = Sx * Sx;

         const double Sy = sx * gx + sy * gy + sz * gz;
         const double Sy2 = Sy * Sy;

         return - k4r4 * (Sx2 * Sx2 - 6.0 * Sx2 * Sy2 + Sy2 * Sy2 );

      }

   } // end of internal namespace

} // end of anisotropy namespace
