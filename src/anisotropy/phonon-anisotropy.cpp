//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: sw766@york.ac.uk
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
      // Function to add induced magnetic moment as anisotropy for spin-phonon parametrisation
      //
      // Induced magnetic moment occurs at frequency f_0 of the E_g phonon mode
      // with amplitude A
      //        this amplitude is a large collection of various factors depending on coherence of the excitation
      //        perhaps it changes with time, space, and the position of the stars
      //
      //    perhaps there is change in the length of the x length vector, or y vector
      //    but the y vector is always nonzero?
      //--------------------------------------------------------------------------------------------------------------
      

      //---------------------------------------------------------------------------------
      // Function to add second order uniaxial anisotropy
      // E = 2/3 * - ku2 (1/2)  * (3sz^2 - 1) == -ku2 sz^2 + const
      //---------------------------------------------------------------------------------
      double uniaxial_second_order_energy(const int atom,
                                          const int mat,
                                          const double sx,
                                          const double sy,
                                          const double sz){

         // get reduced anisotropy constant ku/mu_s (Tesla)
         const double ku2 = internal::ku2[mat];

         const double ex = internal::ku_vector[mat].x;
         const double ey = internal::ku_vector[mat].y;
         const double ez = internal::ku_vector[mat].z;

         const double sdote = (sx*ex + sy*ey + sz*ez);

         return -ku2*(sdote*sdote);

      }

   } // end of internal namespace

} // end of anisotropy namespace
