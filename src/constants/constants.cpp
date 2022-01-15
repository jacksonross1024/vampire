//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Vampire headers (need this for definition of external linkage of const variables)
#include "constants.hpp"

//--------------------------------------------------------------------------------
// Namespace for program constants
//--------------------------------------------------------------------------------
namespace constants{

   // fundamental constants
   const double muB = 9.27400999e-24; // Bohr Magneton (Joules / Tesla)
   const double kB  = 1.380649e-23;  // Boltzmann constant (Joules / Kelvin) NIST standard Reference
   const double m_e = 9.1093837015e-31; //electron mass (kg) NIST standard reference
   const double h   = 6.62607015e-34; // Plank constant (m^2 kg / s)  <or> (AJ /Hz) NIST standard reference
   const double e = 1.602176634e-19; //C
   const double esp_0 = 8.8541878128e-12; //permitivity of free space
  
   // derived constants
   const double J_au = 2.293712658357900e+17; //Joules in a.u.
   const double e_A = 16.02176634; //AC  (Angstrom Coulombs)
   const double m_e_r = 9.1093837015e-1; //kg reduced for future use in fs conversion
   const double m_e_r_i = 1 / m_e_r; //inverse reduced e- mass
   const double eps_0_A = 8.854187812800000000000000000; //reduced permittivity of free space (F/m) <or> (m^-3 kg^-1 s^2 C^2)
   const double K = e*e / (4*M_PI*M_PI *esp_0);
   const double K_A = 1e2*e_A*e_A / (4 * M_PI * M_PI * eps_0_A); // K for Coulomb's law in Angstroms
  
} // end of exchange namespace
