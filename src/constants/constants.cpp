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
   const double h   = 6.62607015e-34; // Plank constant (m^2 kg / s)  <or> (J /Hz) NIST standard reference
   const double J_au = 2.293712658357900e+17;
   const double e = 1.602176634e-19; //C 
   const double eps_0 = 8.8541878128e-12; //permittivity of free space, F/m
   const double K = e*e / (4 * M_PI * eps_0); // K for Coulomb's law
   // derived constants

} // end of exchange namespace
