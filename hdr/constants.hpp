//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <math.h>
//--------------------------------------------------------------------------------
// Namespace for variables and functions for exchange module
//--------------------------------------------------------------------------------
namespace constants{

   extern const double muB; // Bohr Magneton (Joules / Tesla)
   extern const double kB;  // Boltzmann constant (Joules / Kelvin)
   extern const double m_e; // electron mass (kg)
   extern const double h;   // Plank's constant (Js/kg)
   extern const double J_au; //Joule to a.u. 
   extern const double e;     // electron charge
   extern const double eps_0;
   extern const double K; //K for Coulomb's law
} // end of exchange namespace

#endif //CONSTANTS_H_
