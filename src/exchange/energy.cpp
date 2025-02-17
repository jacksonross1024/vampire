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

// C++ standard library headers

// Vampire headers
#include "atoms.hpp" // for exchange list type defs
#include "exchange.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

   //---------------------------------------------------------------------------
   // Calculate isotropic exchange energy for a single spin
   //---------------------------------------------------------------------------
   double spin_exchange_energy_isotropic(const int atom, const double sx, const double sy, const double sz){

   	// energy
   	double energy=0.0;

   	// Loop over neighbouring spins to calculate exchange
   	for(int nn = atoms::neighbour_list_start_index[atom]; nn <= atoms::neighbour_list_end_index[atom]; ++nn){

   		const int natom = atoms::neighbour_list_array[nn];
   		const double Jij = atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij;

         // note: sum over j only (not sum over i for j) leads to a silent factor 1/2 in exchange energy value
         //       - must be normalised in statistics to account for double sum
   		energy -= Jij * (atoms::x_spin_array[natom] * sx + atoms::y_spin_array[natom] * sy + atoms::z_spin_array[natom] * sz);

   	}

   	return energy;

   }

   //---------------------------------------------------------------------------
   // Calculate vectorial exchange energy for a single spin
   //---------------------------------------------------------------------------
   double spin_exchange_energy_vectorial(const int atom, const double sx, const double sy, const double sz){

   	// energy
   	double energy=0.0;

      // Loop over neighbouring spins to calculate exchange
   	for(int nn = atoms::neighbour_list_start_index[atom]; nn <= atoms::neighbour_list_end_index[atom]; ++nn){

   		const int natom = atoms::neighbour_list_array[nn];
   		const double Jij[3]={atoms::v_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[0],
   									atoms::v_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[1],
   									atoms::v_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[2]};

         // note: sum over j only (not sum over i for j) leads to a silent factor 1/2 in exchange energy value
         //       - must be normalised in statistics to account for double sum
   		energy -= ( Jij[0] * atoms::x_spin_array[natom] * sx +
                     Jij[1] * atoms::y_spin_array[natom] * sy +
                     Jij[2] * atoms::z_spin_array[natom] * sz);

   	}

   	return energy;

   }

   //---------------------------------------------------------------------------
   // Calculate tensorial exchange energy for a single spin
   //---------------------------------------------------------------------------
   double spin_exchange_energy_tensorial(const int atom, const double sx, const double sy, const double sz){

   	// energy
   	double energy=0.0;

      // Loop over neighbouring spins to calculate exchange
   	for(int nn = atoms::neighbour_list_start_index[atom]; nn <= atoms::neighbour_list_end_index[atom]; ++nn){

   		const int natom = atoms::neighbour_list_array[nn];
   		const double Jij[9]={atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[0],
   										 atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[1],
   										 atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[2],

   										atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[3],
   										 atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[4],
   										 atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[5],

   										atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[6],
   										 atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[7],
   										 atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[8]};

   		const double S[3]={atoms::x_spin_array[natom],atoms::y_spin_array[natom],atoms::z_spin_array[natom]};

         // note: sum over j only (not sum over i for j) leads to a silent factor 1/2 in exchange energy value
         //       - must be normalised in statistics to account for double sum
   		energy -= (Jij[0] * S[0] * sx + Jij[1] * S[1] * sx + Jij[2] * S[2] * sx +
   					  Jij[3] * S[0] * sy + Jij[4] * S[1] * sy + Jij[5] * S[2] * sy +
   					  Jij[6] * S[0] * sz + Jij[7] * S[1] * sz + Jij[8] * S[2] * sz);

   	}

   	return energy;

   }

   //---------------------------------------------------------------------------
   // Calculate  exchange energy for single spin selecting the correct type
   //---------------------------------------------------------------------------
   double single_spin_energy(const int atom, const double sx, const double sy, const double sz){

      // select calculation based on exchange type
      switch(internal::exchange_type){

   		case exchange::isotropic:
            return spin_exchange_energy_isotropic(atom, sx, sy, sz);
            break;


         case exchange::vectorial:
            return spin_exchange_energy_vectorial(atom, sx, sy, sz);
            break;


         case exchange::tensorial:
            return spin_exchange_energy_tensorial(atom, sx, sy, sz);
            break;

   	}

      return 0.0;

   }

} // end of exchange namespace
