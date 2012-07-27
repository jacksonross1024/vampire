///
/// @file
/// @brief Contains functions to calculate energy for a spin/system
///
/// @details None
///
/// @section notes Implementation Notes
/// This is a list of other notes, not related to functionality but rather to implementation. 
/// Also include references, formulae and other notes here.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///

// Standard Libraries
#include <algorithm>
#include <cmath>
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "material.hpp"
#include "errors.hpp"
#include "demag.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

namespace sim{

/// @brief Calculates the exchange energy for a single spin (isotropic).
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @param[in] int atom number 
/// @param[in] imaterial material of local atom
/// @param[in] Sx x-spin of local atom  
/// @param[in] Sy y-spin of local atom 
/// @param[in] Sz z-spin of local atom 
/// @return exchange energy
///
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///
inline double spin_exchange_energy_isotropic(const int atom, const double Sx, const double Sy, const double Sz){
	
	// energy
	double energy=0.0;
	
	// Loop over neighbouring spins to calculate exchange
	for(int nn=atoms::neighbour_list_start_index[atom];nn<=atoms::neighbour_list_end_index[atom];nn++){
			
		const int natom = atoms::neighbour_list_array[nn];
		const double Jij=atoms::i_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij;

		energy+=Jij*(atoms::x_spin_array[natom]*Sx + atoms::y_spin_array[natom]*Sy + atoms::z_spin_array[natom]*Sz);
	}
		
	return energy;
	
}

/// @brief Calculates the exchange energy for a single spin (vector).
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @param[in] atom atom number 
/// @param[in] Sx x-spin of local atom  
/// @param[in] Sy y-spin of local atom 
/// @param[in] Sz z-spin of local atom 
/// @return exchange energy
///
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///
inline double spin_exchange_energy_vector(const int atom, const double Sx, const double Sy, const double Sz){
	
	// energy
	double energy=0.0;
	
	// Loop over neighbouring spins to calculate exchange
	for(int nn=atoms::neighbour_list_start_index[atom];nn<=atoms::neighbour_list_end_index[atom];nn++){
			
		const int natom = atoms::neighbour_list_array[nn];
		const double Jij[3]={atoms::v_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[0],
									atoms::v_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[1],
									atoms::v_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[2]};

		energy+=(Jij[0]*atoms::x_spin_array[natom]*Sx + Jij[1]*atoms::y_spin_array[natom]*Sy + Jij[2]*atoms::z_spin_array[natom]*Sz);
	}
		
	return energy;
	
}

/// @brief Calculates the exchange energy for a single spin (tensor).
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @param[in] atom atom number 
/// @param[in] Sx x-spin of local atom  
/// @param[in] Sy y-spin of local atom 
/// @param[in] Sz z-spin of local atom 
/// @return exchange energy
///
/// @internal
///	Created:		27/07/2012
///	Revision:	  ---
///=====================================================================================
///
inline double spin_exchange_energy_tensor(const int atom, const double Sx, const double Sy, const double Sz){
	
	// energy
	double energy=0.0;
	
	// Loop over neighbouring spins to calculate exchange
	for(int nn=atoms::neighbour_list_start_index[atom];nn<=atoms::neighbour_list_end_index[atom];nn++){
			
		const int natom = atoms::neighbour_list_array[nn];
		const double Jij[3][3]={atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[0][0],
										atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[0][1],
										atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[0][2],

										atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[1][0],
										atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[1][1],
										atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[1][2],

										atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[2][0],
										atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[2][1],
										atoms::t_exchange_list[atoms::neighbour_interaction_type_array[nn]].Jij[2][2],};
				
		const double S[3]={atoms::x_spin_array[natom],atoms::y_spin_array[natom],atoms::z_spin_array[natom]};
		
		energy+=(Jij[0][0]*S[0]*Sx + Jij[0][1]*S[1]*Sx +Jij[0][2]*S[2]*Sx +
					Jij[1][0]*S[0]*Sy + Jij[1][1]*S[1]*Sy +Jij[1][2]*S[2]*Sy +
					Jij[2][0]*S[0]*Sz + Jij[2][1]*S[1]*Sz +Jij[2][2]*S[2]*Sz);

	}
		
	return energy;
	
}

/// @brief Calculates the uniaxial anisotropy energy for a single spin.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @param[in] atom atom number 
/// @param[in] imaterial material of local atom
/// @param[in] Sx x-spin of local atom  
/// @param[in] Sy y-spin of local atom 
/// @param[in] Sz z-spin of local atom 
/// @return exchange energy
///
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///
inline double spin_uniaxial_energy(const int atom, const int imaterial, const double Sx, const double Sy, const double Sz){
	
	return mp::material[imaterial].Ku*Sz*Sz;
	
}

/// @brief Calculates the applied field energy for a single spin.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @param[in] atom atom number 
/// @param[in] imaterial material of local atom
/// @param[in] Sx x-spin of local atom  
/// @param[in] Sy y-spin of local atom 
/// @param[in] Sz z-spin of local atom 
/// @return applied field energy
///
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///
inline double spin_applied_field_energy(const int atom, const int imaterial, const double Sx, const double Sy, const double Sz){;

	return -sim::H_applied*(sim::H_vec[0]*Sx + sim::H_vec[1]*Sy + sim::H_vec[2]*Sz);

}

/// @brief Calculates the surface anisotropy energy for a single spin.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @param[in] atom atom number 
/// @param[in] imaterial material of local atom
/// @param[in] Sx x-spin of local atom  
/// @param[in] Sy y-spin of local atom 
/// @param[in] Sz z-spin of local atom 
/// @return exchange energy
///
/// @internal
///	Created:		13/09/2011
///	Revision:	  ---
///=====================================================================================
///
inline double spin_surface_anisotropy_energy(const int atom, const int imaterial, const double Sx, const double Sy, const double Sz){
	
	double energy=0.0;

	if(atoms::surface_array[atom]==true && sim::surface_anisotropy==true){
		const double Ks=mp::material[imaterial].Ks;
		for(int nn=atoms::nearest_neighbour_list_si[atom];nn<atoms::nearest_neighbour_list_ei[atom];nn++){
			const double si_dot_eij=(Sx*atoms::eijx[nn]+Sy*atoms::eijy[nn]+Sz*atoms::eijz[nn]);
			energy+=Ks*si_dot_eij*si_dot_eij;
		}
	}
	
	return energy;
}

/// @brief Calculates the magnetostatic energy for a single spin.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2012. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    08/06/2012
///
/// @param[in] atom atom number 
/// @param[in] imaterial material of local atom
/// @param[in] Sx x-spin of local atom  
/// @param[in] Sy y-spin of local atom 
/// @param[in] Sz z-spin of local atom 
/// @return magnetostatic energy
///
/// @internal
///	Created:		08/06/2012
///	Revision:	  ---
///=====================================================================================
///
inline double spin_magnetostatic_energy(const int atom, const int imaterial, const double Sx, const double Sy, const double Sz){
	
	return -1.0*(atoms::x_dipolar_field_array[atom]*Sx+atoms::y_dipolar_field_array[atom]*Sy+atoms::z_dipolar_field_array[atom]*Sz);
}

/// @brief Calculates the total energy for a single spin.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, rfle500@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @param[in] atom atom number 
/// @return total spin energy
///
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///
double calculate_spin_energy(const int atom, const int AtomExchangeType){
	
	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "calculate_spin_energy has been called" << std::endl;
	
	// Local spin value
	const double Sx=atoms::x_spin_array[atom];
	const double Sy=atoms::y_spin_array[atom];
	const double Sz=atoms::z_spin_array[atom];
	
	// Determine neighbour material
	const int imaterial=atoms::type_array[atom];
	
	// Initialise energy to zero
	double energy=0.0;
	
	// Calculate total spin energy
	switch(AtomExchangeType){
		case 0: energy+=spin_exchange_energy_isotropic(atom, Sx, Sy, Sz); break;
		case 1: energy+=spin_exchange_energy_vector(atom, Sx, Sy, Sz); break;
		case 2: energy+=spin_exchange_energy_tensor(atom, Sx, Sy, Sz); break;
		default: zlog << zTs() << "Error. atoms::exchange_type has value " << AtomExchangeType << " which is outside of valid range 0-2. Exiting." << std::endl; err::vexit();  
	}
	energy+=spin_applied_field_energy(atom, imaterial, Sx, Sy, Sz);
	energy+=spin_uniaxial_energy(atom, imaterial, Sx, Sy, Sz);
	energy+=spin_surface_anisotropy_energy(atom, imaterial, Sx, Sy, Sz);
	energy+=spin_magnetostatic_energy(atom, imaterial, Sx, Sy, Sz);
	
	return energy; // Tesla
}

} // end of namespace sim

