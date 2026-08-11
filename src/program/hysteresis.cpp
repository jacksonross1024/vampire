//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo and Richard Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Standard Libraries
#include <cstdlib>
#include <iostream>
#include <string>

// Vampire Header files
#include "vmath.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

namespace{

void hysteresis_warning(const std::string& msg){
	if(vmpi::my_rank==0){
		terminaltextcolor(RED);
		std::cout << "Warning: " << msg << std::endl;
		terminaltextcolor(WHITE);
		zlog << zTs() << "Warning: " << msg << std::endl;
	}
}

int64_t hysteresis_leg_increment(int64_t Hstart, int64_t Hend, int64_t base_inc){
	if(Hstart==Hend) return 0;

	int64_t inc = base_inc;
	const bool need_increase = (Hend > Hstart);
	const bool inc_increases = (inc > 0);
	if(need_increase != inc_increases){
		hysteresis_warning("applied-field-strength-increment sign reversed to avoid infinite loop on hysteresis leg.");
		inc = -inc;
	}
	return inc;
}

} // namespace

namespace program{

/// @brief Function to calculate the hysteresis loop
///
/// @callgraph
/// @callergraph
///
/// @details Consists of a sequence of sub-calculations of fixed temperature. The system is initialised
/// ordered. After initialisation a whole hysteresis loop of the system and coercivity are calculated.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Weijia Fan, wf507@york.ac.uk
/// @version 1.0
/// @date    27/01/2010
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		27/01/2010
///	Revision:	  ---
///=====================================================================================
///
int hysteresis(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::hysteresis has been called" << std::endl;}

	// Setup field values and increment (uT)
	const int64_t iHeq = vmath::iround64(double(sim::equilibrium_H_field)*1.0E6);
	const int64_t iHmax = vmath::iround64(double(sim::Hmax)*1.0E6);
	const int64_t iHmin = vmath::iround64(double(sim::Hmin)*1.0E6);

	// Positive increment: Heq -> Hmax -> Hmin -> Heq; negative increment reverses path.
	const bool reverse_path = (sim::Hinc < 0.0);
	int64_t iHinc = vmath::iround64(fabs(sim::Hinc)*1.0E6);
	if(iHinc==0){
		hysteresis_warning("applied-field-strength-increment is zero; using 1 uT for hysteresis loop.");
		iHinc = 1;
	}
	const int64_t base_inc = reverse_path ? -iHinc : iHinc;

	uint64_t start_time;

	// Equilibrate system at equilibrium field strength
	sim::actual_H_field = sim::equilibrium_H_field;
	sim::actual_H_vector[0] = sim::equilibrium_H_vector[0];
	sim::actual_H_vector[1] = sim::equilibrium_H_vector[1];
	sim::actual_H_vector[2] = sim::equilibrium_H_vector[2];

	// Initialise sim::integrate only if it not a checkpoint
	if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
	else sim::integrate(sim::equilibration_time);

	// Field loop legs always use named Hmax/Hmin regardless of min>max ordering.
	int64_t leg_start[3];
	int64_t leg_end[3];
	if(!reverse_path){
		leg_start[0] = iHeq; leg_end[0] = iHmax;
		leg_start[1] = iHmax; leg_end[1] = iHmin;
		leg_start[2] = iHmin; leg_end[2] = iHeq;
	}
	else{
		leg_start[0] = iHeq; leg_end[0] = iHmin;
		leg_start[1] = iHmin; leg_end[1] = iHmax;
		leg_start[2] = iHmax; leg_end[2] = iHeq;
	}

	bool resume_checkpoint = sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag;
	int64_t leg = resume_checkpoint ? sim::parity : 0;

	// Perform field loop over legs
	for(; leg < 3; leg++){

		const int64_t Hstart = leg_start[leg];
		const int64_t Hend = leg_end[leg];
		const int64_t inc = hysteresis_leg_increment(Hstart, Hend, base_inc);

		int64_t Hfield = Hstart;
		if(resume_checkpoint){
			Hfield = static_cast<int64_t>(sim::iH);
			resume_checkpoint = false;
		}
		if(Hstart == Hend) Hfield = Hend;

		sim::parity = leg;

		// Perform field loop for current leg
		while(true){

			// Set field strength (Tesla); direction set at equilibration
			sim::actual_H_field = double(Hfield)*1.0e-6;

			// Reset start time
			start_time=sim::time;

			// Reset mean magnetisation counters
			stats::reset();

			// Integrate system
			while(sim::time<sim::loop_time+start_time){

				// Integrate system
				sim::integrate(sim::partial_time);

				// Calculate mag_m, mag
				stats::update();

			} // End of integration loop

			// Save field for checkpointing (uT)
			sim::iH=Hfield;

			// Output to screen and file after each field
			vout::data();

			if(Hfield==Hend) break;

			// Advance toward leg end, clamping to Hend on overshoot or checkpoint resume past end
			if(Hstart > Hend){
				if(Hfield <= Hend) Hfield = Hend;
				else{
					Hfield += inc;
					if(Hfield < Hend) Hfield = Hend;
				}
			}
			else if(Hstart < Hend){
				if(Hfield >= Hend) Hfield = Hend;
				else{
					Hfield += inc;
					if(Hfield > Hend) Hfield = Hend;
				}
			}

		} // End of field loop

	} // End of leg loop

	return EXIT_SUCCESS;

} // End of hysteresis program

}//end of namespace program
