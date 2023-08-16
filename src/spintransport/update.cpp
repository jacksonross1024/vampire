//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <fstream>

// Vampire headers
#include "spintransport.hpp"
#include "sim.hpp"
// spintransport module headers
#include "internal.hpp"
#include "vmpi.hpp"
namespace spin_transport{

//---------------------------------------------------------------------------------------------------------
// Function to update spin torque field
//---------------------------------------------------------------------------------------------------------
void update(const unsigned int num_local_atoms,            // number of local atoms
            const std::vector<double>& atoms_x_spin_array, // x-spin vector of atoms
            const std::vector<double>& atoms_y_spin_array, // y-spin vector of atoms
            const std::vector<double>& atoms_z_spin_array, // z-spin-vector of atoms
            const std::vector<double>& atoms_m_spin_array  // moment of atoms
   ){

   //-------------------------------------------------------------------------
   // check that module is needed - if not do nothing
   //-------------------------------------------------------------------------
   if( st::internal::enabled == false ) return;

   //-------------------------------------------------------------------------
   // check that it is time to update
   //-------------------------------------------------------------------------
   if(st::internal::time_counter != st::internal::update_rate){
      // if not increment counter and do nothing
      st::internal::time_counter++;
      return;
   }
   // otherwise reset counter and continue
   else{ st::internal::time_counter = 1; }

   //---------------------------------------------------------------------------------------------------------
   // update cell magnetizations
   //---------------------------------------------------------------------------------------------------------
   st::internal::calculate_cell_magnetization(num_local_atoms, atoms_x_spin_array, atoms_y_spin_array,
                                                               atoms_z_spin_array, atoms_m_spin_array);

   //---------------------------------------------------------------------------------------------------------
   // calculate magnetoresistance
   //---------------------------------------------------------------------------------------------------------
   st::internal::calculate_magnetoresistance();

   //---------------------------------------------------------------------------------------------------------
   // test output of cell-level spin transport data
   //---------------------------------------------------------------------------------------------------------
 if(vmpi::my_rank == 0) {
    if(sim::time % 500 == 0) {
      std::ofstream ofile("spin-transfer/stdata-" + std::to_string(sim::time) + ".txt");
   ofile << "x_i  y_i   z_i   mx_i  my_i  mz_i  Tx_i  Ty_i  Tz_i  A  R_i   R_s_i" << std::endl;
   for(int i =0; i< st::internal::total_num_cells; i++){
      const double isat = st::internal::cell_isaturation[i];
      ofile << 
               st::internal::cell_position[3*i+0] << "\t" <<
               st::internal::cell_position[3*i+1] << "\t" <<
               st::internal::cell_position[3*i+2] << "\t" <<
               st::internal::cell_magnetization[3*i+0] * isat << "\t" <<
               st::internal::cell_magnetization[3*i+1] * isat << "\t" <<
               st::internal::cell_magnetization[3*i+2] * isat << "\t" <<
               st::internal::cell_spin_torque_fields[3*i+0] << "\t" <<
               st::internal::cell_spin_torque_fields[3*i+1] << "\t" <<
               st::internal::cell_spin_torque_fields[3*i+2] << "\t" <<
               st::total_current*1e20/(cs::system_dimensions[0]*cs::system_dimensions[1]) << "\t" << //A/A^2 -> 1e10 A/m
               st::internal::cell_resistance[i] << "\t" <<
               st::internal::cell_spin_resistance[i] << "\t" <<
               st::total_resistance << 
                std::endl;
   }
   ofile.close();
   }
 }
   return;

}

} // end of spin_transport namespace
