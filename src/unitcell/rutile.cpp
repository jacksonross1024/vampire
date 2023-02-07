//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) J L Ross Richard F L Evans 2023. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>

// Vampire headers
#include "unitcell.hpp"

// unitcell module headers
#include "internal.hpp"

namespace unitcell{
namespace internal{

void build_rutile(unitcell::unit_cell_t& unit_cell){

   // Set basic unit cell properties
   unit_cell.dimensions[0] = 1.0;
   unit_cell.dimensions[1] = 1.0;
   unit_cell.dimensions[2] = 1.0;

   unit_cell.shape[0][0]=1.0;
   unit_cell.shape[0][1]=0.0;
   unit_cell.shape[0][2]=0.0;

   unit_cell.shape[1][0]=0.0;
   unit_cell.shape[1][1]=1.0;
   unit_cell.shape[1][2]=0.0;

   unit_cell.shape[2][0]=0.0;
   unit_cell.shape[2][1]=0.0;
   unit_cell.shape[2][2]=1.0;

   unit_cell.lcsize=2;
   unit_cell.hcsize=2;
   unit_cell.interaction_range=2; 
   unit_cell.atom.resize(2);
   unit_cell.surface_threshold=14;
   //-----------------------------
   // nmagnetic atoms only for 0,0,0 and 0.5,0.5,0.5 for now
   unit_cell.atom[0].x=0.0;
   unit_cell.atom[0].y=0.0;
   unit_cell.atom[0].z=0.0;
   unit_cell.atom[0].mat = uc::internal::sublattice_materials ? 0 : 0; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[0].lc=0;
   unit_cell.atom[0].hc=0;
   unit_cell.atom[0].ni=14;
   //-----------------------------
   unit_cell.atom[1].x=0.5;
   unit_cell.atom[1].y=0.5;
   unit_cell.atom[1].z=0.5;
   unit_cell.atom[1].mat = uc::internal::sublattice_materials ? 1 : 0; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[1].lc=1;
   unit_cell.atom[1].hc=1;
   unit_cell.atom[1].ni=14;
   //-----------------------------

   unit_cell.cutoff_radius = 1.0; // extended interaction // normalised to unit cell size

   uc::internal::calculate_interactions(unit_cell);

   // Set actual unit cell size after calculating interactions
   unit_cell.dimensions[0] *= unitcell::internal::unit_cell_size_x;
   unit_cell.dimensions[1] *= unitcell::internal::unit_cell_size_y;
   unit_cell.dimensions[2] *= unitcell::internal::unit_cell_size_z;

   return;

}

} // end of internal namespace
} // end of unitcell namespace
