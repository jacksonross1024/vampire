//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2016. All rights reserved.
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

void build_mn3sn(unitcell::unit_cell_t& unit_cell){

   // rescale y and z unit cell size to preserve unit cell dimensions
   unit_cell.dimensions[0]=0.5746444617738683;
   unit_cell.dimensions[1]=1.0;
   unit_cell.dimensions[2]=0.46002015128929274;

   unit_cell.shape[0][0]=1.0;
   unit_cell.shape[0][1]=0.0;
   unit_cell.shape[0][2]=0.0;

   unit_cell.shape[1][0]=0.0;
   unit_cell.shape[1][1]=1.0;
   unit_cell.shape[1][2]=0.0;

   unit_cell.shape[2][0]=0.0;
   unit_cell.shape[2][1]=0.0;
   unit_cell.shape[2][2]=1.0;

   unit_cell.lcsize=16;
   unit_cell.hcsize=2;
   unit_cell.interaction_range=1;
   unit_cell.atom.resize(16);
   unit_cell.surface_threshold=3;
   //-----------------------------
   unit_cell.atom[0].x = 0;  
   unit_cell.atom[0].y = 0.08294278;
   unit_cell.atom[0].z = 0.75;
   unit_cell.atom[0].mat = uc::internal::sublattice_materials ? 0 : 0; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[0].lc = 0;   // running index 
   unit_cell.atom[0].hc = 1;   // height category
   unit_cell.atom[0].ni = 20;  // (est.) of no. of interactions
   //-----------------------------
   unit_cell.atom[1].x = 0.25;
   unit_cell.atom[1].y = 0;
   unit_cell.atom[1].z = 0.25;
   unit_cell.atom[1].mat = uc::internal::sublattice_materials ? 1 : 1; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[1].lc = 1;
   unit_cell.atom[1].hc = 0;
   unit_cell.atom[1].ni = 20;
   //-----------------------------
   unit_cell.atom[2].x = 0.75;
   unit_cell.atom[2].y = 0;
   unit_cell.atom[2].z = 0.25;
   unit_cell.atom[2].mat = uc::internal::sublattice_materials ? 2 : 2; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[2].lc = 2;
   unit_cell.atom[2].hc = 0;
   unit_cell.atom[2].ni = 20;
   //-----------------------------
   unit_cell.atom[3].x = 0.5;
   unit_cell.atom[3].y = 0.08294278;
   unit_cell.atom[3].z = 0.75;
   unit_cell.atom[3].mat = uc::internal::sublattice_materials ? 3 : 3; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[3].lc = 3;
   unit_cell.atom[3].hc = 1;
   unit_cell.atom[3].ni = 0;
   unit_cell.atom[3].nm = true;
   //-----------------------------
   unit_cell.atom[4].x = 0.5;  // replace by relative coordinates
   unit_cell.atom[4].y = 0.25;
   unit_cell.atom[4].z = 0.25;
   unit_cell.atom[4].mat = uc::internal::sublattice_materials ? 0 : 0; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[4].lc = 4;   // running index 
   unit_cell.atom[4].hc = 0;   // height category
   unit_cell.atom[4].ni = 20;  // (est.) of no. of interactions
   //-----------------------------
   unit_cell.atom[5].x = 0.25;
   unit_cell.atom[5].y = 0.33294278;
   unit_cell.atom[5].z = 0.75;
   unit_cell.atom[5].mat = uc::internal::sublattice_materials ? 2 : 1; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[5].lc = 5;
   unit_cell.atom[5].hc = 1;
   unit_cell.atom[5].ni = 20;
   //-----------------------------
   unit_cell.atom[6].x = 0.75;
   unit_cell.atom[6].y = 0.33294278;
   unit_cell.atom[6].z = 0.75;
   unit_cell.atom[6].mat = uc::internal::sublattice_materials ? 1 : 2; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[6].lc = 6;
   unit_cell.atom[6].hc = 1;
   unit_cell.atom[6].ni = 20;
   //-----------------------------
   unit_cell.atom[7].x = 0.0;
   unit_cell.atom[7].y = 0.25;
   unit_cell.atom[7].z = 0.25;
   unit_cell.atom[7].mat = uc::internal::sublattice_materials ? 3 : 3; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[7].lc = 7;
   unit_cell.atom[7].hc = 0;
   unit_cell.atom[7].ni = 0;
   unit_cell.atom[7].nm = true;
      //-----------------------------
   unit_cell.atom[8].x = 0.5;
   unit_cell.atom[8].y = 0.58294278;
   unit_cell.atom[8].z = 0.75;
   unit_cell.atom[8].mat = uc::internal::sublattice_materials ? 0 : 0; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[8].lc = 8;
   unit_cell.atom[8].hc = 1;
   unit_cell.atom[8].ni = 20;
         //-----------------------------
   unit_cell.atom[9].x = 0.75;
   unit_cell.atom[9].y = 0.5;
   unit_cell.atom[9].z = 0.25;
   unit_cell.atom[9].mat = uc::internal::sublattice_materials ? 1 : 1; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[9].lc = 9;
   unit_cell.atom[9].hc = 0;
   unit_cell.atom[9].ni = 20;
            //-----------------------------
   unit_cell.atom[10].x = 0.25;
   unit_cell.atom[10].y = 0.5;
   unit_cell.atom[10].z = 0.25;
   unit_cell.atom[10].mat = uc::internal::sublattice_materials ? 2 : 2; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[10].lc = 10;
   unit_cell.atom[10].hc = 0;
   unit_cell.atom[10].ni = 20;
           //-----------------------------
   unit_cell.atom[11].x = 0.0;
   unit_cell.atom[11].y = 0.58294278;
   unit_cell.atom[11].z = 0.75;
   unit_cell.atom[11].mat = uc::internal::sublattice_materials ? 3 : 3; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[11].lc = 11;
   unit_cell.atom[11].hc = 1;
   unit_cell.atom[11].ni = 0;
   unit_cell.atom[11].nm = true;
              //-----------------------------
   unit_cell.atom[12].x = 0.0;
   unit_cell.atom[12].y = 0.75;
   unit_cell.atom[12].z = 0.25;
   unit_cell.atom[12].mat = uc::internal::sublattice_materials ? 0 : 0; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[12].lc = 12;
   unit_cell.atom[12].hc = 0;
   unit_cell.atom[12].ni = 20;
                 //-----------------------------
   unit_cell.atom[13].x = 0.75;
   unit_cell.atom[13].y = 0.83294278;
   unit_cell.atom[13].z = 0.75;
   unit_cell.atom[13].mat = uc::internal::sublattice_materials ? 2 : 1; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[13].lc = 13;
   unit_cell.atom[13].hc = 1;
   unit_cell.atom[13].ni = 20;
                    //-----------------------------
   unit_cell.atom[14].x = 0.25;
   unit_cell.atom[14].y = 0.83294278;
   unit_cell.atom[14].z = 0.75;
   unit_cell.atom[14].mat = uc::internal::sublattice_materials ? 1 : 2; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[14].lc = 14;
   unit_cell.atom[14].hc = 1;
   unit_cell.atom[14].ni = 20;
                       //-----------------------------
   unit_cell.atom[15].x = 0.5;
   unit_cell.atom[15].y = 0.75;
   unit_cell.atom[15].z = 0.25;
   unit_cell.atom[15].mat = uc::internal::sublattice_materials ? 3 : 3; // if sublattice material is defined, then identify at as same as ID
   unit_cell.atom[15].lc = 15;
   unit_cell.atom[15].hc = 0;
   unit_cell.atom[15].ni = 0;
   unit_cell.atom[15].nm = true;

   unit_cell.cutoff_radius = 1; // normalised to unit cell size

   // Set actual unit cell size after calculating interactions
   unit_cell.dimensions[0] *= unitcell::internal::unit_cell_size_x;
   unit_cell.dimensions[1] *= unitcell::internal::unit_cell_size_y;
   unit_cell.dimensions[2] *= unitcell::internal::unit_cell_size_z;

   uc::internal::calculate_interactions(unit_cell);

   return;

}

} // end of internal namespace
} // end of unitcell namespace
