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

void build_fe3gete2(unitcell::unit_cell_t& unit_cell){

   // Set basic unit cell properties
   unit_cell.dimensions[0] = 0.49121942;
   unit_cell.dimensions[1] = 0.4254113;
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

   unit_cell.lcsize=4;
   unit_cell.hcsize=2;
   unit_cell.interaction_range=0.25;
   unit_cell.atom.resize(24);
   unit_cell.surface_threshold=12;
   //-----------------------------
   double x1 = 0.125;
   double x2 = 0.375;
   double x3 = 0.625;
   double x4 = 0.875;

   double y1 = 0.0833574739281576;
   double y2 = 0.250072421784473;
   double y3 = 0.416642526071842;
   double y4 = 0.583357473928158;
   double y5 = 0.750072421784473;
   double y6 = 0.916642526071842;

   double z1 = 0.176967157557459;
   double z2 = 0.249984595477232;
   double z3 = 0.323063651488077;
   double z4 = 0.676936348511923;
   double z5 = 0.750015404522768;
   double z6 = 0.823032842442541;

  unit_cell.atom.at(0).x=x1;
  unit_cell.atom.at(0).y=y2;
  unit_cell.atom.at(0).z=z5;
  unit_cell.atom.at(0).mat=1;
  unit_cell.atom.at(0).hc=4;
  unit_cell.atom.at(0).lc=6;
  
  unit_cell.atom.at(1).x=x2;
  unit_cell.atom.at(1).y=y1;
  unit_cell.atom.at(1).z=z1;
  unit_cell.atom.at(1).mat=0;
  unit_cell.atom.at(1).hc=1;
  unit_cell.atom.at(1).lc=1;
  
  unit_cell.atom.at(2).x=x2;
  unit_cell.atom.at(2).y=y3;
  unit_cell.atom.at(2).z=y2;
  unit_cell.atom.at(2).mat=1;
  unit_cell.atom.at(2).hc=1;
  unit_cell.atom.at(2).lc=3;
  
  unit_cell.atom.at(3).x=x2;
  unit_cell.atom.at(3).y=y1;
  unit_cell.atom.at(3).z=z3;
  unit_cell.atom.at(3).mat=0;
  unit_cell.atom.at(3).hc=1;
  unit_cell.atom.at(3).lc=2;
  
  unit_cell.atom.at(4).x=x2;
  unit_cell.atom.at(4).y=y1;
  unit_cell.atom.at(4).z=z4;
  unit_cell.atom.at(4).mat=0;
  unit_cell.atom.at(4).hc=1;
  unit_cell.atom.at(4).lc=4;

  unit_cell.atom.at(5).x=x3;
  unit_cell.atom.at(5).y=y2;
  unit_cell.atom.at(5).z=z5;
  unit_cell.atom.at(5).mat=1;
  unit_cell.atom.at(5).hc=4;
  unit_cell.atom.at(5).lc=6;

  unit_cell.atom.at(6).x=x2;
  unit_cell.atom.at(6).y=y1;
  unit_cell.atom.at(6).z=z6;
  unit_cell.atom.at(6).mat=0;
  unit_cell.atom.at(6).hc=1;
  unit_cell.atom.at(6).lc=5;

  unit_cell.atom.at(7).x=x1;
  unit_cell.atom.at(7).y=y4;
  unit_cell.atom.at(7).z=z1;
  unit_cell.atom.at(7).mat=0;
  unit_cell.atom.at(7).hc=1;
  unit_cell.atom.at(7).lc=1;

  unit_cell.atom.at(8).x=x1;
  unit_cell.atom.at(8).y=y6;
  unit_cell.atom.at(8).z=y2;
  unit_cell.atom.at(8).mat=1;
  unit_cell.atom.at(8).hc=1;
  unit_cell.atom.at(8).lc=3;

  unit_cell.atom.at(9).x=x1;
  unit_cell.atom.at(9).y=y4;
  unit_cell.atom.at(9).z=z3;
  unit_cell.atom.at(9).mat=0;
  unit_cell.atom.at(9).hc=1;
  unit_cell.atom.at(9).lc=2;

  unit_cell.atom.at(10).x=x1;
  unit_cell.atom.at(10).y=y4;
  unit_cell.atom.at(10).z=z4;
  unit_cell.atom.at(10).mat=0;
  unit_cell.atom.at(10).hc=1;
  unit_cell.atom.at(10).lc=4;

  unit_cell.atom.at(11).x=x2;
  unit_cell.atom.at(11).y=z5;
  unit_cell.atom.at(11).z=z5;
  unit_cell.atom.at(11).mat=1;
  unit_cell.atom.at(11).hc=1;
  unit_cell.atom.at(11).lc=6;

  unit_cell.atom.at(12).x=x1;
  unit_cell.atom.at(12).y=y4;
  unit_cell.atom.at(12).z=z6;
  unit_cell.atom.at(12).mat=0;
  unit_cell.atom.at(12).hc=1;
  unit_cell.atom.at(12).lc=5;

  unit_cell.atom.at(13).x=x4;
  unit_cell.atom.at(13).y=y1;
  unit_cell.atom.at(13).z=z1;
  unit_cell.atom.at(13).mat=0;
  unit_cell.atom.at(13).hc=1;
  unit_cell.atom.at(13).lc=1;

  unit_cell.atom.at(14).x=x4;
  unit_cell.atom.at(14).y=y3;
  unit_cell.atom.at(14).z=y2;
  unit_cell.atom.at(14).mat=1;
  unit_cell.atom.at(14).hc=1;
  unit_cell.atom.at(14).lc=3;

  unit_cell.atom.at(15).x=x4;
  unit_cell.atom.at(15).y=y1;
  unit_cell.atom.at(15).z=z3;
  unit_cell.atom.at(15).mat=0;
  unit_cell.atom.at(15).hc=1;
  unit_cell.atom.at(15).lc=2;

  unit_cell.atom.at(16).x=x4;
  unit_cell.atom.at(16).y=y1;
  unit_cell.atom.at(16).z=z4;
  unit_cell.atom.at(16).mat=0;
  unit_cell.atom.at(16).hc=1;
  unit_cell.atom.at(16).lc=4;

  unit_cell.atom.at(17).x=x4;
  unit_cell.atom.at(17).y=y1;
  unit_cell.atom.at(17).z=z6;
  unit_cell.atom.at(17).mat=0;
  unit_cell.atom.at(17).hc=1;
  unit_cell.atom.at(17).lc=5;

  unit_cell.atom.at(18).x=x3;
  unit_cell.atom.at(18).y=y4;
  unit_cell.atom.at(18).z=z1;
  unit_cell.atom.at(18).mat=0;
  unit_cell.atom.at(18).hc=1;
  unit_cell.atom.at(18).lc=1;

  unit_cell.atom.at(19).x=x3;
  unit_cell.atom.at(19).y=y6;
  unit_cell.atom.at(19).z=y2;
  unit_cell.atom.at(19).mat=1;
  unit_cell.atom.at(19).hc=1;
  unit_cell.atom.at(19).lc=3;

  unit_cell.atom.at(20).x=x3;
  unit_cell.atom.at(20).y=y4;
  unit_cell.atom.at(20).z=z3;
  unit_cell.atom.at(20).mat=0;
  unit_cell.atom.at(20).hc=1;
  unit_cell.atom.at(20).lc=2;

  unit_cell.atom.at(21).x=x3;
  unit_cell.atom.at(21).y=y4;
  unit_cell.atom.at(21).z=z4;
  unit_cell.atom.at(21).mat=0;
  unit_cell.atom.at(21).hc=1;
  unit_cell.atom.at(21).lc=4;
  
  unit_cell.atom.at(22).x=x4;
  unit_cell.atom.at(22).y=z5;
  unit_cell.atom.at(22).z=z5;
  unit_cell.atom.at(22).mat=1;
  unit_cell.atom.at(22).hc=1;
  unit_cell.atom.at(22).lc=6;

  unit_cell.atom.at(23).x=x3;
  unit_cell.atom.at(23).y=y4;
  unit_cell.atom.at(23).z=z6;
  unit_cell.atom.at(23).mat=0;
  unit_cell.atom.at(23).hc=1;
  unit_cell.atom.at(23).lc=5;

   unit_cell.cutoff_radius = 0.15; // normalised to unit cell size

   uc::internal::calculate_interactions(unit_cell);

   // Set actual unit cell size after calculating interactions
   unit_cell.dimensions[0] *= unitcell::internal::unit_cell_size_x;
   unit_cell.dimensions[1] *= unitcell::internal::unit_cell_size_y;
   unit_cell.dimensions[2] *= unitcell::internal::unit_cell_size_z;

   return;

}

} // end of internal namespace
} // end of unitcell namespace
