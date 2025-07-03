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

void build_fe5gete2(unitcell::unit_cell_t& unit_cell){

   // Set basic unit cell properties
   unit_cell.dimensions[0] = 1.0;
   unit_cell.dimensions[1] = 1.7321695761;
   unit_cell.dimensions[2] = 4.6608478803;

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
   unit_cell.interaction_range=1.0;
   unit_cell.atom.resize(16);
   unit_cell.surface_threshold=12;

      //16 atoms--monolayer
   //0   0.4010  0.6946  1.8690 90.0000 90.0000 90.0000
   //Te   0.000000(x1)  0.0(y1)  0.328718(z1) -0.079984  1.000000  0.100000  0.100000  0.100000
   //Fe   0.000000(x1)  0.0(y1)  0.462499(z3)  1.499891  4.000000  0.100000  0.100000  0.100000
   //Fe   0.000000(x1)  0.333(y3)  0.542088(z5)  2.080354  3.000000  0.100000  0.100000  0.100000
   //Fe   0.000000(x1)  0.333(y3)  0.405582(z2)  2.599483  5.000000  0.100000  0.100000  0.100000
   //Fe   0.000000(x1)  0.667(y4)  0.622811(z7) -0.029894  1.000000  0.100000  0.100000  0.100000
   //Ge   0.000000(x1)  0.667(y4)  0.491196(z4) -0.031399  1.000000  0.100000  0.100000  0.100000
   //Fe   0.000000(x1)  0.0(y1)  0.590562(z6)  2.336375  2.000000  0.100000  0.100000  0.100000
   //Te   0.5(x2)  0.83331(y5)  0.687683(z8)  0.010792  2.000000  0.100000  0.100000  0.100000
   //Fe   0.500000(x2)  0.1667(y2)  0.622811(z7) -0.029894  1.000000  0.100000  0.100000  0.100000
   //Ge   0.500000(x2)  0.1667(y2)  0.491196(z4) -0.031399  1.000000  0.100000  0.100000  0.100000
   // Fe   0.500000(x2)  0.5(y6)  0.590562(z6)  2.336375  2.000000  0.100000  0.100000  0.100000
   // Te   0.500000(x2)  0.5(y6)  0.328718(z1) -0.079984  1.000000  0.100000  0.100000  0.100000
   // Fe   0.500000(x2)  0.5(y6)  0.462499(z3)  1.499891  4.000000  0.100000  0.100000  0.100000
   // Fe   0.500000(x2)  0.83333(y5) 0.542088(z5)  2.080354  3.000000  0.100000  0.100000  0.100000
   // Fe   0.500000(x2)  0.8333(y5)  0.405582(z2)  2.599483  5.000000  0.100000  0.100000  0.100000
   // Te   0.0(x1)  0.333(y3) 0.687683(z8)  0.010792  2.000000  0.100000  0.100000  0.100000

   //-----------------------------
   double x1 = 0.0;
   double x2 = 0.5;

   double y1 = 0.0;
   double y2 = 0.1667;
   double y3 = 0.333;
   double y4 = 0.667;
   double y5 = 0.8333;
   double y6 = 0.5;

   double z1 = 0.328718;
   double z2 = 0.405582;
   double z3 = 0.462499;
   double z4 = 0.491196;
   double z5 = 0.542088;
   double z6 = 0.590562;
   double z7 = 0.622811;
   double z8 = 0.687683;

      //Te1-nm
  unit_cell.atom.at(0).x=x1;
  unit_cell.atom.at(0).y=y1;
  unit_cell.atom.at(0).z=z1;
  unit_cell.atom.at(0).mat=5;
  unit_cell.atom.at(0).hc=4;
  unit_cell.atom.at(0).lc=6;
  unit_cell.atom[0].nm    = true;

  //Fe-4
  unit_cell.atom.at(1).x=x1;
  unit_cell.atom.at(1).y=y1;
  unit_cell.atom.at(1).z=z3;
  unit_cell.atom.at(1).mat=3;
  unit_cell.atom.at(1).hc=1;
  unit_cell.atom.at(1).lc=1;
  
  //Fe-3
  unit_cell.atom.at(2).x=x1;
  unit_cell.atom.at(2).y=y3;
  unit_cell.atom.at(2).z=z5;
  unit_cell.atom.at(2).mat=2;
  unit_cell.atom.at(2).hc=1;
  unit_cell.atom.at(2).lc=3;
  
  //Fe-5
  unit_cell.atom.at(3).x=x1;
  unit_cell.atom.at(3).y=y3;
  unit_cell.atom.at(3).z=z2;
  unit_cell.atom.at(3).mat=4;
  unit_cell.atom.at(3).hc=1;
  unit_cell.atom.at(3).lc=2;
  
  //Fe-1
  unit_cell.atom.at(4).x=x1;
  unit_cell.atom.at(4).y=y4;
  unit_cell.atom.at(4).z=z7;
  unit_cell.atom.at(4).mat=0;
  unit_cell.atom.at(4).hc=1;
  unit_cell.atom.at(4).lc=4;

   //Ge-1
  unit_cell.atom.at(5).x=x1;
  unit_cell.atom.at(5).y=y4;
  unit_cell.atom.at(5).z=z4;
  unit_cell.atom.at(5).mat=6;
  unit_cell.atom.at(5).hc=4;
  unit_cell.atom.at(5).lc=6;
  unit_cell.atom[5].nm    = true;

   //Fe-2
  unit_cell.atom.at(6).x=x1;
  unit_cell.atom.at(6).y=y1;
  unit_cell.atom.at(6).z=z6;
  unit_cell.atom.at(6).mat=1;
  unit_cell.atom.at(6).hc=1;
  unit_cell.atom.at(6).lc=5;

  //Te-2
  unit_cell.atom.at(7).x=x2;
  unit_cell.atom.at(7).y=y5;
  unit_cell.atom.at(7).z=z8;
  unit_cell.atom.at(7).mat=5;
  unit_cell.atom.at(7).hc=1;
  unit_cell.atom.at(7).lc=1;
  unit_cell.atom[7].nm    = true;

   //Fe-1
  unit_cell.atom.at(8).x=x2;
  unit_cell.atom.at(8).y=y2;
  unit_cell.atom.at(8).z=z7;
  unit_cell.atom.at(8).mat=0;
  unit_cell.atom.at(8).hc=1;
  unit_cell.atom.at(8).lc=3;

   //Ge-1
  unit_cell.atom.at(9).x=x2;
  unit_cell.atom.at(9).y=y2;
  unit_cell.atom.at(9).z=z4;
  unit_cell.atom.at(9).mat=6;
  unit_cell.atom.at(9).hc=1;
  unit_cell.atom.at(9).lc=2;
  unit_cell.atom[9].nm    = true;

   //Fe-2
  unit_cell.atom.at(10).x=x2;
  unit_cell.atom.at(10).y=y6;
  unit_cell.atom.at(10).z=z6;
  unit_cell.atom.at(10).mat=1;
  unit_cell.atom.at(10).hc=1;
  unit_cell.atom.at(10).lc=4;

   //Te-1
  unit_cell.atom.at(11).x=x2;
  unit_cell.atom.at(11).y=y6;
  unit_cell.atom.at(11).z=z1;
  unit_cell.atom.at(11).mat=5;
  unit_cell.atom.at(11).hc=1;
  unit_cell.atom.at(11).lc=6;
  unit_cell.atom[11].nm    = true;

   //Fe-4
  unit_cell.atom.at(12).x=x2;
  unit_cell.atom.at(12).y=y6;
  unit_cell.atom.at(12).z=z3;
  unit_cell.atom.at(12).mat=3;
  unit_cell.atom.at(12).hc=1;
  unit_cell.atom.at(12).lc=5;

   //Fe-3
  unit_cell.atom.at(13).x=x2;
  unit_cell.atom.at(13).y=y5;
  unit_cell.atom.at(13).z=z5;
  unit_cell.atom.at(13).mat=2;
  unit_cell.atom.at(13).hc=1;
  unit_cell.atom.at(13).lc=1;

   //Fe-5
  unit_cell.atom.at(14).x=x2;
  unit_cell.atom.at(14).y=y5;
  unit_cell.atom.at(14).z=z2;
  unit_cell.atom.at(14).mat=4;
  unit_cell.atom.at(14).hc=1;
  unit_cell.atom.at(14).lc=3;

 //Te-2
  unit_cell.atom.at(15).x=x1;
  unit_cell.atom.at(15).y=y3;
  unit_cell.atom.at(15).z=z8;
  unit_cell.atom.at(15).mat=5;
  unit_cell.atom.at(15).hc=1;
  unit_cell.atom.at(15).lc=2;

   unit_cell.cutoff_radius = 0.5; // normalised to unit cell size

   uc::internal::calculate_interactions(unit_cell);

   // Set actual unit cell size after calculating interactions
   unit_cell.dimensions[0] *= unitcell::internal::unit_cell_size_x;
   unit_cell.dimensions[1] *= unitcell::internal::unit_cell_size_y;
   unit_cell.dimensions[2] *= unitcell::internal::unit_cell_size_z;

   return;

}

} // end of internal namespace
} // end of unitcell namespace
