#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "exchange.hpp"
#include "positions.hpp"
#include "initialise.hpp"

double system_size_x;
double system_size_y;
double system_size_z;
int number_of_unit_cells_z;

double twist_angle;
double twist_loction;

double dmi12 = 0.0; // DMI constant between layers 1-2
double dmi23 = 0.0; // DMI constant between layers 2-3
double dmi34 = 0.0; // DMI constant between layers 3-4

double dmi_decay = 1.0; // distance-dependent DMI

double exchange12 = 0.0; // exchange constant between layers 1-2
double exchange23 = 0.0; // exchange constant between layers 2-3
double exchange34 = 0.0; // exchange constant between layers 3-4

double separation = 0.0; // distance between layers 2-3

std::vector < std::vector < std::vector <int> > > unit_cell_shifts;
std::vector < spin > atom(num_atoms);
std::vector < spin > nm_atom(num_nm_atoms);
std::vector < spin > row1;
std::vector < spin > row2;
std::vector < spin > row3;
std::vector < spin > row4;
std::vector < spin > all_m_atoms;
std::vector < spin > all_nm_atoms;


std::ofstream outfile4 ("interactions.ucf");

void resize_arrays(std::vector < std::vector < double > > &A, int sizex, int sizey){
   A.resize(sizex);
      for (int i = 0; i < sizex; ++i)
          A[i].resize(sizey, 0.0);
}

void initialise_variables(){

   number_of_unit_cells_x = ceil(system_size_x/a0x);
   number_of_unit_cells_y = ceil(system_size_y/a1y);
   system_size_z = number_of_unit_cells_z*c0;

   std::cout << "Creating base lattice " << number_of_unit_cells_x << " by " << number_of_unit_cells_y << " rhombic unit cells" << std::endl;
   twist_angle = twist_angle*3.14159265359/180.0;
   resize_arrays(Jinter, 201,201);
   resize_arrays(Jintra1, 201,201);
   resize_arrays(Jintra2, 201,201);
   
   D_inter.resize(201);
   D_intra.resize(201);
   for(int i = 0; i < D_inter.size(); i++) {
      D_inter[i].resize(201);
      D_intra[i].resize(201);
      for(int j = 0; j < D_inter[i].size(); j++) {
         D_inter[i][j].resize(3,0);
         D_intra[i][j].resize(9,0);
      }

   }
   // resize_arrays(Dx_inter, 201,201);
   // resize_arrays(Dy_inter, 201,201);
   // resize_arrays(Dz_inter, 201,201);
   // resize_arrays(Dx_intra, 201,201);
   // resize_arrays(Dy_intra, 201,201);
   // resize_arrays(Dz_intra, 201,201);

   int estimated_system_spins = round(system_size_x*system_size_y*8.0/41.8);

   all_m_atoms.reserve(estimated_system_spins);
   row1.reserve(int(estimated_system_spins/4));
   row2.reserve(int(estimated_system_spins/4));
   row3.reserve(int(estimated_system_spins/4));
   row4.reserve(int(estimated_system_spins/4));

      

   unit_cell_shifts.resize(number_of_unit_cells_x);
   for(int i = 0; i < number_of_unit_cells_x; i++) {
      unit_cell_shifts[i].resize(number_of_unit_cells_y);
      for(int j = 0; j < number_of_unit_cells_y; j++) {
         unit_cell_shifts[i][j].resize(3,0);
      }
   }
}
