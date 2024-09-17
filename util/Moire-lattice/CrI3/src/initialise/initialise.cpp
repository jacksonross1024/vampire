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

   number_of_unit_cells_x = system_size_x/a0x;
   number_of_unit_cells_y = system_size_y/a1y;
   system_size_z = number_of_unit_cells_z*c0;

   twist_angle = twist_angle*3.14159265359/180.0;
   resize_arrays(Jinter, 201,201);
   resize_arrays(Jintra1, 201,201);
   resize_arrays(Jintra2, 201,201);
   
   resize_arrays(Dx_inter, 201,201);
   resize_arrays(Dy_inter, 201,201);
   resize_arrays(Dz_inter, 201,201);
   resize_arrays(Dx_intra, 201,201);
   resize_arrays(Dy_intra, 201,201);
   resize_arrays(Dz_intra, 201,201);

   int estimated_system_spins = round(system_size_x*system_size_y*8.0/41.8);

   all_m_atoms.reserve(estimated_system_spins);
   row1.reserve(int(estimated_system_spins/4));
   row2.reserve(int(estimated_system_spins/4));
   row3.reserve(int(estimated_system_spins/4));
   row4.reserve(int(estimated_system_spins/4));
}
