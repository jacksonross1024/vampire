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
bool DMI = false;

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
std::vector < spin > all_m_atoms_offset;
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
   system_size_x = number_of_unit_cells_x*a0x;
   system_size_y = number_of_unit_cells_y*a1y;
   // number_of_unit_cells_x++;
   // number_of_unit_cells_y++;
   system_size_z = number_of_unit_cells_z*c0;

   std::cout << "Creating base lattice " << number_of_unit_cells_x << " by " << number_of_unit_cells_y << " rhombic unit cells" << std::endl;
   twist_angle = twist_angle*M_PI/180.0;
   resize_arrays(Einter_Cr1, 895,6);
   resize_arrays(Einter_Cr2, 895,6);
   resize_arrays(Einter_Cr3, 895,6);
   resize_arrays(Einter_Cr4, 895,6);

   Eintra_Cr1_1NN.resize(10);
   Eintra_Cr2_1NN.resize(10);
   Eintra_Cr3_1NN.resize(10);
   Eintra_Cr4_1NN.resize(10);

   Eintra_Cr1_2NN.resize(10);
   Eintra_Cr2_2NN.resize(10);
   Eintra_Cr3_2NN.resize(10);
   Eintra_Cr4_2NN.resize(10);

   Eintra_Cr1_3NN.resize(10);
   Eintra_Cr2_3NN.resize(10);
   Eintra_Cr3_3NN.resize(10);
   Eintra_Cr4_3NN.resize(10);

   for(int i = 0; i < Eintra_Cr1_1NN.size(); i++) {
      Eintra_Cr1_1NN[i].resize(11);
      Eintra_Cr2_1NN[i].resize(11);
      Eintra_Cr3_1NN[i].resize(11);
      Eintra_Cr4_1NN[i].resize(11);

      Eintra_Cr1_2NN[i].resize(11);
      Eintra_Cr2_2NN[i].resize(11);
      Eintra_Cr3_2NN[i].resize(11);
      Eintra_Cr4_2NN[i].resize(11);

      Eintra_Cr1_3NN[i].resize(11);
      Eintra_Cr2_3NN[i].resize(11);
      Eintra_Cr3_3NN[i].resize(11);
      Eintra_Cr4_3NN[i].resize(11);
      for(int j = 0; j < Eintra_Cr1_1NN[i].size(); j++) {
         Eintra_Cr1_1NN[i][j].resize(3);
         Eintra_Cr2_1NN[i][j].resize(3);//,std::vector<std::array<std::array<double,4>, 3> >(11));
         Eintra_Cr3_1NN[i][j].resize(3);//,std::vector<std::array<std::array<double,4>, 3> >(11));
         Eintra_Cr4_1NN[i][j].resize(3);

         Eintra_Cr1_2NN[i][j].resize(6);
         Eintra_Cr2_2NN[i][j].resize(6);//,std::vector<std::array<std::array<double,4>, 3> >(11));
         Eintra_Cr3_2NN[i][j].resize(6);//,std::vector<std::array<std::array<double,4>, 3> >(11));
         Eintra_Cr4_2NN[i][j].resize(6);

         Eintra_Cr1_3NN[i][j].resize(3);
         Eintra_Cr2_3NN[i][j].resize(3);//,std::vector<std::array<std::array<double,4>, 3> >(11));
         Eintra_Cr3_3NN[i][j].resize(3);//,std::vector<std::array<std::array<double,4>, 3> >(11));
         Eintra_Cr4_3NN[i][j].resize(3);
         for(int k = 0; k < Eintra_Cr1_1NN[i][j].size(); k++ ) {
            Eintra_Cr1_1NN[i][j][k].resize(4, -11.0);
            Eintra_Cr2_1NN[i][j][k].resize(4, -12.0);
            Eintra_Cr3_1NN[i][j][k].resize(4, -13.0);
            Eintra_Cr4_1NN[i][j][k].resize(4, -14.0);
            Eintra_Cr1_3NN[i][j][k].resize(4, -31.0);
            Eintra_Cr2_3NN[i][j][k].resize(4, -32.0);
            Eintra_Cr3_3NN[i][j][k].resize(4, -33.0);
            Eintra_Cr4_3NN[i][j][k].resize(4, -34.0);
         }
         for(int k = 0; k < Eintra_Cr1_2NN[i][j].size(); k++ ) {
            Eintra_Cr1_2NN[i][j][k].resize(4, -21.0);
            Eintra_Cr2_2NN[i][j][k].resize(4, -22.0);
            Eintra_Cr3_2NN[i][j][k].resize(4, -23.0);
            Eintra_Cr4_2NN[i][j][k].resize(4, -24.0);
         }
      }
   }
   
   // resize_arrays(Dx_inter, 201,201);
   // resize_arrays(Dy_inter, 201,201);
   // resize_arrays(Dz_inter, 201,201);
   // resize_arrays(Dx_intra, 201,201);
   // resize_arrays(Dy_intra, 201,201);
   // resize_arrays(Dz_intra, 201,201);
   // resize_arrays(Dx_intra2, 201,201);
   // resize_arrays(Dy_intra2, 201,201);
   // resize_arrays(Dz_intra2, 201,201);

   int estimated_system_spins = round(system_size_x*system_size_y*8.0/41.8);

   all_m_atoms.reserve(estimated_system_spins);
   row1.reserve(int(estimated_system_spins/2));
   // row2.reserve(int(estimated_system_spins/4));
   row3.reserve(int(estimated_system_spins/2));
   // row4.reserve(int(estimated_system_spins/4));

   
   unit_cell_shifts.resize(number_of_unit_cells_x);
   config_energy.resize(number_of_unit_cells_x);
   for(int i = 0; i < number_of_unit_cells_x; i++) {
      unit_cell_shifts[i].resize(number_of_unit_cells_y);
      config_energy[i].resize(number_of_unit_cells_y);
      for(int j = 0; j < number_of_unit_cells_y; j++) {
         unit_cell_shifts[i][j].resize(3,0);
         config_energy[i][j].resize(30,0.0);
      }
   }
   // std::vector<std::vector<std::vector<double> > > config_energy(number_of_unit_cells_x, std::vector<std::vector<double> >(number_of_unit_cells_y, std::vector<double>(30,0.0)));
}
