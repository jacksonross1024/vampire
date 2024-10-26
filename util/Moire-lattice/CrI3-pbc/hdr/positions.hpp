#ifndef POSITION_HPP
#define POSITION_HPP

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
// #include <array>

   extern double twist_angle;

   extern int number_of_unit_cells_x;
   extern int number_of_unit_cells_y;

   extern double a0x;
   extern double a0y;
   extern double a1x;
   extern double a1y;
   extern double c0;
   extern double a0z;

   extern int num_atoms;
   extern int num_nm_atoms;
   extern double twist_loction;
   extern int total_atoms;
   extern int total_nm_atoms;
   extern int num_above_atoms;
   extern int num_below_atoms;

   class spin {
      public:
         double x;
         double y;
         double z;
         int Gx = 0;
         int Gy = 0;
         int Gz = 0;
         double S;
         int id = 0;
         int l_id;
         int h_id;
         int unit_x;
         int unit_y;
         int dx = 0;
         int dy = 0;
         int inter1 = 0;
         int inter2 = 0;
         int inter3 = 0;
         int intra1 = 0;
         int intra2 = 0;
         int intra3 = 0;
   };

   class interaction {
      public:
         int id_i = -1;
         int id_j = -1;
         double J= 0.0;
         double Dx= 0.0;
         double Dy= 0.0;
         double Dz= 0.0;
         int pbc_x = 0;
         int pbc_y = 0;
         int pbc_z = 0;
   };

   extern std::vector < spin > atom;
   extern std::vector < spin > nm_atom;
   extern std::vector < spin > row1;
   extern std::vector < spin > row2;
   extern std::vector < spin > row3;
   extern std::vector < spin > row4;
   extern std::vector < spin > all_m_atoms;
   extern std::vector < spin > all_nm_atoms;
   extern std::vector < spin > all_m_atoms_offset;
   extern std::vector < spin > new_moire_lattice;
   extern std::vector < std::vector < std::vector <int> > > unit_cell_shifts;

   bool inside_system(double sx, double sy, double x, double y, double offset);
   bool inside_system(double x, double y, double offset);
   void read_in_atoms(std::string filename, int n_atoms, std::vector <spin > &atom2);
   void read_in_exchange(std::string filename, std::vector<std::vector<double> > &Jij);
   void read_in_inter_exchanges(std::string filename, std::vector<std::vector<double> > &Eij);
   void read_in_intra_exchanges(std::string filename, std::vector<std::vector<std::vector<std::vector<double> > > > &Eij_1NN, \
                                                   std::vector<std::vector< std::vector< std::vector<double> > > > &Eij_2NN, \
                                                   std::vector<std::vector< std::vector< std::vector<double> > > > &Eij_3NN );
   void read_in_dmi(std::string filename, std::vector < std::vector < double > > &Dx, std::vector < std::vector < double > > &Dy, std::vector < std::vector < double > > &Dz, bool reflect);
   void read_in_ucf(std::ifstream &ucf_file);
   void print_header();
   void create_magnetic_atom_list(std::string filename);
   void create_magnetic_atom_list_moire_unit(std::string filename, \
                  double Moire_a0x, double Moire_a0y, double Moire_a1x, double Moire_a1y, \
                  double Moire_abs_x, double Moire_abs_y, int Moire_atom_size);

   void create_magnetic_atom_list_moire_unit(std::string filename);
   void create_nm_atom_list();

#endif