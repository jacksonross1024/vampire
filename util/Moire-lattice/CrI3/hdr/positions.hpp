#ifndef POSITION_HPP
#define POSITION_HPP

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>


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
         double S;
         int id;
         int l_id;
         int unit_x;
         int unit_y;
         int dx = 0;
         int dy = 0;
         int inter = 0;
         int intra1 = 0;
         int intra2 = 0;
   };

   class interaction {
      public:
         double xx= 10.0;
         double xy= 0.0;
         double xz= 0.0;
         double yx= 0.0;
         double yy= 0.0;
         double yz= 0.0;
         double zx= 0.0;
         double zy= 0.0;
         double zz= 0.0;
   };

   extern std::vector < spin > atom;
   extern std::vector < spin > nm_atom;
   extern std::vector < spin > row1;
   extern std::vector < spin > row2;
   extern std::vector < spin > row3;
   extern std::vector < spin > row4;
   extern std::vector < spin > all_m_atoms;
   extern std::vector < spin > all_nm_atoms;
   extern std::vector < std::vector < std::vector <int> > > unit_cell_shifts;
   void read_in_atoms(std::string filename, int n_atoms, std::vector <spin > &atom2);
   void read_in_exchange(std::string filename, std::vector<std::vector<double> > &Jij);
   
   void read_in_dmi(std::string filename, std::vector < std::vector < double > > &Dx, std::vector < std::vector < double > > &Dy, std::vector < std::vector < double > > &Dz, bool reflect);

   void print_header();
   void create_magnetic_atom_list(std::string filename);
   void create_nm_atom_list();

#endif