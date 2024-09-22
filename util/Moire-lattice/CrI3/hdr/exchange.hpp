#ifndef EXCHANGE_HPP
#define EXCHANGE_HPP

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "positions.hpp"

   extern double nn_dist_1;
   extern double nn_dist_2;
   extern double nn_dist_3;
   extern double eVtoJ;
   extern double J_intra_1;
   extern double J_intra_2;
   extern double J_intra_3;
   extern double J_constant;
   extern int number_of_interactions;
   
   extern std::vector < std::vector < double > > Jint;
   extern std::vector < std::vector < double > > Jinter;
   extern std::vector < std::vector < double > > Jintra1;
   extern std::vector < std::vector < double > > Jintra2;

   extern double Jinter1_AB;
   extern double Jinter2_AB;
   extern double Jinter2_AB_prime;
   extern double Jintra1_AB;
   extern double Jintra2_AB;
   extern double Jintra2_ABprime;
   extern std::vector < std::vector < double > > Dx_inter;
   extern std::vector < std::vector < double > > Dy_inter;
   extern std::vector < std::vector < double > > Dz_inter;
   extern std::vector < std::vector < double > > Dx_intra;
   extern std::vector < std::vector < double > > Dy_intra;
   extern std::vector < std::vector < double > > Dz_intra;

   // extern std::vector < std::vector < std::vector<double> > > D_intra;
   // extern std::vector < std::vector < std::vector<double> > > D_inter;


   extern std::ofstream outfile4;

   void read_in_exchange(std::string filename);
   void read_in_dft(std::string filename);
   void read_in_dmi(std::string filename, std::vector < std::vector < double > > &Dx, std::vector < std::vector < double > > &Dy, std::vector < std::vector < double > > &Dz);
   void calc_in_plane_exchange(std::vector < spin > atom_list_1);
   interaction calculate_intra_Jani(spin &atom_i, spin &atom_j, double distance, double angle);
   interaction calculate_inter_Jani(spin &atom_i, spin &atom_j, double distance);
   void calc_interactions();
   void calc_out_of_plane_exchange(std::vector < spin > atom_list_1,std::vector < spin > atom_list_2);
   void print_interaction_header();

#endif
