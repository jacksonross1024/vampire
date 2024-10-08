#ifndef EXCHANGE_HPP
#define EXCHANGE_HPP

#include <iostream>
#include <vector>
#include <array>
#include <sstream>
#include <fstream>
#include <cmath>
#include "positions.hpp"
// ##include <omp.h>

   extern double nn_dist_1;
   extern double nn_dist_2;
   extern double nn_dist_3;
   extern double eVtoJ;
   extern double J_intra_1;
   extern double J_intra_2;
   extern double J_intra_3;
   extern double J_constant;
   extern double max_range;
   extern int number_of_interactions;
   
   extern std::vector < std::vector < double > > Jint;
   extern std::vector < std::vector < double > > Jinter;
   extern std::vector < std::vector < double > > Einter_Cr1;
   extern std::vector < std::vector < double > > Einter_Cr2;
   extern std::vector < std::vector < double > > Einter_Cr3;
   extern std::vector < std::vector < double > > Einter_Cr4;

   extern std::vector < std::vector < std::vector<std::vector<double > > > > Eintra_Cr1_1NN;
   extern std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr2_1NN;
   extern std::vector < std::vector < std::vector<std::vector<double > > > > Eintra_Cr3_1NN;
   extern std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr4_1NN;

   extern std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr1_2NN;
   extern std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr2_2NN;
   extern std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr3_2NN;
   extern std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr4_2NN;

   extern std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr1_3NN;
   extern std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr2_3NN;
   extern std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr3_3NN;
   extern std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr4_3NN;

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
   extern std::vector < std::vector < double > > Dx_intra2;
   extern std::vector < std::vector < double > > Dy_intra2;
   extern std::vector < std::vector < double > > Dz_intra2;

   // extern std::vector < std::vector < std::vector<double> > > D_intra;
   // extern std::vector < std::vector < std::vector<double> > > D_inter;


   extern std::ofstream outfile4;

   void read_in_exchange(std::string filename);
   void read_in_dft(std::string filename);
   void read_in_dmi(std::string filename, std::vector < std::vector < double > > &Dx, std::vector < std::vector < double > > &Dy, std::vector < std::vector < double > > &Dz);
   void calc_in_plane_exchange(std::vector < spin > atom_list_1);

   std::vector<double> match_inter_exchange(double dx, double dy, std::vector<std::vector<double> > &Eij);
   std::vector<double> match_intra_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom, std::vector<std::vector< std::vector< std::vector<double> >  > > &Eij);
   std::vector<double> match_intra2_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom, std::vector<std::vector< std::vector< std::vector<double> >  > > &Eij);
   interaction calculate_intra_Jani(spin &atom_i, spin &atom_j, double distance, double angle);
   interaction calculate_inter_Jani(spin &atom_i, spin &atom_j, double distance, double angle);
   void calc_interactions();

   void calc_out_of_plane_exchange(std::vector < spin > atom_list_1,std::vector < spin > atom_list_2);
   void print_interaction_header();

#endif
