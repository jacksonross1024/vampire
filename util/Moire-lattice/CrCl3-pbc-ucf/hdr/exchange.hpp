#ifndef EXCHANGE_HPP
#define EXCHANGE_HPP

#include <iostream>
#include <vector>
#include <array>
#include <sstream>
#include <fstream>
#include <cmath>
#include "positions.hpp"

#include <unistd.h>
#include <omp.h>

   extern double nn_dist_1;
   extern double nn_dist_2;
   extern double nn_dist_3;
   extern double eVtoJ;
   extern double J_intra_1;
   extern double J_intra_2;
   extern double J_intra_3;
   extern double J_constant;
   extern double max_range;
   extern uint64_t number_of_interactions;
   extern uint64_t number_of_bq_interactions;
   extern double J_bq;
   extern double Jz_1NN_aniso;
   
   extern std::vector < std::vector < double > > Jint;
   extern std::vector < std::vector < double > > Jinter;
   extern std::vector < std::vector < double > > Einter_Cr1;
   extern std::vector < std::vector < double > > Einter_Cr2;
   extern std::vector < std::vector < double > > Einter_Cr3;
   extern std::vector < std::vector < double > > Einter_Cr4;

   extern std::vector < std::vector < double  > > Eintra_Cr1_1NN;
   extern std::vector < std::vector < double  > > Eintra_Cr2_1NN;
   extern std::vector < std::vector < double  > > Eintra_Cr3_1NN;
   extern std::vector < std::vector < double  > > Eintra_Cr4_1NN;

   extern std::vector < std::vector < double  > > Eintra_Cr1_2NN;
   extern std::vector < std::vector < double  > > Eintra_Cr2_2NN;
   extern std::vector < std::vector < double  > > Eintra_Cr3_2NN;
   extern std::vector < std::vector < double  > > Eintra_Cr4_2NN;

   extern std::vector < std::vector < double  > > Eintra_Cr1_3NN;
   extern std::vector < std::vector < double  > > Eintra_Cr2_3NN;
   extern std::vector < std::vector < double  > > Eintra_Cr3_3NN;
   extern std::vector < std::vector < double  > > Eintra_Cr4_3NN;

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
   extern std::ofstream outfile_bq;

   // Average both endpoint species maps so i→j and j→i share J / opposite DMI
   // even when Cr1–4 tables differ by imperfect reflections.
   std::array<float,4> match_inter_exchange(int atomi_id, int nn_id, double dx, double dy, double dr,
      std::vector<std::vector<double> > &Eij_i, std::vector<std::vector<double> > &Eij_j);
   // Intra: dual species maps (A/B × bottom/top). Eij_c = central atom's Cr map,
   // Eij_n = neighbour's. 1NN/3NN |θ| bins share θ and θ+π; partner D must come
   // from the other sublattice table (Cr2=R_y(Cr1), Cr3=−Cr1, Cr4=−R_y).
   std::array<float,4> match_intra1_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom,
      std::vector<std::vector<double > > &Eij_c, std::vector<std::vector<double > > &Eij_n);
   std::array<float,4> match_intra2_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom,
      std::vector<std::vector<double > > &Eij_c, std::vector<std::vector<double > > &Eij_n);
   std::array<float,4> match_intra3_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom,
      std::vector<std::vector<double > > &Eij_c, std::vector<std::vector<double > > &Eij_n);
   void spin_config_energy(spin & atom_i, double dr2, spin & atom_j, std::array<float, 4> &exchange, std::vector<std::vector<std::vector<double> > > & local_config_energy);

   void calc_interactions();

   // Compute unit_cell_shifts from an atom list (same algorithm as calc_interactions lines 180-375).
   // Output is written to out_shifts; caller must ensure it has correct dimensions (microcell_Nx+1 x microcell_Ny+1 x 3).
   void compute_unit_cell_shifts_from_atoms(const std::vector<spin>& atoms,
      std::vector<std::vector<std::vector<int>>>& out_shifts);

   void print_interaction_header();

#endif
