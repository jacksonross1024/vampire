#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "positions.hpp"
#include "exchange.hpp"
#include "initialise.hpp"


int main(){

   twist_angle = 1.41; // 1.1
   system_size_x = 1000; //4000
   system_size_y = 1000; //4000
   number_of_unit_cells_z = 1; //2

   // set up new material constants
   dmi12 = 1.0; // DMI constant between layers 1-2
   dmi23 = 1.0; // DMI constant between layers 2-3
   dmi34 = 1.0; // DMI constant between layers 3-4

   dmi_decay = 1.0; // distance-dependent DMI

   exchange12 = 1.0; // exchange constant between layers 1-2
   exchange23 = 1.0; // exchange constant between layers 2-3
   exchange34 = 1.0; // exchange constant between layers 3-4

   separation = 0.0; // distance between layers 2-3

   initialise_variables();

   twist_loction = 2*system_size_z/5 -0.01;
    std::cout << "twisting at: " << twist_loction << std::endl;
    read_in_atoms("files/atom_list_aa_rhombic", num_atoms, atom);
    //  read_in_atoms("files/nm_atoms", num_nm_atoms, nm_atom);
    read_in_exchange("files/interpolated_array-inter-AA", Jinter);
    read_in_exchange("Maps/interpolated_J1_intra-AA.txt", Jintra1);
    read_in_exchange("Maps/interpolated_J2_intra-AA.txt", Jintra2);

    read_in_dmi("Maps/interpolated_1st_Dij_intra", Dx_intra, Dy_intra, Dz_intra);
    read_in_dmi("Maps/interpolated_Dij_inter", Dx_inter, Dy_inter, Dz_inter);


    std::ofstream Jinter_out;
    std::ofstream Jintra1_out;
    std::ofstream Jintra2_out;
    std::ofstream Dinter_out;
    std::ofstream Dintra1_out;
    // std::ofstream Dintra2_out;
    Jinter_out.open("Jinter_interpolation_out.txt");
    Jintra1_out.open("Jintra1_interpolation_out.txt");
    Jintra2_out.open("Jintra2_interpolation_out.txt");
    Dinter_out.open("Dinter_interpolation_out.txt");
    Dintra1_out.open("Dintra1_interpolation_out.txt");
    // Dintra2_out.open("Dintra2_interpolation_out.txt");
    for(int i = 0; i < 201; i++) {
        for(int j = 0; j < 201; j++) {
            Jinter_out << i << ", " << j << ", " << Jinter[i][j] << "\n";
            Jintra1_out << i << ", " << j << ", " << Jintra1[i][j] << "\n";
            Jintra2_out << i << ", " << j << ", " << Jintra2[i][j] << "\n";
            Dinter_out << i << ", " << j << ", " << Dx_inter[i][j] << ", " << Dy_inter[i][j] << ", " << Dz_inter[i][j] << ", " << sqrt(Dx_inter[i][j]*Dx_inter[i][j] \
                                                                                                                                     + Dy_inter[i][j]*Dy_inter[i][j] \
                                                                                                                                     + Dz_inter[i][j]*Dz_inter[i][j]) << "\n";
            Dintra1_out << i << ", " << j << ", " << Dx_intra[i][j] << ", " << Dy_intra[i][j] << ", " << Dz_intra[i][j] <<  ", " << sqrt(Dx_intra[i][j]*Dx_intra[i][j] \
                                                                                                                                     + Dy_intra[i][j]*Dy_intra[i][j] \
                                                                                                                                     + Dz_intra[i][j]*Dz_intra[i][j]) << "\n";
            // Dintra2_out << i << ", " << j << ", " << Jintra2[i][j] << "\n";
        }
    }
    Jinter_out.close();
    Jintra1_out.close();
    Jintra2_out.close();
    Dinter_out.close();
    Dintra1_out.close();

    create_magnetic_atom_list("atom_positions.ucf");
    //  create_nm_atom_list();

    print_header();
    //  calc_in_plane_exchange(row1);
    //  calc_in_plane_exchange(row2);
    //  calc_in_plane_exchange(row3);
    //  calc_in_plane_exchange(row4);

    //  // The order of these exchange calculations is important, as dx,dy shifts are NOT symmetric
    //  calc_out_of_plane_exchange(row4,row3);
    //  calc_out_of_plane_exchange(row3,row2);
    //  calc_out_of_plane_exchange(row2,row1);
    calc_interactions();
    print_interaction_header();

 }
