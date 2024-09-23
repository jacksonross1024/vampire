#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "positions.hpp"
#include "exchange.hpp"
#include "initialise.hpp"


int main(){

   twist_angle = 0.0; // 1.1
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
   //  read_in_dft("files/criteria.txt");
    //  read_in_atoms("files/nm_atoms", num_nm_atoms, nm_atom);
    read_in_exchange("files/interpolated_J_Inter.txt", Jinter);
    read_in_exchange("files/interpolated_J1_Intra.txt", Jintra1);
    read_in_exchange("files/interpolated_J2_Intra.txt", Jintra2);

    read_in_dmi("files/interpolated_1st_Dij_Intra.txt", Dx_intra, Dy_intra, Dz_intra);
   //  read_in_dmi("files/interpolated_1st_Dij_intra.txt", Dx_intra2, Dy_intra2, Dz_intra2);
    read_in_dmi("files/interpolated_Dij_Inter.txt", Dx_inter, Dy_inter, Dz_inter);

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
