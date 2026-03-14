#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "positions.hpp"
#include "exchange.hpp"
#include "initialise.hpp"
#include <array>



int main(int argc, char* argv[]){

    std::string dmi_check = "--nodmi";
    DMI = true;
     std::cout << " with DMI " << std::endl;
    if(argc < 2) {std::cerr << "need twist angle even if zero" << std::endl; exit(1);}
    for(int a = 1; a < argc; a++) {
        if(a == 1) {twist_angle = atof(argv[a]); // 1.1 
        std::cout << "twist angle: "<< twist_angle << std::endl;}
        if(a == 2) {max_range = atof(argv[a]);
        std::cout << "max inter exchange range: " << max_range << std::endl;}
        if(a == 3) {
            if(argv[a]== dmi_check) {DMI = false;
             J_constant = 1.0;//eVtoJ/1000.0; //1 meV
            std::cout << " without DMI " << std::endl;}
        }
        if(a == 4) {
            J_inter_scaling = atof(argv[a]);
            std::cout << " inter exchange scaling: "<< -J_inter_scaling << std::endl;
        }
        if(a == 5) {
            DMI_inter_scaling = atof(argv[a]);
            std::cout << " inter DMI scaling: " << DMI_inter_scaling << std::endl;
        }
        if(a == 6) {
            J_twist_reduction = atof(argv[a]);
            std::cout << " J twist exchange scaling: " << J_twist_reduction << std::endl;
        }
        if(a == 7) {
            J_prist_reduction = atof(argv[a]);
            std::cout << " J prist exchange scaling: " << J_prist_reduction << std::endl;
        }
        if(a == 8) {
            J_intra_reduction = atof(argv[a]);
            std::cout << " J intra exchange scaling: " << J_intra_reduction << std::endl;
        }
        if(a == 9) {
            DMI_sub_scaling = atof(argv[a]);
            std::cout << " DMI substrate exchange scaling: " << DMI_sub_scaling << std::endl;
        }
        if(a == 10) {
            int basis = atof(argv[a]);
            if(basis == 1) {
                DMI_sub_vector_x = 1;
                DMI_sub_vector_y = 0;
                DMI_sub_vector_z = 0;
            } else if(basis == 2) {
                DMI_sub_vector_x = 0;
                DMI_sub_vector_y = 1;
                DMI_sub_vector_z = 0;
            } else if(basis == 3) {
                DMI_sub_vector_x = 0;
                DMI_sub_vector_y = 0;
                DMI_sub_vector_z = 1;
            } else {
                std::cout << "DMI sub vector read error: exiting" << std::endl;
                exit(1);
            }
            std::cout << " DMI_sub_vector: < " << DMI_sub_vector_x << ", " << DMI_sub_vector_y << ", " << DMI_sub_vector_z << "> " << std::endl;
        }

    }


   system_size_x = 6500;//  25.00
   system_size_y = 4500; //4000
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
   // double_bilayer = false;
   // if(double_bilayer) pristine_bilayer_type = "baab";

   initialise_variables();

   twist_loction = 2*system_size_z/5 -0.01;
    std::cout << "twisting at: " << twist_loction << std::endl;
    read_in_atoms("files/atom_list_abprimebprimea_rhombic", num_atoms, atom);
    // read_in_atoms("files/atom_list_aa_rhombic", num_ref_atoms, ref_atom);
   //  read_in_dft("files/criteria.txt");
    //  read_in_atoms("files/nm_atoms", num_nm_atoms, nm_atom);
   //  read_in_exchange("files/Interpolated_J_Inter", Jinter);





#pragma omp parallel sections num_threads(4)
{  
    #pragma omp section
    {
        read_in_intra_exchanges("bilayer_sliding/Cr1_intra_data.txt", Eintra_Cr1_1NN, Eintra_Cr1_2NN, Eintra_Cr1_3NN);
           read_in_inter_exchanges("bilayer_sliding/Cr1_inter_map_2.txt",\
                           "bilayer_sliding/Cr1_Dx_inter_map_2_avg.txt",\
                           "bilayer_sliding/Cr1_Dy_inter_map_2_avg.txt",\
                           "bilayer_sliding/Cr1_Dz_inter_map_2_avg.txt", Einter_Cr1);
    }
    #pragma omp section
    {
        read_in_intra_exchanges("bilayer_sliding/Cr2_intra_data.txt", Eintra_Cr2_1NN, Eintra_Cr2_2NN, Eintra_Cr2_3NN);
            read_in_inter_exchanges("bilayer_sliding/Cr2_inter_map_2.txt",\
                           "bilayer_sliding/Cr2_Dx_inter_map_2_avg.txt",\
                           "bilayer_sliding/Cr2_Dy_inter_map_2_avg.txt",\
                           "bilayer_sliding/Cr2_Dz_inter_map_2_avg.txt", Einter_Cr2);
    }
    #pragma omp section
    {
        read_in_intra_exchanges("bilayer_sliding/Cr3_intra_data.txt", Eintra_Cr3_1NN, Eintra_Cr3_2NN, Eintra_Cr3_3NN);
            read_in_inter_exchanges("bilayer_sliding/Cr3_inter_map_2.txt",\
                           "bilayer_sliding/Cr3_Dx_inter_map_2_avg.txt",\
                           "bilayer_sliding/Cr3_Dy_inter_map_2_avg.txt",\
                           "bilayer_sliding/Cr3_Dz_inter_map_2_avg.txt", Einter_Cr3);
    }
    #pragma omp section
    {
        read_in_intra_exchanges("bilayer_sliding/Cr4_intra_data.txt", Eintra_Cr4_1NN, Eintra_Cr4_2NN, Eintra_Cr4_3NN);
            read_in_inter_exchanges("bilayer_sliding/Cr4_inter_map_2.txt",\
                           "bilayer_sliding/Cr4_Dx_inter_map_2_avg.txt",\
                           "bilayer_sliding/Cr4_Dy_inter_map_2_avg.txt",\
                           "bilayer_sliding/Cr4_Dz_inter_map_2_avg.txt", Einter_Cr4);
    }
}


    create_magnetic_atom_list_central("atom_positions.ucf");

   // calc_interactions computes unit_cell_shifts, then calls build_spin_moire_cell() to determine
   // period and crop to the primary cell; exchange constants use that PBC-applied list.
   calc_interactions();

   // Write atom list (cropped primary cell) and header after crop so UCF is consistent.
//    write_atom_positions_ucf();
//    print_header();
   print_interaction_header();

}
