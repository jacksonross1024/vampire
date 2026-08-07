#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include "positions.hpp"
#include "exchange.hpp"
#include "initialise.hpp"
#include <array>



int main(int argc, char* argv[]){

    // Defaults; flags scanned by name so order after twist/range is flexible:
    //   ./main <twist_deg> <max_range> [--dmi|--nodmi] [--bq] [--bake-only] [--config-atoms]
    //          [--candidate N] [--moire-nn-scale S] [--moire-export MODE]
    //          [J_inter DMI_inter J_twist J_prist J_intra [DMI_sub basis]
    //           [moire_nn_scale] [moire_export]]
    // Legacy: `--dmi` is a no-op placeholder (DMI on by default) so
    //   ./main 0.3 9.9 --dmi 0.085 1 1 5.5 1.45
    // still maps 0.085→J_inter … 1.45→J_intra (old fixed-index CLI).
    DMI = true;
    BQ = false;
    std::cout << " with DMI " << std::endl;
    if(argc < 2) {std::cerr << "need twist angle even if zero" << std::endl; exit(1);}

    int npos = 0; // non-flag positional index
    for(int a = 1; a < argc; a++) {
        const std::string arg = argv[a];
        if(arg == "--bq") {
            BQ = true;
            std::cout << " with biquadratic exchange (1NN intralayer, J_bq=0.44 meV)" << std::endl;
            continue;
        }
        if(arg == "--dmi") {
            DMI = true;
            std::cout << " --dmi (map DMI on; legacy placeholder)" << std::endl;
            continue;
        }
        if(arg == "--nodmi") {
            DMI = false;
            J_constant = 1.0;//eVtoJ/1000.0; //1 meV
            std::cout << " without DMI " << std::endl;
            continue;
        }
        if(arg == "--bake-only") {
            bake_only = true;
            std::cout << " bake-only: detect moiré period and exit (no exchange)" << std::endl;
            continue;
        }
        if(arg == "--config-atoms") {
            write_config_output = true;
            std::cout << " config-atoms: will write config_energy_atomic/cells.txt" << std::endl;
            continue;
        }
        if(arg == "--candidate" || arg.rfind("--candidate=", 0) == 0) {
            std::string val;
            if(arg == "--candidate") {
                if(a + 1 >= argc) { std::cerr << "--candidate needs an index\n"; exit(1); }
                val = argv[++a];
            } else {
                val = arg.substr(std::string("--candidate=").size());
            }
            moire_force_candidate = atoi(val.c_str());
            std::cout << " force AA candidate index: " << moire_force_candidate << std::endl;
            continue;
        }
        if(arg == "--moire-nn-scale" || arg.rfind("--moire-nn-scale=", 0) == 0) {
            std::string val;
            if(arg == "--moire-nn-scale") {
                if(a + 1 >= argc) { std::cerr << "--moire-nn-scale needs a value\n"; exit(1); }
                val = argv[++a];
            } else {
                val = arg.substr(std::string("--moire-nn-scale=").size());
            }
            moire_coarse_nn_tol_scale = atof(val.c_str());
            std::cout << " moire coarse nn_tol scale (x min cell width): "
                      << moire_coarse_nn_tol_scale << std::endl;
            continue;
        }
        if(arg == "--moire-export" || arg.rfind("--moire-export=", 0) == 0) {
            std::string opt;
            if(arg == "--moire-export") {
                if(a + 1 >= argc) { std::cerr << "--moire-export needs a mode\n"; exit(1); }
                opt = argv[++a];
            } else {
                opt = arg.substr(std::string("--moire-export=").size());
            }
            if (opt == "both" || opt == "all") {
                moire_coarse_write_twisted_bilayer = true;
                moire_coarse_write_twisted_double_bilayer = true;
                std::cout << " moire coarse export: twisted-bilayer + twisted-double-bilayer"
                          << std::endl;
            } else if (opt == "bilayer_only") {
                moire_coarse_write_twisted_bilayer = true;
                moire_coarse_write_twisted_double_bilayer = false;
                std::cout << " moire coarse export: twisted-bilayer only (legacy)"
                          << std::endl;
            } else if (opt == "double_only" || opt == "default") {
                moire_coarse_write_twisted_bilayer = false;
                moire_coarse_write_twisted_double_bilayer = true;
                std::cout << " moire coarse export: twisted-double-bilayer only (default)"
                          << std::endl;
            } else {
                std::cerr << " --moire-export: expected both | bilayer_only | "
                             "double_only | default, got: "
                          << opt << std::endl;
                exit(1);
            }
            continue;
        }

        if(npos == 0) {
            twist_angle = atof(argv[a]);
            std::cout << "twist angle: "<< twist_angle << std::endl;
        } else if(npos == 1) {
            max_range = atof(argv[a]);
            std::cout << "max inter exchange range: " << max_range << std::endl;
        } else if(npos == 2) {
            J_inter_scaling = atof(argv[a]);
            std::cout << " inter exchange scaling: "<< -J_inter_scaling << std::endl;
        } else if(npos == 3) {
            DMI_inter_scaling = atof(argv[a]);
            std::cout << " inter DMI scaling: " << DMI_inter_scaling << std::endl;
        } else if(npos == 4) {
            J_twist_reduction = atof(argv[a]);
            std::cout << " J twist exchange scaling: " << J_twist_reduction << std::endl;
        } else if(npos == 5) {
            J_prist_reduction = atof(argv[a]);
            std::cout << " J prist exchange scaling: " << J_prist_reduction << std::endl;
        } else if(npos == 6) {
            J_intra_reduction = atof(argv[a]);
            std::cout << " J intra exchange scaling: " << J_intra_reduction << std::endl;
        } else if(npos == 7) {
            DMI_sub_scaling = atof(argv[a]);
            std::cout << " DMI substrate exchange scaling: " << DMI_sub_scaling << std::endl;
        } else if(npos == 8) {
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
        } else if(npos == 9) {
            moire_coarse_nn_tol_scale = atof(argv[a]);
            std::cout << " moire coarse nn_tol scale (x min cell width): "
                      << moire_coarse_nn_tol_scale << std::endl;
        } else if(npos == 10) {
            std::string opt(argv[a]);
            if (opt == "both" || opt == "all") {
                moire_coarse_write_twisted_bilayer = true;
                moire_coarse_write_twisted_double_bilayer = true;
                std::cout << " moire coarse export: twisted-bilayer + twisted-double-bilayer"
                          << std::endl;
            } else if (opt == "bilayer_only") {
                moire_coarse_write_twisted_bilayer = true;
                moire_coarse_write_twisted_double_bilayer = false;
                std::cout << " moire coarse export: twisted-bilayer only (legacy)"
                          << std::endl;
            } else if (opt == "double_only" || opt == "default") {
                moire_coarse_write_twisted_bilayer = false;
                moire_coarse_write_twisted_double_bilayer = true;
                std::cout << " moire coarse export: twisted-double-bilayer only (default)"
                          << std::endl;
            } else {
                std::cerr << " moire coarse export: expected both | bilayer_only | "
                             "double_only | default, got: "
                          << opt << std::endl;
                exit(1);
            }
        }
        npos++;
    }


   system_size_x = 6500;//  25.00
   system_size_y = 4500; //4000
   number_of_unit_cells_z = 1; //2

   // Per-angle microcell scale (try 2 first; set to 1 if period detection fails).
   // twist_angle is still in degrees here (converted to rad in initialise_variables).
   // θ ≲ 0.8°: scale 1 (large moiré needs fine AA sampling). θ ≳ 0.8°: scale 2.
   // Optional override: --microcell-scale N
   int microcell_scale_override = -1;
   for(int a = 1; a < argc; a++) {
      if(std::string(argv[a]) == "--microcell-scale" && a + 1 < argc) {
         microcell_scale_override = atoi(argv[a + 1]);
      }
   }
   microcell_scale_x = 2;
   microcell_scale_y = 2;
   {
      const double th = twist_angle;
      auto near = [&](double t){ return std::fabs(th - t) < 0.05; };
      // 0°: no moiré — use a small orthogonal working window (~25 tiles of the
      // minimal rectangular cell a×√3 a) so tiling/exchange stay cheap.
      if(near(0.0)) {
         microcell_scale_x = 1; microcell_scale_y = 1;
         system_size_x = 250.0;
         system_size_y = 250.0;
      } else if(th < 0.15) {
         // 0.1°: moiré ~3400 Å — need a larger window for AA–AA edge checks.
         microcell_scale_x = 1; microcell_scale_y = 1;
         system_size_x = 9000.0;
         system_size_y = 9000.0;
      } else if(th < 0.8 || th >= 1.65) {
         // Fine sampling for large moiré (θ≲0.8) and for small cells (θ≳1.7)
         // where scale-2 AA descriptors fail edge checks.
         microcell_scale_x = 1; microcell_scale_y = 1;
      } else {
         microcell_scale_x = 2; microcell_scale_y = 2;
      }
      if(microcell_scale_override > 0) {
         microcell_scale_x = microcell_scale_override;
         microcell_scale_y = microcell_scale_override;
      }
      std::cout << "microcell_scale: " << microcell_scale_x << " x " << microcell_scale_y << std::endl;
   }

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
    // CrI3 bilayer_sliding maps:
    //   Cr{1-4}_inter_map_2.txt + Dx/Dy/Dz_inter_map_2_avg.txt, Cr{1-4}_intra_data.txt

if(!bake_only) {
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
}


    create_magnetic_atom_list_central("atom_positions.ucf");

   // calc_interactions computes unit_cell_shifts, then calls build_spin_moire_cell() to determine
   // period and crop to the primary cell; exchange constants use that PBC-applied list.
   calc_interactions();

   // Write atom list (cropped primary cell) and header after crop so UCF is consistent.
//    write_atom_positions_ucf();
//    print_header();
   if(!bake_only) print_interaction_header();

}
