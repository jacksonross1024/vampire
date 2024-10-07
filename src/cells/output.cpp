

//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Andrea Meo 2016. All rights reserved.
//
//   Email: am1808@york.ac.uk
//
//------------------------------------------------------------------------------


// Vampire headers
#include "cells.hpp"
#include <fstream>

// cells module internal headers
#include "internal.hpp"

namespace cells{


        std::ofstream output_file;

        void open_output_data() {
            output_file.open("microcells.txt");

         if(!output_file.is_open()) {
				std::cerr << "Fatal file directory error for microcell output" << std::endl;
			}
            return;
        }

        void output_cell_positions() {
            std::ofstream output_cell_positions;
            output_cell_positions.open("cells_positions");
            if(!output_cell_positions.is_open()) {
				std::cerr << "Fatal file directory error for microcell positions" << std::endl;
			}
            return;

        }

        void output_data() {

            if(vmpi::my_rank == 0) {
                output_file << output_counter << "\t";

                for(int cell = 0; cell < cells::mag_array_m.size(); cell++) {
                    if(cells::num_atoms_in_cell_global[cell] == 0) continue;
                    // double mod = sqrt(mag_array_x[cell]*mag_array_x[cell]+ mag_array_y[cell]*mag_array_y[cell] + mag_array_z[cell]* mag_array_z[cell]);
                    // if(mod > 0.0)  mod = 1.0/mod;   
                    output_file << mag_array_x[cell] << "\t" << mag_array_y[cell] << "\t" << mag_array_z[cell] << "\t" << mag_array_m[cell] << "\t";
                    
                }
                output_file << std::endl;
                output_counter++;
                
            }
            return;
        }
    
}