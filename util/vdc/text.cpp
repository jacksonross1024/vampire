//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

// program header
#include "vdc.hpp"

// openmp header
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

namespace vdc{

//------------------------------------------------------------------------------
// Binary row: cx, cy, cz, sx, sy, sz (centred coords; NM spins are zero)
//------------------------------------------------------------------------------
struct txt_binary_ctx {
   size_t n_mag;
};

void fill_txt_binary_row(size_t i, double* out, void* vctx){

   const size_t n_mag = static_cast<txt_binary_ctx*>(vctx)->n_mag;

   if(i < n_mag){
      const unsigned int atom = vdc::sliced_atoms_list[i];
      out[0] = coordinates[3*atom+0] - vdc::system_centre[0];
      out[1] = coordinates[3*atom+1] - vdc::system_centre[1];
      out[2] = coordinates[3*atom+2] - vdc::system_centre[2];
      out[3] = spins[3*atom+0];
      out[4] = spins[3*atom+1];
      out[5] = spins[3*atom+2];
   }
   else{
      const unsigned int atom = vdc::sliced_nm_atoms_list[i - n_mag];
      out[0] = nm_coordinates[3*atom+0] - vdc::system_centre[0];
      out[1] = nm_coordinates[3*atom+1] - vdc::system_centre[1];
      out[2] = nm_coordinates[3*atom+2] - vdc::system_centre[2];
      out[3] = 0.0;
      out[4] = 0.0;
      out[5] = 0.0;
   }

}

//------------------------------------------------------------------------------
// Function to output spins-xxxxxxxx.txt file in plaintext format
//------------------------------------------------------------------------------
//
// Writes a single file for each snapshot ID formatted as:
//  cx   cy   cz   sx   sy   sz
//
//------------------------------------------------------------------------------
void output_txt_file(unsigned int spin_file_id){

   // Open Povray Include File
	std::stringstream txt_file_sstr;
	txt_file_sstr << "spins-";
	txt_file_sstr << std::setfill('0') << std::setw(8) << spin_file_id;
	txt_file_sstr << ".txt";
	std::string txt_file = txt_file_sstr.str();

   if(vdc::binary_dump){

      const std::string bin_name = vdc::binary_filename(txt_file);
      if(vdc::verbose) std::cout << "   Writing text file " << bin_name << "..." << std::flush;

      txt_binary_ctx ctx;
      ctx.n_mag = vdc::sliced_atoms_list.size();
      const uint64_t n_rows = static_cast<uint64_t>(ctx.n_mag + vdc::sliced_nm_atoms_list.size());

      binary_dump_spec spec;
      spec.kind = "spins-text";
      spec.notes = "Rows are sliced magnetic atoms followed by sliced non-magnetic atoms. Coordinates are relative to system_centre. NM spins are (0,0,0). Equivalent to --text/--spins .txt columns cx cy cz sx sy sz.";
      spec.columns = {
         {"cx", "Angstrom", "x relative to system_centre"},
         {"cy", "Angstrom", "y relative to system_centre"},
         {"cz", "Angstrom", "z relative to system_centre"},
         {"sx", "1", "spin x (NM: 0)"},
         {"sy", "1", "spin y (NM: 0)"},
         {"sz", "1", "spin z (NM: 0)"}
      };

      vdc::output_binary_file(bin_name, n_rows, 6, fill_txt_binary_row, &ctx, spec);

      if(vdc::verbose) std::cout << "done!" << std::endl;
      return;
   }

   // output informative message to user
   if(vdc::verbose) std::cout << "   Writing text file " << txt_file << "..." << std::flush;

   // open incfile
   std::ofstream txtfile;
   txtfile.open(txt_file.c_str());

   //---------------------------------------------------------------------------
   // parallelise stream formatting for better performance
   // step 1: parallel formatted output to stringstream in memory
   // step 2: binary write of formatted text to output file (awesomely fast!)
   //---------------------------------------------------------------------------
   #pragma omp parallel
   {

      std::stringstream otext;

      // write to output text stream in parallel
      #pragma omp for
      for(unsigned int atom = 0; atom < vdc::num_atoms; atom++){

         // format text for plain text file
         otext << coordinates[3*atom+0]-vdc::system_centre[0] << "\t" << coordinates[3*atom+1]-vdc::system_centre[1] << "\t" << coordinates[3*atom+2]-vdc::system_centre[2] << "\t" <<
                  spins[3*atom+0] << "\t" << spins[3*atom+1] << "\t" << spins[3*atom+2] << "\n";

      } // end of parallel for

      // force each thread to write to file in order
      #pragma omp critical
      txtfile << otext.str();

   } // end of parallel region

   //---------------------------------------------------------------------------
   // write non-magnetic atoms to txt file
   //---------------------------------------------------------------------------
   // parallelise stream formatting for better performance
   //---------------------------------------------------------------------------
   #pragma omp parallel
   {

      std::stringstream otext;

      // write to output text stream in parallel
      #pragma omp for
      for(unsigned int atom = 0; atom < vdc::num_nm_atoms; atom++){

         // format text for text file
         otext << nm_coordinates[3*atom+0]-vdc::system_centre[0] << "\t" << nm_coordinates[3*atom+1]-vdc::system_centre[1] << "\t" << nm_coordinates[3*atom+2]-vdc::system_centre[2] << "\t" <<
                  0.0 << "\t" << 0.0 << "\t" << 0.0 << "\n";

      } // end of parallel for

      // force each thread to write to file in order
      #pragma omp critical
      txtfile << otext.str();

   } // end of parallel region

   // flush data to include file and close
   txtfile << std::flush;
   txtfile.close();

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   return;

}

}
