//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

// program header
#include "vdc.hpp"

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

namespace vdc{

//------------------------------------------------------------------------------
// Binary row: type_id, x, y, z (lab-frame coords, same as crystal.xyz)
//------------------------------------------------------------------------------
struct xyz_binary_ctx {
   size_t n_mag;
};

void fill_xyz_binary_row(size_t i, double* out, void* vctx){

   const size_t n_mag = static_cast<xyz_binary_ctx*>(vctx)->n_mag;

   if(i < n_mag){
      const unsigned int atom = vdc::sliced_atoms_list[i];
      out[0] = static_cast<double>(vdc::type[atom]);
      out[1] = vdc::coordinates[3*atom + 0];
      out[2] = vdc::coordinates[3*atom + 1];
      out[3] = vdc::coordinates[3*atom + 2];
   }
   else{
      const unsigned int atom = vdc::sliced_nm_atoms_list[i - n_mag];
      out[0] = static_cast<double>(vdc::nm_type[atom]);
      out[1] = vdc::nm_coordinates[3*atom + 0];
      out[2] = vdc::nm_coordinates[3*atom + 1];
      out[3] = vdc::nm_coordinates[3*atom + 2];
   }

}

//------------------------------------------------------------------------------
// Function to output crystal.xyz file compatible with rasmol
//------------------------------------------------------------------------------
void output_xyz_file(){

   if(vdc::binary_dump){

      const std::string bin_name = vdc::binary_filename("crystal.xyz");
      if(vdc::verbose) std::cout << "Writing xyz file... " << std::flush;

      xyz_binary_ctx ctx;
      ctx.n_mag = vdc::sliced_atoms_list.size();
      const uint64_t n_rows = static_cast<uint64_t>(ctx.n_mag + vdc::sliced_nm_atoms_list.size());

      binary_dump_spec spec;
      spec.kind = "xyz";
      spec.notes = "Lab-frame coordinates (not centred). type is 0-based index into materials[]. Magnetic atoms first, then non-magnetic. Element symbols are not stored; map type -> materials[type].element.";
      spec.columns = {
         {"type", "1", "0-based material index"},
         {"x", "Angstrom", "lab-frame x"},
         {"y", "Angstrom", "lab-frame y"},
         {"z", "Angstrom", "lab-frame z"}
      };

      vdc::output_binary_file(bin_name, n_rows, 4, fill_xyz_binary_row, &ctx, spec);

      if(vdc::verbose) std::cout << "done!" << std::endl;
      return;
   }

   // output informative message to user
   if(vdc::verbose) std::cout << "Writing xyz file... " << std::flush;

   // output xyz file
   std::ofstream ofile;
   ofile.open("crystal.xyz");

   // output number of atoms
   ofile << vdc::sliced_atoms_list.size() + vdc::sliced_nm_atoms_list.size() << "\n\n";

   #pragma omp parallel
   {

      std::stringstream otext;

      // write magnetic atoms to output text stream in parallel
      #pragma omp for
      for(size_t i=0; i < vdc::sliced_atoms_list.size(); i++){

         // get atom ID
         unsigned int atom = vdc::sliced_atoms_list[i];

         // get atom type
         int type_id = vdc::type[atom];

         otext << materials[type_id].element << "\t" <<
                  vdc::coordinates[3*atom + 0] << "\t" <<
                  vdc::coordinates[3*atom + 1] << "\t" <<
                  vdc::coordinates[3*atom + 2] << "\n";

      } // end of parallel for

      // write non-magnetic atoms
      #pragma omp for
      for(size_t i=0; i < vdc::sliced_nm_atoms_list.size(); i++){

         // get atom ID
         unsigned int atom = vdc::sliced_nm_atoms_list[i];

         // get atom type
         int type_id = vdc::nm_type[atom];

         otext << materials[type_id].element << "\t" <<
                  vdc::nm_coordinates[3*atom + 0] << "\t" <<
                  vdc::nm_coordinates[3*atom + 1] << "\t" <<
                  vdc::nm_coordinates[3*atom + 2] << "\n";

      } // end of parallel for

      // force each thread to write to file in order
      #pragma omp critical
      ofile << otext.str();

   } // end of parallel region

   ofile << std::flush;
   ofile.close();

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   return;

}

}
