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

// C++ standard library headers
#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>

// Vampire headers
#include "cells.hpp"
#include "random.hpp"
#include "errors.hpp"
#include "vmpi.hpp"
#include "create.hpp"
#include "micromagnetic.hpp"

#include "atoms.hpp"

// cells module internal headers
#include "internal.hpp"

namespace cells{

   //-----------------------------------------------------------------------------
   // Function for calculate magnetisation in cells
   //-----------------------------------------------------------------------------
   //int mag(const double time_from_start){
   int mag(){

      if(micromagnetic::discretisation_type != 1){
     // check calling of routine if error checking is activated
      if(err::check==true) std::cout << "cells::mag has been called" << std::endl;

      for(int i=0; i<cells::num_cells; ++i) {
         cells::mag_array_x[i] = 0.0;
         cells::mag_array_y[i] = 0.0;
         cells::mag_array_z[i] = 0.0;
         cells::mag_array_m[i] = 0.0;
      }

      #ifdef MPICF
         int num_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
      #else
         int num_local_atoms = cells::internal::num_atoms;
      #endif

      // calulate total moment in each cell
      for(int i=0;i<num_local_atoms;++i) {
         int cell = cells::atom_cell_id_array[i];
         int type = cells::internal::atom_type_array[i];
            // Consider only magnetic elements
            if(mp::material[type].non_magnetic==0){
               double mm = atoms::m_spin_array[i];
               cells::mag_array_x[cell] += atoms::x_spin_array[i]*mm;//*mus;
               cells::mag_array_y[cell] += atoms::y_spin_array[i]*mm;//*mus;
               cells::mag_array_z[cell] += atoms::z_spin_array[i]*mm;//*mus;
               cells::mag_array_m[cell] += mm;//*mus;
            }
         //}
      }

      #ifdef MPICF
      // Reduce magnetisation on all nodes
      MPI_Allreduce(MPI_IN_PLACE, &cells::mag_array_x[0],   cells::mag_array_x.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &cells::mag_array_y[0],   cells::mag_array_y.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &cells::mag_array_z[0],   cells::mag_array_z.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &cells::mag_array_m[0],   cells::mag_array_m.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
      #endif
      

      // Calculate magnetisation length and normalize
      for(int cell =0; cell < cells::mag_array_m.size(); ++cell){
         if(cells::num_atoms_in_cell[cell] == 0 ) continue;
         double msat = cells::mag_array_m[cell];
         double magm = sqrt(cells::mag_array_x[cell]*cells::mag_array_x[cell] +
                            cells::mag_array_y[cell]*cells::mag_array_y[cell] +
                            cells::mag_array_z[cell]*cells::mag_array_z[cell]);

         cells::mag_array_x[cell] = cells::mag_array_x[cell]/magm; 
         cells::mag_array_y[cell] = cells::mag_array_y[cell]/magm;             
         cells::mag_array_z[cell] = cells::mag_array_z[cell]/magm;               
         cells::mag_array_m[cell] = magm/msat;                     
      }

      if(output_microcells) {
         cells::output_data();
      }
      return EXIT_SUCCESS;
      }
   }
} // end of cells namespace
