//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>

// Vampire headers
#include "spintorque.hpp"
#include "vmpi.hpp"
#include "vio.hpp"
#include "sim.hpp"
#include "material.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace st{
   namespace internal{
      //-----------------------------------------------------------------------------
      // Function to output base microcell properties
      //-----------------------------------------------------------------------------
      void output_base_microcell_data(){

         using st::internal::beta_cond;
         using st::internal::beta_diff;
         using st::internal::sa_infinity;
         using st::internal::lambda_sdl;
         using st::internal::pos;

         const int num_cells = beta_cond.size();

         // only output on root process
         if(vmpi::my_rank==0){
            if(sim::time%(ST_output_rate) ==0){
               zlog << zTs() << "Outputting ST base microcell data" << std::endl;
               std::ofstream ofile;
               ofile.open("spin-acc/st-microcells-base.cfg");
               ofile << num_cells << std::endl;
               ofile << num_stacks << std::endl;
               for(int cell=0; cell < num_cells; ++cell){
                    if( (st::internal::cell_stack_index[cell]-1)%3 == 0) continue;
	               ofile << cell_stack_index[cell] << "\t" << "\t" << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2];
	               ofile << "\t" << beta_cond[cell] << "\t" << beta_diff[cell] << "\t" << sa_infinity[cell] << "\t" << lambda_sdl[cell] << std::endl;
               }

		         ofile.close();
	         }

         }

      return;
      }

      //-----------------------------------------------------------------------------
      // Function to output base microcell properties
      //-----------------------------------------------------------------------------
      void output_microcell_data(){

      
         // only output on root process
         #ifdef MPICF 
           MPI_Barrier(MPI_COMM_WORLD);

            if(sim::time%(ST_output_rate) ==0){
             

                  using st::internal::m;
         using st::internal::sa;
         using st::internal::j;
         using st::internal::coeff_ast;
         using st::internal::coeff_nast;
         using st::internal::ast;
         using st::internal::nast;
         using st::internal::cell_natom;

         const int size = m.size();
         const int num_cells = size/3;
         

               // determine file name
               std::stringstream filename;
               filename << "spin-acc/" << config_file_counter;
          
               // #ifdef MPICF 
                 MPI_Reduce(&sa[0], &sa_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                  MPI_Reduce(&sa_sot[0], &sa_sot_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                  //MPI_Reduce(&m[0], &m_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                  MPI_Reduce(&j[0], &j_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                  MPI_Reduce(&coeff_ast[0], &coeff_ast_sum[0], num_cells, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                  MPI_Reduce(&coeff_nast[0], &coeff_nast_sum[0], num_cells, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                  MPI_Reduce(&ast[0], &ast_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                  MPI_Reduce(&nast[0], &nast_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                  MPI_Reduce(&total_ST[0], &total_ST_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                  MPI_Reduce(&cell_natom[0], &cell_natom_sum[0], num_cells, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
                
               // #endif 
               // zlog << zTs() << "Outputting ST microcell data " << filename.str() << std::endl;

              
            if(vmpi::my_rank == 0) {
                     std::ofstream ofile;
               ofile.open(std::string(filename.str()).c_str());
	          //  ofile<<"Time:"<< "\t" << sim::time*mp::dt_SI<< std::endl;
             ofile << "pos_x \t pos_y \t pos_z \t m_x \t m_y \t m_z \t spin_acc_x \t spin_acc_y \t spin_acc_z \t " << \
             "j_x \t j_y \t j_z \t ast_x \t ast_y \t ast_z \t nast_x \t nast_y \t nast_z \t torque_x \t torque_y \t torque_z \t num_atom" << std::endl;
               for(int cell=0; cell<num_cells; ++cell){
                  //   if( (st::internal::cell_stack_index[cell]-1)%3 == 0) continue;
                  if(cell_natom[cell] == 0) continue;
                  ofile << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] << "\t";
                  ofile << m[3*cell+0] << "\t" << m[3*cell+1] << "\t" << m[3*cell+2] << "\t";
                  ofile << sa_sum[3*cell+0] << "\t" << sa_sum[3*cell+1] << "\t" << sa_sum[3*cell+2] << "\t";
                  ofile << sa_sot_sum[3*cell+0] << "\t" << sa_sot_sum[3*cell+1] << "\t" << sa_sot_sum[3*cell+2] << "\t";
                  ofile << j_sum[3*cell+0] << "\t" << j_sum[3*cell+1] << "\t" << j_sum[3*cell+2] << "\t";
                  ofile << coeff_ast_sum[cell] << "\t";
                  ofile << coeff_nast_sum[cell] << "\t";
                  ofile << ast_sum[3*cell+0] << "\t" << ast_sum[3*cell+1] << "\t" << ast_sum[3*cell+2] << "\t";
                  ofile << nast_sum[3*cell+0] << "\t" << nast_sum[3*cell+1] << "\t" << nast_sum[3*cell+2] << "\t";
                  ofile << total_ST_sum[3*cell+0] << "\t" << total_ST_sum[3*cell+1] << "\t" << total_ST_sum[3*cell+2];
                  ofile << "\t" << cell_natom_sum[cell] << "\n";

                  sa_sum[3*cell+0] = sa_sum[3*cell+1] = sa_sum[3*cell+2] = 0.0;
                  sa_sot_sum[3*cell+0] = sa_sot_sum[3*cell+1] = sa_sot_sum[3*cell+2] = 0.0;
                  j_sum[3*cell+0] = j_sum[3*cell+1] = j_sum[3*cell+2] = 0.0;
                  coeff_ast_sum[cell] = 0.0;
                  coeff_nast_sum[cell]  = 0.0;
                  ast_sum[3*cell+0] = ast_sum[3*cell+1] = ast_sum[3*cell+2] = 0.0;
                  nast_sum[3*cell+0] = nast_sum[3*cell+1] = nast_sum[3*cell+2] = 0.0;
                  total_ST_sum[3*cell+0] = total_ST_sum[3*cell+1] = total_ST_sum[3*cell+2] = 0.0;
                  cell_natom_sum[cell] = 0;
                  // ofile << std::endl;
                  // ofile << spin_torque[3*cell+0] << "\t" << spin_torque[3*cell+1] << "\t" << spin_torque[3*cell+2] << std::endl;
               }

            ofile.close();
       
            // update config_file_counter
          

           }
         config_file_counter++;
         }

           
  #endif 
         return;
      }

   } // end of internal namespace
} // end of st namespace
