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
               ofile << num_stacks_y << std::endl;
               for(int cell=0; cell < num_cells; ++cell){
                  // if(cell_natom[cell] == 0) continue;
                  if(sot_sa) {
                        ofile << cell_stack_index[cell] << "\t" << "\t" << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] \
                        << "\t" << beta_cond[cell] << "\t" << beta_diff[cell] << "\t" << sa_infinity[cell] << "\t" << lambda_sdl[cell] << "\t" << \
                        st::internal::a[cell] << "\t" << st::internal::b[cell] << "\t" \
                        << "\t" << st::internal::sot_beta_cond[cell] << "\t" << st::internal::sot_beta_diff[cell] << "\t" << st::internal::sot_sa_infinity[cell] << "\t" << st::internal::sot_lambda_sdl[cell] << "\t" \
                        << st::internal::sot_a[cell] << "\t" << st::internal::sot_b[cell] << "\t" << st::internal::spin_acc_sign[cell] << "\t" << st::internal::sot_sa_source[cell] << "\t" << st::internal::cell_natom[cell] << std::endl;
                  } else {
                        ofile << cell_stack_index[cell] << "\t" << "\t" << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] \
                        << "\t" << beta_cond[cell] << "\t" << beta_diff[cell] << "\t" << sa_infinity[cell] << "\t" << lambda_sdl[cell] << "\t" << 
                        st::internal::a[cell] << "\t" << st::internal::b[cell] << "\t" <<  std::endl;
                     }
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

               const int size = sa_final.size();
               const int num_cells = size/3;
               // determine file name
               std::stringstream filename;
               filename << "spin-acc/" << config_file_counter;
          
               MPI_Reduce(&sa_final[0], &sa_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&j_final_up_x[0], &j_final_up_x_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&j_final_up_y[0], &j_final_up_y_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&j_final_down_y[0], &j_final_down_y_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&coeff_ast[0], &coeff_ast_sum[0], num_cells, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&coeff_nast[0], &coeff_nast_sum[0], num_cells, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&ast[0], &ast_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&nast[0], &nast_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&total_ST[0], &total_ST_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                          
            if(vmpi::my_rank == 0) {
               std::ofstream ofile;
               ofile.open(std::string(filename.str()).c_str());
	          
               ofile << "pos_x \t pos_y \t pos_z \t m_x \t m_y \t m_z \t spin_acc_x \t spin_acc_y \t spin_acc_z \t " << \
               "jup_x \t jup_y \t jup_z \t jdown_x \t jdown_y \t jdown_z \t ast_x \t ast_y \t ast_z \t nast_x \t nast_y \t nast_z \t torque_x \t torque_y \t torque_z \t num_atom" << std::endl;
               for(int cell=0; cell<num_cells; ++cell){
                  //   if( (st::internal::cell_stack_index[cell]-1)%3 == 0) continue;
                  if(cell_natom[cell] == 0) continue;
                  double mag = sqrt(m[3*cell+0]*m[3*cell+0] + m[3*cell+1]*m[3*cell+1] + m[3*cell+2]*m[3*cell+2]);
                  mag = (mag == 0.0) ? 0.0: 1/mag;
                  ofile << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] << "\t";
                  ofile << m[3*cell+0]*mag << "\t" << m[3*cell+1]*mag << "\t" << m[3*cell+2]*mag << "\t";
                 // if(st::internal::sot_check) ofile << (sa_sum[3*cell+0]-sa_infinity[cell]*m[3*cell]*mag)/sa_infinity[cell] << "\t" << (sa_sum[3*cell+1]-sa_infinity[cell]*m[3*cell+1]*mag)/sa_infinity[cell] << "\t" << (sa_sum[3*cell+2]-m[3*cell+2]*mag)/sa_infinity[cell] << "\t";
                  //else 
                  ofile << std::fixed << std::setprecision(10);
                  ofile << sa_sum[3*cell+0] << "\t" << sa_sum[3*cell+1] << "\t" << sa_sum[3*cell+2] << "\t";
                  ofile << std::setprecision(6);
                  ofile << j_final_up_x_sum[3*cell+0] << "\t" << j_final_up_x_sum[3*cell+1] << "\t" << j_final_up_x_sum[3*cell+2] << "\t";
                  ofile << j_final_up_y_sum[3*cell+0] << "\t" << j_final_up_y_sum[3*cell+1] << "\t" << j_final_up_y_sum[3*cell+2] << "\t";
                  ofile << j_final_down_y_sum[3*cell+0] << "\t" << j_final_down_y_sum[3*cell+1] << "\t" << j_final_down_y_sum[3*cell+2] << "\t";
                  ofile << coeff_ast_sum[cell] << "\t";
                  ofile << coeff_nast_sum[cell] << "\t";
                  ofile << ast_sum[3*cell+0] << "\t" << ast_sum[3*cell+1] << "\t" << ast_sum[3*cell+2] << "\t";
                  ofile << nast_sum[3*cell+0] << "\t" << nast_sum[3*cell+1] << "\t" << nast_sum[3*cell+2] << "\t";
                  ofile << total_ST_sum[3*cell+0] << "\t" << total_ST_sum[3*cell+1] << "\t" << total_ST_sum[3*cell+2] << "\t" << sqrt(total_ST_sum[3*cell+0]*total_ST_sum[3*cell+0] + total_ST_sum[3*cell+1]*total_ST_sum[3*cell+1] + total_ST_sum[3*cell+2]*total_ST_sum[3*cell+2]);
                  ofile << "\t" << cell_natom[cell] << "\n";

                  sa_sum[3*cell+0] = 0.0;
                  sa_sum[3*cell+1] = 0.0;
                  sa_sum[3*cell+2] = 0.0;
                  j_final_up_x_sum[3*cell+0] = j_final_up_x_sum[3*cell+1] = j_final_up_x_sum[3*cell+2] = 0.0;
                  j_final_up_y_sum[3*cell+0] = j_final_up_y_sum[3*cell+1] = j_final_up_y_sum[3*cell+2] = 0.0;
                  j_final_down_y_sum[3*cell+0] = j_final_down_y_sum[3*cell+1] = j_final_down_y_sum[3*cell+2] = 0.0;
                  coeff_ast_sum[cell] = 0.0;
                  coeff_nast_sum[cell]  = 0.0;
                  ast_sum[3*cell+0] = ast_sum[3*cell+1] = ast_sum[3*cell+2] = 0.0;
                  nast_sum[3*cell+0] = nast_sum[3*cell+1] = nast_sum[3*cell+2] = 0.0;
                  total_ST_sum[3*cell+0] = total_ST_sum[3*cell+1] = total_ST_sum[3*cell+2] = 0.0;
           
               }

            ofile.close();
       
            // update config_file_counter
           }
         config_file_counter++;
         }
         #else
      if(sim::time%(ST_output_rate) ==0){ 

      // using namespace st::internal;
         const int size = m.size();
         const int num_cells = size/3;

            // determine file name
            std::stringstream filename;
            filename << "spin-acc/" << config_file_counter;
         
         
                  std::ofstream ofile;
            ofile.open(std::string(filename.str()).c_str());
            //  ofile<<"Time:"<< "\t" << sim::time*mp::dt_SI<< std::endl;
            ofile << "pos_x \t pos_y \t pos_z \t m_x \t m_y \t m_z \t spin_acc_x \t spin_acc_y \t spin_acc_z \t " << \
            "j_x \t j_y \t j_z \t ast_x \t ast_y \t ast_z \t nast_x \t nast_y \t nast_z \t torque_x \t torque_y \t torque_z \t num_atom" << std::endl;
            for(int cell=0; cell<num_cells; ++cell){
               //   if( (st::internal::cell_stack_index[cell]-1)%3 == 0) continue;
               if(cell_natom[cell] == 0) continue;
               double mag = sqrt(m[3*cell+0]*m[3*cell+0] + m[3*cell+1]*m[3*cell+1] + m[3*cell+2]*m[3*cell+2]);
               mag = (mag == 0.0) ? 0.0: 1/mag;
               ofile << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] << "\t";
               ofile << m[3*cell+0] << "\t" << m[3*cell+1] << "\t" << m[3*cell+2] << "\t";
               if(st::internal::sot_check) ofile << (sa_final[3*cell+0]-sa_infinity[cell]*m[3*cell]*mag)/sa_infinity[cell] << "\t" << (sa_final[3*cell+1]-sa_infinity[cell]*m[3*cell+1]*mag)/sa_infinity[cell] << "\t" << (sa_final[3*cell+2]-m[3*cell+2]*mag)/sa_infinity[cell] << "\t";
               else ofile << sa_final[3*cell+0] << "\t" << sa_final[3*cell+1] << "\t" << sa_final[3*cell+2] << "\t";
               ofile << j_final_up_x[3*cell+0] << "\t" << j_final_up_x[3*cell+1] << "\t" << j_final_up_x[3*cell+2] << "\t";
               ofile << j_final_up_y[3*cell+0] << "\t" << j_final_up_y[3*cell+1] << "\t" << j_final_up_y[3*cell+2] << "\t";
               ofile << j_final_down_y[3*cell+0] << "\t" << j_final_down_y[3*cell+1] << "\t" << j_final_down_y[3*cell+2] << "\t";
               ofile << coeff_ast[cell] << "\t";
               ofile << coeff_nast[cell] << "\t";
               ofile << ast[3*cell+0] << "\t" << ast[3*cell+1] << "\t" << ast[3*cell+2] << "\t";
               ofile << nast[3*cell+0] << "\t" << nast[3*cell+1] << "\t" << nast[3*cell+2] << "\t";
               ofile << total_ST[3*cell+0] << "\t" << total_ST[3*cell+1] << "\t" << total_ST[3*cell+2];
               ofile << "\t" << cell_natom[cell] << "\n";
               
            }

            ofile.close();
       
            // update config_file_counter
         config_file_counter++;
      }
           
   #endif 
   
   return;
   }

   } // end of internal namespace
} // end of st namespace
