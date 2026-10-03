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
#include "ltmp.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace st{
   namespace internal{
      //-----------------------------------------------------------------------------
      // Function to output 1D spin currents data
      //-----------------------------------------------------------------------------
      void output_sc1d_data(){
         if(!sc1d_enable) return;

         #ifdef MPICF
            MPI_Barrier(MPI_COMM_WORLD);
         #endif

         // relax: before the LLG, sim::time stays 0. eq: LLG equilibration.
         // drive: sim::time past equilibration-time-steps. The index is the
         // step count inside that stage. ST-output-rate applies in every stage.
         const int rate = std::max(ST_output_rate, 1);
         std::string stage;
         unsigned long index = 0;
         if(sc1d_step_counter <= sc1d_relax_steps){
            stage = "relax";
            index = sc1d_step_counter;
         } else if(sim::time <= sim::equilibration_time){
            stage = "eq";
            index = static_cast<unsigned long>(sim::time);
         } else {
            stage = "drive";
            index = static_cast<unsigned long>(sim::time - sim::equilibration_time);
         }
         const bool due = (index > 0UL) && (index % static_cast<unsigned long>(rate) == 0UL);
         if(due){
            const int size = ns_final.size(); // number of cells
            const int size3 = js_final.size();

            // determine file name
            std::stringstream filename;
            filename << "spin-acc/" << stage << "-"
                     << std::setw(7) << std::setfill('0') << index;

            const int size_sa = sa_final.size(); // spin accumulation size (3*cells)
            #ifdef MPICF
               // Cell fields were summed onto every rank in broadcast_cell_fields.
               // Only rank 0 writes the file.
               if(vmpi::my_rank == 0) {
                  ns_sum = ns_final;
                  ne_sum = ne_final;
                  jc_sum = jc_final;
                  js_sum = js_final;
                  sa_sum = sa_final;
               }
            #else
               // serial copy
               ns_sum = ns_final;
               ne_sum = ne_final;
               jc_sum = jc_final;
               js_sum = js_final;
               sa_sum = sa_final;
            #endif

            if(vmpi::my_rank == 0) {
               std::ofstream ofile;
               ofile.open(std::string(filename.str()).c_str());
               // Header: positions in Angstroms, densities in SI, currents in A/m^2, spin accumulation in C/m^3
               ofile << "# stage " << stage << " index " << index
                     << " step " << sc1d_step_counter << " sim_time " << sim::time << std::endl;
               ofile << "# x(A)\ty(A)\tz(A)\tns(m^-3)\tne(m^-3)\tJc(A/m^2)\tJs_x(A/m^2)\tJs_y(A/m^2)\tJs_z(A/m^2)\tSa_x(C/m^3)\tSa_y(C/m^3)\tSa_z(C/m^3)\tSa_mag(C/m^3)\tm_x\tm_y\tm_z" << std::endl;
               
               // Get cell size for converting indices to physical positions
               const double cell_size_xy = micro_cell_size[0]; // Angstroms
               const double cell_size_z  = micro_cell_thickness; // Angstroms
               
               for(int cell=0; cell<size; ++cell){
                  if(cell_natom[cell] == 0) continue;
                  
                  // Convert cell indices to physical position (center of cell, in Angstroms)
                  const double x_A = (pos[3*cell+0] + 0.5) * cell_size_xy;
                  const double y_A = (pos[3*cell+1] + 0.5) * cell_size_xy;
                  const double z_A = (pos[3*cell+2] + 0.5) * cell_size_z;
                  const double mag = sqrt(m[3*cell+0]*m[3*cell+0] + m[3*cell+1]*m[3*cell+1] + m[3*cell+2]*m[3*cell+2]);
                  const double mag_norm = (mag <= 1.0e-4) ? 0.0: 1/mag;
                  const double sa_mag = sqrt(sa_sum[3*cell+0]*sa_sum[3*cell+0] + sa_sum[3*cell+1]*sa_sum[3*cell+1] + sa_sum[3*cell+2]*sa_sum[3*cell+2]);
                  ofile << std::fixed << std::setprecision(2);
                  ofile << x_A << "\t" << y_A << "\t" << z_A << "\t";
                  ofile << std::scientific << std::setprecision(6);
                  ofile << ns_sum[cell] << "\t" << ne_sum[cell] << "\t" << jc_sum[cell] << "\t";
                  ofile << js_sum[3*cell+0] << "\t" << js_sum[3*cell+1] << "\t" << js_sum[3*cell+2] << "\t";
                  ofile << sa_sum[3*cell+0] << "\t" << sa_sum[3*cell+1] << "\t" << sa_sum[3*cell+2] << "\t" << sa_mag << "\t";
                  
                  ofile << m[3*cell+0]*mag_norm << "\t" << m[3*cell+1]*mag_norm << "\t" << m[3*cell+2]*mag_norm << "\n";
               }
               ofile.close();
            }
            config_file_counter++;
         }

      }

      //-----------------------------------------------------------------------------
      // One-line cumulative timing for the 1D solver. Same output cadence as
      // output_sc1d_data. Rank 0 only; other ranks keep their own accumulators.
      //-----------------------------------------------------------------------------
#if ST_SC1D_TIMINGS
      static void append_sc1d_bucket(std::ostringstream& line,
                                     const char* name,
                                     const double seconds,
                                     const double inv_pct){
         std::ostringstream sec;
         sec << std::scientific << std::setprecision(3) << seconds;
         std::ostringstream pcts;
         pcts << std::fixed << std::setprecision(1) << (seconds * inv_pct);
         line << " " << name << "=" << sec.str() << " (" << pcts.str() << "%)";
      }
#endif

      void report_sc1d_timing(){
#if !ST_SC1D_TIMINGS
         return;
#else
         if(!sc1d_enable) return;

         const int rate = std::max(ST_output_rate, 1);
         unsigned long index = 0;
         if(sc1d_step_counter <= sc1d_relax_steps){
            index = sc1d_step_counter;
         } else if(sim::time <= sim::equilibration_time){
            index = static_cast<unsigned long>(sim::time);
         } else {
            index = static_cast<unsigned long>(sim::time - sim::equilibration_time);
         }
         const bool due = (index > 0UL) && (index % static_cast<unsigned long>(rate) == 0UL);
         if(!due) return;

         double total = sc1d_time_total;
         double io = sc1d_time_io;
         double charge = sc1d_time_charge;
         double spin_current = sc1d_time_spin_current;
         double spin_acc = sc1d_time_spin_acc;
         double interp = sc1d_time_interp;
         double bcast = sc1d_time_bcast;
         double other = sc1d_time_other;
#ifdef MPICF
         // Print the rank that did the most column work. Idle ranks still
         // sit in the broadcast, so a plain rank-0 sample can miss the solve.
         double local_times[8] = {total, io, charge, spin_current, spin_acc, interp, bcast, other};
         const int nrank = std::max(1, vmpi::num_processors);
         std::vector<double> all_times((std::size_t)nrank * 8u, 0.0);
         MPI_Allgather(local_times, 8, MPI_DOUBLE, all_times.data(), 8, MPI_DOUBLE, MPI_COMM_WORLD);
         if(vmpi::my_rank != 0) return;
         int busiest = 0;
         double busiest_work = -1.0;
         for(int r=0; r<nrank; ++r){
            const double* row = all_times.data() + (std::size_t)r * 8u;
            const double work = row[2] + row[3] + row[4] + row[5];
            if(work > busiest_work){
               busiest_work = work;
               busiest = r;
            }
         }
         const double* best = all_times.data() + (std::size_t)busiest * 8u;
         total = best[0];
         io = best[1];
         charge = best[2];
         spin_current = best[3];
         spin_acc = best[4];
         interp = best[5];
         bcast = best[6];
         other = best[7];
#endif
         const double inv_pct = (total > 0.0) ? (100.0 / total) : 0.0;

         std::ostringstream line;
         std::ostringstream tot;
         tot << std::scientific << std::setprecision(3) << total;
         line << "sc1d timing [s] total=" << tot.str();
         append_sc1d_bucket(line, "io", io, inv_pct);
         append_sc1d_bucket(line, "charge", charge, inv_pct);
         append_sc1d_bucket(line, "spin_current", spin_current, inv_pct);
         append_sc1d_bucket(line, "spin_acc", spin_acc, inv_pct);
         append_sc1d_bucket(line, "interp_coarse", interp, inv_pct);
         append_sc1d_bucket(line, "bcast", bcast, inv_pct);
         append_sc1d_bucket(line, "other", other, inv_pct);

         std::cout << line.str() << std::endl;
         zlog << zTs() << line.str() << std::endl;
#endif
      }


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
                        st::internal::a[cell] << "\t" << st::internal::b[cell] << "\t" << st::internal::cell_natom[cell] << std::endl;
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
      void output_microcell_sa_data(){
         
         // If 1D spin currents solver is enabled, skip non-1D output
         if(sc1d_enable) return;

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

      void output_microcell_data(){

         // If 1D spin currents solver is enabled, only output 1D data and return
         if(sc1d_enable){
            output_sc1d_data();
            return;
         }

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
               // MPI_Reduce(&j_final_up_y[0], &j_final_up_y_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               // MPI_Reduce(&j_final_down_y[0], &j_final_down_y_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
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
                  mag = 1;//(mag == 0.0) ? 0.0: 1/mag;
                  ofile << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] << "\t";
                  ofile << m[3*cell+0]*mag << "\t" << m[3*cell+1]*mag << "\t" << m[3*cell+2]*mag << "\t";
                 // if(st::internal::sot_check) ofile << (sa_sum[3*cell+0]-sa_infinity[cell]*m[3*cell]*mag)/sa_infinity[cell] << "\t" << (sa_sum[3*cell+1]-sa_infinity[cell]*m[3*cell+1]*mag)/sa_infinity[cell] << "\t" << (sa_sum[3*cell+2]-m[3*cell+2]*mag)/sa_infinity[cell] << "\t";
                  //else 
                  ofile << std::fixed << std::setprecision(10);
                  ofile << sa_sum[3*cell+0] << "\t" << sa_sum[3*cell+1] << "\t" << sa_sum[3*cell+2] << "\t";
                  ofile << std::setprecision(6);
                  ofile << j_final_up_x_sum[3*cell+0] << "\t" << j_final_up_x_sum[3*cell+1] << "\t" << j_final_up_x_sum[3*cell+2] << "\t";
                  // ofile << j_final_up_y_sum[3*cell+0] << "\t" << j_final_up_y_sum[3*cell+1] << "\t" << j_final_up_y_sum[3*cell+2] << "\t";
                  // ofile << j_final_down_y_sum[3*cell+0] << "\t" << j_final_down_y_sum[3*cell+1] << "\t" << j_final_down_y_sum[3*cell+2] << "\t";
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
                  // j_final_up_y_sum[3*cell+0] = j_final_up_y_sum[3*cell+1] = j_final_up_y_sum[3*cell+2] = 0.0;
                  // j_final_down_y_sum[3*cell+0] = j_final_down_y_sum[3*cell+1] = j_final_down_y_sum[3*cell+2] = 0.0;
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
               // if(st::internal::sot_check) ofile << (sa_final[3*cell+0]-sa_infinity[cell]*m[3*cell]*mag)/sa_infinity[cell] << "\t" << (sa_final[3*cell+1]-sa_infinity[cell]*m[3*cell+1]*mag)/sa_infinity[cell] << "\t" << (sa_final[3*cell+2]-m[3*cell+2]*mag)/sa_infinity[cell] << "\t";
               // else 
               ofile << sa_final[3*cell+0] << "\t" << sa_final[3*cell+1] << "\t" << sa_final[3*cell+2] << "\t";
               ofile << j_final_up_x[3*cell+0] << "\t" << j_final_up_x[3*cell+1] << "\t" << j_final_up_x[3*cell+2] << "\t";
               // ofile << j_final_up_y[3*cell+0] << "\t" << j_final_up_y[3*cell+1] << "\t" << j_final_up_y[3*cell+2] << "\t";
               // ofile << j_final_down_y[3*cell+0] << "\t" << j_final_down_y[3*cell+1] << "\t" << j_final_down_y[3*cell+2] << "\t";
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

   //-----------------------------------------------------------------------------
   // Thermal gradient output functions
   //-----------------------------------------------------------------------------
   
   // Local variables for thermal gradient output
   std::ofstream thermal_microcell_file;
   std::ofstream thermal_temperature_file;
   int thermal_output_counter = 0;
   bool thermal_output_initialised = false;

   //-----------------------------------------------------------------------------
   // Function to output thermal gradient microcell configuration
   //-----------------------------------------------------------------------------
   void output_thermal_microcell_data(){
      if(!sc1d_thermal_gradients_enable || !ltmp::is_enabled()) return;
      
      // only output on root process
      #ifdef MPICF
      if(vmpi::my_rank==0){
      #else
      {
      #endif
         std::ofstream ofile;
         ofile.open("spin-acc/thermal_microcell_config.cfg");
         
         // Header: cell position (Angstroms), attenuation, material constants
         ofile << "# z(A)\tattenuation\tC_e(J/m^3/K)\tC_p(J/m^3/K)\tG_ep(J/s/m^3/K)\tk_e(J/s/m/K)\tk_p(J/s/m/K)\tT_Debye(K)\n";
         
         const int ncz = num_microcells_per_stack;
         const double dz = micro_cell_thickness; // Angstroms
         const double d = std::max(sc1d_optical_absorption_length * 1e10, 1e-8); // m to Angstroms
         const double z_max = (double)ncz * dz; // total stack height in Angstroms
         
         for(int cell = 0; cell < ncz; ++cell){
            const double z = (cell + 0.5) * dz; // cell center position
            const double distance_from_top = z_max - z; // distance from top surface
            const double attenuation = std::exp(-distance_from_top / d);
            
            ofile << std::fixed << std::setprecision(2);
            ofile << z << "\t";
            ofile << std::scientific << std::setprecision(6);
            ofile << attenuation << "\t";
            if(cell < (int)sc1d_Ce.size()){
               ofile << sc1d_Ce[cell] << "\t";
               ofile << sc1d_Cp[cell] << "\t";
               ofile << sc1d_G[cell] << "\t";
               ofile << sc1d_kappa_e[cell] << "\t";
               ofile << sc1d_kappa_p[cell] << "\t";
               ofile << sc1d_T_Debye[cell] << "\n";
            } else {
               ofile << "0\t0\t0\t0\t0\t0\n";
            }
         }
         
         ofile.close();
      }
      
      return;
   }

   //-----------------------------------------------------------------------------
   // Function to open temperature profile output file
   //-----------------------------------------------------------------------------
   void open_thermal_temperature_profile_file(){
      if(!sc1d_thermal_gradients_enable || !ltmp::is_enabled()) return;
      
      thermal_output_counter = 0;
      thermal_output_initialised = true;
      
      #ifdef MPICF
      if(vmpi::my_rank==0){
      #else
      {
      #endif
         thermal_temperature_file.open("spin-acc/thermal_temperature_profile.dat");
         // Header: time step, then Te and Tp for each cell (stack 0 only)
         thermal_temperature_file << "# step\t";
         const int ncz = num_microcells_per_stack;
         for(int cell = 0; cell < ncz; ++cell){
            thermal_temperature_file << "Te_" << cell << "\t" << "Tp_" << cell << "\t";
         }
         thermal_temperature_file << "\n";
      }
      
      return;
   }

   //-----------------------------------------------------------------------------
   // Function to write temperature profile data
   //-----------------------------------------------------------------------------
   void write_thermal_temperature_data(){
      if(!sc1d_thermal_gradients_enable || !ltmp::is_enabled()) return;
      if(!thermal_output_initialised) return;
      
      if(sim::time % ST_output_rate != 0) return;
      
      #ifdef MPICF
      if(vmpi::my_rank==0){
      #else
      {
      #endif
         const int ncz = num_microcells_per_stack;
         const double dz = micro_cell_thickness;
         
         thermal_temperature_file << thermal_output_counter << "\t";
         
         for(int cell = 0; cell < ncz; ++cell){
            const double zA = (cell + 0.5) * dz;
            const double Te = ltmp::get_electron_temperature_at_z(zA);
            const double Tp = ltmp::get_phonon_temperature_at_z(zA);
            thermal_temperature_file << std::fixed << std::setprecision(1);
            thermal_temperature_file << Te << "\t" << Tp << "\t";
         }
         
         thermal_temperature_file << std::endl;
      }
      
      thermal_output_counter++;
      
      return;
   }

   } // end of internal namespace
} // end of st namespace
