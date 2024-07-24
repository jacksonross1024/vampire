//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <fstream>
#include <iostream>

// Vampire headers
#include "stopwatch.hpp"
#include "spintorque.hpp"
#include "vmpi.hpp"
#include "material.hpp"
#include "create.hpp"
#include <complex>
#include "program.hpp"
#include "sim.hpp"
// Spin Torque headers
#include "internal.hpp"

namespace st{
   namespace internal{

      //-----------------------------------------------------------------------------
      // Funtion to calculate the spin accumulation and spin torque
      //-----------------------------------------------------------------------------
   //  void  calculate_spin_accumulation(){}
      void calculate_spin_accumulation(){
         if(sot_sa) std::cout << "standard spin acc and sot spin acc activated " << std::endl;
         stopwatch_t stopwatch;
         stopwatch.start();

         // Zero all shared arrays (essential for parallelisation to work)
         std::fill (spin_torque.begin(),spin_torque.end(),0.0);
         std::fill (sa_final.begin(),sa_final.end(),0.0);
         // std::fill (sa_sot.begin(),sa_sot.end(),0.0);
         // std::fill (j.begin(),j.end(),0.0);
         // Declare resuable temporary variables
         matrix_t itm; // inverse transformation matrix
         matrix_t M; // general matrix
         three_vector_t V(0.0,0.0,0.0); // general 3-vector

         // reference basis vectors
         const three_vector_t bx(1.0,0.0,0.0);
         const three_vector_t by(0.0,1.0,0.0);
         const three_vector_t bz(0.0,0.0,1.0);

         // local basis vectors
         three_vector_t b1(1.0,0.0,0.0);
         three_vector_t b2(0.0,1.0,0.0);
         three_vector_t b3(0.0,0.0,1.0);

         // set local constants
         double je_eff = program::fractional_electric_field_strength* je_eff; // current (C/s) 

         //---------------------------------------------------------------------------------------------------
         //set parameters for TMR calculation
         if(TMRenable == true){
            std::cout << "wrong" << std::endl;
            int FL =	free_layer;
            int RL =	reference_layer;
            double dot = magx_mat[RL]*magx_mat[FL]+
                         magy_mat[RL]*magy_mat[FL]+
                         magz_mat[RL]*magz_mat[FL];

            // RE - this code is not general! Needs to be fixed. Placeholder added in the meantime
            // AM (2020) - Code fixed, but still codes refers to MTJ RL/barrier/FL specifically
            double MgO_thickness = (create::get_material_height_min(FL)-create::get_material_height_max(RL))*cs::system_dimensions[2]*1.0e-10;

            //calculate the relative angle of two FMs
            rel_angle = acos(dot);
            double plus_cos = 1.0+cos(rel_angle);
            double minus_cos = 1.0-cos(rel_angle);
            double exp_t = exp(-MgO_thickness/0.25e-9);

            double jtunnel = je_eff*0.5*(plus_cos+0.5*minus_cos)*exp_t;
//            std::cout << "t_MgO=( " << create::get_material_height_min(FL) << " - " << create::get_material_height_max(RL) << " ) = " << MgO_thickness << "\tje_eff\t" << je_eff << "\tje_eff_tun\t" << jtunnel << std::endl;

            //set the current je_eff and spin poralisation parameters
            je_eff = jtunnel;
            // AM (2020) - I think the default parameters should be rescaled by same factor as tunnelling current and not changed using those of material 0 arbitrarily
            // default_properties.beta_cond *= mp[0].beta_cond*0.5*(plus_cos+0.5*minus_cos)*exp_t;
            // default_properties.beta_diff *= mp[0].beta_diff*0.5*(plus_cos+0.5*minus_cos)*exp_t;

            // Calculate spin torque parameters
            for(size_t cell=0; cell<beta_cond.size(); ++cell){

               // check for zero atoms in cell
               if(cell_natom[cell] <= 0.0001){
                  beta_cond[cell]   = default_properties.beta_cond;
                  beta_diff[cell]   = default_properties.beta_diff;

                  const double hbar = 1.05457162e-34;
                  const double B  = beta_cond[cell];
                  const double Bp = beta_diff[cell];
                  // const double lambda_sdl = lambda_sdl[cell];
                  const double Do = diffusion[cell];
                  const double Jsd = sd_exchange[cell];

                  const double BBp = 1.0/sqrt(1.0-B*Bp);
                  const double lambda_sf = lambda_sdl[cell]*BBp;
                  const double lambda_j = sqrt(2.0*hbar*Do/Jsd); // Angstroms
                  const double lambda_sf2 = lambda_sf*lambda_sf;
                  const double lambda_j2 = lambda_j*lambda_j;

                  std::complex<double> inside (1.0/lambda_sf2, -1.0/lambda_j2);
                  std::complex<double> inv_lplus = sqrt(inside);

                  a[cell] =  real(inv_lplus);
                  b[cell] = -imag(inv_lplus);
               }
            }
            output_base_microcell_data();
         }

         //---------------------------------------------------------------------------------------------------

         const double i_muB = 1.0/9.274e-24; // J/T
         const double i_e = 1.0/1.60217662e-19; // electronic charge (Coulombs)
         const double microcell_volume = (micro_cell_size[stx] *
                                          micro_cell_size[sty] *
                                          micro_cell_thickness)*1.e-30; // m^3
         const double atomcell_volume = 15.7624e-30;
         // loop over all 1D stacks (in parallel)
         int int_stacks;
         #ifdef MPICH
            int_stacks = mpi_stack_list_y.size();
         #else
            int_stacks = stack_index_y.size();
         #endif
         // std::vector<int> stacks_list;
         // stacks_list = mpi_stack_list           
        

         int stack = 0;
         // std::cout << int_stacks << ", " << mpi_stack_list.size() << std::endl;
         for(int s=0; s < int_stacks; ++s) {
            #ifdef MPICH
               stack = mpi_stack_list_y.at(s);
            #else 
               stack = stack_index_x.at(s);
            #endif
            // if(stack % 3 == 0) continue;
               // std::cout << stack << std::endl;
               // std::cout << spin_acc_sign[stack_index[stack]] << ", " << stack << std::endl;
      
            // std::cout << vmpi::my_rank << ", " << stack << std::endl;
            // determine starting cell in stack
            const int idx = stack_index_y[stack] + 1;
           // std::cout << stack << ", " << idx << std::endl;
            // set initial values
           if(fbc) { 
               sa_final[3*idx+0] = 0.0;
               sa_final[3*idx+1] = 0.0;
               sa_final[3*idx+2] = 0.0;
               j_final_up_x [3*idx+0] = initial_beta*je_eff*initial_m[0];
               j_final_up_x [3*idx+1] = initial_beta*je_eff*initial_m[1];
               j_final_up_x [3*idx+2] = initial_beta*je_eff*initial_m[2];
           }  else if (!fbc) {
               const double mod = 1.0/sqrt(m [3*idx+0]*m [3*idx+0] + m [3*idx+1]*m [3*idx+1] + m [3*idx+2]*m[3*idx+2]);
               sa_final[3*idx+0] = sa_infinity[idx]*m [3*idx+0]*mod;
               sa_final[3*idx+1] = sa_infinity[idx]*m [3*idx+1]*mod;
               sa_final[3*idx+2] = sa_infinity[idx]*m [3*idx+2]*mod;

               j_final_up_x [3*idx+0] = initial_beta*je_eff*m [3*idx+0]*mod;
               j_final_up_x [3*idx+1] = initial_beta*je_eff*m [3*idx+1]*mod;
               j_final_up_x [3*idx+2] = initial_beta*je_eff*m [3*idx+2]*mod;
           }

          //  if(sim::time%(ST_output_rate) ==0) std::cout<< stack << "\t" << default_properties.sa_infinity << "\t" << init_stack_mag[((stack)%6)*3 + 0]/sqrt(2.0) << "\t" << init_stack_mag[((stack)%6)*3 + 1]/sqrt(2.0) << "\t" << m[3*idx+0] << " \t" <<  m[3*idx+1] << "\t" << sa[3*idx+0] << "\t" << sa[3*idx+1] << std::endl;

            // loop over all cells in stack after first (idx+1)
            for(int cell=idx+1; cell<idx+num_microcells_per_stack; ++cell) {

               // calculate cell id's
               const int cellx = 3*cell+0;
               const int celly = 3*cell+1;
               const int cellz = 3*cell+2;

               // calculate previous cell id's
               const int pcellx = 3*(cell-1)+0;
               const int pcelly = 3*(cell-1)+1;
               const int pcellz = 3*(cell-1)+2;

               three_vector_t  m_local(0.0,  0.0,  0.0);
               // copy array values to temporaries for readability
       
                  m_local.x = m[cellx];
                  m_local.y = m[celly];
                  m_local.z = m[cellz]; // current cell magnetisations
               
               three_vector_t pm_local(0.0,0.0,0.0);
         
                  pm_local.x = m[pcellx];
                  pm_local.y = m[pcelly];
                  pm_local.z = m[pcellz]; // current cell magnetisations
               

               const double modm = sqrt(m_local.x*m_local.x + m_local.y*m_local.y + m_local.z*m_local.z);
               const double pmodm = sqrt(pm_local.x*pm_local.x + pm_local.y*pm_local.y + pm_local.z*pm_local.z);

               // Check for zero magnetization in normalization
               if(modm > 1.e-11){
                  m_local.x = m_local.x/modm;
                  m_local.y = m_local.y/modm;
                  m_local.z = m_local.z/modm;
               }
               else{
                  m_local.x = 0.0;
                  m_local.y = 0.0;
                  m_local.z = 0.0;
               }
               if(pmodm > 1.e-11){
                  pm_local.x = pm_local.x/pmodm;
                  pm_local.y = pm_local.y/pmodm;
                  pm_local.z = pm_local.z/pmodm;
               }
               else{
                  pm_local.x = 0.0;
                  pm_local.y = 0.0;
                  pm_local.z = 0.0;
               }

               //---------------------------------------------------------------------
               // Step 1 calculate coordinate transformation m -> m' for current cell
               //---------------------------------------------------------------------

               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m_local, itm);

                  // Calculate basis vectors
                  b1 = transform_vector(bz, itm);
                  b2 = transform_vector(bx, itm);
                  b3 = transform_vector(by, itm);
               }
               else{
                  // set default rotational frame
                  b1 = bz;
                  b2 = bx;
                  b3 = by;
               }
               //---------------------------------------------------------------------
               // Step 2 determine coefficients for the spin accumulation
               //---------------------------------------------------------------------

               // Initialise temporary constants
               const double Bc = beta_cond[cell]; // beta
               const double Bd = beta_diff[cell]; // beta_prime
               const double Do = diffusion[cell];
               three_vector_t jm0(j_final_up_x[pcellx],j_final_up_x[pcelly],j_final_up_x[pcellz]);

               //  Calculate gradient dsacc/dx
               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m_local.x*m_local.x - twoDo;     M.xy = preD*m_local.y*m_local.x;             M.xz = preD*m_local.z*m_local.x;
               M.yx = preD*m_local.x*m_local.y;             M.yy = preD*m_local.y*m_local.y - twoDo;     M.yz = preD*m_local.z*m_local.y;
               M.zx = preD*m_local.x*m_local.z;             M.zy = preD*m_local.y*m_local.z;             M.zz = preD*m_local.z*m_local.z - twoDo;

                  V.x = jm0.x - Bc*je_eff*m_local.x;
                  V.y = jm0.y - Bc*je_eff*m_local.y;
                  V.z = jm0.z - Bc*je_eff*m_local.z;
               
               const three_vector_t divm_0 = gaussian_elimination(M, V);

               // Calculate mp(0), c and d
               const double i_lsdl = 1.0/lambda_sdl[cell];
               double mp_inf = sa_infinity[cell];
                
               const double a_local = a[cell];
               const double b_local = b[cell];
               const double two_a = 2.0*a_local;
               const double two_b = 2.0*b_local;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               const three_vector_t coeff = gaussian_elimination(M, V);

               const double mp_0 = coeff.x;
               const double c    = coeff.y;
               const double d    = coeff.z;

               //------------------------------------
               // Step 3 calculate spin accumulation
               //------------------------------------

               const double x = micro_cell_thickness*1.0e-10; // Convert to metres
               const double cos_bx = cos(b_local*x);
               const double sin_bx = sin(b_local*x);
               const double e_xsdl = exp(-x*i_lsdl);
               const double e_ax   = exp(-a_local*x);
               const double prefac = (2.0*e_ax);

               const double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               const double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               const double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);

               //convert mp and m_perp
                double sax = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
                double say = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
                double saz = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;
              
               //--------------------------------------------
               // Step 4 calculate the spin current (jm_end)
               //--------------------------------------------
               const double ac_bd = a_local*c + b_local*d; 
               const double ad_bc = a_local*d - b_local*c;

               const double divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl;
               const double divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx);
               const double divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx);

               const double divsax = b1.x*divsa_para + b2.x*divsa_perp2 + b3.x*divsa_perp3;
               const double divsay = b1.y*divsa_para + b2.y*divsa_perp2 + b3.y*divsa_perp3;
               const double divsaz = b1.z*divsa_para + b2.z*divsa_perp2 + b3.z*divsa_perp3;

               const double dot = m_local.x*divsax + m_local.y*divsay + m_local.z*divsaz;

               const double pre_jmx = divsax - Bc*Bd*m_local.x*dot;
               const double pre_jmy = divsay - Bc*Bd*m_local.y*dot;
               const double pre_jmz = divsaz - Bc*Bd*m_local.z*dot;

                double jmx = Bc*je_eff*m_local.x - twoDo*pre_jmx;
                double jmy = Bc*je_eff*m_local.y - twoDo*pre_jmy;
                double jmz = Bc*je_eff*m_local.z - twoDo*pre_jmz;

          
               if(cell_natom[cell]>0) {
      
                  sa_final[cellx] = sax;
                  sa_final[celly] = say;
                  sa_final[cellz] = saz;

                  // Save values for the spin current
                  j_final_up_x[cellx] = jmx;
                  j_final_up_x[celly] = jmy;
                  j_final_up_x[cellz] = jmz;

                  // }  else {// Calculate spin torque energy for cell (Joules)
                  spin_torque[cellx] += atomcell_volume * sd_exchange[cell] * (sax) * i_e * i_muB;
                  spin_torque[celly] += atomcell_volume * sd_exchange[cell] * (say) * i_e * i_muB;
                  spin_torque[cellz] += atomcell_volume * sd_exchange[cell] * (saz) * i_e * i_muB; 
                  // }
               } 
               else{
                  // Save values for the spin accumulation
                  sa_final[cellx] = sa_final[pcellx];
                  sa_final[celly] = sa_final[pcelly];
                  sa_final[cellz] = sa_final[pcellz];

                  // Save values for the spin current
                  j_final_up_x[cellx] = j_final_up_x[pcellx];
                  j_final_up_x[celly] = j_final_up_x[pcelly];
                  j_final_up_x[cellz] = j_final_up_x[pcellz];

                  // Calculate spin torque energy for cell (Joules)
                  spin_torque[cellx] = spin_torque[pcellx];
                  spin_torque[celly] = spin_torque[pcelly];
                  spin_torque[cellz] = spin_torque[pcellz];
                  // }
               }
               if(output_torque_data == "final") {
               //--------------------------------------------
               // Step 5 calculate the spin torque of each cell
               //--------------------------------------------

               //convert M of previous cell into basis b1, b2, b3

               M.xx = b1.x;    M.xy = b2.x;    M.xz = b3.x;
               M.yx = b1.y;    M.yy = b2.y;    M.yz = b3.y;
               M.zx = b1.z;    M.zy = b2.z;    M.zz = b3.z;

               V.x = pm_local.x;
               V.y = pm_local.y;
               V.z = pm_local.z;

               const three_vector_t pm_basis = gaussian_elimination(M, V);

               //const double pm_b1 = pm_basis.x; // unused variable
               const double pm_b2 = pm_basis.y;
               const double pm_b3 = pm_basis.z;

               // Calculate the spin torque coefficients describing ast and nast
               const double prefac_sc = microcell_volume * sd_exchange[cell] * i_e * i_muB;
               const double plus_perp =  (pm_b2*pm_b2 + pm_b3*pm_b3);
               // const double minus_perp = (pm_b2*pm_b2 - pm_b3*pm_b3); // unused variable

	            double aj; // the ST parameter describing Slonczewski torque
               double bj; // the ST parameter describing field-like torque


                if( ( plus_perp <= 1.0e-11 ) ){
                    aj = 0.0;
                    bj = 0.0; }

                else{
                    aj  = prefac_sc*(sa_perp2*pm_b3 - sa_perp3*pm_b2)/plus_perp;
                    bj  = prefac_sc*(sa_perp2*pm_b2 + sa_perp3*pm_b3)/plus_perp;
                }

                double SxSp[3], SxSxSp[3];
                coeff_ast[cell]  = aj;
                coeff_nast[cell] = bj;

                SxSp[0]=(m_local.y*pm_local.z-m_local.z*pm_local.y);
                SxSp[1]=(m_local.z*pm_local.x-m_local.x*pm_local.z);
                SxSp[2]=(m_local.x*pm_local.y-m_local.y*pm_local.x);

                SxSxSp[0]= (m_local.y*SxSp[2]-m_local.z*SxSp[1]);
                SxSxSp[1]= (m_local.z*SxSp[0]-m_local.x*SxSp[2]);
                SxSxSp[2]= (m_local.x*SxSp[1]-m_local.y*SxSp[0]);

                //calculate directly from J(Sxm)
                total_ST[cellx] = prefac_sc*(m_local.y*saz-m_local.z*say);
                total_ST[celly] = prefac_sc*(m_local.z*sax-m_local.x*saz);
                total_ST[cellz] = prefac_sc*(m_local.x*say-m_local.y*sax);

                ast[cellx] = -aj*SxSxSp[0];
                ast[celly] = -aj*SxSxSp[1];
                ast[cellz] = -aj*SxSxSp[2];

                nast[cellx] = bj*SxSp[0];
                nast[celly] = bj*SxSp[1];
                nast[cellz] = bj*SxSp[2];
               }
            } // end of cell loop


            //reverse process
         } // end of stack loop

         // Reduce all microcell spin torques on all nodes
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &spin_torque[0],spin_torque.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            // if(output_torque_data) MPI_Allreduce(MPI_IN_PLACE, &total_ST[0],total_ST.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif
         output_microcell_data();

         st::spin_acc_time += stopwatch.elapsed_seconds();

         return;
      }

      void calculate_sot_accumulation(){

         stopwatch_t stopwatch;
         stopwatch.start();

         // Zero all shared arrays (essential for parallelisation to work)
         std::fill (spin_torque.begin(),spin_torque.end(),0.0);
         std::fill (sa_final.begin(),sa_final.end(),0.0);
         std::fill (sa_int.begin(),sa_int.end(),0.0);
         
         std::fill (j_init_up_y.begin(),j_init_up_y.end(),0.0);
         std::fill (j_init_down_y.begin(),j_init_down_y.end(),0.0);
         std::fill (j_int_up_y.begin(),j_int_up_y.end(),0.0);
         std::fill (j_int_down_y.begin(),j_int_down_y.end(),0.0);
         
         matrix_t itm; // inverse transformation matrix
         matrix_t M; // general matrix
         three_vector_t V(0.0,0.0,0.0); // general 3-vector

         // reference basis vectors
         const three_vector_t bx(1.0,0.0,0.0);
         const three_vector_t by(0.0,1.0,0.0);
         const three_vector_t bz(0.0,0.0,1.0);

         // local basis vectors
         three_vector_t b1(1.0,0.0,0.0);
         three_vector_t b2(0.0,1.0,0.0);
         three_vector_t b3(0.0,0.0,1.0);

         // set local constants
         const double je_eff = program::fractional_electric_field_strength* je; // current (C/s)
        
         //---------------------------------------------------------------------------------------------------
         //no TMR calculation
         //---------------------------------------------------------------------------------------------------

         const double i_muB = 1.0/9.274e-24; // J/T
         const double i_e = 1.0/1.60217662e-19; // electronic charge (Coulombs)
         const double microcell_volume = (micro_cell_size[stx] *
                                          micro_cell_size[sty] *
                                          micro_cell_thickness)*1.e-30; // m^3
         const double atomcell_volume = 15.7624e-30; //hard code for Mn2Au atom volume for now

         // need mpi run for now
         int   int_stacks = mpi_stack_list_x.size();
         int stack = 0;

         //STT first
         for(int s=0; s < int_stacks; ++s) {

            stack = mpi_stack_list_x.at(s);          
            const int idx = stack_index_x[stack];
   
           if(fbc) { //fbc uses initial beta for spin polarisation from injection
               j_final_up_x[3*idx+0] = initial_beta*je_eff*initial_m[0];
               j_final_up_x[3*idx+1] = initial_beta*je_eff*initial_m[1];
               j_final_up_x[3*idx+2] = initial_beta*je_eff*initial_m[2];
           } else if (!fbc) { //assume infinite magnetisation precedding. 
                double mod = sqrt(m[3*idx+0]*m[3*idx+0] + m[3*idx+1]*m[3*idx+1] + m[3*idx+2]*m[3*idx+2]);
                if(mod > 1e-11) mod = 1.0/mod;
               sa_int[3*idx+0] = sa_infinity[idx]*m[3*idx+0]*mod;
               sa_int[3*idx+1] = sa_infinity[idx]*m[3*idx+1]*mod;
               sa_int[3*idx+2] = sa_infinity[idx]*m[3*idx+2]*mod;

               j_final_up_x[3*idx+0] = initial_beta*je_eff*m[3*idx+0]*mod; //need to remove beta_initial dependency
               j_final_up_x[3*idx+1] = initial_beta*je_eff*m[3*idx+1]*mod;
               j_final_up_x[3*idx+2] = initial_beta*je_eff*m[3*idx+2]*mod;

           }
         
            //standard stt method in x plane
            std::vector<int> cell_container_list = cell_index_x[stack];
            for(int xc = 1; xc < cell_container_list.size(); ++xc) {
              

               int cell = cell_index_x[stack][xc];
               
               // calculate cell id's
               const int cellx = 3*cell+0;
               const int celly = 3*cell+1;
               const int cellz = 3*cell+2;

               // calculate previous cell id's
               const int pcellx = 3*cell_index_x[stack][xc-1]+0;
               const int pcelly = 3*cell_index_x[stack][xc-1]+1;
               const int pcellz = 3*cell_index_x[stack][xc-1]+2;

               three_vector_t  m_local(0.0,  0.0,  0.0);
               // copy array values to temporaries for readability
             
                  m_local.x = m[cellx];
                  m_local.y = m[celly];
                  m_local.z = m[cellz]; // current cell new SA

               const double modm = sqrt(m_local.x*m_local.x + m_local.y*m_local.y + m_local.z*m_local.z);
               
               // Check for zero magnetization in normalization
               if(modm > 1.e-11){
                  m_local.x = m_local.x/modm;
                  m_local.y = m_local.y/modm;
                  m_local.z = m_local.z/modm;
               }
               else{
              //    std::cout << "problem " << std::endl;
                  m_local.x = 0.0;
                  m_local.y = 0.0;
                  m_local.z = 0.0;
               }

               //---------------------------------------------------------------------
               // Step 1 calculate coordinate transformation m -> m' for current cell
               //---------------------------------------------------------------------

               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m_local, itm);

                  // Calculate basis vectors
                  b1 = transform_vector(bz, itm);
                  b2 = transform_vector(bx, itm);
                  b3 = transform_vector(by, itm);
               }
               else{
                  // set default rotational frame
                  b1 = bz;
                  b2 = bx;
                  b3 = by;
               }
               //---------------------------------------------------------------------
               // Step 2 determine coefficients for the spin accumulation
               //---------------------------------------------------------------------

               // Initialise temporary constants
               const double Bc = beta_cond[cell]; // beta
               const double Bd = beta_diff[cell]; // beta_prime
               const double Do = diffusion[cell];
               three_vector_t jm0(j_final_up_x[pcellx],j_final_up_x[pcelly],j_final_up_x[pcellz]);

               //  Calculate gradient dsacc/dx
               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m_local.x*m_local.x - twoDo;     M.xy = preD*m_local.y*m_local.x;             M.xz = preD*m_local.z*m_local.x;
               M.yx = preD*m_local.x*m_local.y;             M.yy = preD*m_local.y*m_local.y - twoDo;     M.yz = preD*m_local.z*m_local.y;
               M.zx = preD*m_local.x*m_local.z;             M.zy = preD*m_local.y*m_local.z;             M.zz = preD*m_local.z*m_local.z - twoDo;

                  V.x = jm0.x - Bc*je_eff*m_local.x;
                  V.y = jm0.y - Bc*je_eff*m_local.y;
                  V.z = jm0.z - Bc*je_eff*m_local.z;
               
               const three_vector_t divm_0 = gaussian_elimination(M, V);

               // Calculate mp(0), c and d
               const double i_lsdl = 1.0/lambda_sdl[cell];
               double mp_inf = sa_infinity[cell];
              
               const double a_local = a[cell];
               const double b_local = b[cell];
               const double two_a = 2.0*a_local;
               const double two_b = 2.0*b_local;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               const three_vector_t coeff = gaussian_elimination(M, V);

               const double mp_0 = coeff.x;
               const double c    = coeff.y;
               const double d    = coeff.z;

               //------------------------------------
               // Step 3 calculate spin accumulation
               //------------------------------------

               const double x = micro_cell_size[0]*1.0e-10; // Convert to metres
               const double cos_bx = cos(b_local*x);
               const double sin_bx = sin(b_local*x);
               const double e_xsdl = exp(-x*i_lsdl);
               const double e_ax   = exp(-a_local*x);
               const double prefac = (2.0*e_ax);

               const double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               const double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               const double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);    

               //--------------------------------------------
               // Step 4 calculate the spin current (jm_end)
               //--------------------------------------------
               const double ac_bd = a_local*c + b_local*d; 
               const double ad_bc = a_local*d - b_local*c;

               const double divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl;
               const double divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx);
               const double divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx);

               const double divsax = b1.x*divsa_para + b2.x*divsa_perp2 + b3.x*divsa_perp3;
               const double divsay = b1.y*divsa_para + b2.y*divsa_perp2 + b3.y*divsa_perp3;
               const double divsaz = b1.z*divsa_para + b2.z*divsa_perp2 + b3.z*divsa_perp3;

               const double dot = m_local.x*divsax + m_local.y*divsay + m_local.z*divsaz;

               const double pre_jmx = divsax - Bc*Bd*m_local.x*dot;
               const double pre_jmy = divsay - Bc*Bd*m_local.y*dot;
               const double pre_jmz = divsaz - Bc*Bd*m_local.z*dot;

               double jmx = Bc*je_eff*m_local.x - twoDo*pre_jmx;
               double jmy = Bc*je_eff*m_local.y - twoDo*pre_jmy;
               double jmz = Bc*je_eff*m_local.z - twoDo*pre_jmz;

               //calculate directly from J(Sxm)
               double sax = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
               double say = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
               double saz = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;

               if(cell_natom[cell]>0) {

                  // Save values for the spin current
                  j_final_up_x[cellx] = jmx;
                  j_final_up_x[celly] = jmy;
                  j_final_up_x[cellz] = jmz;

                  //spin polarisation in cell becomes up and down drift current for z axis
                  j_init_up_y[cellx] = twoDo*pre_jmx;
                  j_init_up_y[celly] = twoDo*pre_jmy;
                  j_init_up_y[cellz] = twoDo*pre_jmz;

                  j_init_down_y[cellx] = twoDo*pre_jmx;
                  j_init_down_y[celly] = twoDo*pre_jmy;
                  j_init_down_y[cellz] = twoDo*pre_jmz;

                  //save intermediate sa for z axis
                  sa_int[cellx] = sax;
                  sa_int[celly] = say;
                  sa_int[cellz] = saz; 

               } 
               else{

                  // Save values for the spin current
                  j_final_up_x[cellx] = jmx;//j_final_up_x[pcellx];
                  j_final_up_x[celly] = jmy;//j_final_up_x[pcelly];
                  j_final_up_x[cellz] = jmz;//j_final_up_x[pcellz];

               }
               if(output_torque_data == "final") {
               //--------------------------------------------
               // Step 5 calculate the spin torque of each cell
               //--------------------------------------------
                three_vector_t  pm_local(0.0,  0.0,  0.0);
               // copy array values to temporaries for readability
             
                  pm_local.x = sa_int[pcellx];
                  pm_local.y = sa_int[pcelly];
                  pm_local.z = sa_int[pcellz]; // current cell magnetisations

               const double pmodm = sqrt(pm_local.x*pm_local.x + pm_local.y*pm_local.y + pm_local.z*pm_local.z);
               // const double pmodm = sqrt(pm_local.x*pm_local.x + pm_local.y*pm_local.y + pm_local.z*pm_local.z);

               // Check for zero magnetization in normalization
               if(pmodm > 1.e-11){
                  pm_local.x = pm_local.x/pmodm;
                  pm_local.y = pm_local.y/pmodm;
                  pm_local.z = pm_local.z/pmodm;
               }
               else{
                  pm_local.x = 0.0;
                  pm_local.y = 0.0;
                  pm_local.z = 0.0;
               }

               // Check for cell magnetization greater than 1e-8 mu_B
               if(pmodm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(pm_local, itm);

                  // Calculate basis vectors
                  b1 = transform_vector(bz, itm);
                  b2 = transform_vector(bx, itm);
                  b3 = transform_vector(by, itm);
               }
               else{
                  // set default rotational frame
                  b1 = bz;
                  b2 = bx;
                  b3 = by;
               }

               //convert M of previous cell into basis b1, b2, b3

               M.xx = b1.x;    M.xy = b2.x;    M.xz = b3.x;
               M.yx = b1.y;    M.yy = b2.y;    M.yz = b3.y;
               M.zx = b1.z;    M.zy = b2.z;    M.zz = b3.z;

               V.x = pm_local.x;
               V.y = pm_local.y;
               V.z = pm_local.z;

               const three_vector_t pm_basis = gaussian_elimination(M, V);

               //const double pm_b1 = pm_basis.x; // unused variable
               const double pm_b2 = pm_basis.y;
               const double pm_b3 = pm_basis.z;

               // Calculate the spin torque coefficients describing ast and nast
               const double prefac_sc = atomcell_volume * sd_exchange[cell] * i_e * i_muB;
               const double plus_perp =  (pm_b2*pm_b2 + pm_b3*pm_b3);
               // const double minus_perp = (pm_b2*pm_b2 - pm_b3*pm_b3); // unused variable

	            double aj; // the ST parameter describing Slonczewski torque
               double bj; // the ST parameter describing field-like torque


                if( ( plus_perp <= 1.0e-11 ) ){
                    aj = 0.0;
                    bj = 0.0; }

                else{
                    aj  = prefac_sc*(sa_perp2*pm_b3 - sa_perp3*pm_b2)/plus_perp;
                    bj  = prefac_sc*(sa_perp2*pm_b2 + sa_perp3*pm_b3)/plus_perp;
                }

               double SxSp[3], SxSxSp[3];
               coeff_ast[cell]  = aj;
               coeff_nast[cell] = bj;

               SxSp[0]=(m_local.y*pm_local.z-m_local.z*pm_local.y);
               SxSp[1]=(m_local.z*pm_local.x-m_local.x*pm_local.z);
               SxSp[2]=(m_local.x*pm_local.y-m_local.y*pm_local.x);

               SxSxSp[0]= (m_local.y*SxSp[2]-m_local.z*SxSp[1]);
               SxSxSp[1]= (m_local.z*SxSp[0]-m_local.x*SxSp[2]);
               SxSxSp[2]= (m_local.x*SxSp[1]-m_local.y*SxSp[0]);

               //calculate directly from J(Sxm)
                double Tx = prefac_sc*(m_local.y*saz-m_local.z*say);
                double Ty = prefac_sc*(m_local.z*sax-m_local.x*saz);
                double Tz = prefac_sc*(m_local.x*say-m_local.y*sax);
                total_ST[cellx] = Ty*m_local.z-Tz*m_local.y;
                total_ST[celly] = Tz*m_local.x-Tx*m_local.z;
                total_ST[cellz] = Tx*m_local.y-Ty*m_local.x;

               ast[cellx] += -aj*SxSxSp[0];
               ast[celly] += -aj*SxSxSp[1];
               ast[cellz] += -aj*SxSxSp[2];

               nast[cellx] += bj*SxSp[0];
               nast[celly] += bj*SxSp[1];
               nast[cellz] += bj*SxSp[2];
               }
            } // end of cell loop
         } // end of stack loop
 
         //reduce for new mpi decomp 
         #ifdef MPICF
            // MPI_Allreduce(MPI_IN_PLACE, &j_final_up_x[0],j_final_up_x.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &sa_int[0],sa_int.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &j_init_up_y[0],j_init_up_y.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &j_init_down_y[0],j_init_down_y.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif

         //new mpi decomp 
         int_stacks = mpi_stack_list_y.size();
         stack = 0;
      
         //sot now
         for(int s=0; s < int_stacks; ++s) {
            stack = mpi_stack_list_y.at(s);
            const int idx = stack_index_y[stack];

            //init process
            for(int cell=idx; cell<idx+num_microcells_per_stack; ++cell) {
               int direction_sign = sot_sa_source[cell] ? -1:1; //direction and sign if Au J_s
               
               // calculate cell id's
               const int cellx = 3*cell+0;
               const int celly = 3*cell+1;
               const int cellz = 3*cell+2;

               const int pcell = cell - 1;
               const int acell = cell + 1;

               const int pcellx = 3*pcell + 0;
               const int pcelly = 3*pcell + 1;
               const int pcellz = 3*pcell + 2;

               const int acellx = 3*acell + 0;
               const int acelly = 3*acell + 1;
               const int acellz = 3*acell + 2;

               three_vector_t  m_local(0.0,  0.0,  0.0);
               // copy array values to temporaries for readability
               if(sot_sa_source[cell]) {
                  m_local.x = 0;
                  m_local.y = 1.0; // Au cell artificial direction for correct basis of spin acc change
                  m_local.z = 0.0; 
               } else {
                  m_local.x = sa_int[cellx];
                  m_local.y = sa_int[celly];
                  m_local.z = sa_int[cellz]; // current cell magnetisations
               }
               const double modm = sqrt(m_local.x*m_local.x + m_local.y*m_local.y + m_local.z*m_local.z);
               // Check for zero magnetization in normalization
               if(modm > 1.e-11){
                  m_local.x = m_local.x/modm;
                  m_local.y = m_local.y/modm;
                  m_local.z = m_local.z/modm;
               }
               else{
                  m_local.x = 0.0;
                  m_local.y = 0.0;
                  m_local.z = 0.0;
               }

               //---------------------------------------------------------------------
               // Step 1 calculate coordinate transformation m -> m' for current cell
               //---------------------------------------------------------------------

               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m_local, itm);

                  // Calculate basis vectors
                  b1 = transform_vector(bz, itm);
                  b2 = transform_vector(bx, itm);
                  b3 = transform_vector(by, itm);
               }
               else{
                  // set default rotational frame
                  b1 = bz;
                  b2 = bx;
                  b3 = by;
               }
               //---------------------------------------------------------------------
               // Step 2 determine coefficients for the spin accumulation
               //---------------------------------------------------------------------

               // Initialise temporary constants
               const double Bc = sot_beta_cond[cell]; // beta
               const double Bd = sot_beta_diff[cell]; // beta_prime
               const double Do = sot_diffusion[cell];
               three_vector_t jm0(0.0,0.0,0.0); //different bounds for top and bottom layers
               if(cell == idx) {jm0 = {j_init_down_y[acellx],\
                                  j_init_down_y[acelly],\
                                  j_init_down_y[acellz]};
               } else if (cell == idx+num_microcells_per_stack - 1 ) {
                                jm0 = {j_init_up_y[pcellx],\
                                  j_init_up_y[pcelly],\
                                  j_init_up_y[pcellz]};
               } else {         jm0 = {j_init_up_y[pcellx]+j_init_down_y[acellx],\
                                  j_init_up_y[pcelly]+j_init_down_y[acelly],\
                                  j_init_up_y[pcellz]+j_init_down_y[acellz]};
               }

               //  Calculate gradient dsacc/dx
               // need to include Roy's new interface dsa contribution for multilayers
               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m_local.x*m_local.x - twoDo;     M.xy = preD*m_local.y*m_local.x;             M.xz = preD*m_local.z*m_local.x;
               M.yx = preD*m_local.x*m_local.y;             M.yy = preD*m_local.y*m_local.y - twoDo;     M.yz = preD*m_local.z*m_local.y;
               M.zx = preD*m_local.x*m_local.z;             M.zy = preD*m_local.y*m_local.z;             M.zz = preD*m_local.z*m_local.z - twoDo;

               //no charge current conversion from out of plane process
                  V.x = jm0.x;// - Bc*je_eff*m_local.x;
                  V.y = jm0.y;// - Bc*je_eff*m_local.x;
                  V.z = jm0.z;// - Bc*je_eff*m_local.z;
              
               const three_vector_t divm_0 = gaussian_elimination(M, V);

               // Calculate mp(0), c and d
               const double i_lsdl = 1.0/sot_lambda_sdl[cell];
               double mp_inf = modm;
               //modify spin acc on Au depending on current. current Au m_inf = 0.0
               if(sot_sa_source[cell]) mp_inf *= program::fractional_electric_field_strength;
               const double a_local = sot_a[cell];
               const double b_local = sot_b[cell];
               const double two_a = 2.0*a_local;
               const double two_b = 2.0*b_local;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               const three_vector_t coeff = gaussian_elimination(M, V);

               const double mp_0 = coeff.x;
               const double c    = coeff.y;
               const double d    = coeff.z;

               //------------------------------------
               // Step 3 calculate spin accumulation
               //------------------------------------

               const double x = micro_cell_thickness*1.0e-10; // Convert to metres
               const double cos_bx = cos(b_local*x);
               const double sin_bx = sin(b_local*x);
               const double e_xsdl = exp(-x*i_lsdl);
               const double e_ax   = exp(-a_local*x);
               const double prefac = (2.0*e_ax);

               const double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               const double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               const double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);
                           
               //--------------------------------------------
               // Step 4 calculate the spin current (jm_end)
               //--------------------------------------------
               const double ac_bd = a_local*c + b_local*d; 
               const double ad_bc = a_local*d - b_local*c;

               const double divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl;
               const double divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx);
               const double divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx);

               const double divsax = b1.x*divsa_para + b2.x*divsa_perp2 + b3.x*divsa_perp3;
               const double divsay = b1.y*divsa_para + b2.y*divsa_perp2 + b3.y*divsa_perp3;
               const double divsaz = b1.z*divsa_para + b2.z*divsa_perp2 + b3.z*divsa_perp3;

               const double dot = m_local.x*divsax + m_local.y*divsay + m_local.z*divsaz;

               const double pre_jmx = divsax - Bc*Bd*m_local.x*dot;
               const double pre_jmy = divsay - Bc*Bd*m_local.y*dot;
               const double pre_jmz = divsaz - Bc*Bd*m_local.z*dot;

                //no charge current conversion
               double jmx =   twoDo*pre_jmx;
               double jmy =   twoDo*pre_jmy;
                  //modify spin current by current density and spin hall theta 
                  if(sot_sa_source[cell]) jmy = -je_eff*initial_theta;
               double jmz =   twoDo*pre_jmz;
               

               //convert mp and m_perp
               double sax = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
               double say = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
               double saz = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;

               // // Save values. assume no 0 magnetisation cells
               sa_int[cellx] = sax;
               sa_int[celly] = say;
               sa_int[cellz] = saz;

               //init value saved in int for fast overwrite
               j_int_up_y[cellx] = jmx;
               j_int_up_y[celly] = jmy;
               j_int_up_y[cellz] = jmz;

               
               //down diffusion for Au spin current negative
               j_int_down_y[cellx] = jmx*direction_sign;
               j_int_down_y[celly] = jmy*direction_sign;
               j_int_down_y[cellz] = jmz*direction_sign;
               
               //conditional spin torque decompositon for output. faster to omit
            if(output_torque_data == "init") {
               //--------------------------------------------
               // Step 5 calculate the spin torque of each cell
               //--------------------------------------------

               //convert M of previous cell into basis b1, b2, b3
               const int pcellx = 3*(cell-1);
               const int pcelly = pcellx +1;
               const int pcellz = pcelly +1;

               three_vector_t  pm_local(0.0,  0.0,  0.0);
               // copy array values to temporaries for readability
               if(sot_sa_source[cell-1]) {
                  pm_local.x = 0;
                  pm_local.y = 1.0; // Au cell artificial direction for correct basis of spin acc change
                  pm_local.z = 0.0; 
               } else {
                  pm_local.x = m[pcellx];
                  pm_local.y = m[pcelly];
                  pm_local.z = m[pcellz]; // current cell magnetisations
               }
               const double pmodm = sqrt(pm_local.x*pm_local.x + pm_local.y*pm_local.y + pm_local.z*pm_local.z);
               // Check for zero magnetization in normalization
               if(pmodm > 1.e-11){
                  pm_local.x = pm_local.x/pmodm;
                  pm_local.y = pm_local.y/pmodm;
                  pm_local.z = pm_local.z/pmodm;
               }
               else{
                  pm_local.x = 0.0;
                  pm_local.y = 0.0;
                  pm_local.z = 0.0;
               }

               M.xx = b1.x;    M.xy = b2.x;    M.xz = b3.x;
               M.yx = b1.y;    M.yy = b2.y;    M.yz = b3.y;
               M.zx = b1.z;    M.zy = b2.z;    M.zz = b3.z;

               V.x = pm_local.x;
               V.y = pm_local.y;
               V.z = pm_local.z;

               const three_vector_t pm_basis = gaussian_elimination(M, V);

               //const double pm_b1 = pm_basis.x; // unused variable
               const double pm_b2 = pm_basis.y;
               const double pm_b3 = pm_basis.z;

               // Calculate the spin torque coefficients describing ast and nast
               const double prefac_sc = atomcell_volume * sd_exchange[cell] * i_e * i_muB;
               const double plus_perp =  (pm_b2*pm_b2 + pm_b3*pm_b3);
               // const double minus_perp = (pm_b2*pm_b2 - pm_b3*pm_b3); // unused variable

	            double aj; // the ST parameter describing Slonczewski torque
               double bj; // the ST parameter describing field-like torque


                if( ( plus_perp <= 1.0e-11 ) ){
                    aj = 0.0;
                    bj = 0.0; }

                else{
                    aj  = prefac_sc*(sa_perp2*pm_b3 - sa_perp3*pm_b2)/plus_perp;
                    bj  = prefac_sc*(sa_perp2*pm_b2 + sa_perp3*pm_b3)/plus_perp;
                }

               double SxSp[3], SxSxSp[3];
               coeff_ast[cell]  = aj;
               coeff_nast[cell] = bj;

               SxSp[0]=(m_local.y*pm_local.z-m_local.z*pm_local.y);
               SxSp[1]=(m_local.z*pm_local.x-m_local.x*pm_local.z);
               SxSp[2]=(m_local.x*pm_local.y-m_local.y*pm_local.x);

               SxSxSp[0]= (m_local.y*SxSp[2]-m_local.z*SxSp[1]);
               SxSxSp[1]= (m_local.z*SxSp[0]-m_local.x*SxSp[2]);
               SxSxSp[2]= (m_local.x*SxSp[1]-m_local.y*SxSp[0]);


               //calculate directly from J(Sxm)
               total_ST[cellx] = prefac_sc*(m_local.y*saz-m_local.z*say);
               total_ST[celly] = prefac_sc*(m_local.z*sax-m_local.x*saz);
               total_ST[cellz] = prefac_sc*(m_local.x*say-m_local.y*sax);

               ast[cellx] = -aj*SxSxSp[0];
               ast[celly] = -aj*SxSxSp[1];
               ast[cellz] = -aj*SxSxSp[2];

               nast[cellx] = bj*SxSp[0];
               nast[celly] = bj*SxSp[1];
               nast[cellz] = bj*SxSp[2];
               } 
            } 
         
               //int process up, 1
            for(int cell=idx; cell<idx+num_microcells_per_stack; cell++) {
               int direction_sign = cell%2==0 ? -1:1; //calculate spin current direction for Au
               
               const int cellx = 3*cell+0;
               const int celly = 3*cell+1;
               const int cellz = 3*cell+2;

               // calculate previous cell id's
               const int pcellx = 3*(cell-1)+0;
               const int pcelly = 3*(cell-1)+1;
               const int pcellz = 3*(cell-1)+2;

               // calculate next cell id's
               const int acellx = 3*(cell+1)+0;
               const int acelly = 3*(cell+1)+1;
               const int acellz = 3*(cell+1)+2;

               three_vector_t  m_local(0.0,  0.0,  0.0);
               // copy array values to temporaries for readability
               if(sot_sa_source[cell]) {
                  m_local.x = 0;
                  m_local.y = direction_sign; // Au cell artificial magnetisations
                  m_local.z = 0.0; 
               } else {
                  m_local.x = sa_int[cellx];
                  m_local.y = sa_int[celly];
                  m_local.z = sa_int[cellz]; // current cell magnetisations
               }
           
               const double modm = sqrt(m_local.x*m_local.x + m_local.y*m_local.y + m_local.z*m_local.z);
      
               // Check for zero magnetization in normalization
               if(modm > 1.e-11){
                  m_local.x = m_local.x/modm;
                  m_local.y = m_local.y/modm;
                  m_local.z = m_local.z/modm;
               }
               else{
                  m_local.x = 0.0;
                  m_local.y = 0.0;
                  m_local.z = 0.0;
               }

               //---------------------------------------------------------------------
               // Step 1 calculate coordinate transformation m -> m' for current cell
               //---------------------------------------------------------------------

               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m_local, itm);

                  // Calculate basis vectors
                  b1 = transform_vector(bz, itm);
                  b2 = transform_vector(bx, itm);
                  b3 = transform_vector(by, itm);
               }
               else{
                  // set default rotational frame
                  b1 = bz;
                  b2 = bx;
                  b3 = by;
               }
               //---------------------------------------------------------------------
               // Step 2 determine coefficients for the spin accumulation
               //---------------------------------------------------------------------

               // Initialise temporary constants
               const double Bc = sot_beta_cond[cell]; // beta
               const double Bd = sot_beta_diff[cell]; // beta_prime
               const double Do = sot_diffusion[cell];
               three_vector_t jm0(0.0,0.0,0.0);
               //grab correct cell's init spin current (int label) bounds for top and bottom cells
               if(direction_sign < 0 && cell != idx+num_microcells_per_stack-1) {jm0.x = j_int_down_y[acellx]; jm0.y = j_int_down_y[acelly]; jm0.z = j_int_down_y[acellz];}
               else if(cell != idx) {jm0.x = j_int_up_y[pcellx]; jm0.y = j_int_up_y[pcelly]; jm0.z = j_int_up_y[pcellz];}
             
               //  Calculate gradient dsacc/dx
               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m_local.x*m_local.x - twoDo;     M.xy = preD*m_local.y*m_local.x;             M.xz = preD*m_local.z*m_local.x;
               M.yx = preD*m_local.x*m_local.y;             M.yy = preD*m_local.y*m_local.y - twoDo;     M.yz = preD*m_local.z*m_local.y;
               M.zx = preD*m_local.x*m_local.z;             M.zy = preD*m_local.y*m_local.z;             M.zz = preD*m_local.z*m_local.z - twoDo;

               //no charge current polarisation
               V.x = jm0.x;// - Bc*je_eff*m_local.x;
               V.y = jm0.y;// + Bc*je_eff*m_local.y;
               V.z = jm0.z;// - Bc*je_eff*m_local.z;
              
               const three_vector_t divm_0 = gaussian_elimination(M, V);

               // Calculate mp(0), c and d
               const double i_lsdl = 1.0/sot_lambda_sdl[cell];
               double mp_inf = modm;//sa_infinity[cell];
                  if(sot_sa_source[cell]) mp_inf *= program::fractional_electric_field_strength;
               const double a_local = sot_a[cell];
               const double b_local = sot_b[cell];
               const double two_a = 2.0*a_local;
               const double two_b = 2.0*b_local;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               const three_vector_t coeff = gaussian_elimination(M, V);

               const double mp_0 = coeff.x;
               const double c    = coeff.y;
               const double d    = coeff.z;

               //------------------------------------
               // Step 3 calculate spin accumulation
               //------------------------------------

               const double x = micro_cell_thickness*1.0e-10; // Convert to metres
               const double cos_bx = cos(b_local*x);
               const double sin_bx = sin(b_local*x);
               const double e_xsdl = exp(-x*i_lsdl);
               const double e_ax   = exp(-a_local*x);
               const double prefac = (2.0*e_ax);

               const double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               const double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               const double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);
              
               //--------------------------------------------
               // Step 4 calculate the spin current (jm_end)
               //--------------------------------------------
               const double ac_bd = a_local*c + b_local*d; 
               const double ad_bc = a_local*d - b_local*c;

               const double divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl;
               const double divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx);
               const double divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx);

               const double divsax = b1.x*divsa_para + b2.x*divsa_perp2 + b3.x*divsa_perp3;
               const double divsay = b1.y*divsa_para + b2.y*divsa_perp2 + b3.y*divsa_perp3;
               const double divsaz = b1.z*divsa_para + b2.z*divsa_perp2 + b3.z*divsa_perp3;

               const double dot = m_local.x*divsax + m_local.y*divsay + m_local.z*divsaz;

               const double pre_jmx = divsax - Bc*Bd*m_local.x*dot;
               const double pre_jmy = divsay - Bc*Bd*m_local.y*dot;
               const double pre_jmz = divsaz - Bc*Bd*m_local.z*dot;

               double jmx =   twoDo*pre_jmx;
               double jmy =   twoDo*pre_jmy;
                  //Au spin source sign direction dependent 
                  if(sot_sa_source[cell]) jmy += -direction_sign*je_eff*initial_theta;
               double jmz =   twoDo*pre_jmz;
              
               if(cell_natom[cell]>0) {
                  //no sa calculation yet
                  // sa_sot_int[cellx] = sax;
                  // sa_sot_int[celly] = say;
                  // sa_sot_int[cellz] = saz;

                  // Save values for the spin current
                  if(direction_sign < 0) {
                     j_int_down_y[cellx] = jmx;
                     j_int_down_y[celly] = jmy;
                     j_int_down_y[cellz] = jmz;
                  } else {
                     j_int_up_y[cellx] = jmx;
                     j_int_up_y[celly] = jmy;
                     j_int_up_y[cellz] = jmz;
                  }
               } 
            if(output_torque_data == "int") { 
               //--------------------------------------------
               // Step 5 calculate the spin torque of each cell
               //--------------------------------------------
               // Note: intermediate spin torque value is currently sum of 
               // up and down process for up and down current. May not be a very physical value
               three_vector_t pm_local(0.0,0.0,0.0);
               double pmodm;
               if(direction_sign < 0){
                  if(sot_sa_source[cell+1]) {
                     pm_local.y = -1;
                     pmodm = 1;
                  } else {
                     pm_local.x = m[acellx];
                     pm_local.y = m[acelly];
                     pm_local.z = m[acellz];
                     pmodm = sqrt(pm_local.x*pm_local.x + pm_local.y*pm_local.y + pm_local.z*pm_local.z);
                  }
                  if(pmodm > 1.e-11){
                     pm_local.x = pm_local.x/pmodm;
                     pm_local.y = pm_local.y/pmodm;
                     pm_local.z = pm_local.z/pmodm;
                  }
                  else{
                     pm_local.x = 0.0;
                     pm_local.y = 0.0;
                     pm_local.z = 0.0;
                  }
               }
               else {
                  if(sot_sa_source[cell-1]) {
                     pm_local.y = 1;
                     pmodm = 1;
                  } else {
                     pm_local.x = m[acellx];
                     pm_local.y = m[acelly];
                     pm_local.z = m[acellz];
                     pmodm = sqrt(pm_local.x*pm_local.x + pm_local.y*pm_local.y + pm_local.z*pm_local.z);
                  }
                  if(pmodm > 1.e-11){
                     pm_local.x = pm_local.x/pmodm;
                     pm_local.y = pm_local.y/pmodm;
                     pm_local.z = pm_local.z/pmodm;
                  }
                  else{
                     pm_local.x = 0.0;
                     pm_local.y = 0.0;
                     pm_local.z = 0.0;
                  }
               }
               
               //convert M of previous cell into basis b1, b2, b3

               M.xx = b1.x;    M.xy = b2.x;    M.xz = b3.x;
               M.yx = b1.y;    M.yy = b2.y;    M.yz = b3.y;
               M.zx = b1.z;    M.zy = b2.z;    M.zz = b3.z;

               V.x = pm_local.x;
               V.y = pm_local.y;
               V.z = pm_local.z;

               const three_vector_t pm_basis = gaussian_elimination(M, V);

               //const double pm_b1 = pm_basis.x; // unused variable
               const double pm_b2 = pm_basis.y;
               const double pm_b3 = pm_basis.z;

               // Calculate the spin torque coefficients describing ast and nast
               const double prefac_sc = atomcell_volume * sd_exchange[cell] * i_e * i_muB;
               const double plus_perp =  (pm_b2*pm_b2 + pm_b3*pm_b3);
               // const double minus_perp = (pm_b2*pm_b2 - pm_b3*pm_b3); // unused variable

	            double aj; // the ST parameter describing Slonczewski torque
               double bj; // the ST parameter describing field-like torque


                if( ( plus_perp <= 1.0e-11 ) ){
                    aj = 0.0;
                    bj = 0.0; }

                else{
                    aj  = prefac_sc*(sa_perp2*pm_b3 - sa_perp3*pm_b2)/plus_perp;
                    bj  = prefac_sc*(sa_perp2*pm_b2 + sa_perp3*pm_b3)/plus_perp;
                }

               double SxSp[3], SxSxSp[3];
               coeff_ast[cell]  += aj;
               coeff_nast[cell] += bj;

               SxSp[0]=(m_local.y*pm_local.z-m_local.z*pm_local.y);
               SxSp[1]=(m_local.z*pm_local.x-m_local.x*pm_local.z);
               SxSp[2]=(m_local.x*pm_local.y-m_local.y*pm_local.x);

               SxSxSp[0]= (m_local.y*SxSp[2]-m_local.z*SxSp[1]);
               SxSxSp[1]= (m_local.z*SxSp[0]-m_local.x*SxSp[2]);
               SxSxSp[2]= (m_local.x*SxSp[1]-m_local.y*SxSp[0]);

               //convert mp and m_perp
               double sax = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
               double say = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
               double saz = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;

               //calculate directly from J(Sxm)
               total_ST[cellx] += prefac_sc*(m_local.y*saz-m_local.z*say);
               total_ST[celly] += prefac_sc*(m_local.z*sax-m_local.x*saz);
               total_ST[cellz] += prefac_sc*(m_local.x*say-m_local.y*sax);

               ast[cellx] += -aj*SxSxSp[0];
               ast[celly] += -aj*SxSxSp[1];
               ast[cellz] += -aj*SxSxSp[2];

               nast[cellx] += bj*SxSp[0];
               nast[celly] += bj*SxSp[1];
               nast[cellz] += bj*SxSp[2];
               } 
            } 
           
            //int process down 1
            for(int cell=idx+num_microcells_per_stack-1; cell >= idx; cell--) {
                  int direction_sign = cell%2==0 ? 1:-1; //sign change for the down process
               
               const int cellx = 3*cell+0;
               const int celly = 3*cell+1;
               const int cellz = 3*cell+2;

               // the sign change in the direction_sign will properly select the previous and next cells
               const int pcellx = 3*(cell-1)+0;
               const int pcelly = 3*(cell-1)+1;
               const int pcellz = 3*(cell-1)+2;

               const int acellx = 3*(cell+1)+0;
               const int acelly = 3*(cell+1)+1;
               const int acellz = 3*(cell+1)+2;

               three_vector_t  m_local(0.0,  0.0,  0.0);
               // copy array values to temporaries for readability
               if(sot_sa_source[cell]) {
                  m_local.x = 0;
                  m_local.y = direction_sign;
                  m_local.z = 0.0; // Au cell artificial magnetisations
               } else {
                  m_local.x = sa_int[cellx];
                  m_local.y = sa_int[celly];
                  m_local.z = sa_int[cellz]; // current cell magnetisations
               }
               
               const double modm = sqrt(m_local.x*m_local.x + m_local.y*m_local.y + m_local.z*m_local.z);
               
               // Check for zero magnetization in normalization
               if(modm > 1.e-11){
                  m_local.x = m_local.x/modm;
                  m_local.y = m_local.y/modm;
                  m_local.z = m_local.z/modm;
               }
               else{
                  m_local.x = 0.0;
                  m_local.y = 0.0;
                  m_local.z = 0.0;
               }

               //--------------------------------------------------------------------
               // Step 1 calculate coordinate transformation m -> m' for current cell
               //---------------------------------------------------------------------

               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m_local, itm);

                  // Calculate basis vectors
                  b1 = transform_vector(bz, itm);
                  b2 = transform_vector(bx, itm);
                  b3 = transform_vector(by, itm);
               }
               else{
                  // set default rotational frame
                  b1 = bz;
                  b2 = bx;
                  b3 = by;
               }
               //---------------------------------------------------------------------
               // Step 2 determine coefficients for the spin accumulation
               //---------------------------------------------------------------------

               // Initialise temporary constants
               const double Bc = sot_beta_cond[cell]; // beta
               const double Bd = sot_beta_diff[cell]; // beta_prime
               const double Do = sot_diffusion[cell];
               three_vector_t jm0(0.0,0.0,0.0);
               if(direction_sign < 0 && cell != idx+num_microcells_per_stack-1) {jm0.x = j_int_down_y[acellx]; jm0.y = j_int_down_y[acelly]; jm0.z = j_int_down_y[acellz];}
               else if(cell != idx) {jm0.x = j_int_up_y[pcellx]; jm0.y = j_int_up_y[pcelly]; jm0.z = j_int_up_y[pcellz];}
               
               //  Calculate gradient dsacc/dx
               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m_local.x*m_local.x - twoDo;     M.xy = preD*m_local.y*m_local.x;             M.xz = preD*m_local.z*m_local.x;
               M.yx = preD*m_local.x*m_local.y;             M.yy = preD*m_local.y*m_local.y - twoDo;     M.yz = preD*m_local.z*m_local.y;
               M.zx = preD*m_local.x*m_local.z;             M.zy = preD*m_local.y*m_local.z;             M.zz = preD*m_local.z*m_local.z - twoDo;

               //no charge current conversion
               V.x = jm0.x;// - Bc*je_eff*m_local.x;
               V.y = jm0.y;// + Bc*je_eff*m_local.y;
               V.z = jm0.z;// - Bc*je_eff*m_local.z;
               
               const three_vector_t divm_0 = gaussian_elimination(M, V);

               // Calculate mp(0), c and d
               const double i_lsdl = 1.0/sot_lambda_sdl[cell];
               double mp_inf = modm;// sa_infinity[cell];
                  if(sot_sa_source[cell]) mp_inf *= program::fractional_electric_field_strength;
               const double a_local = sot_a[cell];
               const double b_local = sot_b[cell];
               const double two_a = 2.0*a_local;
               const double two_b = 2.0*b_local;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               const three_vector_t coeff = gaussian_elimination(M, V);

               const double mp_0 = coeff.x;
               const double c    = coeff.y;
               const double d    = coeff.z;

               //------------------------------------
               // Step 3 calculate spin accumulation
               //------------------------------------

               const double x = micro_cell_thickness*1.0e-10; // Convert to metres
               const double cos_bx = cos(b_local*x);
               const double sin_bx = sin(b_local*x);
               const double e_xsdl = exp(-x*i_lsdl);
               const double e_ax   = exp(-a_local*x);
               const double prefac = (2.0*e_ax);

               const double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               const double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               const double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);
           
               //--------------------------------------------
               // Step 4 calculate the spin current (jm_end)
               //--------------------------------------------
               const double ac_bd = a_local*c + b_local*d; 
               const double ad_bc = a_local*d - b_local*c;

               const double divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl;
               const double divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx);
               const double divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx);

               const double divsax = b1.x*divsa_para + b2.x*divsa_perp2 + b3.x*divsa_perp3;
               const double divsay = b1.y*divsa_para + b2.y*divsa_perp2 + b3.y*divsa_perp3;
               const double divsaz = b1.z*divsa_para + b2.z*divsa_perp2 + b3.z*divsa_perp3;

               const double dot = m_local.x*divsax + m_local.y*divsay + m_local.z*divsaz;

               const double pre_jmx = divsax - Bc*Bd*m_local.x*dot;
               const double pre_jmy = divsay - Bc*Bd*m_local.y*dot;
               const double pre_jmz = divsaz - Bc*Bd*m_local.z*dot;

               double jmx =   twoDo*pre_jmx;
               double jmy =   twoDo*pre_jmy;
                  if(sot_sa_source[cell]) jmy = -direction_sign*je_eff*initial_theta;
               double jmz =   twoDo*pre_jmz;
               
               if(cell_natom[cell]>0) {
                  //no sa calculation yet
                  // sa_sot_final[cellx] = sax;
                  // sa_sot_final[celly] = say;
                  // sa_sot_final[cellz] = saz;

                  // Save values for the spin current
                   if(direction_sign < 0) {
                     j_int_down_y[cellx] = jmx;
                     j_int_down_y[celly] = jmy;
                     j_int_down_y[cellz] = jmz;
                  } else {
                     j_int_up_y[cellx] = jmx;
                     j_int_up_y[celly] = jmy;
                     j_int_up_y[cellz] = jmz;
                  }

               } 
               
            if(output_torque_data == "int") {
               //--------------------------------------------
               // Step 5 calculate the spin torque of each cell
               //--------------------------------------------
               three_vector_t pm_local(0.0,0.0,0.0);
               double pmodm;
               if(direction_sign < 0){
                  if(sot_sa_source[cell+1]) {
                     pm_local.y = -1;
                     pmodm = 1;
                  } else {
                     pm_local.x = m[acellx];
                     pm_local.y = m[acelly];
                     pm_local.z = m[acellz];
                     pmodm = sqrt(pm_local.x*pm_local.x + pm_local.y*pm_local.y + pm_local.z*pm_local.z);
                  }
                  if(pmodm > 1.e-11){
                     pm_local.x = pm_local.x/pmodm;
                     pm_local.y = pm_local.y/pmodm;
                     pm_local.z = pm_local.z/pmodm;
                  }
                  else{
                     pm_local.x = 0.0;
                     pm_local.y = 0.0;
                     pm_local.z = 0.0;
                  }
               }
               else {
                  if(sot_sa_source[cell-1]) {
                     pm_local.y = 1;
                     pmodm = 1;
                  } else {
                     pm_local.x = m[acellx];
                     pm_local.y = m[acelly];
                     pm_local.z = m[acellz];
                     pmodm = sqrt(pm_local.x*pm_local.x + pm_local.y*pm_local.y + pm_local.z*pm_local.z);
                  }
                  if(pmodm > 1.e-11){
                     pm_local.x = pm_local.x/pmodm;
                     pm_local.y = pm_local.y/pmodm;
                     pm_local.z = pm_local.z/pmodm;
                  }
                  else{
                     pm_local.x = 0.0;
                     pm_local.y = 0.0;
                     pm_local.z = 0.0;
                  }
               }
               //convert M of previous cell into basis b1, b2, b3

               M.xx = b1.x;    M.xy = b2.x;    M.xz = b3.x;
               M.yx = b1.y;    M.yy = b2.y;    M.yz = b3.y;
               M.zx = b1.z;    M.zy = b2.z;    M.zz = b3.z;

               V.x = pm_local.x;
               V.y = pm_local.y;
               V.z = pm_local.z;

               const three_vector_t pm_basis = gaussian_elimination(M, V);

               //const double pm_b1 = pm_basis.x; // unused variable
               const double pm_b2 = pm_basis.y;
               const double pm_b3 = pm_basis.z;

               // Calculate the spin torque coefficients describing ast and nast
               const double prefac_sc = microcell_volume * sd_exchange[cell] * i_e * i_muB;
               const double plus_perp =  (pm_b2*pm_b2 + pm_b3*pm_b3);
               // const double minus_perp = (pm_b2*pm_b2 - pm_b3*pm_b3); // unused variable

	            double aj; // the ST parameter describing Slonczewski torque
               double bj; // the ST parameter describing field-like torque


                if( ( plus_perp <= 1.0e-11 ) ){
                    aj = 0.0;
                    bj = 0.0; }

                else{
                    aj  = prefac_sc*(sa_perp2*pm_b3 - sa_perp3*pm_b2)/plus_perp;
                    bj  = prefac_sc*(sa_perp2*pm_b2 + sa_perp3*pm_b3)/plus_perp;
                }

               double SxSp[3], SxSxSp[3];
               coeff_ast[cell]  += aj;
               coeff_nast[cell] += bj;

               SxSp[0]=(m_local.y*pm_local.z-m_local.z*pm_local.y);
               SxSp[1]=(m_local.z*pm_local.x-m_local.x*pm_local.z);
               SxSp[2]=(m_local.x*pm_local.y-m_local.y*pm_local.x);

               SxSxSp[0]= (m_local.y*SxSp[2]-m_local.z*SxSp[1]);
               SxSxSp[1]= (m_local.z*SxSp[0]-m_local.x*SxSp[2]);
               SxSxSp[2]= (m_local.x*SxSp[1]-m_local.y*SxSp[0]);


               //convert mp and m_perp
               double sax = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
               double say = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
               double saz = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;

               //calculate directly from J(Sxm)
               total_ST[cellx] += prefac_sc*(m_local.y*saz-m_local.z*say);
               total_ST[celly] += prefac_sc*(m_local.z*sax-m_local.x*saz);
               total_ST[cellz] += prefac_sc*(m_local.x*say-m_local.y*sax);

               ast[cellx] += -aj*SxSxSp[0];
               ast[celly] += -aj*SxSxSp[1];
               ast[cellz] += -aj*SxSxSp[2];

               nast[cellx] += bj*SxSp[0];
               nast[celly] += bj*SxSp[1];
               nast[cellz] += bj*SxSp[2];
               } 
            }
       
            //int process up, 2
            // second up sweep to use the newly calculated down sweep values
            for(int cell=idx; cell<idx+num_microcells_per_stack; cell++) {
                  int direction_sign = cell%2==0 ? -1:1;
               
               const int cellx = 3*cell+0;
               const int celly = 3*cell+1;
               const int cellz = 3*cell+2;

               // calculate previous cell id's
               const int pcellx = 3*(cell-1)+0;
               const int pcelly = 3*(cell-1)+1;
               const int pcellz = 3*(cell-1)+2;

               const int acellx = 3*(cell+1)+0;
               const int acelly = 3*(cell+1)+1;
               const int acellz = 3*(cell+1)+2;

               three_vector_t  m_local(0.0,  0.0,  0.0);
               // copy array values to temporaries for readability
               if(sot_sa_source[cell]) {
                  m_local.x = 0;
                  m_local.y = direction_sign;
                  m_local.z = 0.0; // Au cell artificial magnetisations
               } else {
                  m_local.x = sa_int[cellx];
                  m_local.y = sa_int[celly];
                  m_local.z = sa_int[cellz]; // current cell magnetisations
               }
           
               const double modm = sqrt(m_local.x*m_local.x + m_local.y*m_local.y + m_local.z*m_local.z);
               
               // Check for zero magnetization in normalization
               if(modm > 1.e-11){
                  m_local.x = m_local.x/modm;
                  m_local.y = m_local.y/modm;
                  m_local.z = m_local.z/modm;
               }
               else{
                  m_local.x = 0.0;
                  m_local.y = 0.0;
                  m_local.z = 0.0;
               }

               //---------------------------------------------------------------------
               // Step 1 calculate coordinate transformation m -> m' for current cell
               //---------------------------------------------------------------------

               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m_local, itm);

                  // Calculate basis vectors
                  b1 = transform_vector(bz, itm);
                  b2 = transform_vector(bx, itm);
                  b3 = transform_vector(by, itm);
               }
               else{
                  // set default rotational frame
                  b1 = bz;
                  b2 = bx;
                  b3 = by;
               }
               //---------------------------------------------------------------------
               // Step 2 determine coefficients for the spin accumulation
               //---------------------------------------------------------------------

               // Initialise temporary constants
               const double Bc = sot_beta_cond[cell]; // beta
               const double Bd = sot_beta_diff[cell]; // beta_prime
               const double Do = sot_diffusion[cell];
               three_vector_t jm0(0.0,0.0,0.0);
               if(direction_sign < 0 && cell != idx+num_microcells_per_stack-1) {jm0.x = j_int_down_y[acellx]; jm0.y = j_int_down_y[acelly]; jm0.z = j_int_down_y[acellz];}
               else if(cell != idx) {jm0.x = j_int_up_y[pcellx]; jm0.y = j_int_up_y[pcelly]; jm0.z = j_int_up_y[pcellz];}
              
               //  Calculate gradient dsacc/dx
               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m_local.x*m_local.x - twoDo;     M.xy = preD*m_local.y*m_local.x;             M.xz = preD*m_local.z*m_local.x;
               M.yx = preD*m_local.x*m_local.y;             M.yy = preD*m_local.y*m_local.y - twoDo;     M.yz = preD*m_local.z*m_local.y;
               M.zx = preD*m_local.x*m_local.z;             M.zy = preD*m_local.y*m_local.z;             M.zz = preD*m_local.z*m_local.z - twoDo;

               // no charge current conversion
                  V.x = jm0.x;// - Bc*je_eff*m_local.x;
                  V.y = jm0.y;// + Bc*je_eff*m_local.y;
                  V.z = jm0.z;// - Bc*je_eff*m_local.z;
              
               const three_vector_t divm_0 = gaussian_elimination(M, V);

               // Calculate mp(0), c and d
               const double i_lsdl = 1.0/sot_lambda_sdl[cell];
               double mp_inf = modm;// sa_infinity[cell];
                  if(sot_sa_source[cell]) mp_inf *= program::fractional_electric_field_strength;
               const double a_local = sot_a[cell];
               const double b_local = sot_b[cell];
               const double two_a = 2.0*a_local;
               const double two_b = 2.0*b_local;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               const three_vector_t coeff = gaussian_elimination(M, V);

               const double mp_0 = coeff.x;
               const double c    = coeff.y;
               const double d    = coeff.z;

               //------------------------------------
               // Step 3 calculate spin accumulation
               //------------------------------------

               const double x = micro_cell_thickness*1.0e-10; // Convert to metres
               const double cos_bx = cos(b_local*x);
               const double sin_bx = sin(b_local*x);
               const double e_xsdl = exp(-x*i_lsdl);
               const double e_ax   = exp(-a_local*x);
               const double prefac = (2.0*e_ax);

               const double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               const double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               const double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);

               //--------------------------------------------
               // Step 4 calculate the spin current (jm_end)
               //--------------------------------------------
               const double ac_bd = a_local*c + b_local*d; 
               const double ad_bc = a_local*d - b_local*c;

               const double divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl;
               const double divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx);
               const double divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx);

               const double divsax = b1.x*divsa_para + b2.x*divsa_perp2 + b3.x*divsa_perp3;
               const double divsay = b1.y*divsa_para + b2.y*divsa_perp2 + b3.y*divsa_perp3;
               const double divsaz = b1.z*divsa_para + b2.z*divsa_perp2 + b3.z*divsa_perp3;

               const double dot = m_local.x*divsax + m_local.y*divsay + m_local.z*divsaz;

               const double pre_jmx = divsax - Bc*Bd*m_local.x*dot;
               const double pre_jmy = divsay - Bc*Bd*m_local.y*dot;
               const double pre_jmz = divsaz - Bc*Bd*m_local.z*dot;

               double jmx =   twoDo*pre_jmx;
               double jmy =   twoDo*pre_jmy;
                  if(sot_sa_source[cell]) jmy += -direction_sign*je_eff*initial_theta;
               double jmz =   twoDo*pre_jmz;

               if(cell_natom[cell]>0) {

                  // sa_sot_int[cellx] = sax;
                  // sa_sot_int[celly] = say;
                  // sa_sot_int[cellz] = saz;

                  // Save values for the spin current
                  if(direction_sign < 0) {
                     j_int_down_y[cellx] = jmx;
                     j_int_down_y[celly] = jmy;
                     j_int_down_y[cellz] = jmz;
                  } else {
                     j_int_up_y[cellx] = jmx;
                     j_int_up_y[celly] = jmy;
                     j_int_up_y[cellz] = jmz;
                  }
               } 
            
            if(output_torque_data == "int") {
               //--------------------------------------------
               // Step 5 calculate the spin torque of each cell
               //--------------------------------------------
               three_vector_t pm_local(0.0,0.0,0.0);
               double pmodm;
               if(direction_sign < 0){
                  if(sot_sa_source[cell+1]) {
                     pm_local.y = -1;
                     pmodm = 1;
                  } else {
                     pm_local.x = m[acellx];
                     pm_local.y = m[acelly];
                     pm_local.z = m[acellz];
                     pmodm = sqrt(pm_local.x*pm_local.x + pm_local.y*pm_local.y + pm_local.z*pm_local.z);
                  }
                  if(pmodm > 1.e-11){
                     pm_local.x = pm_local.x/pmodm;
                     pm_local.y = pm_local.y/pmodm;
                     pm_local.z = pm_local.z/pmodm;
                  }
                  else{
                     pm_local.x = 0.0;
                     pm_local.y = 0.0;
                     pm_local.z = 0.0;
                  }
               }
               else {
                  if(sot_sa_source[cell-1]) {
                     pm_local.y = 1;
                     pmodm = 1;
                  } else {
                     pm_local.x = m[acellx];
                     pm_local.y = m[acelly];
                     pm_local.z = m[acellz];
                     pmodm = sqrt(pm_local.x*pm_local.x + pm_local.y*pm_local.y + pm_local.z*pm_local.z);
                  }
                  if(pmodm > 1.e-11){
                     pm_local.x = pm_local.x/pmodm;
                     pm_local.y = pm_local.y/pmodm;
                     pm_local.z = pm_local.z/pmodm;
                  }
                  else{
                     pm_local.x = 0.0;
                     pm_local.y = 0.0;
                     pm_local.z = 0.0;
                  }
               }
               //convert M of previous cell into basis b1, b2, b3

               M.xx = b1.x;    M.xy = b2.x;    M.xz = b3.x;
               M.yx = b1.y;    M.yy = b2.y;    M.yz = b3.y;
               M.zx = b1.z;    M.zy = b2.z;    M.zz = b3.z;

               V.x = pm_local.x;
               V.y = pm_local.y;
               V.z = pm_local.z;

               const three_vector_t pm_basis = gaussian_elimination(M, V);

               //const double pm_b1 = pm_basis.x; // unused variable
               const double pm_b2 = pm_basis.y;
               const double pm_b3 = pm_basis.z;

               // Calculate the spin torque coefficients describing ast and nast
               const double prefac_sc = atomcell_volume * sd_exchange[cell] * i_e * i_muB;
               const double plus_perp =  (pm_b2*pm_b2 + pm_b3*pm_b3);
               // const double minus_perp = (pm_b2*pm_b2 - pm_b3*pm_b3); // unused variable

	            double aj; // the ST parameter describing Slonczewski torque
               double bj; // the ST parameter describing field-like torque


                if( ( plus_perp <= 1.0e-11 ) ){
                    aj = 0.0;
                    bj = 0.0; }

                else{
                    aj  = prefac_sc*(sa_perp2*pm_b3 - sa_perp3*pm_b2)/plus_perp;
                    bj  = prefac_sc*(sa_perp2*pm_b2 + sa_perp3*pm_b3)/plus_perp;
                }

               double SxSp[3], SxSxSp[3];
               coeff_ast[cell]  += aj;
               coeff_nast[cell] += bj;

               SxSp[0]=(m_local.y*pm_local.z-m_local.z*pm_local.y);
               SxSp[1]=(m_local.z*pm_local.x-m_local.x*pm_local.z);
               SxSp[2]=(m_local.x*pm_local.y-m_local.y*pm_local.x);

               SxSxSp[0]= (m_local.y*SxSp[2]-m_local.z*SxSp[1]);
               SxSxSp[1]= (m_local.z*SxSp[0]-m_local.x*SxSp[2]);
               SxSxSp[2]= (m_local.x*SxSp[1]-m_local.y*SxSp[0]);


               //convert mp and m_perp 
               double sax = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
               double say = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
               double saz = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;
         
               //calculate directly from J(Sxm)
               total_ST[cellx] += prefac_sc*(m_local.y*saz-m_local.z*say);
               total_ST[celly] += prefac_sc*(m_local.z*sax-m_local.x*saz);
               total_ST[cellz] += prefac_sc*(m_local.x*say-m_local.y*sax);

               ast[cellx] += -aj*SxSxSp[0];
               ast[celly] += -aj*SxSxSp[1];
               ast[cellz] += -aj*SxSxSp[2];

               nast[cellx] += bj*SxSp[0];
               nast[celly] += bj*SxSp[1];
               nast[cellz] += bj*SxSp[2];
               } 
            } 
            //summation step
            for(int cell=idx; cell<idx+num_microcells_per_stack; cell++) {
               
                //direction sign for Au outgoing spin current for data output
               int direction_sign = sot_sa_source[cell] ? -1:1;
               int cellx = 3*cell;
               int celly = 3*cell+1;
               int cellz = 3*cell+2;

               three_vector_t  m_local(sa_int[cellx],  sa_int[celly],  sa_int[cellz]);
               double modm = sqrt(m_local.x*m_local.x + m_local.y*m_local.y + m_local.z*m_local.z);
               if(modm > 1.e-11){
                  m_local.x = m_local.x/modm;
                  m_local.y = m_local.y/modm;
                  m_local.z = m_local.z/modm;
               }
               else{
                  m_local.x = 0.0;
                  m_local.y = 1.0; //default to +y for Au, get alternating sa in output but no physical change
                  m_local.z = 0.0;
                  modm = 1.0;
               }
               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m_local, itm);

                  // Calculate basis vectors
                  b1 = transform_vector(bz, itm);
                  b2 = transform_vector(bx, itm);
                  b3 = transform_vector(by, itm);
               }
               else{
                  // set default rotational frame
                  b1 = bz;
                  b2 = bx;
                  b3 = by;
               }

                // Initialise temporary constants
               const double Bc = sot_beta_cond[cell]; // beta
               const double Bd = sot_beta_diff[cell]; // beta_prime
               const double Do = sot_diffusion[cell];                   // no contribution from x current since it was used to calculate the sa_int value at beginning of z axis steps
               three_vector_t jm0(j_int_up_y[cellx]+j_int_down_y[cellx],//+j_final_up_x[cellx], 
                                  j_int_up_y[celly]+j_int_down_y[celly],//+j_final_up_x[cellx],
                                  j_int_up_y[cellz]+j_int_down_y[cellz]);//+j_final_up_x[cellx]);                     

               
               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m_local.x*m_local.x - twoDo;     M.xy = preD*m_local.y*m_local.x;             M.xz = preD*m_local.z*m_local.x;
               M.yx = preD*m_local.x*m_local.y;             M.yy = preD*m_local.y*m_local.y - twoDo;     M.yz = preD*m_local.z*m_local.y;
               M.zx = preD*m_local.x*m_local.z;             M.zy = preD*m_local.y*m_local.z;             M.zz = preD*m_local.z*m_local.z - twoDo;

               V.x = jm0.x;// - Bc*je_eff*m_local.x;
               V.y = jm0.y;//- Bc*je_eff*m_local.y;
               V.z = jm0.z;// - Bc*je_eff*m_local.z;

               const three_vector_t divm_0 = gaussian_elimination(M, V);

               const double i_lsdl = 1.0/sot_lambda_sdl[cell];
               double mp_inf = modm;// sa_infinity[cell];
                  if(sot_sa_source[cell]) mp_inf *= program::fractional_electric_field_strength;
               const double a_local = sot_a[cell];
               const double b_local = sot_b[cell];
               const double two_a = 2.0*a_local;
               const double two_b = 2.0*b_local;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               const three_vector_t coeff = gaussian_elimination(M, V);

               const double mp_0 = coeff.x;
               const double c    = coeff.y;
               const double d    = coeff.z;

               const double x = micro_cell_thickness*1.0e-10; // Convert to metres
               const double cos_bx = cos(b_local*x);
               const double sin_bx = sin(b_local*x);
               const double e_xsdl = exp(-x*i_lsdl);
               const double e_ax   = exp(-a_local*x);
               const double prefac = (2.0*e_ax);

               const double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               const double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               const double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);

               //convert mp and m_perp
               double sax = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
               double say = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
               double saz = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;

               const double ac_bd = a_local*c + b_local*d; 
               const double ad_bc = a_local*d - b_local*c;

               const double divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl;
               const double divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx);
               const double divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx);

               const double divsax = b1.x*divsa_para + b2.x*divsa_perp2 + b3.x*divsa_perp3;
               const double divsay = b1.y*divsa_para + b2.y*divsa_perp2 + b3.y*divsa_perp3;
               const double divsaz = b1.z*divsa_para + b2.z*divsa_perp2 + b3.z*divsa_perp3;

               const double dot = m_local.x*divsax + m_local.y*divsay + m_local.z*divsaz;

               const double pre_jmx = divsax - Bc*Bd*m_local.x*dot;
               const double pre_jmy = divsay - Bc*Bd*m_local.y*dot;
               const double pre_jmz = divsaz - Bc*Bd*m_local.z*dot;

               double jmx =   twoDo*pre_jmx;
               double jmy =   twoDo*pre_jmy;
                  if(sot_sa_source[cell]) jmy = je_eff*initial_theta; //depreciated
               double jmz =   twoDo*pre_jmz;

               sa_final[cellx] = sax;
               sa_final[celly] = say;
               sa_final[cellz] = saz;

               //sot check output flag will run spin acc program without torque to calculate and output values
               if(!sot_check){
                  spin_torque[cellx] = atomcell_volume * sd_exchange[cell] * (sax) * i_e * i_muB;
                  spin_torque[celly] = atomcell_volume * sd_exchange[cell] * (say) * i_e * i_muB;
                  spin_torque[cellz] = atomcell_volume * sd_exchange[cell] * (saz) * i_e * i_muB; 
               }
               j_final_up_y[cellx] = jmx;
               j_final_up_y[celly] = jmy;
               j_final_up_y[cellz] = jmz;

               j_final_down_y[cellx] = jmx*direction_sign;
               j_final_down_y[celly] = jmy*direction_sign;
               j_final_down_y[cellz] = jmz*direction_sign;

               if(output_torque_data == "final") {
               //--------------------------------------------
               // Step 5 calculate the spin torque of each cell
               //--------------------------------------------
               three_vector_t pm_local(0.0,0.0,0.0);
               double pmodm;
               
             
                     pm_local.x = 0.0;
                     pm_local.y = 1.0;
                     pm_local.z = 0.0;
                  
                set_inverse_transformation_matrix(pm_local, itm);

                  // Calculate basis vectors
                  b1 = transform_vector(bz, itm);
                  b2 = transform_vector(bx, itm);
                  b3 = transform_vector(by, itm);
               //convert M of previous cell into basis b1, b2, b3

               M.xx = b1.x;    M.xy = b2.x;    M.xz = b3.x;
               M.yx = b1.y;    M.yy = b2.y;    M.yz = b3.y;
               M.zx = b1.z;    M.zy = b2.z;    M.zz = b3.z;

               V.x = pm_local.x;
               V.y = pm_local.y;
               V.z = pm_local.z;

               const three_vector_t pm_basis = gaussian_elimination(M, V);

               //const double pm_b1 = pm_basis.x; // unused variable
               const double pm_b2 = pm_basis.y;
               const double pm_b3 = pm_basis.z;

               // Calculate the spin torque coefficients describing ast and nast
               const double prefac_sc = atomcell_volume * sot_sd_exchange[cell] * i_e * i_muB;
               const double plus_perp =  (pm_b2*pm_b2 + pm_b3*pm_b3);
               // const double minus_perp = (pm_b2*pm_b2 - pm_b3*pm_b3); // unused variable

	            double aj; // the ST parameter describing Slonczewski torque
               double bj; // the ST parameter describing field-like torque


                if( ( plus_perp <= 1.0e-11 ) ){
                    aj = 0.0;
                    bj = 0.0; }

                else{
                    aj  = prefac_sc*(sa_perp2*pm_b3 - sa_perp3*pm_b2)/plus_perp;
                    bj  = prefac_sc*(sa_perp2*pm_b2 + sa_perp3*pm_b3)/plus_perp;
                }

                double SxSp[3], SxSxSp[3];
                coeff_ast[cell]  = aj;
                coeff_nast[cell] = bj;

                SxSp[0]=(m_local.y*pm_local.z-m_local.z*pm_local.y);
                SxSp[1]=(m_local.z*pm_local.x-m_local.x*pm_local.z);
                SxSp[2]=(m_local.x*pm_local.y-m_local.y*pm_local.x);

                SxSxSp[0]= (m_local.y*SxSp[2]-m_local.z*SxSp[1]);
                SxSxSp[1]= (m_local.z*SxSp[0]-m_local.x*SxSp[2]);
                SxSxSp[2]= (m_local.x*SxSp[1]-m_local.y*SxSp[0]);

                //calculate directly from J(Sxm)
                double mlocal[3] = {m[3*cell], m[3*cell+1], m[3*cell+2]};
                double mmod = sqrt(mlocal[0]*mlocal[0] + mlocal[1]*mlocal[1] + mlocal[2]*mlocal[2]);
                if(mmod > 0.0) {mlocal[0] /= mmod; mlocal[1] /= mmod; mlocal[2] /= mmod;}
                double Tx = prefac_sc*(mlocal[1]*saz-mlocal[2]*say);
                double Ty = prefac_sc*(mlocal[2]*sax-mlocal[0]*saz);
                double Tz = prefac_sc*(mlocal[0]*say-mlocal[1]*sax);
                //Torque field in T. hard code in output to scale with m_s 3.72. Need to change
                total_ST[cellx] = Ty*mlocal[2]-Tz*mlocal[1];
                total_ST[celly] = Tz*mlocal[0]-Tx*mlocal[2];
                total_ST[cellz] = Tx*mlocal[1]-Ty*mlocal[0];

                ast[cellx] = -aj*SxSxSp[0];
                ast[celly] = -aj*SxSxSp[1];
                ast[cellz] = -aj*SxSxSp[2];

                nast[cellx] = bj*SxSp[0];
                nast[celly] = bj*SxSp[1];
                nast[cellz] = bj*SxSp[2];
               }
            }
         } // end of stack loop

        
         // #ifdef MPICH
            // int_stacks = mpi_stack_list_y.size();
         // #else 
            // int_stacks = stack_index_y.size();
         // #endif
         // std::cout << "finilasing sot simulation." << std::endl;
         // for(int s = 0; s < int_stacks; s++) {
         //    // #ifdef MPICH
         //       stack = mpi_stack_list_y.at(s);
         //    // #else 
         //       // stack = stack_index_y.at(s);
         //    // #endif
         //    const int idx = stack_index_y[stack] + 1;
         //   //final loop to calculate spin acc and torque from summed spin currents 
         // }
         

         // Reduce all atomcell spin torques on all nodes
         //only do reduction for this constant since used in LLG step. Others are reduced in output step
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &spin_torque[0],spin_torque.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif
         output_microcell_data();

         st::spin_acc_time += stopwatch.elapsed_seconds();

         return;
      }

   } // end of namespace internal
} // end of namespace st
