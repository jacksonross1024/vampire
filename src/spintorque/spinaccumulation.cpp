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
      void calculate_spin_accumulation(){
         if(st::internal::sot_sa) std::cout << "standard spin acc and sot spin acc activated " << std::endl;
         stopwatch_t stopwatch;
         stopwatch.start();

         // Zero all shared arrays (essential for parallelisation to work)
         std::fill (st::internal::spin_torque.begin(),st::internal::spin_torque.end(),0.0);
         std::fill (st::internal::sa.begin(),st::internal::sa.end(),0.0);
         // std::fill (st::internal::sa_sot.begin(),st::internal::sa_sot.end(),0.0);
         // std::fill (st::internal::j.begin(),st::internal::j.end(),0.0);
         // Declare resuable temporary variables
         st::internal::matrix_t itm; // inverse transformation matrix
         st::internal::matrix_t M; // general matrix
         st::internal::three_vector_t V(0.0,0.0,0.0); // general 3-vector

         // reference basis vectors
         const st::internal::three_vector_t bx(1.0,0.0,0.0);
         const st::internal::three_vector_t by(0.0,1.0,0.0);
         const st::internal::three_vector_t bz(0.0,0.0,1.0);

         // local basis vectors
         st::internal::three_vector_t b1(1.0,0.0,0.0);
         st::internal::three_vector_t b2(0.0,1.0,0.0);
         st::internal::three_vector_t b3(0.0,0.0,1.0);

         // set local constants
         double je = program::fractional_electric_field_strength* st::internal::je; // current (C/s) 

         //---------------------------------------------------------------------------------------------------
         //set parameters for TMR calculation
         if(st::internal::TMRenable == true){
            std::cout << "wrong" << std::endl;
            int FL =	st::internal::free_layer;
            int RL =	st::internal::reference_layer;
            double dot = st::internal::magx_mat[RL]*st::internal::magx_mat[FL]+
                         st::internal::magy_mat[RL]*st::internal::magy_mat[FL]+
                         st::internal::magz_mat[RL]*st::internal::magz_mat[FL];

            // RE - this code is not general! Needs to be fixed. Placeholder added in the meantime
            // AM (2020) - Code fixed, but still codes refers to MTJ RL/barrier/FL specifically
            double MgO_thickness = (create::get_material_height_min(FL)-create::get_material_height_max(RL))*cs::system_dimensions[2]*1.0e-10;

            //calculate the relative angle of two FMs
            st::internal::rel_angle = acos(dot);
            double plus_cos = 1.0+cos(st::internal::rel_angle);
            double minus_cos = 1.0-cos(st::internal::rel_angle);
            double exp_t = exp(-MgO_thickness/0.25e-9);

            double jtunnel = st::internal::je*0.5*(plus_cos+0.5*minus_cos)*exp_t;
//            std::cout << "t_MgO=( " << create::get_material_height_min(FL) << " - " << create::get_material_height_max(RL) << " ) = " << MgO_thickness << "\tJe\t" << st::internal::je << "\tJe_tun\t" << jtunnel << std::endl;

            //set the current je and spin poralisation parameters
            je = jtunnel;
            // AM (2020) - I think the default parameters should be rescaled by same factor as tunnelling current and not changed using those of material 0 arbitrarily
            st::internal::default_properties.beta_cond *= /*st::internal::mp[0].beta_cond**/0.5*(plus_cos+0.5*minus_cos)*exp_t;
            st::internal::default_properties.beta_diff *= /*st::internal::mp[0].beta_diff**/0.5*(plus_cos+0.5*minus_cos)*exp_t;

            // Calculate spin torque parameters
            for(size_t cell=0; cell<st::internal::beta_cond.size(); ++cell){

               // check for zero atoms in cell
               if(st::internal::cell_natom[cell] <= 0.0001){
                  st::internal::beta_cond[cell]   = st::internal::default_properties.beta_cond;
                  st::internal::beta_diff[cell]   = st::internal::default_properties.beta_diff;

                  const double hbar = 1.05457162e-34;
                  const double B  = st::internal::beta_cond[cell];
                  const double Bp = st::internal::beta_diff[cell];
                  const double lambda_sdl = st::internal::lambda_sdl[cell];
                  const double Do = st::internal::diffusion[cell];
                  const double Jsd = st::internal::sd_exchange[cell];

                  const double BBp = 1.0/sqrt(1.0-B*Bp);
                  const double lambda_sf = lambda_sdl*BBp;
                  const double lambda_j = sqrt(2.0*hbar*Do/Jsd); // Angstroms
                  const double lambda_sf2 = lambda_sf*lambda_sf;
                  const double lambda_j2 = lambda_j*lambda_j;

                  std::complex<double> inside (1.0/lambda_sf2, -1.0/lambda_j2);
                  std::complex<double> inv_lplus = sqrt(inside);

                  st::internal::a[cell] =  real(inv_lplus);
                  st::internal::b[cell] = -imag(inv_lplus);
               }
            }
            st::internal::output_base_microcell_data();
         }

         //---------------------------------------------------------------------------------------------------

         const double i_muB = 1.0/9.274e-24; // J/T
         const double i_e = 1.0/1.60217662e-19; // electronic charge (Coulombs)
         const double microcell_volume = (st::internal::micro_cell_size[st::internal::stx] *
                                          st::internal::micro_cell_size[st::internal::sty] *
                                          st::internal::micro_cell_thickness)*1.e-30; // m^3

         // loop over all 1D stacks (in parallel)
         int int_stacks = st::internal::mpi_stack_list.size();
         // std::vector<int> stacks_list;
         // stacks_list = st::internal::mpi_stack_list           
        

         int stack = 0;
         // std::cout << int_stacks << ", " << st::internal::mpi_stack_list.size() << std::endl;
         for(int s=0; s < int_stacks; ++s) {
            stack = st::internal::mpi_stack_list.at(s);
        
            // if(stack % 3 == 0) continue;
               // std::cout << stack << std::endl;
               // std::cout << st::internal::spin_acc_sign[stack_index[stack]] << ", " << stack << std::endl;
      
            // std::cout << vmpi::my_rank << ", " << stack << std::endl;
            // determine starting cell in stack
            const int idx = stack_index[stack] + 1;
           // std::cout << stack << ", " << idx << std::endl;
            // set initial values
           if(st::internal::fbc) { 
               st::internal::sa[3*idx+0] = 0.0;
               st::internal::sa[3*idx+1] = 0.0;
               st::internal::sa[3*idx+2] = 0.0;
               st::internal::j [3*idx+0] = st::internal::initial_beta*je*st::internal::initial_m[0];
               st::internal::j [3*idx+1] = st::internal::initial_beta*je*st::internal::initial_m[1];
               st::internal::j [3*idx+2] = st::internal::initial_beta*je*st::internal::initial_m[2];
           }  else if (!st::internal::fbc) {
               const double mod = 1.0/sqrt(st::internal::m [3*idx+0]*st::internal::m [3*idx+0] + st::internal::m [3*idx+1]*st::internal::m [3*idx+1] + st::internal::m [3*idx+2]*st::internal::m[3*idx+2]);
               st::internal::sa[3*idx+0] = st::internal::sa_infinity[idx]*st::internal::m [3*idx+0]*mod;
               st::internal::sa[3*idx+1] = st::internal::sa_infinity[idx]*st::internal::m [3*idx+1]*mod;
               st::internal::sa[3*idx+2] = st::internal::sa_infinity[idx]*st::internal::m [3*idx+2]*mod;

               st::internal::j [3*idx+0] = st::internal::initial_beta*je*st::internal::m [3*idx+0];
               st::internal::j [3*idx+1] = st::internal::initial_beta*je*st::internal::m [3*idx+1];
               st::internal::j [3*idx+2] = st::internal::initial_beta*je*st::internal::m [3*idx+2];
           }

           if(st::internal::sot_sa && st::internal::sot_sa_source[idx]){
               st::internal::sa[3*idx+0] = st::internal::sa_infinity[idx]*0.0;
               st::internal::sa[3*idx+1] = program::fractional_electric_field_strength*st::internal::sa_infinity[idx];
               st::internal::sa[3*idx+2] = st::internal::sa_infinity[idx]*0.0;

               // std::cout << st::internal::sa[3*idx+0] << ", " << st::internal::sa[3*idx+1] << ", " << st::internal::sa[3*idx+2] << std::endl;

               st::internal::j [3*idx+0] = st::internal::initial_beta*je*0.0;
               st::internal::j [3*idx+1] = 0.0;//st::internal::initial_beta*je*1.0;
               st::internal::j [3*idx+2] = st::internal::initial_beta*je*0.0;
           }

          //  if(sim::time%(ST_output_rate) ==0) std::cout<< stack << "\t" << st::internal::default_properties.sa_infinity << "\t" << st::internal::init_stack_mag[((stack)%6)*3 + 0]/sqrt(2.0) << "\t" << st::internal::init_stack_mag[((stack)%6)*3 + 1]/sqrt(2.0) << "\t" << st::internal::m[3*idx+0] << " \t" <<  st::internal::m[3*idx+1] << "\t" << st::internal::sa[3*idx+0] << "\t" << st::internal::sa[3*idx+1] << std::endl;

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

               st::internal::three_vector_t  m(0.0,  0.0,  0.0);
               // copy array values to temporaries for readability
               if(st::internal::sot_sa && st::internal::sot_sa_source[cell]) {
                  m.x = 0;
                  m.y = -1.0;
                  m.z = 0.0; // Au cell artificial magnetisations
               } else {
                  m.x =st::internal::m[cellx];
                  m.y = st::internal::m[celly];
                  m.z = st::internal::m[cellz]; // current cell magnetisations
               }
               st::internal::three_vector_t pm(0.0,0.0,0.0);
               if(st::internal::sot_sa && st::internal::sot_sa_source[cell-1]) {
                  pm.x = 0;
                  pm.y = 1.0;
                  pm.z = 0.0; // Au cell artificial magnetisations
               } else {
                  pm.x =st::internal::m[pcellx];
                  pm.y = st::internal::m[pcelly];
                  pm.z = st::internal::m[pcellz]; // current cell magnetisations
               }

               const double modm = sqrt(m.x*m.x + m.y*m.y + m.z*m.z);
               const double pmodm = sqrt(pm.x*pm.x + pm.y*pm.y + pm.z*pm.z);

               // Check for zero magnetization in normalization
               if(modm > 1.e-11){
                  m.x = m.x/modm;
                  m.y = m.y/modm;
                  m.z = m.z/modm;
               }
               else{
                  m.x = 0.0;
                  m.y = 0.0;
                  m.z = 0.0;
               }
               if(pmodm > 1.e-11){
                  pm.x = pm.x/pmodm;
                  pm.y = pm.y/pmodm;
                  pm.z = pm.z/pmodm;
               }
               else{
                  pm.x = 0.0;
                  pm.y = 0.0;
                  pm.z = 0.0;
               }

               //---------------------------------------------------------------------
               // Step 1 calculate coordinate transformation m -> m' for current cell
               //---------------------------------------------------------------------

               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m, itm);

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
               const double Bc = st::internal::beta_cond[cell]; // beta
               const double Bd = st::internal::beta_diff[cell]; // beta_prime
               const double Do = st::internal::diffusion[cell];
                st::internal::three_vector_t jm0(st::internal::j[pcellx],st::internal::j[pcelly],st::internal::j[pcellz]);
               if(st::internal::sot_sa_source[cell-1]) {
                  jm0.y += st::internal::initial_beta*je;
               }

               //  Calculate gradient dsacc/dx
               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m.x*m.x - twoDo;     M.xy = preD*m.y*m.x;             M.xz = preD*m.z*m.x;
               M.yx = preD*m.x*m.y;             M.yy = preD*m.y*m.y - twoDo;     M.yz = preD*m.z*m.y;
               M.zx = preD*m.x*m.z;             M.zy = preD*m.y*m.z;             M.zz = preD*m.z*m.z - twoDo;

               if(st::internal::sot_sa ) {
                  V.x = jm0.x;// - Bc*je*m.x;
                  V.y = jm0.y;// + Bc*je*m.y;
                  V.z = jm0.z;// - Bc*je*m.z;
               } else {
                  V.x = jm0.x - Bc*je*m.x;
                  V.y = jm0.y - Bc*je*m.y;
                  V.z = jm0.z - Bc*je*m.z;
               }
               const st::internal::three_vector_t divm_0 = gaussian_elimination(M, V);

               // Calculate mp(0), c and d
               const double i_lsdl = 1.0/st::internal::lambda_sdl[cell];
                double mp_inf = st::internal::sa_infinity[cell];
                 if(st::internal::sot_sa && st::internal::sot_sa_source[cell]) mp_inf *= program::fractional_electric_field_strength;
               const double a = st::internal::a[cell];
               const double b = st::internal::b[cell];
               const double two_a = 2.0*a;
               const double two_b = 2.0*b;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               const st::internal::three_vector_t coeff = gaussian_elimination(M, V);

               const double mp_0 = coeff.x;
               const double c    = coeff.y;
               const double d    = coeff.z;

               //------------------------------------
               // Step 3 calculate spin accumulation
               //------------------------------------

               const double x = st::internal::micro_cell_thickness*1.0e-10; // Convert to metres
               const double cos_bx = cos(b*x);
               const double sin_bx = sin(b*x);
               const double e_xsdl = exp(-x*i_lsdl);
               const double e_ax   = exp(-a*x);
               const double prefac = (2.0*e_ax);

               const double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               const double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               const double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);

               //convert mp and m_perp
                double sax = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
                double say = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
                double saz = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;
              
               //add SOT contribution
               double theta = atan2(m.y,m.x);
                  if(theta != theta) theta = 0.0;
               double phi = acos(m.z);
                  // if(st::internal::spin_acc_sign.at(cell) != 1.0 && st::internal::spin_acc_sign.at(cell) != -1.0) std::cout << st::internal::spin_acc_sign.at(cell) << ", " << stack << ", " << cell << std::endl;
               double say_sot = 0.0;// sin(phi)*(1.4+1.2*sin(theta)*sin(theta))*st::internal::sa_infinity[cell]*je*1e-11*st::internal::spin_acc_sign[cell]; //j_sd / mu_b dS/S 8 10^8 / 10^7 A/m^2 = T
               double sax_sot = 0.0;//sin(phi)*sin(2*theta)*st::internal::sa_infinity[cell]*je*1e-11*st::internal::spin_acc_sign[cell];
               double saz_sot = 0.0;// -cos(theta)*sin(phi)*0.12*st::internal::sa_infinity[cell]*je*1e-11*st::internal::spin_acc_sign[cell];

               // double sax_sot = sin(phi)*sin(2.0*theta)*0.5*1.48e7*je*1e-11; //j_sd / mu_b dS/S 8 10^8 / 10^7 A/m^2 = T
               // double say_sot = sin(phi)*(1.2*sin(theta)+1.4)*1.48e7*je*1e-11;
               // double saz_sot = -cos(theta)*sin(phi)*0.12*1.48e7*je*1e-11;
              //  std::cout << sax << ", " << say << ", " << saz << ", " \
                                          <<sax_sot << ", " << say_sot << ", " << saz_sot << std::endl;
               //--------------------------------------------
               // Step 4 calculate the spin current (jm_end)
               //--------------------------------------------
               const double ac_bd = a*c + b*d; 
               const double ad_bc = a*d - b*c;

               const double divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl;
               const double divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx);
               const double divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx);

               const double divsax = b1.x*divsa_para + b2.x*divsa_perp2 + b3.x*divsa_perp3;
               const double divsay = b1.y*divsa_para + b2.y*divsa_perp2 + b3.y*divsa_perp3;
               const double divsaz = b1.z*divsa_para + b2.z*divsa_perp2 + b3.z*divsa_perp3;

               const double dot = m.x*divsax + m.y*divsay + m.z*divsaz;

               const double pre_jmx = divsax - Bc*Bd*m.x*dot;
               const double pre_jmy = divsay - Bc*Bd*m.y*dot;
               const double pre_jmz = divsaz - Bc*Bd*m.z*dot;

                double jmx = Bc*je*m.x - twoDo*pre_jmx;
                double jmy = Bc*je*m.y - twoDo*pre_jmy;
                double jmz = Bc*je*m.z - twoDo*pre_jmz;

               if(st::internal::sot_sa) {
                   jmx =  - twoDo*pre_jmx;
                   jmy =  - twoDo*pre_jmy;
                   jmz =  - twoDo*pre_jmz;
               }
               if(st::internal::cell_natom[cell]>0) {
               //    if(st::internal::sot_sa && st::internal::sot_sa_source[cell]){
               // //       st::internal::sa[cellx] += st::internal::sa_infinity[cell]*0.0;
               // //       st::internal::sa[celly] += st::internal::sa_infinity[cell]*1.0;
               // //       st::internal::sa[cellz] += st::internal::sa_infinity[cell]*0.0;

               // //    // std::cout << st::internal::sa[3*idx+0] << ", " << st::internal::sa[3*idx+1] << ", " << st::internal::sa[3*idx+2] << std::endl;

               //       st::internal::j [cellx] += st::internal::initial_beta*je*0.0;
               //       st::internal::j [celly] += st::internal::initial_beta*je*1.0;
               //       st::internal::j [cellz] += st::internal::initial_beta*je*0.0;
               // }
                  // Save values for the spin accumulation
                  st::internal::sa[cellx] = sax;
                  st::internal::sa[celly] = say;
                  st::internal::sa[cellz] = saz;

                  // st::internal::sa_sot[cellx] += sax_sot;
                  // st::internal::sa_sot[celly] += say_sot;
                  // st::internal::sa_sot[cellz] += saz_sot;

                  // Save values for the spin current
                  st::internal::j[cellx] = jmx;
                  st::internal::j[celly] = jmy;
                  st::internal::j[cellz] = jmz;

                  // }  else {// Calculate spin torque energy for cell (Joules)
                  st::internal::spin_torque[cellx] += microcell_volume * st::internal::sd_exchange[cell] * (sax+sax_sot) * i_e * i_muB;
                  st::internal::spin_torque[celly] += microcell_volume * st::internal::sd_exchange[cell] * (say+say_sot) * i_e * i_muB;
                  st::internal::spin_torque[cellz] += microcell_volume * st::internal::sd_exchange[cell] * (saz+saz_sot) * i_e * i_muB; 
                  // }
               } 
               else{
                  // Save values for the spin accumulation
                  st::internal::sa[cellx] = st::internal::sa[pcellx];
                  st::internal::sa[celly] = st::internal::sa[pcelly];
                  st::internal::sa[cellz] = st::internal::sa[pcellz];

                  // st::internal::sa_sot[cellx] += st::internal::sa_sot[pcellx];
                  // st::internal::sa_sot[celly] += st::internal::sa_sot[pcelly];
                  // st::internal::sa_sot[cellz] += st::internal::sa_sot[pcellz];

                  // Save values for the spin current
                  st::internal::j[cellx] = st::internal::j[pcellx];
                  st::internal::j[celly] = st::internal::j[pcelly];
                  st::internal::j[cellz] = st::internal::j[pcellz];

                  // Calculate spin torque energy for cell (Joules)
                  st::internal::spin_torque[cellx] = st::internal::spin_torque[pcellx];
                  st::internal::spin_torque[celly] = st::internal::spin_torque[pcelly];
                  st::internal::spin_torque[cellz] = st::internal::spin_torque[pcellz];
                  // }
               }
               if(st::internal::output_torque_data) {
               //--------------------------------------------
               // Step 5 calculate the spin torque of each cell
               //--------------------------------------------

               //convert M of previous cell into basis b1, b2, b3

               M.xx = b1.x;    M.xy = b2.x;    M.xz = b3.x;
               M.yx = b1.y;    M.yy = b2.y;    M.yz = b3.y;
               M.zx = b1.z;    M.zy = b2.z;    M.zz = b3.z;

               V.x = pm.x;
               V.y = pm.y;
               V.z = pm.z;

               const st::internal::three_vector_t pm_basis = gaussian_elimination(M, V);

               //const double pm_b1 = pm_basis.x; // unused variable
               const double pm_b2 = pm_basis.y;
               const double pm_b3 = pm_basis.z;

               // Calculate the spin torque coefficients describing ast and nast
               const double prefac_sc = microcell_volume * st::internal::sd_exchange[cell] * i_e * i_muB;
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
                st::internal::coeff_ast[cell]  = aj;
                st::internal::coeff_nast[cell] = bj;

                SxSp[0]=(m.y*pm.z-m.z*pm.y);
                SxSp[1]=(m.z*pm.x-m.x*pm.z);
                SxSp[2]=(m.x*pm.y-m.y*pm.x);

                SxSxSp[0]= (m.y*SxSp[2]-m.z*SxSp[1]);
                SxSxSp[1]= (m.z*SxSp[0]-m.x*SxSp[2]);
                SxSxSp[2]= (m.x*SxSp[1]-m.y*SxSp[0]);

                //calculate directly from J(Sxm)
                st::internal::total_ST[cellx] = prefac_sc*(m.y*saz-m.z*say);
                st::internal::total_ST[celly] = prefac_sc*(m.z*sax-m.x*saz);
                st::internal::total_ST[cellz] = prefac_sc*(m.x*say-m.y*sax);

                st::internal::ast[cellx] = -aj*SxSxSp[0];
                st::internal::ast[celly] = -aj*SxSxSp[1];
                st::internal::ast[cellz] = -aj*SxSxSp[2];

                st::internal::nast[cellx] = bj*SxSp[0];
                st::internal::nast[celly] = bj*SxSp[1];
                st::internal::nast[cellz] = bj*SxSp[2];
               }
            } // end of cell loop


            //reverse process
         } // end of stack loop

         // Reduce all microcell spin torques on all nodes
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::spin_torque[0],st::internal::spin_torque.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            // if(st::internal::output_torque_data) MPI_Allreduce(MPI_IN_PLACE, &st::internal::total_ST[0],st::internal::total_ST.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif
         st::internal::output_microcell_data();

         st::spin_acc_time += stopwatch.elapsed_seconds();

         return;
      }

       void calculate_sot_accumulation(){

         stopwatch_t stopwatch;
         stopwatch.start();

         // Zero all shared arrays (essential for parallelisation to work)
         std::fill (st::internal::spin_torque.begin(),st::internal::spin_torque.end(),0.0);
         std::fill (st::internal::sa.begin(),st::internal::sa.end(),0.0);
         std::fill (st::internal::sa_sot_up.begin(),st::internal::sa_sot_up.end(),0.0);
         std::fill (st::internal::sa_sot_down.begin(),st::internal::sa_sot_down.end(),0.0);
          std::fill (st::internal::j.begin(),st::internal::j.end(),0.0);
         std::fill (st::internal::j_up.begin(),st::internal::j_up.end(),0.0);
         std::fill (st::internal::j_down.begin(),st::internal::j_down.end(),0.0);
         // Declare resuable temporary variables
         st::internal::matrix_t itm; // inverse transformation matrix
         st::internal::matrix_t M; // general matrix
         st::internal::three_vector_t V(0.0,0.0,0.0); // general 3-vector

         // reference basis vectors
         const st::internal::three_vector_t bx(1.0,0.0,0.0);
         const st::internal::three_vector_t by(0.0,1.0,0.0);
         const st::internal::three_vector_t bz(0.0,0.0,1.0);

         // local basis vectors
         st::internal::three_vector_t b1(1.0,0.0,0.0);
         st::internal::three_vector_t b2(0.0,1.0,0.0);
         st::internal::three_vector_t b3(0.0,0.0,1.0);

         // set local constants
         double je = program::fractional_electric_field_strength* st::internal::je; // current (C/s)

         //---------------------------------------------------------------------------------------------------
         //set parameters for TMR calculation
         if(st::internal::TMRenable == true){
            std::cout << "wrong" << std::endl;
            int FL =	st::internal::free_layer;
            int RL =	st::internal::reference_layer;
            double dot = st::internal::magx_mat[RL]*st::internal::magx_mat[FL]+
                         st::internal::magy_mat[RL]*st::internal::magy_mat[FL]+
                         st::internal::magz_mat[RL]*st::internal::magz_mat[FL];

            // RE - this code is not general! Needs to be fixed. Placeholder added in the meantime
            // AM (2020) - Code fixed, but still codes refers to MTJ RL/barrier/FL specifically
            double MgO_thickness = (create::get_material_height_min(FL)-create::get_material_height_max(RL))*cs::system_dimensions[2]*1.0e-10;

            //calculate the relative angle of two FMs
            st::internal::rel_angle = acos(dot);
            double plus_cos = 1.0+cos(st::internal::rel_angle);
            double minus_cos = 1.0-cos(st::internal::rel_angle);
            double exp_t = exp(-MgO_thickness/0.25e-9);

            double jtunnel = st::internal::je*0.5*(plus_cos+0.5*minus_cos)*exp_t;
//            std::cout << "t_MgO=( " << create::get_material_height_min(FL) << " - " << create::get_material_height_max(RL) << " ) = " << MgO_thickness << "\tJe\t" << st::internal::je << "\tJe_tun\t" << jtunnel << std::endl;

            //set the current je and spin poralisation parameters
            je = jtunnel;
            // AM (2020) - I think the default parameters should be rescaled by same factor as tunnelling current and not changed using those of material 0 arbitrarily
            st::internal::default_properties.beta_cond *= /*st::internal::mp[0].beta_cond**/0.5*(plus_cos+0.5*minus_cos)*exp_t;
            st::internal::default_properties.beta_diff *= /*st::internal::mp[0].beta_diff**/0.5*(plus_cos+0.5*minus_cos)*exp_t;

            // Calculate spin torque parameters
            for(size_t cell=0; cell<st::internal::beta_cond.size(); ++cell){

               // check for zero atoms in cell
               if(st::internal::cell_natom[cell] <= 0.0001){
                  st::internal::beta_cond[cell]   = st::internal::default_properties.beta_cond;
                  st::internal::beta_diff[cell]   = st::internal::default_properties.beta_diff;

                  const double hbar = 1.05457162e-34;
                  const double B  = st::internal::beta_cond[cell];
                  const double Bp = st::internal::beta_diff[cell];
                  const double lambda_sdl = st::internal::lambda_sdl[cell];
                  const double Do = st::internal::diffusion[cell];
                  const double Jsd = st::internal::sd_exchange[cell];

                  const double BBp = 1.0/sqrt(1.0-B*Bp);
                  const double lambda_sf = lambda_sdl*BBp;
                  const double lambda_j = sqrt(2.0*hbar*Do/Jsd); // Angstroms
                  const double lambda_sf2 = lambda_sf*lambda_sf;
                  const double lambda_j2 = lambda_j*lambda_j;

                  std::complex<double> inside (1.0/lambda_sf2, -1.0/lambda_j2);
                  std::complex<double> inv_lplus = sqrt(inside);

                  st::internal::a[cell] =  real(inv_lplus);
                  st::internal::b[cell] = -imag(inv_lplus);
               }
            }
            st::internal::output_base_microcell_data();
         }

         //---------------------------------------------------------------------------------------------------

         const double i_muB = 1.0/9.274e-24; // J/T
         const double i_e = 1.0/1.60217662e-19; // electronic charge (Coulombs)
         const double microcell_volume = (st::internal::micro_cell_size[st::internal::stx] *
                                          st::internal::micro_cell_size[st::internal::sty] *
                                          st::internal::micro_cell_thickness)*1.e-30; // m^3
         const double atomcell_volume = 15.7624e-30;
         // loop over all 1D stacks (in parallel)
         int int_stacks = st::internal::mpi_stack_list.size();
         // std::vector<int> stacks_list;
         // stacks_list = st::internal::mpi_stack_list           
        

         int stack = 0;
         const double sot_sd_exchange_multiple = 1;// 3e2;
         // std::cout << int_stacks << ", " << st::internal::mpi_stack_list.size() << std::endl;
         for(int s=0; s < int_stacks; ++s) {
            stack = st::internal::mpi_stack_list.at(s);
        
            // if(stack % 3 == 0) continue;
               // std::cout << stack << std::endl;
               // std::cout << st::internal::spin_acc_sign[stack_index[stack]] << ", " << stack << std::endl;
      
            // std::cout << vmpi::my_rank << ", " << stack << std::endl;
            // determine starting cell in stack
            const int idx = stack_index[stack] + 1;
           // std::cout << stack << ", " << idx << std::endl;
            // set initial values
         //   if(st::internal::fbc) { 
         //       st::internal::sa[3*idx+0] = 0.0;
         //       st::internal::sa[3*idx+1] = 0.0;
         //       st::internal::sa[3*idx+2] = 0.0;
         //       st::internal::j [3*idx+0] = st::internal::initial_beta*je*st::internal::initial_m[0];
         //       st::internal::j [3*idx+1] = st::internal::initial_beta*je*st::internal::initial_m[1];
         //       st::internal::j [3*idx+2] = st::internal::initial_beta*je*st::internal::initial_m[2];
         //   }  else if (!st::internal::fbc) {
         //       const double mod = 1.0/sqrt(st::internal::m [3*idx+0]*st::internal::m [3*idx+0] + st::internal::m [3*idx+1]*st::internal::m [3*idx+1] + st::internal::m [3*idx+2]*st::internal::m[3*idx+2]);
         //       st::internal::sa[3*idx+0] = st::internal::sa_infinity[idx]*st::internal::m [3*idx+0]*mod;
         //       st::internal::sa[3*idx+1] = st::internal::sa_infinity[idx]*st::internal::m [3*idx+1]*mod;
         //       st::internal::sa[3*idx+2] = st::internal::sa_infinity[idx]*st::internal::m [3*idx+2]*mod;
         //
         //       st::internal::j [3*idx+0] = st::internal::initial_beta*je*st::internal::m [3*idx+0];
         //       st::internal::j [3*idx+1] = st::internal::initial_beta*je*st::internal::m [3*idx+1];
         //       st::internal::j [3*idx+2] = st::internal::initial_beta*je*st::internal::m [3*idx+2];
         //   }
         //
         //   if(st::internal::sot_sa && st::internal::sot_sa_source[idx]){
         //       st::internal::sa[3*idx+0] = st::internal::sa_infinity[idx]*0.0;
         //       st::internal::sa[3*idx+1] = program::fractional_electric_field_strength*st::internal::sa_infinity[idx];
         //       st::internal::sa[3*idx+2] = st::internal::sa_infinity[idx]*0.0;
         //
         //       // std::cout << st::internal::sa[3*idx+0] << ", " << st::internal::sa[3*idx+1] << ", " << st::internal::sa[3*idx+2] << std::endl;
         //
         //       st::internal::j [3*idx+0] = st::internal::initial_beta*je*0.0;
         //       st::internal::j [3*idx+1] = 0.0;//st::internal::initial_beta*je*1.0;
         //       st::internal::j [3*idx+2] = st::internal::initial_beta*je*0.0;
         //   }
         //
          //  if(sim::time%(ST_output_rate) ==0) std::cout<< stack << "\t" << st::internal::default_properties.sa_infinity << "\t" << st::internal::init_stack_mag[((stack)%6)*3 + 0]/sqrt(2.0) << "\t" << st::internal::init_stack_mag[((stack)%6)*3 + 1]/sqrt(2.0) << "\t" << st::internal::m[3*idx+0] << " \t" <<  st::internal::m[3*idx+1] << "\t" << st::internal::sa[3*idx+0] << "\t" << st::internal::sa[3*idx+1] << std::endl;

            // Mn1 process
            for(int cell=idx+1; cell<idx+num_microcells_per_stack; ++cell) {

               if(st::internal::spin_acc_sign[cell] != 1) continue;
               // calculate cell id's
               const int cellx = 3*cell+0;
               const int celly = 3*cell+1;
               const int cellz = 3*cell+2;

               // calculate previous cell id's
               const int pcellx = 3*(cell-1)+0;
               const int pcelly = 3*(cell-1)+1;
               const int pcellz = 3*(cell-1)+2;

               st::internal::three_vector_t  m(0.0,  0.0,  0.0);
               // copy array values to temporaries for readability
               // if(st::internal::sot_sa && st::internal::sot_sa_source[cell]) {
               //    m.x = 0;
               //    m.y = -1.0;
               //    m.z = 0.0; // Au cell artificial magnetisations
               // } else {
                  m.x =st::internal::m[cellx];
                  m.y = st::internal::m[celly];
                  m.z = st::internal::m[cellz]; // current cell magnetisations
               // }
               st::internal::three_vector_t pm(0.0,0.0,0.0);
               // if(st::internal::sot_sa && st::internal::sot_sa_source[cell-1]) {
                  pm.x = 0;
                  pm.y = (m.y >= 0.0) ?  1.0:-1.0;
                  if(m.x == 1.0) pm.y = -1.0;
                  else if (m.x == -1.0)  pm.y = 1.0;
               // } else {
               //    pm.x =st::internal::m[pcellx];
               //    pm.y = st::internal::m[pcelly];
               //    pm.z = st::internal::m[pcellz]; // current cell magnetisations
               // }

               const double modm = sqrt(m.x*m.x + m.y*m.y + m.z*m.z);
               // const double pmodm = sqrt(pm.x*pm.x + pm.y*pm.y + pm.z*pm.z);

               // Check for zero magnetization in normalization
               if(modm > 1.e-11){
                  m.x = m.x/modm;
                  m.y = m.y/modm;
                  m.z = m.z/modm;
               }
               else{
                  m.x = 0.0;
                  m.y = 0.0;
                  m.z = 0.0;
               }
               // if(pmodm > 1.e-11){
               //    pm.x = pm.x/pmodm;
               //    pm.y = pm.y/pmodm;
               //    pm.z = pm.z/pmodm;
               // }
               // else{
               //    pm.x = 0.0;
               //    pm.y = 0.0;
               //    pm.z = 0.0;
               // }

               //---------------------------------------------------------------------
               // Step 1 calculate coordinate transformation m -> m' for current cell
               //---------------------------------------------------------------------

               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m, itm);

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
               const double Bc = st::internal::sot_beta_cond; //sot beta
               const double Bd = st::internal::sot_beta_diff; // beta_prime
               const double Do = st::internal::sot_diffusion;
               // st::internal::three_vector_t jm0(st::internal::initial_theta*je/sqrt(2.0),-st::internal::initial_theta*je/sqrt(2.0),0.0);
               st::internal::three_vector_t jm0(0.0,pm.y*st::internal::initial_theta*je,0.0);
               // if(st::internal::sot_sa_source[cell-1]) {
                  // jm0.y = st::internal::initial_theta*je;
               // }

               //  Calculate gradient dsacc/dx
               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m.x*m.x - twoDo;     M.xy = preD*m.y*m.x;             M.xz = preD*m.z*m.x;
               M.yx = preD*m.x*m.y;             M.yy = preD*m.y*m.y - twoDo;     M.yz = preD*m.z*m.y;
               M.zx = preD*m.x*m.z;             M.zy = preD*m.y*m.z;             M.zz = preD*m.z*m.z - twoDo;

               // if(st::internal::sot_sa ) {
                  V.x = jm0.x;// - Bc*je*m.x;
                  V.y = jm0.y;// + Bc*je*m.y;
                  V.z = jm0.z;// - Bc*je*m.z;
               // } else {
               //    V.x = jm0.x - Bc*je*m.x;
               //    V.y = jm0.y - Bc*je*m.y;
               //    V.z = jm0.z - Bc*je*m.z;
               // }
               const st::internal::three_vector_t divm_0 = gaussian_elimination(M, V);

               // Calculate mp(0), c and d
               double i_lsdl = 1.0/st::internal::sot_lambda_sdl;
               double mp_inf = st::internal::sot_sa_infinity*program::fractional_electric_field_strength;
               //   if(st::internal::sot_sa && st::internal::sot_sa_source[cell]) mp_inf *= program::fractional_electric_field_strength;
               double a = st::internal::sot_a;
               double b = st::internal::sot_b;
               double two_a = 2.0*a;
               double two_b = 2.0*b;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               st::internal::three_vector_t coeff = gaussian_elimination(M, V);

               double mp_0 = coeff.x;
               double c    = coeff.y;
               double d    = coeff.z;

               //------------------------------------
               // Step 3 calculate spin accumulation
               //------------------------------------

               const double x = st::internal::micro_cell_thickness*1.0e-10; // Convert to metres
               double cos_bx = cos(b*x);
               double sin_bx = sin(b*x);
               double e_xsdl = exp(-x*i_lsdl);
               double e_ax   = exp(-a*x);
               double prefac = (2.0*e_ax);

               double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);

               //convert mp and m_perp
               double sax_sot = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
               double say_sot = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
               double saz_sot = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;
                             
               //--------------------------------------------
               // Step 4 calculate the spin current with cell constants (jm_end)
               //--------------------------------------------

               i_lsdl = 1.0/st::internal::lambda_sdl[cell];
               mp_inf = st::internal::sa_infinity[cell];
               a = st::internal::a[cell];
               b = st::internal::b[cell];
               two_a = 2.0*a;
               two_b = 2.0*b;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               coeff = gaussian_elimination(M, V);

               mp_0 = coeff.x;
               c    = coeff.y;
               d    = coeff.z;

               cos_bx = cos(b*x);
               sin_bx = sin(b*x);
               e_xsdl = exp(-x*i_lsdl);
               e_ax   = exp(-a*x);
               prefac = (2.0*e_ax);

               sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               sa_perp3 = prefac*(c*sin_bx + d*cos_bx);

               const double ac_bd = a*c + b*d; 
               const double ad_bc = a*d - b*c;

               const double divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl;
               const double divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx);
               const double divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx);

               const double divsax = b1.x*divsa_para + b2.x*divsa_perp2 + b3.x*divsa_perp3;
               const double divsay = b1.y*divsa_para + b2.y*divsa_perp2 + b3.y*divsa_perp3;
               const double divsaz = b1.z*divsa_para + b2.z*divsa_perp2 + b3.z*divsa_perp3;

               const double dot = m.x*divsax + m.y*divsay + m.z*divsaz;

               const double pre_jmx = divsax - Bc*Bd*m.x*dot;
               const double pre_jmy = divsay - Bc*Bd*m.y*dot;
               const double pre_jmz = divsaz - Bc*Bd*m.z*dot;

               //  double jmx = Bc*je*m.x - twoDo*pre_jmx;
               //  double jmy = Bc*je*m.y - twoDo*pre_jmy;
               //  double jmz = Bc*je*m.z - twoDo*pre_jmz;

               // if(st::internal::sot_sa) {
                  double jmx =  - twoDo*pre_jmx;
                  double jmy =  - twoDo*pre_jmy;
                  double jmz =  - twoDo*pre_jmz;
               // }
               if(st::internal::cell_natom[cell]>0) {
               //    if(st::internal::sot_sa && st::internal::sot_sa_source[cell]){
               // //       st::internal::sa[cellx] += st::internal::sa_infinity[cell]*0.0;
               // //       st::internal::sa[celly] += st::internal::sa_infinity[cell]*1.0;
               // //       st::internal::sa[cellz] += st::internal::sa_infinity[cell]*0.0;
               // //    // std::cout << st::internal::sa[3*idx+0] << ", " << st::internal::sa[3*idx+1] << ", " << st::internal::sa[3*idx+2] << std::endl;
               //       st::internal::j [cellx] += st::internal::initial_beta*je*0.0;
               //       st::internal::j [celly] += st::internal::initial_beta*je*1.0;
               //       st::internal::j [cellz] += st::internal::initial_beta*je*0.0;
               // }
                  // Save values for the spin accumulation
                  // st::internal::sa[cellx] = sax;//(sax-mp_inf*m.x)/(mp_inf*m.x);
                  // st::internal::sa[celly] = say;//(say-mp_inf*m.y)/(mp_inf*m.y);
                  // st::internal::sa[cellz] = saz;//(saz-mp_inf*m.z)/(mp_inf*m.z);

                  st::internal::sa_sot_up[cellx] = sax_sot;
                  st::internal::sa_sot_up[celly] = say_sot;
                  st::internal::sa_sot_up[cellz] = saz_sot;

                  // Save values for the spin current
                  st::internal::j_up[cellx] = jmx;
                  st::internal::j_up[celly] = jmy;
                  st::internal::j_up[cellz] = jmz;

                  // }  else {// Calculate spin torque energy for cell (Joules)
                  st::internal::spin_torque[cellx] = sot_sd_exchange_multiple*atomcell_volume * st::internal::sot_sd_exchange * (sax_sot) * i_e * i_muB;
                  st::internal::spin_torque[celly] = sot_sd_exchange_multiple*atomcell_volume * st::internal::sot_sd_exchange * (say_sot) * i_e * i_muB;
                  st::internal::spin_torque[cellz] = sot_sd_exchange_multiple*atomcell_volume * st::internal::sot_sd_exchange * (saz_sot) * i_e * i_muB; 
                  
                  // st::internal::spin_torque[cellx] = 0.0;// microcell_volume * st::internal::sd_exchange[cell] * (sax+sax_sot) * i_e * i_muB;
                  // st::internal::spin_torque[celly] = 0.0;//microcell_volume * st::internal::sd_exchange[cell] * (say+say_sot) * i_e * i_muB;
                  // st::internal::spin_torque[cellz] = 0.0;//microcell_volume * st::internal::sd_exchange[cell] * (saz+saz_sot) * i_e * i_muB; 
                  
               } 
               else{
                  // Save values for the spin accumulation
                  st::internal::sa_sot_up[cellx] = st::internal::sa_sot_up[pcellx];
                  st::internal::sa_sot_up[celly] = st::internal::sa_sot_up[pcelly];
                  st::internal::sa_sot_up[cellz] = st::internal::sa_sot_up[pcellz];

                  // st::internal::sa_sot_up[cellx] += st::internal::sa_sot[pcellx];
                  // st::internal::sa_sot_up[celly] += st::internal::sa_sot[pcelly];
                  // st::internal::sa_sot_up[cellz] += st::internal::sa_sot[pcellz];

                  // Save values for the spin current
                  st::internal::j_up[cellx] = st::internal::j_up[pcellx];
                  st::internal::j_up[celly] = st::internal::j_up[pcelly];
                  st::internal::j_up[cellz] = st::internal::j_up[pcellz];

                  // Calculate spin torque energy for cell (Joules)
                  st::internal::spin_torque[cellx] = st::internal::spin_torque[pcellx];
                  st::internal::spin_torque[celly] = st::internal::spin_torque[pcelly];
                  st::internal::spin_torque[cellz] = st::internal::spin_torque[pcellz];
                  // }
               }
               if(st::internal::output_torque_data) {
               //--------------------------------------------
               // Step 5 calculate the spin torque of each cell
               //--------------------------------------------

               //convert M of previous cell into basis b1, b2, b3

               M.xx = b1.x;    M.xy = b2.x;    M.xz = b3.x;
               M.yx = b1.y;    M.yy = b2.y;    M.yz = b3.y;
               M.zx = b1.z;    M.zy = b2.z;    M.zz = b3.z;

               V.x = pm.x;
               V.y = pm.y;
               V.z = pm.z;

               const st::internal::three_vector_t pm_basis = gaussian_elimination(M, V);

               //const double pm_b1 = pm_basis.x; // unused variable
               const double pm_b2 = pm_basis.y;
               const double pm_b3 = pm_basis.z;

               // Calculate the spin torque coefficients describing ast and nast
               const double prefac_sc = microcell_volume * st::internal::sot_sd_exchange * i_e * i_muB;
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
                st::internal::coeff_ast[cell]  = aj;
                st::internal::coeff_nast[cell] = bj;

                SxSp[0]=(m.y*pm.z-m.z*pm.y);
                SxSp[1]=(m.z*pm.x-m.x*pm.z);
                SxSp[2]=(m.x*pm.y-m.y*pm.x);

                SxSxSp[0]= (m.y*SxSp[2]-m.z*SxSp[1]);
                SxSxSp[1]= (m.z*SxSp[0]-m.x*SxSp[2]);
                SxSxSp[2]= (m.x*SxSp[1]-m.y*SxSp[0]);

                //calculate directly from J(Sxm)
                st::internal::total_ST[cellx] = prefac_sc*(m.y*saz_sot-m.z*say_sot);
                st::internal::total_ST[celly] = prefac_sc*(m.z*sax_sot-m.x*saz_sot);
                st::internal::total_ST[cellz] = prefac_sc*(m.x*say_sot-m.y*sax_sot);

                st::internal::ast[cellx] = -aj*SxSxSp[0];
                st::internal::ast[celly] = -aj*SxSxSp[1];
                st::internal::ast[cellz] = -aj*SxSxSp[2];

                st::internal::nast[cellx] = bj*SxSp[0];
                st::internal::nast[celly] = bj*SxSp[1];
                st::internal::nast[cellz] = bj*SxSp[2];
               }
            } // end of cell loop
            
            
            //Mn2 process    
            for(int cell=idx+num_microcells_per_stack-2; cell > idx; --cell) {
               
               if(st::internal::spin_acc_sign[cell] != -1) continue;
               // calculate cell id's
               const int cellx = 3*cell+0;
               const int celly = 3*cell+1;
               const int cellz = 3*cell+2;

               // calculate previous cell id's
               const int pcellx = 3*(cell+1)+0;
               const int pcelly = 3*(cell+1)+1;
               const int pcellz = 3*(cell+1)+2;

               // copy array values to temporaries for readability
               st::internal::three_vector_t  m(0.0,  0.0,  0.0);
               // copy array values to temporaries for readability
               // if(st::internal::sot_sa && st::internal::sot_sa_source[cell]) {
               //    m.x = 0;
               //    m.y = 1.0;
               //    m.z = 0.0; // Au cell artificial magnetisations
               // } else {
                  m.x = st::internal::m[cellx];
                  m.y = st::internal::m[celly];
                  m.z = st::internal::m[cellz]; // current cell magnetisations
               // }
               st::internal::three_vector_t pm(0.0,0.0,0.0);
               // if(st::internal::sot_sa && st::internal::sot_sa_source[cell+1]) {
                  pm.x = 0;
                  pm.y = (m.y >= 0.0) ?  1.0:-1.0;
                  if(m.x == 1.0) pm.y = 1.0;
                  else if (m.x == -1.0)  pm.y = -1.0;
                  pm.z = 0.0; // Au cell artificial magnetisations
               // } else {
               //    pm.x =st::internal::m[pcellx];
               //    pm.y = st::internal::m[pcelly];
               //    pm.z = st::internal::m[pcellz]; // current cell magnetisations
               // }

               const double modm = sqrt(m.x*m.x + m.y*m.y + m.z*m.z);
               // const double pmodm = sqrt(pm.x*pm.x + pm.y*pm.y + pm.z*pm.z);

               // Check for zero magnetization in normalization
               if(modm > 1.e-11){
                  m.x = m.x/modm;
                  m.y = m.y/modm;
                  m.z = m.z/modm;
               }
               else{
                  m.x = 0.0;
                  m.y = 0.0;
                  m.z = 0.0;
               }
               // if(pmodm > 1.e-11){
               //    pm.x = pm.x/pmodm;
               //    pm.y = pm.y/pmodm;
               //    pm.z = pm.z/pmodm;
               // }
               // else{
               //    pm.x = 0.0;
               //    pm.y = 0.0;
               //    pm.z = 0.0;
               // }

               //---------------------------------------------------------------------
               // Step 1 calculate coordinate transformation m -> m' for current cell
               //---------------------------------------------------------------------

               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m, itm);

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
               const double Bc = st::internal::sot_beta_cond; //sot beta
               const double Bd = st::internal::sot_beta_diff; // beta_prime
               const double Do = st::internal::sot_diffusion;
               //  st::internal::three_vector_t jm0(-st::internal::initial_theta*je/sqrt(2.0),st::internal::initial_theta*je/sqrt(2.0),0.0);
               st::internal::three_vector_t jm0(0.0,pm.y*st::internal::initial_theta*je,0.0);
               // if(st::internal::sot_sa_source[cell-1]) {
                  // jm0.y = -st::internal::initial_theta*je;
               // }

               //  Calculate gradient dsacc/dx
               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m.x*m.x - twoDo;     M.xy = preD*m.y*m.x;             M.xz = preD*m.z*m.x;
               M.yx = preD*m.x*m.y;             M.yy = preD*m.y*m.y - twoDo;     M.yz = preD*m.z*m.y;
               M.zx = preD*m.x*m.z;             M.zy = preD*m.y*m.z;             M.zz = preD*m.z*m.z - twoDo;

               // if(st::internal::sot_sa) {
                  V.x = jm0.x;// - Bc*je*m.x;
                  V.y = jm0.y;// - Bc*je*m.y;
                  V.z = jm0.z;// - Bc*je*m.z;
               // } else {
               //    V.x = jm0.x - Bc*je*m.x;
               //    V.y = jm0.y - Bc*je*m.y;
               //    V.z = jm0.z - Bc*je*m.z;
               // }

               const st::internal::three_vector_t divm_0 = gaussian_elimination(M, V);

               // Calculate mp(0), c and d
                double i_lsdl = 1.0/st::internal::sot_lambda_sdl;
               double mp_inf = st::internal::sot_sa_infinity*program::fractional_electric_field_strength;
               //   if(st::internal::sot_sa && st::internal::sot_sa_source[cell]) mp_inf *= program::fractional_electric_field_strength;
                double a = st::internal::sot_a;
                double b = st::internal::sot_b;
                double two_a = 2.0*a;
                double two_b = 2.0*b;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               st::internal::three_vector_t coeff = gaussian_elimination(M, V);

                double mp_0 = coeff.x;
                double c    = coeff.y;
                double d    = coeff.z;

               //------------------------------------
               // Step 3 calculate spin accumulation
               //------------------------------------

               const double x = st::internal::micro_cell_thickness*1.0e-10; // Convert to metres
               double cos_bx = cos(b*x);
                double sin_bx = sin(b*x);
                double e_xsdl = exp(-x*i_lsdl);
                double e_ax   = exp(-a*x);
                double prefac = (2.0*e_ax);

                double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
                double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
                double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);

               //convert mp and m_perp
                double sax_sot = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
                double say_sot = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
                double saz_sot = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;
              
               //add SOT contribution
               // double theta = atan2(m.y,m.x);
               //    if(theta != theta) theta = 0.0;
               // double phi = acos(m.z);
                  // if(st::internal::spin_acc_sign.at(cell) != 1.0 && st::internal::spin_acc_sign.at(cell) != -1.0) std::cout << st::internal::spin_acc_sign.at(cell) << ", " << stack << ", " << cell << std::endl;
               // double say_sot = 0.0;// sin(phi)*(1.4+1.2*sin(theta)*sin(theta))*st::internal::sa_infinity[cell]*je*1e-11*st::internal::spin_acc_sign[cell]; //j_sd / mu_b dS/S 8 10^8 / 10^7 A/m^2 = T
               // double sax_sot = 0.0;//sin(phi)*sin(2*theta)*st::internal::sa_infinity[cell]*je*1e-11*st::internal::spin_acc_sign[cell];
               // double saz_sot = 0.0;// -cos(theta)*sin(phi)*0.12*st::internal::sa_infinity[cell]*je*1e-11*st::internal::spin_acc_sign[cell];

               // double sax_sot = sin(phi)*sin(2.0*theta)*0.5*1.48e7*je*1e-11; //j_sd / mu_b dS/S 8 10^8 / 10^7 A/m^2 = T
               // double say_sot = sin(phi)*(1.2*sin(theta)+1.4)*1.48e7*je*1e-11;
               // double saz_sot = -cos(theta)*sin(phi)*0.12*1.48e7*je*1e-11;
              //  std::cout << sax << ", " << say << ", " << saz << ", " \
                                          <<sax_sot << ", " << say_sot << ", " << saz_sot << std::endl;
            
               //--------------------------------------------
               // Step 4 calculate the spin current with cell constants (jm_end)
               //--------------------------------------------

               i_lsdl = 1.0/st::internal::lambda_sdl[cell];
               mp_inf = st::internal::sa_infinity[cell];
               a = st::internal::a[cell];
               b = st::internal::b[cell];
               two_a = 2.0*a;
               two_b = 2.0*b;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               coeff = gaussian_elimination(M, V);

               mp_0 = coeff.x;
               c    = coeff.y;
               d    = coeff.z;

               cos_bx = cos(b*x);
               sin_bx = sin(b*x);
               e_xsdl = exp(-x*i_lsdl);
               e_ax   = exp(-a*x);
               prefac = (2.0*e_ax);

               sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               sa_perp3 = prefac*(c*sin_bx + d*cos_bx);

               const double ac_bd = a*c + b*d; 
               const double ad_bc = a*d - b*c;

               const double divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl;
               const double divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx);
               const double divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx);

               const double divsax = b1.x*divsa_para + b2.x*divsa_perp2 + b3.x*divsa_perp3;
               const double divsay = b1.y*divsa_para + b2.y*divsa_perp2 + b3.y*divsa_perp3;
               const double divsaz = b1.z*divsa_para + b2.z*divsa_perp2 + b3.z*divsa_perp3;

               const double dot = m.x*divsax + m.y*divsay + m.z*divsaz;

               const double pre_jmx = divsax - Bc*Bd*m.x*dot;
               const double pre_jmy = divsay - Bc*Bd*m.y*dot;
               const double pre_jmz = divsaz - Bc*Bd*m.z*dot;

               //  double jmx = Bc*je*m.x - twoDo*pre_jmx;
               //  double jmy = Bc*je*m.y - twoDo*pre_jmy;
               //  double jmz = Bc*je*m.z - twoDo*pre_jmz;

               // if(st::internal::sot_sa) {
               double jmx =  - twoDo*pre_jmx;
               double jmy =  - twoDo*pre_jmy;
               double jmz =  - twoDo*pre_jmz;

               if(st::internal::cell_natom[cell]>0) {

               //    if(st::internal::sot_sa && st::internal::sot_sa_source[cell]){
             //        
               // //       st::internal::sa[3*cell+0] += st::internal::sa_infinity[cell]*0.0;
               // //       st::internal::sa[3*cell+1] += -st::internal::sa_infinity[cell]*1.0;
               // //       st::internal::sa[3*cell+2] += st::internal::sa_infinity[cell]*0.0;
//
               // //       // std::cout << st::internal::sa[3*idx+0] << ", " << st::internal::sa[3*idx+1] << ", " << st::internal::sa[3*idx+2] << std::endl;
//
               //       st::internal::j [3*cell+0] += st::internal::initial_beta*je*0.0;
               //       st::internal::j [3*cell+1] += -st::internal::initial_beta*je*1.0;
               //       st::internal::j [3*cell+2] += st::internal::initial_beta*je*0.0;
               //    }
                  // Save values for the spin accumulation
                  // st::internal::sa[cellx] = sax;//(sax-mp_inf*m.x)/(mp_inf);
                  // st::internal::sa[celly] = say;//(say-mp_inf*m.y)/(mp_inf);
                  // st::internal::sa[cellz] = saz;//(saz-mp_inf*m.z)/(mp_inf);

                  st::internal::sa_sot_down[cellx] = sax_sot;
                  st::internal::sa_sot_down[celly] = say_sot;
                  st::internal::sa_sot_down[cellz] = saz_sot;

                  // Save values for the spin current
                  st::internal::j_down[cellx] = jmx;
                  st::internal::j_down[celly] = jmy;
                  st::internal::j_down[cellz] = jmz;

                  // }  else {// Calculate spin torque energy for cell (Joules)
                  st::internal::spin_torque[cellx] = sot_sd_exchange_multiple*atomcell_volume * st::internal::sot_sd_exchange * (sax_sot) * i_e * i_muB;
                  st::internal::spin_torque[celly] = sot_sd_exchange_multiple*atomcell_volume * st::internal::sot_sd_exchange * (say_sot) * i_e * i_muB;
                  st::internal::spin_torque[cellz] = sot_sd_exchange_multiple*atomcell_volume * st::internal::sot_sd_exchange * (saz_sot) * i_e * i_muB; 
                  
                  // st::internal::spin_torque[cellx] = 0.0;// microcell_volume * st::internal::sd_exchange[cell] * (sax+sax_sot) * i_e * i_muB;
                  // st::internal::spin_torque[celly] = 0.0;//microcell_volume * st::internal::sd_exchange[cell] * (say+say_sot) * i_e * i_muB;
                  // st::internal::spin_torque[cellz] = 0.0;//microcell_volume * st::internal::sd_exchange[cell] * (saz+saz_sot) * i_e * i_muB; 
                                    
               } 
               else{
                  // Save values for the spin accumulation
                  // st::internal::sa[cellx] += st::internal::sa[pcellx];
                  // st::internal::sa[celly] += st::internal::sa[pcelly];
                  // st::internal::sa[cellz] += st::internal::sa[pcellz];

                  st::internal::sa_sot_down[cellx] = st::internal::sa_sot_down[pcellx];
                  st::internal::sa_sot_down[celly] = st::internal::sa_sot_down[pcelly];
                  st::internal::sa_sot_down[cellz] = st::internal::sa_sot_down[pcellz];

                  // Save values for the spin current
                  st::internal::j_down[cellx] = st::internal::j_down[pcellx];
                  st::internal::j_down[celly] = st::internal::j_down[pcelly];
                  st::internal::j_down[cellz] = st::internal::j_down[pcellz];

                  // Calculate spin torque energy for cell (Joules)
                  st::internal::spin_torque[cellx] = st::internal::spin_torque[pcellx];
                  st::internal::spin_torque[celly] = st::internal::spin_torque[pcelly];
                  st::internal::spin_torque[cellz] = st::internal::spin_torque[pcellz];
                  // }
               }
               if(st::internal::output_torque_data) {
               //--------------------------------------------
               // Step 5 calculate the spin torque of each cell
               //--------------------------------------------

               //convert M of previous cell into basis b1, b2, b3

               M.xx = b1.x;    M.xy = b2.x;    M.xz = b3.x;
               M.yx = b1.y;    M.yy = b2.y;    M.yz = b3.y;
               M.zx = b1.z;    M.zy = b2.z;    M.zz = b3.z;

               V.x = pm.x;
               V.y = pm.y;
               V.z = pm.z;

               const st::internal::three_vector_t pm_basis = gaussian_elimination(M, V);

               //const double pm_b1 = pm_basis.x; // unused variable
               const double pm_b2 = pm_basis.y;
               const double pm_b3 = pm_basis.z;

               // Calculate the spin torque coefficients describing ast and nast
               const double prefac_sc = microcell_volume * st::internal::sot_sd_exchange * i_e * i_muB;
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
                st::internal::coeff_ast[cell]  = aj;
                st::internal::coeff_nast[cell] = bj;

                SxSp[0]=(m.y*pm.z-m.z*pm.y);
                SxSp[1]=(m.z*pm.x-m.x*pm.z);
                SxSp[2]=(m.x*pm.y-m.y*pm.x);

                SxSxSp[0]= (m.y*SxSp[2]-m.z*SxSp[1]);
                SxSxSp[1]= (m.z*SxSp[0]-m.x*SxSp[2]);
                SxSxSp[2]= (m.x*SxSp[1]-m.y*SxSp[0]);

                //calculate directly from J(Sxm)
                st::internal::total_ST[cellx] = prefac_sc*(m.y*saz_sot-m.z*say_sot);
                st::internal::total_ST[celly] = prefac_sc*(m.z*sax_sot-m.x*saz_sot);
                st::internal::total_ST[cellz] = prefac_sc*(m.x*say_sot-m.y*sax_sot);

                st::internal::ast[cellx] = -aj*SxSxSp[0];
                st::internal::ast[celly] = -aj*SxSxSp[1];
                st::internal::ast[cellz] = -aj*SxSxSp[2];

                st::internal::nast[cellx] = bj*SxSp[0];
                st::internal::nast[celly] = bj*SxSp[1];
                st::internal::nast[cellz] = bj*SxSp[2];
               }
            } // end of Mn2 sot
                      
            //up process
            for(int cell=idx+1; cell<idx+num_microcells_per_stack; ++cell) {

               // if(st::internal::spin_acc_sign[cell] != 1) continue;
               // calculate cell id's
               const int cellx = 3*cell+0;
               const int celly = 3*cell+1;
               const int cellz = 3*cell+2;

               // calculate previous cell id's
               const int pcellx = 3*(cell-1)+0;
               const int pcelly = 3*(cell-1)+1;
               const int pcellz = 3*(cell-1)+2;

               st::internal::three_vector_t  m(0.0,  0.0,  0.0);
               // copy array values to temporaries for readability
               if(st::internal::sot_sa_source[cell]) {
                  m.x = 0;
                  m.y = -1.0;
                  m.z = 0.0; // Au cell artificial magnetisations
               } else {
                  m.x = st::internal::m[cellx];
                  m.y = st::internal::m[celly];
                  m.z = st::internal::m[cellz]; // current cell magnetisations
               }
               st::internal::three_vector_t pm(0.0,0.0,0.0);
               if(st::internal::sot_sa_source[cell-1]) {
                  pm.x = 0;
                  pm.y = 1.0;
                  pm.z = 0.0; // Au cell artificial magnetisations
               } else {
                  pm.x =st::internal::m[pcellx];
                  pm.y = st::internal::m[pcelly];
                  pm.z = st::internal::m[pcellz]; // current cell magnetisations
               }

               const double modm = sqrt(m.x*m.x + m.y*m.y + m.z*m.z);
               const double pmodm = sqrt(pm.x*pm.x + pm.y*pm.y + pm.z*pm.z);

               // Check for zero magnetization in normalization
               if(modm > 1.e-11){
                  m.x = m.x/modm;
                  m.y = m.y/modm;
                  m.z = m.z/modm;
               }
               else{
                  m.x = 0.0;
                  m.y = 0.0;
                  m.z = 0.0;
               }
               if(pmodm > 1.e-11){
                  pm.x = pm.x/pmodm;
                  pm.y = pm.y/pmodm;
                  pm.z = pm.z/pmodm;
               }
               else{
                  pm.x = 0.0;
                  pm.y = 0.0;
                  pm.z = 0.0;
               }

               //---------------------------------------------------------------------
               // Step 1 calculate coordinate transformation m -> m' for current cell
               //---------------------------------------------------------------------

               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m, itm);

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
               const double Bc = st::internal::beta_cond[cell]; // beta
               const double Bd = st::internal::beta_diff[cell]; // beta_prime
               const double Do = st::internal::diffusion[cell];
                st::internal::three_vector_t jm0(st::internal::j_up[pcellx],st::internal::j_up[pcelly],st::internal::j_up[pcellz]);
               // // if(st::internal::sot_sa_source[cell-1]) {
               //    jm0.y += st::internal::initial_theta*je;
               // // }

               //  Calculate gradient dsacc/dx
               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m.x*m.x - twoDo;     M.xy = preD*m.y*m.x;             M.xz = preD*m.z*m.x;
               M.yx = preD*m.x*m.y;             M.yy = preD*m.y*m.y - twoDo;     M.yz = preD*m.z*m.y;
               M.zx = preD*m.x*m.z;             M.zy = preD*m.y*m.z;             M.zz = preD*m.z*m.z - twoDo;

               // if(st::internal::sot_sa ) {
                  V.x = jm0.x;// - Bc*je*m.x;
                  V.y = jm0.y;// + Bc*je*m.y;
                  V.z = jm0.z;// - Bc*je*m.z;
               // } else {
               //    V.x = jm0.x - Bc*je*m.x;
               //    V.y = jm0.y - Bc*je*m.y;
               //    V.z = jm0.z - Bc*je*m.z;
               // }
               const st::internal::three_vector_t divm_0 = gaussian_elimination(M, V);

               // Calculate mp(0), c and d
               const double i_lsdl = 1.0/st::internal::lambda_sdl[cell];
               double mp_inf = st::internal::sa_infinity[cell];
               //   if(st::internal::sot_sa && st::internal::sot_sa_source[cell]) mp_inf *= program::fractional_electric_field_strength;
               const double a = st::internal::a[cell];
               const double b = st::internal::b[cell];
               const double two_a = 2.0*a;
               const double two_b = 2.0*b;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               const st::internal::three_vector_t coeff = gaussian_elimination(M, V);

               const double mp_0 = coeff.x;
               const double c    = coeff.y;
               const double d    = coeff.z;

               //------------------------------------
               // Step 3 calculate spin accumulation
               //------------------------------------

               const double x = st::internal::micro_cell_thickness*1.0e-10; // Convert to metres
               const double cos_bx = cos(b*x);
               const double sin_bx = sin(b*x);
               const double e_xsdl = exp(-x*i_lsdl);
               const double e_ax   = exp(-a*x);
               const double prefac = (2.0*e_ax);

               const double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               const double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               const double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);

               //convert mp and m_perp
               //  double sax = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
               //  double say = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
               //  double saz = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;
              
               //add SOT contribution
               // double theta = atan2(m.y,m.x);
               //    if(theta != theta) theta = 0.0;
               // double phi = acos(m.z);
                  // if(st::internal::spin_acc_sign.at(cell) != 1.0 && st::internal::spin_acc_sign.at(cell) != -1.0) std::cout << st::internal::spin_acc_sign.at(cell) << ", " << stack << ", " << cell << std::endl;
               // double say_sot = 0.0;// sin(phi)*(1.4+1.2*sin(theta)*sin(theta))*st::internal::sa_infinity[cell]*je*1e-11*st::internal::spin_acc_sign[cell]; //j_sd / mu_b dS/S 8 10^8 / 10^7 A/m^2 = T
               // double sax_sot = 0.0;//sin(phi)*sin(2*theta)*st::internal::sa_infinity[cell]*je*1e-11*st::internal::spin_acc_sign[cell];
               // double saz_sot = 0.0;// -cos(theta)*sin(phi)*0.12*st::internal::sa_infinity[cell]*je*1e-11*st::internal::spin_acc_sign[cell];

               // double sax_sot = sin(phi)*sin(2.0*theta)*0.5*1.48e7*je*1e-11; //j_sd / mu_b dS/S 8 10^8 / 10^7 A/m^2 = T
               // double say_sot = sin(phi)*(1.2*sin(theta)+1.4)*1.48e7*je*1e-11;
               // double saz_sot = -cos(theta)*sin(phi)*0.12*1.48e7*je*1e-11;
              //  std::cout << sax << ", " << say << ", " << saz << ", " \
                                          <<sax_sot << ", " << say_sot << ", " << saz_sot << std::endl;
               //--------------------------------------------
               // Step 4 calculate the spin current (jm_end)
               //--------------------------------------------
               const double ac_bd = a*c + b*d; 
               const double ad_bc = a*d - b*c;

               const double divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl;
               const double divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx);
               const double divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx);

               const double divsax = b1.x*divsa_para + b2.x*divsa_perp2 + b3.x*divsa_perp3;
               const double divsay = b1.y*divsa_para + b2.y*divsa_perp2 + b3.y*divsa_perp3;
               const double divsaz = b1.z*divsa_para + b2.z*divsa_perp2 + b3.z*divsa_perp3;

               const double dot = m.x*divsax + m.y*divsay + m.z*divsaz;

               const double pre_jmx = divsax - Bc*Bd*m.x*dot;
               const double pre_jmy = divsay - Bc*Bd*m.y*dot;
               const double pre_jmz = divsaz - Bc*Bd*m.z*dot;

               //  double jmx = Bc*je*m.x - twoDo*pre_jmx;
               //  double jmy = Bc*je*m.y - twoDo*pre_jmy;
               //  double jmz = Bc*je*m.z - twoDo*pre_jmz;

               // if(st::internal::sot_sa) {
                  double jmx =  - twoDo*pre_jmx;
                  double jmy =  - twoDo*pre_jmy;
                  double jmz =  - twoDo*pre_jmz;
               // }
               if(st::internal::cell_natom[cell]>0) {
               //    if(st::internal::sot_sa && st::internal::sot_sa_source[cell]){
               // //       st::internal::sa[cellx] += st::internal::sa_infinity[cell]*0.0;
               // //       st::internal::sa[celly] += st::internal::sa_infinity[cell]*1.0;
               // //       st::internal::sa[cellz] += st::internal::sa_infinity[cell]*0.0;

               // //    // std::cout << st::internal::sa[3*idx+0] << ", " << st::internal::sa[3*idx+1] << ", " << st::internal::sa[3*idx+2] << std::endl;

               //       st::internal::j [cellx] += st::internal::initial_beta*je*0.0;
               //       st::internal::j [celly] += st::internal::initial_beta*je*1.0;
               //       st::internal::j [cellz] += st::internal::initial_beta*je*0.0;
               // }
                  // Save values for the spin accumulation
                  // st::internal::sa[cellx] = sax;//(sax-mp_inf*m.x)/(mp_inf*m.x);
                  // st::internal::sa[celly] = say;//(say-mp_inf*m.y)/(mp_inf*m.y);
                  // st::internal::sa[cellz] = saz;//(saz-mp_inf*m.z)/(mp_inf*m.z);

                  // st::internal::sa_sot_up[cellx] += sax;
                  // st::internal::sa_sot_up[celly] += sax;
                  // st::internal::sa_sot_up[cellz] += say;

                  // Save values for the spin current
                  st::internal::j_up[cellx] += jmx;
                  st::internal::j_up[celly] += jmy;
                  st::internal::j_up[cellz] += jmz;

                  // }  else {// Calculate spin torque energy for cell (Joules)
                  // st::internal::spin_torque[cellx] += microcell_volume * st::internal::sd_exchange[cell] * (sax) * i_e * i_muB;
                  // st::internal::spin_torque[celly] += microcell_volume * st::internal::sd_exchange[cell] * (say) * i_e * i_muB;
                  // st::internal::spin_torque[cellz] += microcell_volume * st::internal::sd_exchange[cell] * (saz) * i_e * i_muB; 
                  
                  // st::internal::spin_torque[cellx] = 0.0;// microcell_volume * st::internal::sd_exchange[cell] * (sax+sax_sot) * i_e * i_muB;
                  // st::internal::spin_torque[celly] = 0.0;//microcell_volume * st::internal::sd_exchange[cell] * (say+say_sot) * i_e * i_muB;
                  // st::internal::spin_torque[cellz] = 0.0;//microcell_volume * st::internal::sd_exchange[cell] * (saz+saz_sot) * i_e * i_muB; 
               } 
               else{
                  // Save values for the spin accumulation
                  // st::internal::sa[cellx] += st::internal::sa[pcellx];
                  // st::internal::sa[celly] += st::internal::sa[pcelly];
                  // st::internal::sa[cellz] += st::internal::sa[pcellz];

                  // st::internal::sa_sot_up[cellx] += st::internal::sa_sot_up[pcellx];
                  // st::internal::sa_sot_up[celly] += st::internal::sa_sot_up[pcelly];
                  // st::internal::sa_sot_up[cellz] += st::internal::sa_sot_up[pcellz];

                  // Save values for the spin current
                  st::internal::j_up[cellx] += st::internal::j_up[pcellx];
                  st::internal::j_up[celly] += st::internal::j_up[pcelly];
                  st::internal::j_up[cellz] += st::internal::j_up[pcellz];

                  // Calculate spin torque energy for cell (Joules)
                  // st::internal::spin_torque[cellx] += st::internal::spin_torque[pcellx];
                  // st::internal::spin_torque[celly] += st::internal::spin_torque[pcelly];
                  // st::internal::spin_torque[cellz] += st::internal::spin_torque[pcellz];
                  // }
               }
               if(st::internal::output_torque_data) {
               //--------------------------------------------
               // Step 5 calculate the spin torque of each cell
               //--------------------------------------------

               //convert M of previous cell into basis b1, b2, b3

               M.xx = b1.x;    M.xy = b2.x;    M.xz = b3.x;
               M.yx = b1.y;    M.yy = b2.y;    M.yz = b3.y;
               M.zx = b1.z;    M.zy = b2.z;    M.zz = b3.z;

               V.x = pm.x;
               V.y = pm.y;
               V.z = pm.z;

               const st::internal::three_vector_t pm_basis = gaussian_elimination(M, V);

               //const double pm_b1 = pm_basis.x; // unused variable
               const double pm_b2 = pm_basis.y;
               const double pm_b3 = pm_basis.z;

               // Calculate the spin torque coefficients describing ast and nast
               const double prefac_sc = atomcell_volume * st::internal::sd_exchange[cell] * i_e * i_muB;
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
                st::internal::coeff_ast[cell]  = aj;
                st::internal::coeff_nast[cell] = bj;

                SxSp[0]=(m.y*pm.z-m.z*pm.y);
                SxSp[1]=(m.z*pm.x-m.x*pm.z);
                SxSp[2]=(m.x*pm.y-m.y*pm.x);

                SxSxSp[0]= (m.y*SxSp[2]-m.z*SxSp[1]);
                SxSxSp[1]= (m.z*SxSp[0]-m.x*SxSp[2]);
                SxSxSp[2]= (m.x*SxSp[1]-m.y*SxSp[0]);

                //calculate directly from J(Sxm)
               //  st::internal::total_ST[cellx] = prefac_sc*(m.y*saz-m.z*say);
               //  st::internal::total_ST[celly] = prefac_sc*(m.z*sax-m.x*saz);
               //  st::internal::total_ST[cellz] = prefac_sc*(m.x*say-m.y*sax);

                st::internal::ast[cellx] = -aj*SxSxSp[0];
                st::internal::ast[celly] = -aj*SxSxSp[1];
                st::internal::ast[cellz] = -aj*SxSxSp[2];

                st::internal::nast[cellx] = bj*SxSp[0];
                st::internal::nast[celly] = bj*SxSp[1];
                st::internal::nast[cellz] = bj*SxSp[2];
               }
            } 
         
               //down process
            for(int cell=idx+num_microcells_per_stack-2; cell > idx; --cell) {

               // if(st::internal::spin_acc_sign[cell] != 1) continue;
               // calculate cell id's
               const int cellx = 3*cell+0;
               const int celly = 3*cell+1;
               const int cellz = 3*cell+2;

               // calculate previous cell id's
               const int pcellx = 3*(cell+1)+0;
               const int pcelly = 3*(cell+1)+1;
               const int pcellz = 3*(cell+1)+2;

               st::internal::three_vector_t  m(0.0,  0.0,  0.0);
               // copy array values to temporaries for readability
               if(st::internal::sot_sa_source[cell]) {
                  m.x = 0;
                  m.y = 1.0;
                  m.z = 0.0; // Au cell artificial magnetisations
               } else {
                  m.x =st::internal::m[cellx];
                  m.y = st::internal::m[celly];
                  m.z = st::internal::m[cellz]; // current cell magnetisations
               }
               st::internal::three_vector_t pm(0.0,0.0,0.0);
               if(st::internal::sot_sa_source[cell+1]) {
                  pm.x = 0;
                  pm.y = -1.0;
                  pm.z = 0.0; // Au cell artificial magnetisations
               } else {
                  pm.x =st::internal::m[pcellx];
                  pm.y = st::internal::m[pcelly];
                  pm.z = st::internal::m[pcellz]; // current cell magnetisations
               }

               const double modm = sqrt(m.x*m.x + m.y*m.y + m.z*m.z);
               const double pmodm = sqrt(pm.x*pm.x + pm.y*pm.y + pm.z*pm.z);

               // Check for zero magnetization in normalization
               if(modm > 1.e-11){
                  m.x = m.x/modm;
                  m.y = m.y/modm;
                  m.z = m.z/modm;
               }
               else{
                  m.x = 0.0;
                  m.y = 0.0;
                  m.z = 0.0;
               }
               if(pmodm > 1.e-11){
                  pm.x = pm.x/pmodm;
                  pm.y = pm.y/pmodm;
                  pm.z = pm.z/pmodm;
               }
               else{
                  pm.x = 0.0;
                  pm.y = 0.0;
                  pm.z = 0.0;
               }

               //---------------------------------------------------------------------
               // Step 1 calculate coordinate transformation m -> m' for current cell
               //---------------------------------------------------------------------

               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m, itm);

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
               const double Bc = st::internal::beta_cond[cell]; // beta
               const double Bd = st::internal::beta_diff[cell]; // beta_prime
               const double Do = st::internal::diffusion[cell];
                st::internal::three_vector_t jm0(st::internal::j_down[pcellx],st::internal::j_down[pcelly],st::internal::j_down[pcellz]);
               // // if(st::internal::sot_sa_source[cell-1]) {
               //    jm0.y += st::internal::initial_theta*je;
               // // }

               //  Calculate gradient dsacc/dx
               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m.x*m.x - twoDo;     M.xy = preD*m.y*m.x;             M.xz = preD*m.z*m.x;
               M.yx = preD*m.x*m.y;             M.yy = preD*m.y*m.y - twoDo;     M.yz = preD*m.z*m.y;
               M.zx = preD*m.x*m.z;             M.zy = preD*m.y*m.z;             M.zz = preD*m.z*m.z - twoDo;

               // if(st::internal::sot_sa ) {
                  V.x = jm0.x;// - Bc*je*m.x;
                  V.y = jm0.y;// + Bc*je*m.y;
                  V.z = jm0.z;// - Bc*je*m.z;
               // } else {
               //    V.x = jm0.x - Bc*je*m.x;
               //    V.y = jm0.y - Bc*je*m.y;
               //    V.z = jm0.z - Bc*je*m.z;
               // }
               const st::internal::three_vector_t divm_0 = gaussian_elimination(M, V);

               // Calculate mp(0), c and d
               const double i_lsdl = 1.0/st::internal::lambda_sdl[cell];
                double mp_inf = st::internal::sa_infinity[cell];
               //   if(st::internal::sot_sa && st::internal::sot_sa_source[cell]) mp_inf *= program::fractional_electric_field_strength;
               const double a = st::internal::a[cell];
               const double b = st::internal::b[cell];
               const double two_a = 2.0*a;
               const double two_b = 2.0*b;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               const st::internal::three_vector_t coeff = gaussian_elimination(M, V);

               const double mp_0 = coeff.x;
               const double c    = coeff.y;
               const double d    = coeff.z;

               //------------------------------------
               // Step 3 calculate spin accumulation
               //------------------------------------

               const double x = st::internal::micro_cell_thickness*1.0e-10; // Convert to metres
               const double cos_bx = cos(b*x);
               const double sin_bx = sin(b*x);
               const double e_xsdl = exp(-x*i_lsdl);
               const double e_ax   = exp(-a*x);
               const double prefac = (2.0*e_ax);

               const double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               const double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               const double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);

               //convert mp and m_perp
               //  double sax = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
               //  double say = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
               //  double saz = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;
              
               //add SOT contribution
               // double theta = atan2(m.y,m.x);
               //    if(theta != theta) theta = 0.0;
               // double phi = acos(m.z);
                  // if(st::internal::spin_acc_sign.at(cell) != 1.0 && st::internal::spin_acc_sign.at(cell) != -1.0) std::cout << st::internal::spin_acc_sign.at(cell) << ", " << stack << ", " << cell << std::endl;
               // double say_sot = 0.0;// sin(phi)*(1.4+1.2*sin(theta)*sin(theta))*st::internal::sa_infinity[cell]*je*1e-11*st::internal::spin_acc_sign[cell]; //j_sd / mu_b dS/S 8 10^8 / 10^7 A/m^2 = T
               // double sax_sot = 0.0;//sin(phi)*sin(2*theta)*st::internal::sa_infinity[cell]*je*1e-11*st::internal::spin_acc_sign[cell];
               // double saz_sot = 0.0;// -cos(theta)*sin(phi)*0.12*st::internal::sa_infinity[cell]*je*1e-11*st::internal::spin_acc_sign[cell];

               // double sax_sot = sin(phi)*sin(2.0*theta)*0.5*1.48e7*je*1e-11; //j_sd / mu_b dS/S 8 10^8 / 10^7 A/m^2 = T
               // double say_sot = sin(phi)*(1.2*sin(theta)+1.4)*1.48e7*je*1e-11;
               // double saz_sot = -cos(theta)*sin(phi)*0.12*1.48e7*je*1e-11;
              //  std::cout << sax << ", " << say << ", " << saz << ", " \
                                          <<sax_sot << ", " << say_sot << ", " << saz_sot << std::endl;
               //--------------------------------------------
               // Step 4 calculate the spin current (jm_end)
               //--------------------------------------------
               const double ac_bd = a*c + b*d; 
               const double ad_bc = a*d - b*c;

               const double divsa_para = (mp_inf - mp_0)*i_lsdl*e_xsdl;
               const double divsa_perp2 = prefac*(-ac_bd*cos_bx + ad_bc*sin_bx);
               const double divsa_perp3 = prefac*(-ac_bd*sin_bx - ad_bc*cos_bx);

               const double divsax = b1.x*divsa_para + b2.x*divsa_perp2 + b3.x*divsa_perp3;
               const double divsay = b1.y*divsa_para + b2.y*divsa_perp2 + b3.y*divsa_perp3;
               const double divsaz = b1.z*divsa_para + b2.z*divsa_perp2 + b3.z*divsa_perp3;

               const double dot = m.x*divsax + m.y*divsay + m.z*divsaz;

               const double pre_jmx = divsax - Bc*Bd*m.x*dot;
               const double pre_jmy = divsay - Bc*Bd*m.y*dot;
               const double pre_jmz = divsaz - Bc*Bd*m.z*dot;

               //  double jmx = Bc*je*m.x - twoDo*pre_jmx;
               //  double jmy = Bc*je*m.y - twoDo*pre_jmy;
               //  double jmz = Bc*je*m.z - twoDo*pre_jmz;

               // if(st::internal::sot_sa) {
                  double jmx =  - twoDo*pre_jmx;
                  double jmy =  - twoDo*pre_jmy;
                  double jmz =  - twoDo*pre_jmz;
               // }
               if(st::internal::cell_natom[cell]>0) {
               //    if(st::internal::sot_sa && st::internal::sot_sa_source[cell]){
               // //       st::internal::sa[cellx] += st::internal::sa_infinity[cell]*0.0;
               // //       st::internal::sa[celly] += st::internal::sa_infinity[cell]*1.0;
               // //       st::internal::sa[cellz] += st::internal::sa_infinity[cell]*0.0;

               // //    // std::cout << st::internal::sa[3*idx+0] << ", " << st::internal::sa[3*idx+1] << ", " << st::internal::sa[3*idx+2] << std::endl;

               //       st::internal::j [cellx] += st::internal::initial_beta*je*0.0;
               //       st::internal::j [celly] += st::internal::initial_beta*je*1.0;
               //       st::internal::j [cellz] += st::internal::initial_beta*je*0.0;
               // }
                  // Save values for the spin accumulation
                  // st::internal::sa[cellx] = sax;//(sax-mp_inf*m.x)/(mp_inf*m.x);
                  // st::internal::sa[celly] = say;//(say-mp_inf*m.y)/(mp_inf*m.y);
                  // st::internal::sa[cellz] = saz;//(saz-mp_inf*m.z)/(mp_inf*m.z);

                  // st::internal::sa_sot_down[cellx] += sax;
                  // st::internal::sa_sot_down[celly] += sax;
                  // st::internal::sa_sot_down[cellz] += say;

                  // Save values for the spin current
                  st::internal::j_down[cellx] += jmx;
                  st::internal::j_down[celly] += jmy;
                  st::internal::j_down[cellz] += jmz;

                  // }  else {// Calculate spin torque energy for cell (Joules)
                  // st::internal::spin_torque[cellx] += microcell_volume * st::internal::sd_exchange[cell] * (sax) * i_e * i_muB;
                  // st::internal::spin_torque[celly] += microcell_volume * st::internal::sd_exchange[cell] * (say) * i_e * i_muB;
                  // st::internal::spin_torque[cellz] += microcell_volume * st::internal::sd_exchange[cell] * (saz) * i_e * i_muB; 
                  
                  // st::internal::spin_torque[cellx] = 0.0;// microcell_volume * st::internal::sd_exchange[cell] * (sax+sax_sot) * i_e * i_muB;
                  // st::internal::spin_torque[celly] = 0.0;//microcell_volume * st::internal::sd_exchange[cell] * (say+say_sot) * i_e * i_muB;
                  // st::internal::spin_torque[cellz] = 0.0;//microcell_volume * st::internal::sd_exchange[cell] * (saz+saz_sot) * i_e * i_muB; 
                  
               } 
               else{
                  // Save values for the spin accumulation
                  // st::internal::sa[cellx] += st::internal::sa[pcellx];
                  // st::internal::sa[celly] += st::internal::sa[pcelly];
                  // st::internal::sa[cellz] += st::internal::sa[pcellz];

                  // st::internal::sa_sot_down[cellx] += st::internal::sa_sot_down[pcellx];
                  // st::internal::sa_sot_down[celly] += st::internal::sa_sot_down[pcelly];
                  // st::internal::sa_sot_down[cellz] += st::internal::sa_sot_down[pcellz];

                  // Save values for the spin current
                  st::internal::j_down[cellx] += st::internal::j_down[pcellx];
                  st::internal::j_down[celly] += st::internal::j_down[pcelly];
                  st::internal::j_down[cellz] += st::internal::j_down[pcellz];

                  // Calculate spin torque energy for cell (Joules)
                  // st::internal::spin_torque[cellx] += st::internal::spin_torque[pcellx];
                  // st::internal::spin_torque[celly] += st::internal::spin_torque[pcelly];
                  // st::internal::spin_torque[cellz] += st::internal::spin_torque[pcellz];
                  // }
               }
               if(st::internal::output_torque_data) {
               //--------------------------------------------
               // Step 5 calculate the spin torque of each cell
               //--------------------------------------------

               //convert M of previous cell into basis b1, b2, b3

               M.xx = b1.x;    M.xy = b2.x;    M.xz = b3.x;
               M.yx = b1.y;    M.yy = b2.y;    M.yz = b3.y;
               M.zx = b1.z;    M.zy = b2.z;    M.zz = b3.z;

               V.x = pm.x;
               V.y = pm.y;
               V.z = pm.z;

               const st::internal::three_vector_t pm_basis = gaussian_elimination(M, V);

               //const double pm_b1 = pm_basis.x; // unused variable
               const double pm_b2 = pm_basis.y;
               const double pm_b3 = pm_basis.z;

               // Calculate the spin torque coefficients describing ast and nast
               const double prefac_sc = microcell_volume * st::internal::sd_exchange[cell] * i_e * i_muB;
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
                st::internal::coeff_ast[cell]  = aj;
                st::internal::coeff_nast[cell] = bj;

                SxSp[0]=(m.y*pm.z-m.z*pm.y);
                SxSp[1]=(m.z*pm.x-m.x*pm.z);
                SxSp[2]=(m.x*pm.y-m.y*pm.x);

                SxSxSp[0]= (m.y*SxSp[2]-m.z*SxSp[1]);
                SxSxSp[1]= (m.z*SxSp[0]-m.x*SxSp[2]);
                SxSxSp[2]= (m.x*SxSp[1]-m.y*SxSp[0]);

                //calculate directly from J(Sxm)
               //  st::internal::total_ST[cellx] = prefac_sc*(m.y*saz-m.z*say);
               //  st::internal::total_ST[celly] = prefac_sc*(m.z*sax-m.x*saz);
               //  st::internal::total_ST[cellz] = prefac_sc*(m.x*say-m.y*sax);

                st::internal::ast[cellx] = -aj*SxSxSp[0];
                st::internal::ast[celly] = -aj*SxSxSp[1];
                st::internal::ast[cellz] = -aj*SxSxSp[2];

                st::internal::nast[cellx] = bj*SxSp[0];
                st::internal::nast[celly] = bj*SxSp[1];
                st::internal::nast[cellz] = bj*SxSp[2];
               }
            } 
         
            for(int cell=idx+1; cell<idx+num_microcells_per_stack-2; ++cell) {
               // if(st::internal::sot_sa_source[cell]) continue;

               int cellx = 3*cell;
               int celly = 3*cell+1;
               int cellz = 3*cell+2;

               int pcellx = 3*(cell-1);
               int pcelly = 3*(cell-1)+1;
               int pcellz = 3*(cell-1)+2;

               int acellx = 3*(cell+1);
               int acelly = 3*(cell+1)+1;
               int acellz = 3*(cell+1)+2;

               st::internal::three_vector_t  m(st::internal::m[cellx],  st::internal::m[celly],  st::internal::m[cellz]);
               double modm = sqrt(m.x*m.x + m.y*m.y + m.z*m.z);
               if(modm > 1.e-11){
                  m.x = m.x/modm;
                  m.y = m.y/modm;
                  m.z = m.z/modm;
               }
               else{
                  m.x = 0.0;
                  m.y = 1.0;
                  m.z = 0.0;
                  modm = 1.0;
               }
               // Check for cell magnetization greater than 1e-8 mu_B
               if(modm > 1.0e-11){
                  // Initialise inverse transformation matrix
                  set_inverse_transformation_matrix(m, itm);

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
               const double Bc = st::internal::beta_cond[cell]; // beta
               const double Bd = st::internal::beta_diff[cell]; // beta_prime
               const double Do = st::internal::diffusion[cell];
                st::internal::three_vector_t jm0(st::internal::j_up[pcellx]+st::internal::j_down[acellx], \
                                                 st::internal::j_up[pcelly]+st::internal::j_down[acelly], \
                                                 st::internal::j_up[pcellz]+st::internal::j_down[acellz]);
               // st::internal::sa[cellx] = st::internal::sa_sot_up[cellx] + st::internal::sa_sot_down[cellx];
               // st::internal::sa[celly] = st::internal::sa_sot_up[celly] + st::internal::sa_sot_down[celly];
               // st::internal::sa[cellz] = st::internal::sa_sot_up[cellz] + st::internal::sa_sot_down[cellz];

               st::internal::j[cellx] = jm0.x;
               st::internal::j[celly] = jm0.y;
               st::internal::j[cellz] = jm0.z;

            if(st::internal::sot_check) {
               st::internal::sa_sot[cellx] = (st::internal::sa_sot_up[cellx] + st::internal::sa_sot_down[cellx] - st::internal::sa_infinity[cell]*m.x)/st::internal::sa_infinity[cell];
               st::internal::sa_sot[celly] = (st::internal::sa_sot_up[celly] + st::internal::sa_sot_down[celly] - st::internal::sa_infinity[cell]*m.y)/st::internal::sa_infinity[cell];
               st::internal::sa_sot[cellz] = (st::internal::sa_sot_up[cellz] + st::internal::sa_sot_down[cellz] - st::internal::sa_infinity[cell]*m.z)/st::internal::sa_infinity[cell];
            
            } else {
               st::internal::sa_sot[cellx] = st::internal::sa_sot_up[cellx] + st::internal::sa_sot_down[cellx];
               st::internal::sa_sot[celly] = st::internal::sa_sot_up[celly] + st::internal::sa_sot_down[celly];
               st::internal::sa_sot[cellz] = st::internal::sa_sot_up[cellz] + st::internal::sa_sot_down[cellz];
            }
            

               const double twoDo = 2.0*Do;
               const double preD = twoDo*Bc*Bd;

               M.xx = preD*m.x*m.x - twoDo;     M.xy = preD*m.y*m.x;             M.xz = preD*m.z*m.x;
               M.yx = preD*m.x*m.y;             M.yy = preD*m.y*m.y - twoDo;     M.yz = preD*m.z*m.y;
               M.zx = preD*m.x*m.z;             M.zy = preD*m.y*m.z;             M.zz = preD*m.z*m.z - twoDo;

               V.x = jm0.x;// - Bc*je*m.x;
               V.y = jm0.y;// + Bc*je*m.y;
               V.z = jm0.z;// - Bc*je*m.z;

               const st::internal::three_vector_t divm_0 = gaussian_elimination(M, V);

               const double i_lsdl = 1.0/st::internal::lambda_sdl[cell];
               double mp_inf = st::internal::sa_infinity[cell];
               //   if(st::internal::sot_sa && st::internal::sot_sa_source[cell]) mp_inf *= program::fractional_electric_field_strength;
               const double a = st::internal::a[cell];
               const double b = st::internal::b[cell];
               const double two_a = 2.0*a;
               const double two_b = 2.0*b;

               M.xx = -b1.x*i_lsdl;    M.xy = (-two_a*b2.x + two_b*b3.x);    M.xz = (-two_b*b2.x - two_a*b3.x);
               M.yx = -b1.y*i_lsdl;    M.yy = (-two_a*b2.y + two_b*b3.y);    M.yz = (-two_b*b2.y - two_a*b3.y);
               M.zx = -b1.z*i_lsdl;    M.zy = (-two_a*b2.z + two_b*b3.z);    M.zz = (-two_b*b2.z - two_a*b3.z);

               V.x = divm_0.x - b1.x*mp_inf*i_lsdl;
               V.y = divm_0.y - b1.y*mp_inf*i_lsdl;
               V.z = divm_0.z - b1.z*mp_inf*i_lsdl;

               const st::internal::three_vector_t coeff = gaussian_elimination(M, V);

               const double mp_0 = coeff.x;
               const double c    = coeff.y;
               const double d    = coeff.z;

               const double x = st::internal::micro_cell_thickness*1.0e-10; // Convert to metres
               const double cos_bx = cos(b*x);
               const double sin_bx = sin(b*x);
               const double e_xsdl = exp(-x*i_lsdl);
               const double e_ax   = exp(-a*x);
               const double prefac = (2.0*e_ax);

               const double sa_para  = mp_inf + (mp_0 - mp_inf)*e_xsdl;
               const double sa_perp2 = prefac*(c*cos_bx - d*sin_bx);
               const double sa_perp3 = prefac*(c*sin_bx + d*cos_bx);

               //convert mp and m_perp
               double sax = b1.x*sa_para + b2.x*sa_perp2 + b3.x*sa_perp3;
               double say = b1.y*sa_para + b2.y*sa_perp2 + b3.y*sa_perp3;
               double saz = b1.z*sa_para + b2.z*sa_perp2 + b3.z*sa_perp3;

               st::internal::sa[cellx] = sax;
               st::internal::sa[celly] = say;
               st::internal::sa[cellz] = saz;

               st::internal::spin_torque[cellx] += atomcell_volume * st::internal::sd_exchange[cell] * (sax) * i_e * i_muB;
               st::internal::spin_torque[celly] += atomcell_volume * st::internal::sd_exchange[cell] * (say) * i_e * i_muB;
               st::internal::spin_torque[cellz] += atomcell_volume * st::internal::sd_exchange[cell] * (saz) * i_e * i_muB; 

            if(st::internal::sot_check) {    
               st::internal::spin_torque[cellx] = 0.0;// += microcell_volume * st::internal::sd_exchange[cell] * (sax) * i_e * i_muB;
               st::internal::spin_torque[celly] = 0.0;//+= microcell_volume * st::internal::sd_exchange[cell] * (say) * i_e * i_muB;
               st::internal::spin_torque[cellz] = 0.0;//+= microcell_volume * st::internal::sd_exchange[cell] * (saz) * i_e * i_muB; 
            }  
            }
         } // end of stack loop

         // Reduce all microcell spin torques on all nodes
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::spin_torque[0],st::internal::spin_torque.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            // if(st::internal::output_torque_data) MPI_Allreduce(MPI_IN_PLACE, &st::internal::total_ST[0],st::internal::total_ST.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif
         st::internal::output_microcell_data();

         st::spin_acc_time += stopwatch.elapsed_seconds();

         return;
      }

   } // end of namespace internal
} // end of namespace st
