//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sarah Jenkins, Andrea Meo and Richard F L Evans 2017.
//       All rights reserved.
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

// Vampire headers
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "atoms.hpp"
#include "vutil.hpp"

// #ifdef FFT
// #include <fftw3.h>
// #endif

//  #ifdef FFTW_OMP
#include <omp.h>
//  #endif

            //-----------------------------------------------------------------------------
            // Function to initialise dipole field calculation using FFT solver
            //-----------------------------------------------------------------------------
namespace stats {

    void  spinwave_statistic_t::initialize(){

                // if( fftw_init_threads() == 0)
                    std::cout << "Error initialising threads for FFTW!" << std::endl;

        //         int Nthreads = 4;
        //         time_step = 0;
        //         N = atoms::num_atoms;
        //         //gamma point (0,0,0)                       pt. 0
        //         //Sigma line:   (v,v,0)                     pt. 1-19
        //         //M point   (1/2,1/2,0)                     pt. 20
        //         //V line: (1/2,1/2,v)                       pt. 21-39
        //         //A point   (1/2,1/2,1/2)                   pt. 40
        //         //S line (v,v,1/2)                          pt. 41-59
        //         //Z point (0,0,1/2)                         pt. 60
        //         //U line (v,0,1/2)                          pt. 61-69
        //         //R point (1/2,0,1/2)                       pt. 70
        //         //W line (1/2,0,v)                          pt. 71-89   
        //         //X point (1/2,0,0)                         pt. 90
        //         //Delta line  (v,0,0)                      pt. 91-99
             
        //         i_atoms = 230;//2+ atoms::num_atoms;        
        //         r_atoms = 230;
        //         K_points = 100;
        //         fft_coefficients = (double*) fftw_malloc( sizeof(double) * 4 * K_points);
        //         r_cutoff = (double*) fftw_malloc( sizeof(double) * K_points);
        //              vutil::vtimer_t timer;

        //         //   start timer
        //         timer.start();

        //        // tetrahedral points/lines
        //         for(int k = 0; k < K_points; k++) {
        //             // if(k == 0)        {fft_coefficients[4*k] = 0.0; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 0.0;}
        //              if (k < 20)  {r_cutoff[k] = 9.5;fft_coefficients[4*k] = M_PI*0.05*k/4.69; fft_coefficients[4*k +1] = M_PI*0.05*k/4.69; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 0.33333;}
        //             else if (k == 20) {r_cutoff[k] = 9.5;fft_coefficients[4*k] = M_PI*1/4.69; fft_coefficients[4*k +1] = M_PI*1/4.69; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 0.33333;}
        //             else if (k < 40)  {r_cutoff[k] = 9.5;fft_coefficients[4*k] = M_PI*1/4.69; fft_coefficients[4*k +1] = M_PI*1/4.69; fft_coefficients[4*k+2] = M_PI*0.05*(k-20)/3.18; fft_coefficients[4*k+3] = 0.33333;}
        //             else if (k == 40) {r_cutoff[k] = 9.5;fft_coefficients[4*k] = M_PI*1/4.69; fft_coefficients[4*k +1] = M_PI*1/4.69; fft_coefficients[4*k+2] = M_PI*1/3.18; fft_coefficients[4*k+3] = 0.33333;}
        //             else if (k < 60)  {r_cutoff[k] = 9.5;fft_coefficients[4*k] = M_PI*0.05*(60-k)/4.69; fft_coefficients[4*k +1] = M_PI*0.05*(60-k)/4.69; fft_coefficients[4*k+2] = M_PI*1/3.18; fft_coefficients[4*k+3] = 0.333333;}
        //             else if (k == 60) {r_cutoff[k] = 9.5;fft_coefficients[4*k] = 0.0; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = M_PI*1/3.18; fft_coefficients[4*k+3] = 0.3333333;}
        //             else if (k < 70)  {r_cutoff[k] = 9.5;fft_coefficients[4*k] = M_PI*0.1*(k-60)/4.69; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = M_PI*1/3.18; fft_coefficients[4*k+3] = 0.3333333;}
        //             else if (k == 70) {r_cutoff[k] = 9.5;fft_coefficients[4*k] = M_PI*1/4.69; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = M_PI*1/3.18; fft_coefficients[4*k+3] = 0.33333;}
        //             else if (k < 90)  {r_cutoff[k] = 9.5;fft_coefficients[4*k] = M_PI*1/4.69; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = M_PI*0.1*(90-k)/3.18; fft_coefficients[4*k+3] = 0.33333 ;}
        //             else if (k == 90)  {r_cutoff[k] = 9.5;fft_coefficients[4*k] = M_PI*1/4.69; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 0.33333 ;}
        //             else if (k < 100)  {r_cutoff[k] = 9.5;fft_coefficients[4*k] = M_PI*0.1*(100-k)/4.69; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 0.33333 ;}
        //         }
        //         // Nthreads = omp_get_max_threads();
        //         // std::cout << "Planning FFT with Nthreads = " << Nthreads << std::endl;       
                
            
        //         // for(int k = 1; k < K_points; k++) {
                    
        //         //     //if(k == 0)        {r_cutoff[k] = 0.001;fft_coefficients[4*k] = 0.0; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 0.0;} // gamma
        //         //     if (k < 20)  {r_cutoff[k] = 7;fft_coefficients[4*k] = M_PI*0.05*k/2.87; fft_coefficients[4*k +1] =M_PI*0.05*k/2.87; fft_coefficients[4*k+2] = M_PI*0.05*k/2.87; fft_coefficients[4*k+3] = 0.333333;} // Delta line
        //         //     else if(k == 20)  {r_cutoff[k] = 7;fft_coefficients[4*k] = M_PI*1/2.87; fft_coefficients[4*k +1] =  M_PI*1/2.87; fft_coefficients[4*k+2] =  M_PI*1/2.87; fft_coefficients[4*k+3] = 0.3333333;} // H 
        //         //     else if (k < 40)  {r_cutoff[k] = 7;fft_coefficients[4*k] = M_PI*0.05*(40-k)/2.87; fft_coefficients[4*k +1] =M_PI*1/2.87; fft_coefficients[4*k+2] = M_PI*0.05*(40-k)/2.87; fft_coefficients[4*k+3] = 0.333333;} // G line
        //         //     else if (k == 40) {r_cutoff[k] = 10;fft_coefficients[4*k] = 0.0; fft_coefficients[4*k +1] =M_PI*1/2.87; fft_coefficients[4*k+2] = 0; fft_coefficients[4*k+3] = 0.33333;} // N
        //         //     else if (k < 50)  {r_cutoff[k] = 7;fft_coefficients[4*k] = 0.0; fft_coefficients[4*k +1] =M_PI*0.1*(50-k)/2.87; fft_coefficients[4*k+2] = 0; fft_coefficients[4*k+3] = 0.333333;} // Sigma line
        //         //     else if(k == 50)  {r_cutoff[k] = 10;fft_coefficients[4*k] = 0.0; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 0.0;} // gamma
        //         //     else if (k < 70)  {r_cutoff[k] = 7;fft_coefficients[4*k] = M_PI*0.025*(k-50)/2.87; fft_coefficients[4*k +1] = M_PI*0.025*(k-50)/2.87; fft_coefficients[4*k+2] = M_PI*0.025*(k-50)/2.87; fft_coefficients[4*k+3] = 0.333333;} // Sigma line
        //         //     else if(k == 70)  {r_cutoff[k] = 7;fft_coefficients[4*k] = M_PI*0.5/2.87; fft_coefficients[4*k +1] = M_PI*0.5/2.87; fft_coefficients[4*k+2] = M_PI*0.5/2.87; fft_coefficients[4*k+3] = 0.33333333;} // Rho 
        //         //     else if(k < 90)   {r_cutoff[k] = 7;fft_coefficients[4*k] = (M_PI*0.025*(k-70)/2.87)+(M_PI*0.5/2.87); fft_coefficients[4*k +1] = (M_PI*0.025*(k-70)/2.87)+(M_PI*0.5/2.87); fft_coefficients[4*k+2] = (M_PI*0.025*(k-70)/2.87)+(M_PI*0.5/2.87); fft_coefficients[4*k+3] = 0.33333333;} // H 
        //         // }
        //         fftw_plan_with_nthreads(Nthreads);

        //         time_range = (sim::total_time)/frequency_step; //dt / dt
        //         // const double freq_hist_range = ((M_PI / frequency_step))*1e4 - ((M_PI / sim::total_time))*1e4;
        //         freq_hist_cutoff[0] = 0.0;//((M_PI / sim::total_time))*1e4;
        //         freq_hist_cutoff[1] = ((M_PI / frequency_step))*1e4;
        //          const double freq_hist_range = freq_hist_cutoff[1] - freq_hist_cutoff[0];
        //        // std::cout <<   \
        //         " spinwave fft width: " << floor(freq_hist_range) << " (samples), for " << 100.0*time_range/floor(freq_hist_range) << " bins. From " << ((M_PI / sim::total_time))*1e4 << " (THz) to  " \
        //                                 << ((M_PI / frequency_step))*1e4 << " (THz) " <<  std::endl;
                
        //         freq_hist_step = 0.01;//freq_hist_range/freq_hist_bins;
        //         freq_hist_bins = int(freq_hist_range/freq_hist_step);
                
        //         S_t = (double*) fftw_malloc( sizeof(double) * time_range * N * i_atoms );
        //         // S_r = (double*) fftw_malloc( sizeof(double) * time_range * N * i_atoms * 3.0);
        //         // S_o = (double*) fftw_malloc( sizeof(double) * time_range * N * i_atoms * 3.0);
        //         r_i = (double*) fftw_malloc( sizeof(double) *2*r_atoms *N* (K_points-1));
        //         r_s = (int*) fftw_malloc( sizeof(int) * i_atoms * N);
        //         // complex K-space
        //         S_i = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * time_range * N  * i_atoms);
                
        //         S_0 = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * N * i_atoms);

        //         // Calculate memory requirements and inform user
        //         const double mem = double(N) * i_atoms  *time_range*( sizeof(double)  + (1)*sizeof(fftw_complex)) / 1.0e6;
        //         zlog << zTs() << "Atomistic FFT dipole field calculation has been enabled and requires " << mem << " MB of RAM for " << time_range << " collections" << std::endl;
        //         std::cout     << "Atomistic FFT dipole field calculation has been enabled and requires " << mem << " MB of RAM for " << time_range << " collections" << std::endl;
        //         std::cout <<   \
        //         " spinwave fft width: " << freq_hist_bins << " (samples) at " << freq_hist_step << " (THz). From " << freq_hist_cutoff[0] << " (THz) to  " \
        //                                 << freq_hist_cutoff[1] << " (THz) " <<  std::endl;
        //         zlog << zTs() <<   \
        //         " spinwave fft width: " << freq_hist_bins << " (samples) at " << freq_hist_step << " (THz). From " << freq_hist_cutoff[0] << " (THz) to  " \
        //                                 << freq_hist_cutoff[1]<< " (THz) " <<  std::endl;
                
    
        //     double total = 0;
        //     double min = i_atoms;
        //     double max = 0;
        //     #pragma omp parallel for reduction(+:total) num_threads(6) reduction(max:max) reduction(min:min)
        //     for(int i = 0; i < atoms::num_atoms; i++ ) {

        //         const double x = atoms::x_coord_array[i];
        //         const double y = atoms::y_coord_array[i];
        //         const double z = atoms::z_coord_array[i];
             
        //         std::vector<int > keep;
        //         keep.resize(atoms::num_atoms*2, 0);
        //         if(keep.size() != atoms::num_atoms*2) std::cout << keep.size() << std::endl;
        //         int count_r = 1;

        //         for(int k = 1; k < K_points; k++) {
        //             double fft_x = fft_coefficients[4*k];
        //             double fft_y = fft_coefficients[4*k+1];
        //             double fft_z = fft_coefficients[4*k+2];
        //             // double normalisation = fft_coefficients[4*k+3];
        //             int count_k = 1;
        //             for(int nn = 0; nn < atoms::num_atoms; ++nn) {
        //                 if(nn == i) continue;
        //                 // const int natom = atoms::neighbour_list_array[nn];
        //                  double d_x = atoms::x_coord_array[nn] - x;
        //                  double d_y = atoms::y_coord_array[nn] - y;
        //                  double d_z = atoms::z_coord_array[nn] - z;

        //                 if(d_x > cs::system_dimensions[0]*0.5)   d_x -=  cs::system_dimensions[0];//cos(theta)*sin(phi);
        //                 else if (d_x < -1.0*cs::system_dimensions[0]*0.5) d_x += cs::system_dimensions[0];
                        
        //                 if(d_y > cs::system_dimensions[1]*0.5)  d_y -= cs::system_dimensions[1];//sin(theta)*sin(phi);
        //                 else if (d_y < -1.0*cs::system_dimensions[1]*0.5) d_y += cs::system_dimensions[1];
                        
        //                 if(d_z > cs::system_dimensions[2]*0.5)   d_z -= cs::system_dimensions[2];//cos(phi);
        //                 else if(d_z < -1.0*cs::system_dimensions[2]*0.5)  d_z += cs::system_dimensions[2];//cos(phi);

        //                 if(std::abs(d_x) > r_cutoff[k]|| \
        //                    std::abs(d_y) > r_cutoff[k] || \
        //                    std::abs(d_z) > r_cutoff[k]) continue;

        //                 double cutoff =  sin(fft_x*std::abs(d_x)/2.0)+sin(fft_y*std::abs(d_y)/2.0)+sin(fft_z*std::abs(d_z)/2.0);
        //                 if(cutoff > 0.3) {
        //                   r_i[2*i*r_atoms*(K_points-1) + 2*(k-1)*r_atoms + 2*count_k + 0] = cutoff;
        //                   count_k++;
        //                   if(count_k > (r_atoms-2)) { std::cout <<k << ", " << r_cutoff[k] << ", " << count_k << std::endl; break; }
        //                   if(keep[nn*2] == 1) {
        //                     r_i[2*i*r_atoms*(K_points-1) + 2*(k-1)*r_atoms + 2*count_k + 1] = keep[nn*2 +1];
        //                   } else {  
        //                   r_s[i*i_atoms + count_r] = nn;
        //                   keep[nn*2] = 1;
        //                   keep[nn*2 +1] = count_r;
        //                   r_i[2*i*r_atoms*(K_points-1) + 2*(k-1)*r_atoms + 2*count_k + 1] = count_r;
        //                   count_r++;
        //                   if(count_r > (i_atoms -2)) { std::cout << "bad " << count_r << std::endl; break; }
        //                   }
        //                 }
        //             } 
        //             r_i[2*i*r_atoms*(K_points-1) + 2*(k-1)*r_atoms + 0] = count_k;
        //         }
        //         r_s[i*i_atoms + 0] = count_r;
        //         if(count_r > max) max = count_r;
        //         if(count_r < min) min = count_r;
        //         total += count_r;
        //         // if(electron > N-1) std::cout << electron << std::endl;
        //     }

        //     timer.stop();
        //    std::cout << timer.elapsed_time() << " for " << total/N  << " out of " << i_atoms << " fft_r; range: " << min << " to " << max << std::endl;

        //     for(int a = 0; a < atoms::num_atoms; a++) {
        //             // const int start = atoms::neighbour_list_start_index[a];
   		// 		    // const int end   = atoms::neighbour_list_end_index[a]+1;
        //             // if(end - start != 14) std::cout << start << ", " << end << std::endl;
        //         int count = 1;
        //         S_0[a*i_atoms+0][0] = atoms::x_spin_array[a]; //zero-gamma point
        //         S_0[a*i_atoms+0][1] = atoms::y_spin_array[a]; //zero-gamma point

        //         for(int nn = 1; nn < r_s[a*i_atoms]; ++nn) {
        //             const int natom = r_s[a*i_atoms+nn];
        //             S_0[a*i_atoms+count][0] =  atoms::x_spin_array[natom];
        //             S_0[a*i_atoms+count][1] =  atoms::y_spin_array[natom];
        //             count++;
        //         }
        //     }

        //         // create FFTW plans to act on the M and H arrays
        //         int n[] = {time_range};
        //         int rank = 1;
        //         int howmany = N*i_atoms;
        //         int idist = 1, odist = 1;
        //         int istride = N*i_atoms, ostride = N*i_atoms;
        //         int *inembed = n, *onembed = n;

        //         // Here we plan the transforms, making use of the inter-leaved memory
        //         // From real space to K-space uses a real to complex
        //         plan_S = fftw_plan_many_dft_c2r( rank, n, howmany,
        //                 S_i, inembed, istride, idist,
        //                 S_t, onembed, ostride, odist,
        //                 FFTW_MEASURE);

              
        //         // std::cout << "plan_S complete; r time" << std::endl;
        //         // int r[] = {N*time_range, N*time_range, N*time_range};
        //         // howmany = i_atoms;
        //         // rank = 3;
        //         // idist = 3; odist = 3;
        //         // istride = i_atoms; ostride = i_atoms;
        //         // int *inembed_r = r, *onembed_r = r;
        //         //  plan_r = fftw_plan_many_r2r( rank, r, howmany,
        //         //         S_r, inembed_r, istride, idist,
        //         //         S_o, onembed_r, ostride, odist,
        //         //         FFTW_MEASURE, FFTW_R2HC);
        //         fftw_free(fft_coefficients);
        //         fftw_free(r_cutoff);
        //         initialised = true;

}

void spinwave_statistic_t::reset(){

        // // int electron = 0;
        // #pragma omp parallel for num_threads(2)
        // for(int a = 0; a < atoms::num_atoms; a++) {
        //     // const int start = atoms::neighbour_list_start_index[a];
   		//     // const int end   = atoms::neighbour_list_end_index[a]+1;
            
        //     int count = 0;
        //     S_0[a*i_atoms+count][0] = atoms::x_spin_array[a];
        //     S_0[a*i_atoms+count][1] = atoms::y_spin_array[a];
        //     count++;
        //     for(int nn = 1; nn < r_s[a*i_atoms]; ++nn) {
        //                 const int natom = r_s[a*i_atoms+nn];
        //         S_0[a*i_atoms+count][0] = atoms::x_spin_array[natom];
        //         S_0[a*i_atoms+count][1] = atoms::y_spin_array[natom];
        //         count++;
        //     }
        //     // electron++;
        // }
        // time_step = 0;
}
    
void spinwave_statistic_t::update() {

        // if(sim::time % frequency_step != 0) return; //update will call with the statistics
        //                                             //only want spinwave update with frequency setting

        // if(!initialised) {
        //     std::cout << "FFT spinwave has been called but not initialised." << std::endl;                    exit(-1);
        // }
        // // int electron = 0;
       
        // #pragma omp parallel for num_threads(2) 
        // for( int atom = 0; atom < atoms::num_atoms; atom++) {
        //     // const int start = atoms::neighbour_list_start_index[atom];
   		//     // const int end   = atoms::neighbour_list_end_index[atom]+1;
            
        //     int count = 0;
        //     const double S_x = atoms::x_spin_array[atom];
        //     const double S_y = atoms::y_spin_array[atom];
        //     spin_correlation(S_i[time_step*N*i_atoms + atom*i_atoms + 0 ], S_x , S_y , S_0[atom*i_atoms + 0]);
        //     count++;
        //     for(int nn = 1; nn < r_s[atom*i_atoms]; ++nn) {
        //                 // const int natom = r_s[atom*i_atoms+nn];
        //         spin_correlation(S_i[time_step*N*i_atoms + atom*i_atoms + count], S_x, S_y, S_0[atom*i_atoms + count]);
        //         count++;
        //     }
        //     // electron++;
        // }

        // time_step++;
        // return;

    }

            //-----------------------------------------------------------------------------
            // Function to finalize FFT solver and release memory
            //-----------------------------------------------------------------------------
void spinwave_statistic_t::finalize() {

         // instantiate timer
                // vutil::vtimer_t timer;

                // //   start timer
                // timer.start();

        // fftw_execute(plan_S);

    // for(int t = 0; t < time_range; t++) {
    //     for(int a = 0; a < N; a++) {
    //         for(int k = 0; k < i_atoms; k++) {
    //             S_r[t*N*i_atoms*3 + a + k + 0] = S_t[t*N*i_atoms + a + k]*r_i[a + k + 0];
    //             S_r[t*N*i_atoms*3 + a + k + 1] = S_t[t*N*i_atoms + a + k]*r_i[a + k + 1];
    //             S_r[t*N*i_atoms*3 + a + k + 2] = S_t[t*N*i_atoms + a + k]*r_i[a + k + 2];
    //         }
    //     }
    // }
    //     fftw_execute(plan_r);
//         timer.stop();
//         std::cout << "Spinwave compute time = " << timer.elapsed_time() << std::endl;
//         zlog << zTs() << "Spinwave compute time = " << timer.elapsed_time()  << std::endl;
//         timer.start();

//         std::ofstream ofile;
//          std::ofstream ofile_1;
//         // std::vector<std::ofstream> ofile_1;
//         // ofile_1.resize(2);
//         std::vector< double > fft_frequency_hist;
//         // fft_frequency_hist.resize(1);

//         // fft_frequency_hist.resize(2);
        
//             //  std::vector<std::vector< double> > fft_integration;
     
//         freq_hist_step = 1e-3;
//         freq_hist_bins = int(50/freq_hist_step);
//         fft_frequency_hist.resize(4*K_points*freq_hist_bins, 0.0);
       
//         // std::cout << freq_hist_bins << std::endl;
//         // for(int k = 0; k < 1; k++) {
//             // fft_frequency_hist.resize(freq_hist_bins, 0);
//         // }

//         ofile.open("atomistic_spinwave_freq");
//         ofile_1.open("atomistic_spinwave_spatial");
//         // ofile_1[0].open("atomistic_spinwave_rawfreq_0");
//         // ofile_1[1].open("atomistic_spinwave_rawfreq_1");
//         // ofile_1[2].open("atomistic_spinwave_rawfreq_2");
//         // ofile_1[3].open("atomistic_spinwave_rawfreq_3");
//         const double normalise = 1e6/sim::total_time/M_PI;//1e4*1.0/sqrt(1.0*M_PI);
       
 
//     //  std::cout << freq_hist_bins << std::endl;

//     fftw_free(S_i);
//     fftw_free(S_0);
    
//     fftw_destroy_plan(plan_S);
//    int h_value = 0;
//     #pragma omp parallel for num_threads(4) private(h_value)
//     for(int t = 0; t < time_range; t++) {
//         int thread =  0;//omp_get_thread_num();
//         for(int a = 0; a < N; a++) {   
//             double value = std::abs(S_t[t*N*i_atoms + a*i_atoms + 0]*normalise);
//             // if(value < 50) {
//              h_value = int(std::min(double(freq_hist_bins)-1, floor(value/freq_hist_step))); 
//             //  if(h_value == 0 && value != 0.0) std::cout << value << std::endl;
//             fft_frequency_hist[thread*K_points*freq_hist_bins + 0*freq_hist_bins + h_value] ++; //zero-gamma point
//             // }
//             for(int r = 1; r < K_points; r++) {

//                 for(int k = 1; k < r_i[2*a*r_atoms*(K_points-1) + 2*(r-1)*r_atoms + 0]; k++) {
//                     int natom = r_i[2*a*(K_points-1)*r_atoms + 2*(r-1)*r_atoms + 2*k + 1];   
//                     value = std::abs(S_t[t*N*i_atoms + a*i_atoms + natom]*normalise);
//                     // if(value > 50.0) continue;
//                     h_value = int(std::min(double(freq_hist_bins)-1, floor(value/freq_hist_step)));
//                     fft_frequency_hist[thread*K_points*freq_hist_bins + r*freq_hist_bins + h_value] += r_i[2*a*(K_points-1)*r_atoms +2*(r-1)*r_atoms + 2*k];   
//                 }  
//             }
//         }
//     }

//      timer.stop();
//         std::cout << "Spinwave freq output time = " << timer.elapsed_time() << std::endl;
//                 zlog << zTs() << "Spinwave freq output time = " << timer.elapsed_time()  << std::endl;
//        timer.start();
  
    
//     for(int f = 0; f < freq_hist_bins; f++) {
//             ofile <<  f*freq_hist_step + freq_hist_cutoff[0] << ", ";
//             ofile_1 <<  f*freq_hist_step + freq_hist_cutoff[0] << ", ";
//         for(int r = 0; r < K_points; r++) {
//             if(r == 0) ofile << double(fft_frequency_hist[0*K_points*freq_hist_bins + r*freq_hist_bins + f]\
//                                       +fft_frequency_hist[1*K_points*freq_hist_bins + r*freq_hist_bins + f]\
//                                       +fft_frequency_hist[2*K_points*freq_hist_bins + r*freq_hist_bins + f]\
//                                       +fft_frequency_hist[3*K_points*freq_hist_bins + r*freq_hist_bins + f]) << ", ";
//             else ofile_1 << double(fft_frequency_hist[0*K_points*freq_hist_bins + r*freq_hist_bins + f]\
//                                   +fft_frequency_hist[1*K_points*freq_hist_bins + r*freq_hist_bins + f]\
//                                   +fft_frequency_hist[2*K_points*freq_hist_bins + r*freq_hist_bins + f]\
//                                   +fft_frequency_hist[3*K_points*freq_hist_bins + r*freq_hist_bins + f])/r_atoms << ", ";
//         } 
//         ofile << "\n";
//         ofile_1 << "\n";
//     }


//         // ofile << S_t[t*N*i_atoms + a + o]*normalise << "\t";   
//         ofile.close();//("/spinwave/atomistic_spinwave" + std::to_string(o) + ".txt");
//         ofile_1.close();
//         // ofile_1[1].close();
//         // ofile_1[2].close();
//         // ofile_1[3].close();

//         //   ofile.open("atomistic_spinwave_spatial");
          
//     // for(int r = 0; r < K_points+1; r++) {
//     //     ofile << r << "\t" << (fft_integration[0][r]+fft_integration[1][r])/N/(i_atoms-1)/sqrt(2.0*M_PI) << "\n";
//     // } 
//     // ofile.close();

//          timer.stop();
                
       
//         std::cout << "Spinwave spatial output time = " << timer.elapsed_time() << std::endl;
//                 zlog << zTs() << "Spinwave spatial output time = " << timer.elapsed_time()  << std::endl;
//         // timer.start();
           
//                 // Free memory from FFT complex variables
               
//                 // fftw_free(r_i);
//                 fftw_free(S_t);
//                 // fftw_free(S_r);
//                  fftw_free(r_s);
//                 fftw_free(r_i);
               
//                 //  fftw_destroy_plan(plan_r);

//                 fftw_cleanup_threads();

//                 return;

    }
}
