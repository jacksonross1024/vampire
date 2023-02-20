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
#include <fftw3.h>
// #endif

// #ifdef FFTW_OMP
#include <omp.h>
// #endif

            //-----------------------------------------------------------------------------
            // Function to initialise dipole field calculation using FFT solver
            //-----------------------------------------------------------------------------
namespace stats {

    void  spinwave_statistic_t::initialize(){

                if( fftw_init_threads() == 0)
                    std::cout << "Error initialising threads for FFTW!" << std::endl;

                int Nthreads = 4;
                time_step = 0;
                N = atoms::num_atoms;
                //gamma point (0,0,0)                       pt. 0
                //delta line: 10 pts.  (v,0,0) {0<v<1/2}    pt. 1-9
                //X point   (1/2,0,0)                       pt. 10
                //Y line: 10 pts. (1/2,v,0) {0<1/2<v}       pt. 11-19
                //M point   (1/2,1/2,0)                     pt. 20
                //V line (1/2,1/2,v)                        pt. 21-29
                //A point (1/2,1/2,1/2)                     pt. 30
                //A line (v,v,1/2)                          pt. 31-39
                //Z point (0,0,1/2)                         pt. 40
                //lambda line (0,0,v)                       pt. 41-49        
                i_atoms = 1500;//2+ atoms::num_atoms;
                double cutoff = 0.225;
                  K_points = 90;
               fft_coefficients = (double*) fftw_malloc( sizeof(double) * 4 * K_points);
               // tetrahedral points/lines
                // for(int k = 0; k < 50; k++) {
                //     if(k == 0)        {fft_coefficients[4*k] = 0.0; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 0.0;}
                //     else if (k < 10)  {fft_coefficients[4*k] = M_PI*0.1*k/4.69; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 1.0;}
                //     else if (k == 10) {fft_coefficients[4*k] = M_PI*1/4.69; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 1.0;}
                //     else if (k < 20)  {fft_coefficients[4*k] = M_PI*1/4.69; fft_coefficients[4*k +1] = M_PI*0.1*(k-10)/4.69; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 0.50;}
                //     else if (k == 20) {fft_coefficients[4*k] = M_PI*1/4.69; fft_coefficients[4*k +1] = 1/4.69; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 0.50;}
                //     else if (k < 30)  {fft_coefficients[4*k] = M_PI*1/4.69; fft_coefficients[4*k +1] = M_PI*1/4.69; fft_coefficients[4*k+2] = M_PI*0.1*(k-20)/3.18; fft_coefficients[4*k+3] = 0.333333;}
                //     else if (k == 30) {fft_coefficients[4*k] = M_PI*1/4.69; fft_coefficients[4*k +1] = M_PI*1/4.69; fft_coefficients[4*k+2] = M_PI*1/3.18; fft_coefficients[4*k+3] = 0.3333333;}
                //     else if (k < 40)  {fft_coefficients[4*k] = M_PI*0.1*(40-k)/4.69; fft_coefficients[4*k +1] = M_PI*0.1*(40-k)/4.69; fft_coefficients[4*k+2] = M_PI*1/3.18; fft_coefficients[4*k+3] = 0.3333333;}
                //     else if (k == 40) {fft_coefficients[4*k] = 0.0; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = M_PI*1/3.18; fft_coefficients[4*k+3] = 1.0;}
                //     else if (k < 50)  {fft_coefficients[4*k] = 0.0; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = M_PI*0.1*(50-k)/3.18; fft_coefficients[4*k+3] = 1.0;}
                // }
                // Nthreads = omp_get_max_threads();
                // std::cout << "Planning FFT with Nthreads = " << Nthreads << std::endl;
              
                for(int k = 0; k < K_points; k++) {
                    if(k == 0)        {fft_coefficients[4*k] = 0.0; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 0.0;} // gamma
                    else if (k < 20)  {fft_coefficients[4*k] = M_PI*0.05*k/2.87; fft_coefficients[4*k +1] =M_PI*0.05*k/2.87; fft_coefficients[4*k+2] = M_PI*0.05*k/2.87; fft_coefficients[4*k+3] = 0.333333;} // Delta line
                    else if(k == 20)  {fft_coefficients[4*k] = M_PI*1/2.87; fft_coefficients[4*k +1] =  M_PI*1/2.87; fft_coefficients[4*k+2] =  M_PI*1/2.87; fft_coefficients[4*k+3] = 0.3333333;} // H 
                    else if (k < 40)  {fft_coefficients[4*k] = M_PI*0.05*(40-k)/2.87; fft_coefficients[4*k +1] =M_PI*1/2.87; fft_coefficients[4*k+2] = M_PI*0.05*(20-k)/2.87; fft_coefficients[4*k+3] = 0.333333;} // G line
                    else if (k == 40) {fft_coefficients[4*k] = 0.0; fft_coefficients[4*k +1] =M_PI*1/2.87; fft_coefficients[4*k+2] = 0; fft_coefficients[4*k+3] = 1.0;} // N
                    else if (k < 50)  {fft_coefficients[4*k] = 0.0; fft_coefficients[4*k +1] =M_PI*0.05*(50-k)/2.87; fft_coefficients[4*k+2] = 0; fft_coefficients[4*k+3] = 1.0;} // Sigma line
                    else if(k == 50)  {fft_coefficients[4*k] = 0.0; fft_coefficients[4*k +1] = 0.0; fft_coefficients[4*k+2] = 0.0; fft_coefficients[4*k+3] = 0.0;} // gamma
                    else if (k < 70)  {fft_coefficients[4*k] = M_PI*0.025*(k-50)/2.87; fft_coefficients[4*k +1] =M_PI*0.025*(k-50)/2.87; fft_coefficients[4*k+2] = M_PI*0.025*(k-50)/2.87; fft_coefficients[4*k+3] = 0.333333;} // Sigma line
                    else if(k == 70)  {fft_coefficients[4*k] = M_PI*0.5/2.87; fft_coefficients[4*k +1] = M_PI*0.5/2.87; fft_coefficients[4*k+2] = M_PI*0.5/2.87; fft_coefficients[4*k+3] = 0.33333333;} // Rho 
                    else if(k < 90)   {fft_coefficients[4*k] = (M_PI*0.025*(k-70)/2.87)+(M_PI*0.5/2.87); fft_coefficients[4*k +1] = (M_PI*0.025*(k-70)/2.87)+(M_PI*0.5/2.87); fft_coefficients[4*k+2] = (M_PI*0.025*(k-70)/2.87)+(M_PI*0.5/2.87); fft_coefficients[4*k+3] = 0.33333333;} // H 
                }
                fftw_plan_with_nthreads(Nthreads);

                time_range = (sim::total_time)/frequency_step; //dt / dt
                // const double freq_hist_range = ((M_PI / frequency_step))*1e4 - ((M_PI / sim::total_time))*1e4;
                freq_hist_cutoff[0] = 0.0;//((M_PI / sim::total_time))*1e4;
                freq_hist_cutoff[1] = ((M_PI / frequency_step))*1e3;
                 const double freq_hist_range = freq_hist_cutoff[1] - freq_hist_cutoff[0];
               // std::cout <<   \
                " spinwave fft width: " << floor(freq_hist_range) << " (samples), for " << 100.0*time_range/floor(freq_hist_range) << " bins. From " << ((M_PI / sim::total_time))*1e4 << " (THz) to  " \
                                        << ((M_PI / frequency_step))*1e4 << " (THz) " <<  std::endl;
                
                freq_hist_step = 0.1;//freq_hist_range/freq_hist_bins;
                freq_hist_bins = int(freq_hist_range/freq_hist_step);
                
                S_t = (double*) fftw_malloc( sizeof(double) * time_range * N * i_atoms  * 1);
                // S_r = (double*) fftw_malloc( sizeof(double) * time_range * N * i_atoms * 3.0);
                // S_o = (double*) fftw_malloc( sizeof(double) * time_range * N * i_atoms * 3.0);
                r_i = (double*) fftw_malloc( sizeof(double) *i_atoms *N* K_points);
                r_s = (int*) fftw_malloc( sizeof(int) * i_atoms * N);
                // complex K-space
                S_i = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * time_range * N  * i_atoms);
                
                S_0 = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * N * i_atoms);

                // Calculate memory requirements and inform user
                const double mem = double(N) * i_atoms  *time_range*( sizeof(double)  + (1)*sizeof(fftw_complex)) / 1.0e6;
                zlog << zTs() << "Atomistic FFT dipole field calculation has been enabled and requires " << mem << " MB of RAM for " << time_range << " collections" << std::endl;
                std::cout     << "Atomistic FFT dipole field calculation has been enabled and requires " << mem << " MB of RAM for " << time_range << " collections" << std::endl;
                std::cout <<   \
                " spinwave fft width: " << freq_hist_bins << " (samples) at " << freq_hist_step << " (THz). From " << freq_hist_cutoff[0] << " (THz) to  " \
                                        << freq_hist_cutoff[1] << " (THz) " <<  std::endl;
                zlog << zTs() <<   \
                " spinwave fft width: " << freq_hist_bins << " (samples) at " << freq_hist_step << " (THz). From " << freq_hist_cutoff[0] << " (THz) to  " \
                                        << freq_hist_cutoff[1]<< " (THz) " <<  std::endl;
                
            int electron = 0;
           // #pragma omp parallel for reduction(+:lost) num_threads(6) reduction(max:mini)
            for(int i = 0; i < atoms::num_atoms; i++ ) {

                const double x = atoms::x_coord_array[i];
                const double y = atoms::y_coord_array[i];
                const double z = atoms::z_coord_array[i];
             
                std::vector<bool> keep;
                keep.resize(atoms::num_atoms, false);
                if(keep.size() != atoms::num_atoms) std::cout << keep.size() << std::endl;
                int count_r = 1;

                for(int k = 1; k < K_points; k++) {
                    double fft_x = fft_coefficients[4*k];
                    double fft_y = fft_coefficients[4*k+1];
                    double fft_z = fft_coefficients[4*k+2];
                    double normalisation = fft_coefficients[4*k+3];
                    int count_k = 1;
                    for(int nn = 0; nn < atoms::num_atoms; ++nn) {
                        if(nn == i) continue;
                        // const int natom = atoms::neighbour_list_array[nn];
                         double d_x = atoms::x_coord_array[nn] - x;
                         double d_y = atoms::y_coord_array[nn] - y;
                         double d_z = atoms::z_coord_array[nn] - z;

                        if(d_x > cs::system_dimensions[0]*0.5)   d_x -=  cs::system_dimensions[0];//cos(theta)*sin(phi);
                        else if (d_x < -1.0*cs::system_dimensions[0]*0.5) d_x += cs::system_dimensions[0];
                        
                        if(d_y > cs::system_dimensions[1]*0.5)  d_y -= cs::system_dimensions[1];//sin(theta)*sin(phi);
                        else if (d_y < -1.0*cs::system_dimensions[1]*0.5) d_y += cs::system_dimensions[1];
                        
                        if(d_z > cs::system_dimensions[2]*0.5)   d_z -= cs::system_dimensions[2];//cos(phi);
                        else if(d_z < -1.0*cs::system_dimensions[2]*0.5)  d_z += cs::system_dimensions[2];//cos(phi);

                        if(std::abs(d_x) > cs::system_dimensions[0]*cutoff|| \
                           std::abs(d_y) > cs::system_dimensions[1]*cutoff || \
                           std::abs(d_z) > cs::system_dimensions[2]*cutoff) continue;

                        double cutoff =  normalisation*sin(fft_x*std::abs(d_x)/2.0) + normalisation*sin(fft_y*std::abs(d_y)/2.0) + normalisation*sin(fft_z*std::abs(d_z)/2.0);
                        if(cutoff > 0.1) {
                          r_i[electron*i_atoms*(K_points-1) + (k-1)*i_atoms + count_k] = cutoff;
                          count_k++;
                          if(count_k > (i_atoms-2)) { std::cout << count_k << std::endl; break; }
                          if(keep[nn]) continue;
                          r_s[electron*i_atoms + count_r] = nn;
                          keep[nn] = true;
                          count_r++;
                          if(count_r > (i_atoms-2)) { std::cout << count_r << std::endl; break; }
                          
                        }
                    } 
                    r_i[electron*i_atoms*(K_points-1) + (k-1)*i_atoms + 0] = count_k;
                }
                r_s[electron*i_atoms + 0] = count_r;
                electron++;
                if(electron > N-1) std::cout << electron << std::endl;
            }

        //    std::cout << lost  << " out of " << mini << std::endl;

            electron = 0;
            for(int a = 0; a < atoms::num_atoms; a++) {
                    // const int start = atoms::neighbour_list_start_index[a];
   				    // const int end   = atoms::neighbour_list_end_index[a]+1;
                    // if(end - start != 14) std::cout << start << ", " << end << std::endl;
                int count = 1;
                S_0[electron*i_atoms+0][0] = atoms::x_spin_array[a]; //zero-gamma point
                S_0[electron*i_atoms+0][1] = atoms::y_spin_array[a]; //zero-gamma point

                for(int nn = 1; nn < r_s[electron*i_atoms]; ++nn) {
                    const int natom = r_s[electron*i_atoms+nn];
                    S_0[electron*i_atoms+count][0] =  atoms::x_spin_array[natom];
                    S_0[electron*i_atoms+count][1] =  atoms::y_spin_array[natom];
                    count++;
                }
                electron++;
            }

                // create FFTW plans to act on the M and H arrays
                int n[] = {time_range};
                int rank = 1;
                int howmany = N*i_atoms;
                int idist = 1, odist = 1;
                int istride = N*i_atoms, ostride = N*i_atoms;
                int *inembed = n, *onembed = n;

                // Here we plan the transforms, making use of the inter-leaved memory
                // From real space to K-space uses a real to complex
                plan_S = fftw_plan_many_dft_c2r( rank, n, howmany,
                        S_i, inembed, istride, idist,
                        S_t, onembed, ostride, odist,
                        FFTW_MEASURE);

              
                std::cout << "plan_S complete; r time" << std::endl;
                // int r[] = {N*time_range, N*time_range, N*time_range};
                // howmany = i_atoms;
                // rank = 3;
                // idist = 3; odist = 3;
                // istride = i_atoms; ostride = i_atoms;
                // int *inembed_r = r, *onembed_r = r;
                //  plan_r = fftw_plan_many_r2r( rank, r, howmany,
                //         S_r, inembed_r, istride, idist,
                //         S_o, onembed_r, ostride, odist,
                //         FFTW_MEASURE, FFTW_R2HC);
                initialised = true;

}

void spinwave_statistic_t::reset(){

        // int electron = 0;
        #pragma omp parallel for num_threads(2)
        for(int a = 0; a < atoms::num_atoms; a++) {
            // const int start = atoms::neighbour_list_start_index[a];
   		    // const int end   = atoms::neighbour_list_end_index[a]+1;
            
            int count = 0;
            S_0[a*i_atoms+count][0] = atoms::x_spin_array[a];
            S_0[a*i_atoms+count][1] = atoms::y_spin_array[a];
            count++;
            for(int nn = 1; nn < r_s[a*i_atoms]; ++nn) {
                        const int natom = r_s[a*i_atoms+nn];
                S_0[a*i_atoms+count][0] = atoms::x_spin_array[natom];
                S_0[a*i_atoms+count][1] = atoms::y_spin_array[natom];
                count++;
            }
            // electron++;
        }
        time_step = 0;
}
    
void spinwave_statistic_t::update() {

        if(sim::time % frequency_step != 0) return; //update will call with the statistics
                                                    //only want spinwave update with frequency setting

        if(!initialised) {
            std::cout << "FFT spinwave has been called but not initialised." << std::endl;                    exit(-1);
        }
        // int electron = 0;
        #pragma omp parallel for num_threads(2)
        for( int atom = 0; atom < atoms::num_atoms; atom++) {
            // const int start = atoms::neighbour_list_start_index[atom];
   		    // const int end   = atoms::neighbour_list_end_index[atom]+1;
            
            int count = 0;
            spin_correlation(S_i[time_step*N*i_atoms + atom*i_atoms + 0 ], atoms::x_spin_array[atom], atoms::y_spin_array[atom], S_0[atom*i_atoms + 0]);
            count++;
            for(int nn = 1; nn < r_s[atom*i_atoms]; ++nn) {
                        const int natom = r_s[atom*i_atoms+nn];
                spin_correlation(S_i[time_step*N*i_atoms + atom*i_atoms + count], atoms::x_spin_array[natom], atoms::y_spin_array[natom], S_0[atom*i_atoms + count]);
                count++;
            }
            // electron++;
        }

        time_step++;
        return;

    }

            //-----------------------------------------------------------------------------
            // Function to finalize FFT solver and release memory
            //-----------------------------------------------------------------------------
void spinwave_statistic_t::finalize() {

         // instantiate timer
                vutil::vtimer_t timer;

                //   start timer
                timer.start();

        fftw_execute(plan_S);
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
        timer.stop();
        std::cout << "Spinwave compute time = " << timer.elapsed_time() << std::endl;
        zlog << zTs() << "Spinwave compute time = " << timer.elapsed_time()  << std::endl;
        timer.start();

        std::ofstream ofile;
        std::vector< std::vector<int > > fft_frequency_hist;
        // fft_frequency_hist.resize(1);

        fft_frequency_hist.resize(4);
        fft_frequency_hist[0].resize(freq_hist_bins, 0);
        fft_frequency_hist[1].resize(freq_hist_bins, 0);
        fft_frequency_hist[2].resize(freq_hist_bins, 0);
        fft_frequency_hist[3].resize(freq_hist_bins, 0);
        // std::cout << freq_hist_bins << std::endl;
        // for(int k = 0; k < 1; k++) {
            // fft_frequency_hist.resize(freq_hist_bins, 0);
        // }

        ofile.open("atomistic_spinwave_freq");
        const double normalise = 1e4*2.0/sim::total_time/M_PI;//1e4*1.0/sqrt(1.0*M_PI);
       
        std::vector<std::vector< double> > fft_integration;
        fft_integration.resize(4);
        fft_integration[0].resize(K_points+1, 0.0);
        fft_integration[1].resize(K_points+1, 0.0);
        fft_integration[2].resize(K_points+1, 0.0);
        fft_integration[3].resize(K_points+1, 0.0);
    //  std::cout << freq_hist_bins << std::endl;

    #pragma omp parallel for num_threads(4)
    for(int t = 0; t < time_range; t++) {
        for(int a = 0; a < N; a++) {
            for(int k = 0; k < r_s[a*i_atoms]-1; k++) {
                const double value = std::abs(S_t[t*N*i_atoms + a*i_atoms + k]*normalise);
                for(int r = 1; r < K_points+1; r++) {
                    if ( k == 0 ) fft_integration[omp_get_thread_num()][0] += value;
                    else if (k < r_i[a*(K_points-1)*(i_atoms) +(r-1)*i_atoms + 0]) fft_integration[omp_get_thread_num()][r] += value*r_i[a*(K_points-1)*(i_atoms) +(r-1)*i_atoms +(k-1)];
                } 
                 int h_value = std::min(freq_hist_bins-1, int(std::max(0.0, floor((value-freq_hist_cutoff[0])/freq_hist_step))));
                // if(h_value > 108) {
                //     std::cout << h_value << ", " << value << ", " << floor((value-freq_hist_cutoff)/freq_hist_step) << ", " << int(std::max(0.0, floor((value-freq_hist_cutoff)/freq_hist_step))) << std::endl;
                //     h_value = 108;
                // } else if (h_value < 0) {
                //     std::cout << h_value << ", " << value << ", " << floor((value-freq_hist_cutoff)/freq_hist_step) << ", " << int(std::max(0.0, floor((value-freq_hist_cutoff)/freq_hist_step))) << std::endl;
                //     h_value = 0;
                // }
                fft_frequency_hist[omp_get_thread_num()][h_value]++;
            }
        }
    }

    for(int f = 0; f < freq_hist_bins; f++) {
        ofile <<  f*freq_hist_step + freq_hist_cutoff[0] << "\t" << double(fft_frequency_hist[0][f]+fft_frequency_hist[1][f]+fft_frequency_hist[2][f]+fft_frequency_hist[3][f])/i_atoms << "\n";
    }


        // ofile << S_t[t*N*i_atoms + a + o]*normalise << "\t";   
        ofile.close();//("/spinwave/atomistic_spinwave" + std::to_string(o) + ".txt");

    timer.stop();
        std::cout << "Spinwave freq output time = " << timer.elapsed_time() << std::endl;
                zlog << zTs() << "Spinwave freq output time = " << timer.elapsed_time()  << std::endl;
        timer.start();

          ofile.open("atomistic_spinwave_spatial");
          
    for(int r = 0; r < K_points+1; r++) {
        ofile << r << "\t" << (fft_integration[0][r]+fft_integration[1][r]+fft_integration[2][r]+fft_integration[3][r])/N/(i_atoms-1) << "\n";
    } 
    ofile.close();

         timer.stop();
                
       
        std::cout << "Spinwave spatial output time = " << timer.elapsed_time() << std::endl;
                zlog << zTs() << "Spinwave spatial output time = " << timer.elapsed_time()  << std::endl;
        // timer.start();
           
                // Free memory from FFT complex variables
                fftw_free(S_i);
                // fftw_free(r_i);
                fftw_free(S_t);
                // fftw_free(S_r);
                fftw_free(S_0);
                fftw_free(r_s);
                fftw_free(r_i);
                fftw_free(fft_coefficients);
                fftw_destroy_plan(plan_S);
                //  fftw_destroy_plan(plan_r);

                fftw_cleanup_threads();

                return;

    }
}
