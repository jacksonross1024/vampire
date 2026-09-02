//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans, Daniel Meilak 2017-2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
// FFT magnetostatic stray field in a vacuum plane above the sample.
//
// Spins are binned with the existing vdc cell algorithm (cell magnetisation,
// not atomistic coordinates). Each occupied z-slice is a sheet of point
// dipoles at cell centres. The in-plane field map at
//
//    z_obs = z_sample_top + --sf-height
//
// is the discrete dipole sum
//
//    B(r) = (μ0/4π) Σ_src [ 3(m·R)R - m R² ] / R⁵
//
// evaluated by 2D FFT convolution of that kernel. Default --sf-pad 2
// zero-pads the in-plane grid so a finite flake is not tiled (pad 1 is
// periodic). Output is B in tesla. Coordinates in the file are Angstroms
// (vdc default); --sf-height is in nanometres (default 30 nm = 300 Å).
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// program header
#include "vdc.hpp"

#ifdef VDC_FFTW
#include <fftw3.h>
#endif

#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

namespace vdc{

#ifndef VDC_FFTW

namespace {

void fftw_missing_message(){
   std::cerr << "Error: --stray-field requires the FFTW3 library.\n"
             << "This vdc binary was compiled without FFTW3.\n"
             << "Install FFTW3 development files (Debian/Ubuntu: sudo apt install libfftw3-dev;\n"
             << "macOS: brew install fftw) and recompile vdc:\n"
             << "  cd util/vdc && make clean && make FFTW=1\n"
             << "If the library is already installed, `make` auto-detects it; use FFTW=1 to force\n"
             << "a build with FFTW, or FFTW=0 to build without the stray-field solver.\n";
}

} // end of anonymous namespace

//------------------------------------------------------------------------------
void require_stray_field_support(){
   fftw_missing_message();
   std::exit(EXIT_FAILURE);
}

//------------------------------------------------------------------------------
void output_stray_field_file(unsigned int spin_file_id){
   (void)spin_file_id;
   fftw_missing_message();
   std::exit(EXIT_FAILURE);
}

#else

namespace {

const double MU0_4PI = 1.0e-7;                 // T m / A  (μ0/4π)
const double MU_B = 9.2740100783e-24;          // J/T = A m^2
const double ANGSTROM_TO_M = 1.0e-10;
const double NM_TO_ANGSTROM = 10.0;

void cplx_fma(double* acc, const double* n, const double* m){
   acc[0] += n[0] * m[0] - n[1] * m[1];
   acc[1] += n[0] * m[1] + n[1] * m[0];
}

int fft_index(int i, int n){
   return (i < (n + 1) / 2) ? i : i - n;
}

struct stray_fft_workspace {
   int nx;
   int ny;
   fftw_plan fwd;
   fftw_plan bwd;
   fftw_complex* px;
   fftw_complex* py;
   fftw_complex* pz;
   fftw_complex* scratch;
   fftw_complex* Hx;
   fftw_complex* Hy;
   fftw_complex* Hz;
};

stray_fft_workspace g_fft = {
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

void free_fft_workspace(){

   if(g_fft.fwd) fftw_destroy_plan(g_fft.fwd);
   if(g_fft.bwd) fftw_destroy_plan(g_fft.bwd);
   if(g_fft.px) fftw_free(g_fft.px);
   if(g_fft.py) fftw_free(g_fft.py);
   if(g_fft.pz) fftw_free(g_fft.pz);
   if(g_fft.scratch) fftw_free(g_fft.scratch);
   if(g_fft.Hx) fftw_free(g_fft.Hx);
   if(g_fft.Hy) fftw_free(g_fft.Hy);
   if(g_fft.Hz) fftw_free(g_fft.Hz);
   g_fft.fwd = 0;
   g_fft.bwd = 0;
   g_fft.px = g_fft.py = g_fft.pz = g_fft.scratch = 0;
   g_fft.Hx = g_fft.Hy = g_fft.Hz = 0;
   g_fft.nx = 0;
   g_fft.ny = 0;

}

void ensure_fft_workspace(int nx, int ny){

   if(g_fft.px != 0 && g_fft.nx == nx && g_fft.ny == ny) return;

   free_fft_workspace();

   const size_t nxy = static_cast<size_t>(nx) * static_cast<size_t>(ny);
   g_fft.nx = nx;
   g_fft.ny = ny;
   g_fft.px = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * nxy));
   g_fft.py = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * nxy));
   g_fft.pz = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * nxy));
   g_fft.scratch = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * nxy));
   g_fft.Hx = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * nxy));
   g_fft.Hy = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * nxy));
   g_fft.Hz = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * nxy));

   if(!g_fft.px || !g_fft.py || !g_fft.pz || !g_fft.scratch || !g_fft.Hx || !g_fft.Hy || !g_fft.Hz){
      std::cerr << "Error - FFTW malloc failed for stray-field grid " << nx << " x " << ny << std::endl;
      std::exit(EXIT_FAILURE);
   }

   g_fft.fwd = fftw_plan_dft_2d(nx, ny, g_fft.px, g_fft.px, FFTW_FORWARD, FFTW_ESTIMATE);
   g_fft.bwd = fftw_plan_dft_2d(nx, ny, g_fft.Hx, g_fft.Hx, FFTW_BACKWARD, FFTW_ESTIMATE);

}

void zero_cplx(fftw_complex* a, size_t n){
   for(size_t i = 0; i < n; i++){
      a[i][0] = 0.0;
      a[i][1] = 0.0;
   }
}

struct stray_binary_ctx {
   const std::vector<double>* x;
   const std::vector<double>* y;
   const std::vector<double>* bx;
   const std::vector<double>* by;
   const std::vector<double>* bz;
   double z_obs;
};

void fill_stray_binary_row(size_t i, double* out, void* vctx){

   const stray_binary_ctx* ctx = static_cast<stray_binary_ctx*>(vctx);
   out[0] = (*ctx->x)[i];
   out[1] = (*ctx->y)[i];
   out[2] = ctx->z_obs;
   out[3] = (*ctx->bx)[i];
   out[4] = (*ctx->by)[i];
   out[5] = (*ctx->bz)[i];

}

} // end of anonymous namespace

//------------------------------------------------------------------------------
void require_stray_field_support(){
}

//------------------------------------------------------------------------------
void output_stray_field_file(unsigned int spin_file_id){

   if(vdc::total_cells == 0){
      std::cerr << "Error - stray-field FFT: no occupied cells." << std::endl;
      std::exit(EXIT_FAILURE);
   }
   if(vdc::cell_ijk.size() != 3ull * vdc::total_cells){
      std::cerr << "Error - stray-field FFT: cell_ijk size does not match occupied cells." << std::endl;
      std::exit(EXIT_FAILURE);
   }

   // FFT grid is the occupied ijk span (not 1+ceil nx_cells, which can add an empty strip)
   int nx = 0, ny = 0, nz = 0;
   for(unsigned int cell = 0; cell < vdc::total_cells; cell++){
      const int ix = vdc::cell_ijk[3 * cell + 0];
      const int iy = vdc::cell_ijk[3 * cell + 1];
      const int iz = vdc::cell_ijk[3 * cell + 2];
      if(ix + 1 > nx) nx = ix + 1;
      if(iy + 1 > ny) ny = iy + 1;
      if(iz + 1 > nz) nz = iz + 1;
   }

   if(nx < 2 || ny < 2){
      std::cerr << "Error - stray-field FFT needs at least a 2x2 in-plane cell grid "
                << "(got " << nx << " x " << ny << "). Reduce cell-size or use a larger sample."
                << std::endl;
      std::exit(EXIT_FAILURE);
   }
   if(nz < 1){
      std::cerr << "Error - stray-field FFT: no z cells." << std::endl;
      std::exit(EXIT_FAILURE);
   }

   const double dx_A = vdc::cell_size[0];
   const double dy_A = vdc::cell_size[1];
   const double dz_A = vdc::cell_size[2];
   const double dx_m = dx_A * ANGSTROM_TO_M;
   const double dy_m = dy_A * ANGSTROM_TO_M;
   const double area_m2 = dx_m * dy_m;

   if(area_m2 <= 0.0){
      std::cerr << "Error - stray-field FFT: cell in-plane area must be positive." << std::endl;
      std::exit(EXIT_FAILURE);
   }

   // Top face of the highest occupied cell (Angstrom)
   double z_top_A = vdc::cell_origin[2];
   for(unsigned int cell = 0; cell < vdc::total_cells; cell++){
      const double z_face = vdc::cell_coords[3 * cell + 2] + dz_A;
      if(z_face > z_top_A) z_top_A = z_face;
   }

   const double z_obs_A = z_top_A + vdc::sf_height_nm * NM_TO_ANGSTROM;
   const double z_obs_m = z_obs_A * ANGSTROM_TO_M;

   const int pad = vdc::sf_pad;
   if(pad < 1 || nx > 100000000 / pad || ny > 100000000 / pad){
      std::cerr << "Error - stray-field FFT: --sf-pad " << pad
                << " is too large for grid " << nx << " x " << ny << "." << std::endl;
      std::exit(EXIT_FAILURE);
   }
   const int nxfft = nx * pad;
   const int nyfft = ny * pad;
   const size_t nxy = static_cast<size_t>(nx) * static_cast<size_t>(ny);
   const size_t nxy_fft = static_cast<size_t>(nxfft) * static_cast<size_t>(nyfft);
   ensure_fft_workspace(nxfft, nyfft);

   zero_cplx(g_fft.Hx, nxy_fft);
   zero_cplx(g_fft.Hy, nxy_fft);
   zero_cplx(g_fft.Hz, nxy_fft);

   const unsigned int tmid = 1u + static_cast<unsigned int>(vdc::materials.size());

   auto fft_id = [&](int ix, int iy) -> size_t {
      return static_cast<size_t>(ix) * static_cast<size_t>(nyfft) + static_cast<size_t>(iy);
   };

   for(int iz = 0; iz < nz; iz++){

      zero_cplx(g_fft.px, nxy_fft);
      zero_cplx(g_fft.py, nxy_fft);
      zero_cplx(g_fft.pz, nxy_fft);

      bool layer_used = false;
      for(unsigned int cell = 0; cell < vdc::total_cells; cell++){
         if(vdc::cell_ijk[3 * cell + 2] != iz) continue;

         const int ix = vdc::cell_ijk[3 * cell + 0];
         const int iy = vdc::cell_ijk[3 * cell + 1];
         if(ix < 0 || iy < 0 || ix >= nx || iy >= ny) continue;

         // Stored (mx,my,mz) already has length |m|. Moment (μ_B) = n * μ_spin * m.
         double mux = 0.0, muy = 0.0, muz = 0.0;
         for(unsigned int m = 0; m + 1 < tmid; m++){
            const double nmu = static_cast<double>(vdc::num_atoms_in_cell[cell * tmid + m])
                             * vdc::materials[m].moment;
            mux += vdc::cell_magnetization[cell][m][0] * nmu;
            muy += vdc::cell_magnetization[cell][m][1] * nmu;
            muz += vdc::cell_magnetization[cell][m][2] * nmu;
         }

         const size_t id = fft_id(ix, iy);
         g_fft.px[id][0] += mux * MU_B;
         g_fft.py[id][0] += muy * MU_B;
         g_fft.pz[id][0] += muz * MU_B;
         layer_used = true;
      }
      if(!layer_used) continue;

      const double z_layer_A = vdc::cell_origin[2] + (static_cast<double>(iz) + 0.5) * dz_A;
      const double dz = z_obs_m - z_layer_A * ANGSTROM_TO_M;

      fftw_execute_dft(g_fft.fwd, g_fft.px, g_fft.px);
      fftw_execute_dft(g_fft.fwd, g_fft.py, g_fft.py);
      fftw_execute_dft(g_fft.fwd, g_fft.pz, g_fft.pz);

      auto apply_kernel = [&](int which){
         // which: 0 Nxx, 1 Nxy, 2 Nxz, 3 Nyy, 4 Nyz, 5 Nzz
         zero_cplx(g_fft.scratch, nxy_fft);
         #pragma omp parallel for
         for(int ix = 0; ix < nxfft; ix++){
            const double x = static_cast<double>(fft_index(ix, nxfft)) * dx_m;
            for(int iy = 0; iy < nyfft; iy++){
               const double y = static_cast<double>(fft_index(iy, nyfft)) * dy_m;
               const double r2 = x * x + y * y + dz * dz;
               const size_t id = fft_id(ix, iy);
               if(r2 < 1.0e-60) continue;
               const double r = std::sqrt(r2);
               const double r5 = r2 * r2 * r;
               const double pref = MU0_4PI / r5;
               double val = 0.0;
               switch(which){
                  case 0: val = pref * (3.0 * x * x - r2); break;
                  case 1: val = pref * (3.0 * x * y); break;
                  case 2: val = pref * (3.0 * x * dz); break;
                  case 3: val = pref * (3.0 * y * y - r2); break;
                  case 4: val = pref * (3.0 * y * dz); break;
                  default: val = pref * (3.0 * dz * dz - r2); break;
               }
               g_fft.scratch[id][0] = val;
            }
         }
         fftw_execute_dft(g_fft.fwd, g_fft.scratch, g_fft.scratch);
         for(size_t id = 0; id < nxy_fft; id++){
            const double* n = g_fft.scratch[id];
            switch(which){
               case 0:
                  cplx_fma(g_fft.Hx[id], n, g_fft.px[id]);
                  break;
               case 1:
                  cplx_fma(g_fft.Hx[id], n, g_fft.py[id]);
                  cplx_fma(g_fft.Hy[id], n, g_fft.px[id]);
                  break;
               case 2:
                  cplx_fma(g_fft.Hx[id], n, g_fft.pz[id]);
                  cplx_fma(g_fft.Hz[id], n, g_fft.px[id]);
                  break;
               case 3:
                  cplx_fma(g_fft.Hy[id], n, g_fft.py[id]);
                  break;
               case 4:
                  cplx_fma(g_fft.Hy[id], n, g_fft.pz[id]);
                  cplx_fma(g_fft.Hz[id], n, g_fft.py[id]);
                  break;
               default:
                  cplx_fma(g_fft.Hz[id], n, g_fft.pz[id]);
                  break;
            }
         }
      };

      for(int which = 0; which < 6; which++) apply_kernel(which);
   }

   fftw_execute_dft(g_fft.bwd, g_fft.Hx, g_fft.Hx);
   fftw_execute_dft(g_fft.bwd, g_fft.Hy, g_fft.Hy);
   fftw_execute_dft(g_fft.bwd, g_fft.Hz, g_fft.Hz);

   const double inv_n = 1.0 / static_cast<double>(nxy_fft);
   std::vector<double> bx(nxy), by(nxy), bz(nxy), xA(nxy), yA(nxy);

   double bxmin = 1.0e300, bxmax = -1.0e300;
   double bymin = 1.0e300, bymax = -1.0e300;
   double bzmin = 1.0e300, bzmax = -1.0e300;

   for(int ix = 0; ix < nx; ix++){
      const double x = vdc::cell_origin[0] + (static_cast<double>(ix) + 0.5) * dx_A;
      for(int iy = 0; iy < ny; iy++){
         const size_t id = static_cast<size_t>(ix) * static_cast<size_t>(ny) + static_cast<size_t>(iy);
         const size_t idf = fft_id(ix, iy);
         const double y = vdc::cell_origin[1] + (static_cast<double>(iy) + 0.5) * dy_A;
         const double Bx = g_fft.Hx[idf][0] * inv_n;
         const double By = g_fft.Hy[idf][0] * inv_n;
         const double Bz = g_fft.Hz[idf][0] * inv_n;
         xA[id] = x;
         yA[id] = y;
         bx[id] = Bx;
         by[id] = By;
         bz[id] = Bz;
         if(Bx < bxmin) bxmin = Bx; if(Bx > bxmax) bxmax = Bx;
         if(By < bymin) bymin = By; if(By > bymax) bymax = By;
         if(Bz < bzmin) bzmin = Bz; if(Bz > bzmax) bzmax = Bz;
      }
   }

   std::stringstream sf_file_sstr;
   sf_file_sstr << "stray-field-";
   sf_file_sstr << std::setfill('0') << std::setw(8) << spin_file_id;
   sf_file_sstr << ".txt";
   const std::string sf_file_name = sf_file_sstr.str();

   if(vdc::binary_dump){

      const std::string bin_name = vdc::binary_filename(sf_file_name);
      std::cout << "   Writing stray-field file " << bin_name << "..." << std::flush;

      stray_binary_ctx ctx;
      ctx.x = &xA;
      ctx.y = &yA;
      ctx.bx = &bx;
      ctx.by = &by;
      ctx.bz = &bz;
      ctx.z_obs = z_obs_A;

      binary_dump_spec spec;
      spec.kind = "stray-field";
      spec.notes = "2D FFT discrete-dipole convolution in a vacuum plane above the sample. One row per in-plane cell on the occupied nx*ny grid. Dipole moment per cell is n_spins * mu_spin * stored (mx,my,mz) in mu_B (stored vector already has length |m|; |m| is not applied again), converted to A m^2. Each z-slice is a sheet of point dipoles at cell centres; B = (mu0/4pi) Sum [3(m.R)R - m R^2]/R^5. In-plane zero-padding (sf_pad, default 2) avoids periodic tiling of a finite flake; pad 1 is periodic. Coordinates are lab-frame cell centres in Angstroms. B is vacuum flux density in tesla.";
      spec.columns = {
         {"x", "Angstrom", "in-plane cell centre x"},
         {"y", "Angstrom", "in-plane cell centre y"},
         {"z", "Angstrom", "observation plane z_obs"},
         {"Bx", "T", "stray flux density x"},
         {"By", "T", "stray flux density y"},
         {"Bz", "T", "stray flux density z"}
      };

      vdc::output_binary_file(bin_name, static_cast<uint64_t>(nxy), 6, fill_stray_binary_row, &ctx, spec);
      std::cout << "done!" << std::endl;

   }
   else{

      std::cout << "   Writing stray-field file " << sf_file_name << "..." << std::flush;

      std::ofstream ofile;
      ofile.open(sf_file_name.c_str());
      ofile << std::setprecision(10);
      ofile << "# VDC stray field (2D FFT discrete-dipole convolution)\n";
      ofile << "# height_above_sample_nm = " << vdc::sf_height_nm << "\n";
      ofile << "# z_sample_top_Angstrom = " << z_top_A << "\n";
      ofile << "# z_obs_Angstrom = " << z_obs_A << "\n";
      ofile << "# cell_size_Angstrom = " << dx_A << " " << dy_A << " " << dz_A << "\n";
      ofile << "# grid_nx ny nz = " << nx << " " << ny << " " << nz << "\n";
      ofile << "# sf_pad = " << pad << " (1 = periodic, >=2 = zero-padded finite flake)\n";
      ofile << "# fft_grid_nx ny = " << nxfft << " " << nyfft << "\n";
      ofile << "# B = (mu0/4pi) Sum [3(m.R)R - m R^2]/R^5 at cell centres; tesla. x,y,z in Angstrom. Blank line after each y-row for gnuplot.\n";
      ofile << "#x\ty\tz\tBx\tBy\tBz\n";

      for(int iy = 0; iy < ny; iy++){
         for(int ix = 0; ix < nx; ix++){
            const size_t id = static_cast<size_t>(ix) * static_cast<size_t>(ny) + static_cast<size_t>(iy);
            ofile << xA[id] << "\t" << yA[id] << "\t" << z_obs_A << "\t"
                  << bx[id] << "\t" << by[id] << "\t" << bz[id] << "\n";
         }
         ofile << "\n";
      }
      ofile.close();
      std::cout << "done!" << std::endl;
   }

   if(vdc::verbose){
      std::cout << "   stray-field grid " << nx << " x " << ny << " x " << nz
                << " (occupied ijk span), fft pad " << pad << " -> " << nxfft << " x " << nyfft
                << ", z_obs = " << z_obs_A << " A ("
                << vdc::sf_height_nm << " nm above sample top)\n";
      std::cout << "   B range (T): x [" << bxmin << ", " << bxmax
                << "], y [" << bymin << ", " << bymax
                << "], z [" << bzmin << ", " << bzmax << "]\n";
   }

}

#endif // VDC_FFTW

} // end of namespace vdc
