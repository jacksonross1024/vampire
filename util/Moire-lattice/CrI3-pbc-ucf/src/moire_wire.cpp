#include "moire_wire.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>

#include "initialise.hpp"
#include "moire_coarse.hpp"
#include "positions.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

std::vector<spin_interaction> g_moire_spin_interactions;

void moire_spin_interactions_clear() { g_moire_spin_interactions.clear(); }

void moire_spin_interactions_merge(std::vector<spin_interaction>& thread_chunk) {
   g_moire_spin_interactions.insert(
       g_moire_spin_interactions.end(), thread_chunk.begin(), thread_chunk.end());
   thread_chunk.clear();
}

extern double twist_angle;

namespace {

constexpr double kSpinMagnitude = 1.5;
constexpr double kBilayerCellThicknessM = 6.6e-10;
constexpr double kDoubleBilayerLayerThicknessM = 6.54e-10;
constexpr double kInterfaceSpacingM = 6.6e-10;
constexpr int32_t kInterfaceCoordination = 1;
constexpr int32_t kAtomsPerInterfaceCell = 2;
constexpr int32_t kTwistedDoubleBilayerNz = 3;
constexpr int32_t kTwistedDoubleBilayerMiddleCell = 1;
constexpr int32_t kTwistedDoubleBilayerMiddleSub1Layer = 3;
constexpr int32_t kTwistedDoubleBilayerMiddleSub2Layer = 2;

void write_binary_v3(
    const std::string& path,
    int nx_m,
    int ny_m,
    double moire_Lx_angstrom,
    double moire_Ly_angstrom,
    const MoireAverages& avg) {
   const char magicb[4] = {'M', 'O', 'M', 'X'};
   const uint32_t ver = 3;
   const int ncell = nx_m * ny_m;
   std::ofstream f(path.c_str(), std::ios::binary);
   if (!f) {
      std::cerr << "moire_wire: could not open " << path << std::endl;
      return;
   }
   f.write(magicb, 4);
   f.write(reinterpret_cast<const char*>(&ver), 4);
   uint32_t u32 = static_cast<uint32_t>(nx_m);
   f.write(reinterpret_cast<const char*>(&u32), 4);
   u32 = static_cast<uint32_t>(ny_m);
   f.write(reinterpret_cast<const char*>(&u32), 4);
   f.write(reinterpret_cast<const char*>(&moire_Lx_angstrom), sizeof(double));
   f.write(reinterpret_cast<const char*>(&moire_Ly_angstrom), sizeof(double));
   double twist_deg = twist_angle * 180.0 / M_PI;
   f.write(reinterpret_cast<const char*>(&twist_deg), sizeof(double));
   const double cell_thickness_m = kBilayerCellThicknessM;
   const double S = kSpinMagnitude;
   const double d_inter_m = kInterfaceSpacingM;
   const double c_mono_m = kBilayerCellThicknessM;
   int32_t z_inter = 1;
   int32_t n_atoms_cell = kAtomsPerInterfaceCell;
   f.write(reinterpret_cast<const char*>(&cell_thickness_m), sizeof(double));
   f.write(reinterpret_cast<const char*>(&S), sizeof(double));
   f.write(reinterpret_cast<const char*>(&d_inter_m), sizeof(double));
   f.write(reinterpret_cast<const char*>(&c_mono_m), sizeof(double));
   f.write(reinterpret_cast<const char*>(&z_inter), sizeof(int32_t));
   f.write(reinterpret_cast<const char*>(&n_atoms_cell), sizeof(int32_t));
   uint64_t ncells_u = static_cast<uint64_t>(ncell);
   f.write(reinterpret_cast<const char*>(&ncells_u), 8);
   uint64_t nbonds = g_moire_spin_interactions.size();
   f.write(reinterpret_cast<const char*>(&nbonds), 8);

   auto write_vec = [&f](const std::vector<double>& v) {
      f.write(reinterpret_cast<const char*>(v.data()),
              static_cast<std::streamsize>(v.size() * sizeof(double)));
   };
   write_vec(avg.J_cell);
   write_vec(avg.J_nn);
   write_vec(avg.J_homo_mean);
   write_vec(avg.J_homo_Jr2);
   write_vec(avg.D_cell_vec);
   write_vec(avg.M_homo);
   write_vec(avg.M_inhomo);
   f.close();
   std::cout << "moire_wire: wrote " << path << " (" << nbonds
             << " bonds → " << nx_m << " x " << ny_m << " cells)" << std::endl;
}

void write_twisted_double_bilayer_binary_v1(
    const std::string& path,
    int nx_m,
    int ny_m,
    double moire_Lx_angstrom,
    double moire_Ly_angstrom,
    const TwistedDoubleBilayerMoireAverages& avg) {
   const char magicb[4] = {'M', 'O', 'T', 'D'};
   const uint32_t ver = 2;
   const int ncell = nx_m * ny_m;
   std::ofstream f(path.c_str(), std::ios::binary);
   if (!f) {
      std::cerr << "moire_wire: could not open " << path << std::endl;
      return;
   }

   const double twist_deg = twist_angle * 180.0 / M_PI;
   const double mumax_cell_thickness_m =
       (4.0 * kDoubleBilayerLayerThicknessM) / 3.0;
   const double layer_thickness_m[4] = {
       kDoubleBilayerLayerThicknessM, kDoubleBilayerLayerThicknessM,
       kDoubleBilayerLayerThicknessM, kDoubleBilayerLayerThicknessM};
   const double interface_spacing_m[3] = {
       kInterfaceSpacingM, kInterfaceSpacingM, kInterfaceSpacingM};
   const int32_t interface_coordination[3] = {
       kInterfaceCoordination, kInterfaceCoordination, kInterfaceCoordination};
   const int32_t atoms_per_interface_cell[3] = {
       kAtomsPerInterfaceCell, kAtomsPerInterfaceCell, kAtomsPerInterfaceCell};
   const int32_t layer_to_cell[4] = {0, 1, 1, 2};
   const int32_t nz_cells = kTwistedDoubleBilayerNz;
   const int32_t middle_cell = kTwistedDoubleBilayerMiddleCell;
   const uint64_t ncells_u = static_cast<uint64_t>(ncell);
   const uint64_t nbonds = g_moire_spin_interactions.size();

   f.write(magicb, 4);
   f.write(reinterpret_cast<const char*>(&ver), 4);
   uint32_t u32 = static_cast<uint32_t>(nx_m);
   f.write(reinterpret_cast<const char*>(&u32), 4);
   u32 = static_cast<uint32_t>(ny_m);
   f.write(reinterpret_cast<const char*>(&u32), 4);
   f.write(reinterpret_cast<const char*>(&moire_Lx_angstrom), sizeof(double));
   f.write(reinterpret_cast<const char*>(&moire_Ly_angstrom), sizeof(double));
   f.write(reinterpret_cast<const char*>(&twist_deg), sizeof(double));
   f.write(reinterpret_cast<const char*>(&kSpinMagnitude), sizeof(double));
   f.write(reinterpret_cast<const char*>(&mumax_cell_thickness_m), sizeof(double));
   f.write(reinterpret_cast<const char*>(layer_thickness_m), sizeof(layer_thickness_m));
   f.write(reinterpret_cast<const char*>(interface_spacing_m),
           sizeof(interface_spacing_m));
   f.write(reinterpret_cast<const char*>(interface_coordination),
           sizeof(interface_coordination));
   f.write(reinterpret_cast<const char*>(atoms_per_interface_cell),
           sizeof(atoms_per_interface_cell));
   f.write(reinterpret_cast<const char*>(layer_to_cell), sizeof(layer_to_cell));
   f.write(reinterpret_cast<const char*>(&nz_cells), sizeof(int32_t));
   f.write(reinterpret_cast<const char*>(&middle_cell), sizeof(int32_t));
   f.write(reinterpret_cast<const char*>(&ncells_u), sizeof(uint64_t));
   f.write(reinterpret_cast<const char*>(&nbonds), sizeof(uint64_t));

   auto write_vec = [&f](const std::vector<double>& v) {
      f.write(reinterpret_cast<const char*>(v.data()),
              static_cast<std::streamsize>(v.size() * sizeof(double)));
   };
   write_vec(avg.J_23_cell);
   write_vec(avg.J_12_nn);
   write_vec(avg.J_34_nn);
   write_vec(avg.J_layer1_homo_mean);
   write_vec(avg.J_layer1_homo_Jr2);
   write_vec(avg.J_layer2_homo_mean);
   write_vec(avg.J_layer2_homo_Jr2);
   write_vec(avg.J_layer3_homo_mean);
   write_vec(avg.J_layer3_homo_Jr2);
   write_vec(avg.J_layer4_homo_mean);
   write_vec(avg.J_layer4_homo_Jr2);
   write_vec(avg.D_23_cell_vec);
   write_vec(avg.M_inhomo_12);
   write_vec(avg.M_inhomo_34);
   write_vec(avg.M_homo_layer1);
   write_vec(avg.M_homo_layer2);
   write_vec(avg.M_homo_layer3);
   write_vec(avg.M_homo_layer4);
   f.close();
   std::cout << "moire_wire: wrote " << path << " (" << nbonds
             << " bonds → " << nx_m << " x " << ny_m
             << " cells, twisted-double-bilayer MOTD v2 + DMI)" << std::endl;
}

void write_twisted_bilayer_meta_txt(const std::string& path, int nx_m, int ny_m,
                                    double Lx_A, double Ly_A, double bin_x_A,
                                    double bin_y_A, double nn_tol_A,
                                    double nn_scale) {
   std::ofstream m(path.c_str());
   if (!m) return;
   m << "nx_moire " << nx_m << "\n";
   m << "ny_moire " << ny_m << "\n";
   m << "moire_Lx_angstrom " << Lx_A << "\n";
   m << "moire_Ly_angstrom " << Ly_A << "\n";
   m << "micromagnetic_bin_ax_angstrom " << bin_x_A << "\n";
   m << "micromagnetic_bin_ay_angstrom " << bin_y_A << "\n";
   m << "nn_tol_angstrom " << nn_tol_A << "\n";
   m << "moire_coarse_nn_tol_scale " << nn_scale << "\n";
   m << "binary moire_coarse_v2.bin format_version 3\n";
   m << "python: python3 twisted_bilayer_moire_coarse_bin_to_npz.py moire_coarse_v2.bin\n";
   m.close();
}

void write_twisted_double_bilayer_meta_txt(
    const std::string& path, int nx_m, int ny_m, double Lx_A, double Ly_A,
    double bin_x_A, double bin_y_A, double nn_tol_A, double nn_scale) {
   std::ofstream m(path.c_str());
   if (!m) return;
   m << "nx_moire " << nx_m << "\n";
   m << "ny_moire " << ny_m << "\n";
   m << "moire_Lx_angstrom " << Lx_A << "\n";
   m << "moire_Ly_angstrom " << Ly_A << "\n";
   m << "micromagnetic_bin_ax_angstrom " << bin_x_A << "\n";
   m << "micromagnetic_bin_ay_angstrom " << bin_y_A << "\n";
   m << "nn_tol_angstrom " << nn_tol_A << "\n";
   m << "moire_coarse_nn_tol_scale " << nn_scale << "\n";
   m << "binary moire_twisted_double_bilayer_coarse.bin magic MOTD format_version 2\n";
   m << "stack_nz " << kTwistedDoubleBilayerNz << "\n";
   m << "middle_cell_index " << kTwistedDoubleBilayerMiddleCell << "\n";
   m << "middle_sub1_layer " << kTwistedDoubleBilayerMiddleSub1Layer << "\n";
   m << "middle_sub2_layer " << kTwistedDoubleBilayerMiddleSub2Layer << "\n";
   m << "layer_to_cell 0 1 1 2\n";
   m << "python: python3 twisted_double_bilayer_moire_coarse_bin_to_npz.py "
        "moire_twisted_double_bilayer_coarse.bin\n";
   m.close();
}

}  // namespace

void moire_spin_interactions_finalize_and_write(const char* cwd,
                                                const double moire_phys[4]) {
   const double Lx = moire_phys[2];
   const double Ly = moire_phys[3];
   // Micromagnetic moiré bins: cell size 2*a0x by 2*a1y (not microcell_Nx/Ny).
   const double bin_x = 2.0 * a0x;
   const double bin_y = 2.0 * a1y;
   if (bin_x <= 0.0 || bin_y <= 0.0 || Lx <= 0.0 || Ly <= 0.0) {
      std::cerr << "moire_wire: skip coarse grain (bad period or lattice constants)\n";
      return;
   }
   const int nx_m = std::max(1, static_cast<int>(std::floor(Lx / bin_x + 1e-9)));
   const int ny_m = std::max(1, static_cast<int>(std::floor(Ly / bin_y + 1e-9)));
   const double cell_wx = Lx / static_cast<double>(nx_m);
   const double cell_wy = Ly / static_cast<double>(ny_m);
   const double nn_tol =
       moire_coarse_nn_tol_scale * std::min(cell_wx, cell_wy);
   std::cout << "moire coarse: nx_m=" << nx_m << " ny_m=" << ny_m
             << " (target bin " << bin_x << " x " << bin_y
             << " Å), cell " << cell_wx << " x " << cell_wy << " Å, nn_tol="
             << nn_tol << " Å (scale=" << moire_coarse_nn_tol_scale
             << " * min(cell_wx,cell_wy))\n";
   const double interlayer_tol = 1.5;
   const double cell_vol_SI =
       (Lx * 1e-10) * (Ly * 1e-10) * 6.6e-10 /
       (static_cast<double>(nx_m) * static_cast<double>(ny_m));

   const std::string dir(cwd);

   if (moire_coarse_write_twisted_bilayer) {
      MoireAverages avg = coarseGrainMoireBonds(
          g_moire_spin_interactions, nx_m, ny_m, Lx, Ly, cell_vol_SI, nn_tol,
          interlayer_tol);
      write_binary_v3(dir + "/moire_coarse_v2.bin", nx_m, ny_m, Lx, Ly, avg);
      write_twisted_bilayer_meta_txt(dir + "/moire_coarse_v2_meta.txt", nx_m,
                                     ny_m, Lx, Ly, bin_x, bin_y, nn_tol,
                                     moire_coarse_nn_tol_scale);
   }

   if (moire_coarse_write_twisted_double_bilayer) {
      TwistedDoubleBilayerMoireAverages double_avg =
          coarseGrainTwistedDoubleBilayerMoireBonds(
              g_moire_spin_interactions, nx_m, ny_m, Lx, Ly, cell_vol_SI,
              nn_tol, interlayer_tol);
      write_twisted_double_bilayer_binary_v1(
          dir + "/moire_twisted_double_bilayer_coarse.bin", nx_m, ny_m, Lx, Ly,
          double_avg);
      write_twisted_double_bilayer_meta_txt(
          dir + "/moire_twisted_double_bilayer_coarse_meta.txt", nx_m, ny_m,
          Lx, Ly, bin_x, bin_y, nn_tol, moire_coarse_nn_tol_scale);
   }

   if (!moire_coarse_write_twisted_bilayer &&
       !moire_coarse_write_twisted_double_bilayer) {
      std::cerr << "moire_wire: both moire export flags are false; nothing "
                   "written\n";
   }
}
