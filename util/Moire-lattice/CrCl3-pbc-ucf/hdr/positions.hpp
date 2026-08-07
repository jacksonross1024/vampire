#ifndef POSITION_HPP
#define POSITION_HPP

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdint>

#include <unistd.h>
// #include <array>

   extern double twist_angle;

   extern int number_of_unit_cells_x;
   extern int number_of_unit_cells_y;
   extern int microcell_scale_x;
   extern int microcell_scale_y;
   extern int microcell_Nx;
   extern int microcell_Ny;

   extern double a0x;
   extern double a0y;
   extern double a1x;
   extern double a1y;
   extern double c0;
   extern double a0z;

   // Effective microcell pitch for binning. 0 → use scale*a0x / scale*a1y.
   // After PBC id, set so Px*ax_eff = |Tx|, Py*ay_eff = |Ty|.
   extern double microcell_ax_eff;
   extern double microcell_ay_eff;

   // Microcell binning with half-cell grid offset: boundaries at (n+1/2)·cell.
   // ux = floor((x + ax/2)/ax), uy = floor((y + ay/2)/ay).
   inline double microcell_ax() {
      if(microcell_ax_eff > 0.0) return microcell_ax_eff;
      return static_cast<double>(microcell_scale_x) * a0x;
   }
   inline double microcell_ay() {
      if(microcell_ay_eff > 0.0) return microcell_ay_eff;
      return static_cast<double>(microcell_scale_y) * a1y;
   }
   inline int microcell_index_x(double x) {
      const double ax = microcell_ax();
      return static_cast<int>(std::floor((x + 0.5 * ax + 1e-9) / ax));
   }
   inline int microcell_index_y(double y) {
      const double ay = microcell_ay();
      return static_cast<int>(std::floor((y + 0.5 * ay + 1e-9) / ay));
   }

   extern int num_atoms;
   extern int num_nm_atoms;
   extern double twist_loction;
   extern uint64_t total_atoms;
   extern int total_nm_atoms;
   extern int num_above_atoms;
   extern int num_below_atoms;

   extern double J_inter_scaling;
   extern double J_twist_reduction;
   extern double J_intra_reduction;
   extern double J_prist_reduction;
   extern double DMI_inter_scaling;
   extern double DMI_sub_scaling;
   extern double DMI_sub_vector_x;
   extern double DMI_sub_vector_y;
   extern double DMI_sub_vector_z;
   
   class spin {
      public:
         double x;
         double y;
         double z;
         int Gx = 0;
         int Gy = 0;
         int Gz = 0;
         double S;
         uint64_t id = 0;
         uint64_t original_id = 0;
         int l_id;
         int h_id;
         int unit_x;
         int unit_y;
         int unit_x_lr;
         int unit_y_lr;
         int dx = 0;
         int dy = 0;
         
         int inter_twist1_count = 0;
         int inter_twist2_count = 0;
         int inter_twist3_count = 0;

         double J_inter_twist1 = 0;
         double J_inter_twist2 = 0;
         double J_inter_twist3 = 0;
         double Dx_inter_twist1 = 0;
         double Dx_inter_twist2 = 0;
         double Dx_inter_twist3 = 0;
         double Dy_inter_twist1 = 0;
         double Dy_inter_twist2 = 0;
         double Dy_inter_twist3 = 0;
         double Dz_inter_twist1 = 0;
         double Dz_inter_twist2 = 0;
         double Dz_inter_twist3 = 0;


         int inter1_count = 0;
         int inter2_count = 0;
         int inter3_count = 0;

         double J_inter1 = 0;
         double J_inter2 = 0;
         double J_inter3 = 0;
         double Dx_inter1 = 0;
         double Dx_inter2 = 0;
         double Dx_inter3 = 0;
         double Dy_inter1 = 0;
         double Dy_inter2 = 0;
         double Dy_inter3 = 0;
         double Dz_inter1 = 0;
         double Dz_inter2 = 0;
         double Dz_inter3 = 0;


         int intra1_count = 0;
         int intra2_count = 0;
         int intra3_count = 0;
         
         double J_intra1 = 0;
         double J_intra2 = 0;
         double J_intra3 = 0;
         double Dx_intra1 = 0;
         double Dx_intra2 = 0;
         double Dx_intra3 = 0;
         double Dy_intra1 = 0;
         double Dy_intra2 = 0;
         double Dy_intra3 = 0;
         double Dz_intra1 = 0;
         double Dz_intra2 = 0;
         double Dz_intra3 = 0;
   };

   inline void rebin_atom_microcells(spin &at) {
      at.unit_x_lr = microcell_index_x(at.x);
      at.unit_y_lr = microcell_index_y(at.y);
      at.unit_x = at.unit_x_lr;
      at.unit_y = at.unit_y_lr;
   }

   class interaction {
      public:
         int id_i = -1;
         int id_j = -1;
         double J= 0.0;
         double Dx= 0.0;
         double Dy= 0.0;
         double Dz= 0.0;
         int pbc_x = 0;
         int pbc_y = 0;
         int pbc_z = 0;
   };

   extern std::vector < spin > atom;
   extern std::vector < spin > nm_atom;
   extern std::vector < spin > row1;
   extern std::vector < spin > row2;
   extern std::vector < spin > row3;
   extern std::vector < spin > row4;
   extern std::vector < spin > all_m_atoms;
   extern std::vector < spin > all_nm_atoms;
   extern std::vector < spin > all_m_atoms_offset;
   extern std::vector < spin > new_moire_lattice;
   extern std::vector < std::vector < std::vector <int> > > unit_cell_shifts;
   extern std::vector < std::vector < std::vector <double> > > config_energy;
   extern std::vector<std::vector<std::vector<double> > > global_config_energy;

   // Primary spin-moiré microcell period (set by build_spin_moire_cell on PBC success).
   // Used to fold image-tile shift lookups so DMI/J registry is periodic with Tx/Ty.
   extern int moire_period_Px;
   extern int moire_period_Py;
   extern int moire_period_i0;
   extern int moire_period_j0;
   extern bool moire_period_valid;

   // Orthogonal (possibly rotated) Tx/Ty + primary G=0 atoms for deferred Vampire UCF export.
   extern double moire_origin_x, moire_origin_y;
   extern double moire_tx_x, moire_tx_y, moire_ty_x, moire_ty_y;
   // Cached R(-φ) from unit Tx: c = tx_x/|Tx|, s = -tx_y/|Tx| (no atan2 in hot path).
   extern double moire_ucf_c, moire_ucf_s;
   extern std::vector<spin> moire_primary_atoms;
   
   bool inside_system(double sx, double sy, double x, double y, double offset);
   void read_in_atoms(std::string filename, int n_atoms, std::vector <spin > &atom2);
   void read_in_inter_exchanges(std::string J, std::string Dx, std::string Dy, std::string Dz, std::vector<std::vector<double> > &Eij);
   void read_in_intra_exchanges(std::string filename, std::vector<std::vector<double > > &Eij_1NN, \
                                                   std::vector<std::vector<double > > &Eij_2NN, \
                                                   std::vector<std::vector<double > > &Eij_3NN );
   void create_magnetic_atom_list_central(std::string filename);

   // Moiré period result: either from find_moire_periodicity() or user-provided to skip detection.
   struct MoirePeriodResult {
      double x0 = 0, y0 = 0, Lx = 0, Ly = 0;
      double tx_x = 0, tx_y = 0, ty_x = 0, ty_y = 0;
      int Px = 0, Py = 0, i0 = 0, j0 = 0;
      bool from_detection = false;  // true if from find_moire_periodicity; false if user-provided
   };

   // Find moiré periodicity from global_config_energy, unit_cell_shifts, all_m_atoms. Returns true on success.
   // If all_out is non-null, fills it with ranked candidate cells (best first) for PBC retry.
   bool find_moire_periodicity(MoirePeriodResult& out, std::vector<MoirePeriodResult>* all_out = nullptr);

   // Build a spin-moiré primary cell and enforce atomic-level PBC.
   // moire_area: {x0,y0,Lx,Ly[,tx_x,tx_y,ty_x,ty_y[,Px,Py,i0,j0]]}. Lx=Ly=0 → detect.
   // Optional [4..7] = orthogonal (possibly rotated) Tx/Ty for tiling; all-zero ⇒ axis-aligned.
   // Optional [8..11] overrides floor()-derived Px,Py,i0,j0.
   // Crop = half-open AA parallelogram; tile; rezero BL→(0,0); recompute bins.
   // Tile/exchange stay in code frame; call write_vampire_ucf_rotated() after interactions.
   void build_spin_moire_cell(double moire_area[12], const MoirePeriodResult* use_given = nullptr);

   // Write Vampire-orthogonal header.ucf + atom_positions.xyz (positions rotated by -φ).
   void write_vampire_ucf_rotated();
   // Same R(-φ) as atom UCF; uses cached moire_ucf_c/s (set once in build_spin_moire_cell).
   inline void rotate_dmi_for_vampire(double& Dx, double& Dy){
      const double Dx_r = moire_ucf_c*Dx + moire_ucf_s*Dy;
      const double Dy_r = -moire_ucf_s*Dx + moire_ucf_c*Dy;
      Dx = Dx_r; Dy = Dy_r;
   }

#endif