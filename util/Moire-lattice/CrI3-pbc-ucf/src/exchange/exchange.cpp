#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <climits>
#include "initialise.hpp"
#include "exchange.hpp"
#include "moire_bond_record.hpp"
#include "moire_wire.hpp"
#include <list>

// System headers
#include <chrono>

// --- Circmean + signed-y binning fix (plan shift_binning_fix Phase D follow-up):
// signed fractional wrap (no fabs(dy)); circular-mean average of 0..99 sliding indices;
// neighbour-fill empty ortho cells; dual i/j map lookup.
// Backup: src_circmean_backup/exchange.cpp (pre-circmean). Circmean+fabs archive: 1.1_circmean/
static void lab_disp_to_shift_indices(double min_x, double min_y, int &sx, int &sy) {
   // Signed fractional wrap into [0,1): maps are NOT y-fold even (Phase A).
   double fy = min_y / a1y;
   double fx = (min_x - fy * a1x) / a0x;
   fx -= std::floor(fx);
   fy -= std::floor(fy);
   sx = static_cast<int>(std::round(fx * 99.0));
   sy = static_cast<int>(std::round(fy * 99.0));
   if(sx == 100) sx = 0;
   if(sy == 100) sy = 0;
   if(sx < 0) sx += 100;
   if(sy < 0) sy += 100;
}

static void circular_mean_shift(double sum_cx, double sum_sx, double sum_cy, double sum_sy,
                                int occ, int &sx, int &sy) {
   if(occ <= 0) { sx = 66; sy = 0; return; }
   double ang_x = std::atan2(sum_sx, sum_cx);
   double ang_y = std::atan2(sum_sy, sum_cy);
   if(ang_x < 0.0) ang_x += 2.0 * M_PI;
   if(ang_y < 0.0) ang_y += 2.0 * M_PI;
   // Period 100 for labels 0..99 (same range as round(frac*99))
   sx = static_cast<int>(std::round(ang_x * 100.0 / (2.0 * M_PI))) % 100;
   sy = static_cast<int>(std::round(ang_y * 100.0 / (2.0 * M_PI))) % 100;
   if(sx < 0) sx += 100;
   if(sy < 0) sy += 100;
}

// Fold absolute microcell indices into the primary spin-moiré period so image
// tiles (Gy≠0) reuse the same stacking/DMI registry as G=0.
static void fold_microcell_to_primary(int &ux, int &uy) {
   if(!moire_period_valid || moire_period_Px <= 0 || moire_period_Py <= 0) return;
   auto posmod = [](int a, int m)->int {
      int r = a % m;
      return (r < 0) ? r + m : r;
   };
   ux = moire_period_i0 + posmod(ux - moire_period_i0, moire_period_Px);
   uy = moire_period_j0 + posmod(uy - moire_period_j0, moire_period_Py);
}

static void intra_shift_at_atom(const spin &at, int &sx, int &sy) {
   int ux = at.unit_x_lr;
   int uy = at.unit_y_lr;
   fold_microcell_to_primary(ux, uy);
   if(ux < 0 || uy < 0 ||
      ux >= (int)unit_cell_shifts.size() ||
      uy >= (int)unit_cell_shifts[ux].size()) {
      sx = 66; sy = 0;
      return;
   }
   sx = unit_cell_shifts[ux][uy][1];
   sy = unit_cell_shifts[ux][uy][2];
}

static void fill_empty_shifts_from_neighbours(
   std::vector<std::vector<std::vector<int>>> &shifts) {
   const int Nx = static_cast<int>(shifts.size());
   if(Nx == 0) return;
   const int Ny = static_cast<int>(shifts[0].size());
   const int Rmax = 4;
   int nfill = 0;
   for(int i = 0; i < Nx; ++i){
      for(int j = 0; j < Ny; ++j){
         if(shifts[i][j][0] > 0) continue;
         int best_d2 = INT_MAX;
         int bsx = 66, bsy = 0;
         for(int di = -Rmax; di <= Rmax; ++di){
            for(int dj = -Rmax; dj <= Rmax; ++dj){
               if(di == 0 && dj == 0) continue;
               const int ni = i + di, nj = j + dj;
               if(ni < 0 || nj < 0 || ni >= Nx || nj >= Ny) continue;
               if(shifts[ni][nj][0] <= 0) continue;
               const int d2 = di * di + dj * dj;
               if(d2 < best_d2){
                  best_d2 = d2;
                  bsx = shifts[ni][nj][1];
                  bsy = shifts[ni][nj][2];
               }
            }
         }
         if(best_d2 < INT_MAX) nfill++;
         shifts[i][j][1] = bsx;
         shifts[i][j][2] = bsy;
      }
   }
   if(nfill > 0){
      std::cerr << "CIRCMEAN: neighbour-filled " << nfill << " empty unit_cell_shifts cells\n";
   }
}

static void lookup_intra_map(const std::vector<std::vector<double>> &Eij,
                             int sx, int sy, int theta,
                             double &J, double &Dx, double &Dy, double &Dz) {
   const auto &row = Eij[sx * 100 + sy];
   const int o = theta * 4;
   J  = row[o];
   Dx = row[o + 1];
   Dy = row[o + 2];
   Dz = row[o + 3];
}

// Inter map grid encoding (must match read_in_inter_exchanges).
static constexpr double INTER_MAP_A0 = 7.276;
static constexpr double INTER_MAP_A1 = 7.402;
static constexpr double INTER_MAP_DX = 0.02 * INTER_MAP_A0;
static constexpr double INTER_MAP_DY = 0.02 * INTER_MAP_A1;

static inline int inter_map_index(double dx, double dy) {
   int j = (int)std::lround((dx + INTER_MAP_A0) / INTER_MAP_DX);
   int i = (int)std::lround(99.0 - (dy + INTER_MAP_A1) / INTER_MAP_DY);
   if(j < 0) j = 0; else if(j > 99) j = 99;
   if(i < 0) i = 0; else if(i > 99) i = 99;
   return i * 100 + j;
}

static double g_c_bot = 1.0, g_s_bot = 0.0; // S=1,2: R(-φ/2)
static double g_c_top = 1.0, g_s_top = 0.0; // S=3,4: R(+φ/2)

// simple class for performing code timing
class vtimer_t{

private:
   std::chrono::high_resolution_clock::time_point start_time;
   std::chrono::high_resolution_clock::time_point end_time;

public:
   // start the timer
   void start(){
      start_time = std::chrono::high_resolution_clock::now();
   }

   // stop the timer
   void stop(){
      end_time = std::chrono::high_resolution_clock::now();
   }

   // get the elapsed time in milliseconds
   double elapsed_time(){

      // get current time
      std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();

      // work out elapsed time
      return 1.e-9*double(std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count());

   }
};


//set nearest neighbour  distances (in plane nn 1,2,3)
// CrI3 honeycomb shells (a≈6.93 Å): 1NN≈4.00, 2NN≈6.93, 3NN≈8.00.
double intra_nn_dist_1 = 5.465; //A  midpoint 1NN(~4.00)–2NN(~6.93)
double intra_nn_dist_2 = 7.465; //A  midpoint 2NN(~6.93)–3NN(~8.00)
double intra_nn_dist_3 = 8.50;  //A  above 3NN

double inter_nn_dist_1 = 7.0;
double inter_nn_dist_2 = 7.77;
double inter_nn_dist_3 = 9.99;

double inter_AB_dist_1 = 7.0;
double inter_AB_dist_2 = 8.30;
double inter_AB_dist_3 = 9.99;
// double nn_dist_3 = a0x*pow(1.3333333333,0.5);
double nn_dist_1;
double nn_dist_2;
double nn_dist_3;

double max_range = 9.99;
//Set exchange interaction values and associated constants
double eVtoJ = 1.602176634e-19;
double J_constant = 1.0;//eVtoJ/1000.0; //1 meV
double J_intra_1=2.5*J_constant;
double J_intra_2=0.75*J_constant;
double J_intra_3=-0.01*J_constant;
// Easy-axis anisotropy on symmetric exchange: added only to Jzz (J33) for 1NN
double Jz_1NN_aniso = 0.003 * J_constant; // meV
// Biquadratic exchange (--bq): constant, 1NN intralayer only
double J_bq = 0.44 * J_constant; // meV

double Jinter1_AB = 0.1379*J_constant;
double Jinter2_AB = -0.10125*J_constant;
double Jinter3_AB = -0.035*J_constant;

double Jintra1_AB = (2.7657+2.7657+2.4827)*J_constant/3.0;
double Jintra2_AB = (0.7092+0.731+0.7652)*J_constant/3.0;
double Jintra3_AB = (-0.1074-0.1074-0.0994)*J_constant/3.0;

double D_intra_x_constant = 0.0*J_constant;
double D_intra_y_constant =  (-0.0222+0.0767+0.0767)*J_constant/3.0;
double D_intra_z_constant = (-0.0399-0.0399-0.0079)*J_constant/3.0;
double D_intra2_x_constant = 0.1076*J_constant;
double D_intra2_y_constant = -0.0093*J_constant;
double D_intra2_z_constant = -0.0766*J_constant;
double D_intra3_x_constant =	0*J_constant;
double D_intra3_y_constant =  0.0176*J_constant;
double D_intra3_z_constant = 0.0659*J_constant;
//set the initial jumber of interactions to zero for counter

// Interfacial substrate DMI amplitude. Set in calc_interactions from
// DMI_sub_scaling (CLI positional #8); default scale is 0 → no substrate DMI.
// When nonzero, S=1/S=4 intra bonds get an extra D ∝ J×r̂ that is intentionally
// not reverse-bond antisymmetric (interface polarity).
double Dx_substrate = 0.0;
double Dy_substrate = 0.0;
uint64_t number_of_interactions = 0;
uint64_t number_of_bq_interactions = 0;

//initialise arrays to store exchage interactions
std::vector < std::vector < double > > Jint;
std::vector < std::vector < double > > Jinter;
std::vector < std::vector < double > > Einter_Cr1;
std::vector < std::vector < double > > Einter_Cr2;
std::vector < std::vector < double > > Einter_Cr3;
std::vector < std::vector < double > > Einter_Cr4;

std::vector < std::vector < double  > > Eintra_Cr1_1NN;
std::vector < std::vector < double  > > Eintra_Cr2_1NN;
std::vector < std::vector < double  > > Eintra_Cr3_1NN;
std::vector < std::vector < double  > > Eintra_Cr4_1NN;

std::vector < std::vector < double  > > Eintra_Cr1_2NN;
std::vector < std::vector < double  > > Eintra_Cr2_2NN;
std::vector < std::vector < double  > > Eintra_Cr3_2NN;
std::vector < std::vector < double  > > Eintra_Cr4_2NN;

std::vector < std::vector < double  > > Eintra_Cr1_3NN;
std::vector < std::vector < double  > > Eintra_Cr2_3NN;
std::vector < std::vector < double  > > Eintra_Cr3_3NN;
std::vector < std::vector < double  > > Eintra_Cr4_3NN;

std::vector < std::vector < double > > Jintra1;
std::vector < std::vector < double > > Jintra2;

std::vector < std::vector < double > > Dx_inter;
std::vector < std::vector < double > > Dy_inter;
std::vector < std::vector < double > > Dz_inter;
std::vector < std::vector < double > > Dx_intra;
std::vector < std::vector < double > > Dy_intra;
std::vector < std::vector < double > > Dz_intra;
std::vector < std::vector < double > > Dx_intra2;
std::vector < std::vector < double > > Dy_intra2;
std::vector < std::vector < double > > Dz_intra2;

std::vector < std::vector < std::vector< double> > > D_intra;
std::vector < std::vector < std::vector< double> > > D_inter;
   // Config/microcell grid size (from initialise: microcell_Nx + 1, microcell_Ny + 1)

// Classify honeycomb 1NN/3NN bond angle in the layer frame.
// Lattice 1NN/3NN outgoing bonds sit on |θ| ∈ {30°, 90°, 150°} (120°-spaced;
// equivalent to a 0°/120°/240° set after a global frame rotation).
// snap_tol is ONLY for floating-point / PBC rounding — not to reassign wrong neighbours.
static int classify_honeycomb_nn_angle(double angle, const char* tag, const spin &ci, const spin &cj){
   angle = std::remainder(angle, 2.0 * M_PI);
   const double theta_deg = std::abs(angle) * 180.0 / M_PI;
   const double cans[3] = {30.0, 90.0, 150.0};
   const int codes[3] = {1, 0, 2}; // Eij table encoding: 30°→1, 90°→0, 150°→2
   const double snap_tol_deg = 2.0; // 1.1° tiles show ~1.9° seam residuals with rotated Tx/Ty
   double best_err = 1e9;
   int best = -1;
   for(int k = 0; k < 3; ++k){
      const double err = std::fabs(theta_deg - cans[k]);
      if(err < best_err){ best_err = err; best = k; }
   }
   if(best_err <= snap_tol_deg) return codes[best];
   #pragma omp critical
   {
      std::cerr << "Fatal: " << tag << " bond angle off honeycomb lattice (wrong neighbour / tile seam).\n"
                << "  |θ|=" << theta_deg << " deg (nearest canonical err=" << best_err
                << " deg; snap≤" << snap_tol_deg << ")\n"
                << "  angle=" << angle << " central id=" << ci.id << " S=" << ci.S
                << " l_id=" << ci.l_id << " G=(" << ci.Gx << "," << ci.Gy << ")"
                << " pos=(" << ci.x << "," << ci.y << ")"
                << "  j id=" << cj.id << " S=" << cj.S << " l_id=" << cj.l_id
                << " G=(" << cj.Gx << "," << cj.Gy << ") original_id=" << cj.original_id
                << " pos=(" << cj.x << "," << cj.y << ")"
                << " dxy=(" << (cj.x-ci.x) << "," << (cj.y-ci.y) << ")" << std::endl;
      std::cerr.flush();
   }
   // exit() from an OpenMP worker is unreliable; abort the process.
   std::exit(1);
   return -1;
}

// Classify honeycomb 2NN bond angle: θ ∈ {0°,60°,120°,180°,240°,300°} (every 60°).
// Index 0..5 matches read_in_intra_exchanges / Eij_2NN layout. Tight snap only.
static int classify_honeycomb_2nn_angle(double angle, const char* tag, const spin &ci, const spin &cj){
   angle = std::remainder(angle, 2.0 * M_PI);
   double theta = angle;
   if(theta < 0.0) theta += 2.0 * M_PI;
   const double theta_deg = theta * 180.0 / M_PI;
   const double snap_tol_deg = 2.0; // match 1NN tolerance for rotated Tx/Ty seams
   double best_err = 1e9;
   int best = -1;
   for(int k = 0; k < 6; ++k){
      const double can = 60.0 * k;
      double err = std::fabs(theta_deg - can);
      if(err > 180.0) err = 360.0 - err;
      if(err < best_err){ best_err = err; best = k; }
   }
   if(best_err <= snap_tol_deg) return best;
   #pragma omp critical
   {
      std::cerr << "Fatal: " << tag << " 2NN bond angle off honeycomb lattice.\n"
                << "  θ=" << theta_deg << " deg (nearest 60°·k err=" << best_err
                << " deg; snap≤" << snap_tol_deg << ")\n"
                << "  angle=" << angle << " central id=" << ci.id << " S=" << ci.S
                << " pos=(" << ci.x << "," << ci.y << ")"
                << "  j id=" << cj.id << " S=" << cj.S
                << " pos=(" << cj.x << "," << cj.y << ")"
                << " dxy=(" << (cj.x-ci.x) << "," << (cj.y-ci.y) << ")" << std::endl;
      std::cerr.flush();
   }
   std::exit(1);
   return -1;
}

void print_interaction_header(){
   std::ofstream outfile3 ("header_interactions.ucf");
   std::string interaction_type;
   if(DMI) interaction_type = "normalised-tensorial";
   else interaction_type = "normalised-isotropic";
   outfile3 << number_of_interactions <<  "\t" << interaction_type << std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;

   if(BQ) {
      std::ofstream bq_hdr("header_interactions_biquadratic.ucf");
      bq_hdr << number_of_bq_interactions << "\t" << "normalised-isotropic" << std::endl;
   }
}

void compute_unit_cell_shifts_from_atoms(const std::vector<spin>& atoms,
   std::vector<std::vector<std::vector<int>>>& out_shifts) {
   const double range = max_range;
   const double bsize = 1.2*range;

   double min[3] = {1.0e8, 1.0e8, 1.0e8};
   double max[3] = {-1.0e8, -1.0e8, -1.0e8};
   for(size_t i = 0; i < atoms.size(); i++){
      if(atoms[i].S == 1 || atoms[i].S == 4) continue;
      double x_i = atoms[i].x;
      double y_i = atoms[i].y;
      double z_i = atoms[i].z;
      if(x_i < min[0]) min[0] = x_i;
      if(y_i < min[1]) min[1] = y_i;
      if(z_i < min[2]) min[2] = z_i;
      if(x_i > max[0]) max[0] = x_i;
      if(y_i > max[1]) max[1] = y_i;
      if(z_i > max[2]) max[2] = z_i;
   }

   int xb = ceil(system_size_x/bsize)+1;
   int yb = ceil(system_size_y/bsize)+1;
   int zb = ceil(system_size_z/bsize)+1;

   std::vector<std::vector<std::vector<std::vector<spin>>>> shift_boxes;
   shift_boxes.resize(xb);
   for(int i = 0; i < xb; i++){
      shift_boxes[i].resize(yb);
      for(int j = 0; j < yb; j++){
         shift_boxes[i][j].resize(zb);
      }
   }

   for(size_t i = 0; i < atoms.size(); i++){
      if(atoms[i].S == 1 || atoms[i].S == 4) continue;
      double x_i = atoms[i].x - min[0];
      double y_i = atoms[i].y - min[1];
      double z_i = atoms[i].z - min[2];
      const double bxi = x_i / bsize;
      const double byi = y_i / bsize;
      const double bzi = z_i / bsize;
      bool x_ok = bxi >= 0 && bxi < xb;
      bool y_ok = byi >= 0 && byi < yb;
      bool z_ok = bzi >= 0 && bzi < zb;
      if(!(x_ok && y_ok && z_ok)) continue;
      shift_boxes[bxi][byi][bzi].push_back(atoms[i]);
   }

   std::vector<std::list<spin>> moire_shift(atoms.size());
   int moire_shift_count = 0;
   for(int i = 0; i < xb; i++){
      for(int j = 0; j < yb; j++){
         for(int k = 0; k < zb; k++){
            for(int dx = -1; dx < 2; dx++){
               for(int dy = -1; dy < 2; dy++){
                  for(int dz = -1; dz < 2; dz++){
                     const int nx = i+dx, ny = j+dy, nz = k+dz;
                     if(nx >= 0 && nx < xb && ny >= 0 && ny < yb && nz >= 0 && nz < zb){
                        for(size_t ai = 0; ai < shift_boxes[i][j][k].size(); ai++){
                           spin atom_i = shift_boxes[i][j][k][ai];
                           if(atom_i.S == 2) continue;
                           const double x_i = atom_i.x, y_i = atom_i.y, z_i = atom_i.z;
                           for(size_t aj = 0; aj < shift_boxes[nx][ny][nz].size(); aj++){
                              spin atom_j = shift_boxes[nx][ny][nz][aj];
                              if(atom_i.S == atom_j.S || atom_i.id == atom_j.id) continue;
                              double adx = atom_j.x - x_i, ady = atom_j.y - y_i;
                              if(adx*adx <= (a0x*a0x) && ady*ady <= (a1y*a1y) && ((atom_i.l_id-2) == atom_j.l_id)){
                                 moire_shift[atom_i.id].push_back(atom_j);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   const double cos_half = cos(0.5 * twist_angle);
   const double sin_half = sin(0.5 * twist_angle);

   out_shifts.resize(microcell_Nx + 1);
   for(int i = 0; i <= microcell_Nx; i++){
      out_shifts[i].resize(microcell_Ny + 1);
      for(int j = 0; j <= microcell_Ny; j++){
         out_shifts[i][j].resize(3, 0);
      }
   }
   std::vector<std::vector<double>> scx(microcell_Nx + 1), ssx(microcell_Nx + 1),
                                    scy(microcell_Nx + 1), ssy(microcell_Nx + 1);
   for(int i = 0; i <= microcell_Nx; i++){
      scx[i].assign(microcell_Ny + 1, 0.0);
      ssx[i].assign(microcell_Ny + 1, 0.0);
      scy[i].assign(microcell_Ny + 1, 0.0);
      ssy[i].assign(microcell_Ny + 1, 0.0);
   }

   for(size_t i = 0; i < moire_shift.size(); i++){
      double min_x = a0x, min_y = a1y, min_r = 100.0;
      spin ref_spin = atoms[i];
      if(ref_spin.S != 3) continue;
      const double ref_lab_x = ref_spin.x * cos_half - ref_spin.y * sin_half;
      const double ref_lab_y = ref_spin.x * sin_half + ref_spin.y * cos_half;

      for(auto shift = moire_shift[i].begin(); shift != moire_shift[i].end(); ++shift){
         spin shift_spin = *shift;
         const double shift_lab_x = shift_spin.x * cos_half - shift_spin.y * sin_half;
         const double shift_lab_y = shift_spin.x * sin_half + shift_spin.y * cos_half;
         double dx = ref_lab_x - shift_lab_x;
         double dy = ref_lab_y - shift_lab_y;
         double dr = (dx*dx + dy*dy);
         if(dr < min_r){
            min_r = dr;
            min_x = dx;
            min_y = dy;
         }
      }

      int dx_cell = ref_spin.unit_x_lr;
      int dy_cell = ref_spin.unit_y_lr;
      int dx = 0, dy = 0;
      lab_disp_to_shift_indices(min_x, min_y, dx, dy);
      if(dx > 99 || dx < 0 || dy > 99 || dy < 0) continue;
      if(dx_cell < 0 || dy_cell < 0 ||
         dx_cell >= (int)out_shifts.size() ||
         dy_cell >= (int)out_shifts[dx_cell].size()) continue;

      out_shifts[dx_cell][dy_cell][0] += 1;
      const double ax = 2.0 * M_PI * dx / 100.0;
      const double ay = 2.0 * M_PI * dy / 100.0;
      scx[dx_cell][dy_cell] += std::cos(ax);
      ssx[dx_cell][dy_cell] += std::sin(ax);
      scy[dx_cell][dy_cell] += std::cos(ay);
      ssy[dx_cell][dy_cell] += std::sin(ay);
   }

   for(int i = 0; i < (int)out_shifts.size(); i++){
      for(int j = 0; j < (int)out_shifts[i].size(); j++){
         int occupancy = out_shifts[i][j][0];
         if(occupancy > 0){
            int sx = 0, sy = 0;
            circular_mean_shift(scx[i][j], ssx[i][j], scy[i][j], ssy[i][j], occupancy, sx, sy);
            out_shifts[i][j][1] = sx;
            out_shifts[i][j][2] = sy;
         } else {
            out_shifts[i][j][1] = 66;
            out_shifts[i][j][2] = 0;
         }
      }
   }
   fill_empty_shifts_from_neighbours(out_shifts);

   char directory[256];
   if(getcwd(directory, sizeof(directory)) == NULL){
      std::cerr << "compute_unit_cell_shifts_from_atoms: getcwd error" << std::endl;
      return;
   }
   std::ofstream shift_file(std::string(directory) + "/new_shifted_constants.txt");
   for(int i = 0; i < (int)out_shifts.size(); i++){
      for(int j = 0; j < (int)out_shifts[i].size(); j++){
         int occupancy = out_shifts[i][j][0];
         int i_shift = out_shifts[i][j][1];
         int j_shift = out_shifts[i][j][2];
         if(i_shift > 99 || j_shift > 99) std::cout << "problems " << i << ", " << j << ", " << occupancy << ", " << i_shift << ", " << j_shift << std::endl;
         shift_file << i << ", " << j << ", " << occupancy << ", " << i_shift << ", " << j_shift << "\n";
      }
   }
   shift_file.close();
}

 
using EijMap = std::vector<std::vector<double>>;
static inline EijMap &intra_map_1nn(int lid) {
   switch(lid) {
      case 2: return Eintra_Cr2_1NN;
      case 3: return Eintra_Cr3_1NN;
      case 4: return Eintra_Cr4_1NN;
      default: return Eintra_Cr1_1NN;
   }
}
static inline EijMap &intra_map_2nn(int lid) {
   switch(lid) {
      case 2: return Eintra_Cr2_2NN;
      case 3: return Eintra_Cr3_2NN;
      case 4: return Eintra_Cr4_2NN;
      default: return Eintra_Cr1_2NN;
   }
}
static inline EijMap &intra_map_3nn(int lid) {
   switch(lid) {
      case 2: return Eintra_Cr2_3NN;
      case 3: return Eintra_Cr3_3NN;
      case 4: return Eintra_Cr4_3NN;
      default: return Eintra_Cr1_3NN;
   }
}
static inline EijMap &inter_map(int lid) {
   switch(lid) {
      case 2: return Einter_Cr2;
      case 3: return Einter_Cr3;
      case 4: return Einter_Cr4;
      default: return Einter_Cr1;
   }
}



void calc_interactions() {

   moire_spin_interactions_clear();

   std::stringstream ss;
             
   // calculate max range squred
   const double range = max_range;//std::max(inter_nn_dist_3, intra_nn_dist_3);
   const double r2 = range*range;
   intra_nn_dist_1 *= intra_nn_dist_1;
   intra_nn_dist_2 *= intra_nn_dist_2;
   intra_nn_dist_3 *= intra_nn_dist_3;

   inter_nn_dist_1 *= inter_nn_dist_1;
   inter_nn_dist_2 *= inter_nn_dist_2;
   inter_nn_dist_3 *= inter_nn_dist_3;

   inter_AB_dist_1 *= inter_AB_dist_1;
   inter_AB_dist_2 *= inter_AB_dist_2;
   inter_AB_dist_3 *= inter_AB_dist_3;

   const double bsize = 1.2*range;

   Jinter1_AB = 0.1379*J_constant;
   Jinter2_AB = -0.10125*J_constant;
   Jinter3_AB = -0.035*J_constant;
   
   Dx_substrate = 0.01819374 * DMI_sub_scaling; // 0.00606458;//  0.009; // 0.01819374;
   std::cout << "Dx_substrate=" << Dx_substrate
             << " (DMI_sub_scaling=" << DMI_sub_scaling
             << "; 0 ⇒ no interfacial substrate DMI)" << std::endl;
   
   std::cout << "Generating Moire unit cell...." << std::flush;
   // calculate min and max xyz
   double min[3] = {1.0e8, 1.0e8, 1.0e8};
   double max[3] = {-1.0e8, -1.0e8, -1.0e8};
   for(int i=0; i < all_m_atoms.size(); i++){
      if(all_m_atoms[i].S == 1 || all_m_atoms[i].S == 4 ) continue;
      double x_i = all_m_atoms[i].x;
      double y_i = all_m_atoms[i].y;
      double z_i = all_m_atoms[i].z;
      if(x_i < min[0]) min[0] = x_i;
      if(y_i < min[1]) min[1] = y_i;
      if(z_i < min[2]) min[2] = z_i;
      if(x_i > max[0]) max[0] = x_i;
      if(y_i > max[1]) max[1] = y_i;
      if(z_i > max[2]) max[2] = z_i;
   }

   // determine number of blocks in x,y,z
    int xb = ceil(system_size_x/bsize)+1;
    int yb = ceil(system_size_y/bsize)+1;
    int zb = ceil(system_size_z/bsize)+1;
   std::cout << "decomposed into <" << xb << ", " << yb << ", " << zb << "> boxes...." << std::flush;
   // create 4D array to generate blocks
   std::vector< std::vector < std::vector < std::vector < spin > > > > shift_boxes;
   shift_boxes.resize(xb);
   for(int i=0; i<xb; i++){
      shift_boxes[i].resize(yb);
      for(int j=0; j<yb; j++){
         shift_boxes[i][j].resize(zb);
      }
   }

   // determine boxid of each atom and save atoms in boxes
   for(int i=0; i < all_m_atoms.size(); i++){
      if(all_m_atoms[i].S == 1 || all_m_atoms[i].S == 4) continue;
      double x_i = all_m_atoms[i].x - min[0];
      double y_i = all_m_atoms[i].y - min[1];
      double z_i = all_m_atoms[i].z - min[2];
      const double bxi = x_i / bsize;
      const double byi = y_i / bsize;
      const double bzi = z_i / bsize;

      // check that boxid is in range
      bool x_ok = bxi >= 0 && bxi < xb;
      bool y_ok = byi >= 0 && byi < yb;
      bool z_ok = bzi >= 0 && bzi < zb;
      if( !(x_ok && y_ok && z_ok) ){
         std::cerr << "Error! Atom " << i << " out of box range " << bxi << "\t" << byi << "\t" << bzi << "\t" << xb << "\t" << yb << "\t" << zb << std::endl;
         std::cout << x_ok << "\t" << y_ok << "\t" << z_ok << "\t" << (x_ok && y_ok && z_ok) << "\t" << !(x_ok && y_ok && z_ok ) << std::endl;
         exit(1);
      }
      // std::cout << "here?" << std::endl;
      // add atom to box list
      shift_boxes[bxi][byi][bzi].push_back(all_m_atoms[i]);
   }

   std::cout << "[complete]" << std::endl;
   
   
   // now calculate neighbour list looping over boxes
   vtimer_t timer;
   timer.start();
   std::vector< std::list< spin > > moire_shift(all_m_atoms.size());
   for(int i=0; i<xb; i++){
      for(int j=0; j< yb; j++){
         for(int k=0; k<zb; k++){

            // loop over offsets
            for(int dx = -1; dx < 2; dx++){
               for(int dy = -1; dy < 2; dy++){
                  for(int dz = -1; dz < 2; dz++){
                     const int nx = i+dx; // neighbour box ids
                     const int ny = j+dy;
                     const int nz = k+dz;
                     
                     const bool x_ok = nx >= 0 && nx < xb;
                     const bool y_ok = ny >= 0 && ny < yb;
                     const bool z_ok = nz >= 0 && nz < zb;
                     
                     if(x_ok && y_ok && z_ok){
                     // only calculate neighbours for all x,y,z indices ok
                        // loop over all atoms in main box
                        for(int ai = 0; ai < shift_boxes[i][j][k].size(); ai++){
                          
                           // get atom number i
                           spin atom_i = shift_boxes[i][j][k][ai];
                           if(atom_i.S == 2) continue;

                           //moire_shift[atom_i.id].push_back(atom_i);
                           const double x_i = atom_i.x;
                           const double y_i = atom_i.y;
                           const double z_i = atom_i.z;

                           for(int aj = 0; aj < shift_boxes[nx][ny][nz].size(); aj++){

                              // get atom number j
                              spin atom_j = shift_boxes[nx][ny][nz][aj];
                              if(atom_i.S == atom_j.S || atom_i.id == atom_j.id) continue;
  
                              // calculate distance
                              const double x_j = atom_j.x;
                              const double y_j = atom_j.y;
                              const double z_j = atom_j.z;
                              double adx = x_j - x_i;
                              double ady = y_j - y_i;

                              if(adx*adx <= (a0x*a0x) && ady*ady <= (a1y*a1y) && ((atom_i.l_id-2) == atom_j.l_id)){
                                 moire_shift[atom_i.id].push_back(atom_j);
                                 
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   // Rotate ref/shift positions by 0.5*twist_angle into orthogonal lab frame so change_x/change_y (unit_cell_shifts) are in lab; twist is applied as bisect so this restores dx,dy to lab.
   const double cos_half = cos(0.5 * twist_angle);
   const double sin_half = sin(0.5 * twist_angle);

   std::vector<std::vector<double>> scx(microcell_Nx + 1), ssx(microcell_Nx + 1),
                                    scy(microcell_Nx + 1), ssy(microcell_Nx + 1);
   for(int i = 0; i <= microcell_Nx; i++){
      scx[i].assign(microcell_Ny + 1, 0.0);
      ssx[i].assign(microcell_Ny + 1, 0.0);
      scy[i].assign(microcell_Ny + 1, 0.0);
      ssy[i].assign(microcell_Ny + 1, 0.0);
   }

   for(int i = 0; i < moire_shift.size(); i++) {
      double min_x = a0x;
      double min_y = a1y;
      double min_r = 100.0;
      spin ref_spin = all_m_atoms[i];
      if(ref_spin.S != 3) continue;
      const double ref_lab_x = ref_spin.x * cos_half - ref_spin.y * sin_half;
      const double ref_lab_y = ref_spin.x * sin_half + ref_spin.y * cos_half;

      std::list<spin> ref_list = moire_shift[i];
      int count = 0;

      for(auto shift = ref_list.begin(); shift != ref_list.end(); ++shift) {
         spin shift_spin = *shift;
         const double shift_lab_x = shift_spin.x * cos_half - shift_spin.y * sin_half;
         const double shift_lab_y = shift_spin.x * sin_half + shift_spin.y * cos_half;

         double dx = ref_lab_x - shift_lab_x;
         double dy = ref_lab_y - shift_lab_y;
         double dr = (dx*dx + dy*dy);
         if(dr < min_r) {
            min_r = dr;
            min_x = dx;
            min_y = dy;
           
         }
         count++;
      }
      
      int dx_cell = ref_spin.unit_x_lr;
      int dy_cell = ref_spin.unit_y_lr;

      int dx = 0, dy = 0;
      lab_disp_to_shift_indices(min_x, min_y, dx, dy);
       if(dx > 99 || dx < 0 || dy > 99 || dy < 0) {
         std::cerr << "shift problem: (" << dx << ", " << dy << ", " << min_x << ", " << min_y << std::endl;
            exit(1);
      }
      if(dx_cell < 0 || dy_cell < 0 ||
         dx_cell >= (int)unit_cell_shifts.size() ||
         dy_cell >= (int)unit_cell_shifts[dx_cell].size()) continue;
      
      unit_cell_shifts[dx_cell][dy_cell][0] += 1;
      const double ax = 2.0 * M_PI * dx / 100.0;
      const double ay = 2.0 * M_PI * dy / 100.0;
      scx[dx_cell][dy_cell] += std::cos(ax);
      ssx[dx_cell][dy_cell] += std::sin(ax);
      scy[dx_cell][dy_cell] += std::cos(ay);
      ssy[dx_cell][dy_cell] += std::sin(ay);
   }

      char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
         std::cerr << "Fatal getcwd error in datalog." << std::endl;
      }
      for(int i = 0; i < (int)unit_cell_shifts.size(); i++){
      for (int j = 0; j < (int)unit_cell_shifts[i].size(); j++) {
         int occupancy = unit_cell_shifts[i][j][0];
         if(occupancy > 0) {
            int sx = 0, sy = 0;
            circular_mean_shift(scx[i][j], ssx[i][j], scy[i][j], ssy[i][j], occupancy, sx, sy);
            unit_cell_shifts[i][j][1] = sx;
            unit_cell_shifts[i][j][2] = sy;
         } else {
            unit_cell_shifts[i][j][1] = 66;
            unit_cell_shifts[i][j][2] = 0;
         }
      }
   }
   fill_empty_shifts_from_neighbours(unit_cell_shifts);
   {
      std::ofstream shift_file(std::string(directory) + "/shifted_constants.txt");
      for(int i = 0; i < (int)unit_cell_shifts.size(); i++){
         for(int j = 0; j < (int)unit_cell_shifts[i].size(); j++){
            int occupancy = unit_cell_shifts[i][j][0];
            int i_shift = unit_cell_shifts[i][j][1];
            int j_shift = unit_cell_shifts[i][j][2];
            if(i_shift > 99 || j_shift > 99) std::cout << "problems " << i << ", " << j << ", " << occupancy << ", " << i_shift << ", " << j_shift << std::endl;
            shift_file << i << ", " << j << ", " << occupancy << ", " << i_shift << ", " << j_shift << "\n";
         }
      }
   }

   // Period determination and crop to primary spin-moiré cell (uses occupancy + unit_cell_shifts).
   // Exchange constants are then computed on the resulting PBC-applied atom list.
   // moire_area = {x0, y0, Lx, Ly, tx_x, tx_y, ty_x, ty_y[, Px, Py, i0, j0]}.
   // Lx=Ly=0 → find_moire_periodicity(). No baked θ-table — periods are re-detected
   // under circmean+signed-y + orthogonal (rotated OK) Tx/Ty (unvalidated until re-bake).
   // Tile/exchange stay in code frame; Vampire UCF via write_vampire_ucf_rotated() after exchange.
   const double theta_deg = twist_angle * 180.0 / M_PI;
   double moire_area[12] = {0.0}; // auto-detect by default
   auto near = [&](double t){ return std::fabs(theta_deg - t) < 0.05; };
   // Analytic untwisted cell only (not a spin-moiré bake).
   if(!bake_only && near(0.0)) {
      // Untwisted: minimal orthogonal cell a × √3 a (= a0x × 2*a1y), 16 Cr sites.
      const double Lx0 = a0x;
      const double Ly0 = 2.0 * a1y;
      moire_area[0] = 0.5 * system_size_x;
      moire_area[1] = 0.5 * system_size_y;
      moire_area[2] = Lx0;
      moire_area[3] = Ly0;
      moire_area[8] = 1;   // Px
      moire_area[9] = 2;   // Py
      moire_area[10] = microcell_index_x(moire_area[0]);
      moire_area[11] = microcell_index_y(moire_area[1]);
   }
   std::cout << "moire_area for " << theta_deg << " deg: {"
             << moire_area[0] << ", " << moire_area[1] << ", "
             << moire_area[2] << ", " << moire_area[3] << "} TxTy tx=("
             << moire_area[4] << "," << moire_area[5] << ") ty=("
             << moire_area[6] << "," << moire_area[7] << ")" << std::endl;
   build_spin_moire_cell(moire_area);
   if(bake_only) {
      std::cout << "bake-only: period detection complete; skipping exchange." << std::endl;
      return;
   }

   // Retile changes atom coordinates; rebuild neighbour boxes from the new bounds.
   min[0] = min[1] = min[2] = 1.0e8;
   max[0] = max[1] = max[2] = -1.0e8;
   for(size_t i = 0; i < all_m_atoms.size(); ++i){
      const double x_i = all_m_atoms[i].x, y_i = all_m_atoms[i].y, z_i = all_m_atoms[i].z;
      if(x_i < min[0]) min[0] = x_i;
      if(y_i < min[1]) min[1] = y_i;
      if(z_i < min[2]) min[2] = z_i;
      if(x_i > max[0]) max[0] = x_i;
      if(y_i > max[1]) max[1] = y_i;
      if(z_i > max[2]) max[2] = z_i;
   }
   xb = static_cast<int>(std::ceil((max[0] - min[0]) / bsize)) + 2;
   yb = static_cast<int>(std::ceil((max[1] - min[1]) / bsize)) + 2;
   zb = static_cast<int>(std::ceil((max[2] - min[2]) / bsize)) + 2;
   if(xb < 1) xb = 1;
   if(yb < 1) yb = 1;
   if(zb < 1) zb = 1;

      std::vector< std::vector < std::vector < std::vector < spin > > > > boxes;
   boxes.resize(xb);
   for(int i=0; i<xb; i++){
      boxes[i].resize(yb);
      for(int j=0; j<yb; j++){
         boxes[i][j].resize(zb);
      }
   }
   // determine boxid of each atom and save atoms in boxes
   for(int i=0; i < all_m_atoms.size(); i++){
      
      double x_i = all_m_atoms[i].x - min[0];
      double y_i = all_m_atoms[i].y - min[1];
      double z_i = all_m_atoms[i].z - min[2];
      const double bxi = x_i / bsize;
      const double byi = y_i / bsize;
      const double bzi = z_i / bsize;

      // check that boxid is in range
      bool x_ok = bxi >= 0 && bxi < xb;
      bool y_ok = byi >= 0 && byi < yb;
      bool z_ok = bzi >= 0 && bzi < zb;
      if( !(x_ok && y_ok && z_ok) ){
         std::cerr << "Error! Atom " << i << " out of box range " << bxi << "\t" << byi << "\t" << bzi << "\t" << xb << "\t" << yb << "\t" << zb << std::endl;
         std::cout << x_ok << "\t" << y_ok << "\t" << z_ok << "\t" << (x_ok && y_ok && z_ok) << "\t" << !(x_ok && y_ok && z_ok ) << std::endl;
         exit(1);
      }
      // std::cout << "here?" << std::endl;
      // add atom to box list
      boxes[bxi][byi][bzi].push_back(all_m_atoms[i]);
   }

   std::cout << "[complete]" << std::endl;
   // exit(1);

   // Rotation to lab frame for adx,ady used in match_inter_exchange (Eij table is in lab frame, like unit_cell_shifts)
   const double cos_half_ex = cos(0.5 * twist_angle);
   const double sin_half_ex = sin(0.5 * twist_angle);
   g_c_bot = cos(-0.5 * twist_angle);
   g_s_bot = sin(-0.5 * twist_angle);
   g_c_top = cos_half_ex;
   g_s_top = sin_half_ex;

      std::ofstream lattice_file(std::string(directory) + "/moire_lattice_vectors.txt");
   // Per-thread text buffers; flush under critical to avoid torn lines and to cap RAM.
   const std::streamsize flush_bytes = 32 * 1024 * 1024; // 32 MiB
   int nthr = omp_get_max_threads();
   if(write_config_output && nthr > 4) nthr = 4;
   #pragma omp parallel num_threads(nthr) reduction(+:number_of_interactions,number_of_bq_interactions)
   {
      #pragma omp single
      std::cout << "preparing Moire exchange with " << omp_get_num_threads() << " omp threads" << std::endl;

      // Full microcell config grids are ~0.5GB/thread — only when --config-atoms.
      std::vector<std::vector<std::vector<double> > > local_config_energy;
      if(write_config_output){
         local_config_energy.resize(microcell_Nx + 1);
         for(int i = 0; i <= microcell_Nx; i++) {
            local_config_energy[i].resize(microcell_Ny + 1);
            for(int j = 0; j <= microcell_Ny; j++) {
               local_config_energy[i][j].resize(4*16,0.0);
            }
         }
      }
   std::stringstream otext;
   otext.precision(6);
   std::stringstream bq_text;
   bq_text.precision(6);

   std::vector<spin_interaction> local_moire_bonds;
      local_moire_bonds.reserve(all_m_atoms.size()*20/8);
   std::stringstream lattice_text;
   auto flush_interaction_buffers = [&](){
      // Snapshot then clear thread-local buffers; write under critical (no torn lines).
      const std::string chunk = otext.str();
      const std::string bq_chunk = BQ ? bq_text.str() : std::string();
      const std::string lat_chunk = lattice_text.str();
      otext.str(std::string());
      otext.clear();
      otext.precision(6);
      bq_text.str(std::string());
      bq_text.clear();
      bq_text.precision(6);
      lattice_text.str(std::string());
      lattice_text.clear();
      if(chunk.empty() && bq_chunk.empty() && lat_chunk.empty()) return;
      #pragma omp critical
      {
         if(!chunk.empty()) outfile4 << chunk;
         if(!bq_chunk.empty()) outfile_bq << bq_chunk;
         if(!lat_chunk.empty()) lattice_file << lat_chunk;
      }
   };

   #pragma omp for schedule(dynamic) collapse(2)
   for(int i=0; i<xb; i++){
      for(int j=0; j< yb; j++){
         if(j == 0 && (i % std::max(1, xb/10) == 0)){
            #pragma omp critical
            std::cout << "." << std::flush;
         }
         for(int k=0; k<zb; k++){

            // loop over offsets
            for(int dx = -1; dx < 2; dx++){
               for(int dy = -1; dy < 2; dy++){
                  for(int dz = -1; dz < 2; dz++){
                     const int nx = i+dx; // neighbour box ids
                     const int ny = j+dy;
                     const int nz = k+dz;
                     
                     const bool x_ok = nx >= 0 && nx < xb;
                     const bool y_ok = ny >= 0 && ny < yb;
                     const bool z_ok = nz >= 0 && nz < zb;

                     if(x_ok && y_ok && z_ok){
                     // only calculate neighbours for all x,y,z indices ok
                        // loop over all atoms in main box
                        const auto &box_i = boxes[i][j][k];
                        const auto &box_j = boxes[nx][ny][nz];
                        for(int ai = 0; ai < (int)box_i.size(); ai++){
                           const spin &atom_i = box_i[ai];
                           if(atom_i.S == 5) continue;
                           if(atom_i.Gx != 0 || atom_i.Gy != 0) continue;
                           const double x_i = atom_i.x;
                           const double y_i = atom_i.y;
                           const double z_i = atom_i.z;

                           for(int aj = 0; aj < (int)box_j.size(); aj++){
                              const spin &atom_j = box_j[aj];
                              if(atom_i.id == atom_j.id) continue;
                              const double adx = atom_j.x - x_i;
                              const double ady = atom_j.y - y_i;
                              const double adz = atom_j.z - z_i;
                              const double dL2 = adx*adx + ady*ady + adz*adz;
                              if(dL2 >= r2) continue;

                              std::array<float, 4> exchange({0.0f,0.0f,0.0f,0.0f});
                              spin &atom_i_mut = const_cast<spin&>(atom_i);
                              spin &atom_j_mut = const_cast<spin&>(atom_j);

                              if(atom_i.S == atom_j.S) {
                                 if(atom_i.l_id < 1 || atom_i.l_id > 4 || atom_j.l_id < 1 || atom_j.l_id > 4) continue;
                                 double angle_i = std::atan2(ady, adx);
                                 double angle_j = angle_i + M_PI;
                                 if(atom_i.l_id == 1 || atom_i.l_id == 2) {
                                    angle_i += 0.5*twist_angle;
                                    angle_j += 0.5*twist_angle;
                                 } else {
                                    angle_i -= 0.5*twist_angle;
                                    angle_j -= 0.5*twist_angle;
                                 }
                                 if(dL2 < intra_nn_dist_1) {
                                    exchange = match_intra1_exchange(angle_i, angle_j, atom_i_mut, atom_j_mut,
                                       intra_map_1nn(atom_i.l_id), intra_map_1nn(atom_j.l_id));
                                 } else if (dL2 < intra_nn_dist_2) {
                                    exchange = match_intra2_exchange(angle_i, angle_j, atom_i_mut, atom_j_mut,
                                       intra_map_2nn(atom_i.l_id), intra_map_2nn(atom_j.l_id));
                                 } else if (dL2 < intra_nn_dist_3) {
                                    exchange = match_intra3_exchange(angle_i, angle_j, atom_i_mut, atom_j_mut,
                                       intra_map_3nn(atom_i.l_id), intra_map_3nn(atom_j.l_id));
                                 } else continue;
                              } else if(atom_i.h_id == atom_j.h_id) {
                                 if(atom_i.l_id < 1 || atom_i.l_id > 4 || atom_j.l_id < 1 || atom_j.l_id > 4) continue;
                                 const double adx_lab = adx * cos_half_ex - ady * sin_half_ex;
                                 const double ady_lab = adx * sin_half_ex + ady * cos_half_ex;
                                 exchange = match_inter_exchange(atom_i.id, atom_j.id, adx_lab, ady_lab, dL2,
                                    inter_map(atom_i.l_id), inter_map(atom_j.l_id));
                                 exchange[0] *= (float)J_twist_reduction;
                              } else {
                                 if(dL2 <= inter_AB_dist_1)
                                    exchange[0] = (float)(Jinter1_AB * J_prist_reduction);
                                 else if(dL2 <= inter_AB_dist_2)
                                    exchange[0] = (float)(Jinter2_AB * J_prist_reduction);
                                 else
                                    exchange[0] = (float)(Jinter3_AB * J_prist_reduction);
                              }

                                 if(write_config_output)
                                    spin_config_energy(atom_i_mut, dL2, atom_j_mut, exchange, local_config_energy);

                                   if(DMI) {
                                      // Jzz = isotropic J + optional 1NN easy-axis anisotropy (not 2NN/3NN)
                                      const double Jzz = exchange[0] + ((dL2 < intra_nn_dist_1) ? Jz_1NN_aniso : 0.0);
                                      // Rotate in-plane DMI into Vampire axes (same -φ as atom UCF).
                                      double Dx_v = exchange[1], Dy_v = exchange[2];
                                      rotate_dmi_for_vampire(Dx_v, Dy_v);
                                      otext << number_of_interactions <<  "\t" << atom_i.original_id << '\t' << atom_j.original_id << '\t' << atom_j.Gx << '\t' << atom_j.Gy << '\t' << 0 << '\t' <<\
                                                //xx                     xy-> Dz                 xz -> -Dy
                                                  exchange[0] << "\t" << exchange[3] << "\t" << -Dy_v << "\t" << \
                                                //yx -> -Dz              yy                      yz -> Dx
                                                 -exchange[3] << "\t" << exchange[0] << "\t" <<  Dx_v << "\t" << \
                                                //zx -> Dy               zy -> -Dx               zz (J33)
                                                  Dy_v << "\t" <<-Dx_v << "\t" <<  Jzz << "\n"; }

                                    // if(DMI) {  otext << number_of_interactions <<  "\t" << atom_i.id << '\t' << atom_j.id << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' <<\
                                    // //xx                     xy-> Dz                 xz -> -Dy
                                    //    exchange[0] << "\t" << exchange[1] << "\t" << exchange[2] << "\t" << \
                                    // //yx -> -Dz              yy                      yz -> Dx
                                    //    exchange[3] <<  "\n"; }
                                    else {otext << number_of_interactions <<  "\t" << atom_i.original_id << '\t' << atom_j.original_id << '\t' << atom_j.Gx << '\t' << atom_j.Gy << '\t' << 0 << '\t' <<\
                                       //xx                     xy-> Dz                 xz -> -Dy
                                         exchange[0] << "\t" << 0.0 << "\t" << 0.0 << "\t" << \
                                       //yx -> -Dz              yy                      yz -> Dx
                                       0.0 << "\t" << exchange[0] << "\t" <<  0.0 << "\t" << \
                                       //zx -> Dy               yz -> -Dx               zz
                                       0.0 << "\t" << 0.0 << "\t" <<  exchange[0]*1.0 << "\n"; }
                                   
                                    number_of_interactions++;
                                    local_moire_bonds.push_back(
                                        make_spin_interaction_from_exchange(
                                            atom_i, atom_j, exchange, DMI));
                                    if(otext.tellp() >= flush_bytes) flush_interaction_buffers();

                                    // Biquadratic (--bq): constant J_bq, 1NN intralayer only
                                    if(BQ && atom_i.S == atom_j.S && dL2 < intra_nn_dist_1) {
                                       bq_text << number_of_bq_interactions << "\t" << atom_i.original_id << '\t' << atom_j.original_id
                                               << '\t' << atom_j.Gx << '\t' << atom_j.Gy << '\t' << 0 << '\t' << J_bq << "\n";
                                       number_of_bq_interactions++;
                                    }
                           } // end of j atom loop

                        } // end of i atom loop

                     } // end of protection statement
                  }
               }
            }// end of offset loops
         }
      }
   } // end collapse(2) i,j
      flush_interaction_buffers();
      #pragma omp critical
      {
         moire_spin_interactions_merge(local_moire_bonds);
         std::cout << omp_get_thread_num() << " calculated " << number_of_interactions << " of interactions" << std::endl;

         if(write_config_output){
            for(int ii = 0; ii < (int)global_config_energy.size(); ii++) {
               for(int jj = 0; jj < (int)global_config_energy[ii].size(); jj++) {
                  for(int kk =0; kk < (int)global_config_energy[ii][jj].size(); kk++){
                     global_config_energy[ii][jj].at(kk) += local_config_energy[ii][jj].at(kk);
                  }
               }
            }
         }
      }
      
   }
   lattice_file << std::flush;
   outfile4 << std::flush;
   outfile4.close();
   if(BQ) {
      outfile_bq << std::flush;
      outfile_bq.close();
   }

   timer.stop();
   // std::cout << "done!  << std::endl;
   std::cout << number_of_interactions << " [completed] [" << timer.elapsed_time() << " s]" << std::endl;
   if(BQ) std::cout << number_of_bq_interactions << " biquadratic interactions" << std::endl;

   // Vampire UCF: rotate G=0 positions into pure (Lx,0)/(0,Ly); DMI already rotated at write.
   write_vampire_ucf_rotated();

   moire_spin_interactions_finalize_and_write(directory, moire_area);

   if(!write_config_output){
      std::cout << "skipping config_energy dumps (pass --config-atoms to enable)" << std::endl;
      return;
   }

   std::cout << "outputting extra file info:" << std::endl;
   

   std::cout << "config atoms started..." << std::flush;
      std::ofstream config_output(std::string(directory) + "/config_energy_atomic.txt");

      if(!config_output.is_open()) {std::cout << "config energy did not open" << std::endl; exit(1);}
      for(int i = 0; i < all_m_atoms.size(); i++) {
         spin atom = all_m_atoms[i];
         // for(int j = 0; j < number_of_unit_cells_y; j++){
            // double bottom_occ = config_energy[i][j][(2-1)*5+0];
            // double top_occ = config_energy[i][j][(3-1)*5+0];
            // if(bottom_occ == 0 && top_occ == 0) continue;
            // if(2000 < atom.x && atom.x < 3400 && 4750 < atom.y  && atom.y < 5900) {
               config_output << atom.S << ", " << atom.l_id << ", " <<  atom.h_id << ", " \
               << atom.x << ", " << atom.y << ", " \
               << atom.intra1_count << ", " << atom.J_intra1 << ", " << atom.Dx_intra1 << ", " << atom.Dy_intra1 << ", " << atom.Dz_intra1 << ", " \
               << atom.intra2_count << ", " << atom.J_intra2 << ", " << atom.Dx_intra2 << ", " << atom.Dy_intra2 << ", " << atom.Dz_intra2 << ", " \
               << atom.intra3_count << ", " << atom.J_intra3 << ", " << atom.Dx_intra3 << ", " << atom.Dy_intra3 << ", " << atom.Dz_intra3 << ", " \
               << atom.inter1_count << ", " << atom.J_inter1 << ", " << atom.Dx_inter1 << ", " << atom.Dy_inter1 << ", " << atom.Dz_inter1 << ", " \
               << atom.inter2_count << ", " << atom.J_inter2 << ", " << atom.Dx_inter2 << ", " << atom.Dy_inter2 << ", " << atom.Dz_inter2 << ", " \
               << atom.inter3_count << ", " << atom.J_inter3 << ", " << atom.Dx_inter3 << ", " << atom.Dy_inter3 << ", " << atom.Dz_inter3 << ", " \
               << atom.inter_twist1_count << ", " << atom.J_inter_twist1 << ", " << atom.Dx_inter_twist1 << ", " << atom.Dy_inter_twist1 << ", " << atom.Dz_inter_twist1 << ", " \
               << atom.inter_twist2_count << ", " << atom.J_inter_twist2 << ", " << atom.Dx_inter_twist2 << ", " << atom.Dy_inter_twist2 << ", " << atom.Dz_inter_twist2 << ", " \
               << atom.inter_twist3_count << ", " << atom.J_inter_twist3 << ", " << atom.Dx_inter_twist3 << ", " << atom.Dy_inter_twist3 << ", " << atom.Dz_inter_twist3 << std::endl;
               // for(int k = 0; k < config_energy[i][j].size(); k++) config_output << ", " << config_energy[i][j][k]; 
               // config_output << "\n";
            // }
      }
      config_output.close();
      std::cout << "config atoms done." << std::endl;
      std::cout << "config cells started..." << std::flush;

      std::ofstream config_output1(std::string(directory) + "/config_energy_cells.txt");
      if(!config_output1.is_open()) {std::cout << "config energy did not open" << std::endl; exit(1);}
      for(int i = 0; i < global_config_energy.size(); i++) {
         for(int j = 0; j < global_config_energy[i].size(); j++){
               config_output1 << i << ", " << j << ", ";
               // config_output << all_m_atoms[i].S << ", " << all_m_atoms[i].l_id << ", " <<  all_m_atoms[i].h_id << ", " << all_m_atoms[i].x << ", " << all_m_atoms[i].y << ", " << all_m_atoms[i].J_inter << ", " << all_m_atoms[i].Dx_inter << ", " <<  all_m_atoms[i].Dy_inter << ", " << all_m_atoms[i].Dz_inter << ", " << all_m_atoms[i].inter_count  << '\n';// << bottom_occ<< ", " << top_occ;
               for(int k = 0; k < global_config_energy[i][j].size(); k++) config_output1 << global_config_energy[i][j][k] << ", "; 
               config_output1 << "\n";
            // config_output1 << "\n";
         }
         config_output1 << "\n";
      }
      config_output1.close();
      std::cout << "config cells done." << std::endl;
      // // std::cout << "interaction counts started..." << std::flush;
      // std::ofstream interaction_counts(std::string(directory) + "/interaction_counts.txt");
      // if(!interaction_counts.is_open()) {std::cout << "interaction counts did not open" << std::endl; exit(1);}
      // // for(int i = 0; i < all_m_atoms.size(); i++){
      //    interaction_counts << all_m_atoms[i].S  << ", " << all_m_atoms[i].x << ", " << all_m_atoms[i].y <<  ", " << all_m_atoms[i].l_id << ", " <<  all_m_atoms[i].h_id << ", " << all_m_atoms[i].inter1 << ", " << all_m_atoms[i].inter2 << ", " << all_m_atoms[i].inter3 \
      //                                              << ", " << all_m_atoms[i].intra1 << ", " << all_m_atoms[i].intra2 << ", " << all_m_atoms[i].intra3 <<"\n";
      // }
      // interaction_counts.close();
      // outfile4 << ss.str();
      // std::cout << "interaction counts done." << std::endl;
      return;
}

static void intra_pair_shifts_and_angles(spin &central_atom, spin &j_atom,
                                        double angle_i, double angle_j,
                                        int &sx_a, int &sy_a, int &sx_b, int &sy_b,
                                        double &ang_ab, bool &d_flip){
   // Canonical unordered pair (lower original_id → higher) so i→j and j→i share J
   // and opposite DMI regardless of tile/image bin noise.
   spin *a = &central_atom, *b = &j_atom;
   ang_ab = angle_i;
   d_flip = false;
   if(j_atom.original_id < central_atom.original_id){
      a = &j_atom; b = &central_atom;
      ang_ab = angle_j; // j→central == lower→higher after swap
      d_flip = true;
   }
   sx_a = 66; sy_a = 0; sx_b = 66; sy_b = 0;
   if(central_atom.h_id == 1){
      intra_shift_at_atom(*a, sx_a, sy_a);
      intra_shift_at_atom(*b, sx_b, sy_b);
   }
}

static std::array<float,4> match_intra_exchange(int shell, double angle_i, double angle_j,
                                          spin &central_atom, spin &j_atom,
                                          std::vector<std::vector< double > > &Eij_c,
                                          std::vector<std::vector< double > > &Eij_n){
   std::array<float,4> exchange({0.0f,0.0f,0.0f,0.0f});

   int i_x_shift = 66, i_y_shift = 0, j_x_shift = 66, j_y_shift = 0;
   double ang_ab = angle_i;
   bool d_flip = false;
   intra_pair_shifts_and_angles(central_atom, j_atom, angle_i, angle_j,
                                i_x_shift, i_y_shift, j_x_shift, j_y_shift, ang_ab, d_flip);
   const double ang_ba = ang_ab + M_PI;
   std::vector<std::vector<double> > &E_a = d_flip ? Eij_n : Eij_c;
   std::vector<std::vector<double> > &E_b = d_flip ? Eij_c : Eij_n;

   int theta_i, theta_j;
   if(shell == 2) {
      theta_i = classify_honeycomb_2nn_angle(ang_ab, "2NN", central_atom, j_atom);
      theta_j = classify_honeycomb_2nn_angle(ang_ba, "2NN", central_atom, j_atom);
   } else {
      const char *tag = (shell == 1) ? "1NN" : "3NN";
      theta_i = classify_honeycomb_nn_angle(ang_ab, tag, central_atom, j_atom);
      theta_j = classify_honeycomb_nn_angle(ang_ba, tag, central_atom, j_atom);
   }

   double Ji, Dxi, Dyi, Dzi, Jj, Dxj, Dyj, Dzj;
   lookup_intra_map(E_a, i_x_shift, i_y_shift, theta_i, Ji, Dxi, Dyi, Dzi);
   lookup_intra_map(E_b, j_x_shift, j_y_shift, theta_j, Jj, Dxj, Dyj, Dzj);
   const double Javg = 0.5 * (Ji + Jj);
   double Dx_map = 0.5 * (Dxi - Dxj);
   double Dy_map = 0.5 * (Dyi - Dyj);
   double Dz_map = 0.5 * (Dzi - Dzj);
   if(d_flip){ Dx_map = -Dx_map; Dy_map = -Dy_map; Dz_map = -Dz_map; }

   exchange[0] = (float)(J_intra_reduction * Javg);
   const int S = central_atom.S;
   const bool bottom = (S == 1 || S == 2);
   const double c = bottom ? g_c_bot : g_c_top;
   const double s = bottom ? g_s_bot : g_s_top;

   double Dx = Dx_map, Dy = Dy_map, Dz = Dz_map;
   if((S == 1 || S == 4) && Dx_substrate != 0.0) {
      const double ex_vec_x = std::cos(angle_i);
      const double ex_vec_y = std::sin(angle_i);
      const double DMI = Dx_substrate * exchange[0];
      if(S == 1) {
         Dx += DMI*(DMI_sub_vector_y*0.0 - ex_vec_y*(-DMI_sub_vector_z));
         Dy += DMI*((-DMI_sub_vector_z)*ex_vec_x - 0.0*DMI_sub_vector_x);
         Dz += DMI*(DMI_sub_vector_x*ex_vec_y - ex_vec_x*DMI_sub_vector_y);
      } else {
         Dx += DMI*(DMI_sub_vector_y*0.0 - ex_vec_y*DMI_sub_vector_z);
         Dy += DMI*(DMI_sub_vector_z*ex_vec_x - 0.0*DMI_sub_vector_x);
         Dz += DMI*(DMI_sub_vector_x*ex_vec_y - ex_vec_x*DMI_sub_vector_y);
      }
   }
   exchange[1] = (float)(Dx*c - Dy*s);
   exchange[2] = (float)(Dx*s + Dy*c);
   exchange[3] = (float)Dz;
   return exchange;
}

std::array<float,4> match_intra1_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom,
                                          std::vector<std::vector< double > > &Eij_c,
                                          std::vector<std::vector< double > > &Eij_n){
   return match_intra_exchange(1, angle_i, angle_j, central_atom, j_atom, Eij_c, Eij_n);
}

std::array<float,4> match_intra2_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom,
                                          std::vector<std::vector<  double > > &Eij_c,
                                          std::vector<std::vector<  double > > &Eij_n){
   return match_intra_exchange(2, angle_i, angle_j, central_atom, j_atom, Eij_c, Eij_n);
}

std::array<float,4> match_intra3_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom,
                                          std::vector<std::vector< double > > &Eij_c,
                                          std::vector<std::vector< double > > &Eij_n){
   return match_intra_exchange(3, angle_i, angle_j, central_atom, j_atom, Eij_c, Eij_n);
}

std::array<float,4> match_inter_exchange(int /*atom_id*/, int /*nn_id*/, double dx, double dy, double /*dr*/,
   std::vector<std::vector<double> > &Eij_i, std::vector<std::vector<double> > &Eij_j){
   std::array<float,4> exchange({0.0f,0.0f,0.0f,0.0f});
   const int i1 = inter_map_index(dx, dy);
   const int i2 = inter_map_index(-dx, -dy);

   auto centrosym = [&](std::vector<std::vector<double> > &Eij,
                        double &J, double &Dx, double &Dy, double &Dz){
      const auto &r1 = Eij[i1];
      const auto &r2 = Eij[i2];
      J  = 0.5 * (r1[2] + r2[2]);
      Dx = 0.5 * (r1[3] - r2[3]);
      Dy = 0.5 * (r1[4] - r2[4]);
      Dz = 0.5 * (r1[5] - r2[5]);
   };

   double Ji, Dxi, Dyi, Dzi, Jj, Dxj, Dyj, Dzj;
   centrosym(Eij_i, Ji, Dxi, Dyi, Dzi);
   centrosym(Eij_j, Jj, Dxj, Dyj, Dzj);
   exchange[0] = (float)(0.5 * (Ji + Jj));
   exchange[1] = (float)(0.5 * (Dxi + Dxj));
   exchange[2] = (float)(0.5 * (Dyi + Dyj));
   exchange[3] = (float)(0.5 * (Dzi + Dzj));
   return exchange;
}

void spin_config_energy(spin & atom_i, double dr2, spin & atom_j, std::array<float, 4> &exchange, std::vector<std::vector<std::vector<double> > > & local_config_energy) {

   double parity = 1.0;
   //if(atom_i.l_id == 1 || atom_i.l_id == 3) parity = -1.0;
   if(atom_i.S == atom_j.S) {
      if(dr2 < intra_nn_dist_1) {
         all_m_atoms[atom_i.id].intra1_count++;
         all_m_atoms[atom_i.id].J_intra1 += exchange[0];
         all_m_atoms[atom_i.id].Dx_intra1 += exchange[1];
         all_m_atoms[atom_i.id].Dy_intra1 += exchange[2];
         all_m_atoms[atom_i.id].Dz_intra1 += exchange[3];
      } else if (dr2 < intra_nn_dist_2) {
         all_m_atoms[atom_i.id].intra2_count++;
         all_m_atoms[atom_i.id].J_intra2 += exchange[0];
         all_m_atoms[atom_i.id].Dx_intra2 += exchange[1];
         all_m_atoms[atom_i.id].Dy_intra2 += exchange[2];
         all_m_atoms[atom_i.id].Dz_intra2 += exchange[3];
      } else  {
         all_m_atoms[atom_i.id].intra3_count++;
         all_m_atoms[atom_i.id].J_intra3 += exchange[0];
         all_m_atoms[atom_i.id].Dx_intra3 += exchange[1];
         all_m_atoms[atom_i.id].Dy_intra3 += exchange[2];
         all_m_atoms[atom_i.id].Dz_intra3 += exchange[3];
      }
      local_config_energy.at(atom_i.unit_x_lr).at(atom_i.unit_y_lr).at((atom_i.S-1)*16 + 1)++;
      local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 4] += exchange[0];
      local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 5] += fabs(exchange[1]);
      local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 6] += fabs(exchange[2]);
      local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 7] += fabs(exchange[3]);
      if(exchange[0] < 0.0) local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 14] -= exchange[0];
      else local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 15] += exchange[0];
   } else {
      if(atom_i.h_id == atom_j.h_id) {
         if(dr2 <= inter_nn_dist_1) {
            all_m_atoms[atom_i.id].inter_twist1_count++;
            all_m_atoms[atom_i.id].J_inter_twist1 += exchange[0];
            all_m_atoms[atom_i.id].Dx_inter_twist1 += exchange[1];
            all_m_atoms[atom_i.id].Dy_inter_twist1 += exchange[2];
            all_m_atoms[atom_i.id].Dz_inter_twist1 += exchange[3];
         } else if (dr2 <= inter_nn_dist_2) {
            all_m_atoms[atom_i.id].inter_twist2_count++;
            all_m_atoms[atom_i.id].J_inter_twist2 += exchange[0];
            all_m_atoms[atom_i.id].Dx_inter_twist2 += exchange[1];
            all_m_atoms[atom_i.id].Dy_inter_twist2 += exchange[2];
            all_m_atoms[atom_i.id].Dz_inter_twist2 += exchange[3];
         } else {
            all_m_atoms[atom_i.id].inter_twist3_count++;
            all_m_atoms[atom_i.id].J_inter_twist3 += exchange[0];
            all_m_atoms[atom_i.id].Dx_inter_twist3 += exchange[1];
            all_m_atoms[atom_i.id].Dy_inter_twist3 += exchange[2];
            all_m_atoms[atom_i.id].Dz_inter_twist3 += exchange[3];
         }
         local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 2]++;
         local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 8] += exchange[0];
         local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 9] += fabs(exchange[1]);
         local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 10] += fabs(exchange[2]);
         local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 11] += fabs(exchange[3]);
         if(exchange[0] < 0.0) local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 12] -= exchange[0];
         else local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 13] += exchange[0];
      } else {
         if(dr2 <= inter_AB_dist_1) {
            all_m_atoms[atom_i.id].inter1_count++;
            all_m_atoms[atom_i.id].J_inter1 += exchange[0];
            all_m_atoms[atom_i.id].Dx_inter1 += exchange[1];
            all_m_atoms[atom_i.id].Dy_inter1 += exchange[2];
            all_m_atoms[atom_i.id].Dz_inter1 += exchange[3];
         } else if (dr2 <= inter_AB_dist_2) {
            all_m_atoms[atom_i.id].inter2_count++;
            all_m_atoms[atom_i.id].J_inter2 += exchange[0];
            all_m_atoms[atom_i.id].Dx_inter2 += exchange[1];
            all_m_atoms[atom_i.id].Dy_inter2 += exchange[2];
            all_m_atoms[atom_i.id].Dz_inter2 += exchange[3];
         } else {
            all_m_atoms[atom_i.id].inter3_count++;
            all_m_atoms[atom_i.id].J_inter3 += exchange[0];
            all_m_atoms[atom_i.id].Dx_inter3 += exchange[1];
            all_m_atoms[atom_i.id].Dy_inter3 += exchange[2];
            all_m_atoms[atom_i.id].Dz_inter3 += exchange[3];
         }
         local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 3]++;
         // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 12] += exchange[0];
      }      
   }
   
         // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 15] += fabs(exchange[3]);
}

