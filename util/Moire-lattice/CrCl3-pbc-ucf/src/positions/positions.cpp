#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <omp.h>
#include "positions.hpp"
#include "initialise.hpp"
#include "exchange.hpp"


// CrCl3 hexagonal lattice (a≈5.94 Å): McGuire et al., Phys. Rev. Materials 1, 014001 (2017)
// / Morosin & Narath, R-3 at 225 K a=5.942 Å, c=17.333 Å. In-plane Cr at 2/3 cell sites.
double a0x = 5.94;
double a0y= 0.0;
double a1x = -2.97;
double a1y = 5.144191;

// Four Cr planes; interlayer ≈ R-3 c/3 = 17.333/3 Å (Morosin & Narath).
// c0 = 4·interlayer ≈ 23.11 Å = UCF Lz (AFM-style 4L); fz = 0, 1/4, 1/2, 3/4.
double c0 = 23.110667;
double a0z = c0/4.0;

int num_atoms = 8;
int num_nm_atoms = 24;

int number_of_unit_cells_x;
int number_of_unit_cells_y;

int num_above_atoms =0;
int num_below_atoms =0;

double J_inter_scaling = 0.0;
double J_twist_reduction = 1.0;
double J_intra_reduction = 1.0;
double J_prist_reduction = 1.0;
double DMI_inter_scaling = 1.0;
double DMI_sub_scaling = 0.0;
double DMI_sub_vector_x = 0.0;
double DMI_sub_vector_y = 0.0;
double DMI_sub_vector_z = 1.0;

uint64_t total_atoms = 0;
int total_nm_atoms = 0;


bool inside_system(double sx, double sy, double x, double y, double offset){

   if ((x >=offset) && (x <= sx-offset) && (y >=offset) && (y <= sy-offset)) return true;
   else return false;
}


// Variant that first builds a central window [-system_size/5, system_size/5] in x,y,
// then recenters coordinates so all atoms lie in [0, system_size/ (something < 1)] with no negative x or y.
// Unit-cell indices (unit_x_lr/unit_y_lr and unit_x/unit_y) and global_config_energy are filled
// only after the offset is applied, so all_m_atoms stores the offset coordinates.
void create_magnetic_atom_list_central(std::string filename){
   std::cout << "Generating lattice structure in central window...." << std::flush;

   const double x_min = -system_size_x / 2.0;
   const double x_max =  system_size_x / 2.0;
   const double y_min = -system_size_y / 2.0;
   const double y_max =  system_size_y / 2.0;

   // Remember starting index so we only post-process atoms created in this call.
   int start_index = 0; //all_m_atoms.size();

   for (int i = -3*number_of_unit_cells_x; i < 3*number_of_unit_cells_x; i++) {
      for (int j = -3*number_of_unit_cells_y; j < 3*number_of_unit_cells_y; j++){
         for (int atom_i = 0; atom_i < num_atoms; atom_i ++){

            double x_j = atom[atom_i].x + i*a0x + j*a1x;
            double y_j = atom[atom_i].y         + j*a1y;
            double z_j = atom[atom_i].z;

            if ( z_j > twist_loction){
               // upper (twisted) layers: rotate by +twist_angle/2 around origin
               double x_new = x_j*cos(twist_angle*0.5) - y_j*sin(twist_angle*0.5);
               double y_new = y_j*cos(twist_angle*0.5) + x_j*sin(twist_angle*0.5);

               // keep only atoms in central window
               if (x_new >= x_min && x_new <= x_max && y_new >= y_min && y_new <= y_max){
                  spin new_atom;
                  new_atom.x = x_new;
                  new_atom.y = y_new;
                  new_atom.z = z_j;
                  new_atom.id = total_atoms;
                  new_atom.l_id = atom[atom_i].l_id;
                  new_atom.h_id = atom[atom_i].h_id;

                  // Layer assignment based on z (epsilon for file vs a0z rounding)
                  const double zeps = 1e-3;
                  if (z_j <= a0z*2 + zeps){
                     new_atom.S = 3;
                     new_atom.dx = 66;
                     new_atom.dy = 0;
                  } else if (z_j <= a0z*3 + zeps){
                     new_atom.S = 4;
                     new_atom.dx = 66;
                     new_atom.dy = 0;
                  } else if (z_j == 600.0) {
                     new_atom.S = 5;
                  } else {
                     std::cerr << "Error! Atom " << total_atoms << " twist layer: " << z_j << " < " << twist_loction << std::endl;
                     exit(1);
                  }

                  // microcell indices and global_config_energy are filled after recentering
                  new_atom.unit_x_lr = 0;
                  new_atom.unit_y_lr = 0;
                  new_atom.unit_x = 0;
                  new_atom.unit_y = 0;

                  total_atoms++;
                  all_m_atoms.push_back(new_atom);
                  num_above_atoms++;
               }
            } else {
               // lower (twisted) layers: rotate by -twist_angle/2 around origin
               double x_new = x_j*cos(-twist_angle*0.5) - y_j*sin(-twist_angle*0.5);
               double y_new = y_j*cos(-twist_angle*0.5) + x_j*sin(-twist_angle*0.5);

               if (x_new > x_min && x_new < x_max && y_new > y_min && y_new < y_max){
                  spin new_atom;
                  new_atom.x = x_new;
                  new_atom.y = y_new;
                  new_atom.z = z_j;
                  new_atom.id = total_atoms;
                  new_atom.l_id = atom[atom_i].l_id;
                  new_atom.h_id = atom[atom_i].h_id;

                  // Layer assignment based on z (epsilon for file vs a0z rounding)
                  const double zeps = 1e-3;
                  if (z_j <= zeps){
                     new_atom.S = 1;
                     new_atom.dx = 66;
                     new_atom.dy = 0;
                  } else if (z_j <= a0z + zeps){
                     new_atom.S = 2;
                     new_atom.dx = 66;
                     new_atom.dy = 0;
                  } else {
                     std::cerr << "Error! Atom " << total_atoms << " twist layer: " << z_j << " > " << twist_loction << std::endl;
                     exit(1);
                  }

                  new_atom.unit_x_lr = 0;
                  new_atom.unit_y_lr = 0;
                  new_atom.unit_x = 0;
                  new_atom.unit_y = 0;

                  total_atoms++;
                  all_m_atoms.push_back(new_atom);
                  num_below_atoms++;
               }
            }
         }
      } // j-loop
   } // i-loop

   // Recentre so all coordinates are non-negative, then fill microcell indices and global_config_energy.
   const double offset_x = system_size_x / 2.0;
   const double offset_y = system_size_y / 2.0;

   for(int idx = 0; idx < all_m_atoms.size(); ++idx){
      spin &at = all_m_atoms[idx];
      at.x += offset_x;
      at.y += offset_y;

      if(at.x < 0.0 || at.y < 0.0) {
         std::cerr << "Atom " << idx << " has negative coordinates: x=" << at.x << " y=" << at.y << std::endl;
         std::exit(1);
      }
      int dx_cell = microcell_index_x(at.x);
      int dy_cell = microcell_index_y(at.y);
      at.unit_x_lr = dx_cell;
      at.unit_y_lr = dy_cell;
      at.unit_x = dx_cell;
      at.unit_y = dy_cell;

      
      if(at.unit_x_lr < 0 || at.unit_x_lr > microcell_Nx || at.unit_y_lr < 0 || at.unit_y_lr > microcell_Ny) {
         std::cerr << "microcell index out of range (central): unit_x_lr=" << at.unit_x_lr << " unit_y_lr=" << at.unit_y_lr
                   << " (max " << microcell_Nx << "," << microcell_Ny << ") x=" << at.x << " y=" << at.y << std::endl;
         std::exit(1);
      }

      // Fill global_config_energy exactly as in original creator
      if(at.S == 1){
         global_config_energy[at.unit_x_lr][at.unit_y_lr][(at.S-1)*16 + 0] += 1;
      } else if(at.S == 2){
         global_config_energy[at.unit_x_lr][at.unit_y_lr][(at.S-1)*16 + 0] += 1;
      } else if(at.S == 3){
         global_config_energy[at.unit_x_lr][at.unit_y_lr][(at.S-1)*16] += 1;
      } else if(at.S == 4){
         global_config_energy[at.unit_x_lr][at.unit_y_lr][(at.S-1)*16] += 1;
      }
   }

   std::cout << " " << (all_m_atoms.size() - start_index) << " central-window atoms; [complete]" << std::endl;
}


// Find moiré periodicity from global_config_energy, unit_cell_shifts, all_m_atoms.
// Fills out and returns true on success; returns false (and logs to stderr) on failure.
bool find_moire_periodicity(MoirePeriodResult& out, std::vector<MoirePeriodResult>* all_out){
   if(global_config_energy.empty() || global_config_energy[0].empty()) {
      std::cerr << "find_moire_periodicity: global_config_energy not initialised" << std::endl;
      return false;
   }

   const int Nx = static_cast<int>(global_config_energy.size());
   const int Ny = static_cast<int>(global_config_energy[0].size());
   const int shifts_Nx = static_cast<int>(unit_cell_shifts.size());
   const int shifts_Ny = (shifts_Nx > 0 && !unit_cell_shifts[0].empty()) ? static_cast<int>(unit_cell_shifts[0].size()) : 0;

   const int shift_tol = 1;  // tolerance for binned shift values (0–99 scale): descriptor match and AA selection (0, 0)

   auto occ = [&](int i, int j, int layer)->double {
      // layer in [1,4], stored at (layer-1)*16
      return global_config_energy[i][j][(layer-1)*16];
   };

   auto descriptors_equal = [&](int i1, int j1, int i2, int j2)->bool {
      // Per-layer occupancy: same spin count in each layer at the two microcells (no averaging across layers).
      // double n1a = occ(i1,j1,1), n2a = occ(i1,j1,2), n3a = occ(i1,j1,3), n4a = occ(i1,j1,4);
      // double n1b = occ(i2,j2,1), n2b = occ(i2,j2,2), n3b = occ(i2,j2,3), n4b = occ(i2,j2,4);
      // const double occ_tol = 0.0;
      // if(!(std::fabs(n1a-n1b) <= occ_tol &&
      //      std::fabs(n2a-n2b) <= occ_tol &&
      //      std::fabs(n3a-n3b) <= occ_tol &&
      //      std::fabs(n4a-n4b) <= occ_tol)) return false;

      // Registry: unit_cell_shifts must match at the two microcells (same binned dx, dy), within shift_tol.
      // Sliding indices are periodic on 0..99 — compare toroidally (0≈99).
      bool in1 = (i1 >= 0 && i1 < shifts_Nx-1 && j1 >= 0 && j1 < shifts_Ny-1);
      bool in2 = (i2 >= 0 && i2 < shifts_Nx-1 && j2 >= 0 && j2 < shifts_Ny-1);
      if(in1 != in2) return false; // one in bounds, one out -> not equal
      if(in1 && in2) {
         if(unit_cell_shifts[i1][j1][0] < 1 || unit_cell_shifts[i2][j2][0] < 1) return false;
         auto tor_d = [](int a, int b){
            int d = std::abs(a - b) % 100;
            return std::min(d, 100 - d);
         };
         const int d1 = tor_d(unit_cell_shifts[i1][j1][1], unit_cell_shifts[i2][j2][1]);
         const int d2 = tor_d(unit_cell_shifts[i1][j1][2], unit_cell_shifts[i2][j2][2]);
         if(d1 <= shift_tol && d2 <= shift_tol) return true;
      }
      return false;
   };

   // Restrict period detection to interior: 2 cells from each edge (indices 2 .. Nx-3, 2 .. Ny-3).
   const int edge_margin = 2;
   const int i_min = edge_margin;
   const int i_max = Nx - 1 - edge_margin;   // Nx-3
   const int j_min = edge_margin;
   const int j_max = Ny - 1 - edge_margin;   // Ny-3

   const double cell_ax = microcell_ax();
   const double cell_ay = microcell_ay();

   // Minimum dx, dy from analytic moiré period L = a/(2*sin(θ/2)) so we target the true AA–AA distance (~360 Å at 1.1°, a=6.9 Å), not a shorter sub-repeat (~70 Å).
   // Guard θ→0 (sin(θ/2)→0): period is formally infinite; bake a finite orthogonal cell instead of detecting.
   const double sin_half = std::sin(twist_angle / 2.0);
   if(std::fabs(sin_half) < 1e-12){
      std::cerr << "find_moire_periodicity: twist angle ~0; analytic moiré period diverges. "
                << "Bake a finite orthogonal cell (moire_area) rather than auto-detecting." << std::endl;
      return false;
   }
   const double L_moire_angstrom = a0x / (2.0 * sin_half);
   const int expected_cells_x = static_cast<int>(std::floor(L_moire_angstrom / cell_ax + 0.5));
   const int expected_cells_y = static_cast<int>(std::floor(L_moire_angstrom / cell_ay + 0.5));
   const int min_dx_period = std::max(2, static_cast<int>(expected_cells_x * 0.25));
   const int min_dy_period = std::max(2, static_cast<int>(expected_cells_y * 0.25));
   // Cap search so we prefer a primary moiré cell; widen later if nothing passes.
   int max_dx_period = std::max(min_dx_period + 1, static_cast<int>(std::ceil(expected_cells_x * 2.5)));
   int max_dy_period = std::max(min_dy_period + 1, static_cast<int>(std::ceil(expected_cells_y * 2.5)));
   std::cout << "build_spin_moire_cell: analytic moiré period L=" << L_moire_angstrom << " Å; min_dx=" << min_dx_period
             << " max_dx=" << max_dx_period << " min_dy=" << min_dy_period << " max_dy=" << max_dy_period
             << " (expected " << expected_cells_x << "x" << expected_cells_y << " cells)." << std::endl;

   // --- Step 1: AA cells from unit_cell_shifts (AA = shift near 0,0; empty/AB default is 66,0) ---
   const int aa_shift_x = 0, aa_shift_y = 0;
   std::vector<std::pair<int,int>> aa_cells;
   for(int i = i_min; i <= i_max; ++i){
      for(int j = j_min; j <= j_max; ++j){
         if(unit_cell_shifts[i][j][0] > 0 &&
            std::abs(unit_cell_shifts[i][j][1] - aa_shift_x) <= shift_tol &&
            std::abs(unit_cell_shifts[i][j][2] - aa_shift_y) <= shift_tol)
            aa_cells.push_back({i, j});
      }
   }
   if(aa_cells.empty()){
      std::cerr << "find_moire_periodicity: no AA-configuration cells found (shift near " << aa_shift_x << "," << aa_shift_y << " ±" << shift_tol << ")." << std::endl;
      return false;
   } 
   std::cout << "build_spin_moire_cell: Step 1 complete: " << aa_cells.size() << " AA cells found. ";

   // Atoms per cell for atom-level check (indexed by unit_x_lr, unit_y_lr).
   std::vector<std::vector<std::vector<const spin*>>> atoms_in_cell(shifts_Nx);
   int atoms_in_cell_count = 0;
   for(int i = 0; i < shifts_Nx; ++i) atoms_in_cell[i].resize(shifts_Ny);
   for(size_t k = 0; k < all_m_atoms.size(); ++k){
      const spin& at = all_m_atoms[k];
      // if(at.S != 2 && at.S != 3) continue;
      int ui = at.unit_x_lr, uj = at.unit_y_lr;
      if(ui >= 2 && ui < shifts_Nx-2 && uj >= 2 && uj < shifts_Ny-2) {
         atoms_in_cell[ui][uj].push_back(&at);
         atoms_in_cell_count++;
      }
   }
   std::cout <<  atoms_in_cell_count << " atoms in cells. " << std::endl;
   // Bake-only: slightly looser atom match helps close quartets at awkward angles (e.g. 1.8°/1.9°).
   const double atom_match_tol = bake_only ? 0.35 : 1e-1;

   int Px = -1, Py = -1, i0 = -1, j0 = -1;
   int cell_count = 0;
   int atom_count = 0;

   struct CandidateCell {
      int i0, j0;
      int Px, Py;
   };
   std::vector<CandidateCell> candidate_cells;
   candidate_cells.reserve(aa_cells.size());

   std::cout << "build_spin_moire_cell: Step 2: scanning " << aa_cells.size() << " AA cells for period (Px,Py) and edge-periodicity check." << std::endl;
   // --- Step 2: Parallel loop over AA cells to build list of descriptor-periodic candidates. ---
   #pragma omp parallel num_threads(4)
   {
      std::vector<CandidateCell> my_candidates;

      #pragma omp for schedule(dynamic) reduction(+:cell_count) reduction(+:atom_count)
      for(size_t idx = 0; idx < aa_cells.size(); ++idx){
         const int i0_cand = aa_cells[idx].first, j0_cand = aa_cells[idx].second;

         // Scan over possible dx,dy periods; for each pair that passes edge checks,
         // add a candidate but keep searching for larger periods as well.
         for(int dx = min_dx_period; dx <= max_dx_period && i0_cand + dx <= i_max; ++dx){
            if(!descriptors_equal(i0_cand, j0_cand, i0_cand + dx, j0_cand)) continue;

            for(int dy = min_dy_period; dy <= max_dy_period && j0_cand + dy <= j_max; ++dy){
               if(!descriptors_equal(i0_cand, j0_cand, i0_cand, j0_cand + dy)) continue;

               const int Px_cand = dx;
               const int Py_cand = dy;

               // Edge-periodicity: descriptor must repeat across the cell boundaries so the crop is periodic.
               // Left edge (column i0_cand) must match the next period's left edge (column i0_cand + Px_cand) for every row j in the cell.
               // Bottom edge (row j0_cand) must match the next period's bottom edge (row j0_cand + Py_cand) for every column i in the cell.
               bool edge_ok = true;
               for(int j = j0_cand; edge_ok && j < j0_cand + Py_cand; ++j)
                  if(!descriptors_equal(i0_cand, j, i0_cand + Px_cand, j)) edge_ok = false;
               for(int i = i0_cand; edge_ok && i < i0_cand + Px_cand; ++i)
                  if(!descriptors_equal(i, j0_cand, i, j0_cand + Py_cand)) edge_ok = false;
               if(!edge_ok) continue;

               // Descriptor-level criteria passed: record this AA cell + (Px_cand,Py_cand) as a candidate.
               cell_count++;
               my_candidates.push_back({i0_cand, j0_cand, Px_cand, Py_cand});
            }
         }
      }

      #pragma omp critical
      {
         candidate_cells.insert(candidate_cells.end(), my_candidates.begin(), my_candidates.end());
      }
   }

   if(candidate_cells.empty()){
      // Widen period window once — some angles only match at multi-period coincidence cells.
      max_dx_period = std::max(max_dx_period + 1, static_cast<int>(std::ceil(expected_cells_x * 6.0)));
      max_dy_period = std::max(max_dy_period + 1, static_cast<int>(std::ceil(expected_cells_y * 6.0)));
      std::cout << "find_moire_periodicity: no primary-scale candidates; widening to max_dx="
                << max_dx_period << " max_dy=" << max_dy_period << std::endl;
      #pragma omp parallel num_threads(4)
      {
         std::vector<CandidateCell> my_candidates;
         #pragma omp for schedule(dynamic) reduction(+:cell_count)
         for(size_t idx = 0; idx < aa_cells.size(); ++idx){
            const int i0_cand = aa_cells[idx].first, j0_cand = aa_cells[idx].second;
            for(int dx = min_dx_period; dx <= max_dx_period && i0_cand + dx <= i_max; ++dx){
               if(!descriptors_equal(i0_cand, j0_cand, i0_cand + dx, j0_cand)) continue;
               for(int dy = min_dy_period; dy <= max_dy_period && j0_cand + dy <= j_max; ++dy){
                  if(!descriptors_equal(i0_cand, j0_cand, i0_cand, j0_cand + dy)) continue;
                  const int Px_cand = dx;
                  const int Py_cand = dy;
                  bool edge_ok = true;
                  for(int j = j0_cand; edge_ok && j < j0_cand + Py_cand; ++j)
                     if(!descriptors_equal(i0_cand, j, i0_cand + Px_cand, j)) edge_ok = false;
                  for(int i = i0_cand; edge_ok && i < i0_cand + Px_cand; ++i)
                     if(!descriptors_equal(i, j0_cand, i, j0_cand + Py_cand)) edge_ok = false;
                  if(!edge_ok) continue;
                  cell_count++;
                  my_candidates.push_back({i0_cand, j0_cand, Px_cand, Py_cand});
               }
            }
         }
         #pragma omp critical
         {
            candidate_cells.insert(candidate_cells.end(), my_candidates.begin(), my_candidates.end());
         }
      }
   }

   if(candidate_cells.empty()){
      std::cerr << "find_moire_periodicity: no AA candidates passed descriptor edge checks (Step 2). Cells passed to diagonal: " << cell_count << std::endl;
      return false;
   }

   // Step 3: atomic quartets defining orthogonal (possibly rotated) translation vectors Tx,Ty.
   // Gate on |cos θ| only — rotated-orthogonal is allowed; true shear is not.
   // Each quartet uses atoms from (i0,j0), (i0+Px,j0), (i0,j0+Py), (i0+Px,j0+Py).
   struct QuartetResult {
      double area;
      double cos_abs;  // |Tx·Ty|/(|Tx||Ty|); ideal 0
      double Lx, Ly;
      double tx_x, tx_y;  // x-period translation vector (atom a -> atom b)
      double ty_x, ty_y;  // y-period translation vector (atom a -> atom c)
      double x0, y0;
      int i0, j0, Px, Py;
   };

   bool have_best_quartet = false;
   QuartetResult best_q{};
   std::vector<QuartetResult> all_quartets;

   // Tight |cos θ| gate: residual Tx/Ty skew shows up as ~1° 1NN angle errors at seams.
   // Bake-only may relax slightly so awkward-angle cells can close; never axis-lock.
   const double orthogonality_tol = bake_only ? 5e-4 : 5e-5; // max |cos(theta)|

   #pragma omp parallel num_threads(4)
   {
      std::vector<QuartetResult> local_quartets;

      #pragma omp for schedule(dynamic)
      for(size_t ci = 0; ci < candidate_cells.size(); ++ci){
         const CandidateCell &c = candidate_cells[ci];
         const int i0_cand = c.i0;
         const int j0_cand = c.j0;
         const int Px_cand = c.Px;
         const int Py_cand = c.Py;

         const int cell_x_right = i0_cand + Px_cand;
         const int cell_y_top   = j0_cand + Py_cand;
         const int diag_i       = cell_x_right;
         const int diag_j       = cell_y_top;

         if(cell_x_right <= 0 || cell_x_right >= shifts_Nx) continue;
         if(cell_y_top   <= 0 || cell_y_top   >= shifts_Ny) continue;

         const auto &base_cells  = atoms_in_cell[i0_cand][j0_cand];
         const auto &right_cells = atoms_in_cell[cell_x_right][j0_cand];
         const auto &top_cells   = atoms_in_cell[i0_cand][cell_y_top];
         const auto &diag_cells  = atoms_in_cell[diag_i][diag_j];
         if(base_cells.empty() || right_cells.empty() || top_cells.empty() || diag_cells.empty()) continue;

         for(const spin* a : base_cells){
            for(const spin* b : right_cells){
               if(b->S != a->S || b->l_id != a->l_id) continue;
               const double vx_x = b->x - a->x;
               const double vx_y = b->y - a->y;

               const double Lx_candidate = std::sqrt(vx_x*vx_x + vx_y*vx_y);
               if(Lx_candidate <= 0.0) continue;

               for(const spin* c_top : top_cells){
                  if(c_top->S != a->S || c_top->l_id != a->l_id) continue;
                  const double vy_x = c_top->x - a->x;
                  const double vy_y = c_top->y - a->y;

                  const double Ly_candidate = std::sqrt(vy_x*vy_x + vy_y*vy_y);
                  if(Ly_candidate <= 0.0) continue;

                  // Positive oriented area (right-handed Tx,Ty); allows rotated frames.
                  const double oriented = vx_x*vy_y - vx_y*vy_x;
                  if(oriented <= 0.0) continue;

                  const double dot = vx_x*vy_x + vx_y*vy_y;
                  const double cos_theta = std::fabs(dot / (Lx_candidate * Ly_candidate));
                  if(cos_theta > orthogonality_tol) continue;

                  const double x_diag_target = a->x + vx_x + vy_x;
                  const double y_diag_target = a->y + vx_y + vy_y;

                  bool has_diag = false;
                  for(const spin* d : diag_cells){
                     if(d->S != a->S || d->l_id != a->l_id) continue;
                     if(std::fabs(d->x - x_diag_target) < atom_match_tol &&
                        std::fabs(d->y - y_diag_target) < atom_match_tol){
                        has_diag = true;
                        break;
                     }
                  }
                  if(!has_diag) continue;

                  QuartetResult q{};
                  q.area = oriented;
                  q.cos_abs = cos_theta;
                  q.Lx   = Lx_candidate;
                  q.Ly   = Ly_candidate;
                  q.tx_x = vx_x;
                  q.tx_y = vx_y;
                  q.ty_x = vy_x;
                  q.ty_y = vy_y;
                  q.x0   = a->x;
                  q.y0   = a->y;
                  q.i0   = i0_cand;
                  q.j0   = j0_cand;
                  q.Px   = Px_cand;
                  q.Py   = Py_cand;
                  local_quartets.push_back(q);
               }
            }
         }
      }

      #pragma omp critical
      {
         all_quartets.insert(all_quartets.end(), local_quartets.begin(), local_quartets.end());
      }
   }

   if(all_quartets.empty()){
      std::cerr << "find_moire_periodicity: no atomic quartets found with orthogonal translations and matching diagonal." << std::endl;
      return false;
   }

   std::sort(all_quartets.begin(), all_quartets.end(),
             [&](const QuartetResult &a, const QuartetResult &b){
                // Prefer highest orthogonality (|cos θ|→0), then smaller area.
                if(std::fabs(a.cos_abs - b.cos_abs) > 1e-12)
                   return a.cos_abs < b.cos_abs;
                if(std::fabs(a.area - b.area) > 1e-6)
                   return a.area < b.area;
                // Tie-break: closer to analytic L_moire.
                auto scale_score = [&](const QuartetResult &q){
                   const double rx = q.Lx / L_moire_angstrom;
                   const double ry = q.Ly / L_moire_angstrom;
                   return (rx - 1.0)*(rx - 1.0) + (ry - 1.0)*(ry - 1.0);
                };
                return scale_score(a) < scale_score(b);
             });

   best_q = all_quartets.front();
   have_best_quartet = true;

   // dump full reduced list for diagnostics
   {
      char directory[256];
      if(getcwd(directory, sizeof(directory)) != NULL){
         std::ofstream qout(std::string(directory) + "/moire_quartets.txt");
         qout << "# i0 j0 Px Py Lx Ly tx_x tx_y ty_x ty_y area cos_abs\n";
         for(const auto &q : all_quartets){
            qout << q.i0 << " " << q.j0 << " "
                 << q.Px << " " << q.Py << " "
                 << q.Lx << " " << q.Ly << " "
                 << q.tx_x << " " << q.tx_y << " "
                 << q.ty_x << " " << q.ty_y << " "
                 << q.area << " " << q.cos_abs << "\n";
         }
      }
   }
   Px = best_q.Px;
   Py = best_q.Py;
   i0 = best_q.i0;
   j0 = best_q.j0;
   atom_count = 1; // at least one quartet found
   
   if(Px < 0){
      std::cerr << "find_moire_periodicity: no periodic unit cell found (no AA candidate passed edge and atom-level checks). cell_count=" << cell_count << ", atom_count=" << atom_count << std::endl;
      return false;
   }
   std::cout << "find_moire_periodicity: Steps 2&3 complete: chosen Px=" << Px << " Py=" << Py
             << " x0=" << best_q.x0 << " y0=" << best_q.y0 << "x1=" << best_q.x0 + best_q.Lx << " y1=" << best_q.y0 + best_q.Ly
             << " Lx=" << best_q.Lx << " Ly=" << best_q.Ly << "." << std::endl;

   out.x0 = best_q.x0;
   out.y0 = best_q.y0;
   out.Lx = best_q.Lx;
   out.Ly = best_q.Ly;
   out.tx_x = best_q.tx_x;
   out.tx_y = best_q.tx_y;
   out.ty_x = best_q.ty_x;
   out.ty_y = best_q.ty_y;
   out.Px = Px;
   out.Py = Py;
   out.i0 = i0;
   out.j0 = j0;
   out.from_detection = true;

   if(all_out){
      all_out->clear();
      // Prefer distinct periods; keep up to 3 AA origins per (Px,Py,~L) so seam/origin
      // can be retried without exploding to thousands of quartets.
      std::vector<int> per_period_count;
      for(const auto &q : all_quartets){
         MoirePeriodResult r;
         r.x0 = q.x0; r.y0 = q.y0; r.Lx = q.Lx; r.Ly = q.Ly;
         r.tx_x = q.tx_x; r.tx_y = q.tx_y; r.ty_x = q.ty_x; r.ty_y = q.ty_y;
         r.Px = q.Px; r.Py = q.Py; r.i0 = q.i0; r.j0 = q.j0;
         r.from_detection = true;
         int bucket = -1;
         for(size_t pi = 0; pi < all_out->size(); ++pi){
            const auto &prev = (*all_out)[pi];
            if(prev.Px == r.Px && prev.Py == r.Py &&
               std::fabs(prev.Lx - r.Lx) < 0.5 && std::fabs(prev.Ly - r.Ly) < 0.5){
               // Exact same AA origin → skip
               if(std::fabs(prev.x0 - r.x0) < 0.5 && std::fabs(prev.y0 - r.y0) < 0.5){
                  bucket = -2; break;
               }
               if(bucket < 0) bucket = (int)pi;
            }
         }
         if(bucket == -2) continue;
         if(bucket >= 0){
            // Count how many AA origins already stored for this period.
            int n_same = 0;
            for(const auto &prev : *all_out){
               if(prev.Px == r.Px && prev.Py == r.Py &&
                  std::fabs(prev.Lx - r.Lx) < 0.5 && std::fabs(prev.Ly - r.Ly) < 0.5)
                  ++n_same;
            }
            if(n_same >= 3) continue;
         }
         all_out->push_back(r);
         if(all_out->size() >= 40) break;
      }
      std::cout << "find_moire_periodicity: " << all_out->size() << " unique candidate cells for PBC retry." << std::endl;
      // Deterministic --candidate indices; prefer larger Py (double-cell) then smaller area.
      std::sort(all_out->begin(), all_out->end(),
                [](const MoirePeriodResult &a, const MoirePeriodResult &b){
                   if(a.Py != b.Py) return a.Py > b.Py;
                   if(a.Px != b.Px) return a.Px < b.Px;
                   const double aa = a.Lx * a.Ly, bb = b.Lx * b.Ly;
                   if(std::fabs(aa - bb) > 1.0) return aa < bb;
                   if(std::fabs(a.x0 - b.x0) > 0.25) return a.x0 < b.x0;
                   return a.y0 < b.y0;
                });
   }
   return true;
}

// Build a spin-moiré primary cell from the large patch and perform an atomic-level PBC check.
// If use_given is non-null and has Lx > 0, Ly > 0, period detection is skipped and the given lattice vectors/cutoffs are used.
void build_spin_moire_cell(double moire_area[12], const MoirePeriodResult* use_given){
   if(global_config_energy.empty() || global_config_energy[0].empty()) {
      std::cerr << "build_spin_moire_cell: global_config_energy not initialised" << std::endl;
      return;
   }

   std::vector<MoirePeriodResult> candidates;
   if(moire_area[2] > 0 && moire_area[3] > 0){
      MoirePeriodResult r;
      r.x0 = moire_area[0];
      r.y0 = moire_area[1];
      r.Lx = moire_area[2];
      r.Ly = moire_area[3];
      // Optional baked orthogonal (possibly rotated) Tx/Ty in moire_area[4..7]. All-zero → axis-aligned.
      const double vec_sum = std::fabs(moire_area[4]) + std::fabs(moire_area[5])
                             + std::fabs(moire_area[6]) + std::fabs(moire_area[7]);
      if(vec_sum > 1e-12){
         r.tx_x = moire_area[4]; r.tx_y = moire_area[5];
         r.ty_x = moire_area[6]; r.ty_y = moire_area[7];
         r.Lx = std::sqrt(r.tx_x*r.tx_x + r.tx_y*r.tx_y);
         r.Ly = std::sqrt(r.ty_x*r.ty_x + r.ty_y*r.ty_y);
      } else {
         r.tx_x = r.Lx; r.tx_y = 0.0;
         r.ty_x = 0.0; r.ty_y = r.Ly;
      }
      // Optional baked (Px,Py,i0,j0) in [8..11]; else derive from floor(L/cell).
      if(moire_area[8] > 0.5 && moire_area[9] > 0.5){
         r.Px = (int)std::lround(moire_area[8]);
         r.Py = (int)std::lround(moire_area[9]);
         r.i0 = (int)std::lround(moire_area[10]);
         r.j0 = (int)std::lround(moire_area[11]);
      } else {
         r.Px = static_cast<int>(std::floor((r.Lx+1e-9)/microcell_ax()));
         r.Py = static_cast<int>(std::floor((r.Ly+1e-9)/microcell_ay()));
         r.i0 = microcell_index_x(r.x0);
         r.j0 = microcell_index_y(r.y0);
      }
      r.from_detection = false;
      candidates.push_back(r);
      std::cout << "Using read-in moire area: x0=" << r.x0 << " y0=" << r.y0 << " Lx=" << r.Lx << " Ly=" << r.Ly
                << " tx=(" << r.tx_x << "," << r.tx_y << ") ty=(" << r.ty_x << "," << r.ty_y << ")"
                << " Px=" << r.Px << " Py=" << r.Py << " i0=" << r.i0 << " j0=" << r.j0 << std::endl;
   } else {
      MoirePeriodResult res;
      if(!find_moire_periodicity(res, &candidates)){
         exit(1);
      }
      if(candidates.empty()) candidates.push_back(res);
   }

   {
      char directory[256];
      if(getcwd(directory, sizeof(directory)) != NULL){
         std::ofstream cout_cands(std::string(directory) + "/aa_candidates.txt");
         cout_cands << "# idx Px Py i0 j0 x0 y0 Lx Ly tx_x tx_y ty_x ty_y\n";
         for(size_t i = 0; i < candidates.size(); ++i){
            const auto &c = candidates[i];
            cout_cands << i << " " << c.Px << " " << c.Py << " " << c.i0 << " " << c.j0 << " "
                       << c.x0 << " " << c.y0 << " " << c.Lx << " " << c.Ly << " "
                       << c.tx_x << " " << c.tx_y << " " << c.ty_x << " " << c.ty_y << "\n";
            std::cout << "AA candidate[" << i << "]: Px=" << c.Px << " Py=" << c.Py
                      << " AA=(" << c.x0 << "," << c.y0 << ") |T|=(" << c.Lx << "," << c.Ly << ")"
                      << " tx=(" << c.tx_x << "," << c.tx_y << ") ty=(" << c.ty_x << "," << c.ty_y << ")"
                      << std::endl;
         }
         std::cout << "Wrote " << candidates.size() << " AA candidates to aa_candidates.txt" << std::endl;
      }
   }
   if(moire_force_candidate >= 0){
      if(moire_force_candidate >= (int)candidates.size()){
         std::cerr << "build_spin_moire_cell: --candidate " << moire_force_candidate
                   << " out of range (n=" << candidates.size() << ")" << std::endl;
         exit(1);
      }
      MoirePeriodResult only = candidates[(size_t)moire_force_candidate];
      candidates.clear();
      candidates.push_back(only);
      std::cout << "Using forced AA candidate index " << moire_force_candidate
                << " only (Px=" << only.Px << " Py=" << only.Py << ")" << std::endl;
   }

   // Tile with detected Tx/Ty. Spins replace sites in-place: keep pre-retile
   // unit_x/y and unit_cell_shifts; only original_id (and G) change.
   const std::vector<spin> atoms_backup = all_m_atoms;
   bool pbc_ok = false;
   double x0=0, y0=0, Lx=0, Ly=0, tx_x=0, tx_y=0, ty_x=0, ty_y=0;
   int Px=0, Py=0, i0=0, j0=0;
   std::vector<spin> moire_patch;

   auto in_parallelogram = [](double x, double y, double ox, double oy,
                              double tx_x, double tx_y, double ty_x, double ty_y,
                              double edge_eps)->bool {
      const double rx = x - ox;
      const double ry = y - oy;
      const double det = tx_x * ty_y - tx_y * ty_x;
      if(std::fabs(det) < 1e-18) return false;
      const double a = ( rx * ty_y - ry * ty_x) / det;
      const double b = (-rx * tx_y + ry * tx_x) / det;
      return (a >= -edge_eps) && (a < 1.0 - edge_eps) && (b >= -edge_eps) && (b < 1.0 - edge_eps);
   };

   for(size_t ci = 0; ci < candidates.size(); ++ci){
      const MoirePeriodResult &res = candidates[ci];
      all_m_atoms = atoms_backup;
      x0 = res.x0; y0 = res.y0;
      tx_x = res.tx_x; tx_y = res.tx_y; ty_x = res.ty_x; ty_y = res.ty_y;
      if(std::fabs(tx_x) + std::fabs(tx_y) < 1e-12){ tx_x = res.Lx; tx_y = 0.0; }
      if(std::fabs(ty_x) + std::fabs(ty_y) < 1e-12){ ty_x = 0.0; ty_y = res.Ly; }
      Lx = std::sqrt(tx_x*tx_x + tx_y*tx_y);
      Ly = std::sqrt(ty_x*ty_x + ty_y*ty_y);
      Px = res.Px; Py = res.Py; i0 = res.i0; j0 = res.j0;

      {
         const double cos_abs = std::fabs(tx_x*ty_x + tx_y*ty_y) / (Lx * Ly);
         if(cos_abs > 5e-5){
            std::cerr << "build_spin_moire_cell: candidate " << ci
                      << " rejected: |cos θ|=" << cos_abs << " (true shear)" << std::endl;
            continue;
         }
      }

      // Resize bin pitch so Px·ax' = |Tx|, Py·ay' = |Ty| (removes Tx vs Px·a0x slip).
      if(Px <= 0 || Py <= 0){
         std::cerr << "build_spin_moire_cell: candidate " << ci << " bad Px/Py" << std::endl;
         continue;
      }
      {
         const double ax0 = static_cast<double>(microcell_scale_x) * a0x;
         const double ay0 = static_cast<double>(microcell_scale_y) * a1y;
         microcell_ax_eff = Lx / static_cast<double>(Px);
         microcell_ay_eff = Ly / static_cast<double>(Py);
         std::cout << "bin resize: ax " << ax0 << " -> " << microcell_ax_eff
                   << " (Px*ax'=|Tx|=" << (microcell_ax_eff * Px) << ")"
                   << "  ay " << ay0 << " -> " << microcell_ay_eff
                   << " (Py*ay'=|Ty|=" << (microcell_ay_eff * Py) << ")" << std::endl;
      }
      // Re-bin open lattice on ax',ay'; recompute shifts; refresh i0,j0.
      for(auto &at : all_m_atoms) rebin_atom_microcells(at);
      {
         std::vector<std::vector<std::vector<int>>> shifts_rebin;
         compute_unit_cell_shifts_from_atoms(all_m_atoms, shifts_rebin);
         unit_cell_shifts.swap(shifts_rebin);
      }
      i0 = microcell_index_x(x0);
      j0 = microcell_index_y(y0);

      std::cout << "Trying candidate " << (ci+1) << "/" << candidates.size()
                << ": Px=" << Px << " Py=" << Py << " AA=(" << x0 << "," << y0 << ")"
                << " |tx|=" << Lx << " |ty|=" << Ly
                << " tx=(" << tx_x << "," << tx_y << ") ty=(" << ty_x << "," << ty_y << ")"
                << " i0=" << i0 << " j0=" << j0
                << " ax'=" << microcell_ax() << " ay'=" << microcell_ay() << std::endl;

      // Dump shifts NOW — after PBC id for this candidate, before retile/recompute.
      // Includes whole-grid Δ vs (i+Px,j+Py), (i+Px,j), (i,j+Py) where both cells exist.
      {
         char directory[256];
         if(getcwd(directory, sizeof(directory)) != NULL){
            const int Nx = (int)unit_cell_shifts.size();
            std::ofstream sf(std::string(directory) + "/unit_cell_shifts_pre_tile.txt");
            sf << "# unit_cell_shifts AFTER PBC id, BEFORE retile/recompute\n";
            sf << "# Px Py i0 j0\n";
            sf << Px << " " << Py << " " << i0 << " " << j0 << "\n";
            sf << "# i j occ sx sy  dsx_pp dsy_pp  dsx_px dsy_px  dsx_py dsy_py\n";
            sf << "# pp=vs(i+Px,j+Py)  px=vs(i+Px,j)  py=vs(i,j+Py); blank if partner missing\n";
            for(int ii = 0; ii < Nx; ++ii){
               const int Ny = (int)unit_cell_shifts[ii].size();
               for(int jj = 0; jj < Ny; ++jj){
                  const int occ = unit_cell_shifts[ii][jj][0];
                  const int sx = unit_cell_shifts[ii][jj][1];
                  const int sy = unit_cell_shifts[ii][jj][2];
                  sf << ii << " " << jj << " " << occ << " " << sx << " " << sy;
                  auto emit_delta = [&](int ip, int jp){
                     if(ip < 0 || jp < 0 || ip >= Nx || jp >= (int)unit_cell_shifts[ip].size()){
                        sf << " nan nan";
                        return;
                     }
                     if(occ < 1 || unit_cell_shifts[ip][jp][0] < 1){
                        sf << " nan nan";
                        return;
                     }
                     sf << " " << (sx - unit_cell_shifts[ip][jp][1])
                        << " " << (sy - unit_cell_shifts[ip][jp][2]);
                  };
                  emit_delta(ii + Px, jj + Py);
                  emit_delta(ii + Px, jj);
                  emit_delta(ii, jj + Py);
                  sf << "\n";
               }
            }
            std::cout << "Wrote unit_cell_shifts_pre_tile.txt (full-grid PBC deltas) before retile" << std::endl;
         }
      }

      const double edge_eps = 1e-4;
      const double aa_x = x0, aa_y = y0;
      moire_patch.clear();
      uint64_t new_id = 0;
      // Primary G=0 = half-open AA parallelogram. Keep unit bins; only remap IDs.
      for(size_t k = 0; k < all_m_atoms.size(); ++k){
         spin at = all_m_atoms[k];
         if(!in_parallelogram(at.x, at.y, aa_x, aa_y, tx_x, tx_y, ty_x, ty_y, edge_eps)) continue;
         at.original_id = new_id++;
         moire_patch.push_back(at);
      }
      std::cout << "moire patch size: " << moire_patch.size() << std::endl;
      if(moire_patch.empty()){
         std::cerr << "build_spin_moire_cell: empty patch for candidate " << ci << std::endl;
         continue;
      }

      // Safety dedup under {0,Tx,Ty,Tx+Ty} (S+l_id+h_id).
      {
         const double dup_tol = 0.05;
         const double inv = 1.0 / std::max(dup_tol, 1e-6);
         struct Key { int ix, iy, S, lid, hid; };
         struct KeyHash {
            size_t operator()(Key const& k) const {
               return (size_t)k.ix * 73856093u ^ (size_t)k.iy * 19349663u
                    ^ (size_t)k.S * 83492791u ^ (size_t)k.lid * 39916801u
                    ^ (size_t)k.hid;
            }
         };
         struct KeyEq {
            bool operator()(Key const& a, Key const& b) const {
               return a.ix==b.ix && a.iy==b.iy && a.S==b.S && a.lid==b.lid && a.hid==b.hid;
            }
         };
         std::unordered_map<Key, std::vector<size_t>, KeyHash, KeyEq> bins;
         bins.reserve(moire_patch.size()*2);
         for(size_t k = 0; k < moire_patch.size(); ++k){
            Key key{(int)std::floor(moire_patch[k].x * inv),
                    (int)std::floor(moire_patch[k].y * inv),
                    (int)moire_patch[k].S, moire_patch[k].l_id, moire_patch[k].h_id};
            bins[key].push_back(k);
         }
         auto bary_sum = [&](const spin& at)->double {
            const double rx = at.x - aa_x, ry = at.y - aa_y;
            const double det = tx_x * ty_y - tx_y * ty_x;
            const double a = ( rx * ty_y - ry * ty_x) / det;
            const double b = (-rx * tx_y + ry * tx_x) / det;
            return a + b;
         };
         std::vector<char> drop(moire_patch.size(), 0);
         auto mark_pair = [&](size_t ia, size_t ib){
            if(drop[ia] || drop[ib]) return;
            if(bary_sum(moire_patch[ia]) <= bary_sum(moire_patch[ib])) drop[ib] = 1;
            else drop[ia] = 1;
         };
         const double shifts[4][2] = {{0,0},{tx_x,tx_y},{ty_x,ty_y},{tx_x+ty_x,tx_y+ty_y}};
         for(size_t ia = 0; ia < moire_patch.size(); ++ia){
            if(drop[ia]) continue;
            for(int s = 0; s < 4; ++s){
               const double tx = moire_patch[ia].x + shifts[s][0];
               const double ty = moire_patch[ia].y + shifts[s][1];
               const int ix0 = (int)std::floor(tx * inv);
               const int iy0 = (int)std::floor(ty * inv);
               for(int dx = -1; dx <= 1; ++dx){
                  for(int dy = -1; dy <= 1; ++dy){
                     Key key{ix0+dx, iy0+dy, (int)moire_patch[ia].S, moire_patch[ia].l_id, moire_patch[ia].h_id};
                     auto it = bins.find(key);
                     if(it == bins.end()) continue;
                     for(size_t ib : it->second){
                        if(ib <= ia || drop[ib]) continue;
                        const double ddx = moire_patch[ib].x - tx;
                        const double ddy = moire_patch[ib].y - ty;
                        if(ddx*ddx + ddy*ddy < dup_tol*dup_tol) mark_pair(ia, ib);
                     }
                  }
               }
            }
         }
         size_t n_drop = 0;
         for(char d : drop) if(d) ++n_drop;
         if(n_drop > 0){
            std::vector<spin> deduped;
            deduped.reserve(moire_patch.size() - n_drop);
            uint64_t nid = 0;
            for(size_t k = 0; k < moire_patch.size(); ++k){
               if(drop[k]) continue;
               spin at = moire_patch[k];
               at.original_id = nid++;
               deduped.push_back(at);
            }
            std::cout << "build_spin_moire_cell: removed " << n_drop
                      << " period-duplicate atoms from patch" << std::endl;
            moire_patch.swap(deduped);
         }
      }

      // UCF/tile frame: AA → (0,0). Do NOT recompute unit bins or unit_cell_shifts.
      for(auto &at : moire_patch){
         at.x -= aa_x;
         at.y -= aa_y;
         at.Gx = 0;
         at.Gy = 0;
         // unit_x_lr / unit_y_lr / unit_x / unit_y unchanged
      }
      x0 = 0.0;
      y0 = 0.0;
      // i0,j0 stay as period-finder microcell origin (for bookkeeping only).

      // Drop any patch atoms outside the microcell grid; compact Vampire ids.
      {
         std::vector<spin> kept;
         kept.reserve(moire_patch.size());
         uint64_t nid = 0;
         for(auto &at : moire_patch){
            if(at.unit_x_lr < 0 || at.unit_x_lr > microcell_Nx ||
               at.unit_y_lr < 0 || at.unit_y_lr > microcell_Ny) continue;
            at.original_id = nid++;
            kept.push_back(at);
         }
         if(kept.empty()){
            std::cerr << "build_spin_moire_cell: candidate " << ci << " primary bins OOB" << std::endl;
            continue;
         }
         moire_patch.swap(kept);
      }

      const double span = std::max(system_size_x, system_size_y);
      int n_tiles_x = static_cast<int>(std::ceil(span / std::max(Lx, 1.0))) + 2;
      int n_tiles_y = static_cast<int>(std::ceil(span / std::max(Ly, 1.0))) + 2;

      std::vector<spin> tiled_lattice;
      tiled_lattice.reserve(moire_patch.size() * (size_t)(n_tiles_x+1) * (size_t)(n_tiles_y+1));
      uint64_t tiled_id = 0;
      for(int gi = -n_tiles_x; gi <= n_tiles_x; ++gi){
         for(int gj = -n_tiles_y; gj <= n_tiles_y; ++gj){
            for(size_t atom_i = 0; atom_i < moire_patch.size(); ++atom_i){
               spin new_atom(moire_patch[atom_i]);
               new_atom.x = moire_patch[atom_i].x + gi*tx_x + gj*ty_x;
               new_atom.y = moire_patch[atom_i].y + gi*tx_y + gj*ty_y;
               new_atom.z = moire_patch[atom_i].z;
               new_atom.Gx = gi;
               new_atom.Gy = gj;
               if(inside_system(system_size_x, system_size_y,
                                new_atom.x + aa_x, new_atom.y + aa_y, 0.0)
                  || inside_system(system_size_x, system_size_y, new_atom.x, new_atom.y, 0.0)){
                  // Same site geometry/shift as open-lattice microcell at
                  // primary_bin + (gi·Px, gj·Py); only ID is the primary-cell id.
                  new_atom.original_id = moire_patch[atom_i].original_id;
                  new_atom.unit_x_lr = moire_patch[atom_i].unit_x_lr + gi * Px;
                  new_atom.unit_y_lr = moire_patch[atom_i].unit_y_lr + gj * Py;
                  new_atom.unit_x = new_atom.unit_x_lr;
                  new_atom.unit_y = new_atom.unit_y_lr;
                  if(new_atom.unit_x_lr < 0 || new_atom.unit_x_lr > microcell_Nx ||
                     new_atom.unit_y_lr < 0 || new_atom.unit_y_lr > microcell_Ny) continue;
                  new_atom.id = tiled_id++;
                  tiled_lattice.push_back(new_atom);
               }
            }
         }
      }
      if(tiled_lattice.empty()){
         std::cerr << "build_spin_moire_cell: candidate " << ci << " empty after tile" << std::endl;
         continue;
      }

      std::cout << "build_spin_moire_cell: retile kept pre-tile bins/shifts; origin=(0,0) i0=" << i0
                << " j0=" << j0 << " tiled=" << tiled_lattice.size() << std::endl;

      // Working-lattice BL → (0,0) for exchange boxing only. Bins/shifts unchanged.
      // Primary UCF patch keeps AA at (0,0) and is NOT shifted here.
      {
         double wx0 = tiled_lattice[0].x, wy0 = tiled_lattice[0].y;
         for(const auto &at : tiled_lattice){
            if(at.x < wx0) wx0 = at.x;
            if(at.y < wy0) wy0 = at.y;
         }
         for(auto &at : tiled_lattice){
            at.x -= wx0;
            at.y -= wy0;
         }
         std::cout << "build_spin_moire_cell: working-lattice BL-=(" << wx0 << "," << wy0
                   << ") for boxing (UCF primary unshifted; bins unchanged)" << std::endl;
      }

      all_m_atoms = std::move(tiled_lattice);

      // Recompute unit_cell_shifts from the retiled lattice (diagnostic + exchange registry).
      {
         std::vector<std::vector<std::vector<int>>> unit_cell_shifts_new;
         compute_unit_cell_shifts_from_atoms(all_m_atoms, unit_cell_shifts_new);
         unit_cell_shifts.swap(unit_cell_shifts_new);

         char directory[256];
         if(getcwd(directory, sizeof(directory)) != NULL){
            std::ofstream sf(std::string(directory) + "/unit_cell_shifts_post_tile.txt");
            sf << "# unit_cell_shifts AFTER retile recompute\n";
            sf << "# Px Py i0 j0\n";
            sf << Px << " " << Py << " " << i0 << " " << j0 << "\n";
            sf << "# i j occupancy sx sy\n";
            for(int ii = 0; ii < (int)unit_cell_shifts.size(); ++ii){
               for(int jj = 0; jj < (int)unit_cell_shifts[ii].size(); ++jj){
                  sf << ii << " " << jj << " "
                     << unit_cell_shifts[ii][jj][0] << " "
                     << unit_cell_shifts[ii][jj][1] << " "
                     << unit_cell_shifts[ii][jj][2] << "\n";
               }
            }
            std::cout << "Wrote unit_cell_shifts_post_tile.txt after retile recompute" << std::endl;
         }
      }

      pbc_ok = true;
      moire_period_Px = Px;
      moire_period_Py = Py;
      moire_period_i0 = i0;
      moire_period_j0 = j0;
      moire_period_valid = true;

      std::cout << "build_spin_moire_cell: constructed periodic spin-moiré cell with Px=" << Px << " Py=" << Py
                << " atoms=" << moire_patch.size()
                << " (Tx/Ty tile; unit_cell_shifts recomputed post-tile)" << std::endl;
      {
         char directory[256];
         if(getcwd(directory, sizeof(directory)) != NULL){
            std::ofstream bout(std::string(directory) + "/moire_period_bake.txt");
            bout.precision(17);
            bout << "# x0 y0 Lx Ly tx_x tx_y ty_x ty_y Px Py i0 j0 natoms ax_eff ay_eff\n";
            bout << "# ax_eff=|Tx|/Px so Px*ax_eff=|Tx|; ay_eff=|Ty|/Py\n";
            bout << x0 << " " << y0 << " " << Lx << " " << Ly << " "
                 << tx_x << " " << tx_y << " " << ty_x << " " << ty_y << " "
                 << Px << " " << Py << " " << i0 << " " << j0 << " "
                 << moire_patch.size() << " "
                 << microcell_ax() << " " << microcell_ay() << "\n";
         }
      }
      break;
   }

   if(!pbc_ok){
      std::cerr << "build_spin_moire_cell: no candidate produced a tiled lattice." << std::endl;
      exit(1);
   }

   moire_area[0] = x0;
   moire_area[1] = y0;
   moire_area[2] = Lx;
   moire_area[3] = Ly;
   moire_area[4] = tx_x;
   moire_area[5] = tx_y;
   moire_area[6] = ty_x;
   moire_area[7] = ty_y;

   // Stash lattice + primary atoms for deferred Vampire UCF export (rotate after exchange).
   moire_origin_x = x0;
   moire_origin_y = y0;
   moire_tx_x = tx_x; moire_tx_y = tx_y;
   moire_ty_x = ty_x; moire_ty_y = ty_y;
   // R(-φ) with φ=atan2(tx_y,tx_x): maps Tx → (Lx,0), orthogonal Ty → (0,±Ly).
   // xp = c·x + s·y, yp = −s·x + c·y  ⇒  c=cosφ=tx_x/Lx, s=sinφ=tx_y/Lx.
   // (Previously s=−tx_y/Lx implemented R(+φ), skewing the UCF and causing fold duplicates.)
   {
      const double invL = 1.0 / Lx;
      moire_ucf_c = tx_x * invL;
      moire_ucf_s = tx_y * invL;
   }
   moire_primary_atoms = moire_patch;

   // Bake-only: no exchange, so write Vampire-orthogonal UCF now.
   // Full runs: write_vampire_ucf_rotated() after interactions so DMI rotates with positions.
   if(bake_only) write_vampire_ucf_rotated();
}

// Rotate code-frame G=0 positions so Vampire UCF is pure (Lx,0)/(0,Ly).
// Uses cached moire_ucf_c/s from build_spin_moire_cell. Maps / exchange stay in code frame.
void write_vampire_ucf_rotated(){
   if(moire_primary_atoms.empty()){
      std::cerr << "write_vampire_ucf_rotated: no primary atoms stored" << std::endl;
      return;
   }
   const double x0 = moire_origin_x, y0 = moire_origin_y;
   const double Lx = std::sqrt(moire_tx_x*moire_tx_x + moire_tx_y*moire_tx_y);
   const double Ly = std::sqrt(moire_ty_x*moire_ty_x + moire_ty_y*moire_ty_y);
   if(Lx < 1e-12 || Ly < 1e-12){
      std::cerr << "write_vampire_ucf_rotated: degenerate Tx/Ty" << std::endl;
      exit(1);
   }
   const double c = moire_ucf_c, s = moire_ucf_s;

   char directory[256];
   if(getcwd(directory, sizeof(directory)) == NULL){
      std::cerr << "Fatal getcwd error in write_vampire_ucf_rotated." << std::endl;
      return;
   }

   // Four Cr planes at z = 0, a0z, 2·a0z, 3·a0z (a0z = interlayer ≈ 5.778 Å).
   // UCF height = 4·interlayer = c0 ≈ 23.11 Å (AFM-style 4L), so fz = 0, 1/4, 1/2, 3/4.
   const double Lz = c0;

   std::ofstream outfile1(std::string(directory) + "/header.ucf");
   outfile1 << " #unit cell size " << std::endl;
   outfile1 << Lx << '\t' << Ly << '\t' << Lz << std::endl;
   outfile1 << " #unit cell vectors" << std::endl;
   outfile1 << " 1     0	   0" << std::endl;
   outfile1 << " 0     1	   0" << std::endl;
   outfile1 << " 0     0	   1" << std::endl;
   outfile1 << " #Atoms" << std::endl;
   outfile1 << moire_primary_atoms.size() << '\t'	<< 4 << std::endl;

   std::ofstream out(std::string(directory) + "/atom_positions.xyz");
   // After R(-φ), half-open parallelogram maps into [0,Lx)×[0,Ly) up to
   // barycentric edge_eps (~1e-4) → O(0.01–0.1) Å overshoot. Clamp that;
   // larger excursions still Fatal (would fmod-wrap into duplicates).
   const double fold_eps = 0.50; // Å (Tx/Ty shear edge overshoot)
   auto clamp01 = [](double v, double L, double eps)->int {
      // returns 0 ok, 1 clamped, 2 fatal outside
      if(v >= -eps && v < 0.0) { return 1; } // caller sets v=0
      if(v >= L && v < L + eps) { return 1; } // caller sets nextafter
      if(v < -eps || v >= L + eps) return 2;
      return 0;
   };
   int n_outside = 0;
   for(const auto &at : moire_primary_atoms){
      const double dx = at.x - x0;
      const double dy = at.y - y0;
      double xp = c*dx + s*dy;
      double yp = -s*dx + c*dy;
      const int cx = clamp01(xp, Lx, fold_eps);
      const int cy = clamp01(yp, Ly, fold_eps);
      if(cx == 2 || cy == 2){
         std::cerr << "Fatal: write_vampire_ucf_rotated: atom " << at.original_id
                   << " outside [0,L) after R(-φ): (" << xp << "," << yp
                   << ") L=(" << Lx << "," << Ly << ")" << std::endl;
         ++n_outside;
         continue;
      }
      if(cx == 1){
         if(xp < 0.0) xp = 0.0;
         else xp = std::nextafter(Lx, 0.0);
      }
      if(cy == 1){
         if(yp < 0.0) yp = 0.0;
         else yp = std::nextafter(Ly, 0.0);
      }
      const double fx = xp / Lx;
      const double fy = yp / Ly;
      double fz = at.z / Lz;
      if(fz < 0.0) fz = 0.0;
      if(fz >= 1.0) fz = std::nextafter(1.0, 0.0);

      out << at.original_id << "\t" << fx << "\t" << fy << "\t" << fz
          << "\t" << static_cast<int>(at.S)-1 << "\t" << at.l_id << "\t" << at.h_id << "\n";
   }
   if(n_outside > 0){
      std::cerr << "Fatal: write_vampire_ucf_rotated: " << n_outside
                << " atoms outside orthogonal cell after R(-φ) (would fold to duplicates)."
                << std::endl;
      exit(1);
   }
   std::cout << "write_vampire_ucf_rotated: wrote " << moire_primary_atoms.size()
             << " atoms, Lx=" << Lx << " Ly=" << Ly << " Lz=" << Lz
             << " c=" << c << " s=" << s << std::endl;
}
