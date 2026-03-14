#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <omp.h>
#include "positions.hpp"
#include "initialise.hpp"
#include "exchange.hpp"


double a0x = 6.93;
double a0y= 0.0;
double a1x = -3.465;
double a1y = 6.002;

double c0 = 26.16;
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

void print_header(){

   char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
         std::cerr << "Fatal getcwd error in datalog." << std::endl;
      }
   
      std::string header = "/header.ucf";
      std::ofstream outfile1(std::string(directory) + header);
   // std::ofstream outfile1 ("header.ucf");

   outfile1 << " #unit cell size " << std::endl;
   outfile1 << system_size_x << '\t' << system_size_y << '\t' << 26.0 << std::endl;
   outfile1 << " #unit cell vectors" << std::endl;
   outfile1 << " 1     0	   0" << std::endl;
   outfile1 << " 0     1	   0" << std::endl;
   outfile1 << " 0     0	   1" << std::endl;
   outfile1 << " #Atoms" << std::endl;
   outfile1 << total_atoms << '\t'	<< 4 << std::endl;

}

void write_atom_positions_ucf(){
   char directory [256];
   if(getcwd(directory, sizeof(directory)) == NULL){
      std::cerr << "Fatal getcwd error in write_atom_positions_ucf." << std::endl;
      return;
   }
   std::ofstream out(std::string(directory) + "/atom_positions.xyz");
   const double Lz = 26.0;
   for(const auto &at : all_m_atoms){
      out << total_atoms << "\t" << at.x/system_size_x << "\t" << at.y/system_size_y << "\t" << at.z/Lz
          << "\t" << static_cast<int>(at.S)-1 << "\t" << at.l_id << "\t" << at.h_id << "\n";
   }
}

bool inside_system(double sx, double sy, double x, double y, double offset){

   if ((x >=offset) && (x <= sx-offset) && (y >=offset) && (y <= sy-offset)) return true;
   else return false;
}

bool inside_system(double x, double y, double offset){

   if (x >=offset && x <= system_size_x-offset && y >=offset && y <= system_size_y -offset) return true;
   else return false;
}


void create_magnetic_atom_list(std::string filename){
   std::cout << "Generating lattice structure...." << std::flush;
   // double normalise_x = 100.0/(a0x*3.0);
   // double normalise_y = 100.0/(a0x*sqrt(3));
   char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
         std::cerr << "Fatal getcwd error in datalog." << std::endl;
      }
   

   // atom_positions.xyz is written after build_spin_moire_cell() via write_atom_positions_ucf()
   for (int i = -2*number_of_unit_cells_x; i < 3*number_of_unit_cells_x; i++) {
         for (int j = -2*number_of_unit_cells_y; j < 3*number_of_unit_cells_y; j++){
            // turn off replication in z to allow for explicit abba/abab stacking
            
               for (int atom_i = 0; atom_i < num_atoms; atom_i ++){

                  double x_j = atom[atom_i].x + i*a0x + j*a1x;
                  double y_j = atom[atom_i].y         + j*a1y;
                  double z_j = atom[atom_i].z; 
                  
                  if ( z_j > twist_loction){
                     // calculate rotated atom positions
                     double x_new = x_j*cos(twist_angle*0.5) - y_j*sin(twist_angle*0.5);
                     double y_new = y_j*cos(twist_angle*0.5) + x_j*sin(twist_angle*0.5);

                     // if atom is in system bounds, then generate it
                     if (inside_system(system_size_x, system_size_y, x_new, y_new, 0.0)){
                        
                        spin new_atom;
                        new_atom.x = x_new;
                        new_atom.y = y_new;
                        new_atom.z = z_j;
                        new_atom.id = total_atoms;
                        new_atom.l_id = atom[atom_i].l_id;
                        new_atom.h_id = atom[atom_i].h_id;
                       
                        new_atom.unit_y_lr = static_cast<int>(std::floor((y_new + 1e-9) / (microcell_scale_y * a1y)));
                        new_atom.unit_x_lr = static_cast<int>(std::floor((x_new + 1e-9) / (microcell_scale_x * a0x)));

                        int dy_cell = new_atom.unit_y_lr;
                        int dx_cell = new_atom.unit_x_lr;

                        // Set layer number
                        new_atom.unit_x = dx_cell;
                        new_atom.unit_y = dy_cell;      
                        
                        if(new_atom.unit_x_lr < 0 || new_atom.unit_x_lr > microcell_Nx || new_atom.unit_y_lr < 0 || new_atom.unit_y_lr > microcell_Ny) {
                           std::cerr << "microcell index out of range: unit_x_lr=" << new_atom.unit_x_lr << " unit_y_lr=" << new_atom.unit_y_lr
                                     << " (max " << microcell_Nx << "," << microcell_Ny << ") x=" << x_j << " y=" << y_j << std::endl;
                           std::exit(1);
                        }
                        if (z_j <= a0z*2){
                         
                           new_atom.S = 3;
                           new_atom.dx = 66;//dx;
                           new_atom.dy = 0;//dy;                          
                           global_config_energy[new_atom.unit_x_lr][new_atom.unit_y_lr][(new_atom.S-1)*16] += 1;                          
                        } else if (z_j <= a0z*3){
                           new_atom.S = 4;
                           new_atom.dx = 66;
                           new_atom.dy = 0;
                           global_config_energy[new_atom.unit_x_lr][new_atom.unit_y_lr][(new_atom.S-1)*16] += 1;
                        } else if (z_j == 600.0) {
                           new_atom.S = 5;
                        }else {
                           std::cerr << "Error! Atom " << total_atoms << " twist layer: " << z_j << " < " << twist_loction << std::endl;
                           
                           exit(1);
                        }    
                           total_atoms++;
                        // }
                        all_m_atoms.push_back(new_atom);
                        num_above_atoms++;
                     }
                  } else{
                     double x_new = x_j*cos(-twist_angle*0.5) - y_j*sin(-twist_angle*0.5);
                     double y_new = y_j*cos(-twist_angle*0.5) + x_j*sin(-twist_angle*0.5);

                     if(inside_system(system_size_x, system_size_y, x_new, y_new, 0.0)) { 
                     


                        spin new_atom;
                        new_atom.x = x_new;
                        new_atom.y = y_new;  
                        new_atom.z = z_j;
                        new_atom.id = total_atoms;
                        new_atom.l_id = atom[atom_i].l_id;
                        new_atom.h_id = atom[atom_i].h_id;

                        new_atom.unit_y_lr = static_cast<int>(std::floor((y_new + 1e-9) / (microcell_scale_y * a1y)));
                        new_atom.unit_x_lr = static_cast<int>(std::floor((x_new + 1e-9) / (microcell_scale_x * a0x)));

                        new_atom.unit_x = new_atom.unit_x_lr;
                        new_atom.unit_y = new_atom.unit_y_lr;

                        if(new_atom.unit_x_lr < 0 || new_atom.unit_x_lr > microcell_Nx || new_atom.unit_y_lr < 0 || new_atom.unit_y_lr > microcell_Ny) {
                           std::cerr << "microcell index out of range: unit_x_lr=" << new_atom.unit_x_lr << " unit_y_lr=" << new_atom.unit_y_lr
                                     << " (max " << microcell_Nx << "," << microcell_Ny << ") x=" << x_j << " y=" << y_j << std::endl;
                           std::exit(1);
                        }  

                        // Set layer number
                        if (z_j == 0.0){
                           new_atom.S = 1;
                           new_atom.dx = 66;
                           new_atom.dy = 0;
                           global_config_energy[new_atom.unit_x_lr][new_atom.unit_y_lr][(new_atom.S-1)*16 + 0] += 1;
                        
                        } else if (z_j <= a0z){
                           new_atom.S = 2;
                           new_atom.dx = 66; // need a dx,dy to take into account the actual stacking!
                           new_atom.dy = 0;
                           global_config_energy[new_atom.unit_x_lr][new_atom.unit_y_lr][(new_atom.S-1)*16 + 0] += 1;
                           
                        } else {
                              std::cerr << "Error! Atom " << total_atoms << " twist layer: " << z_j << " > " << twist_loction << std::endl;
                              exit(1);
                        }
                           total_atoms++;
                        
                        all_m_atoms.push_back(new_atom);       
                        num_below_atoms++;
                     }
                  } 
               
               }  
         } // j-loop
   } // i-loop

   std::cout << total_atoms << " atoms; [complete]" << std::endl;
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

                  // Layer assignment based on z, same as original
                  if (z_j <= a0z*2){
                     new_atom.S = 3;
                     new_atom.dx = 66;
                     new_atom.dy = 0;
                  } else if (z_j <= a0z*3){
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

                  // Layer assignment based on z, same as original
                  if (z_j == 0.0){
                     new_atom.S = 1;
                     new_atom.dx = 66;
                     new_atom.dy = 0;
                  } else if (z_j <= a0z){
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
      int dx_cell = static_cast<int>(std::floor((at.x + 1e-9) / (microcell_scale_x * a0x)));
      int dy_cell = static_cast<int>(std::floor((at.y + 1e-9) / (microcell_scale_y * a1y)));
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

void create_magnetic_atom_list_moire_unit(std::string filename, \
                  double Moire_aix, double Moire_aiy, double Moire_ajx, double Moire_ajy, \
                  double Moire_abs_x, double Moire_abs_y, int Moire_atom_size){
   std::cout << "Generating lattice structure from Moire unit cell file...." << std::flush;

   std::ofstream shift_file;
   shift_file.open("shifted_constants_Moire_ucf.txt");

   std::ofstream new_moire_output(filename);
   if(!new_moire_output.is_open()) {std::cout << "New moire did not open" << std::endl; exit(1);}

   int number_of_Moire_unit_cells_x = 1;
   int number_of_Moire_unit_cells_y = 1;
   std::cout << "reconstructing lattice using <" << number_of_Moire_unit_cells_x << ", " << number_of_Moire_unit_cells_y << "> unit cells " << std::endl;

   int new_lattice_atoms = 0;
   for (int i = -1*number_of_Moire_unit_cells_x; i < 2*number_of_Moire_unit_cells_x; i++){
         for (int j = -1*number_of_Moire_unit_cells_y; j < 2*number_of_Moire_unit_cells_y; j++){
               for (int atom_i = 0; atom_i < all_m_atoms_offset.size(); atom_i ++){

                  double x_j = all_m_atoms_offset[atom_i].x + i*Moire_aix + j*Moire_ajx;
                  double y_j = all_m_atoms_offset[atom_i].y + i*Moire_aiy + j*Moire_ajy;
                  double z_j = all_m_atoms_offset[atom_i].z;

                  if(inside_system(system_size_x, system_size_y, x_j, y_j, 0.0)){

                     spin new_atom(all_m_atoms_offset[atom_i]);
                     new_atom.x = x_j;
                     new_atom.y = y_j;
                     new_atom.Gx = i;
                     new_atom.Gy = j;

                     int dx_cell = static_cast<int>(std::floor((x_j + 1e-9) / (microcell_scale_x * a0x)));
                     int dy_cell = static_cast<int>(std::floor((y_j + 1e-9) / (microcell_scale_y * a1y)));
                     new_atom.unit_x_lr = dx_cell;
                     new_atom.unit_y_lr = dy_cell;
                     new_atom.unit_x = dx_cell;
                     new_atom.unit_y = dy_cell;

                     // new_moire_output << new_lattice_atoms << "\t" << x_j/(Moire_abs_x) << '\t' <<  y_j/(Moire_abs_y) <<  "\t" << z_j/system_size_z << "\t" << new_atom.S-1 << "\t" << new_atom.l_id << "\t" << new_atom.h_id << "\n"; 
                     new_lattice_atoms++;

                     new_moire_lattice.push_back(new_atom);              
                  }
               }  
         } // j-loop
   } // i-loop

   for(int i = 0; i < unit_cell_shifts.size(); i++){
      for (int j = 0; j < unit_cell_shifts[i].size(); j++) {
         int occupancy = std::max(1,unit_cell_shifts[i][j][0]);
         unit_cell_shifts[i][j][1] = round(unit_cell_shifts[i][j][1]/occupancy);
         unit_cell_shifts[i][j][2] = round(unit_cell_shifts[i][j][2]/occupancy);
         int i_shift = unit_cell_shifts[i][j][1];
         int j_shift = unit_cell_shifts[i][j][2];
         
         shift_file << i << ", " << j << ", " << occupancy << ", " << i_shift << ", " << j_shift << "\n";// << 

      }
   }
   shift_file.close();
   std::cout << new_lattice_atoms << " atoms; [complete]" << std::endl;
}


// Helper key for approximate atomic positions (for deduplication / PBC checks)
namespace {
   struct AtomKey {
      long ix;
      long iy;
      long iz;
      int S;
      int l_id;
      int h_id;
   };

   struct AtomKeyHash {
      std::size_t operator()(const AtomKey &k) const noexcept {
         // Simple mixed hash
         std::size_t h = std::hash<long>()(k.ix);
         h ^= std::hash<long>()(k.iy + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2));
         h ^= std::hash<long>()(k.iz + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2));
         h ^= std::hash<int>()(k.S + 0x9e3779b9 + int(h<<1));
         h ^= std::hash<int>()(k.l_id + 0x7f4a7c15 + int(h>>1));
         h ^= std::hash<int>()(k.h_id + 0x6a09e667 + int(h<<3));
         return h;
      }
   };

   struct AtomKeyEq {
      bool operator()(const AtomKey &a, const AtomKey &b) const noexcept {
         return a.ix == b.ix && a.iy == b.iy && a.iz == b.iz &&
                a.S == b.S && a.l_id == b.l_id && a.h_id == b.h_id;
      }
   };
}

// Find moiré periodicity from global_config_energy, unit_cell_shifts, all_m_atoms.
// Fills out and returns true on success; returns false (and logs to stderr) on failure.
bool find_moire_periodicity(MoirePeriodResult& out){
   if(global_config_energy.empty() || global_config_energy[0].empty()) {
      std::cerr << "find_moire_periodicity: global_config_energy not initialised" << std::endl;
      return false;
   }

   const int Nx = static_cast<int>(global_config_energy.size());
   const int Ny = static_cast<int>(global_config_energy[0].size());
   const int shifts_Nx = static_cast<int>(unit_cell_shifts.size());
   const int shifts_Ny = (shifts_Nx > 0 && !unit_cell_shifts[0].empty()) ? static_cast<int>(unit_cell_shifts[0].size()) : 0;

   const int shift_tol = 1;  // tolerance for binned shift values (0–99 scale): descriptor match and AA selection (66, 0)

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
      bool in1 = (i1 >= 0 && i1 < shifts_Nx-1 && j1 >= 0 && j1 < shifts_Ny-1);
      bool in2 = (i2 >= 0 && i2 < shifts_Nx-1 && j2 >= 0 && j2 < shifts_Ny-1);
      if(in1 != in2) return false; // one in bounds, one out -> not equal
      if(in1 && in2) {
         if(unit_cell_shifts[i1][j1][0] < 1 || unit_cell_shifts[i2][j2][0] < 1) return false;
         const int d1 = std::abs(unit_cell_shifts[i1][j1][1] - unit_cell_shifts[i2][j2][1]);
         const int d2 = std::abs(unit_cell_shifts[i1][j1][2] - unit_cell_shifts[i2][j2][2]);
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

   const double cell_ax = static_cast<double>(microcell_scale_x) * a0x;
   const double cell_ay = static_cast<double>(microcell_scale_y) * a1y;

   // Minimum dx, dy from analytic moiré period L = a/(2*sin(θ/2)) so we target the true AA–AA distance (~360 Å at 1.1°, a=6.9 Å), not a shorter sub-repeat (~70 Å).
   const double L_moire_angstrom = a0x / (2.0 * std::sin(twist_angle / 2.0));
   const int expected_cells_x = static_cast<int>(std::floor(L_moire_angstrom / cell_ax + 0.5));
   const int expected_cells_y = static_cast<int>(std::floor(L_moire_angstrom / cell_ay + 0.5));
   const int min_dx_period = std::max(2, static_cast<int>(expected_cells_x * 0.25));
   const int min_dy_period = std::max(2, static_cast<int>(expected_cells_y * 0.25));
   std::cout << "build_spin_moire_cell: analytic moiré period L=" << L_moire_angstrom << " Å; min_dx=" << min_dx_period << " min_dy=" << min_dy_period << " (expected " << expected_cells_x << "x" << expected_cells_y << " cells)." << std::endl;

   // --- Step 1: AA cells from unit_cell_shifts (shift near 66, 0 = AA configuration), within shift_tol ---
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
   const double atom_match_tol = 1e-1;

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
   #pragma omp parallel num_threads(16)
   {
      std::vector<CandidateCell> my_candidates;

      #pragma omp for schedule(dynamic) reduction(+:cell_count) reduction(+:atom_count)
      for(size_t idx = 0; idx < aa_cells.size(); ++idx){
         const int i0_cand = aa_cells[idx].first, j0_cand = aa_cells[idx].second;

         // Scan over possible dx,dy periods; for each pair that passes edge checks,
         // add a candidate but keep searching for larger periods as well.
         for(int dx = min_dx_period; i0_cand + dx <= i_max; ++dx){
            if(!descriptors_equal(i0_cand, j0_cand, i0_cand + dx, j0_cand)) continue;

            for(int dy = min_dy_period; j0_cand + dy <= j_max; ++dy){
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
      std::cerr << "find_moire_periodicity: no AA candidates passed descriptor edge checks (Step 2). Cells passed to diagonal: " << cell_count << std::endl;
      return false;
   }

   // Step 3: brute-force search, over all descriptor-accepted candidates, for atomic quartets defining orthogonal
   // translation vectors. Each quartet uses atoms from the four cells:
   // (i0,j0), (i0+Px,j0), (i0,j0+Py), (i0+Px,j0+Py).
   struct QuartetResult {
      double area;
      double mse_ortho_shear;  // MSE of orthogonality, horiz_shear, vert_shear (ideal 0 for all)
      double Lx, Ly;
      double tx_x, tx_y;  // x-period translation vector (atom a -> atom b)
      double ty_x, ty_y;  // y-period translation vector (atom a -> atom c)
      double x0, y0;
      int i0, j0, Px, Py;
   };

   bool have_best_quartet = false;
   QuartetResult best_q{};
   std::vector<QuartetResult> all_quartets;

   const double orthogonality_tol = 0.00015; // max |cos(theta)| for orthogonality
   const double horiz_shear_tol   = 0.005*cell_ay; // allow small vertical drift in x-translation
   const double vert_shear_tol   = 0.005*cell_ax; // allow small horizontal drift in y-translation

   #pragma omp parallel num_threads(16)
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
               if(vx_x <= 0.0) continue;
               if(std::fabs(vx_y) > horiz_shear_tol) continue;

               const double Lx_candidate = std::sqrt(vx_x*vx_x + vx_y*vx_y);
               if(Lx_candidate <= 0.0) continue;

               for(const spin* c_top : top_cells){
                  if(c_top->S != a->S || c_top->l_id != a->l_id) continue;
                  const double vy_x = c_top->x - a->x;
                  const double vy_y = c_top->y - a->y;
                  if(vy_y <= 0.0) continue;
                  if(std::fabs(vy_x) > vert_shear_tol) continue;

                  const double Ly_candidate = std::sqrt(vy_x*vy_x + vy_y*vy_y);
                  if(Ly_candidate <= 0.0) continue;

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

                  const double area = std::fabs(vx_x*vy_y - vx_y*vy_x);
                  if(area <= 0.0) continue;

                  // MSE of orthogonality, horiz_shear, vert_shear (ideal 0 for all); normalize by tolerances for comparable scales
                  const double ortho_val = cos_theta / (orthogonality_tol > 1e-12 ? orthogonality_tol : 1.0);
                  const double horiz_val = std::fabs(vx_y) / (horiz_shear_tol > 1e-12 ? horiz_shear_tol : 1.0);
                  const double vert_val  = std::fabs(vy_x) / (vert_shear_tol > 1e-12 ? vert_shear_tol : 1.0);
                  const double mse = (ortho_val*ortho_val + horiz_val*horiz_val + vert_val*vert_val) / 3.0;

                  QuartetResult q{};
                  q.area = area;
                  q.mse_ortho_shear = mse;
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
             [](const QuartetResult &a, const QuartetResult &b){
                return a.mse_ortho_shear < b.mse_ortho_shear;
             });

   // choose best_q as lowest-MSE quartet from the reduced list
   best_q = all_quartets.front();
   have_best_quartet = true;

   // dump full reduced list for diagnostics
   {
      char directory[256];
      if(getcwd(directory, sizeof(directory)) != NULL){
         std::ofstream qout(std::string(directory) + "/moire_quartets.txt");
         qout << "# i0 j0 Px Py Lx Ly tx_x tx_y ty_x ty_y area mse\n";
         for(const auto &q : all_quartets){
            qout << q.i0 << " " << q.j0 << " "
                 << q.Px << " " << q.Py << " "
                 << q.Lx << " " << q.Ly << " "
                 << q.tx_x << " " << q.tx_y << " "
                 << q.ty_x << " " << q.ty_y << " "
                 << q.area << " " << q.mse_ortho_shear << "\n";
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
   return true;
}

// Build a spin-moiré primary cell from the large patch and perform an atomic-level PBC check.
// If use_given is non-null and has Lx > 0, Ly > 0, period detection is skipped and the given lattice vectors/cutoffs are used.
void build_spin_moire_cell(double moire_area[4], const MoirePeriodResult* use_given){
   if(global_config_energy.empty() || global_config_energy[0].empty()) {
      std::cerr << "build_spin_moire_cell: global_config_energy not initialised" << std::endl;
      return;
   }

   double x0, y0, Lx, Ly, tx_x, tx_y, ty_x, ty_y;
   int Px = 0, Py = 0, i0 = 0, j0 = 0;
   bool from_detection = false;

   if(moire_area[2] > 0 && moire_area[3] > 0){
      x0 = moire_area[0];
      y0 = moire_area[1];
      Lx = moire_area[2];
      Ly = moire_area[3];
      tx_x = Lx;
      tx_y = 0;
      ty_x = 0;
      ty_y = Ly;
      Px = static_cast<int>(std::floor((Lx+1e-9)/(microcell_scale_x * a0x)));
      Py = static_cast<int>(std::floor((Ly+1e-9)/(microcell_scale_y * a1y)));
      i0 = static_cast<int>(std::floor((x0 + 1e-9) / (microcell_scale_x * a0x)));
      j0 = static_cast<int>(std::floor((y0 + 1e-9) / (microcell_scale_y * a1y)));
      from_detection = true;

      std::cout << "Using read-in moire area: x0=" << x0 << " y0=" << y0 << " Lx=" << Lx << " Ly=" << Ly << std::endl;
      std::cout << "Px=" << Px << " Py=" << Py << " i0=" << i0 << " j0=" << j0 << std::endl;
      std::cout << "tx_x=" << tx_x << " tx_y=" << tx_y << " ty_x=" << ty_x << " ty_y=" << ty_y << std::endl;
      std::cout << "from_detection=" << from_detection << std::endl;
   } else {
      MoirePeriodResult res;
      if(!find_moire_periodicity(res)){
         exit(1);
      }
      x0 = res.x0;
      y0 = res.y0;
      Lx = res.Lx;
      Ly = res.Ly;
      tx_x = res.tx_x;
      tx_y = res.tx_y;
      ty_x = res.ty_x;
      ty_y = res.ty_y;
      Px = res.Px;
      Py = res.Py;
      i0 = res.i0;
      j0 = res.j0;
      from_detection = true;
   }

   const int shift_tol = 1;

   std::cout << "moire lattice vectors: tx_x=" << tx_x << " tx_y=" << tx_y << " ty_x=" << ty_x << " ty_y=" << ty_y << std::endl;
   std::vector<spin> moire_patch;
   uint64_t new_id = 0;
   //use exact atomic positions to determine the moire patch
   for(size_t k = 0; k < all_m_atoms.size(); ++k){
      spin at = all_m_atoms[k];
      if(at.x >= x0 - 0.0016869 && at.x < x0 + Lx - 0.0016869 &&
         at.y >= y0 - 0.0032932 && at.y < y0 + Ly - 0.00329322){
         at.original_id = new_id;
         new_id++;
         moire_patch.push_back(at);
      }
   }

   std::cout << "moire patch size: " << moire_patch.size() << std::endl;

   int n_tiles_x = static_cast<int>(std::ceil(system_size_x / Lx)) + 2;
   int n_tiles_y = static_cast<int>(std::ceil(system_size_y / Ly)) + 2;

   std::vector<spin> tiled_lattice;
   tiled_lattice.reserve(moire_patch.size() * (n_tiles_x+1) * (n_tiles_y+1));
   uint64_t tiled_id = 0;
   for(int gi = -n_tiles_x; gi <= n_tiles_x; ++gi){
      for(int gj = -n_tiles_y; gj <= n_tiles_y; ++gj){
         for(size_t atom_i = 0; atom_i < moire_patch.size(); ++atom_i){
            spin new_atom(moire_patch[atom_i]);
            new_atom.x = moire_patch[atom_i].x + gi*tx_x + 0.0*gj*ty_x;
            new_atom.y = moire_patch[atom_i].y + 0.0*gi*tx_y + gj*ty_y;
            new_atom.z = moire_patch[atom_i].z;
            new_atom.Gx = gi;
            new_atom.Gy = gj;
            if(inside_system(system_size_x, system_size_y, new_atom.x, new_atom.y, 0.0)){
               new_atom.original_id = moire_patch[atom_i].original_id;
               new_atom.id = tiled_id++;
               new_atom.unit_x_lr = static_cast<int>(std::floor((new_atom.x + 1e-9) / (microcell_scale_x * a0x)));
               new_atom.unit_y_lr = static_cast<int>(std::floor((new_atom.y + 1e-9) / (microcell_scale_y * a1y)));
               new_atom.unit_x = new_atom.unit_x_lr;
               new_atom.unit_y = new_atom.unit_y_lr;
               tiled_lattice.push_back(new_atom);
            }
         }
      }
   }

   // std::vector<std::vector<std::vector<int>>> unit_cell_shifts_old = unit_cell_shifts;
   all_m_atoms = std::move(tiled_lattice);

   std::vector<std::vector<std::vector<int>>> unit_cell_shifts_new;
   compute_unit_cell_shifts_from_atoms(all_m_atoms, unit_cell_shifts_new);

   if(from_detection){
      if(unit_cell_shifts.size() != unit_cell_shifts_new.size()){
         std::cerr << "build_spin_moire_cell: unit_cell_shifts mismatch after tiling: size " << unit_cell_shifts.size() << " vs " << unit_cell_shifts_new.size() << std::endl;
         exit(1);
      }
      for(size_t i = i0; i < i0+Px; i++){
         if(unit_cell_shifts[i].size() != unit_cell_shifts_new[i].size()){
            std::cerr << "build_spin_moire_cell: unit_cell_shifts mismatch at i=" << i << std::endl;
            exit(1);
         }
         for(size_t j = j0; j < j0+Py; j++){
            for(int c = 0; c < 3; ++c){
               if( std::abs(unit_cell_shifts[i][j][c] - unit_cell_shifts_new[i][j][c]) > shift_tol*1 ) {
                  std::cerr << "build_spin_moire_cell: unit_cell_shifts mismatch at (" << i << "," << j << ",[" << c << "]): old=" << unit_cell_shifts[i][j][c] << " new=" << unit_cell_shifts_new[i][j][c] << std::endl;
                  exit(1);
               }
            }
         }
      }
   }
   std::cout << "build_spin_moire_cell: constructed periodic spin-moiré cell with Px=" << Px << " Py=" << Py << " atoms=" << moire_patch.size() << std::endl;
   moire_area[0] = x0;
   moire_area[1] = y0;
   moire_area[2] = Lx;
   moire_area[3] = Ly;

   
   char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
         std::cerr << "Fatal getcwd error in datalog." << std::endl;
      }
   
      std::string header = "/header.ucf";
      std::ofstream outfile1(std::string(directory) + header);
   // std::ofstream outfile1 ("header.ucf");

   outfile1 << " #unit cell size " << std::endl;
   outfile1 << Lx << '\t' << Ly << '\t' << 26.0 << std::endl;
   outfile1 << " #unit cell vectors" << std::endl;
   outfile1 << " 1     0	   0" << std::endl;
   outfile1 << " 0     1	   0" << std::endl;
   outfile1 << " 0     0	   1" << std::endl;
   outfile1 << " #Atoms" << std::endl;
   outfile1 << moire_patch.size() << '\t'	<< 4 << std::endl;


   std::ofstream out(std::string(directory) + "/atom_positions.xyz");
   const double Lz = 26.0;
   double x = 0.0;
   double y = 0.0;
   for(const auto &at : moire_patch){
      x = at.x- x0- 0.0016869;
      y = at.y- y0- 0.0032932;
      if(x < 0 ) x = 0.0;
      if(y < 0 ) y = 0.0;

      if(x/Lx > 1 || y/Ly > 1) {
         std::cout << "atom " << at.original_id << " is outside the moire patch, x=" << x << " y=" << y << std::endl;
         exit(1);
      }
      
      out << at.original_id << "\t" << x/Lx << "\t" << y/Ly << "\t" << at.z/Lz
          << "\t" << static_cast<int>(at.S)-1 << "\t" << at.l_id << "\t" << at.h_id << "\n";
   }
}
