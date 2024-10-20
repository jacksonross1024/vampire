#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "positions.hpp"
#include "initialise.hpp"
#include "exchange.hpp"


double a0x = 6.93;
double a0y= 0.0;
double a1x = -3.465;
double a1y = 6.002;

double c0 = 26.16/2;
double a0z = c0/2.0;

int num_atoms = 4;
int num_nm_atoms = 24;

int number_of_unit_cells_x;
int number_of_unit_cells_y;

int num_above_atoms =0;
int num_below_atoms =0;


int total_atoms = 1;
int total_nm_atoms = 0;

void print_header(){

   std::ofstream outfile1 ("header.ucf");

   outfile1 << " #unit cell size " << std::endl;
   outfile1 << system_size_x << '\t' << system_size_y << '\t' << system_size_z << std::endl;
   outfile1 << " #unit cell vectors" << std::endl;
   outfile1 << " 1     0	   0" << std::endl;
   outfile1 << " 0     1	   0" << std::endl;
   outfile1 << " 0     0	   1" << std::endl;
   outfile1 << " #Atoms" << std::endl;
   outfile1 << total_atoms << '\t'	<< 4 << std::endl;

}

bool inside_system(double sx, double sy, double x, double y, double offset){

   if (x >=offset && x <= sx-offset && y >=offset && y <= sy -offset) return true;
   else return false;
}

bool inside_system(double x, double y, double offset){

   if (x >=offset && x <= system_size_x-offset && y >=offset && y <= system_size_y -offset) return true;
   else return false;
}


int calc_dxy(const double x_new, const double x, const double normalise){

   double change = fabs(x_new - x);
   int val_x = sqrt(change*change)*normalise;
   const int hundreds = int(val_x/100);

   //while (val_x > 99){
   //   val_x = val_x - 100;
   //}
   return val_x - hundreds*100;
}

void create_magnetic_atom_list(std::string filename){
   std::cout << "Generating lattice structure...." << std::flush;
   // double normalise_x = 100.0/(a0x*3.0);
   // double normalise_y = 100.0/(a0x*sqrt(3));
   std::ofstream shift_file;
   shift_file.open("shifted_constants.txt");

   // resize_arrays(unit_cell_shifts, number_of_unit_cells_x, number_of_unit_cells_y);
   // int total_atoms_kept = 1;
   for (int i = -1*number_of_unit_cells_x; i < 2*number_of_unit_cells_x; i++){
         for (int j = -1*number_of_unit_cells_y; j < 2*number_of_unit_cells_y; j++){
            // turn off replication in z to allow for explicit abba/abab stacking
            //for (int k = 0; k < number_of_unit_cells_z; k++){
               for (int atom_i = 0; atom_i < num_atoms; atom_i ++){

                  double x_j = atom[atom_i].x + i*a0x + j*a1x;
                  double y_j = atom[atom_i].y         + j*a1y;

                  double z_j = atom[atom_i].z; // + k*c0;
                  // std::cout << x_j << ", " << y_j << ", " << z_j << std::endl;
                  //std::cout << z_j << '\t' << twist_loction <<std::endl;
                  // only twist the top two layers z_j > twist
                  if ( z_j > twist_loction){
                     // calculate rotated atom positions
                     double x_new = x_j*cos(twist_angle*0.5) - y_j*sin(twist_angle*0.5);
                     double y_new = y_j*cos(twist_angle*0.5) + x_j*sin(twist_angle*0.5);
                     if( !inside_system(system_size_x, system_size_y, x_new, y_new, 0.0) && (inside_system(system_size_x, system_size_y, x_new, y_new, -0.000)) ){
                        if(x_new < 0.0) x_new = 0.0;
                        else if(x_new > system_size_x) x_new = system_size_x-0.0001;

                        if(y_new < 0.0) y_new = 0.0;
                        else if (y_new > system_size_y) y_new = system_size_y-0.0001;
                     }
                     // if atom is in system bounds, then generate it
                     if (inside_system(system_size_x, system_size_y, x_new, y_new, 0.0)){
                        
                        spin new_atom;
                        new_atom.x = x_new;
                        new_atom.y = y_new;
                        new_atom.z = z_j;
                        new_atom.id = total_atoms;
                        new_atom.l_id = atom[atom_i].l_id;
                        new_atom.h_id = atom[atom_i].h_id;
                        // new_atom.S = 1;
                        
                        // double changex = std::abs(x_new - x_j);
                        // double changey = std::abs(y_new - y_j);
                        
                        int dy_cell = floor((y_new +0.0000001)/ a1y);
                        // changex += dy_cell*std::abs(a1x);
                        int dx_cell = floor((x_new +0.0000001)/ a0x);
                        // double unit_x = dx_cell*a0x + dy_cell*a1x + atom[atom_i].x;
                        // double unit_y = dy_cell*a1y + atom[atom_i].y;
                        // int changex = round(-100.0*remainder(x_new - x_j,a0x)/a0x)+100;
                        // int changey = round(-100.0*remainder(y_new - y_j,a1y)/a0x)+100;
                        double x_eff = x_j*cos(twist_angle) - y_j*sin(twist_angle);
                        double y_eff = y_j*cos(twist_angle) + x_j*sin(twist_angle);
                        int changey = int(round(10*(fmod(std::abs(y_eff-y_j) , a1y)/a1y)));
                        int changex = int(round(9*(fmod(std::abs(x_eff-x_j +changey*a1y/11.0) , a0x)/a0x)));
                        
                        if(changex > 9 || changex < 0 || changey > 10 || changey < 0) {
                           std::cerr << "shift problem: (" << x_new << ", " << x_j << ") in cell: [" << dx_cell << ", " << dy_cell << "] indexing " << changex << ", " << changey  << std::endl;
                            exit(1);
                        }
                       
                        // Set layer number
                        new_atom.unit_x = dx_cell;
                        new_atom.unit_y = dy_cell;

                        if (z_j <= a0z*2){
                           new_atom.S = 3;
                           new_atom.dx = changex;
                           new_atom.dy = changey;
                            unit_cell_shifts.at(dx_cell).at(dy_cell)[0] += 1;
                           unit_cell_shifts[dx_cell][dy_cell][1] += changex;
                           unit_cell_shifts[dx_cell][dy_cell][2] += changey;
                           row3.push_back(new_atom);
                        } else if (z_j <= a0z*3){
                           new_atom.S = 3;
                           // new_atom.dx = changex;
                           // new_atom.dy = changey;
                           // row4.push_back(new_atom);
                        } else {
                           std::cerr << "Error! Atom " << total_atoms << " twist layer: " << z_j << " < " << twist_loction << std::endl;
                           
                           exit(1);
                        }    
                           // outfile2 << total_atoms << "\t" << x_new/(system_size_x) << '\t' <<  y_new/(system_size_y) <<  "\t" << z_j/system_size_z << "\t" << new_atom.S-1 << "\t" << new_atom.l_id << "\t" << new_atom.h_id << "\n"; 
                           total_atoms++;
                        // }
                        all_m_atoms.push_back(new_atom);
                        num_above_atoms++;
                     }
                  } else if (inside_system(system_size_x, system_size_y, x_j, y_j, 0.0)){  // not twisted layer
                     double x_new = x_j*cos(-twist_angle*0.5) - y_j*sin(-twist_angle*0.5);
                     double y_new = y_j*cos(-twist_angle*0.5) + x_j*sin(-twist_angle*0.5);
                     if( !inside_system(system_size_x, system_size_y, x_j, y_j, 0.0) && (inside_system(system_size_x, system_size_y, x_j, y_j, -0.00)) ){
                        if(x_new < 0.0) x_new = 0.0;
                        else if(x_new > system_size_x) x_new = system_size_x-0.0001;

                        if(y_new < 0.0) y_new = 0.0;
                        else if (y_new > system_size_y) y_new = system_size_y-0.0001;
                     }

                     spin new_atom;
                     new_atom.x = x_j;
                     new_atom.y = y_j;
                     new_atom.z = z_j;
                     new_atom.id = total_atoms;
                     new_atom.l_id = atom[atom_i].l_id;
                     new_atom.h_id = atom[atom_i].h_id;
                     // new_atom.S = 0;
                     new_atom.unit_y = int(floor((y_j +0.00001)/ a1y));
                        // changex += dy_cell*std::abs(a1x);
                     new_atom.unit_x = int(floor((x_j +0.00001)/ a0x));
                     if(new_atom.unit_x  > number_of_unit_cells_x || new_atom.unit_y > number_of_unit_cells_y) {
                        std::cerr << new_atom.unit_x  << ", " << new_atom.unit_y  << ", " << x_j << ", " << y_j << ", " << \
                         int(floor(y_j / a1y)) << ", " <<  int(floor(x_j / a0x)) << std::endl;
                         std::exit(1);
                     }
                     // Set layer number
                     if (z_j == 0.0){
                        new_atom.S = 2;
                        // new_atom.dx = 0; // need a dx,dy to take into account the actual stacking!
                        // new_atom.dy = 0;
                        row1.push_back(new_atom);
                        //std::cout << total_atoms << "\t" << new_atom.S << "\t" << new_atom.dx << "\t" << new_atom.dy << "\t" << Jint[new_atom.dx][new_atom.dy] << std::endl;
                     } else if (z_j <= a0z){
                        new_atom.S = 2;
                        // new_atom.dx = 0; // need a dx,dy to take into account the actual stacking!
                        // new_atom.dy = 0;
                        // row2.push_back(new_atom);
                        //std::cout << total_atoms << "\t" << new_atom.S << "\t" << new_atom.dx << "\t" << new_atom.dy << "\t" << Jint[new_atom.dx][new_atom.dy] << std::endl;
                     } else {
                           std::cerr << "Error! Atom " << total_atoms << " twist layer: " << z_j << " > " << twist_loction << std::endl;
                           exit(1);
                     }
                        // outfile2 << total_atoms << "\t" << x_j/(system_size_x) << '\t' <<  y_j/(system_size_y) <<  "\t" << z_j/system_size_z << "\t" << new_atom.S-1 << "\t" << new_atom.l_id << "\t" << new_atom.h_id << "\n"; 
                        total_atoms++;
                     // }
                     all_m_atoms.push_back(new_atom);       
                     num_below_atoms++;
                  } 
               }  
         } // j-loop
   } // i-loop

  
   // if(row1.size() != row4.size() || row2.size() != row3.size() || total_atoms != all_m_atoms.size()) {
   //    std::cout << row1.size() << "\t" << row2.size() << "\t" << row3.size() << "\t" << row4.size() << std::endl;
   //    // exit(1);
   // }

   for(int i = 0; i < unit_cell_shifts.size(); i++){
      for (int j = 0; j < unit_cell_shifts[i].size(); j++) {
         double occupancy = std::max(1,unit_cell_shifts[i][j][0]);
         unit_cell_shifts[i][j][1] = round(unit_cell_shifts[i][j][1]/occupancy);
         unit_cell_shifts[i][j][2] = round(unit_cell_shifts[i][j][2]/occupancy);
         int i_shift = unit_cell_shifts[i][j][1];
         int j_shift = unit_cell_shifts[i][j][2];
         //  std::cout << "problems " << unit_cell_shifts[i][j][2] << ", " << j_shift << ", " << occupancy << std::endl;
         shift_file << i << ", " << j << ", " << occupancy << ", " << i_shift << ", " << j_shift << "\n";// << 
                        // Einter_Cr1.at(i_shift).at(j_shift)[2]  << ", " <<\
                        // Einter_Cr1.at(i_shift).at(j_shift) << ", " << \
                        // Einter_Cr1.at(i_shift).at(j_shift) <<  ", " << \
                        // Dx_inter.at(i_shift).at(j_shift) << ", " << \
                        // Dy_inter.at(i_shift).at(j_shift) << ", " << \
                        // Dz_inter.at(i_shift).at(j_shift) << ", " << \
                        // Dx_intra.at(i_shift).at(j_shift) << ", " << \
                        // Dy_intra.at(i_shift).at(j_shift) << ", " << \
                        // Dz_intra.at(i_shift).at(j_shift) << "\n";
      }
   }
   shift_file.close();
   std::cout << total_atoms << " atoms; [complete]" << std::endl;
}

void create_magnetic_atom_list_moire_unit(std::string filename, \
                  double Moire_a0x, double Moire_a0y, double Moire_a1x, double Moire_a1y, \
                  double Moire_abs_x, double Moire_abs_y, int Moire_atom_size){
   std::cout << "Generating lattice structure from Moire unit cell file...." << std::flush;
   // double normalise_x = 100.0/(a0x*3.0);
   // double normalise_y = 100.0/(a0x*sqrt(3));
   std::ofstream shift_file;
   shift_file.open("shifted_constants_Moire_ucf.txt");

   std::ofstream new_moire_output(filename);
   if(!new_moire_output.is_open()) {std::cout << "New moire did not open" << std::endl; exit(1);}

   int number_of_Moire_unit_cells_x = ceil(system_size_x/Moire_abs_x);
   int number_of_Moire_unit_cells_y = ceil(system_size_x/Moire_abs_y);
   // resize_arrays(unit_cell_shifts, number_of_unit_cells_x, number_of_unit_cells_y);
   // int total_atoms_kept = 1;
   int new_lattice_atoms = 0;
   for (int i = -1*number_of_Moire_unit_cells_x; i < 2*number_of_Moire_unit_cells_x; i++){
         for (int j = -1*number_of_Moire_unit_cells_y; j < 2*number_of_Moire_unit_cells_y; j++){
            // turn off replication in z to allow for explicit abba/abab stacking
            //for (int k = 0; k < number_of_unit_cells_z; k++){
               for (int atom_i = 0; atom_i < all_m_atoms_offset.size(); atom_i ++){

                  double x_j = all_m_atoms_offset[atom_i].x + i*Moire_a0x + j*Moire_a0y;
                  double y_j = all_m_atoms_offset[atom_i].y + i*Moire_a1x + j*Moire_a1y;
                  double z_j = all_m_atoms_offset[atom_i].z;

                  spin new_atom(all_m_atoms_offset[atom_i]);
                  new_atom.x = x_j;
                  new_atom.y = y_j;
                  new_atom.Gx = i;
                  new_atom.Gy = j;
                  new_atom.id = new_lattice_atoms;

                  int dy_cell = floor((y_j +0.0000001)/ a1y);
                  int dx_cell = floor((x_j +0.0000001)/ a0x);
                  new_atom.unit_x = dx_cell;
                  new_atom.unit_y = dy_cell;

                  if(new_atom.z >= a0z) {
                        unit_cell_shifts.at(dx_cell).at(dy_cell)[0] += 1;
                        unit_cell_shifts[dx_cell][dy_cell][1] += new_atom.dx;
                        unit_cell_shifts[dx_cell][dy_cell][2] += new_atom.dy;
                        // row3.push_back(new_atom);
                  }  
   
                  new_moire_output << new_lattice_atoms << "\t" << x_j/(Moire_abs_x) << '\t' <<  y_j/(Moire_abs_y) <<  "\t" << z_j/system_size_z << "\t" << new_atom.S-1 << "\t" << new_atom.l_id << "\t" << new_atom.h_id << "\n"; 
                  new_lattice_atoms++;

                     new_moire_lattice.push_back(new_atom);              
                  
               }  
         } // j-loop
   } // i-loop

   for(int i = 0; i < unit_cell_shifts.size(); i++){
      for (int j = 0; j < unit_cell_shifts[i].size(); j++) {
         double occupancy = std::max(1,unit_cell_shifts[i][j][0]);
         unit_cell_shifts[i][j][1] = round(unit_cell_shifts[i][j][1]/occupancy);
         unit_cell_shifts[i][j][2] = round(unit_cell_shifts[i][j][2]/occupancy);
         int i_shift = unit_cell_shifts[i][j][1];
         int j_shift = unit_cell_shifts[i][j][2];
         //  std::cout << "problems " << unit_cell_shifts[i][j][2] << ", " << j_shift << ", " << occupancy << std::endl;
         shift_file << i << ", " << j << ", " << occupancy << ", " << i_shift << ", " << j_shift << "\n";// << 
                        // Einter_Cr1.at(i_shift).at(j_shift)[2]  << ", " <<\
                        // Einter_Cr1.at(i_shift).at(j_shift) << ", " << \
                        // Einter_Cr1.at(i_shift).at(j_shift) <<  ", " << \
                        // Dx_inter.at(i_shift).at(j_shift) << ", " << \
                        // Dy_inter.at(i_shift).at(j_shift) << ", " << \
                        // Dz_inter.at(i_shift).at(j_shift) << ", " << \
                        // Dx_intra.at(i_shift).at(j_shift) << ", " << \
                        // Dy_intra.at(i_shift).at(j_shift) << ", " << \
                        // Dz_intra.at(i_shift).at(j_shift) << "\n";
      }
   }
   shift_file.close();
   std::cout << total_atoms << " atoms; [complete]" << std::endl;
   
}


void create_nm_atom_list(){


   for (int i = -number_of_unit_cells_x; i < 2*number_of_unit_cells_x; i++){
         for (int j = -number_of_unit_cells_y; j < 2*number_of_unit_cells_y; j++){
            for (int k = 0; k < number_of_unit_cells_z; k++){
               for (int atom_i = 0; atom_i < num_nm_atoms; atom_i ++){

                  double x_j = nm_atom[atom_i].x*a0x + i*a0x;
                  double y_j = nm_atom[atom_i].y*a0y + j*a0y;
                  double z_j = nm_atom[atom_i].z*c0 + k*c0;
                  //std::cout << z_j << '\t' << twist_loction <<std::endl;
                  if ( z_j > twist_loction){
                     double x_new = x_j*cos(twist_angle) - y_j*sin(twist_angle);
                     double y_new = y_j*cos(twist_angle) + x_j*sin(twist_angle);
                     x_j = x_new;
                     y_j = y_new;
                  }

                  if (inside_system(x_j, y_j, 0.0)){
                     spin new_atom;
                     new_atom.x = x_j;
                     new_atom.y = y_j;
                     new_atom.z = z_j;
                     new_atom.S = 1;
                     new_atom.id = total_atoms;
                     all_nm_atoms.push_back(new_atom);
                  //   std::cout << total_nm_atoms << "\t" << x_j/system_size_x << '\t' <<  y_j/system_size_y <<  "\t" << z_j/system_size_z << "\t" << 1 << "\t" << 0 << "\t" << 0 << std::endl;
                     total_nm_atoms ++;
                  }
            }
         }
      }
   }
   std::cout << total_nm_atoms <<std::endl;// "\t" << row2.size() << "\t" << row3.size() << "\t" << row4.size() << std::endl;
}
