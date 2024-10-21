#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "initialise.hpp"
#include "exchange.hpp"

// System headers
#include <chrono>

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

// simple exchange selection function (for brevity)
double exchange_factor(int layer1, int layer2){

   if(layer1 == 0 && layer2 == 1) return exchange12;
   if(layer1 == 1 && layer2 == 0) return exchange12;

   if(layer1 == 1 && layer2 == 2) return exchange23;
   if(layer1 == 2 && layer2 == 1) return exchange23;

   if(layer1 == 2 && layer2 == 3) return exchange34;
   if(layer1 == 3 && layer2 == 2) return exchange34;

   return 0.0;
}

// simple exchange selection function (for brevity)
double dmi_factor(int layer1, int layer2){

   if(layer1 == 0 && layer2 == 1) return dmi12;
   if(layer1 == 1 && layer2 == 0) return dmi12;

   if(layer1 == 1 && layer2 == 2) return dmi23;
   if(layer1 == 2 && layer2 == 1) return dmi23;

   if(layer1 == 2 && layer2 == 3) return dmi34;
   if(layer1 == 3 && layer2 == 2) return dmi34;

   return 0.0;
}

//-------------------------------------------
// function to create neighbour list i <-> j
//-------------------------------------------
std::vector<std::vector <int> > generate_neighbours(const double range, std::vector < spin >& atom_list_1){

   // calculate max range squred
   const double r2 = range*range;

   // data structure to hold neighbour list
   std::vector<std::vector <int> > nn_list;
   nn_list.resize(atom_list_1.size());

   // set block size as 1.5*max range
   const double bsize = 1.5*range;

   // calculate min and max xyz
   double min[3] = {1.0e6, 1.0e6, 1.0e6};
   double max[3] = {-1.0e6, -1.0e6, -1.0e6};
   for(int i=0; i < atom_list_1.size(); i++){
      double x_i = atom_list_1[i].x;
      double y_i = atom_list_1[i].y;
      double z_i = atom_list_1[i].z;
      if(x_i < min[0]) min[0] = x_i;
      if(y_i < min[1]) min[1] = y_i;
      if(z_i < min[2]) min[2] = z_i;
      if(x_i > max[0]) max[0] = x_i;
      if(y_i > max[1]) max[1] = y_i;
      if(z_i > max[2]) max[2] = z_i;
   }

   //std::cout << "min: " << min[0] << "\t" << min[1] << "\t" << min[2] << std::endl;
   //std::cout << "max: " << max[0] << "\t" << max[1] << "\t" << max[2] << std::endl;

   // determine number of blocks in x,y,z
   const int xb = ceil((max[0]-min[0])/bsize)+1;
   const int yb = ceil((max[1]-min[1])/bsize)+1;
   const int zb = ceil((max[2]-min[2])/bsize)+1;

   // create 4D array to generate blocks
   std::vector< std::vector < std::vector < std::vector <int> > > > boxes;
   boxes.resize(xb);
   for(int i=0; i<xb; i++){
      boxes[i].resize(yb);
      for(int j=0; j<yb; j++){
         boxes[i][j].resize(zb);
      }
   }

   // determine boxid of each atom and save atoms in boxes
   for(int i=0; i < atom_list_1.size(); i++){
      double x_i = atom_list_1[i].x - min[0];
      double y_i = atom_list_1[i].y - min[1];
      double z_i = atom_list_1[i].z - min[2];
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

      // add atom to box list
      boxes[bxi][byi][bzi].push_back(i);

   }

   std::cout << "Generating neighbours..." << std::flush;
   // std::vector <int> exchange_count(atom_list_1.size(), 0);
   // int atom_index = 0;
   // now calculate neighbour list looping over boxes
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
                     // only calculate neighbours for all x,y,z indices ok
                     if(x_ok && y_ok && z_ok){
                        // loop over all atoms in main box
                        for(int ai = 0; ai < boxes[i][j][k].size(); ai++){
                           // atom_index++;
                           // get atom number i
                           const int atom_i = boxes[i][j][k][ai];
                           const double x_i = atom_list_1[atom_i].x;
                           const double y_i = atom_list_1[atom_i].y;
                           const double z_i = atom_list_1[atom_i].z;

                           // loop over all atoms in neighbour box
                           for(int aj = 0; aj < boxes[nx][ny][nz].size(); aj++){

                              // get atom number j
                              const int atom_j = boxes[nx][ny][nz][aj];

                              // calculate distance
                              const double x_j = atom_list_1[atom_j].x;
                              const double y_j = atom_list_1[atom_j].y;
                              const double z_j = atom_list_1[atom_j].z;
                              const double adx = x_i - x_j;
                              const double ady = y_i - y_j;
                              const double adz = z_i - z_j;
                              double dL2 = adx*adx + ady*ady+adz*adz;
                              // check for atoms in interaction range, if so add to neighbour list
                              if(dL2 < r2 && atom_i != atom_j){
                                 nn_list[atom_i].push_back(atom_j);
                                 // exchange_count[atom_index]++;
                              }

                           } // end of j atom loop

                        } // end of i atom loop

                     } // end of protection statement

                  }
               }
            }// end of offset loops

         }
      }
   } // end of box loops
   // for(int i = 1; i < atom_list_1.size(); i++){
   //    if(exchange_count[i-1] != exchange_count[i] ) std::cout << i << " has " << exchange_count[i] << " vs " << i-1 << " with " << exchange_count[i-1] << std::endl;
   // }
   std::cout << "done!" << std::endl;

   return nn_list;

}

//-----------------------------------------------
// function to create neighbour list i -> j only
//-----------------------------------------------
std::vector<std::vector <int> > generate_neighbours(const double range, std::vector < spin >& atom_list_1, std::vector < spin >& atom_list_2){

   // calculate max range squred
   const double r2 = range*range;

   // data structure to hold neighbour list
   std::vector<std::vector <int> > nn_list;
   nn_list.resize(atom_list_1.size());

   // set block size as 1.5*max range
   const double bsize = 1.5*range;

   // calculate min and max xyz for both lists
   double min[3] = {1.0e6, 1.0e6, 1.0e6};
   double max[3] = {-1.0e6, -1.0e6, -1.0e6};
   for(int i=0; i < atom_list_1.size(); i++){
      double x_i = atom_list_1[i].x;
      double y_i = atom_list_1[i].y;
      double z_i = atom_list_1[i].z;
      if(x_i < min[0]) min[0] = x_i;
      if(y_i < min[1]) min[1] = y_i;
      if(z_i < min[2]) min[2] = z_i;
      if(x_i > max[0]) max[0] = x_i;
      if(y_i > max[1]) max[1] = y_i;
      if(z_i > max[2]) max[2] = z_i;
   }
   for(int i=0; i < atom_list_2.size(); i++){
      double x_i = atom_list_2[i].x;
      double y_i = atom_list_2[i].y;
      double z_i = atom_list_2[i].z;
      if(x_i < min[0]) min[0] = x_i;
      if(y_i < min[1]) min[1] = y_i;
      if(z_i < min[2]) min[2] = z_i;
      if(x_i > max[0]) max[0] = x_i;
      if(y_i > max[1]) max[1] = y_i;
      if(z_i > max[2]) max[2] = z_i;
   }

   //std::cout << "min: " << min[0] << "\t" << min[1] << "\t" << min[2] << std::endl;
   //std::cout << "max: " << max[0] << "\t" << max[1] << "\t" << max[2] << std::endl;

   // determine number of blocks in x,y,z
   const int xb = ceil((max[0]-min[0])/bsize)+1;
   const int yb = ceil((max[1]-min[1])/bsize)+1;
   const int zb = ceil((max[2]-min[2])/bsize)+1;

   // create 4D array to generate blocks
   std::vector< std::vector < std::vector < std::vector <int> > > > boxes1, boxes2;
   boxes1.resize(xb);
   boxes2.resize(xb);
   for(int i=0; i<xb; i++){
      boxes1[i].resize(yb);
      boxes2[i].resize(yb);
      for(int j=0; j<yb; j++){
         boxes1[i][j].resize(zb);
         boxes2[i][j].resize(zb);
      }
   }

   // determine boxid of each atom and save atoms in boxes
   for(int i=0; i < atom_list_1.size(); i++){
      double x_i = atom_list_1[i].x - min[0];
      double y_i = atom_list_1[i].y - min[1];
      double z_i = atom_list_1[i].z - min[2];
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

      // add atom to box list
      boxes1[bxi][byi][bzi].push_back(i);

   }
   for(int i=0; i < atom_list_2.size(); i++){
      double x_i = atom_list_2[i].x - min[0];
      double y_i = atom_list_2[i].y - min[1];
      double z_i = atom_list_2[i].z - min[2];
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

      // add atom to box list
      boxes2[bxi][byi][bzi].push_back(i);

   }

   std::cout << "Generating neighbours..." << std::flush;

   // now calculate neighbour list looping over boxes
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
                     // only calculate neighbours for all x,y,z indices ok
                     if(x_ok && y_ok && z_ok){
                        // loop over all atoms in box 1
                        for(int ai = 0; ai < boxes1[i][j][k].size(); ai++){

                           // get atom number i
                           const int atom_i = boxes1[i][j][k][ai];
                           const double x_i = atom_list_1[atom_i].x;
                           const double y_i = atom_list_1[atom_i].y;
                           const double z_i = atom_list_1[atom_i].z;

                           // loop over all atoms in neighbour box 2
                           for(int aj = 0; aj < boxes2[nx][ny][nz].size(); aj++){

                              // get atom number j
                              const int atom_j = boxes2[nx][ny][nz][aj];

                              // calculate distance
                              const double x_j = atom_list_2[atom_j].x;
                              const double y_j = atom_list_2[atom_j].y;
                              const double z_j = atom_list_2[atom_j].z;
                              const double adx = x_i - x_j;
                              const double ady = y_i - y_j;
                              const double adz = z_i - z_j;
                              double dL2 = adx*adx + ady*ady+adz*adz;
                              // check for atoms in interaction range, if so add to neighbour list
                              if(dL2 < r2){
                                 nn_list[atom_i].push_back(atom_j);
                              }

                           } // end of j atom loop


                        } // end of i atom loop

                     } // end of protection statement

                  }
               }
            }// end of offset loops

         }
      }
   } // end of box loops

   std::cout << "done!" << std::endl;

   return nn_list;

}

//set nearest neighbour  distances (in plane nn 1,2,3)
double intra_nn_dist_1 = 4.01; //A
double intra_nn_dist_2 = 6.96; //A
double intra_nn_dist_3 = 8.005; //A

double inter_nn_dist_1 = 7.0;
double inter_nn_dist_2 = 7.77;
double inter_nn_dist_3 = 9.9;

double inter_AB_dist_1 = 6.541;
double inter_AB_dist_2 = 7.77;
double inter_AB_dist_3 = 7.77;
// double nn_dist_3 = a0x*pow(1.3333333333,0.5);
double nn_dist_1;
double nn_dist_2;
double nn_dist_3;

double max_range = 9.9;
//Set exchange interaction values and associated constants
double eVtoJ = 1.602176634e-19;
double J_constant = 1.0*eVtoJ/1000.0; //1 meV
double J_intra_1=2.5*J_constant;
double J_intra_2=0.75*J_constant;
double J_intra_3=-0.01*J_constant;

double Jinter1_AB = 0.135*J_constant;
double Jinter2_AB = 0.185*J_constant;
double Jinter2_AB_prime = 0.15*J_constant;
double Jintra1_AB = 1.98*J_constant;
double Jintra2_AB = 0.275*J_constant;
double Jintra2_ABprime = 0.305*J_constant;

double D_intra_x_constant = 0.025*J_constant;
double D_intra_y_constant =  0.03*J_constant;
double D_intra_z_constant = 0.01*J_constant;
double D_intra2_x_constant = 0.0125*J_constant;
double D_intra2_y_constant = 0.015*J_constant;
double D_intra2_z_constant = 0.015*J_constant;
//set the initial jumber of interactions to zero for counter
int number_of_interactions = 0;

//initialise arrays to store exchage interactions
std::vector < std::vector < double > > Jint;
std::vector < std::vector < double > > Jinter;
std::vector < std::vector < double > > Einter_Cr1;
std::vector < std::vector < double > > Einter_Cr2;
std::vector < std::vector < double > > Einter_Cr3;
std::vector < std::vector < double > > Einter_Cr4;

std::vector < std::vector < std::vector<std::vector<double > > > > Eintra_Cr1_1NN;
std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr2_1NN;
std::vector < std::vector < std::vector<std::vector<double > > > > Eintra_Cr3_1NN;
std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr4_1NN;

std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr1_2NN;
std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr2_2NN;
std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr3_2NN;
std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr4_2NN;

std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr1_3NN;
std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr2_3NN;
std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr3_3NN;
std::vector < std::vector < std::vector<std::vector<double> > > > Eintra_Cr4_3NN;


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
std::vector<std::vector<std::vector<double> > > config_energy(number_of_unit_cells_x, std::vector<std::vector<double> >(number_of_unit_cells_y, std::vector<double>(30,0.0)));//(number_of_unit_cells_x, std::vector<std::array<double, 30> >(number_of_unit_cells_y, {0.0}));

std::vector <double > crossProduct(std::vector <double >A, std::vector <double > B){
   std::vector <double > P(3,0.0);
    P[0] = A[1] * B[2] - A[2] * B[1];
    P[1] = A[2] * B[0] - A[0] * B[2];
    P[2] = A[0] * B[1] - A[1] * B[0];
    return P;
}

void print_interaction_header(){
   std::ofstream outfile3 ("header_interactions.ucf");
   outfile3 << number_of_interactions <<  "\t" << "tensorial"<< std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
   
}

// double match_inter_exchange(double dx, double dy, std::vector<std::vector<double> > &Eij){
//    double new_shift_error = 10.0;
//    double old_shift_error = 10.0;
//    int min_index = -1;
//    #pragma omp parallel for num_threads(2) reduction(min, min_index)
//    for(int i = 0; i < Eij.size(); i++) {
//       new_shift_error = sqrt((Eij[i][0] - dx)*(Eij[i][0] - dx) + (Eij[i][1] - dy)*(Eij[i][1] - dy));
//       if(new_shift_error < old_shift_error) {
//          old_shift_error = new_shift_error;
//          min_index = i;
//       }
//    }
//    return Eij.at(min_index)[2];
// }


std::vector <double >  calculate_dmi_vector(double x_i, double y_i, double z_i,double x_j, double y_j, double z_j, double dx, double dy, double dz){

   std::vector <double > Duv(3,0.0);
    for (int nm = 0; nm < total_nm_atoms; nm ++){
       double x_nm = all_nm_atoms[nm].x;
       double y_nm = all_nm_atoms[nm].y;
       double z_nm = all_nm_atoms[nm].z;
       double dl1 = (x_nm - x_i)*(x_nm - x_i) + (y_nm - y_i)*(y_nm - y_i) + (z_nm - z_i)*(z_nm - z_i);
       double dl2 = (x_nm - x_j)*(x_nm - x_j) + (y_nm - y_j)*(y_nm - y_j) + (z_nm - z_j)*(z_nm - z_j);
      // std::cout << dl1 << '\t' << dl2 << "\t" << nn_dist_1*nn_dist_1 << std::endl;
       if (dl1  <nn_dist_2*nn_dist_2 && dl2 < nn_dist_2*nn_dist_2){
          std::vector <double > uij(3,0.0);
          std::vector <double > zij(3,0.0);
          std::vector <double > Dij(3,0.0);
         // std::cout << "a" << std::endl;
         double norm = sqrt(dx*dx + dy*dy + dz*dz);
         uij[0] = dx/norm;
         uij[1] = dy/norm;
         uij[2] = dz/norm;
         double midx = ((x_i + x_j)/2);
         double midy = ((y_i + y_j)/2);
         double midz = ((z_i + z_j)/2);
         //std::cout << atom << "\t" << atomi << "\t"<< midx << "\t" <<  midy << "\t" <<  midz/c0 << "\t" <<  x_i << "\t" <<  y_i << "\t" <<  z_i << "\t" << x_j << "\t" <<  y_j << "\t" <<  z_j<< "\t" <<  std::endl;

         zij[0] = midx - x_nm;
         zij[1] = midy - y_nm;
         zij[2] = midz - z_nm;
         double norm2 = sqrt(zij[0]*zij[0] + zij[1]*zij[1] + zij[2]*zij[2]);
         zij[0] = zij[0]/norm2;
         zij[1] = zij[1]/norm2;
         zij[2] = zij[2]/norm2;
         Dij = crossProduct(uij,zij);
         Duv[0] = Duv[0] + Dij[0];
         Duv[1] = Duv[1] + Dij[1];
         Duv[2] = Duv[2] + Dij[2];
         //std::cout << uij[0] << "\t" <<  uij[1] << "\t" <<  uij[2] << "\t" <<  zij[0] << "\t" <<  zij[1] << "\t" <<  zij[2] << "\t" <<  Dij[0] << "\t" <<  Dij[1] << "\t" <<  Dij[2] << "\t" <<  std::endl;
         //std::cout << x_i << "\t" << atomi << "\t" << nm <<  "\t" << Dij[0] << '\t' << Dij[1] << "\t" << Dij[2] << std::endl;
         //std::cout << atom << "\t" << atomi << "\t" << nm <<  "\t" << sqrt(dl1) << '\t' << sqrt(dl2) << std::endl;
      }
   }
   double Dl = sqrt(Duv[0]*Duv[0] + Duv[1]*Duv[1] + Duv[2]*Duv[2]);
   Duv[0] = Duv[0]/Dl;
   Duv[1] = Duv[1]/Dl;
   Duv[2] = Duv[2]/Dl;
      return Duv;

}

// calculate DMI using explit atom lists
std::vector <double >  calculate_dmi_vector(int atom_i, int atom_j,
                                            double x_i, double y_i, double z_i,
                                            double x_j, double y_j, double z_j,
                                            double dx, double dy, double dz,
                                            std::vector< std::vector <int> >& nm_list1,
                                            std::vector< std::vector <int> >& nm_list2){

   // calculate a list of mutual nm atoms shared between i and j
   std::vector<int> nm_atom_list;
   for( int i=0 ; i < nm_list1[atom_i].size() ; i++){
      int ai = nm_list1[atom_i][i];
      for( int j=0 ; j < nm_list2[atom_j].size() ; j++){
         int aj = nm_list2[atom_j][j];
         if(ai == aj) nm_atom_list.push_back(ai);
      }
   }

   std::vector <double > Duv(3,0.0);
   //for (int nm = 0; nm < total_nm_atoms; nm ++){
   for(int i=0; i<nm_atom_list.size(); i++){

      int nm = nm_atom_list[i];

       double x_nm = all_nm_atoms[nm].x;
       double y_nm = all_nm_atoms[nm].y;
       double z_nm = all_nm_atoms[nm].z;
       double dl1 = (x_nm - x_i)*(x_nm - x_i) + (y_nm - y_i)*(y_nm - y_i) + (z_nm - z_i)*(z_nm - z_i);
       double dl2 = (x_nm - x_j)*(x_nm - x_j) + (y_nm - y_j)*(y_nm - y_j) + (z_nm - z_j)*(z_nm - z_j);
      // std::cout << dl1 << '\t' << dl2 << "\t" << nn_dist_1*nn_dist_1 << std::endl;
       if (dl1  <nn_dist_2*nn_dist_2 && dl2 < nn_dist_2*nn_dist_2){
          std::vector <double > uij(3,0.0);
          std::vector <double > zij(3,0.0);
          std::vector <double > Dij(3,0.0);
         // std::cout << "a" << std::endl;
         double norm = sqrt(dx*dx + dy*dy + dz*dz);
         uij[0] = dx/norm;
         uij[1] = dy/norm;
         uij[2] = dz/norm;
         double midx = ((x_i + x_j)/2);
         double midy = ((y_i + y_j)/2);
         double midz = ((z_i + z_j)/2);
         //std::cout << atom << "\t" << atomi << "\t"<< midx << "\t" <<  midy << "\t" <<  midz/c0 << "\t" <<  x_i << "\t" <<  y_i << "\t" <<  z_i << "\t" << x_j << "\t" <<  y_j << "\t" <<  z_j<< "\t" <<  std::endl;

         zij[0] = midx - x_nm;
         zij[1] = midy - y_nm;
         zij[2] = midz - z_nm;
         double norm2 = sqrt(zij[0]*zij[0] + zij[1]*zij[1] + zij[2]*zij[2]);
         zij[0] = zij[0]/norm2;
         zij[1] = zij[1]/norm2;
         zij[2] = zij[2]/norm2;
         Dij = crossProduct(uij,zij);
         Duv[0] = Duv[0] + Dij[0];
         Duv[1] = Duv[1] + Dij[1];
         Duv[2] = Duv[2] + Dij[2];
         //std::cout << uij[0] << "\t" <<  uij[1] << "\t" <<  uij[2] << "\t" <<  zij[0] << "\t" <<  zij[1] << "\t" <<  zij[2] << "\t" <<  Dij[0] << "\t" <<  Dij[1] << "\t" <<  Dij[2] << "\t" <<  std::endl;
         //std::cout << x_i << "\t" << atomi << "\t" << nm <<  "\t" << Dij[0] << '\t' << Dij[1] << "\t" << Dij[2] << std::endl;
         //std::cout << atom << "\t" << atomi << "\t" << nm <<  "\t" << sqrt(dl1) << '\t' << sqrt(dl2) << std::endl;
      }
   }
   double Dl = sqrt(Duv[0]*Duv[0] + Duv[1]*Duv[1] + Duv[2]*Duv[2]);
   Duv[0] = Duv[0]/Dl;
   Duv[1] = Duv[1]/Dl;
   Duv[2] = Duv[2]/Dl;
      return Duv;

}

void calc_in_plane_exchange(std::vector < spin > atom_list_1){

   int num_atoms_1 = atom_list_1.size();
   //std::cout << num_atoms_1 <<std::endl;

   std::stringstream ss;

   // calculate neighbour lists
   std::vector< std::vector <int> > nn_list = generate_neighbours(nn_dist_3 + 0.01, atom_list_1);
   // generate non-magnetic lists
   // std::vector< std::vector <int> > nm_list1 = generate_neighbours(nn_dist_2 + 0.01, atom_list_1, all_nm_atoms);

   for (int atom_i = 0; atom_i < num_atoms_1; atom_i++){


      double x_i = atom_list_1[atom_i].x;
      double y_i = atom_list_1[atom_i].y;
      double z_i = atom_list_1[atom_i].z;
      int id_i = atom_list_1[atom_i].id;
      double changex = atom_list_1[atom_i].dx;
      double changey = atom_list_1[atom_i].dy;
      //std::cout << x_i << '\t' << y_i << "\t" << z_i << "\t" << id_i << std::endl;
      double Dx = -Dx_intra[changex][changey]*J_constant;
      double Dy = -Dy_intra[changex][changey]*J_constant;
      double Dz = -Dz_intra[changex][changey]*J_constant;
      int nn = 0;
      int nnn = 0;
      int nnnn = 0;

      // magic loop over neighbour list
      for(int nni = 0; nni < nn_list[atom_i].size(); nni++){
      //for (int atom_j = 0; atom_j < num_atoms_1; atom_j ++){

         int atom_j = nn_list[atom_i][nni];
         int id_j = atom_list_1[atom_j].id;

         //if (id_i != id_j){
            double x_j = atom_list_1[atom_j].x;
            double y_j = atom_list_1[atom_j].y;
            double z_j = atom_list_1[atom_j].z;
            double dx = x_i - x_j;
            double dy = y_i - y_j;
            double dz = z_i - z_j;
            double dL = sqrt(dx*dx + dy*dy+dz*dz);
            double D_v = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
            D_v = 0.0; // turn off IP DMI
            if (dL < nn_dist_1 + 0.01){
            //            std::cout << "JIJ" << "\t" <<x_j << '\t' << y_j << "\t" << z_j << "\t" << id_j << std::endl;
            // std::vector <double > Dij = calculate_dmi_vector(atom_i, atom_j, x_i, y_i,z_i, x_j,y_j,z_j,dx,dy,dz, nm_list1, nm_list1);
            //std::cout << Dij[0] << '\t' << Dij[1] << '\t' << Dij[2] << std::endl;
            ss << number_of_interactions <<  "\t" << id_i << '\t' << id_j <<" 0 0 0 "<<J_intra_1 << "\t" << Dz << "\t" << -Dy << "\t" << -Dz << "\t" << J_intra_1 << "\t" << Dx << "\t" << Dy << "\t" << -Dx << "\t" << J_intra_1 <<  std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;

            //outfile4 << number_of_interactions <<  "\t" << id_i << '\t' << id_j << " 0 0 0 "<<J_intra_1 << "\t" << 0 << "\t" << 0<< "\t" << 0 << "\t" << J_intra_1 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << J_intra_1 <<  std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
               nn++;
               number_of_interactions++;
            }
            else if (dL < nn_dist_2+0.01){
               ss << number_of_interactions <<  "\t" << id_i << '\t' << id_j << " 0 0 0 "<<J_intra_2 << "\t" << 0 << "\t" << 0<< "\t" << 0 << "\t" << J_intra_2 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << J_intra_2 <<  std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
               nnn++ ;
               number_of_interactions++;
            }
            else if (dL < nn_dist_3 + 0.01){
               ss << number_of_interactions <<  "\t" << id_i << '\t' << id_j << " 0 0 0 "<<J_intra_3 << "\t" << 0 << "\t" << 0<< "\t" << 0 << "\t" << J_intra_3 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << J_intra_3 <<  std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
               nnnn++;
               number_of_interactions++;
            }
         //} // end of i !=j if
      }
   //   std::cout << x_i << '\t' << y_i << '\t' << z_i << "\t" << nn << '\t' << nnn << '\t' << nnnn << std::endl;
   }

   // write data to file from stringstream
   vtimer_t timer;
   std::cout << "Writing data to file..." << std::flush;
   timer.start();
   outfile4 << ss.str();
   timer.stop();
   std::cout << "done! [" << timer.elapsed_time() << std::endl;

   return;

}


void calc_out_of_plane_exchange(std::vector < spin > atom_list_1,std::vector < spin > atom_list_2){

   int num_atoms_1 = atom_list_1.size();
   int num_atoms_2 = atom_list_2.size();

   std::stringstream ss;

   // calculate neighbour list
   std::vector< std::vector <int> > nn_list = generate_neighbours(nn_dist_3 + 0.01, atom_list_1, atom_list_2);
   // generate non-magnetic lists
   // std::vector< std::vector <int> > nm_list1 = generate_neighbours(nn_dist_2 + 0.01, atom_list_1, all_nm_atoms);
   // std::vector< std::vector <int> > nm_list2 = generate_neighbours(nn_dist_2 + 0.01, atom_list_2, all_nm_atoms);

   for (int atom_i = 0; atom_i < num_atoms_1; atom_i ++){

      double x_i = atom_list_1[atom_i].x;
      double y_i = atom_list_1[atom_i].y;
      double z_i = atom_list_1[atom_i].z;
      int id_i = atom_list_1[atom_i].id;
      double changex = atom_list_1[atom_i].dx;
      double changey = atom_list_1[atom_i].dy;
      //std::cout << x_i << "\t" << y_i << "\t" << changex << '\t' << changey << '\t' << Jint[changex][changey] << std::endl;
      double Jj = Jint[changex][changey]*J_constant;
      double Dx = -Dx_inter[changex][changey]*J_constant;
      double Dy = -Dy_inter[changex][changey]*J_constant;
      double Dz = -Dz_inter[changex][changey]*J_constant;
      // double D_v = sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
      //std::cout << atom_i << "\t" << x_i << '\t' << y_i << "\t" << z_i << '\t' << changex << '\t' << changey << '\t' << Jj << std::endl;
      int nn = 0;

      // std::vector<double> Dij(3);
      // magic loop over neighbour list
      for(int nni = 0; nni < nn_list[atom_i].size(); nni++){
      //for (int atom_j = 0; atom_j < num_atoms_2; atom_j ++){

         int atom_j = nn_list[atom_i][nni];

         double x_j = atom_list_2[atom_j].x;
         double y_j = atom_list_2[atom_j].y;
         double z_j = atom_list_2[atom_j].z;
         int id_j = atom_list_2[atom_j].id;

         double dx = x_i - x_j;
         double dy = y_i - y_j;
         double dz = z_i - z_j;
         double dL = sqrt(dx*dx + dy*dy+dz*dz);
         //std::cout << dL << '\t' << nn_dist_3 <<std::endl;
         if (dL < nn_dist_3){
            nn++;

            // calculate exchange factor
            const double Jf = exchange_factor(atom_list_2[atom_j].S, atom_list_1[atom_i].S);
            // const double Df = dmi_factor(atom_list_2[atom_j].S, atom_list_1[atom_i].S);
            //std::cout << id_i << std::endl;
//            std::cout << x_i << "\t" << y_i <<'\t' << x_j << "\t" << y_j << "\t" << Dx*Dij[0] << "\t" <<  Dy*Dij[1] << "\t" <<  Dz*Dij[2]<< "\t" <<  Jj << std::endl;
            //std::vector <double > Dij = calculate_dmi_vector(x_i, y_i,z_i, x_j,y_j,z_j,dx,dy,dz); // original version
            // std::vector <double > Dij = calculate_dmi_vector(atom_i, atom_j, x_i, y_i,z_i, x_j,y_j,z_j,dx,dy,dz, nm_list1, nm_list2);
         //   std::cout << id_i << "\t" << id_j << '\t' << x_i << "\t" << y_i <<'\t' << x_j << "\t" << y_j << "\t" << Dx*Dij[0] << "\t" <<  Dy*Dij[1] << "\t" <<  Dz*Dij[2]<< "\t" <<  Jj << std::endl;//"\t" << "\t" << Dij[0] << "\t" << Dij[1] << "\t" << Dij[2] << std::endl;
            //std::vector <double > Dij(3,0.0);
         // std::cout << x_i << "\t" << y_i <<'\t' << x_j << "\t" << y_j << "\t" << Dx << "\t" <<  Dy << "\t" <<  Dz<< "\t" <<  Dij[0] << '\t' << Dij[1] << '\t' << Dij[2]  << std::endl;
         //std::cout << x_i << "\t" << y_i <<'\t' << x_j << "\t" << y_j << "\t" << Dx*Dij[0] << "\t" <<  Dy*Dij[1] << "\t" <<  Dz*Dij[2]<< "\t" <<  Jj << std::endl;

            //std::cout << x_i << "\t" << y_i << "\t" << z_i << "\t" << z_j << "\t" << Jj <<'\t' << Dx << '\t' << Dy << '\t' << Dz << "\t" << sqrt(Dx*Dx + Dy*Dy + Dz*Dz) << std::endl;
            ss << number_of_interactions <<  "\t" << id_i << '\t' << id_j <<" 0 0 0 "<< Jf*Jj << "\t" << Dz << "\t" << -Dy << "\t" << -Dz << "\t" << Jf*Jj << "\t" << Dx << "\t" << Dy << "\t" << -Dx << "\t" << Jf*Jj <<  std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
            number_of_interactions++;

            // Dij = calculate_dmi_vector(atom_i, atom_j, x_j, y_j,z_j, x_i,y_i,z_i,-dx,-dy,-dz, nm_list1, nm_list2);
            // Dij = calculate_dmi_vector(x_j, y_j,z_j, x_i,y_i,z_i,-dx,-dy,-dz); // original version
            ss << number_of_interactions <<  "\t" << id_i << '\t' << id_j <<" 0 0 0 "<< Jf*Jj << "\t" << Dz << "\t" << -Dy << "\t" << -Dz << "\t" << Jf*Jj << "\t" << Dx << "\t" << Dy << "\t" << -Dx << "\t" << Jf*Jj <<  std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
            //std::cout << id_i << "\t" << id_j << '\t' << x_i << "\t" << y_i <<'\t' << x_j << "\t" << y_j << "\t" << Dx*Dij[0] << "\t" <<  Dy*Dij[1] << "\t" <<  Dz*Dij[2]<< "\t" <<  Jj << std::endl;//"\t" << "\t" << Dij[0] << "\t" << Dij[1] << "\t" << Dij[2] << std::endl;

            //outfile4 << number_of_interactions <<  "\t" << id_j << '\t' << id_i <<" 0 0 0 "<<Jj << "\t" << Dz*Dij[2] << "\t" << -Dy*Dij[1]<< "\t" << -Dz*Dij[2] << "\t" << Jj << "\t" << Dx*Dij[0] << "\t" << Dy*Dij[1] << "\t" << -Dx*Dij[0] << "\t" << Jj <<  std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
            number_of_interactions++;
         }
      }
      //std::cout << nn << std::endl;
   }

   // write data to file from stringstream
   // write data to file from stringstream
   vtimer_t timer;
   std::cout << "Writing data to file..." << std::flush;
   timer.start();
   outfile4 << ss.str();
   timer.stop();
   std::cout << "done! [" << timer.elapsed_time() << "]" << std::endl;

   return;

}

 
void calc_interactions() {

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

   const double bsize = 1.5*range;
   std::cout << "Generating Moire unit cell...." << std::flush;
   // calculate min and max xyz
   // double min[3] = {1.0e8, 1.0e8, 1.0e8};
   // double max[3] = {-1.0e8, -1.0e8, -1.0e8};
   // for(int i=0; i < all_m_atoms.size(); i++){
   //    double x_i = all_m_atoms[i].x;
   //    double y_i = all_m_atoms[i].y;
   //    double z_i = all_m_atoms[i].z;
   //    if(x_i < min[0]) min[0] = x_i;
   //    if(y_i < min[1]) min[1] = y_i;
   //    if(z_i < min[2]) min[2] = z_i;
   //    if(x_i > max[0]) max[0] = x_i;
   //    if(y_i > max[1]) max[1] = y_i;
   //    if(z_i > max[2]) max[2] = z_i;
   // }

   // determine number of blocks in x,y,z
    int xb = ceil(system_size_x/bsize)+1;
    int yb = ceil(system_size_y/bsize)+1;
    int zb = ceil(system_size_z/bsize)+1;
   // std::/cout << "decomposed into <" << xb << ", " << yb << ", " << zb << "> boxes...." << std::flush;
   // create 4D array to generate blocks
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
      double x_i = all_m_atoms[i].x;// - min[0];
      double y_i = all_m_atoms[i].y;// - min[1];
      double z_i = all_m_atoms[i].z;// - min[2];
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
   
   // std::vector< int> interactions_list;
   std::ofstream correlation_file;
   correlation_file.open("moire-lattice-constants.txt");
   std::vector<std::array<double, 3> > zero_correlation;
  
   // now calculate neighbour list looping over boxes
   vtimer_t timer;
      timer.start();
   
   for(int i=0; i<xb; i++){
      // if(i%10 == 0) std::cout << "." << std::flush;
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
                     // int i_index = nx;
                     // int j_index = ny;
                     // int k_index = nz;
                     // int pbc_x = 0;
                     // int pbc_y = 0;
                     // int pbc_z = 0;

                     // if(nx < 0) {i_index = xb-1; pbc_x = -1;}
                     // else if (nx >= xb) {i_index = 0; pbc_x = 1;}

                     // if (ny < 0) {j_index = yb-1; pbc_y = -1;}
                     // else if (ny >= yb) {j_index = 0; pbc_y = 1;}

                     // if(nz < 0 || nz >= zb) continue;
                     if(x_ok && y_ok && z_ok){
                     // only calculate neighbours for all x,y,z indices ok
                        // loop over all atoms in main box
                        for(int ai = 0; ai < boxes[i][j][k].size(); ai++){
                           // atom_index++;
                           // get atom number i
                           spin atom_i = boxes[i][j][k][ai];
                           const double x_i = atom_i.x;
                           const double y_i = atom_i.y;
                           const double z_i = atom_i.z;
                           // if(atom_i.id == 0) continue;
                           // loop over all atoms in neighbour box
                           for(int aj = 0; aj < boxes[nx][ny][nz].size(); aj++){

                              // get atom number j
                              spin atom_j = boxes[nx][ny][nz][aj];
                              if((atom_i.id == atom_j.id)) continue;
                              // if(interactions_list[atom_i.id*11 + interactions_list[atom_j.id*11]])
                              // calculate distance
                              const double x_j = atom_j.x;
                              const double y_j = atom_j.y;
                              const double z_j = atom_j.z;
                              double adx = x_j - x_i;
                              double ady = y_j - y_i;
                              // if(adx < -1*range) adx += system_size_x;
                              // else if(adx > 1*range) adx -= system_size_x;

                              // if(ady < -1*range) ady += system_size_y;
                              // else if(ady > 1*range) ady -= system_size_y;

                              const double adz = z_j - z_i;
                              double dL2 = adx*adx + ady*ady + adz*adz;

                              // check for atoms in interaction range, if so add to neighbour list
                              if(dL2 < 6.6*6.6 ){
                                 if(std::abs(adx) < 1e-2 && std::abs(ady) < 1e-2) {
                                    //correlation_file << std::abs(adx) << ", " << std::abs(ady) << ", " << atom_i.id << ", " << atom_i.x << ", " << atom_i.y << ", " << atom_j.id << ", " << atom_j.x << ", " << atom_j.y << std::endl;
                                    zero_correlation.push_back({atom_i.x, atom_i.y, double(atom_i.l_id)});
                                 }
                                 // std::cout << dL2 << ", " << r2 << ", " << x_i << ", " << y_i << ", " << z_i << ", " << x_j << ", " << y_j << ", " << z_j << std::endl;
                              //    double angle_i = atan2(ady,adx);// - twist_angle;// - M_PI*0.5;
                              //    double angle_j = atan2(-ady,-adx);// - twist_angle;
                              //    std::array<double, 4> exchange({0.0,0.0,0.0,0.0});
                              //    if(atom_i.S == atom_j.S) {
                              //       if(atom_i.l_id == 1) {
                              //          if(dL2 < intra_nn_dist_1) {exchange = match_intra_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr1_1NN ); all_m_atoms[atom_i.id].intra1++;}
                              //          else if (dL2 < intra_nn_dist_2) {exchange = match_intra2_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr1_2NN ); all_m_atoms[atom_i.id].intra2++;}
                              //          else if (dL2 < intra_nn_dist_3) {exchange = match_intra_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr1_3NN ); all_m_atoms[atom_i.id].intra3++;}
                              //          else continue;                                 
                              //       } else if (atom_i.l_id == 2) {
                              //          if(dL2 < intra_nn_dist_1) {exchange = match_intra_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr2_1NN ); all_m_atoms[atom_i.id].intra1++;}
                              //          else if (dL2 < intra_nn_dist_2) {exchange = match_intra2_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr2_2NN ); all_m_atoms[atom_i.id].intra2++;}
                              //          else if (dL2 < intra_nn_dist_3) {exchange = match_intra_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr2_3NN ); all_m_atoms[atom_i.id].intra3++;}
                              //          else continue;
                              //       } else if (atom_i.l_id == 3) {
                              //          angle_i += twist_angle;
                              //          angle_j += twist_angle;
                              //          if(dL2 < intra_nn_dist_1) {exchange = match_intra_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr3_1NN ); all_m_atoms[atom_i.id].intra1++;}
                              //          else if (dL2 < intra_nn_dist_2) {exchange = match_intra2_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr3_2NN ); all_m_atoms[atom_i.id].intra2++;}
                              //          else if (dL2 < intra_nn_dist_3) {exchange = match_intra_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr3_3NN ); all_m_atoms[atom_i.id].intra3++;}
                              //          else continue;
                              //          double r_Dx = exchange[1]*cos(-twist_angle)-exchange[2]*sin(-twist_angle);
                              //          double r_Dy = exchange[1]*sin(-twist_angle)+exchange[2]*cos(-twist_angle);
                              //          exchange[1] = r_Dx;
                              //          exchange[2] = r_Dy;
                              //       } else if (atom_i.l_id == 4) {
                              //          angle_i += twist_angle;
                              //          angle_j += twist_angle;
                              //          if(dL2 < intra_nn_dist_1) {exchange = match_intra_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr4_1NN ); all_m_atoms[atom_i.id].intra1++;}
                              //          else if (dL2 < intra_nn_dist_2) {exchange = match_intra2_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr4_2NN );  all_m_atoms[atom_i.id].intra2++;}
                              //          else if (dL2 < intra_nn_dist_3) {exchange = match_intra_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr4_3NN ); all_m_atoms[atom_i.id].intra3++;}
                              //          else continue;
                              //          double r_Dx = exchange[1]*cos(-twist_angle)-exchange[2]*sin(-twist_angle);
                              //          double r_Dy = exchange[1]*sin(-twist_angle)+exchange[2]*cos(-twist_angle);
                              //          exchange[1] = r_Dx;
                              //          exchange[2] = r_Dy;
                              //       } else if(atom_i.S == 1) {
                              //          exchange = calculate_intra_Jani(atom_i, atom_j, dL2, angle_i);
                              //          exchange[3] *= -1;
                              //          // std::cout << atom_i.id << ", " << exchange[3]/J_constant << std::endl;
                              //       } else if(atom_i.S == 4) {
                              //          angle_i -= twist_angle;
                              //          exchange = calculate_intra_Jani(atom_i, atom_j, dL2, angle_i);
                              //       } else continue;
                              //    } else {
                              //       if (atom_j.h_id == 0 || atom_i.h_id == 0) exchange = calculate_inter_Jani(atom_i, atom_j, dL2, angle_i);
                              //       else if(atom_i.l_id == 1) {  
                              //          if(dL2 <= inter_nn_dist_1) {exchange = match_inter_exchange(adx, ady, Einter_Cr1);
                              //          all_m_atoms[atom_i.id].inter1++;}
                              //          else if(dL2 <= inter_nn_dist_2) {exchange = match_inter_exchange(adx, ady, Einter_Cr1);
                              //          all_m_atoms[atom_i.id].inter2++;}
                              //          else if(dL2 <= inter_nn_dist_3) {exchange = match_inter_exchange(adx, ady, Einter_Cr1);
                              //          all_m_atoms[atom_i.id].inter3++;}
                              //          else continue;
                              //       } else if(atom_i.l_id == 2) {
                              //          if(dL2 <= inter_nn_dist_1) {exchange = match_inter_exchange(adx, ady, Einter_Cr2);
                              //          all_m_atoms[atom_i.id].inter1++;}
                              //          else if(dL2 <= inter_nn_dist_2) {exchange = match_inter_exchange(adx, ady, Einter_Cr2);
                              //          all_m_atoms[atom_i.id].inter2++;}
                              //          else if(dL2 <= inter_nn_dist_3) {exchange = match_inter_exchange(adx, ady, Einter_Cr2);
                              //          all_m_atoms[atom_i.id].inter3++;}
                              //          else continue;
                              //       } else if(atom_i.l_id == 3) {
                              //          adx = adx*cos(-twist_angle)-ady*sin(-twist_angle);
                              //          ady = adx*sin(-twist_angle)+ady*cos(-twist_angle);
                              //          if(dL2 <= inter_nn_dist_1) {exchange = match_inter_exchange(adx, ady, Einter_Cr3);
                              //          all_m_atoms[atom_i.id].inter1++;}
                              //          else if(dL2 <= inter_nn_dist_2) {exchange = match_inter_exchange(adx, ady, Einter_Cr3);
                              //          all_m_atoms[atom_i.id].inter2++;}
                              //          else if(dL2 <= inter_nn_dist_3) {exchange = match_inter_exchange(adx, ady, Einter_Cr3);
                              //          all_m_atoms[atom_i.id].inter3++;}
                              //          else continue;
                              //       } else if(atom_i.l_id == 4) {
                              //          adx = adx*cos(-twist_angle)-ady*sin(-twist_angle);
                              //          ady = adx*sin(-twist_angle)+ady*cos(-twist_angle);
                              //          if(dL2 < inter_nn_dist_1) {exchange = match_inter_exchange(adx, ady, Einter_Cr4);
                              //          all_m_atoms[atom_i.id].inter1++;}
                              //          else if(dL2 < inter_nn_dist_2) {exchange = match_inter_exchange(adx, ady, Einter_Cr4);
                              //          all_m_atoms[atom_i.id].inter2++;}
                              //          else if(dL2 < inter_nn_dist_3) {exchange = match_inter_exchange(adx, ady, Einter_Cr4);
                              //          all_m_atoms[atom_i.id].inter3++;}
                              //          else continue;                                      
                              //       } else continue;
                              //    }
                              //       // bond_avg  = calculate_intra_Jani(atom_i, atom_j, dL2, angle);
                              //    // else bond_avg  = calculate_inter_Jani(atom_i, atom_j, dL2, angle);
                              //    // interactions_list[atom_j.id*11 + interactions_list[atom_j.id*11]] = 0;
                              //    // interactions_list[atom_j.id*11]++;
                              //    // exchange_count[atom_index]++;
                              //    if(exchange[0] == 0.0 ) continue;

                              //    config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+0] += 1.0;
                              //    config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+1] += exchange[0]/J_constant;
                              //    config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+2] += exchange[1]/J_constant;
                              //    config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+3] += exchange[2]/J_constant;
                              //    config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+4] += exchange[3]/J_constant;
                              
                              //     interaction new_interaction;
                              //     new_interaction.id_i = atom_i.id;
                              //     new_interaction.id_j = atom_j.id;
                              //     new_interaction.J = exchange[0];
                              //    if(DMI) {
                              //       new_interaction.Dx = exchange[1];
                              //       new_interaction.Dx = exchange[2];
                              //       new_interaction.Dx = exchange[3];
                              //    }
                              //    // interaction_list.push_back(new_interaction);
                              // // if(DMI) {  outfile4 << number_of_interactions <<  "\t" << atom_i.id << '\t' << atom_j.id << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' <<\
                              // //                   //xx                     xy-> Dz                 xz -> -Dy
                              // //                     exchange[0] << "\t" << exchange[3] << "\t" << -exchange[2] << "\t" << \
                              // //                   //yx -> -Dz              yy                      yz -> Dx
                              // //                    -exchange[3] << "\t" << exchange[0] << "\t" <<  exchange[1] << "\t" << \
                              // //                   //zx -> Dy               yz -> -Dx               zz
                              // //                     exchange[2] << "\t" <<-exchange[1] << "\t" <<  exchange[0] << "\n"; }
                              // // else {   outfile4 << number_of_interactions <<  "\t" << atom_i.id << '\t' << atom_j.id << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' <<\
                              // // //xx                     xy-> Dz                 xz -> -Dy
                              // //    exchange[0] << "\t" << 0.0 << "\t" << 0.0 << "\t" << \
                              // // //yx -> -Dz              yy                      yz -> Dx
                              // //    0.0 << "\t" << exchange[0] << "\t" <<  0.0 << "\t" << \
                              // // //zx -> Dy               yz -> -Dx               zz
                              // //    0.0 << "\t" << 0.0 << "\t" <<  exchange[0] << "\n"; }
                              
                              //    // std::cout << number_of_interactions <<  "\t" << adx << '\t' << ady <<" 0 0 0 "<<\
                              //    //                // xx                     xy-> Dz                 xz -> -Dy 
                              //    //                  exchange[0]/J_constant << "\t" << exchange[3]/J_constant << "\t" << -exchange[2]/J_constant << "\t" << \
                              //    //                // yx -> -Dz              yy                      yz -> Dx
                              //    //                 -exchange[3]/J_constant << "\t" << exchange[0]/J_constant << "\t" <<  exchange[1]/J_constant << "\t" << \
                              //    //                // zx -> Dy               yz -> -Dx               zz
                              //    //                  exchange[2]/J_constant << "\t" <<-exchange[1]/J_constant << "\t" <<  exchange[0]/J_constant << std::endl;}
                              //    number_of_interactions++;                          
                              }
                           } // end of j atom loop

                        } // end of i atom loop

                     } // end of protection statement
                  }
               }
            }// end of offset loops

         }
      }
   }
 
 /*     
   if(ucf_file.is_open()) {
      std::cout << "ucf file add on has been selected. reading secondary file..." << std::endl;
      std::cout << "Warning. DMI rotation has not been designed yet." << std::endl;
      std::cout << "Warning. Lattice basis unity has not been ensured yet." << std::endl;
      std::string line;
      // keep record of current line
      unsigned int line_counter=0;
      unsigned int line_id=0;

      std::string exchange_type_string; // string defining exchange type

      // defaults for interaction list
      unsigned int interaction_range = 1; // assume +-1 unit cell as default

      // Loop over all lines
      while (! ucf_file.eof() ){
         line_counter++;
         // read in whole line
         std::string line;
         getline(ucf_file,line);
         //std::cout << line.c_str() << std::endl;

         // ignore blank lines
         std::string empty="";
         if(line==empty) continue;

         // set character triggers
         const char* hash="#";	// Comment identifier

         bool has_hash=false;
         // Determine if line is a comment line
         for(unsigned int i=0;i<line.length();i++){
            char c=line.at(i);

            if(c== *hash){
                  has_hash=true;
                  break;
            }
         }
         // if hash character found then read next line
         if(has_hash==true) continue;

         // convert line to string stream
         std::istringstream iss(line,std::istringstream::in);

         // non-comment line found - check for line number
         switch(line_id){
            case 0:
               // iss >> unit_cell.dimensions[0] >> unit_cell.dimensions[1] >> unit_cell.dimensions[2];
               break;
            case 1:
               // iss >> unit_cell.shape[0][0] >> unit_cell.shape[0][1] >> unit_cell.shape[0][2];
               break;
            case 2:
               // iss >> unit_cell.shape[1][0] >> unit_cell.shape[1][1] >> unit_cell.shape[1][2];
               break;
            case 3:
               // iss >> unit_cell.shape[2][0] >> unit_cell.shape[2][1] >> unit_cell.shape[2][2];
               break;
            case 4:
               int num_uc_atoms;
               iss >> num_uc_atoms;
               //std::cout << "Reading in " << num_uc_atoms << " atoms" << std::endl;
               // resize unit_cell.atom array if within allowable bounds
               if( (num_uc_atoms >0) && (num_uc_atoms <= 1000000)) all_m_atoms.reserve(num_uc_atoms + all_m_atoms.size());
               else {
                  // terminaltextcolor(RED);
                  std::cerr << "Error! Requested number of atoms " << num_uc_atoms << " on line " << line_counter
                  << " of unit cell input file is outside of valid range 1-1,000,000. Exiting" << std::endl; exit(1);
                  // terminaltextcolor(WHITE);
               }

               std::cout << "\nProcessing data for " << all_m_atoms.size() << " atoms..." << std::flush;
               // zlog << zTs() << "\t" << "Processing data for " << unit_cell.atom.size() << " unit cell atoms..." << std::endl;

            // loop over all atoms and read into class
            for(unsigned int i = 0; i < all_m_atoms.size(); i++){

                line_counter++;

                // declare safe temporaries for atom input
                int id=i;
                double cx=2.0, cy=2.0,cz=2.0; // coordinates - default will give an error
                int mat_id=0, lcat_id=0, hcat_id=0; // sensible defaults if omitted
                // get line
                std::string atom_line;
                getline(ucf_file,atom_line);
                std::istringstream atom_iss(atom_line,std::istringstream::in);
                atom_iss >> id >> cx >> cy >> cz >> mat_id >> lcat_id >> hcat_id;
                
                spin new_atom;
                id += all_m_atoms.size();
                if(cx>=0.0 && cx <=1.0) new_atom.x=cx;
                else{
                    // terminaltextcolor(RED);
                    std::cerr << "Error! atom x-coordinate for atom " << id << " on line " << line_counter
                                    << " of unit cell input file is outside of valid range 0.0-1.0. Exiting" << std::endl;
                    // terminaltextcolor(WHITE);
                    // zlog << zTs() << "Error! atom x-coordinate for atom " << id << " on line " << line_counter
                                //  << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
                    exit(1);
                }
                if(cy>=0.0 && cy <=1.0) new_atom.y=cy;
                else{
                    // terminaltextcolor(RED);
                    std::cerr << "Error! atom y-coordinate for atom " << id << " on line " << line_counter
                                    << " of unit cell input file is outside of valid range 0.0-1.0. Exiting" << std::endl;
                    // terminaltextcolor(WHITE);
                    // zlog << zTs() << "Error! atom y-coordinate for atom " << id << " on line " << line_counter
                    // 			     << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
                    exit(1);
                }
                if(cz>=0.0 && cz <=1.0) new_atom.z=cz;
                else{
                    // terminaltextcolor(RED);
                    std::cerr << "Error! atom z-coordinate for atom " << id << " on line " << line_counter
                    << " of unit cell input file is outside of valid range 0.0-1.0. Exiting" << std::endl;
                    // terminaltextcolor(WHITE);
                    // zlog << zTs() << "Error! atom z-coordinate for atom " << id << " on line " << line_counter
                    // 				  << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
                    exit(1);
                }
                new_atom.unit_x = floor((cx +0.0000001)/ a1y);
                // changex += dy_cell*std::abs(a1x);
                new_atom.unit_y = floor((cy +0.0000001)/ a0x);
                new_atom.S = 5;
                new_atom.id = id;
                new_atom.l_id =lcat_id;
                new_atom.h_id = hcat_id;
                all_m_atoms.push_back(new_atom);
                //std::cout << i << "\t" << id << "\t" << cx << "\t" << cy << "\t" << cz << "\t" << mat_id << "\t" << lcat_id << "\t" << hcat_id << std::endl;
            }
            break;
			case 5:{

            int num_interactions = 0; // assume no interactions
            std::string exchange_type_string; // string defining exchange type

            // get number of exchange types
            iss >> num_interactions >> exchange_type_string;

            // loop over all interactions and read into class
            for (int i=0; i<num_interactions; i++){

                // Output progress counter to screen for large interaction counts
                if( (i % (num_interactions/10 + 1)) == 0 && num_interactions > 10000) std::cout << "." << std::flush;

                // declare safe temporaries for interaction input
                int id=i;
                int iatom=-1,jatom=-1; // atom pairs
                int dx=0, dy=0,dz=0; // relative unit cell coordinates
                // get line
                std::string int_line;
                getline(ucf_file,int_line);
                //std::cout << int_line.c_str() << std::endl;
                std::istringstream int_iss(int_line,std::istringstream::in);
                int_iss >> id >> iatom >> jatom >> dx >> dy >> dz;
                //inputfile >> id >> iatom >> jatom >> dx >> dy >> dz;
                line_counter++;
                // check for sane input
                id += num_interactions;
                iatom += all_m_atoms.size();
                jatom += all_m_atoms.size();
                // check for long range interactions
                if(dx*dx+dy*dy+dz*dz > r2) continue;
                interaction ucf_interaction;
                            //xx                     xy-> Dz             xz -> -Dy
                int_iss >> ucf_interaction.xx >> ucf_interaction.xy >> ucf_interaction.xz;
                            //yx -> -Dz              yy                  yz -> Dx
                int_iss >> ucf_interaction.yx >> ucf_interaction.yy >> ucf_interaction.yz;
                            //zx -> Dy               yz -> -Dx           zz
                int_iss >> ucf_interaction.zx >> ucf_interaction.zy >> ucf_interaction.zz;
        
                spin atom_i = all_m_atoms[id];
                spin atom_j = all_m_atoms[id];
                config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+0] += 1.0;
                config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+1] += ucf_interaction.xx/J_constant;
                config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+2] += 0.5*(ucf_interaction.yz-ucf_interaction.zy)/J_constant;
                config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+3] += 0.5*(ucf_interaction.zx-ucf_interaction.xz)/J_constant;
                config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+4] += 0.5*(ucf_interaction.xy-ucf_interaction.yx)/J_constant;
            
                if(DMI) {  outfile4 << number_of_interactions <<  "\t" << atom_i.id << '\t' << atom_j.id <<" 0 0 0 "<<\
                    //xx                     xy-> Dz                 xz -> -Dy
                        ucf_interaction.xx << "\t" << ucf_interaction.xy << "\t" << ucf_interaction.xz << "\t" << \
                    //yx -> -Dz              yy                      yz -> Dx
                        ucf_interaction.yx << "\t" << ucf_interaction.yy << "\t" <<  ucf_interaction.yz << "\t" << \
                    //zx -> Dy               yz -> -Dx               zz
                        ucf_interaction.zx << "\t" << ucf_interaction.yz << "\t" <<  ucf_interaction.zz << "\n"; }
                else {   outfile4 << number_of_interactions <<  "\t" << atom_i.id << '\t' << atom_j.id <<" 0 0 0 "<<\
                    //xx                     xy-> Dz                 xz -> -Dy
                        ucf_interaction.xx << "\t" << 0.0 << "\t" << 0.0 << "\t" << \
                    //yx -> -Dz              yy                      yz -> Dx
                        0.0 << "\t" << ucf_interaction.yy << "\t" <<  0.0 << "\t" << \
                    //zx -> Dy               yz -> -Dx               zz
                        0.0 << "\t" << 0.0 << "\t" <<  ucf_interaction.zz << "\n"; }

                    number_of_interactions++;
                }
            }
	    }
	    line_id++;
    }
    */
      int max_symmetry = 0;
      for(int i = 0; i < zero_correlation.size(); i+=1) {
         double min_x = 1.0;
         double min_y = 1.0;
         double x_vec0 = 0.0;
         double y_vec0 = 0.0;
         double x_vec1 = 0.0;
         double y_vec1 = 0.0;
         double x00 = -1;
         double y00 = -1;
         double x01 = -1;
         double y01 = -1;
         double x10 = -1;
         double y10 = -1;
         double x11 = -1;
         double y11 = -1;
         int max_local_symm = 0;
         double min_r = 1e6;
         bool vector_id0 = false;// = (zero_correlation[j][2] == zero_correlation[i][2]) ? (true) : (false);
         bool vector_id1 = false;// (zero_correlation[j][2] == zero_correlation[i][2]) ? (true) : (false);
         std::vector<std::array<double, 7 > > orthogonal_set;
         for(int j = 0; j < zero_correlation.size(); j+=1) {
            if(i == j) continue;
            double dx = zero_correlation[j][0] - zero_correlation[i][0];
            double dy = zero_correlation[j][1] - zero_correlation[i][1];
            double r = sqrt(dx*dx+dy*dy);
            dx /= r;
            dy /= r;
            if(r < min_r) min_r = r;

            // if( 1.0-std::abs(dx) < min_x) {
               min_x = 1.0-std::abs(dx);
               x_vec0 = dx;
               y_vec0 = dy;
               x00 = zero_correlation[i][0];
               y00 = zero_correlation[i][1];
               x01 = zero_correlation[j][0];
               y01 = zero_correlation[j][1];
               vector_id0 = (zero_correlation[j][2] == zero_correlation[i][2]) ? (true) : (false);
            // }
            // if( 1.0-std::abs(dy) < min_y) {
               min_y = 1.0-std::abs(dy);
               x_vec1 = dx;
               y_vec1 = dy;
               x10 = zero_correlation[i][0];
               y10 = zero_correlation[i][1];
               x11 = zero_correlation[j][0];
               y11 = zero_correlation[j][1];
               vector_id1 = (zero_correlation[j][2] == zero_correlation[i][2]) ? (true) : (false);
            // }
               // double r_vec = sqrt((x_vec0-x_vec1)*(x_vec0-x_vec1) + (y_vec1-y_vec0)*(y_vec1-y_vec0));

               if( vector_id0 && vector_id1) {
                 
                  orthogonal_set.push_back({x_vec0, y_vec0, r, x00, y00, x01, y01});
                  // correlation_file << i << ", " << j << ", " << x_vec0 << ", " << y_vec0 << ", " << x00 << ", " << y00 << ", " << x01 << ", " << y01 << '\n';
                  // correlation_file << i << ", " << j << ", " << x_vec1 << ", " << y_vec1 << ", " << x10 << ", " << y10 << ", " << x11 << ", " << y11 << '\n';
               }
         }
         for(int k = 0; k < orthogonal_set.size(); k++) {
            for(int kprime = k+1; kprime < orthogonal_set.size(); kprime++) {
               if(k == kprime) continue;
               if( std::abs(orthogonal_set[k][0]*orthogonal_set[kprime][0]+orthogonal_set[k][1]*orthogonal_set[kprime][1]) < 1e-4){
                  max_local_symm++;
                  correlation_file << i;
                  for(int o = 0; o < orthogonal_set[k].size(); o++) correlation_file << ", " <<  orthogonal_set[k][o];
                  for(int o = 0; o < orthogonal_set[kprime].size(); o++) correlation_file << ", " << orthogonal_set[kprime][o];
                  correlation_file << "\n";
                  // correlation_file << orthogonal_set[k][0] << ", " << orthogonal_set[k][1] << ", " << orthogonal_set[kprime][0] << ", " << orthogonal_set[kprime][1] << '\n';
               }
            }
         }
      if(max_local_symm > max_symmetry) {
         std::cout << i << " with symm " << max_local_symm << std::endl;
         max_symmetry = max_local_symm;
      }
      }
      correlation_file.close();
   // double x00 =  3662.51;// A 3662.51, 1594.52 
   // double y00 =  1594.52;// 
   // double x01 =  3659.04;// B 3659.04, 3185.05
   // double y01 =  3185.05;// 
   // double x10 =  6417.19;// A;  6417.19, 1600.54 
   // double y10 =  1600.53;// 
   // double x11 =  6413.71;// B' 6413.71, 3191.07
   // double y11 =  3191.06;//

   double x00 =  3439.04;// A 3662.51, 1594.52 
   double y00 =  2779.95;// 
   double x01 =  3439.04;// B 3659.04, 3185.05
   double y01 =  5162.77;// 
   double x10 =  4814.66;// A;  6417.19, 1600.54 
   double y10 =  2779.95;// 
   double x11 =  4814.66;// B' 6413.71, 3191.07
   double y11 =  5162.77;//

   double rAAprime = sqrt((x10-x00)*(x10-x00)+(y10-y00)*(y10-y00));
   double rAB      = sqrt((x01-x00)*(x01-x00)+(y01-y00)*(y01-y00));
   double rAprimeBprime = sqrt((x10-x11)*(x10-x11)+(y10-y11)*(y10-y11));
   double rBBprime = sqrt((x11-x01)*(x11-x01)+(y11-y01)*(y11-y01));

   std::cout << rAAprime << " == " << rBBprime << " != " << rAB << " == " << rAprimeBprime << std::endl;

   double Moire_aix = x10-x00;
   double Moire_aiy = y10-y00;
   double Moire_ajx = x01-x00; 
   double Moire_ajy = y01-y00;

   double Moire_aix_prime = x11-x01;
   double Moire_aiy_prime = y11-y01;
   double Moire_ajx_prime = x11-x10;
   double Moire_ajy_prime = y11-y10;


   std::cout << "Moire lattice shift " << Moire_aix << " == " << Moire_aix_prime << ", " \
                                       << Moire_aiy << " == " << Moire_aiy_prime << ", " \
                                       << Moire_ajx << " == " << Moire_ajx_prime << ", " \
                                       << Moire_ajy << " == " << Moire_ajy_prime <<  std::endl; 
                                      
   Moire_aix = 0.5*(Moire_aix+Moire_aix_prime);
   Moire_aiy = 0.5*(Moire_aiy+Moire_aiy_prime);
   Moire_ajx = 0.5*(Moire_ajx+Moire_ajx_prime);
   Moire_ajy = 0.5*(Moire_ajy+Moire_ajy_prime);

   x10 =  x00 + Moire_aix;//6417.19;// A;  6417.19, 1600.54 
   y10 =  y00 + Moire_aiy;//1600.53;// 
   x11 =  x01 + Moire_aix;//6413.71;// B' 6413.71, 3191.07
   y11 =  y01 + Moire_aiy;//3191.06;//

   std::cout << "Shifting unit cell range [" << x00 << ", " << y00 << "], [" << x11 << ", " << y11 << "]"  << std::endl;
   auto zero_zero_to_zero_one = [&](double x, double y){
      // double x0 = 3669.43;
      // double x1 = 3659.04;
      // double y0 = 10.0089;
      // double y1 = 3185.05;

      double y_prime = 0.00+y00 + ((y01-y00)/(x01-x00))*(x - x00);
      if(y >= y_prime) return true;
      else return false;
   };
   auto zero_zero_to_one_zero = [&](double x, double y){
      // double x0 = 3669.43;
      // double x1 = 6417.19;
      // double y0 = 10.0089;
      // double y1 = 1600.54;

      double y_prime = 0.00+y00 + ((y10-y00)/(x10-x00))*(x - x00);
      if(y >= y_prime) return true;
      else return false;
   };
   auto zero_one_to_one_one = [&](double x, double y){
      // double x0 = 3659.04;
      // double x1 = 7778.92;
      // double y0 = 3185.05;
      // double y1 = 7164.39;

      double y_prime = -0.00 + y01 + ((y11-y01)/(x11-x01))*(x - x01);
      if(y < y_prime) return true;
      else return false;
   };
   auto one_zero_to_one_one = [&](double x, double y){
      // double x0 = 6417.19;
      // double x1 = 7778.92;
      // double y0 = 1600.54;
      // double y1 = 7164.39;

      double y_prime = -0.00 + y10 + ((y11-y10)/(x11-x10))*(x - x10);
      if(y < y_prime) return true;
      else return false;
   };

   double x_offset = 3439.04;// std::min(x00, x01);//-a0x*0.5-0.001;
   double y_offset = 2779.95;//std::min(y00, y10);
   double max_x = 4814.66;//std::max(x10, x11);// - 0.0001;
   double max_y = 5162.77;//std::max(y01, y11);//  - 0.0001; 

   const double new_system_size_x = max_x-x_offset;
   const double new_system_size_y = max_y-y_offset;

   int interaction_estimate = all_m_atoms.size()*new_system_size_x*new_system_size_y/system_size_x/system_size_y;
   all_m_atoms_offset.reserve(interaction_estimate);
   int new_atom_count = 0;
   std::ofstream outfile2;
   outfile2.open("atom_positions.xyz");
   for(int i = 0; i < all_m_atoms.size(); i++) {
      spin offset_atom(all_m_atoms[i]);
      double x = offset_atom.x;
      double y = offset_atom.y;
      bool A_B = zero_zero_to_zero_one(x,y);
      bool A_Aprime = zero_zero_to_one_zero(x,y);
      bool B_Bprime = zero_one_to_one_one(x,y);
      bool Aprime_Bprime = one_zero_to_one_one(x,y);
      bool x_to_xprime = (x >= x_offset && x<(max_x-0.1));
      bool y_to_yprime = (y >= y_offset && y<(max_y-0.1));
      // && A_B && A_Aprime && B_Bprime && Aprime_Bprime
      if(x_to_xprime && y_to_yprime ) {
         offset_atom.id = new_atom_count;
         // offset_atom.x -= x_offset;
         // offset_atom.y -= y_offset;
         all_m_atoms_offset.push_back(offset_atom);
         outfile2 << new_atom_count << "\t" << (offset_atom.x-x_offset)/new_system_size_x << '\t' <<  (offset_atom.y-y_offset)/new_system_size_y <<  "\t" << offset_atom.z/system_size_z << "\t" << offset_atom.S-1 << "\t" << offset_atom.l_id << "\t" << offset_atom.h_id << "\n"; 
         new_atom_count++;
      }
   }
   outfile2.close();

   std::ofstream outfile1 ("header.ucf");

   outfile1 << " #unit cell size " << std::endl;
   outfile1 << new_system_size_x << '\t' << new_system_size_y << '\t' << system_size_z << std::endl;
   outfile1 << " #unit cell vectors" << std::endl;
   outfile1 << " " <<  Moire_aix/new_system_size_x << '\t' << Moire_aiy/new_system_size_x << '\t' <<  "0" << std::endl;
   outfile1 << " " <<  Moire_ajx/new_system_size_y << '\t' << Moire_ajy/new_system_size_y << '\t' <<  "0" << std::endl;
   outfile1 << " 0     0	   1" << std::endl;
   outfile1 << " #Atoms" << std::endl;
   outfile1 << new_atom_count << '\t'	<< 4 << std::endl;

   outfile1.close();

   new_moire_lattice.reserve(total_atoms);

   create_magnetic_atom_list_moire_unit("atom_positions.ucf", Moire_aix,Moire_aiy,\
                                                              Moire_ajx, Moire_ajy, \
                                                              new_system_size_x, new_system_size_y, new_atom_count);

     // determine number of blocks in x,y,z
   // xb = ceil(new_system_size_x/bsize)+1;
   // yb = ceil(new_system_size_y/bsize)+1;
   //   zb = ceil(system_size_z/bsize)+1;
   std::cout << "decomposed into <" << xb << ", " << yb << ", " << zb << "> boxes." << std::endl;
   // create 4D array to generate blocks
   std::vector< std::vector < std::vector < std::vector < spin > > > > new_boxes;
   new_boxes.resize(xb);
   for(int i=0; i<xb; i++){
      new_boxes[i].resize(yb);
      for(int j=0; j<yb; j++){
         new_boxes[i][j].resize(zb);
      }
   }
   int boxed = 0;
   // determine boxid of each atom and save atoms in boxes
   for(int i=0; i < new_moire_lattice.size(); i++){
      double x_i = new_moire_lattice[i].x;// - min[0];
      double y_i = new_moire_lattice[i].y;// - min[1];
      double z_i = new_moire_lattice[i].z;// - min[2];
      const double bxi = x_i / bsize;
      const double byi = y_i / bsize;
      const double bzi = z_i / bsize;

      // if(bxi >= xb) bxi = xb-1;
      // if(byi >= yb) byi = yb-1;

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
      boxed++;
      new_boxes[bxi][byi][bzi].push_back(new_moire_lattice[i]);
   }
   if(boxed != new_moire_lattice.size()) {
      std::cout << "atoms missed in supercell offset " << std::endl;
      exit(1);
   }
   
   interaction_estimate = all_m_atoms_offset.size()*22;
   // std::vector<interaction> interaction_list;
   // interaction_list.reserve(interaction_estimate);


   std::cout << "Generating estimated " << interaction_estimate << " interactions for remaining " << all_m_atoms_offset.size() << " atoms " << std::endl;
      // exit(1);

   for(int i=0; i<xb; i++){
      if(i%10 == 0) std::cout << "." << std::flush;
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
                     int i_index = nx;
                     int j_index = ny;
                     int k_index = nz;
                     int pbc_x = 0;
                     int pbc_y = 0;
                     int pbc_z = 0;

                     // if(nx < 0) {i_index = xb-1; pbc_x = -1;}
                     // else if (nx >= xb) {i_index = 0; pbc_x = 1;}
                     // if (ny < 0) {j_index = yb-1; pbc_y = -1;}
                     // else if (ny >= yb) {j_index = 0; pbc_y = 1;}
                     // if(nz < 0 || nz >= zb) continue;
                     if(x_ok && y_ok && z_ok){
                     // only calculate neighbours for all x,y,z indices ok
                        // loop over all atoms in main box
                        for(int ai = 0; ai < new_boxes[i][j][k].size(); ai++){
                           // atom_index++;
                           // get atom number i
                           spin atom_i = new_boxes[i][j][k][ai];
                           // if(atom_i.Gx != 0 || atom_i.Gy != 0) continue;
                           const double x_i = atom_i.x;
                           const double y_i = atom_i.y;
                           const double z_i = atom_i.z;
                           // if(atom_i.id == 0) continue;
                           // if(pbc_x != 0 || pbc_y != 0 ) std::cout << i_index << ", " <<  j_index << ", " << k_index  << std::endl;
                           // loop over all atoms in neighbour box
                           for(int aj = 0; aj < new_boxes[i_index][j_index][k_index].size(); aj++){

                              // get atom number j
                              spin atom_j = new_boxes[i_index][j_index][k_index][aj];
                              if(atom_i.id == atom_j.id) continue;
                              // if(interactions_list[atom_i.id*11 + interactions_list[atom_j.id*11]])
                              // calculate distance
                              const double x_j = atom_j.x;
                              const double y_j = atom_j.y;
                              const double z_j = atom_j.z;

                              double adx = x_j - x_i;
                              double ady = y_j - y_i;
                              // if(adx < -2*bsize) adx += new_system_size_x;
                              // else if(adx > 2*bsize) adx -= new_system_size_x;

                              // if(ady < -2*bsize) ady += new_system_size_y;
                              // else if(ady > 2*bsize) ady -= new_system_size_y;

                              // if(pbc_x != 0 || pbc_y != 0 ) std::cout << x_j - x_i << ", " << y_j - y_i << ", " <<  adx << ", " << ady << std::endl;
                              const double adz = z_j - z_i;
                              double dL2 = adx*adx + ady*ady + adz*adz;

                              // check for atoms in interaction range, if so add to neighbour list
                              if(dL2 < r2 ){
                                 // if(std::abs(adx) < 1e-2 && std::abs(ady) < 1e-2) {
                                 //    //correlation_file << std::abs(adx) << ", " << std::abs(ady) << ", " << atom_i.id << ", " << atom_i.x << ", " << atom_i.y << ", " << atom_j.id << ", " << atom_j.x << ", " << atom_j.y << std::endl;
                                 //    zero_correlation.push_back({atom_i.x, atom_i.y});
                                 // }
                                 // std::cout << dL2 << ", " << r2 << ", " << x_i << ", " << y_i << ", " << z_i << ", " << x_j << ", " << y_j << ", " << z_j << std::endl;
                                 double angle_i = atan2(ady,adx);// - twist_angle;// - M_PI*0.5;
                                 double angle_j = atan2(-ady,-adx);// - twist_angle;
                                 std::array<double, 4> exchange({0.0,0.0,0.0,0.0});
                                 if(atom_i.S == atom_j.S) {
                                    if(atom_i.l_id == 1) {
                                       angle_i -= twist_angle;
                                       angle_j -= twist_angle;
                                       if(dL2 < intra_nn_dist_1) exchange = match_intra1_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr1_1NN );
                                       else if (dL2 < intra_nn_dist_2) exchange = match_intra2_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr1_2NN ); 
                                       else if (dL2 < intra_nn_dist_3) exchange = match_intra3_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr1_3NN ); 
                                       else continue;     
                                       double r_Dx = exchange[1]*cos(-0.5*twist_angle)-exchange[2]*sin(-0.5*twist_angle);
                                       double r_Dy = exchange[1]*sin(-0.5*twist_angle)+exchange[2]*cos(-0.5*twist_angle);
                                       exchange[1] = r_Dx;
                                       exchange[2] = r_Dy;                            
                                    } else if (atom_i.l_id == 2) {
                                       angle_i -= twist_angle;
                                       angle_j -= twist_angle;
                                       if(dL2 < intra_nn_dist_1) exchange = match_intra1_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr2_1NN ); 
                                       else if (dL2 < intra_nn_dist_2) exchange = match_intra2_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr2_2NN ); 
                                       else if (dL2 < intra_nn_dist_3) exchange = match_intra3_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr2_3NN ); 
                                       else continue;
                                       double r_Dx = exchange[1]*cos(-0.5*twist_angle)-exchange[2]*sin(-0.5*twist_angle);
                                       double r_Dy = exchange[1]*sin(-0.5*twist_angle)+exchange[2]*cos(-0.5*twist_angle);
                                       exchange[1] = r_Dx;
                                       exchange[2] = r_Dy;
                                    } else if (atom_i.l_id == 3) {
                                       angle_i += twist_angle;
                                       angle_j += twist_angle;
                                       if(dL2 < intra_nn_dist_1) exchange = match_intra1_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr3_1NN ); 
                                       else if (dL2 < intra_nn_dist_2) exchange = match_intra2_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr3_2NN ); 
                                       else if (dL2 < intra_nn_dist_3) exchange = match_intra3_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr3_3NN ); 
                                       else continue;
                                       double r_Dx = exchange[1]*cos(0.5*twist_angle)-exchange[2]*sin(0.5*twist_angle);
                                       double r_Dy = exchange[1]*sin(0.5*twist_angle)+exchange[2]*cos(0.5*twist_angle);
                                       exchange[1] = r_Dx;
                                       exchange[2] = r_Dy;
                                    } else if (atom_i.l_id == 4) {
                                       angle_i += twist_angle;
                                       angle_j += twist_angle;
                                       if(dL2 < intra_nn_dist_1) exchange = match_intra1_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr4_1NN ); 
                                       else if (dL2 < intra_nn_dist_2) exchange = match_intra2_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr4_2NN ); 
                                       else if (dL2 < intra_nn_dist_3) exchange = match_intra3_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr4_3NN ); 
                                       else continue;
                                       double r_Dx = exchange[1]*cos(0.5*twist_angle)-exchange[2]*sin(0.5*twist_angle);
                                       double r_Dy = exchange[1]*sin(0.5*twist_angle)+exchange[2]*cos(0.5*twist_angle);
                                       exchange[1] = r_Dx;
                                       exchange[2] = r_Dy;
                                    } else if(atom_i.S == 1) {
                                       exchange = calculate_intra_Jani(atom_i, atom_j, dL2, angle_i);
                                       exchange[3] *= -1;
                                       // std::cout << atom_i.id << ", " << exchange[3]/J_constant << std::endl;
                                    } else if(atom_i.S == 4) {
                                       angle_i -= 0.5*twist_angle;
                                       exchange = calculate_intra_Jani(atom_i, atom_j, dL2, angle_i);
                                    } else continue;
                                 } else {
                                    if (atom_j.h_id == 0 || atom_i.h_id == 0) exchange = calculate_inter_Jani(atom_i, atom_j, dL2, angle_i);
                                    else if(atom_i.l_id == 1) {  
                                       if(dL2 <= inter_nn_dist_1) {exchange = match_inter_exchange(adx, ady, Einter_Cr1);
                                       all_m_atoms_offset[atom_i.id].inter1++;}
                                       else if(dL2 <= inter_nn_dist_2) {exchange = match_inter_exchange(adx, ady, Einter_Cr1);
                                       all_m_atoms_offset[atom_i.id].inter2++;}
                                       else if(dL2 <= inter_nn_dist_3) {exchange = match_inter_exchange(adx, ady, Einter_Cr1);
                                       all_m_atoms_offset[atom_i.id].inter3++;}
                                       else continue;
                                    } else if(atom_i.l_id == 2) {
                                       //  adx = adx*cos(twist_angle)-ady*sin(twist_angle);
                                       // ady = adx*sin(-twist_angle)+ady*cos(twist_angle);
                                       if(dL2 <= inter_nn_dist_1) {exchange = match_inter_exchange(adx, ady, Einter_Cr2);
                                       all_m_atoms_offset[atom_i.id].inter1++;}
                                       else if(dL2 <= inter_nn_dist_2) {exchange = match_inter_exchange(adx, ady, Einter_Cr2);
                                       all_m_atoms_offset[atom_i.id].inter2++;}
                                       else if(dL2 <= inter_nn_dist_3) {exchange = match_inter_exchange(adx, ady, Einter_Cr2);
                                       all_m_atoms_offset[atom_i.id].inter3++;}
                                       else continue;
                                    } else if(atom_i.l_id == 3) {
                                       // adx = adx*cos(-twist_angle)-ady*sin(-twist_angle);
                                       // ady = adx*sin(-twist_angle)+ady*cos(-twist_angle);
                                       if(dL2 <= inter_nn_dist_1) {exchange = match_inter_exchange(adx, ady, Einter_Cr3);
                                       all_m_atoms_offset[atom_i.id].inter1++;}
                                       else if(dL2 <= inter_nn_dist_2) {exchange = match_inter_exchange(adx, ady, Einter_Cr3);
                                       all_m_atoms_offset[atom_i.id].inter2++;}
                                       else if(dL2 <= inter_nn_dist_3) {exchange = match_inter_exchange(adx, ady, Einter_Cr3);
                                       all_m_atoms_offset[atom_i.id].inter3++;}
                                       else continue;
                                    } else if(atom_i.l_id == 4) {
                                       // adx = adx*cos(-twist_angle)-ady*sin(-twist_angle);
                                       // ady = adx*sin(-twist_angle)+ady*cos(-twist_angle);
                                       if(dL2 < inter_nn_dist_1) {exchange = match_inter_exchange(adx, ady, Einter_Cr4);
                                       all_m_atoms_offset[atom_i.id].inter1++;}
                                       else if(dL2 < inter_nn_dist_2) {exchange = match_inter_exchange(adx, ady, Einter_Cr4);
                                       all_m_atoms_offset[atom_i.id].inter2++;}
                                       else if(dL2 < inter_nn_dist_3) {exchange = match_inter_exchange(adx, ady, Einter_Cr4);
                                       all_m_atoms_offset[atom_i.id].inter3++;}
                                       else continue;                                      
                                    } else continue;
                                 }
   
                                 if(exchange[0] == 0.0 ) continue;
                                 config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+0] += 1.0;
                                 config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+1] += exchange[0]/J_constant;
                                 config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+2] += exchange[1]/J_constant;
                                 config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+3] += exchange[2]/J_constant;
                                 config_energy[atom_i.unit_x][atom_i.unit_y][(atom_i.S-1)*5+4] += exchange[3]/J_constant;
                              
                                 //  interaction new_interaction;
                                 //  new_interaction.id_i = atom_i.id;
                                 //  new_interaction.id_j = atom_j.id;
                                 //  new_interaction.J = exchange[0];
                                 // if(DMI) {
                                 //    new_interaction.Dx = exchange[1];
                                 //    new_interaction.Dx = exchange[2];
                                 //    new_interaction.Dx = exchange[3];
                                 // }
                                 // interaction_list.push_back(new_interaction);
                              if(DMI) {  outfile4 << number_of_interactions <<  "\t" << atom_i.id << '\t' << atom_j.id << '\t' << atom_j.Gx << '\t' << atom_j.Gy << '\t' << atom_j.Gz << '\t' <<\
                                                //xx                     xy-> Dz                 xz -> -Dy
                                                  exchange[0] << "\t" << exchange[3] << "\t" << -exchange[2] << "\t" << \
                                                //yx -> -Dz              yy                      yz -> Dx
                                                 -exchange[3] << "\t" << exchange[0] << "\t" <<  exchange[1] << "\t" << \
                                                //zx -> Dy               yz -> -Dx               zz
                                                  exchange[2] << "\t" <<-exchange[1] << "\t" <<  exchange[0] << "\n"; }
                              else {   outfile4 << number_of_interactions <<  "\t" << atom_i.id << '\t' << atom_j.id << '\t' << atom_j.Gx << '\t' << atom_j.Gy << '\t' << atom_j.Gz << '\t' <<\
                              //xx                     xy-> Dz                 xz -> -Dy
                                 exchange[0] << "\t" << 0.0 << "\t" << 0.0 << "\t" << \
                              //yx -> -Dz              yy                      yz -> Dx
                                 0.0 << "\t" << exchange[0] << "\t" <<  0.0 << "\t" << \
                              //zx -> Dy               yz -> -Dx               zz
                                 0.0 << "\t" << 0.0 << "\t" <<  exchange[0] << "\n"; }
                              
                                 // std::cout << number_of_interactions <<  "\t" << adx << '\t' << ady <<" 0 0 0 "<<\
                                 //                // xx                     xy-> Dz                 xz -> -Dy 
                                 //                  exchange[0]/J_constant << "\t" << exchange[3]/J_constant << "\t" << -exchange[2]/J_constant << "\t" << \
                                 //                // yx -> -Dz              yy                      yz -> Dx
                                 //                 -exchange[3]/J_constant << "\t" << exchange[0]/J_constant << "\t" <<  exchange[1]/J_constant << "\t" << \
                                 //                // zx -> Dy               yz -> -Dx               zz
                                 //                  exchange[2]/J_constant << "\t" <<-exchange[1]/J_constant << "\t" <<  exchange[0]/J_constant << std::endl;}
                                 number_of_interactions++;                          
                              }
                           } // end of j atom loop
                        } // end of i atom loop
                     } // end of protection statement
                  }
               }
            }// end of offset loops
         }
      }
   }

      // std::cout << "Writing data to file..." << std::flush;
      std::ofstream config_output;
      config_output.open("config_energy.txt");
      for(int i = 0; i < number_of_unit_cells_x; i++) {
         for(int j = 0; j < number_of_unit_cells_y; j++){
            double bottom_occ = config_energy[i][j][0];
            double top_occ = config_energy[i][j][1];
            if(bottom_occ == 0 && top_occ == 0) continue;
            config_output << i << ", " << j << ", " << bottom_occ<< ", " << top_occ;
            for(int k = 2; k < config_energy[i][j].size(); k++) config_output << ", " << config_energy[i][j][k]; 
            config_output << "\n";
         }
      }
      config_output.close();

      std::ofstream interaction_counts;
      interaction_counts.open("interaction_counts.txt");
      for(int i = 0; i < all_m_atoms_offset.size(); i++){
         interaction_counts << all_m_atoms_offset[i].S  << ", " << all_m_atoms_offset[i].x << ", " << all_m_atoms_offset[i].y <<  ", " <<  all_m_atoms_offset[i].l_id << ", " << all_m_atoms_offset[i].inter1 << ", " << all_m_atoms_offset[i].inter2 << ", " << all_m_atoms_offset[i].inter3 \
                                                   << ", " << all_m_atoms_offset[i].intra1 << ", " << all_m_atoms_offset[i].intra2 << ", " << all_m_atoms_offset[i].intra3 <<"\n";
      }
      interaction_counts.close();
      // outfile4 << ss.str();
      timer.stop();
      // std::cout << "done!  << std::endl;
      std::cout << number_of_interactions << " [completed] [" << timer.elapsed_time() << " s]" << std::endl;
      return;
}

std::array<double,4> match_intra1_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom, std::vector<std::vector< std::vector< std::vector<double> >  > > &Eij){
   std::array<double,4> exchange({0.0});

   int i_x_shift = (unit_cell_shifts[central_atom.unit_x][central_atom.unit_y][1]);
   int i_y_shift = (unit_cell_shifts[central_atom.unit_x][central_atom.unit_y][2]);

   int j_x_shift = (unit_cell_shifts[j_atom.unit_x][j_atom.unit_y][1]);
   int j_y_shift = (unit_cell_shifts[j_atom.unit_x][j_atom.unit_y][2]);

   double theta =  std::abs(angle_i)*180.0/M_PI;
   int theta_i = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0)== 1.0) ? (1) : ( round(theta/30.0) == 3.0 ? (0):-1) );
   theta =  std::abs(angle_j)*180.0/M_PI;
   int theta_j = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0)== 1.0) ? (1) : (round(theta/30.0) == 3.0 ? (0):-1) );
   if(theta_i == -1 || theta_j == -1) {
      std::cout << "\n " << round(std::abs(angle_i)*180.0/M_PI/30.0) << ", " << round(std::abs(angle_j)*180.0/M_PI/30.0) << ", " << angle_i << ", " << angle_j << ", " << central_atom.Gx << ", " << central_atom.Gy << ", " << j_atom.Gx << ", " << j_atom.Gy <<  std::endl;
      exit(1);
      return exchange;
   }
   exchange[0] = 0.5*(Eij.at(i_x_shift).at(i_y_shift)[theta_i][0]+Eij[j_x_shift][j_y_shift][theta_j][0]);
   exchange[1] = 0.5*(Eij[i_x_shift][i_y_shift][theta_i][1]-Eij[j_x_shift][j_y_shift][theta_j][1]);
   exchange[2] = 0.5*(Eij[i_x_shift][i_y_shift][theta_i][2]-Eij[j_x_shift][j_y_shift][theta_j][2]);
   exchange[3] = 0.5*(Eij[i_x_shift][i_y_shift][theta_i][3]-Eij[j_x_shift][j_y_shift][theta_j][3]);
   
   all_m_atoms_offset[central_atom.id].intra1++;
   // std::cout << i_x_shift << ", " << i_y_shift << ", " << j_x_shift << ", " << j_y_shift << ", " << 0.5*(J_i+J_j) << ", " << Dx*cos(twist_angle)-Dy*sin(twist_angle) << ", " <<  Dx*sin(twist_angle)+Dy*cos(twist_angle) << ", " << Dz << std::endl;
   return exchange;
}

std::array<double,4> match_intra2_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom, std::vector<std::vector< std::vector< std::vector<double> >  > > &Eij){
   std::array<double,4> exchange({0.0});
   int i_x_shift = (unit_cell_shifts[central_atom.unit_x][central_atom.unit_y][1]);
   int i_y_shift = (unit_cell_shifts[central_atom.unit_x][central_atom.unit_y][2]);

   int j_x_shift = (unit_cell_shifts[j_atom.unit_x][j_atom.unit_y][1]);
   int j_y_shift = (unit_cell_shifts[j_atom.unit_x][j_atom.unit_y][2]);

   // double theta =  std::abs(angle)*180.0/M_PI;
   int theta_i = int(round(((angle_i < 0.0) ? (angle_i+=2*M_PI) : (angle_i)) *180.0/M_PI/60.0));
   // theta =  std::abs(angle-180.0)*180.0/M_PI;
   int theta_j = int(round(((angle_j < 0.0) ? (angle_j+=2*M_PI) : (angle_j)) *180.0/M_PI/60.0));

   if( theta_i > 5 || theta_j > 5 || theta_i < 0 || theta_j < 0 ) {
      std::cout << "problem: " << theta_i << " , " << theta_j << ", " << angle_i << ", " << angle_j << std::endl;
      exit(1);
      return exchange;
   }
   exchange[0] = 0.5*(Eij.at(i_x_shift).at(i_y_shift).at(theta_i)[0]+Eij[j_x_shift][j_y_shift].at(theta_j)[0]);
   exchange[1] = 0.5*(Eij[i_x_shift][i_y_shift][theta_i][1]-Eij[j_x_shift][j_y_shift][theta_j][1]);
   exchange[2] = 0.5*(Eij[i_x_shift][i_y_shift][theta_i][2]-Eij[j_x_shift][j_y_shift][theta_j][2]);
   exchange[3] = 0.5*(Eij[i_x_shift][i_y_shift][theta_i][3]-Eij[j_x_shift][j_y_shift][theta_j][3]);

   all_m_atoms_offset[central_atom.id].intra2++;

   return exchange;
}

std::array<double,4> match_intra3_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom, std::vector<std::vector< std::vector< std::vector<double> >  > > &Eij){
   std::array<double,4> exchange({0.0});

   int i_x_shift = (unit_cell_shifts[central_atom.unit_x][central_atom.unit_y][1]);
   int i_y_shift = (unit_cell_shifts[central_atom.unit_x][central_atom.unit_y][2]);

   int j_x_shift = (unit_cell_shifts[j_atom.unit_x][j_atom.unit_y][1]);
   int j_y_shift = (unit_cell_shifts[j_atom.unit_x][j_atom.unit_y][2]);

   double theta =  std::abs(angle_i)*180.0/M_PI;
   int theta_i = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0)== 1.0) ? (1) : ( round(theta/30.0) == 3.0 ? (0):-1) );
   theta =  std::abs(angle_j)*180.0/M_PI;
   int theta_j = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0)== 1.0) ? (1) : (round(theta/30.0) == 3.0 ? (0):-1) );
   if(theta_i == -1 || theta_j == -1) {
      // std::cout << "\n " << round(std::abs(angle_i)*180.0/M_PI/30.0) << ", " << round(std::abs(angle_j)*180.0/M_PI/30.0) << ", " << angle_i << ", " << angle_j << ", " << central_atom.Gx << ", " << central_atom.Gy << ", " << j_atom.Gx << ", " << j_atom.Gy <<  std::endl;
      // exit(1);
      return exchange;
   }
   exchange[0] = 0.5*(Eij.at(i_x_shift).at(i_y_shift)[theta_i][0]+Eij[j_x_shift][j_y_shift][theta_j][0]);
   exchange[1] = 0.5*(Eij[i_x_shift][i_y_shift][theta_i][1]-Eij[j_x_shift][j_y_shift][theta_j][1]);
   exchange[2] = 0.5*(Eij[i_x_shift][i_y_shift][theta_i][2]-Eij[j_x_shift][j_y_shift][theta_j][2]);
   exchange[3] = 0.5*(Eij[i_x_shift][i_y_shift][theta_i][3]-Eij[j_x_shift][j_y_shift][theta_j][3]);
   
   all_m_atoms_offset[central_atom.id].intra3++;
   // std::cout << i_x_shift << ", " << i_y_shift << ", " << j_x_shift << ", " << j_y_shift << ", " << 0.5*(J_i+J_j) << ", " << Dx*cos(twist_angle)-Dy*sin(twist_angle) << ", " <<  Dx*sin(twist_angle)+Dy*cos(twist_angle) << ", " << Dz << std::endl;
   return exchange;
}

std::array<double,4> match_inter_exchange(double dx, double dy, std::vector<std::vector<double> > &Eij){
   std::array<double,4> exchange({-20.0});
   double new_shift_error = 10.0;
   double old_shift_error = 10.0;
   int min_index = -1;

   // #pragma omp parallel for num_threads(2) reduction(min, min_index)
   for(int i = 0; i < Eij.size(); i++) {
      new_shift_error = sqrt((Eij[i][0] - dx)*(Eij[i][0] - dx) + (Eij[i][1] - dy)*(Eij[i][1] - dy));
      if(new_shift_error < old_shift_error) {
         old_shift_error = new_shift_error;
         min_index = i;
      }
   }
   
   exchange[0] = Eij.at(min_index)[2];
   exchange[1] = Eij[min_index][3];
   exchange[2] = Eij[min_index][4];
   exchange[3] = Eij[min_index][5];
   return exchange;
}

std::array<double,4> calculate_intra_Jani(spin &atom_i, spin &atom_j, double distance, double angle){

   std::array<double,4> exchange({0.0});
   if (distance <= intra_nn_dist_1) {
      exchange[0] = Jintra1_AB;
      exchange[3] = D_intra_z_constant;
      exchange[2] = (D_intra_x_constant*sin(angle)+D_intra_y_constant*cos(angle)); //-D_y
      exchange[1] = (D_intra_x_constant*cos(angle)-D_intra_y_constant*sin(angle)); //D_x
      all_m_atoms[atom_i.id].intra1++;
   } else if (distance <= intra_nn_dist_2) {
      if(atom_i.l_id == 5) exchange[0] = Jintra2_ABprime;
      else exchange[0] = Jintra2_AB;
      all_m_atoms[atom_i.id].intra2++;
   }
   return exchange;
}

std::array<double,4> calculate_inter_Jani(spin &atom_i, spin &atom_j, double distance, double angle) {
   std::array<double,4> exchange({0.0});

      if (distance <= inter_AB_dist_1) {
         exchange[0] = Jinter1_AB;
         all_m_atoms[atom_i.id].inter1++;
      }  else if  (distance <= inter_AB_dist_2 && (atom_i.l_id == 6 || atom_j.l_id == 6)) {
         exchange[0] = Jinter2_AB;
         all_m_atoms[atom_i.id].inter2++;
      }  else if (distance <= inter_AB_dist_3 && (atom_i.l_id == 5 || atom_j.l_id == 5)) {
         exchange[0] = Jinter2_AB_prime;
         all_m_atoms[atom_i.id].inter2++;
      }
      return exchange;
}