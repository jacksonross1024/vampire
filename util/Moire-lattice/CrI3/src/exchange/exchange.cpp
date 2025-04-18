#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "initialise.hpp"
#include "exchange.hpp"
#include <list>

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
double intra_nn_dist_3 = 8.01; //A

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
double J_constant = eVtoJ/1000.0; //1 meV
double J_intra_1=2.5*J_constant;
double J_intra_2=0.75*J_constant;
double J_intra_3=-0.01*J_constant;

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

double Dx_substrate = 0.025;
double Dy_substrate = 0.0;
uint64_t number_of_interactions = 0;

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
int config_cells_x = (number_of_unit_cells_x/2)+1;
int config_cells_y = (number_of_unit_cells_y/2)+1;

std::vector <double > crossProduct(std::vector <double >A, std::vector <double > B){
   std::vector <double > P(3,0.0);
    P[0] = A[1] * B[2] - A[2] * B[1];
    P[1] = A[2] * B[0] - A[0] * B[2];
    P[2] = A[0] * B[1] - A[1] * B[0];
    return P;
}

void print_interaction_header(){
   std::ofstream outfile3 ("header_interactions.ucf");
   std::string interaction_type;
   if(DMI) interaction_type = "tensorial";
   else interaction_type = "normalised-isotropic";
   outfile3 << number_of_interactions <<  "\t" << interaction_type << std::endl;//<<" 0 0 0 "<<J3S2_x << " 0 0 0 " << J3S2_z << std::endl;
   
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

   const double bsize = 1.2*range;

   Jinter1_AB = 0.1379*J_constant;
   Jinter2_AB = -0.10125*J_constant;
   Jinter3_AB = -0.035*J_constant;
   
   Dx_substrate = 0.01819374 * DMI_sub_scaling; // 0.00606458;//  0.009; // 0.01819374; 
   
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
   
   // // std::vector< int> interactions_list;
   // std::ofstream correlation_file;
   // correlation_file.open("moire-lattice-constants.txt");
   // std::vector<std::array<double, 3> > zero_correlation;
   // now calculate neighbour list looping over boxes
   vtimer_t timer;
      timer.start();
   //(J-J_inter_scaling*std::abs(J))*J_constant;
   // Jinter1_AB = Jinter1_AB - std::abs(Jinter1_AB)*J_inter_scaling*2;
   // Jinter2_AB = Jinter2_AB - std::abs(Jinter2_AB)*J_inter_scaling*2;
   // Jinter3_AB = Jinter3_AB - std::abs(Jinter3_AB)*J_inter_scaling*2;
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
                        for(int ai = 0; ai < shift_boxes[i][j][k].size(); ai++){
                           // atom_index++;
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

                              if(adx*adx <= (6.93*6.93) && ady*ady <= (6.003*6.003) && ((atom_i.l_id-2) == atom_j.l_id)){
                                 moire_shift[atom_i.id].push_back(atom_j);
                                 // if(adx == 0.0 && ady == 0.0) std::cout << atom_i.l_id << ", " << atom_j.l_id << ", " << x_i << ", " << y_i << ", " << x_j << ", " << y_j << std::endl;
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

   for(int i = 0; i < moire_shift.size(); i++) {
      double min_x = a0x;
      double min_y = a1y;
      double min_r = 100.0;
      spin ref_spin = all_m_atoms[i];
      if(ref_spin.S != 3) continue;
      std::list<spin> ref_list = moire_shift[i];
      int count = 0;

      for(auto shift = ref_list.begin(); shift != ref_list.end(); ++shift) {
         spin shift_spin = *shift;
         
         double dx = ref_spin.x - shift_spin.x ;
         double dy = ref_spin.y - shift_spin.y ;
         double dr = (dx*dx + dy*dy);
         if(dr < min_r) {
            min_r = dr;
            min_x = dx;
            min_y = dy;
            // if(dx == 0.0 && dy == 0.0) std::cout << i << ", " <<  shift_spin.id << ", " << shift_spin.x << ", " << ref_spin.x << std::endl;
         }
         count++;
      }
      // if(count == 0 ) continue;
      // if( count > 3 ) std::cout << i << ", " <<  count << ", " << min_x << ", " << min_y << std::endl;
      int dx_cell = ref_spin.unit_x_lr;
      int dy_cell = ref_spin.unit_y_lr;
      
      double change_x = min_x;
      double change_y = min_y;
      change_y = fabs(change_y)/a1y;
      

      change_x = (change_x - change_y*a1x)/a0x;
      if(change_x < 0.0) change_x += 1.0;
      else if (change_x >= 1.0) change_x -= 1.0;

      int dx = int(round(change_x*99.0));
      int dy = int(round(change_y*99.0));
       if(dx > 99 || dx < 0 || dy > 99 || dy < 0) {
         std::cerr << "shift problem: (" <<  dx << ", " << dy  << ", " << change_x << ", " << change_y << ", " << min_x << ", " << min_y << std::endl;
            exit(1);
      }
      
      unit_cell_shifts.at(dx_cell).at(dy_cell)[0] += 1;
      unit_cell_shifts[dx_cell][dy_cell][1] += dx;
      unit_cell_shifts[dx_cell][dy_cell][2] += dy;
   }

      char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
         std::cerr << "Fatal getcwd error in datalog." << std::endl;
      }
      std::ofstream shift_file(std::string(directory) + "/shifted_constants.txt");
      for(int i = 0; i < unit_cell_shifts.size(); i++){
      for (int j = 0; j < unit_cell_shifts[i].size(); j++) {
         int occupancy = unit_cell_shifts[i][j][0];
         if(occupancy > 0) {
            unit_cell_shifts[i][j][1] = round(unit_cell_shifts[i][j][1]/occupancy);
            unit_cell_shifts[i][j][2] = round(unit_cell_shifts[i][j][2]/occupancy);
         } else {
            unit_cell_shifts[i][j][1] = 66;
         }
        
         int i_shift = unit_cell_shifts[i][j][1];
         int j_shift = unit_cell_shifts[i][j][2];
         if(i_shift > 99 || j_shift > 99) std::cout << "problems " << i << ", " << j << ", " << occupancy << ", " << i_shift << ", " << j_shift << std::endl;
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
   // exit(1);

      std::ofstream lattice_file(std::string(directory) + "/moire_lattice_vectors.txt");
   #pragma omp parallel num_threads(16) reduction(+:number_of_interactions) 
   {
      #pragma omp single 
      std::cout << "preparing Moire exchange with " << omp_get_num_threads() << " omp threads" << std::endl;

      std::vector<std::vector<std::vector<double> > > local_config_energy;
      local_config_energy.resize(number_of_unit_cells_x/2 + 1);
      for(int i = 0; i < number_of_unit_cells_x/2 + 1; i++) {
     
         local_config_energy[i].resize(number_of_unit_cells_y/2 + 1);
         for(int j = 0; j < number_of_unit_cells_y/2 + 1; j++) {
         
            local_config_energy[i][j].resize(4*16,0.0);
         }
      }
   std::stringstream otext;
   otext.precision(6);

   std::stringstream lattice_text;
   for(int i=0; i<xb; i++){

       #pragma omp single nowait
      if(i%int(xb*0.1) == 0) std::cout << "." << std::flush;
      // if(number_of_interactions > 0) std::cout << 100*number_of_interactions/interaction_estimate << "%..." << std::flush;
      #pragma omp for schedule(dynamic)
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
                           if(atom_i.S == 5) continue;
                           const double x_i = atom_i.x;
                           const double y_i = atom_i.y;
                           const double z_i = atom_i.z;
                           // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 0] += 1;
                           // if(atom_i.unit_x_lr == 20 && atom_i.unit_y_lr == 20) std::cout << atom_i.id << ", " << x_i << ", " << y_i << std::endl;
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
                              if(dL2 < r2 ){
                                 // if(std::abs(adx) < 1e-2 && std::abs(ady) < 1e-2) {
                                 //    //correlation_file << std::abs(adx) << ", " << std::abs(ady) << ", " << atom_i.id << ", " << atom_i.x << ", " << atom_i.y << ", " << atom_j.id << ", " << atom_j.x << ", " << atom_j.y << std::endl;
                                 //    zero_correlation.push_back({atom_i.x, atom_i.y, double(atom_i.l_id)});
                                 // }
                                 // std::cout << dL2 << ", " << r2 << ", " << x_i << ", " << y_i << ", " << z_i << ", " << x_j << ", " << y_j << ", " << z_j << std::endl;
                                 double angle_i = atan2(ady,adx);// - twist_angle;// - M_PI*0.5;
                                 double angle_j = atan2(-ady,-adx);// - twist_angle;
                                 std::array<double, 4> exchange({0.0,0.0,0.0,0.0});
                                 if(atom_i.S == atom_j.S) {
                                 
                                    if(atom_i.l_id == 1) {
                                       // angle_i += twist_angle;
                                       // angle_j += twist_angle;
                                       if(dL2 < intra_nn_dist_1) {
                                          exchange = match_intra1_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr1_1NN );
                                          // all_m_atoms[atom_i.id].intra1_count++;
                                          // all_m_atoms[atom_i.id].J_intra1 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_intra1 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_intra1 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_intra1 += exchange[3]/J_constant;
                                          
                                          // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 1] += exchange[0];
                                       } else if (dL2 < intra_nn_dist_2) {
                                          exchange = match_intra2_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr1_2NN );
                                          // all_m_atoms[atom_i.id].intra2_count++;
                                          // all_m_atoms[atom_i.id].J_intra2 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_intra2 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_intra2 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_intra2 += exchange[3]/J_constant;
                                          
                                          // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 2] += exchange[0];
                                       } else if (dL2 < intra_nn_dist_3) {
                                          exchange = match_intra3_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr1_3NN );
                                          // all_m_atoms[atom_i.id].intra3_count++;
                                          // all_m_atoms[atom_i.id].J_intra3 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_intra3 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_intra3 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_intra3 += exchange[3]/J_constant;
                                          
                                          // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 3] += exchange[0];
                                       } else continue;   
                                       
                                       
                                    } if(atom_i.l_id == 2) {
                                       // angle_i += twist_angle;
                                       // angle_j += twist_angle;
                                       if(dL2 < intra_nn_dist_1) {
                                          exchange = match_intra1_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr2_1NN );
                                          // all_m_atoms[atom_i.id].intra1_count++;
                                          // all_m_atoms[atom_i.id].J_intra1 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_intra1 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_intra1 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_intra1 += exchange[3]/J_constant;
                                          lattice_text << x_i << "\t" << y_i << "\t" << x_j << "\t" << y_j << "\n";
                                          // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 1] += exchange[0];
                                       } else if (dL2 < intra_nn_dist_2) {
                                          exchange = match_intra2_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr2_2NN );
                                          // all_m_atoms[atom_i.id].intra2_count++;
                                          // all_m_atoms[atom_i.id].J_intra2 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_intra2 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_intra2 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_intra2 += exchange[3]/J_constant;
                                          
                                          // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 2] += exchange[0];
                                       } else if (dL2 < intra_nn_dist_3) {
                                          exchange = match_intra3_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr2_3NN );
                                          // all_m_atoms[atom_i.id].intra3_count++;
                                          // all_m_atoms[atom_i.id].J_intra3 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_intra3 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_intra3 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_intra3 += exchange[3]/J_constant;

                                          // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 3] += exchange[0];
                                       } else continue;                              
                                    } if(atom_i.l_id == 3) {
                                       angle_i -= twist_angle*0.99;
                                       angle_j -= twist_angle*0.99;
                                       if(dL2 < intra_nn_dist_1) {
                                          exchange = match_intra1_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr3_1NN );
                                          lattice_text << x_i << "\t" << y_i << "\t" << x_j << "\t" << y_j << "\n";
                                          // all_m_atoms[atom_i.id].intra1_count++;
                                          // all_m_atoms[atom_i.id].J_intra1 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_intra1 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_intra1 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_intra1 += exchange[3]/J_constant;

                                          // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 1] += exchange[0];
                                       } else if (dL2 < intra_nn_dist_2) {
                                          exchange = match_intra2_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr3_2NN );
                                          

                                         // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 2] += exchange[0];
                                       } else if (dL2 < intra_nn_dist_3) {
                                          exchange = match_intra3_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr3_3NN );
                                          // all_m_atoms[atom_i.id].intra3_count++;
                                          // all_m_atoms[atom_i.id].J_intra3 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_intra3 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_intra3 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_intra3 += exchange[3]/J_constant;

                                          // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 3] += exchange[0];
                                       } else continue;                              
                                    }if(atom_i.l_id == 4) {
                                       angle_i -= twist_angle*0.99;
                                       angle_j -= twist_angle*0.99;
                                       if(dL2 < intra_nn_dist_1) {
                                          exchange = match_intra1_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr4_1NN );
                                          // all_m_atoms[atom_i.id].intra1_count++;
                                          // all_m_atoms[atom_i.id].J_intra1 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_intra1 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_intra1 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_intra1 += exchange[3]/J_constant;

                                          // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 1] += exchange[0];

                                       } else if (dL2 < intra_nn_dist_2) {
                                          exchange = match_intra2_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr4_2NN );
                                          // all_m_atoms[atom_i.id].intra2_count++;
                                          // all_m_atoms[atom_i.id].J_intra2 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_intra2 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_intra2 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_intra2 += exchange[3]/J_constant;

                                          // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 2] += exchange[0];
                                       } else if (dL2 < intra_nn_dist_3) {
                                          exchange = match_intra3_exchange(angle_i, angle_j, atom_i, atom_j, Eintra_Cr4_3NN );
                                          // all_m_atoms[atom_i.id].intra3_count++;
                                          // all_m_atoms[atom_i.id].J_intra3 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_intra3 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_intra3 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_intra3 += exchange[3]/J_constant;

                                          // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 3] += exchange[0];
                                       } else continue;                          
                                    }
                               
                                 } else {
                                    if(atom_i.h_id == atom_j.h_id) {
                                       if(atom_i.l_id == 1) {exchange = match_inter_exchange(atom_i.id, atom_j.id, adx, ady, dL2, Einter_Cr1);}
                                       
                                       else if(atom_i.l_id == 2) {exchange = match_inter_exchange(atom_i.id, atom_j.id, adx, ady,  dL2,Einter_Cr2);}
                                          
                                       else if(atom_i.l_id == 3) {exchange = match_inter_exchange(atom_i.id, atom_j.id, adx, ady,  dL2,Einter_Cr3);}
                                          
                                       else if(atom_i.l_id == 4) {exchange = match_inter_exchange(atom_i.id, atom_j.id, adx, ady, dL2, Einter_Cr4);}
                                       else continue;
                                       exchange[0] *= J_twist_reduction;

                                    } else {
                                       if(dL2 <= inter_AB_dist_1) {
                                          exchange[0] = Jinter1_AB;
                                          exchange[0] *= J_prist_reduction;
                                          // all_m_atoms[atom_i.id].inter1_count++;
                                          // all_m_atoms[atom_i.id].J_inter1 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_inter1 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_inter1 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_inter1 += exchange[3]/J_constant;
                                       } else if (dL2 <= inter_AB_dist_2) {
                                          exchange[0] = Jinter2_AB;
                                          exchange[0] *= J_prist_reduction;
                                          // all_m_atoms[atom_i.id].inter2_count++;
                                          // all_m_atoms[atom_i.id].J_inter2 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_inter2 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_inter2 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_inter2 += exchange[3]/J_constant;
                                       } else {
                                          exchange[0] = Jinter3_AB;
                                          exchange[0] *= J_prist_reduction;
                                          // all_m_atoms[atom_i.id].inter3_count++;
                                          // all_m_atoms[atom_i.id].J_inter3 += exchange[0]/J_constant;
                                          // all_m_atoms[atom_i.id].Dx_inter3 += exchange[1]/J_constant;
                                          // all_m_atoms[atom_i.id].Dy_inter3 += exchange[2]/J_constant;
                                          // all_m_atoms[atom_i.id].Dz_inter3 += exchange[3]/J_constant;
                                       }
                                       
                                    }
                                    // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 1] += 1;
                                    // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 3] += exchange[0];
                                    // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 4] += exchange[1];
                                    // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 5] += exchange[2];
                                    // local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 6] += exchange[3];
                                 }
                                                              
                                 // if(exchange[0] == -60) continue;

                              // ##pragma omp critical 
                                 // {
                                 spin_config_energy(atom_i, dL2, atom_j, exchange, local_config_energy);

                                   if(DMI) {  otext << number_of_interactions <<  "\t" << atom_i.id << '\t' << atom_j.id << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' <<\
                                                //xx                     xy-> Dz                 xz -> -Dy
                                                  exchange[0] << "\t" << exchange[3] << "\t" << -exchange[2] << "\t" << \
                                                //yx -> -Dz              yy                      yz -> Dx
                                                 -exchange[3] << "\t" << exchange[0] << "\t" <<  exchange[1] << "\t" << \
                                                //zx -> Dy               yz -> -Dx               zz
                                                  exchange[2] << "\t" <<-exchange[1] << "\t" <<  exchange[0]*1.0 << "\n"; }
                                 else {   otext << number_of_interactions <<  "\t" << atom_i.id << '\t' << atom_j.id << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' <<\
                              //xx                     xy-> Dz                 xz -> -Dy
                                 exchange[0] << "\n";// << 0.0 << "\t" << 0.0 << "\t" << \
                              //yx -> -Dz              yy                      yz -> Dx
                             //    0.0 << "\t" << exchange[0] << "\t" <<  0.0 << "\t" << \
                              //zx -> Dy               yz -> -Dx               zz
                              //   0.0 << "\t" << 0.0 << "\t" <<  exchange[0] << "\n"; 
                                 }
                                    // config_energy.at(atom_i.unit_x_lr).at(atom_i.unit_y_lr).at((atom_i.S-1)*5+0) += 1.0;
                                    // config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*5+1] += exchange[0]/J_constant;
                                    // config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*5+2] += exchange[1]/J_constant;
                                    // config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*5+3] += exchange[2]/J_constant;
                                    // config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr].at((atom_i.S-1)*5+4) += exchange[3]/J_constant;

                                    
                                    number_of_interactions++;
                                 // }                           
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
      #pragma omp critical
      {
         lattice_file << lattice_text.str();
         outfile4 << otext.str();
         std::cout << omp_get_thread_num() << " calculated " << number_of_interactions << " of interactions" << std::endl;

         //double J_const_inv = 1.0/J_constant;
         
         for(int i = 0; i < global_config_energy.size(); i++) {
            for(int j = 0; j < global_config_energy[i].size(); j++) {
               for(int k =0; k < global_config_energy[i][j].size(); k++){
                     
                  global_config_energy[i][j].at(k) += local_config_energy[i][j].at(k); 
                  //else global_config_energy[i][j].at(k) += local_config_energy[i][j].at(k)*J_const_inv;
                  
               }
            }
         }
      }
      
   }
   lattice_file << std::flush;
   outfile4 << std::flush;
   outfile4.close();

   timer.stop();
   // std::cout << "done!  << std::endl;
   std::cout << number_of_interactions << " [completed] [" << timer.elapsed_time() << " s]" << std::endl;
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

std::array<double,4> match_intra1_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom, std::vector<std::vector< double > > &Eij){
   std::array<double,4> exchange({0.0});

   int i_x_shift = 66;
   int i_y_shift = 0;

   if(central_atom.h_id == 1) {
      i_x_shift = (unit_cell_shifts[central_atom.unit_x][central_atom.unit_y][1]);
      i_y_shift = (unit_cell_shifts[central_atom.unit_x][central_atom.unit_y][2]);

      // j_x_shift = (unit_cell_shifts[j_atom.unit_x][j_atom.unit_y][1]);
      // j_y_shift = (unit_cell_shifts[j_atom.unit_x][j_atom.unit_y][2]);
   }
   

   double theta =  std::abs(angle_i)*180.0/M_PI;
   int theta_i = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0)== 1.0) ? (1) : ( round(theta/30.0) == 3.0 ? (0):-1) );
   theta =  std::abs(angle_j)*180.0/M_PI;
   int theta_j = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0)== 1.0) ? (1) : (round(theta/30.0) == 3.0 ? (0):-1) );

    if(theta_i == -1 || theta_j == -1) {
      std::cout << "\n " << round(std::abs(angle_i)*180.0/M_PI/30.0) << ", " << round(std::abs(angle_j)*180.0/M_PI/30.0) << ", " << angle_i << ", " << angle_j << ", " << theta_i << ", " << theta_j << ", " << j_atom.Gx << ", " << j_atom.Gy <<  std::endl;
     exit(1);
    //  return exchange;
   }
      // double DMI_sub_vector_x = 0.0;
      // double DMI_sub_vector_y = 0.0;
      // double DMI_sub_vector_z = 1.0;
      double ex_vec_x = cos(angle_i);
      double ex_vec_y = sin(angle_i);
      double ex_vec_z = 0.0;
   if(central_atom.S == 1) {
      exchange[0] = J_intra_reduction*0.5*(Eij.at(i_x_shift*100 +i_y_shift).at(theta_i*4 +0))*2;//+Eij.at(j_x_shift*100 + j_y_shift).at(theta_j*4 + 0));

      double DMI = Dx_substrate*exchange[0];
      exchange[1] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +1])*2 + DMI*(DMI_sub_vector_y*ex_vec_z - ex_vec_y*(-DMI_sub_vector_z));//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 1]);
      exchange[2] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +2])*2 + DMI*((-DMI_sub_vector_z)*ex_vec_x - ex_vec_z*DMI_sub_vector_x);
      exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2 + DMI*(DMI_sub_vector_x*ex_vec_y - ex_vec_x*DMI_sub_vector_y);

    //  std::cout << DMI << ", " << ex_vec_x << ", " << ex_vec_y << ", " << (DMI_sub_vector_x*ex_vec_y - ex_vec_x*DMI_sub_vector_y) << std::endl;
   //   std::cout << 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +1])*2/J_constant << " + " << Dx_substrate*exchange[0]/J_constant << " = " << exchange[1]/J_constant << std::endl;
   } else if(central_atom.S == 2) {
      exchange[0] = J_intra_reduction*0.5*(Eij.at(i_x_shift*100 +i_y_shift).at(theta_i*4 +0))*2;//+Eij.at(j_x_shift*100 + j_y_shift).at(theta_j*4 + 0));
      exchange[1] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +1])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 1]);
      exchange[2] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +2])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 2]);
      exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 3]);
   } else if(central_atom.S == 3) {
      exchange[0] = J_intra_reduction*0.5*(Eij.at(i_x_shift*100 +i_y_shift).at(theta_i*4 +0))*2;//+Eij.at(j_x_shift*100 + j_y_shift).at(theta_j*4 + 0));
      double Dx = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +1])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 1]);
      double Dy = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +2])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 2]);
      exchange[1] = Dx*cos(twist_angle) - Dy*sin(twist_angle);
      exchange[2] = Dx*sin(twist_angle) + Dy*cos(twist_angle);
      exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 3]);
   } else if(central_atom.S == 4) {
      exchange[0] = J_intra_reduction*0.5*(Eij.at(i_x_shift*100 +i_y_shift).at(theta_i*4 +0))*2;

      double DMI = Dx_substrate*exchange[0];
      double Dx = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +1])*2 + DMI*(DMI_sub_vector_y*ex_vec_z - ex_vec_y*DMI_sub_vector_z);//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 1]);
      double Dy = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +2])*2 + DMI*(DMI_sub_vector_z*ex_vec_x - ex_vec_z*DMI_sub_vector_x);//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 2]);
      exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2 + DMI*(DMI_sub_vector_x*ex_vec_y - ex_vec_x*DMI_sub_vector_y);

      exchange[1] = Dx*cos(twist_angle) - Dy*sin(twist_angle);
      exchange[2] = Dx*sin(twist_angle) + Dy*cos(twist_angle);
    
   }

   return exchange;
}

std::array<double,4> match_intra2_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom, std::vector<std::vector<  double > > &Eij){
   std::array<double,4> exchange({0.0});

      int i_x_shift = 66;
   int i_y_shift = 0;
   
   if(central_atom.h_id == 1 ) {
      i_x_shift = (unit_cell_shifts[central_atom.unit_x][central_atom.unit_y][1]);
      i_y_shift = (unit_cell_shifts[central_atom.unit_x][central_atom.unit_y][2]);

      // j_x_shift = (unit_cell_shifts[j_atom.unit_x][j_atom.unit_y][1]);
      // j_y_shift = (unit_cell_shifts[j_atom.unit_x][j_atom.unit_y][2]);
   } 
   
   // double theta =  std::abs(angle)*180.0/M_PI;
   int theta_i = int(round(((angle_i < 0.0) ? (angle_i+=2*M_PI) : (angle_i)) *179.0/M_PI/60.0));
   // theta =  std::abs(angle-180.0)*180.0/M_PI;
   int theta_j = int(round(((angle_j < 0.0) ? (angle_j+=2*M_PI) : (angle_j)) *179.0/M_PI/60.0));

   if( theta_i > 5 || theta_j > 5 || theta_i < 0 || theta_j < 0 ) {
      std::cout << "problem: " << theta_i << " , " << theta_j << ", " << angle_i << ", " << angle_j << std::endl;
      exit(1);
      return exchange;
   }
      // double DMI_sub_vector_x = 0.0;
      // double DMI_sub_vector_y = 0.0;
      // double DMI_sub_vector_z = 1.0;
      double ex_vec_x = cos(angle_i);
      double ex_vec_y = sin(angle_i);
      double ex_vec_z = 0.0;

   if(central_atom.S == 1) {
      exchange[0] = J_intra_reduction*0.5*(Eij.at(i_x_shift*100 +i_y_shift).at(theta_i*4 +0))*2;//+Eij.at(j_x_shift*100 + j_y_shift).at(theta_j*4 + 0));

      double DMI = Dx_substrate*exchange[0];

      exchange[1] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +1])*2 + DMI*(DMI_sub_vector_y*ex_vec_z - ex_vec_y*(-DMI_sub_vector_z));//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 1]);
      exchange[2] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +2])*2 + DMI*((-DMI_sub_vector_z)*ex_vec_x - ex_vec_z*DMI_sub_vector_x);
      exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2 + DMI*(DMI_sub_vector_x*ex_vec_y - ex_vec_x*DMI_sub_vector_y);
   } else if(central_atom.S == 2) {
      exchange[0] = J_intra_reduction*0.5*(Eij.at(i_x_shift*100 +i_y_shift).at(theta_i*4 +0))*2;//+Eij.at(j_x_shift*100 + j_y_shift).at(theta_j*4 + 0));
      exchange[1] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +1])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 1]);
      exchange[2] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +2])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 2]);
      exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 3]);
   } else if(central_atom.S == 3) {
      exchange[0] = J_intra_reduction*0.5*(Eij.at(i_x_shift*100 +i_y_shift).at(theta_i*4 +0))*2;//+Eij.at(j_x_shift*100 + j_y_shift).at(theta_j*4 + 0));
      double Dx = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +1])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 1]);
      double Dy = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +2])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 2]);
      exchange[1] = Dx*cos(twist_angle) - Dy*sin(twist_angle);
      exchange[2] = Dx*sin(twist_angle) + Dy*cos(twist_angle);
      exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 3]);
   } else if(central_atom.S == 4) {
      exchange[0] = J_intra_reduction*0.5*(Eij.at(i_x_shift*100 +i_y_shift).at(theta_i*4 +0))*2;

      double DMI = Dx_substrate*exchange[0];
      double Dx = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +1])*2 + DMI*(DMI_sub_vector_y*ex_vec_z - ex_vec_y*DMI_sub_vector_z);//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 1]);
      double Dy = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +2])*2 + DMI*(DMI_sub_vector_z*ex_vec_x - ex_vec_z*DMI_sub_vector_x);//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 2]);
      exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2 + DMI*(DMI_sub_vector_x*ex_vec_y - ex_vec_x*DMI_sub_vector_y);

      exchange[1] = Dx*cos(twist_angle) - Dy*sin(twist_angle);
      exchange[2] = Dx*sin(twist_angle) + Dy*cos(twist_angle);
      //exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 3]);
   }


   // exchange[0] *= J_intra_reduction;
   return exchange;
}

std::array<double,4> match_intra3_exchange(double angle_i, double angle_j, spin &central_atom, spin &j_atom, std::vector<std::vector< double > > &Eij){
   std::array<double,4> exchange({0.0});

      int i_x_shift = 66;
   int i_y_shift = 0;
   
   if(central_atom.h_id == 1 ) {
      i_x_shift = (unit_cell_shifts[central_atom.unit_x][central_atom.unit_y][1]);
      i_y_shift = (unit_cell_shifts[central_atom.unit_x][central_atom.unit_y][2]);

      // j_x_shift = (unit_cell_shifts[j_atom.unit_x][j_atom.unit_y][1]);
      // j_y_shift = (unit_cell_shifts[j_atom.unit_x][j_atom.unit_y][2]);
   } 
   

   double theta =  std::abs(angle_i)*180.0/M_PI;
   int theta_i = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0)== 1.0) ? (1) : ( round(theta/30.0) == 3.0 ? (0):-1) );
   theta =  std::abs(angle_j)*180.0/M_PI;
   int theta_j = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0)== 1.0) ? (1) : (round(theta/30.0) == 3.0 ? (0):-1) );
   if(theta_i == -1 || theta_j == -1) {
      std::cout << "\n " << round(std::abs(angle_i)*180.0/M_PI/30.0) << ", " << round(std::abs(angle_j)*180.0/M_PI/30.0) << ", " << angle_i << ", " << angle_j << ", " << central_atom.Gx << ", " << central_atom.Gy << ", " << j_atom.Gx << ", " << j_atom.Gy <<  std::endl;
      exit(1);
      return exchange;
   }
      // double DMI_sub_vector_x = 0.0;
      // double DMI_sub_vector_y = 0.0;
      // double DMI_sub_vector_z = 1.0;
      double ex_vec_x = cos(angle_i);
      double ex_vec_y = sin(angle_i);
      double ex_vec_z = 0.0;
   if(central_atom.S == 1) {
      exchange[0] = J_intra_reduction*0.5*(Eij.at(i_x_shift*100 +i_y_shift).at(theta_i*4 +0))*2;//+Eij.at(j_x_shift*100 + j_y_shift).at(theta_j*4 + 0));

      double DMI = Dx_substrate*exchange[0];
      exchange[1] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +1])*2 + DMI*(DMI_sub_vector_y*ex_vec_z - ex_vec_y*(-DMI_sub_vector_z));//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 1]);
      exchange[2] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +2])*2 + DMI*((-DMI_sub_vector_z)*ex_vec_x - ex_vec_z*DMI_sub_vector_x);
      exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2 + DMI*(DMI_sub_vector_x*ex_vec_y - ex_vec_x*DMI_sub_vector_y);
   } else if(central_atom.S == 2) {
      exchange[0] = J_intra_reduction*0.5*(Eij.at(i_x_shift*100 +i_y_shift).at(theta_i*4 +0))*2;//+Eij.at(j_x_shift*100 + j_y_shift).at(theta_j*4 + 0));
      exchange[1] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +1])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 1]);
      exchange[2] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +2])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 2]);
      exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 3]);
   } else if(central_atom.S == 3) {
      exchange[0] = J_intra_reduction*0.5*(Eij.at(i_x_shift*100 +i_y_shift).at(theta_i*4 +0))*2;//+Eij.at(j_x_shift*100 + j_y_shift).at(theta_j*4 + 0));
      double Dx = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +1])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 1]);
      double Dy = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +2])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 2]);
      exchange[1] = Dx*cos(twist_angle) - Dy*sin(twist_angle);
      exchange[2] = Dx*sin(twist_angle) + Dy*cos(twist_angle);
      exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 3]);
   } else if(central_atom.S == 4) {
      exchange[0] = J_intra_reduction*0.5*(Eij.at(i_x_shift*100 +i_y_shift).at(theta_i*4 +0))*2;

      double DMI = Dx_substrate*exchange[0];
      double Dx = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +1])*2 + DMI*(DMI_sub_vector_y*ex_vec_z - ex_vec_y*DMI_sub_vector_z);//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 1]);
      double Dy = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +2])*2 + DMI*(DMI_sub_vector_z*ex_vec_x - ex_vec_z*DMI_sub_vector_x);//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 2]);
      exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2 + DMI*(DMI_sub_vector_x*ex_vec_y - ex_vec_x*DMI_sub_vector_y);

      exchange[1] = Dx*cos(twist_angle) - Dy*sin(twist_angle);
      exchange[2] = Dx*sin(twist_angle) + Dy*cos(twist_angle);
      //exchange[3] = 0.5*(Eij[i_x_shift*100 +i_y_shift][theta_i*4 +3])*2;//+Eij[j_x_shift*100 + j_y_shift][theta_j*4 + 3]);
   }
   // exchange[0] *= J_intra_reduction;
   // all_m_atoms[central_atom.id].intra3++;
   // std::cout << i_x_shift << ", " << i_y_shift << ", " << j_x_shift << ", " << j_y_shift << ", " << 0.5*(J_i+J_j) << ", " << Dx*cos(twist_angle)-Dy*sin(twist_angle) << ", " <<  Dx*sin(twist_angle)+Dy*cos(twist_angle) << ", " << Dz << std::endl;
   return exchange;
}

std::array<double,4> match_inter_exchange(int atom_id, int nn_id, double dx, double dy, double dr, std::vector<std::vector<double> > &Eij){
   std::array<double,4> exchange({0.0});
   double new_shift_error = 10000.0;
   double old_shift_error = 10000.0;
   int min_index = -1;

   // #pragma omp parallel for num_threads(2) reduction(min, min_index)
   for(int i = 0; i < Eij.size(); i++) {
      new_shift_error = ((Eij[i][0] - dx)*(Eij[i][0] - dx) + (Eij[i][1] - dy)*(Eij[i][1] - dy));
      if(new_shift_error < old_shift_error) {
         old_shift_error = new_shift_error;
         min_index = i;
      }
   }
   if(min_index < 0 || min_index > 10000) std::cout << min_index << std::endl;
   exchange[0] = Eij.at(min_index)[2];
   exchange[1] = Eij[min_index][3];
   exchange[2] = Eij[min_index][4];
   exchange[3] = Eij[min_index][5];
 
   return exchange;
}

void spin_config_energy(spin & atom_i, double dr2, spin & atom_j, std::array<double, 4> &exchange, std::vector<std::vector<std::vector<double> > > & local_config_energy) {

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
         local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 12] += exchange[0];
         local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 13] += fabs(exchange[1]);
         local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 14] += fabs(exchange[2]);
         local_config_energy[atom_i.unit_x_lr][atom_i.unit_y_lr][(atom_i.S-1)*16 + 15] += fabs(exchange[3]);
      }    
   }
}
