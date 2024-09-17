#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "positions.hpp"
#include "exchange.hpp"
#include "initialise.hpp"


void read_in_atoms(std::string filename, int n_atoms, std::vector <spin > &atom2){
   std::ifstream ifile(filename);
   std::string line;
   int temp;
   for (int i = 0; i < n_atoms; i ++){
      getline(ifile,line);
     std::stringstream liness(line.c_str());
     liness >> atom2[i].id >> atom2[i].x >> atom2[i].y >> atom2[i].z >> atom2[i].S >> atom2[i].l_id >> temp;
     //   outfile <<atom_id[i] << "\t" << x_in[i] << '\t' << y_in[i] << "\t" << z_in[i] << "\t" << S_in[i] << "\t" << temp << "\t" << temp << std::endl;

   }
   ifile.close();
}

void read_in_exchange(std::string filename, std::vector<std::vector<double>> & Jij){
   std::ifstream ifile2(filename);
   std::string line;

   //optional control
//    double x_step = 0.05;
//    double y_step = 0.0866;
    for(int i=0; i<10000; i++){
        getline(ifile2,line);
        std::stringstream liness(line.c_str());
        double ii;
        double ij;
        liness >> ii >> ij >> Jij[ii*2][ij*2];
        // std::cout << ii << ", " << ii*2 << std::endl;
   //   outfile <<atom_id[i] << "\t" << x_in[i] << '\t' << y_in[i] << "\t" << z_in[i] << "\t" << S_in[i] << "\t" << temp << "\t" << temp << std::endl;
   //   std::cout <<ii << "\t" << ij << "\t" <<  Jint[ii][ij] << std::endl;
 }

    for(int i = 0; i < 201; i+=2) {
        for(int j = 1; j < 200; j+=2) {
            Jij[i][j] = 0.5*(Jij[i][j-1]+Jij[i][j+1]);
        }
    }
    for(int j = 0; j < 201; j++) {
        for(int i = 1; i < 201; i+=2) {
            Jij[i][j] = 0.5*(Jij[i-1][j]+Jij[i+1][j]);
        }
    }
    ifile2.close();

}


void read_in_dmi(std::string filename, std::vector < std::vector < double > > &Dx, std::vector < std::vector < double > > &Dy, std::vector < std::vector < double > > &Dz){
   std::ifstream ifile2(filename);
   std::string line;

   for(int i=0; i<40000; i++){
     getline(ifile2,line);
     std::stringstream liness(line.c_str());
     double ii;
     double ij;
     liness >> ii >> ij >> Dx[ii][ij] >> Dy[ii][ij] >> Dz[ii][ij];
   //   outfile <<atom_id[i] << "\t" << x_in[i] << '\t' << y_in[i] << "\t" << z_in[i] << "\t" << S_in[i] << "\t" << temp << "\t" << temp << std::endl;
   // std::cout <<ii << "\t" << ij << "\t" <<  Dx[ii][ij]<< "\t" <<  Dy[ii][ij]<< "\t" <<  Dz[ii][ij] << std::endl;
 }
   ifile2.close();
}
