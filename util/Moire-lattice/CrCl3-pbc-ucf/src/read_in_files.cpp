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
   if(!ifile.is_open()) {std::cerr << filename << " is not open" << std::endl; exit(1);}
   std::string line;
   int temp;
   for (int i = 0; i < n_atoms; i ++){
      getline(ifile,line);
     std::stringstream liness(line.c_str());
     liness >> atom2[i].id >> atom2[i].x >> atom2[i].y >> atom2[i].z >> atom2[i].S >> atom2[i].l_id >> atom2[i].h_id;
     //   outfile <<atom_id[i] << "\t" << x_in[i] << '\t' << y_in[i] << "\t" << z_in[i] << "\t" << S_in[i] << "\t" << temp << "\t" << temp << std::endl;

   }
   ifile.close();
}


void read_in_inter_exchanges(std::string J, std::string Dx, std::string Dy, std::string Dz,  std::vector<std::vector<double> > &Eij ) {
    std::ifstream J_file(J);
    std::string Jline;
    if(!J_file.is_open()) {std::cerr << J << " is not open" << std::endl; exit(1);}
    
    std::ifstream Dx_file(Dx);
    std::string Dxline;
    if(!Dx_file.is_open()) {std::cerr << Dx << " is not open" << std::endl; exit(1);}

    std::ifstream Dy_file(Dy);
    std::string Dyline;
    if(!Dy_file.is_open()) {std::cerr << Dy << " is not open" << std::endl; exit(1);}

    std::ifstream Dz_file(Dz);
    std::string Dzline;
    if(!Dz_file.is_open()) {std::cerr << Dz << " is not open" << std::endl; exit(1);}

    const double a_0 = 7.276;
    const double a_1 = 7.402;
    int i = 0;
    int total = 0;
    while(i < 100){
        int j  = 0;
        while( j < 100) {
            getline(J_file, Jline); 
            getline(Dx_file, Dxline); 
            getline(Dy_file, Dyline); 
            getline(Dz_file, Dzline);
            std::stringstream Jliness(Jline.c_str());
            std::stringstream Dxliness(Dxline.c_str());
            std::stringstream Dyliness(Dyline.c_str());
            std::stringstream Dzliness(Dzline.c_str());
            double J;
            double Dx;
            double Dy;
            double Dz;
            Jliness >> J; 
            Dxliness >> Dx; 
            Dyliness >> Dy; 
            Dzliness >> Dz;
            Eij.at(i*100+j)[0] = j*0.02*a_0 - a_0; //x pos
            Eij[i*100+j][1] = (99-i)*0.02*a_1 - a_1;
            //if(std::abs(Dz) > 0.2) std::cout << i << ", " << j << ", "  << j*0.01*a_0 - a_0 << ", " << (199-i)*0.01*a_1 - a_1 << ", " <<  Dz << std::endl;
            Eij[i*100+j][2] = (J > 0.0) ? (J-J_inter_scaling)*J_constant*1.0 : (J-J_inter_scaling)*J_constant;
            Eij[i*100+j][3] = DMI_inter_scaling*Dx*J_constant;
            Eij[i*100+j][4] = DMI_inter_scaling*Dy*J_constant;
            Eij[i*100+j][5] = DMI_inter_scaling*Dz*J_constant;
            j++;
            total++;
        }
        i++;
    }
    std::cout << total << std::endl;
    
    std::ofstream inter_out("Inter_exchange_out.txt");
    for(i = 0; i < Eij.size(); i++) {
        inter_out << Eij[i][0] << ", " << Eij[i][1] << ", " << Eij[i][2] << ", " << Eij[i][3] << ", " << Eij[i][4] << ", " << Eij[i][5] << "\n";
    }
    inter_out.close();
    
}

void read_in_intra_exchanges(std::string filename, std::vector<std::vector<double > > &Eij_1NN, \
                                                   std::vector<std::vector<double > > &Eij_2NN, \
                                                   std::vector<std::vector<double > > &Eij_3NN ) {

    std::ifstream ifile1(filename);
    std::string line;
    if(!ifile1.is_open()) {std::cerr  << " is not open" << std::endl; exit(1);}
    
    int count = 0;
    for(int k = 0; k < 3; k++) {
        
        for(int i = 0; i < 100; i++) {
           
           for(int j = 0; j < 100; j++) {
        
                getline(ifile1,line);
                std::stringstream liness(line.c_str());
                
                double dx;
                double dy;
                
                double J;
                double Dx;
                double Dy;
                double Dz;

                double theta;
                liness >> dx >> dy >> J >> Dx >> Dy >> Dz;

                // Honeycomb 1NN: |θ| ∈ {30,90,150} → codes {1,0,2} (must match classify_honeycomb_nn_angle)
                theta = std::abs(atan2(dy,dx))*180.0/M_PI;
                int theta_i = -1;
                {
                   const double cans[3] = {30.0, 90.0, 150.0};
                   const int codes[3] = {1, 0, 2};
                   double best_err = 1e9;
                   for(int t = 0; t < 3; ++t){
                      const double err = std::fabs(theta - cans[t]);
                      if(err < best_err){ best_err = err; theta_i = codes[t]; }
                   }
                   if(best_err > 1.0) theta_i = -1;
                }
                if(theta_i < 0) std::cout << k << ", " << i << ", " << j << ", " << theta << ", " << dx << ", " << dy << std::endl;
                Eij_1NN.at(i*100 + j).at(theta_i*4 + 0) = J*J_constant;
                Eij_1NN[i*100 + j][theta_i*4 + 1] = DMI_inter_scaling*Dx*J_constant;
                Eij_1NN[i*100 + j][theta_i*4 + 2] = DMI_inter_scaling*Dy*J_constant;
                Eij_1NN[i*100 + j][theta_i*4 + 3] = DMI_inter_scaling*Dz*J_constant;
            //  std::cout  <<  Eintra_Cr1_1NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr1_1NN.at(int_x).at(int_y).at(theta_i*4 + 1) << ", " << Eintra_Cr1_1NN.at(int_x).at(int_y).at(theta_i*4 + 2) << ", " << Eintra_Cr1_1NN.at(int_x).at(int_y).at(theta_i*4 + 3) << std::endl;
                count++;
            }
       
        }
    }
 std::cout << "1NN by thread " << omp_get_thread_num() << " count: " << count << std::endl;
    
    for(int k = 0; k < 6; k++) {
       
        for(int i = 0; i < 100; i++) {
           
           for(int j = 0; j < 100; j++) {
            getline(ifile1,line);
            std::stringstream liness(line.c_str());
                double dx;
                double dy;
                
                double J;
                double Dx;
                double Dy;
                double Dz;
                
                double theta;
                liness >> dx >> dy >> J >> Dx >> Dy >> Dz;
            // Honeycomb 2NN: θ ∈ {0,60,...,300} → index 0..5 (must match classify_honeycomb_2nn_angle)
            theta = atan2(dy,dx);
            if(theta < 0.0) theta += 2.0 * M_PI;
            if(theta >= 2.0 * M_PI) theta -= 2.0 * M_PI;
            int theta_i = -1;
            {
               const double theta_deg = theta * 180.0 / M_PI;
               double best_err = 1e9;
               int best = -1;
               for(int t = 0; t < 6; ++t){
                  const double can = 60.0 * t;
                  double err = std::fabs(theta_deg - can);
                  if(err > 180.0) err = 360.0 - err;
                  if(err < best_err){ best_err = err; best = t; }
               }
               if(best_err <= 1.0) theta_i = best;
            }
            if(theta_i < 0 || theta_i > 5) std::cout << k << ", " << i << ", " << j << ", " << theta << ", " << theta_i << ", " << dx << ", " << dy << std::endl;
            // int theta_i = (theta/30.0 == 5.0) ? (2) : ((theta/30.0 == 4.0) ? (1) : (theta/30.0 == 3.0 ? (0):-1) );
            Eij_2NN.at(i*100 + j).at(theta_i*4 + 0) = J*J_constant;
            Eij_2NN[i*100 + j][theta_i*4 + 1] = DMI_inter_scaling*Dx*J_constant;
            Eij_2NN[i*100 + j][theta_i*4 + 2] = DMI_inter_scaling*Dy*J_constant;
            Eij_2NN[i*100 + j][theta_i*4 + 3] = DMI_inter_scaling*Dz*J_constant;
            //  std::cout  <<  Eintra_Cr1_2NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr1_2NN.at(int_x).at(int_y).at(theta_i*4 + 1) << ", " << Eintra_Cr1_2NN.at(int_x).at(int_y).at(theta_i*4 + 2) << ", " << Eintra_Cr1_2NN.at(int_x).at(int_y).at(theta_i*4 + 3) << std::endl;
            count++;
            }
        }
        
    }
     std::cout << "2NN by thread " << omp_get_thread_num() << " count: " << count << std::endl;
    
    for(int k = 0; k < 3; k++) {
        for(int i = 0; i < 100; i++) {
           
           for(int j = 0; j < 100; j++) {
                getline(ifile1,line);
                std::stringstream liness(line.c_str());
                double dx;
                double dy;
                
                double J;
                double Dx;
                double Dy;
                double Dz;
  
                double theta;
                liness >> dx >> dy >> J >> Dx >> Dy >> Dz;

                // Honeycomb 3NN: same |θ| set / codes as 1NN
                theta = std::abs(atan2(dy,dx))*180.0/M_PI;
                int theta_i = -1;
                {
                   const double cans[3] = {30.0, 90.0, 150.0};
                   const int codes[3] = {1, 0, 2};
                   double best_err = 1e9;
                   for(int t = 0; t < 3; ++t){
                      const double err = std::fabs(theta - cans[t]);
                      if(err < best_err){ best_err = err; theta_i = codes[t]; }
                   }
                   if(best_err > 1.0) theta_i = -1;
                }
                if(theta_i < 0) std::cout << k << ", " << i << ", " << j << ", " << theta << ", " << dx << ", " << dy << std::endl;
                Eij_3NN.at(i*100 + j).at(theta_i*4 + 0) = J*J_constant;
                Eij_3NN[i*100 + j][theta_i*4 + 1] = DMI_inter_scaling*Dx*J_constant;
                Eij_3NN[i*100 + j][theta_i*4 + 2] = DMI_inter_scaling*Dy*J_constant;
                Eij_3NN[i*100 + j][theta_i*4 + 3] = DMI_inter_scaling*Dz*J_constant;
            // std::cout  <<  Eintra_Cr1_3NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr1_3NN.at(int_x).at(int_y).at(theta_i*4 + 1) << ", " << Eintra_Cr1_3NN.at(int_x).at(int_y).at(theta_i*4 + 2) << ", " << Eintra_Cr1_3NN.at(int_x).at(int_y).at(theta_i*4 + 3) << std::endl;
                count++;
            }
        }
    }
     std::cout << "3NN by thread " << omp_get_thread_num() << " count: " << count << std::endl;
    // if(count != 100*100*12) std::cout << count  << std::endl;
    std::ofstream intra_out("Intra_exchange_out_" + std::to_string(omp_get_thread_num()) + ".txt");
    for(int i = 0; i < Eij_3NN.size(); i++) {
        intra_out << int(i/100) << ", " << i%100  ;
        for(int j = 0; j < Eij_3NN[0].size(); j++) {intra_out << ", " <<  Eij_3NN[i][j]; }
        intra_out << "\n";
    }
    intra_out.close();

   ifile1.close();
}


