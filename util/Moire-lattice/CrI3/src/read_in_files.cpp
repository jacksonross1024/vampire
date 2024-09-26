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
    if(!ifile2.is_open()) {std::cerr << filename << " is not open" << std::endl; exit(1);}
   //optional control
//    double x_step = 0.05;
//    double y_step = 0.0866;
    for(int i=0; i<40000; i++){
        getline(ifile2,line);
        std::stringstream liness(line.c_str());
        double ii;
        double ij;
        liness >> ij >> ii >> Jij[int(ii)+100][int(ij)+100];
        // std::cout << ii << ", " << ii*2 << std::endl;
   //   outfile <<atom_id[i] << "\t" << x_in[i] << '\t' << y_in[i] << "\t" << z_in[i] << "\t" << S_in[i] << "\t" << temp << "\t" << temp << std::endl;
   //   std::cout <<ii << "\t" << ij << "\t" <<  Jint[ii][ij] << std::endl;
 }

    // for(int i = 0; i < 201; i+=2) {
    //     for(int j = 1; j < 200; j+=2) {
    //         Jij[i][j] = 0.5*(Jij[i][j-1]+Jij[i][j+1]);
    //     }
    // }
    // for(int j = 0; j < 201; j++) {
    //     for(int i = 1; i < 201; i+=2) {
    //         Jij[i][j] = 0.5*(Jij[i-1][j]+Jij[i+1][j]);
    //     }
    // }
    ifile2.close();
    // double 
    // double rad_120 = 120.0*M_PI/180.0;
    // double rad_210 = -120.0*M_PI/180.0;
    std::ofstream out_file(filename + "_out.txt");
   if(!out_file.is_open()) {std::cerr << filename + "_out.txt is not open" << std::endl; exit(1);}
   for(int i = 0; i < Jij.size(); i++) {
    for(int j = 0; j < Jij[i].size(); j++) {
        if((j > 100) && (i*0.0693 > (6.93+3.465-j*0.03465)) ) {
            // out_file << i << "\t" << j << "\t" << Jij[i][j] << "\n";
            double x = i*0.0693 - a0x;
            double y = j*0.0693- a0x;
            double theta = atan2(y,x);
            double r = sqrt(x*x+y*y);
            double dx = r*cos(theta+2*M_PI/3.0);//-y*sin(rad_120);
            double dy = r*sin(theta+2*M_PI/3.0);//+x*sin(rad_120);
            double dx2 = r*cos(theta + 4*M_PI/3.0);//-y*sin(rad_210);
            double dy2 = r*sin(theta + 4*M_PI/3.0);//+x*sin(rad_210);
            int dx_int = round(100*dx/a0x)+100;
            int dy_int = round(100*dy/a0x)+100;
            int dx2_int = round(100*dx2/a0x)+100;
            int dy2_int = round(100*dy2/a0x)+100;

            if(dx_int > 200 || dx_int < 0 || dy_int > 200 || dy_int < 0 || dx2_int > 200 || dx2_int < 0 || dy2_int > 200 || dy2_int < 0) continue;
            out_file << i << "\t" << j << "\t" << Jij[i][j] << "\n";
            Jij[dx_int][dy_int] = Jij[i][j];
            Jij[dx2_int][dy2_int] = Jij[i][j];

            
            out_file << dx_int << "\t" << dy_int << "\t" << Jij[dx_int][dy_int] << "\n";
            out_file << dx2_int << "\t" << dy2_int << "\t" << Jij[dx2_int][dy2_int] << "\n";
            // out_file << dx_int << "\t" << dy_int << "\t" << Jij[i][j] << "\n";
            // out_file << dx2_int << "\t" << dy2_int << "\t" << Jij[i][j] << "\n";
            
        }
        // else if((i*0.0693-6.93-3.465) < (0- 0.03465*j) && i<100) out_file << i << "\t" << j << "\t" << Jij[Jij.size()-i][Jij[i].size()-j] << "\n";
        // else  out_file << i << "\t" << j << "\t" << Jij[200+(200-i)][200+(200-j)] << "\n";
    }
   }
   out_file.close();
}


void read_in_dmi(std::string filename, std::vector < std::vector < double > > &Dx, std::vector < std::vector < double > > &Dy, std::vector < std::vector < double > > &Dz){
   std::ifstream ifile2(filename);
   std::string line;
    if(!ifile2.is_open()) {std::cerr << filename << " is not open" << std::endl; exit(1);}
   for(int i=0; i<40000; i++){
     getline(ifile2,line);
     std::stringstream liness(line.c_str());
     double ii;
     double ij;
     liness >> ij >> ii >> Dx[ii+100][ij+100] >> Dy[ii+100][ij+100] >> Dz[ii+100][ij+100];
   //   outfile <<atom_id[i] << "\t" << x_in[i] << '\t' << y_in[i] << "\t" << z_in[i] << "\t" << S_in[i] << "\t" << temp << "\t" << temp << std::endl;
   // std::cout <<ii << "\t" << ij << "\t" <<  Dx[ii][ij]<< "\t" <<  Dy[ii][ij]<< "\t" <<  Dz[ii][ij] << std::endl;
 }
   ifile2.close();

   std::ofstream out_file(filename + "_out.txt");
   if(!out_file.is_open()) {std::cerr << filename + "_out.txt is not open" << std::endl; exit(1);}
   for(int i = 0; i < Dx.size(); i++){
    for(int j = 0; j < Dx[i].size(); j++) {
        out_file << i << "\t" << j << "\t" << Dx[i][j] << "\t" << Dy[i][j] << "\t" << Dz[i][j] << "\n";
    }
   }
   out_file.close();
}

void read_in_dft(std::string filename) {
    std::ifstream ifile2(filename);
    std::string line;

    double x_step = 0.05;
    double y_step = 0.0866;
    
    for(int i=0; i<331; i++){
        getline(ifile2,line);
        std::stringstream liness(line.c_str());
        double ii;
        double ij;
        double Jij;
        double Jii;
        double Jii2;
        double Dij_x;
        double Dij_y;
        double Dij_z;
        double Dii90_x;
        double Dii90_y;
        double Dii90_z;
        double Dii210_x;
        double Dii210_y;
        double Dii210_z;
        double Dii330_x;
        double Dii330_y;
        double Dii330_z;
        
        liness >> ii >> ij >> Jii >> Jii2 >> Jij >> Dij_x >> Dij_y >> Dij_z >> Dii90_x >> Dii90_y >> Dii90_z >> Dii210_x >> Dii210_y >> Dii210_z >> Dii330_x >> Dii330_y >> Dii330_z;
        std::vector<double> Dvector_ij(9);
        Dvector_ij = {Dij_x, Dij_y, Dij_z};
        // std::cout << ii << " , " << ij << ", " << int(round(ii/x_step))+20 << ", " << int(round(ij/y_step))+10 << std::endl;
        ii = int(round(ii/x_step))+20;
        ij = int(round(ij/y_step))+10;
        
        Jinter.at(ii*5).at(ij*10) = Jij;
        Jintra1.at(ii*5).at(ij*10) = Jii;
        Jintra2.at(ii*5).at(ij*10) = Jii2;
        // D_inter.at(ii*5).at(ij*10) = Dvector_ij;
        // std::vector<double> Dvector_ii(9);
        // Dvector_ii = {Dii90_x, Dii90_y, Dii90_z,\
        //                         Dii210_x, Dii210_y, Dii210_z, \
        //                         Dii330_x, Dii330_y, Dii330_z };
        // D_intra.at(ii*5).at(ij*10) = Dvector_ii;
    }

    for(int i = 0; i < Jinter.size()-1; i+=5){
        for(int j = 0; j < Jinter.size()-1; j+=10){
            double Jinter_x_1 = Jinter[i][j];
            double Jinter_x_2 = Jinter[i+5][j];
            double Jinter_dx = (Jinter_x_2 - Jinter_x_1)/5.0;
            double Jinter_y_1 = Jinter[i][j];
            double Jinter_y_2 = Jinter[i][j+10];
            double Jinter_dy = (Jinter_y_2-Jinter_y_1)/10.0;

            double Jintra_x_1 = Jintra1[i][j];
            double Jintra_x_2 = Jintra1[i+5][j];
            double Jintra_dx = (Jintra_x_2 - Jintra_x_1)/5.0;
            double Jintra_y_1 = Jintra1[i][j];
            double Jintra_y_2 = Jintra1[i][j+10];
            double Jintra_dy = (Jintra_y_2-Jintra_y_1)/10.0;

            double Jintra2_x_1 = Jintra2[i][j];
            double Jintra2_x_2 = Jintra2[i+5][j];
            double Jintra2_dx = (Jintra_x_2 - Jintra_x_1)/5.0;
            double Jintra2_y_1 = Jintra2[i][j];
            double Jintra2_y_2 = Jintra2[i][j+10];
            double Jintra2_dy = (Jintra_y_2-Jintra_y_1)/10.0;

            for(int iprime = 0; iprime < 5; iprime++){
                for(int jprime = 0; jprime < 10; jprime++) {
                    double x_weight = iprime/5.0;
                    double y_weight = jprime/10.0;
                    double norm = (x_weight == 0 && y_weight==0) ? (1.0) : (x_weight+y_weight);
                    Jinter[i+iprime][j+jprime] = (x_weight*Jinter_dx*iprime+Jinter_dy*jprime*y_weight)/norm+Jinter[i][j];
                    Jintra1[i+iprime][j+jprime] = 0.5*(Jintra_dx*iprime+Jintra_dy*jprime)+Jintra_x_1;
                    Jintra2[i+iprime][j+jprime] = 0.5*(Jintra2_dx*iprime+Jintra2_dy*jprime)+Jintra2_x_1;
                }
            }
        }
    }

    std::ofstream Jinter_out;
    std::ofstream Jintra1_out;
    std::ofstream Jintra2_out;
    // std::ofstream Dinter_out;
    // std::ofstream Dintra1_out;
    // std::ofstream Dintra2_out;
    Jinter_out.open("Jinter_interpolation_out.txt");
    Jintra1_out.open("Jintra1_interpolation_out.txt");
    Jintra2_out.open("Jintra2_interpolation_out.txt");
    // Dinter_out.open("Dinter_interpolation_out.txt");
    // Dintra1_out.open("Dintra1_interpolation_out.txt");
    // Dintra2_out.open("Dintra2_interpolation_out.txt");
    for(int i = 0; i < 201; i++) {
        for(int j = 0; j < 201; j++) {
            Jinter_out << i << ", " << j << ", " << Jinter[i][j] << "\n";
            Jintra1_out << i << ", " << j << ", " << Jintra1[i][j] << "\n";
            Jintra2_out << i << ", " << j << ", " << Jintra2[i][j] << "\n";
            // Dinter_out << i << ", " << j << ", " << Dx_inter[i][j] << ", " << Dy_inter[i][j] << ", " << Dz_inter[i][j] << ", " << sqrt(Dx_inter[i][j]*Dx_inter[i][j] \
            //                                                                                                                          + Dy_inter[i][j]*Dy_inter[i][j] \
            //                                                                                                                          + Dz_inter[i][j]*Dz_inter[i][j]) << "\n";
            // Dintra1_out << i << ", " << j << ", " << Dx_intra[i][j] << ", " << Dy_intra[i][j] << ", " << Dz_intra[i][j] <<  ", " << sqrt(Dx_intra[i][j]*Dx_intra[i][j] \
                                                                                                                                     + Dy_intra[i][j]*Dy_intra[i][j] \
                                                                                                                                     + Dz_intra[i][j]*Dz_intra[i][j]) << "\n";
            // Dintra2_out << i << ", " << j << ", " << Jintra2[i][j] << "\n";
        }
    }
    Jinter_out.close();
    Jintra1_out.close();
    Jintra2_out.close();
    // Dinter_out.close();
    // Dintra1_out.close();   
   //   outfile <<atom_id[i] << "\t" << x_in[i] << '\t' << y_in[i] << "\t" << z_in[i] << "\t" << S_in[i] << "\t" << temp << "\t" << temp << std::endl;
   // std::cout <<ii << "\t" << ij << "\t" <<  Dx[ii][ij]<< "\t" <<  Dy[ii][ij]<< "\t" <<  Dz[ii][ij] << std::endl;
    
    ifile2.close();
    exit(1);
}
