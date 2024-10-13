#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include "positions.hpp"
#include "exchange.hpp"
#include "initialise.hpp"
#include <array>

int main(int argc, char* argv[]){

    std::string dmi_check = "--dmi";
    if(argc < 2) {std::cerr << "need twist angle even if zero" << std::endl; exit(1);}
    for(int a = 1; a < argc; a++) {
        if(a == 1) {twist_angle = atof(argv[a]); // 1.1 
        std::cout << "twist angle: "<< twist_angle << std::endl;}
        if(a == 2) {max_range = atof(argv[a]);
        std::cout << "max inter exchange range: " << max_range << std::endl;}
        if(a == 3) {
            if(argv[a]== dmi_check) {DMI = true;
            std::cout << " with DMI " << std::endl;}
        }
    }

   system_size_x = 2000; //4000
   system_size_y = 2000; //4000
   number_of_unit_cells_z = 1; //2

   // set up new material constants
   dmi12 = 1.0; // DMI constant between layers 1-2
   dmi23 = 1.0; // DMI constant between layers 2-3
   dmi34 = 1.0; // DMI constant between layers 3-4
   
   dmi_decay = 1.0; // distance-dependent DMI

   exchange12 = 1.0; // exchange constant between layers 1-2
   exchange23 = 1.0; // exchange constant between layers 2-3
   exchange34 = 1.0; // exchange constant between layers 3-4

   separation = 0.0; // distance between layers 2-3
   // double_bilayer = false;
   // if(double_bilayer) pristine_bilayer_type = "baab";

   initialise_variables();

   twist_loction = 2*system_size_z/5 -0.01;
    std::cout << "twisting at: " << twist_loction << std::endl;
    read_in_atoms("files/atom_list_aa_rhombic", num_atoms, atom);
   //  read_in_dft("files/criteria.txt");
    //  read_in_atoms("files/nm_atoms", num_nm_atoms, nm_atom);
   //  read_in_exchange("files/Interpolated_J_Inter", Jinter);
   read_in_inter_exchanges("bilayer_sliding/Cr1_inter.txt", Einter_Cr1);
   read_in_inter_exchanges("bilayer_sliding/Cr2_inter.txt", Einter_Cr2);
   read_in_inter_exchanges("bilayer_sliding/Cr3_inter.txt", Einter_Cr3);
   read_in_inter_exchanges("bilayer_sliding/Cr4_inter.txt", Einter_Cr4);

   // read_in_intra_exchanges("bilayer_sliding/Cr1_intra.txt", Eintra_Cr1_1NN, Eintra_Cr1_2NN, Eintra_Cr1_3NN);
   // read_in_intra_exchanges("bilayer_sliding/Cr2_intra.txt", Eintra_Cr2_1NN, Eintra_Cr2_2NN, Eintra_Cr2_3NN);
   // read_in_intra_exchanges("bilayer_sliding/Cr3_intra.txt", Eintra_Cr3_1NN, Eintra_Cr3_2NN, Eintra_Cr3_3NN);
   // read_in_intra_exchanges("bilayer_sliding/Cr4_intra.txt", Eintra_Cr4_1NN, Eintra_Cr4_2NN, Eintra_Cr4_3NN);
   std::ifstream ifile1("bilayer_sliding/Cr1_intra.txt");
    std::string line;
    if(!ifile1.is_open()) {std::cerr  << " is not open" << std::endl; exit(1);}
   for(int i=0; i<Eintra_Cr1_1NN.size()*11; i++){
        for(int j = 0; j < 3; j++) {
            getline(ifile1,line);
            std::stringstream liness(line.c_str());
            double r;
            double dx;
            double dy;
            double dz;
            double J;
            double Dx;
            double Dy;
            double Dz;
            double sx;
            double sy;
            double theta;
            liness >> r >> dx >> dy >> dz >> J >> Dx >> Dy >> Dz >> sx >> sy;

            theta = std::abs(atan2(dy,dx))*180.0/M_PI;

            int int_x = int(sx*10);
            int int_y = int(sy*10);
                        // std::cout << int_x << ", " << int_y << std::endl;
            int theta_i = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0)== 1.0) ? (1) : (round(theta/30.0 == 3.0) ? (0):-1) );
            if(theta_i < 0) std::cout << round(theta/30.0) << ", " << dx << ", " << dy << std::endl;
            Eintra_Cr1_1NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eintra_Cr1_1NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eintra_Cr1_1NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eintra_Cr1_1NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
            //  std::cout  <<  Eintra_Cr1_1NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr1_1NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eintra_Cr1_1NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eintra_Cr1_1NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
        }
        for(int j = 0; j < 6; j++) {
            getline(ifile1,line);
            std::stringstream liness(line.c_str());
            double r;
            double dx;
            double dy;
            double dz;
            double J;
            double Dx;
            double Dy;
            double Dz;
            double sx;
            double sy;
            double theta;

            liness >> r >> dx >> dy >> dz >> J >> Dx >> Dy >> Dz >> sx >> sy;
            theta = atan2(dy,dx);//*180.0/M_PI;
            int theta_i = int(round(((theta < 0.0) ? (theta+=2*M_PI) : (theta)) *180.0/M_PI/60.0));
           
            // if(theta_i == 0) std::cout <<  theta << ", " << dx << ", " << dy << std::endl;
            int int_x = int(sx*10);
            int int_y = int(sy*10);
            // int theta_i = (theta/30.0 == 5.0) ? (2) : ((theta/30.0 == 4.0) ? (1) : (theta/30.0 == 3.0 ? (0):-1) );
            Eintra_Cr1_2NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eintra_Cr1_2NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eintra_Cr1_2NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eintra_Cr1_2NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
            //  std::cout  <<  Eintra_Cr1_2NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr1_2NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eintra_Cr1_2NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eintra_Cr1_2NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
        }
        for(int j = 0; j < 3; j++) {
            getline(ifile1,line);
            std::stringstream liness(line.c_str());
            double r;
            double dx;
            double dy;
            double dz;
            double J;
            double Dx;
            double Dy;
            double Dz;
            double sx;
            double sy;
            double theta;
            liness >> r >> dx >> dy >> dz >> J >> Dx >> Dy >> Dz >> sx >> sy;
            theta = std::abs(atan2(dy,dx))*180.0/M_PI;
            // std::cout << round(theta/30.0) << ", " << dx << ", " << dy << std::endl;
            int int_x = int(sx*10);
            int int_y = int(sy*10);
            int theta_i = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0) == 1.0) ? (1) : (round(theta/30.0) == 3.0 ? (0):-1) );
            Eintra_Cr1_3NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eintra_Cr1_3NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eintra_Cr1_3NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eintra_Cr1_3NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
            // std::cout  <<  Eintra_Cr1_3NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr1_3NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eintra_Cr1_3NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eintra_Cr1_3NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
        }
    }

   ifile1.close();

   std::ifstream ifile2("bilayer_sliding/Cr2_intra.txt");
   //  std::string line;
    if(!ifile2.is_open()) {std::cerr  << " is not open" << std::endl; exit(1);}
   for(int i=0; i<Eintra_Cr2_1NN.size()*11; i++){
        for(int j = 0; j < 3; j++) {
            getline(ifile2,line);
            std::stringstream liness(line.c_str());
            double r;
            double dx;
            double dy;
            double dz;
            double J;
            double Dx;
            double Dy;
            double Dz;
            double sx;
            double sy;
            double theta;
            liness >> r >> dx >> dy >> dz >> J >> Dx >> Dy >> Dz >> sx >> sy;

            theta = std::abs(atan2(dy,dx))*180.0/M_PI;
            // std::cout << round(theta/30.0) << ", " << dx << ", " << dy << std::endl;
            int int_x = int(sx*10);
            int int_y = int(sy*10);
            int theta_i = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0)== 1.0) ? (1) : (round(theta/30.0 == 3.0) ? (0):-1) );
            if(theta_i < 0) std::cout << round(theta/30.0) << ", " << dx << ", " << dy << std::endl;
            Eintra_Cr2_1NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eintra_Cr2_1NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eintra_Cr2_1NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eintra_Cr2_1NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
            //  std::cout  <<  Eintra_Cr2_1NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr2_1NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eintra_Cr2_1NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eintra_Cr2_1NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
        }
        for(int j = 0; j < 6; j++) {
            getline(ifile2,line);
            std::stringstream liness(line.c_str());
            double r;
            double dx;
            double dy;
            double dz;
            double J;
            double Dx;
            double Dy;
            double Dz;
            double sx;
            double sy;
            double theta;

            liness >> r >> dx >> dy >> dz >> J >> Dx >> Dy >> Dz >> sx >> sy;
            theta = atan2(dy,dx);//*180.0/M_PI;
            int theta_i = int(round(((theta < 0.0) ? (theta+=2*M_PI) : (theta)) *180.0/M_PI/60.0));
           
            // if(theta_i == 0) std::cout <<  theta << ", " << dx << ", " << dy << std::endl;
            int int_x = int(sx*10);
            int int_y = int(sy*10);
            // int theta_i = (theta/30.0 == 5.0) ? (2) : ((theta/30.0 == 4.0) ? (1) : (theta/30.0 == 3.0 ? (0):-1) );
            Eintra_Cr2_2NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eintra_Cr2_2NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eintra_Cr2_2NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eintra_Cr2_2NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
            //  std::cout  <<  Eintra_Cr2_2NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr2_2NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eintra_Cr2_2NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eintra_Cr2_2NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
        }
        for(int j = 0; j < 3; j++) {
            getline(ifile2,line);
            std::stringstream liness(line.c_str());
            double r;
            double dx;
            double dy;
            double dz;
            double J;
            double Dx;
            double Dy;
            double Dz;
            double sx;
            double sy;
            double theta;
            liness >> r >> dx >> dy >> dz >> J >> Dx >> Dy >> Dz >> sx >> sy;
            theta = std::abs(atan2(dy,dx))*180.0/M_PI;
            // std::cout << round(theta/30.0) << ", " << dx << ", " << dy << std::endl;
            int int_x = int(sx*10);
            int int_y = int(sy*10);
            int theta_i = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0) == 1.0) ? (1) : (round(theta/30.0) == 3.0 ? (0):-1) );
            Eintra_Cr2_3NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eintra_Cr2_3NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eintra_Cr2_3NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eintra_Cr2_3NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
            // std::cout  <<  Eintra_Cr2_3NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr2_3NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eintra_Cr2_3NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eintra_Cr2_3NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
        }
    }
    ifile2.close();

    std::ifstream ifile3("bilayer_sliding/Cr3_intra.txt");
   //  std::string line;
    if(!ifile3.is_open()) {std::cerr  << " is not open" << std::endl; exit(1);}
   for(int i=0; i<Eintra_Cr3_1NN.size()*11; i++){
        for(int j = 0; j < 3; j++) {
            getline(ifile3,line);
            std::stringstream liness(line.c_str());
            double r;
            double dx;
            double dy;
            double dz;
            double J;
            double Dx;
            double Dy;
            double Dz;
            double sx;
            double sy;
            double theta;
            liness >> r >> dx >> dy >> dz >> J >> Dx >> Dy >> Dz >> sx >> sy;

            theta = std::abs(atan2(dy,dx))*180.0/M_PI;
            // std::cout << round(theta/30.0) << ", " << dx << ", " << dy << std::endl;
            int int_x = int(sx*10);
            int int_y = int(sy*10);
            int theta_i = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0)== 1.0) ? (1) : (round(theta/30.0 == 3.0) ? (0):-1) );
            if(theta_i < 0) std::cout << round(theta/30.0) << ", " << dx << ", " << dy << std::endl;
            Eintra_Cr3_1NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eintra_Cr3_1NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eintra_Cr3_1NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eintra_Cr3_1NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
            // std::cout  <<  Eintra_Cr3_1NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr3_1NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eintra_Cr3_1NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eintra_Cr3_1NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
        }
        for(int j = 0; j < 6; j++) {
            getline(ifile3,line);
            std::stringstream liness(line.c_str());
            double r;
            double dx;
            double dy;
            double dz;
            double J;
            double Dx;
            double Dy;
            double Dz;
            double sx;
            double sy;
            double theta;

            liness >> r >> dx >> dy >> dz >> J >> Dx >> Dy >> Dz >> sx >> sy;
            theta = atan2(dy,dx);//*180.0/M_PI;
            int theta_i = int(round(((theta < 0.0) ? (theta+=2*M_PI) : (theta)) *180.0/M_PI/60.0));
           
            // if(theta_i == 0) std::cout <<  theta << ", " << dx << ", " << dy << std::endl;
            int int_x = int(sx*10);
            int int_y = int(sy*10);
            // int theta_i = (theta/30.0 == 5.0) ? (2) : ((theta/30.0 == 4.0) ? (1) : (theta/30.0 == 3.0 ? (0):-1) );
            Eintra_Cr3_2NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eintra_Cr3_2NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eintra_Cr3_2NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eintra_Cr3_2NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
            //  std::cout  <<  Eintra_Cr3_2NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr3_2NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eintra_Cr3_2NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eintra_Cr3_2NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
        }
        for(int j = 0; j < 3; j++) {
            getline(ifile3,line);
            std::stringstream liness(line.c_str());
            double r;
            double dx;
            double dy;
            double dz;
            double J;
            double Dx;
            double Dy;
            double Dz;
            double sx;
            double sy;
            double theta;
            liness >> r >> dx >> dy >> dz >> J >> Dx >> Dy >> Dz >> sx >> sy;
            theta = std::abs(atan2(dy,dx))*180.0/M_PI;
            // std::cout << round(theta/30.0) << ", " << dx << ", " << dy << std::endl;
            int int_x = int(sx*10);
            int int_y = int(sy*10);
            int theta_i = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0) == 1.0) ? (1) : (round(theta/30.0) == 3.0 ? (0):-1) );
            Eintra_Cr3_3NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eintra_Cr3_3NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eintra_Cr3_3NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eintra_Cr3_3NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
            // std::cout  <<  Eintra_Cr3_3NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr3_3NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eintra_Cr3_3NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eintra_Cr3_3NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
        }
    }
    ifile3.close();

    std::ifstream ifile4("bilayer_sliding/Cr4_intra.txt");
   //  std::string line;
    if(!ifile4.is_open()) {std::cerr  << " is not open" << std::endl; exit(1);}
   for(int i=0; i<Eintra_Cr4_1NN.size()*11; i++){
        for(int j = 0; j < 3; j++) {
            getline(ifile4,line);
            std::stringstream liness(line.c_str());
            double r;
            double dx;
            double dy;
            double dz;
            double J;
            double Dx;
            double Dy;
            double Dz;
            double sx;
            double sy;
            double theta;
            liness >> r >> dx >> dy >> dz >> J >> Dx >> Dy >> Dz >> sx >> sy;

            theta = std::abs(atan2(dy,dx))*180.0/M_PI;
            // std::cout << round(theta/30.0) << ", " << dx << ", " << dy << std::endl;
            int int_x = int(sx*10);
            int int_y = int(sy*10);
            int theta_i = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0)== 1.0) ? (1) : (round(theta/30.0 == 3.0) ? (0):-1) );
            if(theta_i < 0) std::cout << round(theta/30.0) << ", " << dx << ", " << dy << std::endl;
            Eintra_Cr4_1NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eintra_Cr4_1NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eintra_Cr4_1NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eintra_Cr4_1NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
            // std::cout  <<  Eintra_Cr4_1NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr4_1NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eintra_Cr4_1NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eintra_Cr4_1NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
        }
        for(int j = 0; j < 6; j++) {
            getline(ifile4,line);
            std::stringstream liness(line.c_str());
            double r;
            double dx;
            double dy;
            double dz;
            double J;
            double Dx;
            double Dy;
            double Dz;
            double sx;
            double sy;
            double theta;

            liness >> r >> dx >> dy >> dz >> J >> Dx >> Dy >> Dz >> sx >> sy;
            theta = atan2(dy,dx);//*180.0/M_PI;
            int theta_i = int(round(((theta < 0.0) ? (theta+=2*M_PI) : (theta)) *180.0/M_PI/60.0));
           
            // if(theta_i == 0) std::cout <<  theta << ", " << dx << ", " << dy << std::endl;
            int int_x = int(sx*10);
            int int_y = int(sy*10);
            // int theta_i = (theta/30.0 == 5.0) ? (2) : ((theta/30.0 == 4.0) ? (1) : (theta/30.0 == 3.0 ? (0):-1) );
            Eintra_Cr4_2NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eintra_Cr4_2NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eintra_Cr4_2NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eintra_Cr4_2NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
            //  std::cout  <<  Eintra_Cr4_2NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr4_2NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eintra_Cr4_2NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eintra_Cr4_2NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
        }
        for(int j = 0; j < 3; j++) {
            getline(ifile4,line);
            std::stringstream liness(line.c_str());
            double r;
            double dx;
            double dy;
            double dz;
            double J;
            double Dx;
            double Dy;
            double Dz;
            double sx;
            double sy;
            double theta;
            liness >> r >> dx >> dy >> dz >> J >> Dx >> Dy >> Dz >> sx >> sy;
            theta = std::abs(atan2(dy,dx))*180.0/M_PI;
            // std::cout << round(theta/30.0) << ", " << dx << ", " << dy << std::endl;
            int int_x = int(sx*10);
            int int_y = int(sy*10);
            int theta_i = (round(theta/30.0) == 5.0) ? (2) : ((round(theta/30.0) == 1.0) ? (1) : (round(theta/30.0) == 3.0 ? (0):-1) );
            Eintra_Cr4_3NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eintra_Cr4_3NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eintra_Cr4_3NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eintra_Cr4_3NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
            // std::cout  <<  Eintra_Cr4_3NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eintra_Cr4_3NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eintra_Cr4_3NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eintra_Cr4_3NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
        }
    }
    ifile4.close();
   // for(int i = 0; i < Eintra_Cr1_1NN.size(); i++) {
   //    for(int j = 0; j < Eintra_Cr1_1NN[i].size(); j++) {
   //       for(int k = 0; k < Eintra_Cr1_1NN[i][j].size(); k++) {
   //          std::cout <<  Eintra_Cr1_1NN[i][j][k][0] << ", " << Eintra_Cr1_1NN[i][j][k][1] << ", " <<  Eintra_Cr1_1NN[i][j][k][2] << ", " << Eintra_Cr1_1NN[i][j][k][3] << std::endl;
   //          std::cout <<  Eintra_Cr2_1NN[i][j][k][0] << ", " << Eintra_Cr2_1NN[i][j][k][1] << ", " <<  Eintra_Cr2_1NN[i][j][k][2] << ", " << Eintra_Cr2_1NN[i][j][k][3] << std::endl;
   //          std::cout <<  Eintra_Cr3_1NN[i][j][k][0] << ", " << Eintra_Cr3_1NN[i][j][k][1] << ", " <<  Eintra_Cr3_1NN[i][j][k][2] << ", " << Eintra_Cr3_1NN[i][j][k][3] << std::endl;
   //          std::cout <<  Eintra_Cr4_1NN[i][j][k][0] << ", " << Eintra_Cr4_1NN[i][j][k][1] << ", " <<  Eintra_Cr4_1NN[i][j][k][2] << ", " << Eintra_Cr4_1NN[i][j][k][3] << std::endl;
   //       }
   //    }
   // }
   // exit(1);
   //  read_in_exchange("files/Interpolated_J1_Intra", Jintra1);
   //  read_in_exchange("files/Interpolated_J2_Intra", Jintra2);
   //  read_in_dmi("files/Interpolated_1st_Dij_Intra", Dx_intra, Dy_intra, Dz_intra, true);
   //  read_in_dmi("files/Interpolated_2nd_Dij_Intra", Dx_intra2, Dy_intra2, Dz_intra2, false);
   //  read_in_dmi("files/interpolated_1st_Dij_intra.txt", Dx_intra2, Dy_intra2, Dz_intra2);
   //  read_in_dmi("files/Interpolated_Dij_Inter", Dx_inter, Dy_inter, Dz_inter, true);

    create_magnetic_atom_list("atom_positions.ucf");
    //  create_nm_atom_list();

    print_header();
   //  //  calc_in_plane_exchange(row1);
   //  //  calc_in_plane_exchange(row2);
   //  //  calc_in_plane_exchange(row3);
   //  //  calc_in_plane_exchange(row4);

   //  //  // The order of these exchange calculations is important, as dx,dy shifts are NOT symmetric
   //  //  calc_out_of_plane_exchange(row4,row3);
   //  //  calc_out_of_plane_exchange(row3,row2);
   //  //  calc_out_of_plane_exchange(row2,row1);
    calc_interactions();
    print_interaction_header();

 }
