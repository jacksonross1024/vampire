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
     liness >> atom2[i].id >> atom2[i].x >> atom2[i].y >> atom2[i].z >> atom2[i].S >> atom2[i].l_id >> atom2[i].h_id;
     //   outfile <<atom_id[i] << "\t" << x_in[i] << '\t' << y_in[i] << "\t" << z_in[i] << "\t" << S_in[i] << "\t" << temp << "\t" << temp << std::endl;

   }
   ifile.close();
}

void read_in_exchange(std::string filename, std::vector<std::vector<double>> & Jij){
   std::ifstream ifile2(filename);
   std::string line;
    if(!ifile2.is_open()) {std::cerr << filename << " is not open" << std::endl; exit(1);}
    
    for(int i=0; i<(Jij.size()*Jij.size()); i++){
        getline(ifile2,line);
        std::stringstream liness(line.c_str());
        double ii;
        double ij;
        liness >> ij >> ii >> Jij[int(ii)+100][int(ij)+100];
    }

    ifile2.close();
    std::ofstream out_file(filename + "_out.txt");
   if(!out_file.is_open()) {std::cerr << filename + "_out.txt is not open" << std::endl; exit(1);}
   for(int i = 0; i < Jij.size(); i++) {
    for(int j = 0; j < Jij[i].size(); j++) {
        // if((j > 100) && (i*0.0693 > (6.93+3.465-j*0.03465)) ) {
        if(j > 100) {
            // out_file << i << "\t" << j << "\t" << Jij[i][j] << "\n";
            // double x = i*0.0693 - a0x;
            // double y = j*0.06002 - a1y;
            // double theta = atan2(y,x);
            // double r = sqrt(x*x+y*y);
            // double dx = -x;// r*cos(theta+2*M_PI/3.0);//-y*sin(rad_120);
            // double dy = j;//r*sin(theta+2*M_PI/3.0);//+x*sin(rad_120);
            // double dx2 = r*cos(theta + 4*M_PI/3.0);//-y*sin(rad_210);
            // double dy2 = r*sin(theta + 4*M_PI/3.0);//+x*sin(rad_210);
            int dx_int = i ;//round(100*dx/a0x)+100;
            int dy_int = 200-j;//round(100*dy/a0x)+100;
            // int dx2_int = round(100*dx2/a0x)+100;
            // int dy2_int = round(100*dy2/a0x)+100;

            if(dx_int > 200 || dx_int < 0 || dy_int > 200 || dy_int < 0 ) continue;
            out_file << i << "\t" << j << "\t" << Jij[i][j] << "\n";
            Jij[dx_int][dy_int] = Jij[i][j];
            // Jij[dx2_int][dy2_int] = Jij[i][j];

            
            out_file << dx_int << "\t" << dy_int << "\t" << Jij[dx_int][dy_int] << "\n";
            // out_file << dx2_int << "\t" << dy2_int << "\t" << Jij[dx2_int][dy2_int] << "\n";
            // out_file << dx_int << "\t" << dy_int << "\t" << Jij[i][j] << "\n";
            // out_file << dx2_int << "\t" << dy2_int << "\t" << Jij[i][j] << "\n";
            
        }
        // else if((i*0.0693-6.93-3.465) < (0- 0.03465*j) && i<100) out_file << i << "\t" << j << "\t" << Jij[Jij.size()-i][Jij[i].size()-j] << "\n";
        // else  out_file << i << "\t" << j << "\t" << Jij[200+(200-i)][200+(200-j)] << "\n";
    }
   }
   out_file.close();
}


void read_in_dmi(std::string filename, std::vector < std::vector < double > > &Dx, std::vector < std::vector < double > > &Dy, std::vector < std::vector < double > > &Dz, bool reflect){
   std::ifstream ifile2(filename);
   std::string line;
    if(!ifile2.is_open()) {std::cerr << filename << " is not open" << std::endl; exit(1);}
   for(int i=0; i<(Dx.size()*Dx.size()); i++){
     getline(ifile2,line);
     std::stringstream liness(line.c_str());
     double ii;
     double ij;
     liness >> ii >> ij >> Dx[ii+100][ij+100] >> Dy[ii+100][ij+100] >> Dz[ii+100][ij+100];
   //   outfile <<atom_id[i] << "\t" << x_in[i] << '\t' << y_in[i] << "\t" << z_in[i] << "\t" << S_in[i] << "\t" << temp << "\t" << temp << std::endl;
   // std::cout <<ii << "\t" << ij << "\t" <<  Dx[ii][ij]<< "\t" <<  Dy[ii][ij]<< "\t" <<  Dz[ii][ij] << std::endl;
 }
   ifile2.close();
 if(reflect){
   std::ofstream out_file(filename + "_out.txt");
   if(!out_file.is_open()) {std::cerr << filename + "_out.txt is not open" << std::endl; exit(1);}
   for(int i = 0; i < Dx.size(); i++) {
    for(int j = 0; j < Dx[i].size(); j++) {
        if( i > 100 ) {

            int dx = 200-i;//r*cos(theta+2*M_PI/3.0);//-y*sin(rad_120);
            int dy = j;//r*sin(theta+2*M_PI/3.0);//+x*sin(rad_120);
            Dx[dx][dy] = Dx[i][j];
            Dy[dx][dy] = -Dy[i][j];
            // Dy[i][j] *= -1;
            Dz[dx][dy] = Dz[i][j];
            out_file << i << "\t" << j << "\t" << Dx[i][j] << "\t" << Dy[i][j] << "\t" << Dz[i][j] << "\n";
            out_file << dx << "\t" << dy << "\t" << Dx[dx][dy] << "\t" << Dy[dx][dy] << "\t" << Dz[dx][dy] <<  "\n";
        }
    }
   }
   out_file.close();
 }
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

void read_in_inter_exchanges(std::string filename, std::vector<std::vector<double> > &Eij) {
    std::ifstream ifile2(filename);
    std::string line;
    if(!ifile2.is_open()) {std::cerr << filename << " is not open" << std::endl; exit(1);}
    
    for(int i=0; i<Eij.size(); i++){
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

        liness >> r >> dx >> dy >> dz >> J >> Dx >> Dy >> Dz;
        Eij[i][0] = dx;
        Eij[i][1] = dy;
        Eij[i][2] = J*J_constant;
        Eij[i][3] = Dx*J_constant;
        Eij[i][4] = Dy*J_constant;
        Eij[i][5] = Dz*J_constant;
    }
}

void read_in_intra_exchanges(std::string filename, std::vector<std::vector<std::vector<std::vector<double> > > > &Eij_1NN, \
                                                   std::vector<std::vector<std::vector<std::vector<double> > > > &Eij_2NN, \
                                                   std::vector<std::vector<std::vector<std::vector<double> > > > &Eij_3NN ) {

    std::ifstream ifile2(filename);
    std::string line;
    if(!ifile2.is_open()) {std::cerr << filename << " is not open" << std::endl; exit(1);}
    
    for(int i=0; i<Eij_1NN.size(); i++){
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
            Eij_1NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eij_1NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eij_1NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eij_1NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
             std::cout  <<  Eij_1NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eij_1NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eij_1NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eij_1NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
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
            Eij_2NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eij_2NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eij_2NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eij_2NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
             std::cout  <<  Eij_2NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eij_2NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eij_2NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eij_2NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
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
            Eij_3NN.at(int_x).at(int_y).at(theta_i)[0] = J*J_constant;
            Eij_3NN.at(int_x).at(int_y).at(theta_i)[1] = Dx*J_constant;
            Eij_3NN.at(int_x).at(int_y).at(theta_i)[2] = Dy*J_constant;
            Eij_3NN.at(int_x).at(int_y).at(theta_i)[3] = Dz*J_constant;
            std::cout  <<  Eij_3NN.at(int_x).at(int_y).at(theta_i)[0] << ", " << Eij_3NN.at(int_x).at(int_y).at(theta_i)[1] << ", " << Eij_3NN.at(int_x).at(int_y).at(theta_i)[2] << ", " << Eij_3NN.at(int_x).at(int_y).at(theta_i)[3] << std::endl;
        }
    }
}

void read_in_ucf(std::ifstream &ucf_file) {
    std::cout << "ucf file add on has been selected. reading secondary file..." << std::endl;
    // std::cout << "Warning. DMI rotation has not been designed yet." << std::endl;
    // std::cout << "Warning. Lattice basis unity has not been ensured yet." << std::endl;
    if(!ucf_file.is_open()) {std::cerr << "ucf_file is not yet open" << std::endl; exit(1);}
    std::string line;
      // keep record of current line
    unsigned int line_counter=0;
    unsigned int line_id=0;
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
        double ucf_cx;
        double ucf_cy;
        double ucf_cz;
        // non-comment line found - check for line number
        switch(line_id){
        case 0:
            iss >> ucf_cx >> ucf_cy >> ucf_cz;
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
            std::cout << "Reading in " << num_uc_atoms << " atoms" << std::endl;
            // resize unit_cell.atom array if within allowable bounds
            if( (num_uc_atoms >0) && (num_uc_atoms <= 1000000)) all_m_atoms.reserve(num_uc_atoms);
            else {
                // terminaltextcolor(RED);
                std::cerr << "Error! Requested number of atoms " << num_uc_atoms << " on line " << line_counter
                << " of unit cell input file is outside of valid range 1-1,000,000. Exiting" << std::endl; exit(1);
                // terminaltextcolor(WHITE);
            }

            // std::cout << "\nProcessing data for " << all_m_atoms.size() << " atoms..." << std::flush;
            // zlog << zTs() << "\t" << "Processing data for " << unit_cell.atom.size() << " unit cell atoms..." << std::endl;

        // loop over all atoms and read into class
        for(unsigned int i = 0; i < num_uc_atoms; i++){

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
            // id += all_m_atoms.size();
            if(cx>=0.0 && cx <=1.0) new_atom.x=cx*ucf_cx;
            else{
                // terminaltextcolor(RED);
                std::cerr << "Error! atom x-coordinate for atom " << id << " on line " << line_counter
                                << " of unit cell input file is outside of valid range 0.0-1.0. Exiting" << std::endl;
                // terminaltextcolor(WHITE);
                // zlog << zTs() << "Error! atom x-coordinate for atom " << id << " on line " << line_counter
                            //  << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
                exit(1);
            }
            if(cy>=0.0 && cy <=1.0) new_atom.y=cy*ucf_cy;
            else{
                // terminaltextcolor(RED);
                std::cerr << "Error! atom y-coordinate for atom " << id << " on line " << line_counter
                                << " of unit cell input file is outside of valid range 0.0-1.0. Exiting" << std::endl;
                // terminaltextcolor(WHITE);
                // zlog << zTs() << "Error! atom y-coordinate for atom " << id << " on line " << line_counter
                // 			     << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
                exit(1);
            }
            if(cz>=0.0 && cz <=1.0) new_atom.z=cz*ucf_cz;
            else{
                // terminaltextcolor(RED);
                std::cerr << "Error! atom z-coordinate for atom " << id << " on line " << line_counter
                << " of unit cell input file is outside of valid range 0.0-1.0. Exiting" << std::endl;
                // terminaltextcolor(WHITE);
                // zlog << zTs() << "Error! atom z-coordinate for atom " << id << " on line " << line_counter
                // 				  << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
                exit(1);
            }
            new_atom.unit_y = floor((new_atom.y +0.0000001)/ a1y);
            // changex += dy_cell*std::abs(a1x);
            new_atom.unit_x = floor((new_atom.x +0.0000001)/ a0x);
            new_atom.S = 5;
            if(cz > 0.5) new_atom.S = 6;
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
            double eVtoJ = 1.602176634e-19;
            double J_constant = 1.0*eVtoJ/1000.0; //1 meV
            // get number of exchange types
            iss >> num_interactions >> exchange_type_string;
            std::cout << "Reading in " << num_interactions  << " interactions " << std::endl;
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
            
                line_counter++;
                // check for sane input
                // id += num_interactions;
                // iatom += all_m_atoms.size();
                // jatom += all_m_atoms.size();
                // check for long range interactions
                // if(dx*dx+dy*dy+dz*dz > r2) continue;
                interaction ucf_interaction;
                            //xx                     xy-> Dz             xz -> -Dy
                int_iss >> ucf_interaction.xx >> ucf_interaction.xy >> ucf_interaction.xz;
                            //yx -> -Dz              yy                  yz -> Dx
                int_iss >> ucf_interaction.yx >> ucf_interaction.yy >> ucf_interaction.yz;
                            //zx -> Dy               yz -> -Dx           zz
                int_iss >> ucf_interaction.zx >> ucf_interaction.zy >> ucf_interaction.zz;

                spin atom_i = all_m_atoms[iatom];
                spin atom_j = all_m_atoms[jatom];
                config_energy.at(atom_i.unit_x).at(atom_i.unit_y).at((atom_i.S-1)*5+0) += 1.0;
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
      // std::cout << "Writing data to file..." << std::flush;
      std::ofstream config_output;
      config_output.open("config_energy.txt");
      for(int i = 0; i < number_of_unit_cells_x; i++) {
         for(int j = 0; j < number_of_unit_cells_y; j++){
            double bottom_occ = config_energy[i][j][0];
            double top_occ = config_energy[i][j][1];
            config_output << i << ", " << j << ", " << bottom_occ<< ", " << top_occ;
            for(int k = 2; k < config_energy[i][j].size(); k++) config_output << ", " << config_energy[i][j][k]; 
            config_output << "\n";
         }
      }
      config_output.close();

      std::ofstream interaction_counts;
      interaction_counts.open("interaction_counts.txt");
      for(int i = 0; i < all_m_atoms.size(); i++){
         interaction_counts << all_m_atoms[i].S  << ", " <<  all_m_atoms[i].l_id << ", " << all_m_atoms[i].inter1 << ", " << all_m_atoms[i].inter2 << ", " << all_m_atoms[i].inter3 \
                                                   << ", " << all_m_atoms[i].intra1 << ", " << all_m_atoms[i].intra2 << ", " << all_m_atoms[i].intra3 <<"\n";
      }
      interaction_counts.close();
      // outfile4 << ss.str();
    //   timer.stop();
      // std::cout << "done!  << std::endl;
    //   std::cout << number_of_interactions << " [completed] [" << timer.elapsed_time() << " s]" << std::endl;
      return;
}

