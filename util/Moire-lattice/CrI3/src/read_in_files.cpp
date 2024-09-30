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

// void read_in_from_ucf(std::string filename) {
//     std::ifstream ifile2(filename);
//     std::string line;
//     if(!ifile2.is_open()) {std::cerr << filename << " is not open" << std::endl; exit(1);}
    
//     unsigned int line_counter=0;
// 	unsigned int line_id=0;

//    std::string exchange_type_string; // string defining exchange type

//    // defaults for interaction list
//    unsigned int interaction_range = 1; // assume +-1 unit cell as default

// 	// Loop over all lines
// 	while (! ifile2.eof() ){
// 		line_counter++;
// 		// read in whole line
// 		std::string line;
// 		getline(ifile2,line);
// 		//std::cout << line.c_str() << std::endl;

// 		// ignore blank lines
// 		std::string empty="";
// 		if(line==empty) continue;

// 		// set character triggers
// 		const char* hash="#";	// Comment identifier

// 		bool has_hash=false;
// 		// Determine if line is a comment line
// 		for(unsigned int i=0;i<line.length();i++){
// 			char c=line.at(i);

// 			if(c== *hash){
// 					has_hash=true;
// 					break;
// 			}
// 		}
// 		// if hash character found then read next line
// 		if(has_hash==true) continue;

// 		// convert line to string stream
// 		std::istringstream iss(line,std::istringstream::in);

// 		// non-comment line found - check for line number
//         double unit_cell_dimensions[3] = {0.0,0.0,0.0};
//         double unit_cell_shape[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
//         double num_atoms = 0.0;
// 		switch(line_id){
// 			case 0:
// 				iss >> unit_cell_dimensions[0] >> unit_cell_dimensions[1] >> unit_cell_dimensions[2];
// 				break;
// 			case 1:
// 				iss >> unit_cell_shape[0][0] >> unit_cell_shape[0][1] >> unit_cell_shape[0][2];
// 				break;
// 			case 2:
// 				iss >> unit_cell_shape[1][0] >> unit_cell_shape[1][1] >> unit_cell_shape[1][2];
// 				break;
// 			case 3:
// 				iss >> unit_cell_shape[2][0] >> unit_cell_shape[2][1] >> unit_cell_shape[2][2];
// 				break;
// 			case 4:
// 				int num_uc_atoms;
// 				iss >> num_uc_atoms;
// 				//std::cout << "Reading in " << num_uc_atoms << " atoms" << std::endl;
// 				// resize unit_cell.atom array if within allowable bounds
// 				if( (num_uc_atoms >0) && (num_uc_atoms <= 1000000)) num_atoms = num_uc_atoms;
// 				else {
					
// 					std::cerr << "Error! Requested number of atoms " << num_uc_atoms << " on line " << line_counter
// 					<< " of unit cell input file " << filename.c_str() << " is outside of valid range 1-1,000,000. Exiting" << std::endl; 
					
// 				}

//            // loop over all atoms and read into class
//             for(unsigned int i = 0; i < num_atoms; i++){

// 					line_counter++;

// 					// declare safe temporaries for atom input
// 					int id=i;
// 					double cx=2.0, cy=2.0,cz=2.0; // coordinates - default will give an error
// 					int mat_id=0, lcat_id=0, hcat_id=0; // sensible defaults if omitted
// 					// get line
// 					std::string atom_line;
// 					getline(ifile2,atom_line);
// 					std::istringstream atom_iss(atom_line,std::istringstream::in);
//                     spin new_atom;

// 					atom_iss >> new_atom.id >> cx >> cy >> cz >> new_atom.S >> new_atom.l_id >> hcat_id;
// 					//std::cout << id << "\t" << cx << "\t" << cy << "\t" << cz<< "\t"  << mat_id << "\t" << lcat_id << "\t" << hcat_id << std::endl;
// 					//inputfile >> id >> cx >> cy >> cz >> mat_id >> lcat_id >> hcat_id;
// 					// now check for mostly sane input
// 					if(cx>=0.0 && cx <=1.0) new_atom.x=cx*unit_cell_dimensions[0];
// 					else{
					
// 						std::cerr << "Error! atom x-coordinate for atom " << id << " on line " << line_counter
// 									 << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
// 						exit(1);
// 					}
// 					if(cy>=0.0 && cy <=1.0) new_atom.y=cy*unit_cell_dimensions[1];
// 					else{
						
// 						std::cerr << "Error! atom y-coordinate for atom " << id << " on line " << line_counter
// 									 << " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
// 						exit(1);
// 					}
// 					if(cz>=0.0 && cz <=1.0) new_atom.z=cz*unit_cell_dimensions[2];
// 					else{
						
// 						std::cerr << "Error! atom z-coordinate for atom " << id << " on line " << line_counter
// 						<< " of unit cell input file " << filename.c_str() << " is outside of valid range 0.0-1.0. Exiting" << std::endl;
// 						exit(1);
// 					}
// 					all_m_atoms.push_back(new_atom);

//                     if(double_bilayer) {
//                         if(pristine_bilayer_type == "baab") {
//                             if(new_atom.z > twist_loction) {
//                                 spin new_atom_2;
//                                 new_atom_2.z = 1.0*unit_cell_dimensions[2];
//                                 if(new_atom.l_id) 
//                             }


//                         }
//                     }
                    
//                }
					
// 				break;
// }