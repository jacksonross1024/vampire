#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string> 
#include <iomanip>


//to compile: 'g++ Topological-charge.cpp'
// executable: a.out
int main(){ 
//  std::cout << M_PIl <<std::endl ;  
  double min[3] = {-2001,-2001,-500};
  double max[3] = { 2001, 2001, 500};


  std::ofstream ofile;
  ofile.open("Top-charge.txt");
for(int f = 0; f < 196; f++) { //which cell files to look at 


  // Determine cell file name
  std::stringstream cell_file_sstr;
  cell_file_sstr << "cells-";
  cell_file_sstr << std::setfill('0') << std::setw(8) << f;
  cell_file_sstr << ".txt";
  std::string cell_file_name = cell_file_sstr.str();

  // output informative message to user
//  std::cout << "   Reading cell file " << cell_file_name << "..." << std::flush;

 
  std::ifstream ifile(cell_file_name); // very-small
  if(!ifile.is_open()) exit(1);


  int y_max = 1480; //
  int y_min = 20;
  int x_max = 1480;
  int x_min = 20;
  int cell_size = 10; //A
  int number_of_lines = (y_max-y_min)*(x_max-x_min)/cell_size/cell_size; // very-small

  //std::cout << "Number of lines in text file: " << number_of_lines << std::endl;  
  int j = 0,m=0,n=0; 
  std::cout << "test0" <<std::endl ;  
  double xyz [number_of_lines][3] ; 
  double S [number_of_lines][3] ;
  double S_up [number_of_lines][3] ;
  double S_down [number_of_lines][3] ;

  // std::cout << "test1 "<< number_of_lines << std::endl ;  
  int line_count = 0;
  while(!ifile.eof()){
    line_count++;
    std::string line;
    getline(ifile, line);
    double nope1,nope2,x,y,z,sx,sy,sz, length;
    std::stringstream if_ss(line);
    
    if_ss >> x >> y >> z >> nope1 >> nope1 >> nope1 >> nope1 >> nope1 \
                         >> nope1 >> nope1 >> nope1 >> nope1 >> nope1 \
                         >> nope1 >> nope1 >> nope1 >> nope1 >> nope1 \
                         >> nope1 >> nope1 >> nope1 >> nope1 >> nope1 \
                         //>> nope1 >> nope1 >> nope1 >> nope1 >> nope1 
                         >> sx >> sy >> sz >> length;

    if(z != 0.0) continue;

    // std::cout << sx << ", " << sy << ", " << sz << ", " << length << std::endl;

    if(x < x_min) continue;
    else if (x >= x_max) continue;

    if(y < y_min) continue;
    else if (y >= y_max) continue;
    xyz[j][0] = x;
    xyz[j][1] = y;
    xyz[j][2] = z;
    S[j][0] = sx*length;
    S[j][1] = sy*length;
    S[j][2] = sz*length;
    j++;
    if (sz>0){
    	S_up[m][0] = sx;
      	S_up[m][1] = sy;
    	S_up[m][2] = sz;
    	m++;
    }
    else{
    	S_down[n][0] = sx;
      	S_down[n][1] = sy;
    	S_down[n][2] = sz;
    	n++;
    }
    
  }




  //std::cout << "test2"<<std::endl ; 
  double OMEGA_tot = 0;
/////////////////////////////////////////////////////////////////
// REMEMBER TO SORT YOUR DATASET WITH RESPECT TO THE x COORDINATE
/////////////////////////////////////////////////////////////////

  for(int i=0;i<number_of_lines; i++){
    for(int j=i+1;j<i+(y_max-y_min); j++){
      if(j >= number_of_lines) break;
      double x1 = xyz[i][0];
      double x2 = xyz[j][0];
      double y1 = xyz[i][1];
      double y2 = xyz[j][1];
      //double d12 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)); //distance S1-S2
      if( fabs(x1-x2) < cell_size+0.01 && fabs(x1-x2) > cell_size-0.01 && (y1-y2) == 0){
        for(int k=i+1; k<i+(y_max-y_min); k++){
          if(k >= number_of_lines) break;
	        double x4 = xyz[k][0];
          double y4 = xyz[k][1];
          double d14 = (x1-x4)*(x1-x4) + (y1-y4)*(y1-y4); //distance S1-S4
	        double d24 = (x2-x4)*(x2-x4) + (y2-y4)*(y2-y4);
          if(d24 < 2*(cell_size*cell_size+0.1) && d24 > 2*(cell_size*cell_size-0.1) && d14 > (cell_size*cell_size-0.1) && d14 < (cell_size*cell_size+0.1) && (y4-y1)>0){
            for(int l=i+1; l<j+(y_max-y_min); l++){
              if(l >= number_of_lines) break;
              double x3 = xyz[l][0];
              double y3 = xyz[l][1];
	            double d13 = (x1-x3)*(x1-x3) + (y1-y3)*(y1-y3);
	            double d23 = (x2-x3)*(x2-x3) + (y2-y3)*(y2-y3);
	            double d43 = (x4-x3)*(x4-x3) + (y4-y3)*(y4-y3);
              if(d23 < (cell_size*cell_size+0.1) && d23 > (cell_size*cell_size-0.1) && d43 > (cell_size*cell_size-0.1) && d43 < (cell_size*cell_size+0.1) && d13 > 2*(cell_size*cell_size-0.1) && d13 < 2*(cell_size*cell_size+0.1) ){
                double Num1 = S[i][0]*(S[j][1]*S[l][2] - S[j][2]*S[l][1]) + S[i][1]*(S[j][2]*S[l][0] - S[j][0]*S[l][2]) + S[i][2]*(S[j][0]*S[l][1] - S[j][1]*S[l][0]);
                double Dem1 = 1 + S[i][0]*S[j][0] + S[i][1]*S[j][1] + S[i][2]*S[j][2] + S[i][0]*S[l][0] + S[i][1]*S[l][1] + S[i][2]*S[l][2] + S[j][0]*S[l][0] + S[j][1]*S[l][1] + S[j][2]*S[l][2];
                double Num2 = S[i][0]*(S[l][1]*S[k][2] - S[l][2]*S[k][1]) + S[i][1]*(S[l][2]*S[k][0] - S[l][0]*S[k][2]) + S[i][2]*(S[l][0]*S[k][1] - S[l][1]*S[k][0]);
                double Dem2 = 1 + S[i][0]*S[l][0] + S[i][1]*S[l][1] + S[i][2]*S[l][2] + S[i][0]*S[k][0] + S[i][1]*S[k][1] + S[i][2]*S[k][2] + S[l][0]*S[k][0] + S[l][1]*S[k][1] + S[l][2]*S[k][2];
                double OMEGA_triangle1 = 2*atan(Num1/Dem1);
                double OMEGA_triangle2 = 2*atan(Num2/Dem2);
	              OMEGA_tot = OMEGA_triangle1 + OMEGA_triangle2 + OMEGA_tot; 
              }
            }
          }
        }
      } 
    }

  }
//////////////////////////////
//COMPUTING POLARITY
//////////////////////////////

  double Smax = 0;
  int p=0;
  for(int i=0;i<number_of_lines; i++){
    if(std::abs(S[i][2])>std::abs(Smax)){
	    Smax = S[i][2];
    }
  }
  if(Smax>0){
    p=1;
  }
  else{
    p=-1;
  }
  std::cout << "Q = "<< OMEGA_tot/(4*M_PIl)<<", p = "<< p << std::endl ;  
  ofile << f << ", " <<  OMEGA_tot/(4*M_PIl)<<", "<< p << std::endl ;  
}
  ofile.close();
  return 0;

}

