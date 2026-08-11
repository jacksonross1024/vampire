//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

// program header
#include "vdc.hpp"

// openmp header
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
#endif

namespace vdc{

//------------------------------------------------------------------------------
// Centroid of sliced atoms (equal weight); used for --spin-cells camera look-at.
//------------------------------------------------------------------------------
static void sliced_atoms_centroid(double com[3]){

   com[0] = 0.0;
   com[1] = 0.0;
   com[2] = 0.0;
   const size_t n = vdc::sliced_atoms_list.size();
   if(n == 0){
      com[0] = vdc::system_centre[0];
      com[1] = vdc::system_centre[1];
      com[2] = vdc::system_centre[2];
      return;
   }
   for(size_t i = 0; i < n; i++){
      const unsigned int a = vdc::sliced_atoms_list[i];
      com[0] += vdc::coordinates[3*a + 0];
      com[1] += vdc::coordinates[3*a + 1];
      com[2] += vdc::coordinates[3*a + 2];
   }
   const double invn = 1.0 / static_cast<double>(n);
   com[0] *= invn;
   com[1] *= invn;
   com[2] *= invn;
}

//------------------------------------------------------------------------------
// Function to output cells.inc file compatible with povray
//------------------------------------------------------------------------------
void output_cells_inc_file(unsigned int spin_file_id){

   // Open Povray Include File
	std::stringstream incpov_file_sstr;
	incpov_file_sstr << "cells-";
	incpov_file_sstr << std::setfill('0') << std::setw(8) << spin_file_id;
	incpov_file_sstr << ".inc";
	std::string incpov_file = incpov_file_sstr.str();

   // output informative message to user
   if(vdc::verbose) std::cout << "   Writing povray file " << incpov_file << "..." << std::flush;

   // open incfile
   std::ofstream incfile;
   incfile.open(incpov_file.c_str());

   const unsigned int tmid = vdc::materials.size();

   //---------------------------------------------------------------------------
   // parallelise stream formatting for better performance
   // step 1: parallel formatted output to stringstream in memory
   // step 2: binary write of formatted text to output file (awesomely fast!)
   //---------------------------------------------------------------------------
   #pragma omp parallel num_threads(4)
   {

      std::stringstream otext;

      const double ux = vdc::cell_size[0]*0.5;
      const double uy = vdc::cell_size[1]*0.5;
      const double uz = vdc::cell_size[2]*0.5;

      // write to output text stream in parallel
      #pragma omp for
      for( unsigned int cell = 0; cell < total_cells; cell++){

         // get magnetization for colour contrast
         double mm = 1.0;//vdc::cell_magnetization[cell][tmid][3];
         double mx = vdc::cell_magnetization[cell][tmid][0]*mm;
         double my = vdc::cell_magnetization[cell][tmid][1]*mm;
         double mz = vdc::cell_magnetization[cell][tmid][2]*mm;

         // temporary thread private variables defining spin colours
         double red=0.0, green=0.0, blue=1.0;

         // calculate rgb components based on spin orientation
         vdc::rgb(mx, my, mz, red, green, blue);

         double cx = vdc::cell_coords[3*cell + 0] - vdc::system_centre[0];
         double cy = vdc::cell_coords[3*cell + 1] - vdc::system_centre[1];
         double cz = vdc::cell_coords[3*cell + 2] - vdc::system_centre[2];
         // if(cz >= -9.8) continue;
         // format text for povray file
         otext << "cell"<< "(" <<
                  cx << "," << cy << "," << cz << "," <<
                  cx-ux << "," << cy-uy << "," << cz-uz << "," <<
                  cx+ux << "," << cy+uy << "," << cz+uz << "," <<
                  mx << "," << my << "," << mz << "," <<
                  red << "," << green << "," << blue << ")\n";

      } // end of parallel for

      // force each thread to write to file in order
      #pragma omp critical
      incfile << otext.str();

   } // end of parallel region

   // flush data to include file and close
   incfile << std::flush;
   incfile.close();

   // output informative message to user
   if(vdc::verbose) std::cout << "done!" << std::endl;

   return;

}

//------------------------------------------------------------------------------
// Function to output spins.pov file compatible with povray
//------------------------------------------------------------------------------
void output_povray_cells_file(){

	std::ofstream pfile;
	pfile.open("cells.pov");

   // Calculate location of camera
   std::vector<double> dim = {vdc::system_size[0]+0.001, vdc::system_size[1]+0.001, vdc::system_size[2]+0.001};
   std::vector<double> vec(3);

   // direction camera looks (normalised)
   // technically only position if lookat is not (0,0,0)
	double size = sqrt(dim[0]*dim[0] + dim[1]*dim[1] + dim[2]*dim[2]);
   if (default_camera_pos){
      vec[0] = (1.0/dim[0]);
      vec[1] = (1.0/dim[1]);
      vec[2] = (1.0/dim[2]);
   }
   else { vec = vdc::camera_pos; }

   // normalise camera position vector
   double mag_vec = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
   vec[0]/=mag_vec;
   vec[1]/=mag_vec;
   vec[2]/=mag_vec;

   pfile << "//-------------------------------------------------------------------------" << std::endl;
   pfile << "// Povray file generated using vampire" << std::endl;
   pfile << "//-------------------------------------------------------------------------" << std::endl;
   pfile << "#version 3.5;"            << std::endl;
	pfile << "#include \"colors.inc\""  << std::endl;
	pfile << "#include \"metals.inc\""	<< std::endl;
	pfile << "#include \"screen.inc\""	<< std::endl;
   // look at position
	pfile << "#declare LX=" << vdc::camera_look_at[0]*dim[0]/2.0 << ";" << std::endl;
	pfile << "#declare LY=" << vdc::camera_look_at[1]*dim[1]/2.0 << ";" << std::endl;
	pfile << "#declare LZ=" << vdc::camera_look_at[2]*dim[2]/2.0 << ";" << std::endl;
   // camera position
	pfile << "#declare CX=" << size*vec[0]*6.0*vdc::camera_zoom << ";" << std::endl;
	pfile << "#declare CY=" << size*vec[1]*6.0*vdc::camera_zoom << ";" << std::endl;
	pfile << "#declare CZ=" << size*vec[2]*6.0*vdc::camera_zoom << ";" << std::endl;
	pfile << "#declare ref=0.05;" << std::endl;
	pfile << "global_settings { assumed_gamma 2.0 }" << std::endl;
   // background colour
	pfile << "background { color " << vdc::background_colour << " }" << std::endl;

	pfile << "Set_Camera(<CX,CY,CZ>, <LX,LY,LZ>, 15)" << std::endl;
	pfile << "Set_Camera_Aspect(4,3)"  << std::endl;
	pfile << "Set_Camera_Sky(<0,0,1>)" << std::endl;
	pfile << "light_source { <2*CX, 2*CY, 2*CZ> color White}" << std::endl;

   pfile << "#declare Initial_Frame = " << vdc::start_file_id << ";" << std::endl;
   pfile << "#declare Final_Frame = "   << vdc::final_file_id << ";" << std::endl;

   // sscale affects the spin arrow
   pfile << "#declare mscale" << "=" << 2.0 << ";" << std::endl;
   // rscale affects sphere(atom) size
   // cscale affects cube size
   pfile << "#declare spx = 0.1; // spacing between cells" << std::endl;
   pfile << "#declare spy = 0.1;" << std::endl;
   pfile << "#declare spz = 0.1;" << std::endl;
   pfile << "#declare cones"   << " = true;" << std::endl;
   pfile << "#declare arrows"  << " = false;" << std::endl;
   pfile << "#declare cubes"   << " = true;" << std::endl;
   pfile << "#declare mcolors" << " = true;" << std::endl;
   pfile << "#declare mcolor"  << " = pigment {color rgb < 0.1 0.1 0.1 >};" << std::endl;
   pfile << "#macro cell"<< "(cx,cy,cz,sx,sy,sz,ex,ey,ez,mx,my,mz,cr,cg,cb)" << std::endl;
   pfile << "union{"      << std::endl;
   pfile << "#if(cubes) box {<sx+spx*0.5,sy+spy*0.5,sz+spz*0.5>,<ex-spx*0.5,ey-spy*0.5,ez-spz*0.5>} #end" << std::endl;
   pfile << "#if(cones) cone {<cx+0.5*mx*mscale,cy+0.5*my*mscale,cz+0.5*mz*mscale>,0.0 <cx-0.5*mx*mscale,cy-0.5*my*mscale,cz-0.5*mz*mscale>,mscale*0.5} #end" << std::endl;
   pfile << "#if(arrows) cylinder {<cx+mx*0.5*mscale,cy+my*0.5*mscale,cz+mz*0.5*mscale>," <<
                                  "<cx-mx*0.5*mscale,cy-my*0.5*mscale,cz-mz*0.5*mscale>,mscale*0.12}";
   pfile << "            cone {<cx+mx*0.5*1.6*mscale, cy+my*0.5*1.6*mscale, cz+mz*0.5*1.6*mscale>,mscale*0.0" <<
                              "<cx+mx*0.5*mscale,     cy+my*0.5*mscale,     cz+mz*0.5*mscale    >,mscale*0.2} #end" << std::endl;
   pfile << "#if(mcolors) texture { pigment {color rgb <cr cg cb>}finish {reflection {ref} diffuse 1 ambient 0}}" << std::endl;
   pfile << "#else texture { mcolor finish {reflection {ref} diffuse 1 ambient 0}} #end" << std::endl;
   pfile << "}"    << std::endl;
   pfile << "#end" << std::endl;

   // frame specific povray output
	pfile << "#include concat(\"cells-\", str(frame_number, -8, 0) \".inc\")" << std::endl;

   // close output file
	pfile.close();

   //---------------------------------------------------------------------------
   // output povray ini file for rendering all files by default
   //---------------------------------------------------------------------------
   std::ofstream pifile;
	pifile.open("cells.ini");

   pifile << "Input_File_Name = \"cells.pov\"" << std::endl;
   pifile << "Width = 800"                     << std::endl;
   pifile << "Height = 600"                    << std::endl;
   pifile << "Antialias = On"                  << std::endl;
   pifile << "Antialias_Threshold = 0.3"       << std::endl;
   pifile << "Output_File_Type = N"            << std::endl;
   pifile << "Initial_Frame = " << vdc::start_file_id << std::endl;
   pifile << "Final_Frame = "   << vdc::final_file_id << std::endl;

   pifile.close();

   return;

}

//------------------------------------------------------------------------------
// Per-frame spins-XXXXXXXX.inc for --spin-cells (cell magnetization as spinmc objects)
//------------------------------------------------------------------------------
void output_spin_cells_inc_file(unsigned int spin_file_id){

   std::stringstream incpov_file_sstr;
   incpov_file_sstr << "spins-";
   incpov_file_sstr << std::setfill('0') << std::setw(8) << spin_file_id;
   incpov_file_sstr << ".inc";
   std::string incpov_file = incpov_file_sstr.str();

   if(vdc::verbose) std::cout << "   Writing povray file " << incpov_file << "..." << std::flush;

   std::ofstream incfile;
   incfile.open(incpov_file.c_str());

   const unsigned int tmid = vdc::materials.size();
   const double max_cell_dim = std::max(vdc::cell_size[0], std::max(vdc::cell_size[1], vdc::cell_size[2]));
   const double arrow_scale = (vdc::arrow_sizes.empty() ? 2.0 : vdc::arrow_sizes[0]);
   const double ss = 0.40 * max_cell_dim * (arrow_scale / 2.0);

   #pragma omp parallel num_threads(4)
   {

      std::stringstream otext;

      #pragma omp for
      for(unsigned int cell = 0; cell < total_cells; cell++){

         double mx = vdc::cell_magnetization[cell][tmid][0];
         double my = vdc::cell_magnetization[cell][tmid][1];
         double mz = vdc::cell_magnetization[cell][tmid][2];

         // Same lab-frame unit vector for PoVRAY arrow geometry and vdc::rgb() (vector-z/x, colourmap, 3D)
         double sx = mx, sy = my, sz = mz;
         const double nm = std::sqrt(sx*sx + sy*sy + sz*sz);
         if(nm > 1.0e-12){
            sx /= nm;
            sy /= nm;
            sz /= nm;
         }
         else{
            sx = 0.0;
            sy = 0.0;
            sz = 1.0;
         }

         double red = 0.0, green = 0.0, blue = 1.0;
         vdc::rgb(sx, sy, sz, red, green, blue);

         double cx = vdc::cell_coords[3*cell + 0] - vdc::system_centre[0];
         double cy = vdc::cell_coords[3*cell + 1] - vdc::system_centre[1];
         double cz = vdc::cell_coords[3*cell + 2] - vdc::system_centre[2];

         otext << "spinmc(" <<
                  cx << "," << cy << "," << cz << "," <<
                  sx << "," << sy << "," << sz << "," <<
                  ss << "," <<
                  red << "," << green << "," << blue << ")\n";

      }

      #pragma omp critical
      incfile << otext.str();

   }

   incfile << std::flush;
   incfile.close();

   if(vdc::verbose) std::cout << "done!" << std::endl;

   return;

}

void output_povray_spin_cells_file(){

	std::ofstream pfile;
	pfile.open("spins.pov");

   const double max_cell_dim = std::max(vdc::cell_size[0], std::max(vdc::cell_size[1], vdc::cell_size[2]));

   // Calculate location of camera
   std::vector<double> dim = {vdc::system_size[0]+0.001, vdc::system_size[1]+0.001, vdc::system_size[2]+0.001};
   std::vector<double> vec(3);

   // direction camera looks (normalised)
   // technically only position if lookat is not (0,0,0)
	double size = sqrt(dim[0]*dim[0] + dim[1]*dim[1] + dim[2]*dim[2]);
   if (default_camera_pos){
      vec[0] = (1.0/dim[0]);
      vec[1] = (1.0/dim[1]);
      vec[2] = (1.0/dim[2]);
   }
   else { vec = vdc::camera_pos; }

   // normalise camera position vector
   double mag_vec = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
   vec[0]/=mag_vec;
   vec[1]/=mag_vec;
   vec[2]/=mag_vec;

   const double camr = size*6.0*vdc::camera_zoom;
   const double camp = 180.0*acos(vec[1]/sqrt(vec[0]*vec[0]+vec[1]*vec[1]))/M_PI;
   const double camt = 180.0*acos(vec[2]/sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]))/M_PI;

   // Scene coords subtract bbox centre (system_centre); aim default look-at at sliced centroid (COM for equal masses)
   double com[3];
   sliced_atoms_centroid(com);
   const double look_com_x = com[0] - vdc::system_centre[0];
   const double look_com_y = com[1] - vdc::system_centre[1];
   const double look_com_z = com[2] - vdc::system_centre[2];

   pfile << "//-------------------------------------------------------------------------" << std::endl;
   pfile << "// Povray file generated using vampire" << std::endl;
   pfile << "//-------------------------------------------------------------------------" << std::endl;
   pfile << "#version 3.5;"            << std::endl;
	pfile << "#include \"colors.inc\""  << std::endl;
	pfile << "#include \"metals.inc\""	<< std::endl;
	pfile << "#include \"screen.inc\""	<< std::endl;
   // look at position
   pfile << "//-----------------------------------------------------------------------------------" << std::endl;
   pfile << "// Camera parameters" << std::endl;
   pfile << "//-----------------------------------------------------------------------------------" << std::endl;
	pfile << "#declare LX = " << look_com_x + vdc::camera_look_at[0]*dim[0]/2.0 << "; // look-at: sliced centroid + user camera-look-at offset" << std::endl;
	pfile << "#declare LY = " << look_com_y + vdc::camera_look_at[1]*dim[1]/2.0 << ";" << std::endl;
	pfile << "#declare LZ = " << look_com_z + vdc::camera_look_at[2]*dim[2]/2.0 << ";" << std::endl;
   pfile << "// camera location" << std::endl;
   pfile << "#declare cam_theta  = " << camt << "; // angle from z in degrees" << std::endl;
   pfile << "#declare cam_phi    = " << camp << "; // angle from x in degrees" << std::endl;
   pfile << "#declare cam_radius = " << camr << "; // distance from origin" << std::endl;
   pfile << "#declare CX = cam_radius * cos(cam_phi*pi/180.0) * sin(cam_theta*pi/180.0);" << std::endl;
   pfile << "#declare CY = cam_radius * sin(cam_phi*pi/180.0) * sin(cam_theta*pi/180.0);" << std::endl;
   pfile << "#declare CZ = cam_radius * cos(cam_theta*pi/180.0);" << std::endl;
   // camera position
	//pfile << "#declare CX = " << size*vec[0]*6.0*vdc::camera_zoom << "; // camera location" << std::endl;
	//pfile << "#declare CY = " << size*vec[1]*6.0*vdc::camera_zoom << ";" << std::endl;
	//pfile << "#declare CZ = " << size*vec[2]*6.0*vdc::camera_zoom << ";" << std::endl;
   pfile << "Set_Camera(<CX,CY,CZ>, <LX,LY,LZ>, 15)" << std::endl;
	pfile << "Set_Camera_Aspect(4,3)"  << std::endl;
	pfile << "Set_Camera_Sky(<0,0,1>)" << std::endl;

   pfile << "//-----------------------------------------------------------------------------------" << std::endl;
   pfile << "// Global constants and appearance" << std::endl;
   pfile << "//-----------------------------------------------------------------------------------" << std::endl;
	pfile << "global_settings { assumed_gamma 2.0 }" << std::endl;
   // background colour
	pfile << "background { color " << vdc::background_colour << " } // background  colour" << std::endl;


	pfile << "light_source { <2*CX, 2*CY, 2*CZ> color White} // lights" << std::endl;
   pfile << "#declare Initial_Frame = " << vdc::start_file_id << ";" << std::endl;
   pfile << "#declare Final_Frame = "   << vdc::final_file_id << ";\n" << std::endl;

   pfile << "//------------------------------------------" << std::endl;
   pfile << "// Global settings for appearance" << std::endl;
   pfile << "//------------------------------------------\n" << std::endl;

   pfile << "#declare ref = 0.05; // reflection level of objects" << std::endl;
   pfile << "#declare dif = 1.0; // diffuse level of objects" << std::endl;
   pfile << "#declare amb = 0.1; // ambient level of objects\n" << std::endl;

   pfile << "// spin scale (default 2.0)" << std::endl;
   pfile << "#declare global_spin_scale = true;" << std::endl;
   pfile << "#declare sscale = 2.0;\n" << std::endl;

   pfile << "// radius / cube scale scaled to microcell size (larger than atomistic --povray)" << std::endl;
   pfile << "#declare global_radius_scale = true;" << std::endl;
   pfile << "#declare rscale = " << 0.30 * max_cell_dim << ";\n" << std::endl;

   pfile << "#declare global_cube_scale = true;" << std::endl;
   pfile << "#declare cscale = " << max_cell_dim << ";\n" << std::endl;

   pfile << "// global objects" << std::endl;
   pfile << "#declare global_cones   = false;" << std::endl;
   pfile << "#declare global_arrows  = true;" << std::endl;
   pfile << "#declare global_spheres = true;" << std::endl;
   pfile << "#declare global_cubes   = false;\n" << std::endl;

   pfile << "// non-magnetic element colours" << std::endl;
   pfile << "#declare global_non_magnet_colour = true;" << std::endl;
   pfile << "#declare nmcolour = pigment {color rgb < 0.2 0.2 0.2 >};\n" << std::endl;

   pfile << "//-----------------------------------------------------------------------------------" << std::endl;
   pfile << "// Colour macro" << std::endl;
   pfile << "// The default is vdc generated and takes rgb values, but can be overridden to" << std::endl;
   pfile << "// apply non-linear colour scales, colour maps etc." << std::endl;
   pfile << "//-----------------------------------------------------------------------------------" << std::endl;
   pfile << "#macro spinrgb(sx, sy, sz, cr, cg, cb)" << std::endl;
   pfile << "   pigment {color rgb <cr cg cb>}" << std::endl;
   pfile << "#end\n" << std::endl;

   pfile << "// alternative colour schemes" << std::endl;
   pfile << "//#include \"util/povray_colours/jet.inc\"" << std::endl;
   pfile << "//#include \"util/povray_colours/purple_white.inc\"" << std::endl;
   pfile << "//#include \"util/povray_colours/blue_gold.inc\"" << std::endl;
   pfile << "//#include \"util/povray_colours/color_wheel.inc\"\n" << std::endl;

   

   pfile << "#if(global_spin_scale) #declare sscalec = sscale;" << std::endl;
   pfile << "#else #declare sscalec = " << 2.0 << "; #end" << std::endl;
   pfile << "#if(global_radius_scale) #declare rscalec = rscale;" << std::endl;
   pfile << "#else #declare rscalec = " << 1.2<< "; #end" << std::endl;
   pfile << "#if(global_cube_scale) #declare cscalec = cscale;" << std::endl;
   pfile << "#else #declare cscalec = 3.54; #end" << std::endl;

   pfile << "#if(global_cones)   #declare conesc = global_cones;   #else #declare conesc = false; #end" << std::endl;
   pfile << "#if(global_arrows)  #declare arrowsc = global_arrows;  #else #declare arrowsc = true;  #end" << std::endl;
   pfile << "#if(global_spheres) #declare spheresc = global_spheres; #else #declare spheresc = true;  #end" << std::endl;
   pfile << "#if(global_cubes)   #declare cubesc = global_cubes;   #else #declare cubesc = false; #end" << std::endl;
   pfile << "#declare spincolorsc = true; // enable colours defined in vdc" << std::endl;
   pfile << "#declare spincolorc  = pigment {color rgb < 0.1 0.1 0.1 >};" << std::endl;
   pfile << "//-------------------------------------" << std::endl;
   pfile << "#macro spinmc(cx,cy,cz,sx,sy,sz,ss, cr,cg,cb)" << std::endl;
   pfile << "#declare sscalec = ss;"      << std::endl;
   pfile << "union{"      << std::endl;
   pfile << "   #if(spheresc) sphere {<cx,cy,cz>,0.5*rscalec} #end" << std::endl;
   pfile << "   #if(cubesc) box {<cx-cscalec*0.5,cy-cscalec*0.5,cz-cscalec*0.5>,<cx+cscalec*0.5,cy+cscalec*0.5,cz+cscalec*0.5>} #end" << std::endl;
   pfile << "   #if(conesc) cone {<cx+0.5*sx*sscalec,cy+0.5*sy*sscalec,cz+0.5*sz*sscalec>,0.0 <cx-0.5*sx*sscalec,cy-0.5*sy*sscalec,cz-0.5*sz*sscalec>,sscalec*0.5} #end" << std::endl;
   pfile << "   #if(arrowsc)" << std::endl;
   pfile << "      cylinder {<cx+sx*0.5*sscalec,    cy+sy*0.5*sscalec,    cz+sz*0.5*sscalec>" << std::endl;
   pfile << "                <cx-sx*0.5*sscalec,    cy-sy*0.5*sscalec,    cz-sz*0.5*sscalec>,sscalec*0.12}" << std::endl;
   pfile << "      cone     {<cx+sx*0.5*1.6*sscalec,cy+sy*0.5*1.6*sscalec,cz+sz*0.5*1.6*sscalec>,sscalec*0.0" << std::endl;
   pfile << "                <cx+sx*0.5*sscalec,    cy+sy*0.5*sscalec,    cz+sz*0.5*sscalec    >,sscalec*0.2}" << std::endl;
   pfile << "   #end" << std::endl;
   pfile << "   #if(spincolorsc) texture { spinrgb(sx,sy,sz,cr,cg,cb) finish {reflection ref diffuse dif ambient amb } }" << std::endl;
   pfile << "   #else texture { spincolorc finish {reflection ref diffuse dif ambient amb }} #end" << std::endl;
   pfile << "}"    << std::endl;
   pfile << "#end\n" << std::endl;
   // frame specific povray output
   pfile << "//----------------------------------------------------------------" << std::endl;
   pfile << "// Include spin data" << std::endl;
   pfile << "//----------------------------------------------------------------" << std::endl;
   pfile << "#include concat(\"spins-\", str(frame_number, -8, 0) \".inc\")" << std::endl;

   // close output file
	pfile.close();

   //---------------------------------------------------------------------------
   // output povray ini file for rendering all files by default
   //---------------------------------------------------------------------------
   std::ofstream pifile;
	pifile.open("spins.ini");

   pifile << "Input_File_Name = \"spins.pov\"" << std::endl;
   pifile << "Width = 800"                     << std::endl;
   pifile << "Height = 600"                    << std::endl;
   pifile << "Antialias = On"                  << std::endl;
   pifile << "Antialias_Threshold = 0.3"       << std::endl;
   pifile << "Output_File_Type = N"            << std::endl;
   pifile << "Initial_Frame = " << vdc::start_file_id << std::endl;
   pifile << "Final_Frame = "   << vdc::final_file_id << std::endl;

   pifile.close();

   return;

}

//------------------------------------------------------------------------------
// Rewrite spins.ini after snapshots processed so animation range matches .inc files
//------------------------------------------------------------------------------
void rewrite_spins_ini_spin_cells(){

   std::ofstream pifile;
   pifile.open("spins.ini");

   pifile << "Input_File_Name = \"spins.pov\"" << std::endl;
   pifile << "Width = 800"                     << std::endl;
   pifile << "Height = 600"                    << std::endl;
   pifile << "Antialias = On"                  << std::endl;
   pifile << "Antialias_Threshold = 0.3"       << std::endl;
   pifile << "Output_File_Type = N"            << std::endl;
   pifile << "Initial_Frame = " << vdc::start_file_id << std::endl;
   pifile << "Final_Frame = "   << vdc::final_file_id << std::endl;

   pifile.close();

   return;
}

}

