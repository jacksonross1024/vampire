//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans, Daniel Meilak 2017-2019. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iostream>
#include <sstream>
#include <cmath>
#include <cctype>
#include <vector>
#include <algorithm>
#include <fstream>

// program header
#include "vdc.hpp"

namespace vdc{

// forward function declarations
// void extract_vector( std::string arg_string, std::vector<double>& arg_vector );
// void extract_materials( std::string arg_string, std::vector<int>& arg_vector);
// void extract_slice_param( std::string arg_string, std::vector<double>& arg_vector, int number_of_param);
void check_arg( int &arg, int argc, char* argv[], std::string &temp_str, std::string error_output );
void init_vector_y();

//------------------------------------------------------------------------------
// Command line parsing function
//------------------------------------------------------------------------------
void command( int argc, char* argv[] ){
// Command line options for utility to be implemented:
//    --xyz - generate xyz file
//    --povray - generate povray files
//    --vector - generate raw xyz vector data
//    --vtk - generate vtk files
//    --cells - collate data into cells
//    --binary - native-endian binary dump of data files
//    --omp-threads N - OpenMP threads for --binary write path
//    --cell-size = x - define cell size
//    --stray-field - FFT stray field above the sample (needs FFTW3)
//    --sf-height = h - observation height in nm (default 30)
//    --sf-pad = n - FFT zero-padding factor (1 = periodic, default 2)
//    --povray-cells - generate renderable cell positions for povray
//    --verbose [= true, false] - set verbose output to screen
//    --colours = default, rwb [red-white-blue], oyb [orange-yellow-blue], jet, gs [grey-scale], cw [colour-wheel], <filename>
//    --objects = spins, cones, spheres, cubes
//    --slice = x,x,y,y,z,z
//    --multiscale = gradient, material, region
//    --track  -generate povray track vectors
   // temporary string for storing command line argument
   std::string temp_str;

   for (int arg = 1; arg < argc; arg++){

      // read prefix
      std::string sw=argv[arg];

      // convert to lowercase
      for (char &c : sw){ c = std::tolower(c); }

      //------------------------------------------------------------------------
      // Check for appropriate data outputs
      //------------------------------------------------------------------------
      if      (sw == "--xyz"   ){ vdc::xyz    = true; } // xyz coordinate file output
      else if (sw == "--povray"){ vdc::povray = true; } // pov coordinate file output
      else if (sw == "--track") { vdc::track  = true; }
      else if (sw == "--vtk"   ){ vdc::vtk    = true; } // vtk coordinate file output
      else if (sw == "--text" || sw == "--spins"){ vdc::txt = true; } // plain text / spin data file output
      else if (sw == "--binary"){ vdc::binary_dump = true; } // native-endian numeric dump of data files
      else if (sw == "--omp-threads"){
         check_arg(arg, argc, argv, temp_str, "Error - number of OpenMP threads required for \'--omp-threads\'." );
         try {
            vdc::omp_threads = std::stoi(temp_str);
         }
         catch(...){
            std::cerr << "Error - invalid value \'" << temp_str << "\' for --omp-threads." << std::endl;
            std::exit(EXIT_FAILURE);
         }
         if(vdc::omp_threads < 1){
            std::cerr << "Error - --omp-threads must be at least 1." << std::endl;
            std::exit(EXIT_FAILURE);
         }
      }
      else if (sw == "--cells" ){ // cell raw data
         vdc::cells  = true; // calculate cell data
         vdc::cellsf = true; // output cell data file
      }
      else if (sw == "--spin-cells" ){ // cell raw data
         vdc::cells  = true; // calculate cell data
         vdc::cellsf = true; // output cell data file
         vdc::spin_cellsf = true;
         vdc::povspincells = true; // enable povray cells output
      }
      else if (sw == "--povray-cells" ){ // povray cell data
         vdc::cells  = true; // calculate cell data
         vdc::povcells = true; // enable povray cells output
      }
      else if (sw == "--stray-field" ){
         vdc::require_stray_field_support();
         vdc::cells = true;         // reuse cell binning of spins / magnetisation
         vdc::stray_field = true;
      }
      else if (sw == "--sf-height"){
         check_arg(arg, argc, argv, temp_str, "Error - height in nm required for \'--sf-height\'." );
         try {
            vdc::sf_height_nm = std::stod(temp_str);
         }
         catch(...){
            std::cerr << "Error - invalid value \'" << temp_str << "\' for --sf-height." << std::endl;
            std::exit(EXIT_FAILURE);
         }
         if(vdc::sf_height_nm < 0.0){
            std::cerr << "Error - --sf-height must be >= 0 nm." << std::endl;
            std::exit(EXIT_FAILURE);
         }
      }
      else if (sw == "--sf-pad"){
         check_arg(arg, argc, argv, temp_str, "Error - integer padding factor required for \'--sf-pad\'." );
         try {
            vdc::sf_pad = std::stoi(temp_str);
         }
         catch(...){
            std::cerr << "Error - invalid value \'" << temp_str << "\' for --sf-pad." << std::endl;
            std::exit(EXIT_FAILURE);
         }
         if(vdc::sf_pad < 1){
            std::cerr << "Error - --sf-pad must be >= 1 (1 = periodic, 2 = finite flake)." << std::endl;
            std::exit(EXIT_FAILURE);
         }
      }
      else if (sw == "--ssc" || sw == "--spin-spin-correlation"){ vdc::ssc = true; }
      //------------------------------------------------------------------------
      // Check for cell size
      //------------------------------------------------------------------------
      /*else if (sw == "--cell-size"){
         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - size of cells in Angstroms." );

         // set cell size
         if ( stof(temp_str) >= 0.5 ){ vdc::cell_size = stof(temp_str); }
         else{
            std::cerr << "Error - cell size must be greater than 0.5 Angstroms." << std::endl;
            std::exit(EXIT_FAILURE);
         }
      }*/
      //------------------------------------------------------------------------
      // Check for verbose output
      //------------------------------------------------------------------------
      else if (sw == "--verbose" || sw == "-v"){ vdc::verbose = true; }
      //------------------------------------------------------------------------
      // Input file name
      //------------------------------------------------------------------------
      else if (sw == "--input-file"){
         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - no file specified for \'--input-file\' command line option." );

         // set input_file name
         vdc::input_file = temp_str;
      }
      //------------------------------------------------------------------------
      // Help information
      //------------------------------------------------------------------------
      else if (sw == "--help" || sw == "-h"){
         // check number of args not exceeded
         check_arg(arg, argc, argv, temp_str, "Error - no parameter specified for "+sw+" command line option." );

         // change temp_str to lowercase for comparison
         for (char &c : temp_str){
            c = std::tolower(c);
         }

         // check if parameter exists
         if (vdc::key_list.find(temp_str) == vdc::key_list.end()){
            std::cerr << "Error - No matching parameter '" << temp_str << "'.\n";
            std::exit(EXIT_FAILURE);
         }

         input_t input;

         // create argument to pass to function wrapper
         input.value = {"-h"};

         // pass to function wrapper
         vdc::key_list.at(temp_str)(input);
      }
      //------------------------------------------------------------------------
      // check for input file parameters passed in command line (parameters must use -- notation)
      //------------------------------------------------------------------------
      else if ( vdc::key_list.count(sw.substr(2)) ){

         input_t input;

         // set key. line number is set to -1 for error checking
         input.key = sw.substr(2);
         input.line_number = -1;

         // extract and sanitize parameter arguments, if any
         if (++arg < argc){
            temp_str = argv[arg];

            std::string delimiters = ",(){}[]:=!";
            std::string temp_val;

            // stop adding arguments when a command line paramter is reached (i.e. something starting with -- or -h)
            while ( !(temp_str[0] == '-' && (temp_str[1] == '-' || std::isalpha(temp_str[1]))) ){

               // work through argument
               for (char &c : temp_str){

                  // if uppercase, make lowercase
                  c = std::tolower(c);

                  // skip delimiters and push last value
                  if (delimiters.find(c) != std::string::npos){

                     // in case delimiters follow each other
                     if (temp_val == ""){ continue; }

                     // add to values
                     input.value.push_back(temp_val);

                     temp_val.clear();
                  }
                  // otherwise add char to temp_val
                  else { temp_val.push_back(c); }

               }

               // push back last value if not empty and clear temp_val
               if (temp_val != ""){
                  input.value.push_back(temp_val);
                  temp_val.clear();
               }

               // if more arguments remain, move to next, else finished
               if (++arg < argc){ temp_str = argv[arg]; }
               else { break; }
            }

            // reduce arg to process next parameter again
            arg--;
         }

         // keep track of parameter to check it is not used agin in input file
         // slices are excluded as multiple can be defined
         if (input.key.find("slice") == std::string::npos){
            vdc::cmdl_parameters.push_back(input.key);
         }


         // set parameters using function wrapper and arguments
         vdc::key_list.at(input.key)(input);
      }
      else {
         std::cerr << "Error - unknown command line parameter \'" << sw << "\'" << std::endl;
         std::exit(EXIT_FAILURE);
      }
   }

   //---------------------------------------------------------------------------
   // Additional checks on command line parameters
   //---------------------------------------------------------------------------

   // check that some kind of data output is requested
   if( !vdc::xyz && !vdc::povray && !vdc::vtk && !vdc::txt && !vdc::ssc && !vdc::cells && !vdc::track){
      std::cerr << "Error! No output data formats requested. Available options are: " << std::endl;
      std::cerr << "\t\t --xyz    Data output in .xyz format for viewing in rasmol/jmol" << std::endl;
      std::cerr << "\t\t --povray Data output in PoVRAY format for rendering" << std::endl;
      std::cerr << "\t\t --track  Data output in PoVRAY vector format for tracking rendering" << std::endl;
      std::cerr << "\t\t --vtk    Data output in VTK format for viewing in Paraview" << std::endl;
      std::cerr << "\t\t --text   Data output in plain text format for plotting in gnuplot/excel etc" << std::endl;
      std::cerr << "\t\t --spins  Same as --text" << std::endl;
      std::cerr << "\t\t --cells  Data output in plain text format in cells" << std::endl;
      std::cerr << "\t\t --spin-cells  PoVRAY cell-scale spin arrows (use same --slice, --camera-*, --colourmap, etc. as --povray)" << std::endl;
      std::cerr << "\t\t --stray-field  FFT stray field above the sample (requires FFTW3); --sf-height <nm> (default 30), --sf-pad <n> (default 2)" << std::endl;
      std::cerr << "\t\t --ssc    Spin-spin correlation data in text format" << std::endl;
      std::cerr << "\t\t --binary Native-endian binary dump of data files; also writes vdc-binary.meta" << std::endl;
      std::cerr << "\t\t --omp-threads N  OpenMP threads for --binary write path (ignored without --binary)" << std::endl;
      std::exit(EXIT_FAILURE);
   }
}

//------------------------------------------------------------------------------
// Check number of command line args not exceeded
//------------------------------------------------------------------------------
void check_arg( int &arg, int argc, char* argv[], std::string &temp_str, std::string error_output ){

   if (++arg < argc){ temp_str = argv[arg]; }
   else {
      std::cerr << error_output << std::endl;
      std::exit(EXIT_FAILURE);
   }
}

} // end of namespace vdc
