//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <fstream>
#include <iostream>

// Vampire headers
#include "spintorque.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "sim.hpp"


// Spin Torque headers
#include "internal.hpp"

namespace st{

bool match_material(string const word, string const value, string const unit, int const line, int const super_index){

   // add prefix string
   std::string prefix="material:";

   // Check for material id > current array size and if so dynamically expand mp array
   if((unsigned int) super_index + 1 > st::internal::mp.size() && super_index + 1 < 101) st::internal::mp.resize(super_index + 1);

   //------------------------------------------------------------
   std::string test="spin-diffusion-length"; 
   /*
      float spin-diffusion-length 
         Details
      */
   if(word==test){
      double lsdl=atof(value.c_str());
      //std::cout << "lsdl" << super_index << "\t" << lsdl << "\t" << value.c_str() << std::endl;
      vin::check_for_valid_value(lsdl, word, line, prefix, unit, "length", 0.01, 1.0e10,"material"," 0.01 - 1e10 Angstroms");
      st::internal::mp[super_index].lambda_sdl=lsdl*1.e-10; // defined in metres
      std::cout << "lsdl" << super_index << "\t" << lsdl << std::endl;
      return true;
   }
   //--------------------------------------------------------------------
   test="spin-polarisation-conductivity"; //beta (eq 2)
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double betac=atof(value.c_str());
      vin::check_for_valid_value(betac, word, line, prefix, unit, "none", 1.0e-3, 1.0e3,"material"," 0.001 - 1000");
      st::internal::mp[super_index].beta_cond=betac;
      return true;
   }
   //--------------------------------------------------------------------
   test="spin-polarisation-diffusion"; //beta'
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double betad=atof(value.c_str());
      vin::check_for_valid_value(betad, word, line, prefix, unit, "none", 1.0e-3, 1.0e3,"material"," 0.001 - 1000");
      st::internal::mp[super_index].beta_diff=betad;
      return true;
   }
   //--------------------------------------------------------------------
   test="spin-accumulation"; //m(infinity)
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double sa=atof(value.c_str());
      vin::check_for_valid_value(sa, word, line, prefix, unit, "none", 0.0, 1.0e10,"material"," 1.0e-6 - 1.0e10");
      st::internal::mp[super_index].sa_infinity=sa;
      return true;
   }
   //--------------------------------------------------------------------
   test="diffusion-constant"; //D_0 (eq 2)
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double dc=atof(value.c_str());
      vin::check_for_valid_value(dc, word, line, prefix, unit, "none", 1.0e-9, 100.0,"material"," 1.0e-9 - 100"); //m^2/s
      st::internal::mp[super_index].diffusion=dc;
      return true;
   }
   //--------------------------------------------------------------------
   test="sd-exchange-constant";// 
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double sd=atof(value.c_str());
      vin::check_for_valid_value(sd, word, line, prefix, unit, "energy", 0.0, 1.e-17,"material"," 1.0e-30 - 1.e-17 J");
      st::internal::mp[super_index].sd_exchange=sd;
      return true;
   }
   //--------------------------------------------------------------------
   test="spin-torque-free-layer"; //
   /*
    * Spin torque free layer flag
    */
   if(word==test){
      st::internal::free_layer = super_index;
      return true;
   }
   //--------------------------------------------------------------------
   test="spin-torque-reference-layer";
   /*
    * Spin torque free layer flag
    */
   if(word==test){
      st::internal::reference_layer = super_index;
      return true;
   }

   //--------------------------------------------------------------------
   // keyword not found
   //--------------------------------------------------------------------
   return false;

}



   //-----------------------------------------------------------------------------
   // Function to process input file parameters for ST module
   //-----------------------------------------------------------------------------
   bool match_input_parameter(string const key, string const word, string const value, string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="spin-torque";


      if(key!=prefix) return false;

      //----------------------------------
      // Now test for all valid options
      //----------------------------------

      std::string test="enable-ST-fields";

       if(word==test){
         st::internal::enabled = true;
         return true;
       }

       //-------------------------------------------------

       test="TMRenable"; //set false

       if(word==test){
         st::internal::TMRenable = true;
         return true;
       }


      //-------------------------------------------------
      test="current-density"; //
       if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e15,"input","0.0 - 1.0e13 A/m2");
         st::internal::je =T;
         return true;
        }

      //-------------------------------------------------

       test="current-direction";
        if(word==test){
         std::string je_dir="x";
             if(value==je_dir){
                st::internal::current_direction=0;
                return true;
             }
             else
               je_dir ="y";
             if(value==je_dir){
              st::internal::current_direction=1;
              return true;
            }
            else{
              st::internal::current_direction=2;
              return true;
            }
        }

      //-------------------------------------------------

       test="micro-cell-size";
       if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         //vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         st::internal::micro_cell_size[0]=u.at(0);
         st::internal::micro_cell_size[1]=u.at(1);
         st::internal::micro_cell_size[2]=u.at(2);
         std::cout << "microcell size " << st::internal::micro_cell_size[0] << ", " << st::internal::micro_cell_size[1] << ", " << st::internal::micro_cell_size[2] << std::endl;
         return true;
      }
      test = "micro-cell-decomp";
      if(word == test){
         std::string decomp_value = "mpi";
         if(decomp_value == value) {
            st::internal::microcell_decomp_type = decomp_value;
            return true;
         }
         decomp_value = "A";
         if(decomp_value == value) {
            st::internal::microcell_decomp_type = decomp_value;
            return true;
         }
         terminaltextcolor(RED);
            std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"mpi\"" << std::endl;
            std::cerr << "\t\"A\"" << std::endl;
         terminaltextcolor(WHITE);
            return false;

      }
      //-------------------------------------------------

      test="micro-cell-thickness";
       if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 20,"input","0.0 - 20 A");
         st::internal::micro_cell_thickness =T;
         return true;
        }

     //-------------------------------------------------

       test="initial-spin-polarisation"; //
       if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 20,"input","0.0 - 20 C/m^3");
         st::internal::initial_beta =T;
         return true;
        }

      test="spin-hall-angle"; //
       if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 10,"input","Spin Hall Angle: 0.0 - 1.0");
         st::internal::initial_theta =T;
         return true;
        }

      test="remove-material-type"; //
       if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 20,"input","material type");
         st::internal::remove_nm =T;
         return true;
        }
    //--------------------------------------------------------------------
      test="initial-mag-direction";//might not matter
      if(word==test){
      std::vector<double> u(3);
      u=vin::doubles_from_string(value);
      vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
      st::internal::initial_m[0]=u.at(0);
      st::internal::initial_m[1]=u.at(1);
      st::internal::initial_m[2]=u.at(2);
      return true;
   }
   test="finite-boundary-condition";
   if(word==test){
      if(value == "true"){
         st::internal::fbc = true;
      }  else if(value == "false"){
         st::internal::fbc = false;
      } else {
         terminaltextcolor(RED);
                std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
                std::cerr << "\t\"true\"" << std::endl;
                std::cerr << "\t\"false\"" << std::endl;
            terminaltextcolor(WHITE);
               return false;
      }

      return true;
   }

   test = "output-torque-data";
   if(word == test) {
      std::string output_type = "init";
      if(value == output_type){
         st::internal::output_torque_data = output_type;
         return true;
      }
      output_type = "int";
      if(value == output_type) {
         st::internal::output_torque_data = output_type;
         return true;
      }
      output_type = "final";
      if(value == output_type) {
         st::internal::output_torque_data = output_type;
         return true;
      }
      
      terminaltextcolor(RED);
      std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
      std::cerr << "\t\"init\"" << std::endl;
      std::cerr << "\t\"int\"" << std::endl;
      std::cerr << "\t\"final\"" << std::endl;
      terminaltextcolor(WHITE);
      return false;
      
   }

   test = "sot-check-only";
   if(word == test) {
      st::internal::sot_check = true;

      return true;
   }
      //-------------------------------------------------

   test="ST-output-rate";
   if(word==test){
      int T=atoi(value.c_str());
      vin::check_for_valid_int(T, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
      st::internal::ST_output_rate =T;
      return true;
   }
   test="SOT-spin-accumulation";
   if(word==test){
      
      st::internal::sot_sa =true;
      return true;
   }

      //--------------------------------------------------------------------
      // input parameter not found here
      return false;
   }











} // end of namespace st
