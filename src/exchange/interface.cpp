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
#include <string>

// Vampire headers
#include "exchange.hpp"
#include "errors.hpp"
#include "vio.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{

   // internal namespacve for exchange module
   namespace internal{

      //-----------------------------------------------------------------------------------------------------------------------
      // general function for reading exchange value from input file and storing in 4D format
      //-----------------------------------------------------------------------------------------------------------------------
      void read_exchange_values(int material_i, int material_j, int neighbour, std::string const word, std::string const prefix, std::string const value, std::string const unit, int const line, exchange_matrix_4D_t& exchange_matrix){

         // extract comma separated values from string
         std::vector<double> Jij = vin::doubles_from_string(value);

         // optional factor 2 correction for ab-initio
         const double ef = internal::exchange_factor;
         for(auto& J : Jij) J *= ef;

         if(Jij.size() == 1){
            vin::check_for_valid_value(Jij[0], word, line, prefix, unit, "energy", -1e18, 1e18,"material"," < +/- 1.0e18");
            // set exchange constants
            exchange_matrix.set_exchange_values(material_i, material_j, neighbour, Jij);
         }
         else if(Jij.size() == 3){
            vin::check_for_valid_vector(Jij, word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e18");
            // set exchange constants
            exchange_matrix.set_exchange_values(material_i, material_j, neighbour, Jij);
            internal::minimum_needed_exchange_type = exchange::vectorial;
         }
         else if(Jij.size() == 9){
            vin::check_for_valid_vector(Jij, word, line, prefix, unit, "energy", -1e-18, 1e-18,"material"," < +/- 1.0e18");
            // set exchange constants
            exchange_matrix.set_exchange_values(material_i, material_j, neighbour, Jij);
            internal::minimum_needed_exchange_type = exchange::tensorial;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error in input file - material[" << material_i << "]:exchange_matrix[" << material_j << "] must have one, three, or nine values." << std::endl;
            terminaltextcolor(WHITE);
            zlog << zTs() << "Error in input file - material[" << material_i << "]:exchange_matrix[" << material_j << "] must have one, three, or nine values." << std::endl;
            err::vexit();
         }

         return;

      }

   } // end of internal namespace

   //---------------------------------------------------------------------------
   // Function to process input file parameters for exchange module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="exchange";
      if(key!=prefix) return false;

      //-------------------------------------------------------------------
      std::string test="dmi-cutoff-range";
      if(word==test){
          double cr = atof(value.c_str());
          // Test for valid range
          vin::check_for_valid_value(cr, word, line, prefix, unit, "length", 0.0, 1.0e3,"input","0.0 - 1e3");
          internal::dmi_cutoff_range = cr;
          return true;
      }
      //-------------------------------------------------------------------
      test="kitaev-cutoff-range";
      if(word==test){
          double cr = atof(value.c_str());
          // Test for valid range
          vin::check_for_valid_value(cr, word, line, prefix, unit, "length", 0.0, 1e3,"input","0.0 - 1e3");
          internal::kitaev_cutoff_range = cr;
          return true;
      }
      //--------------------------------------------------------------------
      if( word == "ab-initio-factor" ){
         double ab = atof(value.c_str());
         vin::check_for_valid_value(ab, word, line, prefix, unit, "scaling", -10.0, 10.0,"input","-10.0 - 10");
         internal::exchange_factor = ab;
         std::cout << "Note: exchange factor set to " << ab << " only for normalised-isotropic exchange " << std::endl;
         
         return true;
      }

       if( word == "meV-interaction-unit" ){
         internal::meV_interactions = true;
         return true;
      }
      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index, const int max_materials){

      // add prefix string
      std::string prefix="material:";

      // Check for empty material parameter array and resize to avoid segmentation fault
      if(internal::mp.size() == 0){
         internal::mp.resize(max_materials);
      }

      //------------------------------------------------------------
      // Check for material properties
      //------------------------------------------------------------
      std::string test = "dmi-constant"; // short form
      std::string test2 = "dzyaloshinskii-moriya-interaction-constant"; // long form
      if( (word == test) || (word == test2) ){
         double dmi = atof(value.c_str());
         vin::check_for_valid_value(dmi, word, line, prefix, unit, "energy", -1e-17, 1e-17,"material"," < +/- 1.0e17");
         internal::mp[super_index].dmi[sub_index] = dmi;
         internal::enable_dmi = true; // Switch on dmi calculation and fully unrolled tensorial anisotropy
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 0, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-1st-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 0, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-2nd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 1, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-3rd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 2, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-4th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 3, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-5th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 4, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-6th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 5, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-7th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 6, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-8th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 7, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-9th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 8, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-10th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 9, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      test = "exchange-matrix-11th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 10, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-12th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 11, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-13th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 12, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-14th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 13, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-15th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 14, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-16th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 15, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-17th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 16, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-18th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 17, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-19th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 18, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-20th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 19, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-21st-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 20, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-22nd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 21, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-23rd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 22, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-24th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 23, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      test = "exchange-matrix-25th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 24, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-26th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 25, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-27th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 26, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-28th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 27, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-29th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 28, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-30th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 29, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      test = "exchange-matrix-31st-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 30, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-32nd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 31, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-33rd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 32, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-34th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 33, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-35th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 34, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-36th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 35, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-37th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 36, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-38th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 37, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-39th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 38, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-40th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 39, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-41st-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 40, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-42nd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 41, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-43rd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 42, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-44th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 43, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-45th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 44, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-46th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 45, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-47th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 46, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-48th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 47, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-49th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 48, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-50th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 49, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-51st-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 50, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-52nd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 51, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-53rd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 52, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "exchange-matrix-44th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 43, word, prefix, value, unit, line, internal::bilinear_exchange_constants);
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 0, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-1st-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 0, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-2nd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 1, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-3rd-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 2, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-4th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 3, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-5th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 4, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-6th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 5, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-7th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 6, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-8th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 7, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-9th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 8, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //------------------------------------------------------------------------
      test = "biquadratic-exchange-matrix-10th-nn";
      if( word == test ){
         read_exchange_values(super_index, sub_index, 9, word, prefix, value, unit, line, internal::biquadratic_exchange_constants);
         exchange::biquadratic = true; // Switch on biquadratic exchange
         return true;
      }
      //--------------------------------------------------------------------
      test = "kitaev-constant"; // short form
      if( word == test ){
         double k = atof(value.c_str());
         vin::check_for_valid_value(k, word, line, prefix, unit, "energy", -1e-17, 1e-17,"material"," < +/- 1.0e17");
         internal::mp[super_index].kitaev[sub_index] = k;
         internal::enable_kitaev = true; // Switch on kitaev calculation and fully unrolled tensorial anisotropy
         return true;
      }
      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of exchange namespace
