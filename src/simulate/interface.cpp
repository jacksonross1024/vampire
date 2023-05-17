//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <sstream>

// Vampire headers
#include "errors.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "stats.hpp"

// Internal sim header
#include "internal.hpp"

namespace sim{

   //-----------------------------------------------------------------------------
   // Function to process input file parameters for sim module
   //-----------------------------------------------------------------------------
   bool match_input_parameter(string const key, string const word, string const value, string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="sim";
      if(key!=prefix) return false;

      // set maximum allowable value for time steps (10^12)
      const uint64_t max_time = 1000000000000;
      const std::string max_time_str = "1,000,000,000,000";

      //----------------------------------
      // Now test for all valid options
      //----------------------------------
      std::string test="spin-transfer-torque-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::stt_polarization_unit_vector = u;
         return true;
      }
      //-------------------------------------------------------------------
      test="spin-orbit-torque-polarization-unit-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         // Test for valid range
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         // save sanitized unit vector
         sim::internal::sot_polarization_unit_vector = u;
         return true;
      }
      //-------------------------------------------------------------------
      test="preconditioning-steps";
      if(word==test){
         int n = atoi(value.c_str());
         // Test for valid range
         vin::check_for_valid_int(n, word, line, prefix, 0, 1000000,"input","0 - 1,000,000");
         sim::num_monte_carlo_preconditioning_steps = n;
         return true;
      }
      //-------------------------------------------------------------------
      test="time-step";
      if(word==test){
         double dt = atof(value.c_str());
         vin::check_for_valid_value(dt, word, line, prefix, unit, "time", 1.0e-20, 1.0e-6,"input","0.01 attosecond - 1 picosecond");
         mp::dt_SI = dt;
         return true;
      }
      //--------------------------------------------------------------------
      test="total-time-steps";
      if(word==test){
         uint64_t tt = vin::str_to_uint64(value); // convert string to uint64_t
         vin::check_for_valid_int(tt, word, line, prefix, 0, max_time,"input","0 - "+max_time_str);
         sim::total_time = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="loop-time-steps";
      if(word==test){
         uint64_t tt = vin::str_to_uint64(value); // convert string to uint64_t
         vin::check_for_valid_int(tt, word, line, prefix, 0, max_time,"input","0 - "+max_time_str);
         sim::loop_time = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="equilibration-time-steps";
      if(word==test){
         uint64_t tt = vin::str_to_uint64(value); // convert string to uint64_t
         vin::check_for_valid_int(tt, word, line, prefix, 0, max_time,"input","0 - "+max_time_str);
         sim::equilibration_time = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="time-steps-increment";
      if(word==test){
         uint64_t tt = vin::str_to_uint64(value); // convert string to uint64_t
         vin::check_for_valid_int(tt, word, line, prefix, 1, max_time,"input","1 - "+max_time_str);
         sim::partial_time = tt;
         return true;
      }
      test = "integrator";
      if( word == test ){
         //--------------------------------------------------------------------
         test="llg-heun";
         if( value == test ){
            sim::integrator = sim::llg_heun;
            return true;
         }
         //--------------------------------------------------------------------
         test="monte-carlo";
         if( value == test ){
            sim::integrator = sim::monte_carlo;
            return true;
         }
         //--------------------------------------------------------------------
         test="llg-midpoint";
         if( value == test ){
            sim::integrator = sim::llg_midpoint;
            return true;
         }
         //--------------------------------------------------------------------
         test="constrained-monte-carlo";
         if( value == test ){
            sim::integrator = sim::cmc;
            return true;
         }
         //--------------------------------------------------------------------
         test="hybrid-constrained-monte-carlo";
         if( value == test ){
            sim::integrator = sim::hybrid_cmc;
            return true;
         }
         //--------------------------------------------------------------------
         test="llg-quantum";
         if( value == test ){
            sim::integrator = sim::llg_quantum;
            return true;
         }
         test = "velocity-verlet";
         if (value == test) {
            sim::integrator = sim::velocity_verlet;
            return true;
         }
         //--------------------------------------------------------------------
         else{
            terminaltextcolor(RED);
               std::cerr << "Error - value for \'sim:" << word << "\' must be one of:" << std::endl;
               std::cerr << "\t\"llg-heun\"" << std::endl;
               std::cerr << "\t\"llg-midpoint\"" << std::endl;
               std::cerr << "\t\"llg-quantum\"" << std::endl;
               std::cerr << "\t\"monte-carlo\"" << std::endl;
               std::cerr << "\t\"constrained-monte-carlo\"" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
          }
      }
      //--------------------------------------------------------------------
      test="domain-wall-axis";
      if(word==test){
         //vin::check_for_valid_int(tt, word, line, prefix, 0, max_time,"input","0 - "+max_time_str);
         if (value == "x") {
         sim::domain_wall_axis = 0;
         }
         else if (value == "y") {
         sim::domain_wall_axis = 1;
         }
         else if (value == "z") {
         sim::domain_wall_axis = 2;
         }
         else {
            std::cout << "domain wall axis must equal x or y or z" <<std::endl;
            return false;
         }
         return true;
      }
      //--------------------------------------------------------------------
      test="domain-wall-discretisation";
      if(word==test){
         double tt = atof(value.c_str()); // convert string to uint64_t
         vin::check_for_valid_value(tt, word, line, prefix, unit, "length", 10, 1000,"input","10 - 1 A");
         sim::domain_wall_discretisation = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="domain-wall-anti-pbc-x";
      if(word==test){
         sim::anti_PBC[0] = true;
         cs::pbc[0]=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="domain-wall-anti-pbc-y";
      if(word==test){
         sim::anti_PBC[1] = true;
         cs::pbc[1]=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="domain-wall-anti-pbc-z";
      if(word==test){
         sim::anti_PBC[2] = true;
         cs::pbc[2]=true;
         return true;
      }
      //--------------------------------------------------------------------
      test="domain-wall-position";
      if(word==test){
         double tt = atof(value.c_str()); // convert string to uint64_t
         vin::check_for_valid_value(tt, word, line, prefix, unit, "none", 0, 1,"input","0 - 1");
         sim::domain_wall_position = tt;
         return true;
      }
      //--------------------------------------------------------------------
      test="domain-wall-width";
      if(word==test){
         double tt = atof(value.c_str()); // convert string to uint64_t
         vin::check_for_valid_value(tt, word, line, prefix, unit, "length", 0, 1000,"input","0 - 1 A");
         sim::domain_wall_width = tt;
         return true;
      }
      test = "convergence-criteria";
      if (word == test) {
         sim::calculate_program_convergence = true;
         stats::calculate_system_energy = true;
         stats::calculate_system_magnetization = true;
         double c = atof(value.c_str());
         vin::check_for_valid_value(c, word, line, prefix, unit, "unitless", 0.0, 1.0, "input", "0 - 1");
         sim::convergence_criteria = c;
         std::cout << "convergence-criteria: " << c << std::endl;
         return true;
      }
      test = "convergence-check";
      if (word == test) {
         sim::calculate_program_convergence = true;
         unsigned int c = atoi(value.c_str());
         vin::check_for_valid_int(c, word, line, prefix, 2, sim::equilibration_time, "input", "2 - equilibration-time-steps");
         sim::convergence_check = c;
     //    std::cout << "convergence-check: " << c << std::endl;
         return true;
      }
      test = "conduction-electrons";
      if (word == test) {
         double c = atof(value.c_str());
         vin::check_for_valid_value(c, word, line, prefix, unit, "e/atom", 0.0, 128, "input", "Number of conduction electrons / atom");
         sim::conduction_electrons = c;
         return true;
      }
      test = "applied-voltage";
      if(word == test) {
         double c = atof(value.c_str());
         vin::check_for_valid_value(c, word, line, prefix, unit, "V", 0.0, 1e20, "input", "0-1e10");
         sim::applied_voltage_sim = true;
         sim::applied_voltage = c;
         return true;
      }

      test = "fluence";
      if(word == test) {
         double c = atof(value.c_str());
         vin::check_for_valid_value(c, word, line, prefix, unit, "J", 0.0, 1e10, "input", "0-1e10");
         sim::heat_pulse_sim = true;
         sim::fluence = c;
         return true;
      }
      test = "photon-energy";
      if(word == test) {
         double c = atof(value.c_str());
         vin::check_for_valid_value(c, word, line, prefix, unit, "eV", 0.0, 1e10, "input", "0-1e10");
         sim::photon_energy = c;
         return true;
      }

      test = "x_vector";
      if(word == test) {
         double c = atof(value.c_str());
         vin::check_for_valid_value(c, word, line, prefix, unit, "unit-value", 0.0, 1.0, "input", "0-1; -1");
         sim::CASTLE_x_vector = c;
         return true;
      }
      test = "y_vector";
      if(word == test) {
         double c = atof(value.c_str());
         vin::check_for_valid_value(c, word, line, prefix, unit, "unit-value", 0.0, 1.0, "input", "0-1; -1");
         sim::CASTLE_y_vector = c;
         return true;
      }
      test = "z_vector";
      if(word == test) {
         double c = atof(value.c_str());
         vin::check_for_valid_value(c, word, line, prefix, unit, "unit-value", 0.0, 1.0, "input", "0-1; -1");
         sim::CASTLE_z_vector = c;
         return true;
      }
      test = "e-e_coupling";
      if(word == test) {
         sim::ee_coupling = true;
         double c = atof(value.c_str());
         vin::check_for_valid_value(c, word, line, prefix, unit, "AJ", -1.0, 1e10, "input", "0 - 1e10; -1");
         sim::ee_coupling_strength = c;
         return true;
      }
      test = "e-a_coupling";
      if(word == test) {
         sim::ea_coupling = true;
         double c = atof(value.c_str());
         vin::check_for_valid_value(c, word, line, prefix, unit, "AJ", -1.0, 1e20, "input", "0 - 1e10; -1");
         sim::ea_coupling_strength = c;
         return true;
      }
      test = "CASTLE-MD-rate";
      if(word == test) {
         int c = atoi(value.c_str());
         vin::check_for_valid_int(c, word, line, prefix, 1, 100000, "input", "1-1e10");
         sim::CASTLE_MD_rate = c;
         return true;
      }
      test = "ee-scattering-angle";
      if(word == test) {
         double a = atof(value.c_str());
         vin::check_for_valid_value(a, word, line, prefix, unit, "2pi radian", 0.0, 1.0, "input", "0 to 1.0, default = 0.1");
         sim::ee_scattering_angle = a;
         return true;
      }
      test = "omp-threads";
      if(word == test) {
         int t = atoi(value.c_str());
         vin::check_for_valid_int(t, word, line, prefix, 1, 64, "input", "1, 8 or 64");
         sim::CASTLE_omp_threads = t;
         return true;
      }
      //--------------------------------------------------------------------
      // input parameter not found here
      return false;
   }

   //----------------------------------------------------------------------------------
   // material parameter match function
   //----------------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index){

      // add prefix string
      std::string prefix="material:";

      // Check for material id > current array size and if so dynamically expand mp array
      if((unsigned int) super_index + 1 > sim::internal::mp.size() && super_index + 1 < 101) sim::internal::mp.resize(super_index + 1);

      //------------------------------------------------------------
      std::string test="domain-wall-second-magnetisation-vector";
      if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         std::cout << sim::domain_wall_second_vector_x.size() << "\t" << super_index << "\t" << u[0] << '\t' << u[1] << '\t' << u[2] <<std::endl;
         sim::domain_wall_second_vector_x[super_index] = u[0];
         sim::domain_wall_second_vector_y[super_index] = u[1];
         sim::domain_wall_second_vector_z[super_index] = u[2];
         return true;
      }
      //------------------------------------------------------------
      test  = "spin-transfer-relaxation-torque";
      std::string test2 = "slonczewski-adiabatic-spin-torque";
      std::string test3 = "spin-transfer-torque";
      std::string test4 = "antidamping-torque";
      // aj parameter for material in slonczewski torque calculation
      if( word==test || word==test2 || word==test3 || word==test4){
         double aj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(aj, word, line, prefix, unit, "field", -1.0e-2, 1.0e2,"input","-100 - 100T");
         sim::internal::mp[super_index].stt_rj.set(aj);
         sim::internal::enable_spin_torque_fields = true;
         return true;
      }
      //------------------------------------------------------------
      test2 = "spin-transfer-precession-torque";
      test  = "slonczewski-non-adiabatic-spin-torque";
      test3 = "field-like-torque";
      test4 = "slonczewski-precession-spin-torque";
      // bj parameter for material in slonczewski torque calculation
      if( word==test || word==test2 || word==test3 ){
         double bj=atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(bj, word, line, prefix, unit, "field", -1.0e-2, 1.0e2,"input","-100 - 100T");
         sim::internal::mp[super_index].stt_pj.set(bj);
         sim::internal::enable_spin_torque_fields = true;
         return true;
      }
      //------------------------------------------------------------
      test = "spin-transfer-torque-asymmetry";
      // damping-like parameter for material in spin orbit torque calculation
      if( word==test ){
         double sttasm = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(sttasm, word, line, prefix, unit, "", 0.0, 1.0e2,"input","0 - 100");
         sim::internal::mp[super_index].stt_asm.set(sttasm);
         return true;
      }
      //------------------------------------------------------------
      // field-like parameter for material in spin orbit torque calculation
      test = "spin-orbit-relaxation-torque";
      test2 = "spin-orbit-anti-damping-torque";
      if( word==test || word==test2 ){
         double aj = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(aj, word, line, prefix, unit, "field", -1.0e2, 1.0e2,"input","-100 - 100T");
         sim::internal::mp[super_index].sot_rj.set(aj);
         sim::internal::enable_spin_torque_fields = true;
         return true;
      }
      //------------------------------------------------------------
      test = "spin-orbit-precession-torque";
      test2 = "spin-orbit-torque";
      test3 = "spin-orbit-field-like-torque";
      // damping-like parameter for material in spin orbit torque calculation
      if( word==test || word==test2 || word==test3 ){
         double bj = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(bj, word, line, prefix, unit, "field", -1.0e2, 1.0e2,"input","-100 - 100T");
         sim::internal::mp[super_index].sot_pj.set(bj);
         sim::internal::enable_spin_torque_fields = true;
         return true;
      }
      //------------------------------------------------------------
      test = "spin-orbit-torque-asymmetry";
      // damping-like parameter for material in spin orbit torque calculation
      if( word==test ){
         double sotasm = atof(value.c_str());
         // Test for valid range
         vin::check_for_valid_value(sotasm, word, line, prefix, unit, "", 0.0, 1.0e2,"input","0 - 100");
         sim::internal::mp[super_index].sot_asm.set(sotasm);
         return true;
      }
      test = "perfect-fermi-gas";
      if (word == test) {
         sim::calculate_fermi_distribution = true;
         return true;
      }
      test = "conduction-electrons";
      

      //--------------------------------------------------------------------
      // keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of namespace sim
