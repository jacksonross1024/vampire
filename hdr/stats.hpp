//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2011-2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
#ifndef STATS_H_
#define STATS_H_

// C++ include files
#include <vector>
#include <string>

#include "atoms.hpp"

namespace stats
//==========================================================
// Namespace statistics
//==========================================================
{
	extern int num_atoms;				/// Number of atoms for statistic purposes
	extern double inv_num_atoms;	///1.0/num_atoms
	extern double max_moment;		/// Total Maximum moment
	extern double data_counter;		/// number of data points for averaging

	// Member Functions
	extern int mag_m();
	extern void mag_m_reset();
	extern double max_torque();

	extern bool calculate_torque;
	extern double total_system_torque[3];
	extern double total_mean_system_torque[3];

	extern std::vector <double> sublattice_mean_torque_x_array;
	extern std::vector <double> sublattice_mean_torque_y_array;
	extern std::vector <double> sublattice_mean_torque_z_array;

	extern double energy_data_counter;
	extern double torque_data_counter;

   extern bool calculate_energy;

   /// Statistics energy types
   enum energy_t { total = 0, exchange = 1, anisotropy = 2, applied_field = 3, magnetostatic = 4};

   /// Statistics types
   enum stat_t { atotal=0, mean=1};

   /// Statistics output functions
   extern void output_energy(std::ostream&, enum energy_t, enum stat_t,bool header);

   //-------------------------------------------------
   // New statistics module functions and variables
   //-------------------------------------------------

   // Control functions
   void initialize(const int num_atoms, const int num_materials,
                   const std::vector<double>& magnetic_moment_array,
                   const std::vector<int>& material_type_array,
                   const std::vector<int>& height_category_array,
                   const std::vector<bool>& non_magnetic_materials_array);

   void update(const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz,
               const std::vector<double>& mm, const std::vector<int>& mat, const double temperature);
   void reset();

   // Statistics control flags (to be moved internally when long-awaited refactoring of vio is done)
   extern bool calculate_system_energy;
   extern bool calculate_material_energy;
   extern bool calculate_system_magnetization;
   extern bool calculate_material_magnetization;
   extern bool calculate_height_magnetization;
   extern bool calculate_material_height_magnetization;
   extern bool calculate_system_specific_heat;
   extern bool calculate_material_specific_heat;
   extern bool calculate_system_standard_deviation;// AJN
   extern bool calculate_system_susceptibility;
   extern bool calculate_material_susceptibility;
   extern bool calculate_system_spin_temperature;
   extern bool calculate_material_spin_temperature;

   


   // forward declaration of friend classes
   class susceptibility_statistic_t;
   class specific_heat_statistic_t;

   class standard_deviation_statistic_t;
   class spin_temperature_statistic_t;
   class convergence_statistic_t;





//--------------------------------
// convergence class definition
//--------------------------------
   class convergence_statistic_t {
      public:
         convergence_statistic_t (std::string n) {
            name = n;
 
    };
    //----
         void initialize(double criteria, unsigned int convergence_check_value);
         void initialize_first_pass();
         void initialize_second_pass();
         bool did_converge();
         bool get_method();
         unsigned int get_converged_counter();
         void is_converged(unsigned int n_steps);
         void update_convergence();
         void update_data(const std::vector<double>& sx, // spin unit vector
               const std::vector<double>& sy,
               const std::vector<double>& sz,
               const std::vector<double>& mm,
               const std::vector<int>& mat,
               const double temperature);
         void update_counter();
         void update_mean_convergence();
         std::string output_convergence_rate(bool header);
         void do_converge();
         void end_convergence();
         void output_convergence();



      private:
         bool initialized;
         bool magnetisation_converged;      // did it converge
         bool energy_converged;
         bool converged;
         bool convergence_method; // do we check the convergence
         unsigned int counter;
         unsigned int converged_counter;
         bool output_convergence_counter;
         std::string name;
         std::vector<double> convergence_magnetisation;
         std::vector<double> convergence_energy;
         double convergence_criteria;  //convergence must reach below this error value
         unsigned int convergence_check; // check convergence every ith time
         std::vector<unsigned int> convergence_output_list;
         

   };

   //----------------------------------
   // Energy class definition
   //----------------------------------
   class energy_statistic_t{

      friend class specific_heat_statistic_t;

   public:
      energy_statistic_t (std::string n):initialized(false){
        name = n;
      };
      bool is_initialized();
      void set_mask(const int in_mask_size, const std::vector<int> in_mask);
      void get_mask(std::vector<int>& out_mask, std::vector<double>& out_normalisation);
      void calculate(const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz,
                     const std::vector<double>& mm, const std::vector<int>& mat, const double temperature);

      void reset_averages();

      void set_total_energy(         std::vector<double>& new_energy, std::vector<double>& new_mean_energy);
      void set_exchange_energy(      std::vector<double>& new_energy, std::vector<double>& new_mean_energy);
      void set_anisotropy_energy(    std::vector<double>& new_energy, std::vector<double>& new_mean_energy);
      void set_applied_field_energy( std::vector<double>& new_energy, std::vector<double>& new_mean_energy);
      void set_magnetostatic_energy( std::vector<double>& new_energy, std::vector<double>& new_mean_energy);

      const std::vector<double>& get_total_energy();
      const double get_total_system_energy();
      const double get_mean_total_system_energy();
      const std::vector<double>& get_exchange_energy();
      const std::vector<double>& get_anisotropy_energy();
      const std::vector<double>& get_applied_field_energy();
      const std::vector<double>& get_magnetostatic_energy();

      void update_mean_counter(long counter);

      std::string output_energy(enum energy_t energy_type, bool header);
      std::string output_mean_energy(enum energy_t energy_type, bool header);

   private:
      bool initialized;
      int num_atoms;
      int mask_size;
      double mean_counter;

      std::vector<int> mask;

      std::vector<double> total_energy;
      std::vector<double> exchange_energy;
      std::vector<double> anisotropy_energy;
      std::vector<double> applied_field_energy;
      std::vector<double> magnetostatic_energy;

      std::vector<double> mean_total_energy;
      std::vector<double> mean_exchange_energy;
      std::vector<double> mean_anisotropy_energy;
      std::vector<double> mean_applied_field_energy;
      std::vector<double> mean_magnetostatic_energy;

      std::vector<int> zero_list;
      std::vector<double> normalisation;

      std::string name;

   };

   //----------------------------------
   // Magnetization Class definition
   //----------------------------------
   class magnetization_statistic_t{

      friend class susceptibility_statistic_t;
      friend class standard_deviation_statistic_t;
      friend class spin_temperature_statistic_t;

      public:
         magnetization_statistic_t (std::string n):initialized(false){
           name = n;
         };
         bool is_initialized();
         void set_mask(const int mask_size, std::vector<int> inmask, const std::vector<double>& mm);
         void get_mask(std::vector<int>& out_mask, std::vector<double>& out_saturation);
         void calculate_magnetization(const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz, const std::vector<double>& mm);
         void set_magnetization(std::vector<double>& magnetization, std::vector<double>& mean_magnetization, long counter);
         void reset_magnetization_averages();
         const std::vector<int>& get_mask_only();
         const int get_mask_size();
         const double get_system_magnetization_length();
         const int get_mask_length();
         const double get_normalized_mean_system_magnetization_length();
         const std::vector<double>& get_magnetization();
         std::string output_magnetization(bool header);
         std::string output_normalized_magnetization(bool header);
         std::string output_normalized_magnetization_length(bool header);
         std::string output_normalized_mean_magnetization(bool header);
         std::string output_normalized_mean_magnetization_length(bool header);
         std::string output_normalized_magnetization_dot_product(const std::vector<double>& vec,bool header);
         std::string output_mean_magnetization_length(bool header);
         std::string output_mean_magnetization(bool header);

      private:
         bool initialized;
         int num_atoms;
         int mask_size;
         double mean_counter;
         std::vector<int> mask;
         std::vector<double> magnetization;
         std::vector<double> mean_magnetization;
         std::vector<int> zero_list;
         std::vector<double> saturation;
         std::string name;

   };

	//----------------------------------
   // Specific Heat Class definition
   //----------------------------------
   class specific_heat_statistic_t{

      public:
         specific_heat_statistic_t (std::string n):initialized(false){
           name = n;
         };
         void initialize(energy_statistic_t& energy_statistic);
         void calculate(const std::vector<double>& energy);
         void reset_averages();
         std::string output_mean_specific_heat(const double temperature,bool header);


      private:
         bool initialized;
         int num_elements;
         double mean_counter;
         std::vector<double> mean_specific_heat;
         std::vector<double> mean_specific_heat_squared;
         std::vector<double> normalisation;

         std:: string name;

   };
    //----------------------------------
   //Spin_Temperature_Class_Definition
   //----------------------------------

   class spin_temperature_statistic_t{

      public:
         spin_temperature_statistic_t(std::string n):initialized(false) {
            name = n;
         };
         void calculate(const std::vector<double>& sx, const std::vector<double>& sy, const std::vector<double>& sz, const std::vector<double>& sm,
         const std::vector<double>& Hx_int, const std::vector<double>& Hy_int,const std::vector<double>& Hz_int,
         const std::vector<double>& Hx_ext, const std::vector<double>& Hy_ext, const std::vector<double>& Hz_ext);
         void set_mask(const int mask_length, const std::vector<int>& mag_mask);
         std::string output_spin_temperature(bool header);
         std::string output_mean_spin_temperature(bool header);
         void reset();

      private:
         bool initialized;
         int mask_size; //how many different materials
         int num_atoms; //total number of atoms
       
         std::string name;
         std::vector<int> zero_list;
         std::vector<int> mask; 
         std::vector <double> total_spin_temperature;
         std::vector<double> spin_temperature_top;
         std::vector<double> spin_temperature_bottom;
         std::vector<double> mean_spin_temperature;
         std::vector<double> mean_spin_counter;

   };

   //----------------------------------
   // Susceptibility Class definition
   //----------------------------------
   class susceptibility_statistic_t{

      public:
         susceptibility_statistic_t (std::string n):initialized(false){
           name = n;
         };
         void initialize(magnetization_statistic_t& mag_stat);
         void calculate(const std::vector<double>& magnetization);
         void reset_averages();
         std::string output_mean_susceptibility(const double temperature,bool header);
         //std::string output_mean_absolute_susceptibility();

      private:
         bool initialized;
         int num_elements;
         double mean_counter;
         std::vector<double> mean_susceptibility;
         std::vector<double> mean_susceptibility_squared;
         std::vector<double> mean_absolute_susceptibility;
         std::vector<double> mean_absolute_susceptibility_squared;
         std::vector<double> saturation;
         std::string name;

   };
   //----------------------------------
   // Standard Deviation of magnetisation in time Class definition
   //----------------------------------
   class standard_deviation_statistic_t{ // AJN

      public:
         standard_deviation_statistic_t (std::string n):initialized(false){
           name=n;
         };
         void initialize(magnetization_statistic_t& mag_stat);
         bool is_initialized();
         void calculate(const std::vector<double>& magnetization);
         void reset_averages();
         std::string output_standard_deviation(bool header);
         double get_standard_deviation_length();
         std::string output_standard_deviation_length(bool header);

      private:
         bool initialized;
         int num_elements;//number of elements in the system
         int idx;// index for looping through directions
         double mean_counter;// counts time steps in loop - factor out for normal time
         double res1; // residuals calculated at each time
         double res2;
         std::vector<double> residual_sq;// running squared residual for each direction
         std::vector<double> mean; // running mean for each direction
         std::string name;
         std::vector<double> stdDv_array;

   };

   // Statistics classes
   extern energy_statistic_t system_energy;
   extern energy_statistic_t material_energy;

   extern magnetization_statistic_t system_magnetization;
   extern magnetization_statistic_t material_magnetization;
   extern magnetization_statistic_t height_magnetization;
   extern magnetization_statistic_t material_height_magnetization;

   extern specific_heat_statistic_t system_specific_heat;
   extern specific_heat_statistic_t material_specific_heat;

   extern susceptibility_statistic_t system_susceptibility;
   extern susceptibility_statistic_t material_susceptibility;
   extern standard_deviation_statistic_t system_standard_deviation;

    extern spin_temperature_statistic_t system_spin_temperature;
    extern spin_temperature_statistic_t material_spin_temperature;

    //convergence classes
    extern convergence_statistic_t program_convergence;
   
}

#endif /*STATS_H_*/
