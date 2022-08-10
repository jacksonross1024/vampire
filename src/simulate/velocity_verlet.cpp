// ==============================================
//   Coupled Atomistic Spin Thermalized Lattice Environment
//
//  =========    ========     ========  ============  ||           =========
// ||          ||        ||  ||              ||       ||          ||
// ||          ||        ||  ||              ||       ||          ||
// ||          ||        ||  ||              ||       ||          ||
// ||          || ====== ||   ========       ||       ||           =========
// ||          ||        ||          ||      ||       ||          ||
// ||          ||        ||          ||      ||       ||          ||
// ||          ||        ||          ||      ||       ||          ||
//  =========                 ========                  =========   =========
//

//
//   This file is part of the VAMPIRE open source package and its subset CASTLE under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans and J L Ross 2021. All rights reserved.
//
//   Email: richard.evans@york.ac.uk and jackson.ross@york.ac.uk
//
//
//==============================================
#include <string>
#include <vector>
#include <math.h>

#include "CASTLE.hpp"

namespace CASTLE {
 //ouble return_phonon_distribution(const double& epsilon, const double& beta );


int velocity_verlet_step(double time_step) {
  
    TEKE = 0.0;
    TLE  = 0.0;
     p_x = 0.0;
       p_y = 0.0;
       p_z = 0.0;
          //  std::cout << "Updating new electron position." << std::endl;
    update_position();

         // std::cout << current_time_step << "; positions updated " << std::endl;
    update_dynamics();
        //  std::cout << "Output mean data" << std::endl;
    
    if (current_time_step % CASTLE_output_rate == 0)   output_data(); //std::cout << "x_flux: " << x_flux / CASTLE_output_rate << "\n"; x_flux = 0;
      //    std::cout << current_time_step << "; output updated " << std::endl;
    //reset integration
    CASTLE_real_time += dt;
    current_time_step++;
   
    return EXIT_SUCCESS;
}

void setup_output() {

    CASTLE_output_data = true;
}


void update_position(){

    if(current_time_step % half_int_var == 0) {
      
      old_cell_integration_lists.swap(cell_integration_lists);
      #pragma simd
      for(int c = 0; c < total_cells; c++) {
        cell_integration_lists.at(c).at(0) = 1;
      }
    }
    omp_set_dynamic(0);
    omp_set_num_threads(omp_threads);
 
 
  #pragma omp parallel reduction(+:x_flux,y_flux,z_flux, p_x,p_y,p_z)
  {
  for (int l = 0 ; l < cells_per_thread; l++) {
    const unsigned int cell = lattice_cells_per_omp.at(omp_get_thread_num()).at(l);
    const unsigned int size = old_cell_integration_lists.at(cell).at(0);
   // cell_integration_lists.at(cell).at(0) = 1;

  //  if((current_time_step % half_int_var) == 0 && ( (l % ((z_omp_cells/2)-1) )== 0 || (l % (z_omp_cells/2)) == 0) ) {
      #pragma omp barrier
   // }
    
    for (int e = 1; e < size; e++) { 
      const unsigned int electron = old_cell_integration_lists.at(cell).at(e);
      const unsigned int array_index = 3*electron;
    //  const int array_index_y = array_index + 1;
     // const int array_index_z = array_index + 2;
  
      double x_pos = electron_position.at(array_index);
      double y_pos = electron_position.at(array_index+1);
      double z_pos = electron_position.at(array_index+2);
    
     // if (x_pos != x_pos || y_pos != y_pos || z_pos != z_pos) std::cout << "access pos " << x_pos << ", " << y_pos << ", " << z_pos << ", " << \
      electron_position.at(array_index) << ", " << electron_position.at(array_index+1) << ", " << electron_position.at(array_index+2) << ", " << electron << ", " << array_index << std::endl;
      
        p_x += electron_velocity.at(array_index);
        p_y += electron_velocity.at(array_index+1) ;
        p_z += electron_velocity.at(array_index+2) ;
        
      if(electron_transport_list.at(electron)) {

        x_pos += electron_velocity.at(array_index)   * dt;// + (electron_force.at(array_index)   * dt * dt * constants::K_A / 2); // x superarray component
        y_pos += electron_velocity.at(array_index+1) * dt;// + (electron_force.at(array_index_y) * dt * dt * constants::K_A / 2); // y superarray component
        z_pos += electron_velocity.at(array_index+2) * dt;// + (electron_force.at(array_index_z) * dt * dt * constants::K_A / 2); // z superarray component

      //  if (x_pos != x_pos || y_pos != y_pos || z_pos != z_pos) std::cout << "update pos " << x_pos << ", " << y_pos << ", " << z_pos << ", " << \
         electron_position.at(array_index) << ", " << electron_position.at(array_index+1) << ", " << electron_position.at(array_index+2) << ", " << electron << ", " << array_index << ", " << dt << std::endl;
      
        if (x_pos < 0.0) {x_pos += lattice_width; x_flux--;}
        else if (x_pos > lattice_width) {x_pos -= lattice_width; x_flux++;}

        if (y_pos < 0.0) {y_pos += lattice_depth; y_flux--;}
        else if (y_pos > lattice_depth) {y_pos -= lattice_depth; y_flux++;}

        if (z_pos < 0.0) {z_pos += lattice_height; z_flux--;}
        else if (z_pos > lattice_height) {z_pos -= lattice_height; z_flux++;}

        //   if (x_pos != x_pos || y_pos != y_pos || z_pos != z_pos) std::cout << "periodic boundary " << x_pos << ", " << y_pos << ", " << z_pos << ", " << \
              electron_position.at(array_index) << ", " << electron_position.at(array_index+1) << ", " << electron_position.at(array_index+2) << ", " << electron << ", " << array_index << ", " << dt << std::endl;
      
        electron_position.at(array_index)   = x_pos;
        electron_position.at(array_index+1) = y_pos;
        electron_position.at(array_index+2) = z_pos;
        }
    
      if(current_time_step % half_int_var == 0) {

           int x_cell = int(floor(x_pos / x_step_size));

           int y_cell = int(floor(y_pos / y_step_size));
    
           int z_cell = int(floor(z_pos / z_step_size));

         //  if(x_cell >= x_omp_cells || y_cell >= y_omp_cells || z_cell >= z_omp_cells) {std::cout << "cell sorting for ee integration exceeds bounds" << \
      // x_cell << ", " << y_cell << ", " << z_cell << std::endl;
      // x_cell --;
      // y_cell --;
      // z_cell --; }
      //     if(x_cell < 0 || y_cell < 0 || z_cell < 0) {std::cout << "cell sorting for ee integration less than zero" << \
      // x_cell << ", " << y_cell << ", " << z_cell << std::endl;
      // x_cell = 0;
      // y_cell = 0;
      // z_cell = 0; }
          const unsigned int omp_cell = lattice_cell_coordinate.at(x_cell).at(y_cell).at(z_cell);
          
          if ((abs(x_cell - cell_lattice_coordinate.at(cell).at(0)) > 1 &&\
              abs(x_cell - cell_lattice_coordinate.at(cell).at(0)) <= (x_omp_cells-2)) ||
             (abs(y_cell - cell_lattice_coordinate.at(cell).at(1)) > 1  &&\
              abs(y_cell - cell_lattice_coordinate.at(cell).at(1)) <= (y_omp_cells-2)) ||\
             (abs(z_cell - cell_lattice_coordinate.at(cell).at(2)) > 1  &&\
              abs(z_cell - cell_lattice_coordinate.at(cell).at(2)) <= (z_omp_cells-2) )) {
            #pragma omp critical(halo)
            {
            escaping_electrons.at(escaping_electrons.at(0)) = electron;
            escaping_electrons.at(0)++;
            }
          } else if (omp_cell != cell) {
              #pragma omp critical(shell)
              {
              cell_integration_lists.at(omp_cell).at(cell_integration_lists.at(omp_cell).at(0)) = electron;
              cell_integration_lists.at(omp_cell).at(0)++;
              }
          } else {
            cell_integration_lists.at(cell).at(cell_integration_lists.at(cell).at(0)) = electron;
            cell_integration_lists.at(cell).at(0)++;
           // if(cell_integration_lists.at(cell).at(0) >= cell_integration_lists.at(cell).size() - 1) std::cerr << "too big? " << cell_integration_lists.at(cell).at(0) << " > " << cell_integration_lists.at(cell).size() << std::endl;
           }
      }
    }
  }
  }
      if(current_time_step % half_int_var == 0) {
    const unsigned int size = escaping_electrons.at(0);
   // std::cout << "escaping electrons " << size-1 << ", " << escaping_electrons.size() << std::endl;
    
    for (int e = 1; e < size; e++) {
      const unsigned int electron = escaping_electrons.at(e);
      const unsigned int array_index = 3*electron;
      int x_cell = int(floor(electron_position.at(array_index) / x_step_size));

      int y_cell = int(floor(electron_position.at(array_index+1) / y_step_size));
    
      int z_cell = int(floor(electron_position.at(array_index+2) / z_step_size));

      //      if(x_cell >= x_omp_cells || y_cell >= y_omp_cells || z_cell >= z_omp_cells) {std::cout << "cell sorting for ee integration exceeds bounds" << \
      // x_cell << ", " << y_cell << ", " << z_cell << std::endl;
      // x_cell --;
      // y_cell --;
      // z_cell --; }
      //     if(x_cell < 0 || y_cell < 0 || z_cell < 0) {std::cout << "cell sorting for ee integration less than zero" << \
      // x_cell << ", " << y_cell << ", " << z_cell << std::endl;
      // x_cell = 0;
      // y_cell = 0;
      // z_cell = 0; }
      
        const int omp_cell = lattice_cell_coordinate.at(x_cell).at(y_cell).at(z_cell);
        cell_integration_lists.at(omp_cell).at(cell_integration_lists.at(omp_cell).at(0)) = electron;
        cell_integration_lists.at(omp_cell).at(0)++;

      //  escaping_electrons.at(e) = 0;
    }
    escaping_electrons.at(0) = 1;
  }
}

void update_dynamics() {
  
  //  int array_index;
    const static double photon_energy = sim::applied_voltage*constants::eV_to_AJ; //AJ/hv
    //AJ/fs/nm**2 -> .at(1e-3volume(A**3)/(AJ/hv)) hv/fs 
    const static double photon_rate = 1e-2*power_density*lattice_width*lattice_depth/photon_energy; //hv/fs
    double photons_at_dt = 0.0; //hv*dt
    double pump = 0.0; // AJ / fs
    double external_potential = 0.0; //AJ/e-
    const static double sigma = dt * 0.1;
    unsigned int count = 0;
   
     std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> test_int(0, conduction_electrons -1);
    std::uniform_real_distribution<double> test_uniform;
    if(!equilibrium_step) {
      if(heat_pulse_sim) {
        //hv(dt)/fs
        photons_at_dt = round(photon_rate*dt*exp(-0.5*sigma*sigma*((double(current_time_step) - ((40.0/dt)+sim::equilibration_time))*(double(current_time_step) - ((40.0 / dt)+sim::equilibration_time))))); // AJ/fs/nm**3
        pump = 1e3*photons_at_dt*photon_energy/(dt*lattice_depth*lattice_height*lattice_width); //AJ/fs/nm**3
        external_potential = photon_energy; //AJ/hv     ;//1e27*pump*dt/ n_f; // AJ / particle
     
        TTMe = d_TTMe;
        TTMp = d_TTMp;
          //AJ/fs/K -> g(T-T)=C/t
        d_TTMe = ((G*(TTMp - TTMe)+pump)*dt*e_heat_capacity_i*300.0/TTMe) + TTMe;
       
        d_TTMp = ( G*(TTMe - TTMp)      *dt*a_heat_capacity_i)      + TTMp;
      }
    }
      omp_set_dynamic(0);
       omp_set_num_threads(omp_threads);
    #pragma omp parallel for schedule(dynamic, 10) 
    for (int e = 0; e < conduction_electrons; e++) {
      const unsigned int array_index = 3*e;        
      
      if(current_time_step % full_int_var == 0) e_e_coulomb(e, array_index);
      else if (current_time_step % half_int_var == 0 && electron_transport_list.at(e)) e_e_coulomb(e, array_index); 
      else neighbor_e_e_coulomb(e, array_index);
      //replace neighbor ee with step delayed approximation for ee scattering inside DoS calculations for scattering

      if((photons_at_dt > count) && (test_int(gen) % conduction_electrons < photons_at_dt) \
      && electron_ea_scattering_list.at(e).at(0) != 1) {
        #pragma omp atomic
        count++;

        electron_thermal_field(e, array_index, external_potential);
      }

      //  if(!equilibrium_step) electron_applied_voltage(e, array_index, external_potential);
  if(!equilibrium_step) ea_scattering(e, array_index);
   } 
  //  TEKE += external_potential;
      // std::cout << " nearest neighbor complete. ee step" << std::endl;
    // std::cout << current_time_step << "; integration and ea scattering updated " << std::endl;
    ee_scattering();
  //    std::cout << current_time_step << "; ee scattering updated " << std::endl;

  Tp = Tp +  a_heat_capacity_i*1e-27*TLE *n_f/lattice_atoms;
  if(Te > 1.0) Te = Te + (e_heat_capacity_i*1e-27*TEKE*n_f*300.0/conduction_electrons/Te) + (e_heat_capacity_i*pump*dt*300.0/Te);
  else         Te = Te + (e_heat_capacity_i*1e-27*TEKE*n_f*300.0/conduction_electrons)    + (e_heat_capacity_i*pump*dt*300.0);
        
        if (err::check) std::cout << "reset scattering." << std::endl;
    
    #pragma simd
    for(int e = 0; e < conduction_electrons; e++) {
      electron_ee_scattering_list.at(e).at(0) = 0;
    }   
}

void electron_thermal_field(const int e, const int array_index, const double EKE) {
        
        //  old_vel += EKE;
   // external_interaction_list.at(e) = false;
    // int array_index_y = array_index + 1;
    // int array_index_z = array_index + 2;
   // double x_vel = electron_velocity.at(array_index);//   + ((electron_force.at(array_index)   + new_electron_force.at(array_index))   * dt  * constants::K_A / 2); 
    //double y_vel = electron_velocity.at(array_index+1);// + ((electron_force.at(array_index_y) + new_electron_force.at(array_index_y)) * dt  * constants::K_A / 2);
    //double z_vel = electron_velocity.at(array_index+2);// + ((electron_force.at(array_index_z) + new_electron_force.at(array_index_z)) * dt  * constants::K_A / 2);
 //   double vel = sqrt((x_vel*x_vel)+(y_vel*y_vel)+(z_vel*z_vel));
  //  double theta = atan(y_vel / x_vel);
   // double phi = acos(z_vel / vel);
    //if(x_vel < 0.0) theta += M_PI;
    const double theta = omp_uniform_random.at(omp_get_thread_num())() * 2.0 * M_PI;
    const double phi   = omp_uniform_random.at(omp_get_thread_num())() * M_PI; 

    electron_potential.at(e) += EKE;
    
    const double vel = sqrt(2.0*electron_potential.at(e)*constants::m_e_r_i);
        if(electron_potential.at(e) < 0.9817*E_f_A) std::cout << electron_potential.at(e) << ", " << EKE << ", " << vel << std::endl;
    electron_velocity.at(array_index)   = vel*cos(theta)*sin(phi);
    electron_velocity.at(array_index+1) = vel*sin(theta)*sin(phi);
    electron_velocity.at(array_index+2) = vel*cos(phi);

  //  electron_ee_scattering_list.at(e).at(0) = 1;
    electron_ea_scattering_list.at(e).at(0) = 1;
}
/*
// void e_a_coulomb(const int& e, const int& array_index){

//     double x_distance;
//     double y_distance;
//     double z_distance;
   
//     double length;
 
//     int array_index_a, nearest_electron_count = 1;
//     int count = 2;

//   for (int a = 0; a < lattice_atoms; a++) {

//         array_index_a = 3*a;
//         x_distance = electron_position.at(array_index)     - atom_position.at(array_index_a);
//         y_distance = electron_position.at(array_index + 1) - atom_position.at(array_index_a + 1);
//         z_distance = electron_position.at(array_index + 2) - atom_position.at(array_index_a + 2); 
       
//         if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
//         else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

//         if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
//         else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

//         if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
//         else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

//         length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
        
//         if(length > e_a_neighbor_cutoff) continue;

//         atomic_nearest_electron_list.at(e).at(nearest_electron_count) = array_index_a;
//         nearest_electron_count++;

//         if (length > e_a_coulomb_cutoff) continue;
     
//         if(mean_radius.at(2*e) > length) mean_radius.at(2*e) = length;

//         if(ea_coupling) {
//             electron_ea_scattering_list.at(e).at(count) = a;
//             count++;
//         }
//   }

//  atomic_nearest_electron_list.at(e).at(0) = nearest_electron_count;
//  electron_ea_scattering_list.at(e).at(1) = count;
 
// }

// void neighbor_e_a_coulomb(const int& e, const int& array_index){
                     
//     double x_distance;
//     double y_distance;
//     double z_distance;
    
//     double length;
//     int size = atomic_nearest_electron_list.at(e).at(0);
//     int count = 2;
  
   
//     for (int a = 1; a < size; a++) {
      
//         int array_index_a = atomic_nearest_electron_list.at(e).at(a);
        
//         x_distance = electron_position.at(array_index)     - atom_position.at(array_index_a);
//         y_distance = electron_position.at(array_index + 1) - atom_position.at(array_index_a + 1);
//         z_distance = electron_position.at(array_index + 2) - atom_position.at(array_index_a + 2); 
      
//         if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
//         else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

//         if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
//         else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

//         if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
//         else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

//       length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
   
//       if (length > e_a_coulomb_cutoff) continue;
     
//       if(mean_radius.at(2*e) > length) mean_radius.at(2*e) = length;
      
//       if(ea_coupling) {
//         electron_ea_scattering_list.at(e).at(count) = array_index_a / 3;
//         count++;
//       }
//     }
//    electron_ea_scattering_list.at(e).at(1) = count;
// }
*/
void e_e_coulomb(const int e, const int array_index) {

     int x_cell = int(floor(electron_position.at(array_index) / x_step_size));
      //if (x_cell < 0 || x_cell > x_omp_cells) std::cout << electron_position.at(array_index) << ", " << x_step_size << ", " << floor(electron_position.at(array_index) / x_step_size) << std::endl;

     int y_cell = int(floor(electron_position.at(array_index+1) / y_step_size));
      //if (y_cell < 0 || y_cell > y_omp_cells) std::cout << electron_position.at(array_index+1) << ", " << y_step_size << ", " << floor(electron_position.at(array_index+1) / y_step_size) << std::endl;
    
     int z_cell = int(floor(electron_position.at(array_index+2) / z_step_size));
      //if (z_cell < 0 || z_cell > z_omp_cells) std::cout << electron_position.at(array_index+2) << ", " << z_step_size << ", " << floor(electron_position.at(array_index+2) / z_step_size) << std::endl;
    //std::cout << x_cell << ", " << y_cell << ", " << z_cell << ", " << electron_position.at(array_index) << ", " << electron_position.at(array_index+1) << electron_position.at(array_index+2)  << std::endl;
      // if(x_cell >= x_omp_cells || y_cell >= y_omp_cells || z_cell >= z_omp_cells) {std::cout << "cell sorting for ee integration exceeds bounds" << \
      // x_cell << ", " << y_cell << ", " << z_cell << std::endl;
      // x_cell --;
      // y_cell --;
      // z_cell --; }
      // if(x_cell < 0 || y_cell < 0 || z_cell < 0) {std::cout << "cell sorting for ee integration less than zero" << \
      // x_cell << ", " << y_cell << ", " << z_cell << std::endl;
      // x_cell = 0;
      // y_cell = 0;
      // z_cell = 0; }
  ///  std::cout << lattice_cell_coordinate.at(0).at(0).at(0) << std::endl;
    const unsigned int cell = lattice_cell_coordinate.at(x_cell).at(y_cell).at(z_cell);
    unsigned int ee_dos_count = 1;
    unsigned int ee_integration_count = 1;
    unsigned int ee_scattering_list = 2;
   
    for(int s = 0; s < 27; s++) {
      const int int_cell = cell_nearest_neighbor_list.at(cell).at(s);
      const int size = cell_integration_lists.at(int_cell).at(0);
    //  std::cout << size << std::endl;
      //if(e==0) std::cout << cell << ", " << int_cell << ", " << size << std::endl;
      for(int i = 1; i < size; i++) {
        if(e == cell_integration_lists.at(int_cell).at(i)) continue;
        const int array_index_i = 3*cell_integration_lists.at(int_cell).at(i);

        double x_distance = electron_position.at(array_index)   - electron_position.at(array_index_i);
        double y_distance = electron_position.at(array_index+1) - electron_position.at(array_index_i + 1);
        double z_distance = electron_position.at(array_index+2) - electron_position.at(array_index_i + 2); 
    
        if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance += lattice_width;
        else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance -= lattice_width;

        if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance += lattice_depth;
        else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance -= lattice_depth;

        if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance += lattice_height;
        else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance -= lattice_height;
        
        const double length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
            
        if(length == 0.0) continue;
       
            if (length != length) std::cout << x_distance << ", " << y_distance << ", " << z_distance <<  ", " << \
        electron_position.at(array_index) << ", " << electron_position.at(array_index_i) << ", " << electron_position.at(array_index+1) << ", " << electron_position.at(array_index_i + 1) << ", " <<\
        electron_position.at(array_index+2) << ", " <<  electron_position.at(array_index_i + 2) << std::endl;
        
        if(length > e_e_integration_cutoff) continue;
        electron_integration_list.at(e).at(ee_integration_count) = array_index_i;
        //  if(ee_integration_count >= electron_integration_list.at(e).size() - 2) {std::cout << e << ", " << ee_integration_count << " > " << electron_integration_list.at(e).size() << ", " << length << ", " << electron_potential.at(e)  << std::endl;
        //   break; }
        ee_integration_count++;

        if(length > e_e_neighbor_cutoff) continue;
        electron_nearest_electron_list.at(e).at(ee_dos_count) = array_index_i/3;
        ee_dos_count++;
       //  std::cout << e << ", " << i << ", " << length << std::endl;
    
        //  if(ee_dos_count >= electron_nearest_electron_list.at(e).size() - 2) {std::cout << e << ", " << ee_dos_count << " > " << electron_nearest_electron_list.at(e).size() << ", " << length << ", " << electron_potential.at(e)  << std::endl;
        //   break; }

        if (length > e_e_coulomb_cutoff) continue; 
        electron_ee_scattering_list.at(e).at(ee_scattering_list) = array_index_i/3;
      // if(ee_scattering_list >= electron_ee_scattering_list.at(e).size() - 2) {std::cout << e << ", " << ee_scattering_list << " > " << electron_ee_scattering_list.at(e).size() << ", " << length << ", " << electron_potential.at(e)  << std::endl;
      //     break; }
        ee_scattering_list++;
      //   std::cout << e << ", " << i << ", " << length << std::endl;
    }
    }
    electron_integration_list.at(e).at(0) = ee_integration_count;
    electron_nearest_electron_list.at(e).at(0) = ee_dos_count;
    electron_ee_scattering_list.at(e).at(1) = ee_scattering_list;
   // if(e % 100000 == 0 ) std::cout << ee_integration_count << ", " << ee_dos_count << ", " << ee_scattering_list << std::endl;
  }
void neighbor_e_e_coulomb(const int e, const int array_index) {
    
    const unsigned int size = electron_integration_list.at(e).at(0);
    unsigned int neighbor_count = 1;
    unsigned int scattering_count = 2;
 
   // std::cout << size << std::endl;
    for (int i = 1; i < size; i++) {
        
        const unsigned int array_index_i = electron_integration_list.at(e).at(i);

        double x_distance = electron_position.at(array_index)   - electron_position.at(array_index_i);
        double y_distance = electron_position.at(array_index+1) - electron_position.at(array_index_i + 1);
        double z_distance = electron_position.at(array_index+2) - electron_position.at(array_index_i + 2); 
    
        if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance = x_distance + lattice_width;
        else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance = x_distance - lattice_width;

        if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance = y_distance + lattice_depth;
        else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance = y_distance - lattice_depth;

        if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance = z_distance + lattice_height;
        else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance = z_distance - lattice_height;
        
        const double length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
        if(length == 0.0) continue;
        
       if (length != length) std::cout << x_distance << ", " << y_distance << ", " << z_distance <<  ", " << \
        electron_position.at(array_index) << ", " << electron_position.at(array_index_i) << ", " << electron_position.at(array_index+1) << ", " << electron_position.at(array_index_i + 1) << ", " <<\
        electron_position.at(array_index+2) << ", " <<  electron_position.at(array_index_i + 2) << std::endl;
        
        if(length > e_e_neighbor_cutoff) continue;
        electron_nearest_electron_list.at(e).at(neighbor_count) = array_index_i/3;
        neighbor_count++;
       //  std::cout << e << ", " << i << ", " << length << std::endl;
    
        //  if(neighbor_count >= electron_nearest_electron_list.at(e).size() - 2) {std::cout << e << ", " << neighbor_count << " > " << electron_nearest_electron_list.at(e).size() << ", " << length << ", " << electron_potential.at(e)  << std::endl;
        //    break; }

        if (length > e_e_coulomb_cutoff) continue; 
        electron_ee_scattering_list.at(e).at(scattering_count) = array_index_i/3;
        scattering_count++;
        // if(scattering_count >= electron_ee_scattering_list.at(e).size() - 1) {std::cout << e << ", " << scattering_count << " > " << electron_ee_scattering_list.at(e).size() << ", " << length << ", " << electron_potential.at(e)  << std::endl;
        //    std::cout << x_distance << ", " << y_distance << ", " << z_distance << std::endl;
        //    break; }
      //   std::cout << e << ", " << i << ", " << length << std::endl;
    }

    electron_nearest_electron_list.at(e).at(0) = neighbor_count;
    electron_ee_scattering_list.at(e).at(1) = scattering_count;
}
/*
// void a_a_coulomb(const int a, const int array_index, \
//                         double& a_x_force, double& a_y_force, double& a_z_force, double& LPE) {
   
//     double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
//  //   int neighbor_count = 1;

//         x_distance = atom_position.at(array_index);
//         y_distance = atom_position.at(array_index+1);
//         z_distance = atom_position.at(array_index+2);
      
//                // if (err::check)  std::cout << "Calculating atomic interactions" << std::endl;

//        // neighbor_count = 1;

//         for (int i = 0; i < lattice_atoms; i++) {
//             if (i == a) continue; //no self repulsion

//             x_distance -= atom_position.at(array_index);
//             y_distance -= atom_position.at(array_index + 1);
//             z_distance -= atom_position.at(array_index + 2);

//             length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
           
//             if(length > 6) continue;

//             //length = sqrt(length);
//          //   force = -2000*(length - 2);
//            // PE = 1000*(length - 2)*(length - 2);
//         //   if(e == 0) std::cout << force << ", " << length <<  std::endl;
 
//             //phi   = acos(z_distance / length);
//             //theta = atan(y_distance / x_distance);
//             if (x_distance < 0) theta += M_PI;

//             atom_position.at(array_index)     += force * cos(theta)*sin(phi) / atomic_mass;
//             atom_position.at(array_index + 1) += force * sin(theta)*sin(phi) / atomic_mass;
//             atom_position.at(array_index + 2) += force * cos(phi) / atomic_mass;
//         }
//    // atomic_nearest_atom_list.at(a).at(0) = neighbor_count;
//     LPE += PE/2;
//    // if(a == 100) std::cout << neighbor_count << std::endl;
// }

// void neighbor_a_a_coulomb(const int a, const int array_index, \
//                         double& a_x_force, double& a_y_force, double& a_z_force, double& LPE) {
    
//     int array_index_i;
//     double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
//     double size = atomic_nearest_atom_list.at(a).at(0);

//         x_distance = atom_position.at(array_index);
//         y_distance = atom_position.at(array_index+1);
//         z_distance = atom_position.at(array_index+2);

//         for (int i = 1; i < size; i++) {
//            // if (i == a) continue; //no self repulsion
           
//             array_index_i = atomic_nearest_atom_list.at(a).at(i);

//             x_distance -= atom_position.at(array_index_i);
//             y_distance -= atom_position.at(array_index_i + 1);
//             z_distance -= atom_position.at(array_index_i + 2); 

//             if (x_distance < -30)     x_distance = x_distance + 40;
//             else if (x_distance > 30) x_distance = x_distance - 40;
//             if (y_distance < -30)     y_distance = y_distance + 40;
//             else if (y_distance > 30) y_distance = y_distance - 40;
//             if (z_distance <  -30)    z_distance = z_distance + 40;
//             else if (z_distance > 30) z_distance = z_distance - 40;

//             length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
             
//             // if (length > a_a_coulomb_cutoff) continue; 

//             length = sqrt(length);
//             force = 4*(-2*exp(4 - 2*length) + 2*exp(length + 2));
//             PE = 4*(exp(4 - 2*length) - (2*exp(2 -length)) + 1);
            
//             phi   = acos(z_distance / length);
//             theta = atan(y_distance / x_distance);
//             if (x_distance < 0) theta += M_PI;

//             a_x_force += force * cos(theta)*sin(phi) / atomic_mass;
//             a_y_force += force * sin(theta)*sin(phi) / atomic_mass;
//             a_z_force += force * cos(phi) / atomic_mass;

//         }
//     LPE += PE/2;
//    // if(a == 100) std::cout << size << std::endl;
// }
*/
double electron_applied_voltage(const int& e, const int& array_index, double& external_potential) {
  

  const static double vel = 5e-11*dt*applied_voltage*constants::e_A*constants::m_e_r_i; //V/m * q = F_x/kg = 0.5*a_x*dt = d_v_x
  double e_energy = electron_potential.at(e);
  double deltaE = 0.5*constants::m_e_r*(((vel+electron_velocity.at(array_index))*(vel+electron_velocity.at(array_index)))+(electron_velocity.at(array_index+1)*electron_velocity.at(array_index+1))+(electron_velocity.at(array_index+2)*electron_velocity.at(array_index+2))) - e_energy;

  double DoS1 = electron_nearest_electron_list.at(e).at(0) - 1;
if(e_energy + deltaE > 0.9817*E_f_A) {
  for (int i = 1; i < electron_nearest_electron_list.at(e).at(0); i++) {
    if((electron_potential.at(electron_nearest_electron_list.at(e).at(i)) < e_energy + fmax(1.2*deltaE, 0.8*deltaE)) \
    && (electron_potential.at(electron_nearest_electron_list.at(e).at(i)) > e_energy + fmin(0.8*deltaE, 1.2*deltaE))) {
      DoS1 -= 1.0;
    }
  }
  DoS1 /= double(electron_nearest_electron_list.at(e).at(0) - 1.0);
  if(omp_uniform_random.at(omp_get_thread_num())() > exp(-3.0*DoS1)) {
    electron_velocity.at(array_index) += vel;
    external_potential += deltaE;
    electron_potential.at(e) += deltaE;
  }
} 
  return 0.0;
}
/*
// void aa_scattering() {
//         if (err::check) std::cout << "aa_scattering." << std::endl;
//   double scattering_prob;
  
//   const static double aa_rate = -10.0*dt/sim::ee_coupling_strength;
//     for(int e = 0; e < conduction_electrons; e++) {
  
//       if(atomic_nearest_atom_list.at(e).at(0)) continue;
      
//       double deltaE = atom_potential.at(e) - E_f_A; - 3.0*constants::kB_r*Tp;

//      // if(deltaE < 8.0) scattering_prob = 8.0;
//        scattering_prob = deltaE*deltaE;

//       if(omp_uniform_random.at(omp_get_thread_num())() > exp(aa_rate*scattering_prob)) {
        
//         deltaE *= 0.5-0.05*exp(-2.0*deltaE/(E_f_A-3.0*constants::kB_r*Tp));
//         int size = atomic_nearest_atom_list.at(e).at(1) - 2;
//         size = 2+ (omp_int_random.at(omp_get_thread_num())() % size);
//       //  if(size >= atomic_nearest_atom_list.at(e).size() -1) std::cout << "bus error: " << size << " < " << atomic_nearest_atom_list.at(e).size() << std::endl; 
//         int a = atomic_nearest_atom_list.at(e).at(size);
      
//        // if(a >= conduction_electrons -1) std::cout << "bus error: " << a << " < " << conduction_electrons-1 << std::endl; 
//         if(atomic_nearest_atom_list.at(a).at(0)) continue;

//         deltaE *= 0.5;
//         atom_potential.at(e) -= deltaE;
//         atom_potential.at(a) += deltaE;  

//         atomic_nearest_atom_list.at(e).at(0) = 1;
//         atomic_nearest_atom_list.at(a).at(0) = 1;
        
//         a_a_scattering_count++;
//       }
//     }
//     if (err::check) std::cout << "aa_scattering done." << std::endl;
// }
*/
void ea_scattering(const int e, const int array_index) {

      std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> test_int(0, conduction_electrons -1);
    std::uniform_real_distribution<double> test_uniform;
    double scattering_velocity = electron_potential.at(e);
    double deltaE = sqrt(phonon_energy*E_f_A);
  
    if(test_uniform(gen) > exp(ea_rate/scattering_velocity)) {
     
      const unsigned int size = electron_nearest_electron_list.at(e).at(0);
      int DoS1 = size - 1;
      int DoS2 = size - 1;
      for (int i = 1; i < size ; i++) {
        if((electron_potential.at(electron_nearest_electron_list.at(e).at(i)) < scattering_velocity + 1.2*deltaE) \
        && (electron_potential.at(electron_nearest_electron_list.at(e).at(i)) > scattering_velocity + 0.8*deltaE)) {
          DoS1 -= 1.0;
        }
        else if ((electron_potential.at(electron_nearest_electron_list.at(e).at(i)) < 1.2*scattering_velocity) \
        && (electron_potential.at(electron_nearest_electron_list.at(e).at(i)) > 0.8*scattering_velocity)) {
          DoS2 -= 1.0;
        }
      }
      if(DoS1 < DoS2) deltaE *= -1.0;
      if(Tp < Te) deltaE *= -1.0;
      if((scattering_velocity + deltaE) < 0.9817*E_f_A) deltaE = 0.0;

      electron_potential.at(e) += deltaE;

      if(electron_potential.at(e) > 0.99*E_f_A) electron_transport_list.at(e) = true;
      else electron_transport_list.at(e) = false;

      const double unit_vec = sqrt(2.0*scattering_velocity*constants::m_e_r_i);
      const double theta = 2.0*M_PI*test_uniform(gen);
      const double phi   = M_PI*test_uniform(gen) ;
     // if(electron_velocity.at(array_index+2) < 0.0) theta += M_PI;

      scattering_velocity = sqrt(2.0*electron_potential.at(e)*constants::m_e_r_i);
      electron_velocity.at(array_index)   = scattering_velocity * cos(theta)*sin(phi);
      electron_velocity.at(array_index+1) = scattering_velocity * sin(theta)*sin(phi);
      electron_velocity.at(array_index+2) = scattering_velocity * cos(phi); 
    
     // electron_ee_scattering_list.at(e).at(0) = 1;

      #pragma omp critical(eascattering)
      {
      TEKE += deltaE;
      TLE -= deltaE;
      e_a_scattering_count++;
      }
    }  
}

void ee_scattering() {
         if (err::check) std::cout << "ee_scattering." << std::endl;
  omp_set_dynamic(0);
       omp_set_num_threads(omp_threads);

  std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> test_int(0, conduction_electrons -1);
    std::uniform_real_distribution<double> test_uniform;
  #pragma omp parallel reduction(+:e_e_scattering_count) 
  {
  for(int l = 0; l < cells_per_thread; l++) {
    const unsigned int cell = lattice_cells_per_omp.at(omp_get_thread_num()).at(l);
    const unsigned int size = cell_integration_lists.at(cell).at(0);
  //  if(size >= cell_integration_lists.at(cell).size()) std::cout << size << ", " << cell_integration_lists.at(cell).at(0) << ", " << cell << ", " << lattice_cells_per_omp.at(omp_get_thread_num()).at(l) << std::endl;
    
  //  #pragma omp critical 
   // if(!equilibrium_step) std::cout << omp_get_thread_num() << ", " << cell << ", " << size << ", " << l << std::endl;

    for(int e = 1; e < size; e++) {
      const unsigned int electron = cell_integration_lists.at(cell).at(e);
      if (electron_ee_scattering_list.at(electron).at(0) == 1) continue;

        // if(electron != electron) std::cout << electron << ", " << cell << ", " << size << ", " << l << std::endl;

      const uint32_t c_size = electron_ee_scattering_list.at(electron).at(1) - 2;

      if(c_size <= 3) continue;

      // #pragma omp critical
      // if(c_size != c_size) std::cout << "c_size selection " << c_size << ", " << electron << ", " << cell << ", " << size << ", " << l << std::endl;

      unsigned int electron_collision;
      int found = 0;
      for(uint32_t a = 0; a < c_size; a++) {
        electron_collision = 2 + (test_int(gen) % c_size);

            // #pragma omp critical
            // if(electron_collision != electron_collision) std::cout << "random collision selection" << electron << ", " << cell << ", " << size << ", " << l << std::endl;
            
            // if(electron_collision >= electron_ee_scattering_list.at(electron).at(1)) std::cout << "random collision error " << electron_collision << ", " << electron_ee_scattering_list.at(electron).at(0) << std::endl;
        
        electron_collision = electron_ee_scattering_list.at(electron).at(electron_collision);
         if(electron_ee_scattering_list.at(electron_collision).at(0) == 1) continue;
         if(!electron_transport_list.at(electron)) {
          if(!electron_transport_list.at(electron_collision)) continue;
         }
         //  if(electron_collision != electron_collision) std::cout << electron << ", " << cell << ", " << size << ", " << l << std::endl;
        const unsigned int array_index_i = 3*electron_collision;
         int x_cell = int(floor(electron_position.at(array_index_i)  / x_step_size));
         int y_cell = int(floor(electron_position.at(array_index_i+1)/ y_step_size));
         int z_cell = int(floor(electron_position.at(array_index_i+2)/ z_step_size));

      //   if(x_cell >= x_omp_cells || y_cell >= y_omp_cells || z_cell >= z_omp_cells) {std::cout << "cell sorting for ee integration exceeds bounds" << \
      // x_cell << ", " << y_cell << ", " << z_cell << std::endl;
      // x_cell --;
      // y_cell --;
      // z_cell --; }
      // if(x_cell < 0 || y_cell < 0 || z_cell < 0) {std::cout << "cell sorting for ee integration less than zero" << \
      // x_cell << ", " << y_cell << ", " << z_cell << std::endl;
      // x_cell = 0;
      // y_cell = 0;
      // z_cell = 0; }

        const int cell_i = lattice_cell_coordinate.at(x_cell).at(y_cell).at(z_cell);

        if (cell_i != cell) continue;
       
        found = 1;
        break;
      }
      if (found == 0) continue;
    //if(electron == electron_collision) { std::cout << "uh oh " << electron << ", " << electron_collision << ", " << cell << ", " << l << std::endl; continue;}
      electron_ee_scattering_list.at(electron).at(0) = 1;
      electron_ee_scattering_list.at(electron_collision).at(0) = 1;

      const double e_energy = electron_potential.at(electron);
      const double d_e_energy = electron_potential.at(electron_collision);
      if(equilibrium_step) {
        double deltaE = test_uniform(gen)*(e_energy - d_e_energy);
    //   if (deltaE != deltaE) std::cout  << deltaE << ", " << e_energy << ", " << d_e_energy << std::endl;
    
        if(test_uniform(gen) < 0.5*exp(abs(deltaE)/(-0.5))) deltaE *= -1.0;
     //    if (deltaE != deltaE) std::cout  << deltaE << ", " << e_energy << ", " << d_e_energy << std::endl;

        int DoS1 = electron_nearest_electron_list.at(electron).at(0) - 1;
        int DoS2 = electron_nearest_electron_list.at(electron_collision).at(0) - 1;
 
        if((e_energy - deltaE < 0.9817*E_f_A) || (d_e_energy + deltaE < 0.9817*E_f_A)) DoS2 = 0;
        else {
          for (int e = 1; e < electron_nearest_electron_list.at(electron).at(0); ++e) {
            if((electron_potential.at(electron_nearest_electron_list.at(electron).at(e)) < (e_energy - fmin(0.8*deltaE, 1.2*deltaE))) \
            && (electron_potential.at(electron_nearest_electron_list.at(electron).at(e)) > (e_energy - fmax(0.8*deltaE, 1.2*deltaE)))) {
              DoS1 -= 1;
            //  if (DoS1 < 0) std::cout << "inelastic DoS1: " << DoS1 << ", " <<  electron_nearest_electron_list.at(electron).at(0) << ", " << e << std::endl;
            }
          }
      
          for (int i = 1; i < electron_nearest_electron_list.at(electron_collision).at(0); ++i) {
            if((electron_potential.at(electron_nearest_electron_list.at(electron_collision).at(i)) < (d_e_energy + fmax(0.8*deltaE, 1.2*deltaE))) \
            && (electron_potential.at(electron_nearest_electron_list.at(electron_collision).at(i)) > (d_e_energy + fmin(0.8*deltaE, 1.2*deltaE)))) {
              DoS2 -= 1;
              //#pragma omp critical
            //  if (DoS2 < 0) std::cout << "inelastic DoS2: " << DoS2 << ", " <<  electron_nearest_electron_list.at(electron_collision).at(0)-1 << ", " << i << std::endl;
            }
          }
        }

        if(test_uniform(gen) > exp(ee_rate*DoS1*DoS2/(0.6666+(deltaE*deltaE)))) {
       
          electron_potential.at(electron) -= deltaE;
          electron_potential.at(electron_collision)   += deltaE;
          electron_transport_list.at(electron) = true;
          electron_transport_list.at(electron_collision) = true;
          if (electron_potential.at(electron) < 0.99*E_f_A)  electron_transport_list.at(electron) = false;
          if (electron_potential.at(electron_collision) < 0.99*E_f_A) electron_transport_list.at(electron_collision) = false;
          ee_inelastic(electron, electron_collision, deltaE);
          e_e_scattering_count += 2;
    
        }
      } else e_e_scattering_count += int(ee_elestic(electron, electron_collision, e_energy, d_e_energy, test_uniform(gen)));
    }
  }
    
  }   if (err::check) std::cout << "ee_scattering done." << std::endl;
}

int ee_inelastic(const int electron, const int electron_collision, const double deltaE){

        std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> test_int(0, conduction_electrons -1);
    std::uniform_real_distribution<double> test_uniform;

      const unsigned int array_index = 3*electron;
      const unsigned int array_index_i = 3*electron_collision;
       double theta = 2.0*M_PI*test_uniform(gen);
       double phi = M_PI*test_uniform(gen);

      const double v_1 = sqrt(2.0*constants::m_e_r_i*electron_potential.at(electron));
      const double v_2 = sqrt(2.0*constants::m_e_r_i*electron_potential.at(electron_collision));

      electron_velocity.at(array_index)   = v_1*cos(theta)*sin(phi);
      electron_velocity.at(array_index+1) = v_1*sin(theta)*sin(phi);
      electron_velocity.at(array_index+2) = v_1*cos(phi);

       theta = 2.0*M_PI*test_uniform(gen);
        phi = M_PI*test_uniform(gen);
      electron_velocity.at(array_index_i)   = v_2*cos(theta)*sin(phi);
      electron_velocity.at(array_index_i+1) = v_2*sin(theta)*sin(phi);
      electron_velocity.at(array_index_i+2) = v_2*cos(phi);

    return 2;
}

int ee_elestic(const int electron, const int electron_collision, const double e_energy, const double d_e_energy, const double probability) {
   const unsigned int array_index = 3*electron;
   const unsigned int array_index_i = 3*electron_collision;

    double x_distance = electron_position.at(array_index)  - electron_position.at(array_index_i);
    double y_distance = electron_position.at(array_index+1)- electron_position.at(array_index_i+1);
    double z_distance = electron_position.at(array_index+2)- electron_position.at(array_index_i+2); 

//    if(x_distance == y_distance == z_distance == 0) std::cout << x_distance << ", " << y_distance << ", " << z_distance <<  ", " << \
            electron_position.at(array_index) << ", " << electron_position.at(array_index_i) << ", " << electron_position.at(array_index+1) << ", " << electron_position.at(array_index_i + 1) << ", " <<\
            electron_position.at(array_index+2) << ", " <<  electron_position.at(array_index_i + 2) << ", " << boundary_conditions_cutoff << ", " << lattice_depth << ", " << lattice_width << ", " << lattice_height << std::endl;

    if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance += lattice_width;
    else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance -= lattice_width;

    if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance += lattice_depth;
    else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance -= lattice_depth;
    
    if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance += lattice_height;
    else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance -= lattice_height;
        
    // #pragma omp critical
    // if(x_distance == 0 && y_distance == 0 && z_distance == 0) std::cout << "distance " << x_distance << ", " << y_distance << ", " << z_distance <<  ", " << \
            electron_position.at(array_index) << ", " << electron_position.at(array_index_i) << ", " << electron_position.at(array_index+1) << ", " << electron_position.at(array_index_i + 1) << ", " <<\
            electron_position.at(array_index+2) << ", " <<  electron_position.at(array_index_i + 2) << ", " << boundary_conditions_cutoff << ", " << lattice_depth << ", " << lattice_width << ", " << lattice_height << std::endl;
            
    const double length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);

    // #pragma omp critical
    // if(length == 0.0) std::cout << "length "<< x_distance << ", " << y_distance << ", " << z_distance <<  ", " << \
            electron_position.at(array_index) << ", " << electron_position.at(array_index_i) << ", " << electron_position.at(array_index+1) << ", " << electron_position.at(array_index_i + 1) << ", " <<\
            electron_position.at(array_index+2) << ", " <<  electron_position.at(array_index_i + 2) << ", " << boundary_conditions_cutoff << ", " << lattice_depth << ", " << lattice_width << ", " << lattice_height << std::endl;
            // #pragma omp critical
          // if (x_distance != x_distance || y_distance != y_distance ||  z_distance != z_distance) std::cout << "nan distance "<< x_distance << ", " << y_distance << ", " << z_distance <<  ", " << \
            electron_position.at(array_index) << ", " << electron_position.at(array_index_i) << ", " << electron_position.at(array_index+1) << ", " << electron_position.at(array_index_i + 1) << ", " <<\
            electron_position.at(array_index+2) << ", " <<  electron_position.at(array_index_i + 2) << ", " << boundary_conditions_cutoff << ", " << lattice_depth << ", " << lattice_width << ", " << lattice_height << std::endl;
    const double v_x_dot_product = ((electron_velocity.at(array_index)  - electron_velocity.at(array_index_i)   )*x_distance\
                                  + (electron_velocity.at(array_index+1)- electron_velocity.at(array_index_i+1) )*y_distance\
                                  + (electron_velocity.at(array_index+2)- electron_velocity.at(array_index_i+2) )*z_distance) /length;

      // #pragma omp critical
      // if(v_x_dot_product != v_x_dot_product) std::cout << "dot product " << v_x_dot_product << ", " << electron_velocity.at(array_index)  << ", " <<  electron_velocity.at(array_index_i) << ", " << x_distance << ", " <<\
                                             electron_velocity.at(array_index+1) << ", " << electron_velocity.at(array_index_i+1) << ", " << y_distance << ", " << \
                                             electron_velocity.at(array_index+2) << ", " << electron_velocity.at(array_index_i+2) << ", " << z_distance << ", " << length << std::endl;

      if(v_x_dot_product == 0.0) return 0;
    const double normalised_dot_product = v_x_dot_product;///sqrt(length*v_x_dot_product*v_x_dot_product);
    //if(v_x_dot_product/sqrt(length*v_x_dot_product*v_x_dot_product) > 1.0) std::cout << " normalised dot product " << v_x_dot_product/sqrt(length*v_x_dot_product*v_x_dot_product);

    const double d_v_x_1 = electron_velocity.at(array_index_i)   + x_distance*normalised_dot_product;
    const double d_v_y_1 = electron_velocity.at(array_index_i+1) + y_distance*normalised_dot_product;
    const double d_v_z_1 = electron_velocity.at(array_index_i+2) + z_distance*normalised_dot_product;
    const double d_v_x_2 = electron_velocity.at(array_index)   - x_distance*normalised_dot_product;
    const double d_v_y_2 = electron_velocity.at(array_index+1) - y_distance*normalised_dot_product;
    const double d_v_z_2 = electron_velocity.at(array_index+2) - z_distance*normalised_dot_product;

    const double deltaE = -1.0*(0.5*constants::m_e_r*(d_v_x_1*d_v_x_1 + d_v_y_1*d_v_y_1 + d_v_z_1*d_v_z_1) -\
                           0.5*constants::m_e_r*(d_v_x_2*d_v_x_2 + d_v_y_2*d_v_y_2 + d_v_z_2*d_v_z_2));
    // #pragma omp critical
    if(deltaE != deltaE) std::cout << "deltaE " << deltaE << ", " << normalised_dot_product << ", " << length << std::endl;
     int DoS1 = electron_nearest_electron_list.at(electron).at(0) - 1;
     int DoS2 = electron_nearest_electron_list.at(electron_collision).at(0) - 1;
 
      if((e_energy - deltaE < 0.9817*E_f_A) || (d_e_energy + deltaE < 0.9817*E_f_A)) DoS2 = 0;
      else {
      for (int i = 1; i < electron_nearest_electron_list.at(electron).at(0); i++) {
        if((electron_potential.at(electron_nearest_electron_list.at(electron).at(i)) < (e_energy - fmin(0.8*deltaE, 1.2*deltaE))) \
        && (electron_potential.at(electron_nearest_electron_list.at(electron).at(i)) > (e_energy - fmax(0.8*deltaE, 1.2*deltaE)))) {
          DoS1 -= 1.0;
        }
      }
      
      for (int i = 1; i < electron_nearest_electron_list.at(electron_collision).at(0); i++) {
        if((electron_potential.at(electron_nearest_electron_list.at(electron_collision).at(i)) < (d_e_energy + fmax(0.8*deltaE, 1.2*deltaE))) \
        && (electron_potential.at(electron_nearest_electron_list.at(electron_collision).at(i)) > (d_e_energy + fmin(0.8*deltaE, 1.2*deltaE)))) {
          DoS2 -= 1.0;
        }
      } 
      // #pragma omp critical
      // if(DoS1 < 0 || DoS2 < 0) {
      //     std::cout << "elastic DoS < 0 " << DoS1 << ", " <<  electron_nearest_electron_list.at(electron).at(0) << ", " << DoS2 << ", " << electron_nearest_electron_list.at(electron_collision).at(0) <<std::endl;
      //     DoS1 = int(fmax(0, DoS1));
      //     DoS2 = int(fmax(0, DoS2));
      //   }
      //   if(DoS1 != DoS1 || DoS2 != DoS2) {
      //     std::cout << "elastic DoS < 0 " << DoS1 << ", " <<  electron_nearest_electron_list.at(electron).at(0) << ", " << DoS2 << ", " << electron_nearest_electron_list.at(electron_collision).at(0) <<std::endl; 
      //   }
      }
      if(probability > exp(ee_rate*DoS1*DoS2/(0.6666+(deltaE*deltaE)))) {
       
        electron_potential.at(electron) -= deltaE;
        electron_potential.at(electron_collision)   += deltaE;
        electron_transport_list.at(electron) = true;
        electron_transport_list.at(electron_collision) = true;
        if (electron_potential.at(electron) < 0.99*E_f_A)  electron_transport_list.at(electron) = false;
        if (electron_potential.at(electron_collision) < 0.99*E_f_A) electron_transport_list.at(electron_collision) = false;
    
        electron_velocity.at(array_index_i)  += x_distance*normalised_dot_product;
        electron_velocity.at(array_index_i+1)+= y_distance*normalised_dot_product; 
        electron_velocity.at(array_index_i+2)+= z_distance*normalised_dot_product;
    
        electron_velocity.at(array_index)  -= x_distance*normalised_dot_product;
        electron_velocity.at(array_index+1)-= y_distance*normalised_dot_product;
        electron_velocity.at(array_index+2)-= z_distance*normalised_dot_product; 
        return 2;
      }
     return 0;
}
   
} //end CASTLE namespace

