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
    const int size = ee_dos_hist[0].size();
    #pragma omp parallel for
    for(int e = 0; e < conduction_electrons; e++) {
      electron_ee_scattering_list[e][0] = 0;
      for(int h = 0; h < size; h++) {
        ee_dos_hist[e][h] = 0;
      }
    }
    //reset integration
    CASTLE_real_time += dt;
    current_time_step++;
   
    return EXIT_SUCCESS;
}

void setup_output() {

    CASTLE_output_data = true;
}


void update_position(){

    omp_set_dynamic(0);
    omp_set_num_threads(omp_threads);
    if(current_time_step % half_int_var == 0) {
      
      old_cell_integration_lists.swap(cell_integration_lists);
      #pragma omp parallel for
      for(int c = 0; c < total_cells; c++) {
        cell_integration_lists[c][0] = 1;
      }
    }
    
  #pragma omp parallel reduction(+:x_flux,y_flux,z_flux, p_x,p_y,p_z)
  {
    const int thread = omp_get_thread_num();
  for (int l = 0 ; l < cells_per_thread; l++) {
    const unsigned int cell = lattice_cells_per_omp[thread][l];
    const unsigned int size = old_cell_integration_lists[cell][0];
   // cell_integration_lists[cell][0] = 1;

  //  if((current_time_step % half_int_var) == 0 && ( (l % ((z_omp_cells/2)-1) )== 0 || (l % (z_omp_cells/2)) == 0) ) {
      #pragma omp barrier
   // }
    
    for (int e = 1; e < size; e++) { 
      const unsigned int electron = old_cell_integration_lists[cell][e];

        const unsigned int array_index = 3*electron;
      
        double x_pos = electron_position[array_index];
        double y_pos = electron_position[array_index+1];
        double z_pos = electron_position[array_index+2];
    
     // if(electron_potential[electron] > E_f_A) {
        const double v_x = electron_velocity[array_index];
        const double v_y = electron_velocity[array_index+1] ;
        const double v_z = electron_velocity[array_index+2] ;

        p_x += v_x;
        p_y += v_y;
        p_z += v_z;
        
        x_pos += v_x * dt;// + (electron_force[array_index]   * dt * dt * constants::K_A / 2); // x superarray component
        y_pos += v_y * dt;// + (electron_force[array_index_y) * dt * dt * constants::K_A / 2); // y superarray component
        z_pos += v_z * dt;// + (electron_force[array_index_z) * dt * dt * constants::K_A / 2); // z superarray component

        if (x_pos < 0.0) {x_pos += lattice_width; x_flux--;}
        else if (x_pos > lattice_width) {x_pos -= lattice_width; x_flux++;}

        if (y_pos < 0.0) {y_pos += lattice_depth; y_flux--;}
        else if (y_pos > lattice_depth) {y_pos -= lattice_depth; y_flux++;}

        if (z_pos < 0.0) {z_pos += lattice_height; z_flux--;}
        else if (z_pos > lattice_height) {z_pos -= lattice_height; z_flux++;}

        electron_position[array_index]   = x_pos;
        electron_position[array_index+1] = y_pos;
        electron_position[array_index+2] = z_pos;
    //  }    
    
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
          const unsigned int omp_cell = lattice_cell_coordinate[x_cell][y_cell][z_cell];
          
          if ((abs(x_cell - cell_lattice_coordinate[cell][0]) > 1 &&\
              abs(x_cell - cell_lattice_coordinate[cell][0]) <= (x_omp_cells-2)) ||
             (abs(y_cell - cell_lattice_coordinate[cell][1]) > 1  &&\
              abs(y_cell - cell_lattice_coordinate[cell][1]) <= (y_omp_cells-2)) ||\
             (abs(z_cell - cell_lattice_coordinate[cell][2]) > 1  &&\
              abs(z_cell - cell_lattice_coordinate[cell][2]) <= (z_omp_cells-2) )) {
            #pragma omp critical(halo)
            {
            escaping_electrons[escaping_electrons[0]] = electron;
            escaping_electrons[0]++;
            }
          } else if (omp_cell != cell) {
              #pragma omp critical(shell)
              {
              cell_integration_lists[omp_cell][cell_integration_lists[omp_cell][0]] = electron;
              cell_integration_lists[omp_cell][0]++;
              }
          } else {
            cell_integration_lists[cell][cell_integration_lists[cell][0]] = electron;
            cell_integration_lists[cell][0]++;
           // if(cell_integration_lists[cell][0] >= cell_integration_lists[cell).size() - 1) std::cerr << "too big? " << cell_integration_lists[cell][0] << " > " << cell_integration_lists[cell).size() << std::endl;
           }
      }
    }
  }
  }
      if(current_time_step % half_int_var == 0) {
    const unsigned int size = escaping_electrons[0];
   // std::cout << "escaping electrons " << size-1 << ", " << escaping_electrons.size() << std::endl;
    
    for (int e = 1; e < size; e++) {
      const unsigned int electron = escaping_electrons[e];
      const unsigned int array_index = 3*electron;
      int x_cell = int(floor(electron_position[array_index] / x_step_size));

      int y_cell = int(floor(electron_position[array_index+1] / y_step_size));
    
      int z_cell = int(floor(electron_position[array_index+2] / z_step_size));

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
      
        const int omp_cell = lattice_cell_coordinate[x_cell][y_cell][z_cell];
        cell_integration_lists[omp_cell][cell_integration_lists[omp_cell][0]] = electron;
        cell_integration_lists[omp_cell][0]++;

      //  escaping_electrons[e] = 0;
    }
    escaping_electrons[0] = 1;
  }
}

void update_dynamics() {
  
  //  int array_index;
    const static double photon_energy = sim::applied_voltage*constants::eV_to_AJ; //AJ/hv
    //AJ/fs/nm**2 -> [1e-3volume(A**3)/(AJ/hv)) hv/fs 
    const static double photon_rate = 1e-2*power_density*lattice_width*lattice_depth/photon_energy; //hv/fs
    int photons_at_dt = 0; //hv*dt
    double pump = 0.0; // AJ / fs
    double external_potential = 0.0; //AJ/e-
    const static double sigma = dt * 0.1;
    unsigned int count = 0;

    if(!equilibrium_step) {
      if(heat_pulse_sim) {
        //hv(dt)/fs
        photons_at_dt = int(round(photon_rate*dt*exp(-0.5*sigma*sigma*((double(current_time_step) - ((10.0/dt)+sim::equilibration_time))*(double(current_time_step) - ((10.0 / dt)+sim::equilibration_time)))))); // AJ/fs/nm**3
        pump = 1e3*photons_at_dt*photon_energy/(dt*lattice_depth*lattice_height*lattice_width); //AJ/fs/nm**3
        external_potential = photon_energy; //AJ/hv     ;//1e27*pump*dt/ n_f; // AJ / particle
     //   std::cout << photons_at_dt << std::endl;
        TTMe = d_TTMe;
        TTMp = d_TTMp;
          //AJ/fs/K -> g(T-T)=C/t
        d_TTMe = ((G*(TTMp - TTMe)+pump)*dt*e_heat_capacity_i/TTMe) + TTMe;
      //  std::cout <<" TTMe: " << pump*dt*e_heat_capacity_i*300.0/TTMe << ", " << pump*dt*e_heat_capacity_i*300.0/Te;
        d_TTMp = ( G*(TTMe - TTMp)      *dt*a_heat_capacity_i)      + TTMp;
      }
    }
      omp_set_dynamic(0);
       omp_set_num_threads(omp_threads);
       std::vector<int> chosen;
       chosen.resize(photons_at_dt);
      while(count < photons_at_dt) {
        int select = int(omp_uniform_random[omp_get_thread_num()]()*2147483647) % (conduction_electrons);
        double e_energy = electron_potential[select];
        if(return_phonon_distribution((e_energy+photon_energy-E_f_A)/E_f_A, constants::kB_r*Te/E_f_A) < (1.0-1e-3)) {
          chosen[count] = select;
          count++;
        }
      } count = 0;
    #pragma omp parallel for schedule(dynamic, 10)
    for (int e = 0; e < conduction_electrons; e++) {
      const unsigned int array_index = 3*e;        
      
      if(current_time_step % half_int_var == 0 ) e_e_coulomb(e, array_index);
     // else if (current_time_step % full_int_var == 0) e_e_coulomb(e, array_index);
       else neighbor_e_e_coulomb(e, array_index);
      

      if(photons_at_dt > 0 && std::end(chosen) != std::find(chosen.begin(), chosen.end(), e)) {
        #pragma omp atomic
        count++;
      
        electron_thermal_field(e, array_index, external_potential, omp_get_thread_num());
      }
      
      //  if(!equilibrium_step) electron_applied_voltage(e, array_index, external_potential);
    if(!equilibrium_step) ea_scattering(e, array_index, omp_get_thread_num());
    }
   if(count != photons_at_dt) std::cout << photons_at_dt << ", " << count<< std::endl;
  
   ee_scattering();
  
  Tp +=  a_heat_capacity_i*1e-27*TLE *n_f/conduction_electrons;
  Te += (e_heat_capacity_i*1e-27*TEKE*n_f/conduction_electrons/Te) + (e_heat_capacity_i*pump*dt/Te);
 
        if (err::check) std::cout << "reset scattering." << std::endl;
       
}

void electron_thermal_field(const int e, const int array_index, const double EKE, const int thread) {
        
        //  old_vel += EKE;
   // external_interaction_list[e] = false;
    // int array_index_y = array_index + 1;
    // int array_index_z = array_index + 2;
   // double x_vel = electron_velocity[array_index];//   + ((electron_force[array_index]   + new_electron_force[array_index])   * dt  * constants::K_A / 2); 
    //double y_vel = electron_velocity[array_index+1];// + ((electron_force[array_index_y) + new_electron_force[array_index_y)) * dt  * constants::K_A / 2);
    //double z_vel = electron_velocity[array_index+2];// + ((electron_force[array_index_z) + new_electron_force[array_index_z)) * dt  * constants::K_A / 2);
 //   double vel = sqrt((x_vel*x_vel)+(y_vel*y_vel)+(z_vel*z_vel));
  //  double theta = atan(y_vel / x_vel);
   // double phi = acos(z_vel / vel);
    //if(x_vel < 0.0) theta += M_PI;
    const double theta = omp_uniform_random[thread]() * 2.0 * M_PI;
    const double phi   = omp_uniform_random[thread]() * M_PI; 

    electron_potential[e] += EKE;
    
    const double vel = sqrt(2.0*electron_potential[e]*constants::m_e_r_i);
      //  if(electron_potential[e] < 0.9817*E_f_A) std::cout << electron_potential[e] << ", " << EKE << ", " << vel << std::endl;
    electron_velocity[array_index]   = vel*cos(theta)*sin(phi);
    electron_velocity[array_index+1] = vel*sin(theta)*sin(phi);
    electron_velocity[array_index+2] = vel*cos(phi);

  //  electron_ee_scattering_list[e][0] = 1;
    electron_ea_scattering_list[e][0] = 1;
}
/*
// void e_a_coulomb(const int& e, const int& array_index]{

//     double x_distance;
//     double y_distance;
//     double z_distance;
   
//     double length;
 
//     int array_index_a, nearest_electron_count = 1;
//     int count = 2;

//   for (int a = 0; a < lattice_atoms; a++) {

//         array_index_a = 3*a;
//         x_distance = electron_position[array_index]     - atom_position[array_index_a);
//         y_distance = electron_position[array_index + 1) - atom_position[array_index_a + 1);
//         z_distance = electron_position[array_index + 2) - atom_position[array_index_a + 2); 
       
//         if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
//         else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

//         if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
//         else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

//         if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
//         else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

//         length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
        
//         if(length > e_a_neighbor_cutoff) continue;

//         atomic_nearest_electron_list[e][nearest_electron_count) = array_index_a;
//         nearest_electron_count++;

//         if (length > e_a_coulomb_cutoff) continue;
     
//         if(mean_radius[2*e) > length) mean_radius[2*e) = length;

//         if(ea_coupling) {
//             electron_ea_scattering_list[e][count) = a;
//             count++;
//         }
//   }

//  atomic_nearest_electron_list[e][0] = nearest_electron_count;
//  electron_ea_scattering_list[e][1] = count;
 
// }

// void neighbor_e_a_coulomb(const int& e, const int& array_index]{
                     
//     double x_distance;
//     double y_distance;
//     double z_distance;
    
//     double length;
//     int size = atomic_nearest_electron_list[e][0];
//     int count = 2;
  
   
//     for (int a = 1; a < size; a++) {
      
//         int array_index_a = atomic_nearest_electron_list[e][a);
        
//         x_distance = electron_position[array_index]     - atom_position[array_index_a);
//         y_distance = electron_position[array_index + 1) - atom_position[array_index_a + 1);
//         z_distance = electron_position[array_index + 2) - atom_position[array_index_a + 2); 
      
//         if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
//         else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

//         if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
//         else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

//         if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
//         else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

//       length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
   
//       if (length > e_a_coulomb_cutoff) continue;
     
//       if(mean_radius[2*e) > length) mean_radius[2*e) = length;
      
//       if(ea_coupling) {
//         electron_ea_scattering_list[e][count) = array_index_a / 3;
//         count++;
//       }
//     }
//    electron_ea_scattering_list[e][1] = count;
// }
*/
void e_e_coulomb(const int e, const int array_index) {

     int x_cell = int(floor(electron_position[array_index] / x_step_size));
      //if (x_cell < 0 || x_cell > x_omp_cells) std::cout << electron_position[array_index] << ", " << x_step_size << ", " << floor(electron_position[array_index] / x_step_size) << std::endl;

     int y_cell = int(floor(electron_position[array_index+1] / y_step_size));
      //if (y_cell < 0 || y_cell > y_omp_cells) std::cout << electron_position[array_index+1] << ", " << y_step_size << ", " << floor(electron_position[array_index+1] / y_step_size) << std::endl;
    
     int z_cell = int(floor(electron_position[array_index+2] / z_step_size));
      //if (z_cell < 0 || z_cell > z_omp_cells) std::cout << electron_position[array_index+2] << ", " << z_step_size << ", " << floor(electron_position[array_index+2] / z_step_size) << std::endl;
    //std::cout << x_cell << ", " << y_cell << ", " << z_cell << ", " << electron_position[array_index] << ", " << electron_position[array_index+1] << electron_position[array_index+2]  << std::endl;
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
  ///  std::cout << lattice_cell_coordinate[0][0][0] << std::endl;
    const unsigned int cell = lattice_cell_coordinate[x_cell][y_cell][z_cell];
    unsigned int ee_dos_count = 1;
    unsigned int ee_integration_count = 1;
    unsigned int ee_scattering_list = 2;
   
    for(int s = 0; s < 27; s++) {
      const int int_cell = cell_nearest_neighbor_list[cell][s];
      const int size = cell_integration_lists[int_cell][0];
    //  std::cout << size << std::endl;
      //if(e==0) std::cout << cell << ", " << int_cell << ", " << size << std::endl;
      for(int i = 1; i < size; i++) {
        const int electron = cell_integration_lists[int_cell][i];
        if(e == electron) continue;
        const int array_index_i = 3*electron;

        double x_distance = electron_position[array_index]   - electron_position[array_index_i];
        double y_distance = electron_position[array_index+1] - electron_position[array_index_i + 1];
        double z_distance = electron_position[array_index+2] - electron_position[array_index_i + 2]; 
    
        if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance += lattice_width;
        else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance -= lattice_width;

        if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance += lattice_depth;
        else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance -= lattice_depth;

        if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance += lattice_height;
        else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance -= lattice_height;
        
        const double length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
            
        if(length == 0.0) continue;
       
            if (length != length) std::cout << x_distance << ", " << y_distance << ", " << z_distance <<  ", " << \
        electron_position[array_index] << ", " << electron_position[array_index_i] << ", " << electron_position[array_index+1] << ", " << electron_position[array_index_i + 1] << ", " <<\
        electron_position[array_index+2] << ", " <<  electron_position[array_index_i + 2] << std::endl;
        
        if(length > e_e_integration_cutoff) continue;
        electron_integration_list[e][ee_integration_count] = array_index_i;
         if(ee_integration_count >= electron_integration_list[e].size() - 2) {std::cout << e << ", " << ee_integration_count << " > " << electron_integration_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
          break; }
        ee_integration_count++;

        if(length > e_e_neighbor_cutoff) continue;
        electron_nearest_electron_list[e][ee_dos_count] = electron;
        ee_dos_count++;
       //  std::cout << e << ", " << i << ", " << length << std::endl;
        
        if(electron_potential[electron] > transport_cutoff) ee_dos_hist[e].at(int(std::min(51.0, std::max((transport_cutoff-core_cutoff)/4.0, (transport_cutoff-core_cutoff)/4.0 + floor((electron_potential[electron]-transport_cutoff)/1)))))++;
        else ee_dos_hist[e].at(int(std::max(0.0, floor((electron_potential[electron]-core_cutoff)/4))))++;
        
        if(ee_dos_count >= electron_nearest_electron_list[e].size() - 2) {std::cout << e << ", " << ee_dos_count << " > " << electron_nearest_electron_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
          break; }

        if (length > e_e_coulomb_cutoff) continue; 
        electron_ee_scattering_list[e][ee_scattering_list] = electron;
      if(ee_scattering_list >= electron_ee_scattering_list[e].size() - 2) {std::cout << e << ", " << ee_scattering_list << " > " << electron_ee_scattering_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
          break; }
        ee_scattering_list++;
      //   std::cout << e << ", " << i << ", " << length << std::endl;
    }
    }
    electron_integration_list[e][0] = ee_integration_count;
    electron_nearest_electron_list[e][0] = ee_dos_count;
    electron_ee_scattering_list[e][1] = ee_scattering_list;
   // if(e % 100000 == 0 ) std::cout << ee_integration_count << ", " << ee_dos_count << ", " << ee_scattering_list << std::endl;
  }
void neighbor_e_e_coulomb(const int e, const int array_index) {
    
    const unsigned int size = electron_integration_list[e][0];
    unsigned int neighbor_count = 1;
    unsigned int scattering_count = 2;
 
    for (int i = 1; i < size; i++) {
        const unsigned int array_index_i = electron_integration_list[e][i];
        const int electron = array_index_i/3;

        double x_distance = electron_position[array_index]   - electron_position[array_index_i];
        double y_distance = electron_position[array_index+1] - electron_position[array_index_i + 1];
        double z_distance = electron_position[array_index+2] - electron_position[array_index_i + 2]; 
    
        if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance = x_distance + lattice_width;
        else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance = x_distance - lattice_width;

        if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance = y_distance + lattice_depth;
        else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance = y_distance - lattice_depth;

        if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance = z_distance + lattice_height;
        else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance = z_distance - lattice_height;
        
        const double length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
        if(length == 0.0) continue;
        
       if (length != length) std::cout << x_distance << ", " << y_distance << ", " << z_distance <<  ", " << \
        electron_position[array_index] << ", " << electron_position[array_index_i] << ", " << electron_position[array_index+1] << ", " << electron_position[array_index_i + 1] << ", " <<\
        electron_position[array_index+2] << ", " <<  electron_position[array_index_i + 2] << std::endl;
        
        if(length > e_e_neighbor_cutoff) continue;
        electron_nearest_electron_list[e][neighbor_count] = array_index_i/3;
        neighbor_count++;
      
        if(electron_potential[electron] > transport_cutoff) ee_dos_hist[e].at(int(std::min(51.0, std::max((transport_cutoff-core_cutoff)/4.0, (transport_cutoff-core_cutoff)/4.0 + floor((electron_potential[electron]-transport_cutoff)/1)))))++;
        else ee_dos_hist[e].at(int(std::max(0.0, floor((electron_potential[electron]-core_cutoff)/4))))++;
      
         if(neighbor_count >= electron_nearest_electron_list[e].size() - 2) {std::cout << e << ", " << neighbor_count << " > " << electron_nearest_electron_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
           break; }

        if (length > e_e_coulomb_cutoff) continue; 
        electron_ee_scattering_list[e][scattering_count] = array_index_i/3;
        scattering_count++;
        if(scattering_count >= electron_ee_scattering_list[e].size() - 1) {std::cout << e << ", " << scattering_count << " > " << electron_ee_scattering_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
           std::cout << x_distance << ", " << y_distance << ", " << z_distance << std::endl;
           break; }
      //   std::cout << e << ", " << i << ", " << length << std::endl;
    }

    electron_nearest_electron_list[e][0] = neighbor_count;
    electron_ee_scattering_list[e][1] = scattering_count;
}
/*
// void a_a_coulomb(const int a, const int array_index, \
//                         double& a_x_force, double& a_y_force, double& a_z_force, double& LPE) {
   
//     double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
//  //   int neighbor_count = 1;

//         x_distance = atom_position[array_index];
//         y_distance = atom_position[array_index+1];
//         z_distance = atom_position[array_index+2];
      
//                // if (err::check)  std::cout << "Calculating atomic interactions" << std::endl;

//        // neighbor_count = 1;

//         for (int i = 0; i < lattice_atoms; i++) {
//             if (i == a) continue; //no self repulsion

//             x_distance -= atom_position[array_index];
//             y_distance -= atom_position[array_index + 1);
//             z_distance -= atom_position[array_index + 2);

//             length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
           
//             if(length > 6) continue;

//             //length = sqrt(length);
//          //   force = -2000*(length - 2);
//            // PE = 1000*(length - 2)*(length - 2);
//         //   if(e == 0) std::cout << force << ", " << length <<  std::endl;
 
//             //phi   = acos(z_distance / length);
//             //theta = atan(y_distance / x_distance);
//             if (x_distance < 0) theta += M_PI;

//             atom_position[array_index]     += force * cos(theta)*sin(phi) / atomic_mass;
//             atom_position[array_index + 1) += force * sin(theta)*sin(phi) / atomic_mass;
//             atom_position[array_index + 2) += force * cos(phi) / atomic_mass;
//         }
//    // atomic_nearest_atom_list[a][0] = neighbor_count;
//     LPE += PE/2;
//    // if(a == 100) std::cout << neighbor_count << std::endl;
// }

// void neighbor_a_a_coulomb(const int a, const int array_index, \
//                         double& a_x_force, double& a_y_force, double& a_z_force, double& LPE) {
    
//     int array_index_i;
//     double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
//     double size = atomic_nearest_atom_list[a][0];

//         x_distance = atom_position[array_index];
//         y_distance = atom_position[array_index+1];
//         z_distance = atom_position[array_index+2];

//         for (int i = 1; i < size; i++) {
//            // if (i == a) continue; //no self repulsion
           
//             array_index_i = atomic_nearest_atom_list[a][i);

//             x_distance -= atom_position[array_index_i];
//             y_distance -= atom_position[array_index_i + 1);
//             z_distance -= atom_position[array_index_i + 2); 

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
  double e_energy = electron_potential[e];
  double deltaE = 0.5*constants::m_e_r*(((vel+electron_velocity[array_index])*(vel+electron_velocity[array_index]))+(electron_velocity[array_index+1]*electron_velocity[array_index+1])+(electron_velocity[array_index+2]*electron_velocity[array_index+2])) - e_energy;

  double DoS1 = electron_nearest_electron_list[e][0] - 1;
if(e_energy + deltaE > 0.9817*E_f_A) {
  for (int i = 1; i < electron_nearest_electron_list[e][0]; i++) {
    if((electron_potential[electron_nearest_electron_list[e][i]] < e_energy + fmax(1.2*deltaE, 0.8*deltaE)) \
    && (electron_potential[electron_nearest_electron_list[e][i]] > e_energy + fmin(0.8*deltaE, 1.2*deltaE))) {
      DoS1 -= 1.0;
    }
  }
  DoS1 /= double(electron_nearest_electron_list[e][0] - 1.0);
  if(omp_uniform_random[omp_get_thread_num()]() > exp(-3.0*DoS1)) {
    electron_velocity[array_index] += vel;
    external_potential += deltaE;
    electron_potential[e] += deltaE;
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
  
//       if(atomic_nearest_atom_list[e][0]) continue;
      
//       double deltaE = atom_potential[e] - E_f_A; - 3.0*constants::kB_r*Tp;

//      // if(deltaE < 8.0) scattering_prob = 8.0;
//        scattering_prob = deltaE*deltaE;

//       if(omp_uniform_random[omp_get_thread_num())() > exp(aa_rate*scattering_prob)) {
        
//         deltaE *= 0.5-0.05*exp(-2.0*deltaE/(E_f_A-3.0*constants::kB_r*Tp));
//         int size = atomic_nearest_atom_list[e][1] - 2;
//         size = 2+ (omp_int_random[omp_get_thread_num())() % size);
//       //  if(size >= atomic_nearest_atom_list[e].size() -1) std::cout << "bus error: " << size << " < " << atomic_nearest_atom_list[e].size() << std::endl; 
//         int a = atomic_nearest_atom_list[e][size);
      
//        // if(a >= conduction_electrons -1) std::cout << "bus error: " << a << " < " << conduction_electrons-1 << std::endl; 
//         if(atomic_nearest_atom_list[a][0]) continue;

//         deltaE *= 0.5;
//         atom_potential[e] -= deltaE;
//         atom_potential[a) += deltaE;  

//         atomic_nearest_atom_list[e][0] = 1;
//         atomic_nearest_atom_list[a][0] = 1;
        
//         a_a_scattering_count++;
//       }
//     }
//     if (err::check) std::cout << "aa_scattering done." << std::endl;
// }
*/
void ea_scattering(const int e, const int array_index, const int thread) {

   // if(!electron_transport_list[e]) return;
    if(Te == Tp) return;
    double e_energy = electron_potential[e];
   
    //if(!electron_transport_list[e]) return;
    if(omp_uniform_random[thread]() > exp(ea_rate/e_energy)) {
       double deltaE = sqrt(phonon_energy*E_f_A);
       if( e_energy + deltaE > (E_f_A+37.0) ) return;
      const unsigned int size = electron_nearest_electron_list[e][0];
      const double DoS_width = 4;//AJ
      const double FD_width = std::min(size - 1.0, ee_density/3.0)/(3*constants::kB_r*300.0+E_f_A - core_cutoff);
      double e_dos; 
      double d_e_dos;
        
        // int hist;
        // if(e_energy < transport_cutoff) {hist = int(std::max(0.0,  floor((e_energy - core_cutoff)/4.0))); e_dos = DoS_width*FD_width;}
        // else {hist = int(std::max(ee_dos_hist[0].size() -1.0, ((transport_cutoff-core_cutoff)/4.0) + floor((e_energy - transport_cutoff)/1.0)));  e_dos = DoS_width*FD_width/4.0;}
        // double e_occupation = double(ee_dos_hist[e].at(hist))/e_dos;

        // if(e_energy + deltaE < transport_cutoff) {hist = int(std::min(0.0,  floor((e_energy + deltaE - core_cutoff)/4.0)));  d_e_dos = DoS_width*FD_width;} 
        // else {hist = int(std::max(ee_dos_hist[0].size() -1.0, ((transport_cutoff-core_cutoff)/4.0) + floor((e_energy + deltaE - transport_cutoff)/1.0)));  d_e_dos = DoS_width*FD_width/4.0;}
        // double d_e_occupation = double(ee_dos_hist[e].at(hist))/d_e_dos;

        // if(e_energy < transport_cutoff) e_occupation = 1.0;
        // if(e_energy + deltaE < transport_cutoff) d_e_occupation = 1.0;
        // if(e_occupation - d_e_occupation == 0.0) return;
      // if(scattering_velocity + deltaE < transport_cutoff) {
      //   DoS1 = int(round(DoS_width*ee_density/3.0/(3*constants::kB_r*300.0+E_f_A - core_cutoff)));
      //   DoS1_n = DoS1;
      //   for (int i = 1; i < size ; i++) {
      //     if( ((scattering_velocity + deltaE -  0.5*DoS_width) < core_cutoff) || ((electron_potential[electron_nearest_electron_list[e][i)) < scattering_velocity + deltaE + 0.5*DoS_width) \
      //     && (electron_potential[electron_nearest_electron_list[e][i)) > scattering_velocity + deltaE - 0.5*DoS_width) ) ) {
      //       DoS1 -= 1.0;
      //       if(DoS1 < 0) {DoS1 = 0; break;}
      //     }
      //   }
      // } else {
      //   DoS1 = int(round(DoS_width*ee_density/3.0/(3*constants::kB_r*300.0+E_f_A - core_cutoff))/6.0);
      //   DoS1_n = DoS1;
      //   for (int i = 1; i < size ; i++) {
      //     if( ((scattering_velocity + deltaE -  0.25) < core_cutoff) || ((electron_potential[electron_nearest_electron_list[e][i)) < scattering_velocity + deltaE + 0.25) \
      //     && (electron_potential[electron_nearest_electron_list[e][i)) > scattering_velocity + deltaE - 0.25) ) ) {
      //       DoS1 -= 1.0;
      //       if(DoS1 < 0) {DoS1 = 0; break;}
      //     }
      //   }
      // }
      // if(scattering_velocity < transport_cutoff) {
      //   DoS2 = int(round(DoS_width*ee_density/3.0/(3*constants::kB_r*300.0+E_f_A - core_cutoff)));
      //   DoS2_n = DoS2;
      //   for (int i = 1; i < size ; i++) {
      //     if ( (scattering_velocity - 0.5*DoS_width) < core_cutoff || ((electron_potential[electron_nearest_electron_list[e][i)) < scattering_velocity + 0.5*DoS_width) \
      //     && (electron_potential[electron_nearest_electron_list[e][i)) > scattering_velocity - 0.5*DoS_width) ) ) {
      //       DoS2 -= 1.0;
      //       if(DoS2 < 0) {DoS2 = 0; break;}
      //     }
      //   }
      // } else {
      //   DoS2 = int(round(DoS_width*ee_density/3.0/(3*constants::kB_r*300.0+E_f_A - core_cutoff))/6.0);
      //   DoS2_n = DoS2;
      //   for (int i = 1; i < size ; i++) {
      //     if ( (scattering_velocity - 0.25) < core_cutoff || ((electron_potential[electron_nearest_electron_list[e][i)) < scattering_velocity + 0.25) \
      //     && (electron_potential[electron_nearest_electron_list[e][i)) > scattering_velocity - 0.25) ) ) {
      //       DoS2 -= 1.0;
      //       if(DoS2 < 0) {DoS2 = 0; break;}
      //     }
      //   }
      // }
     // if(e_occupation > 180.0/195.0 && d_e_occupation > 180.0/195.0) return;
      if(Tp <= Te) deltaE *= -1.0;
      if(return_phonon_distribution(e_energy/E_f_A, constants::kB_r*Te/E_f_A) < return_phonon_distribution((e_energy+deltaE)/E_f_A, constants::kB_r*Te/E_f_A)) deltaE *= -1.0;
      if((return_phonon_distribution(e_energy/E_f_A, constants::kB_r*Te/E_f_A) - return_phonon_distribution((e_energy+deltaE)/E_f_A, constants::kB_r*Te/E_f_A)) , 1e-3) return;
      if(e_energy + deltaE < core_cutoff ) return;

      electron_potential[e] += deltaE;

      electron_transport_list[e] = true;
      if(electron_potential[e] < transport_cutoff) electron_transport_list[e] = false;

      const double theta = 2.0*M_PI*omp_uniform_random[thread]();
      const double phi   = M_PI*omp_uniform_random[thread]();
     // if(electron_velocity[array_index+2] < 0.0) theta += M_PI;

      const double scattering_velocity = sqrt(2.0*electron_potential[e]*constants::m_e_r_i);
      electron_velocity[array_index]   = scattering_velocity * cos(theta)*sin(phi);
      electron_velocity[array_index+1] = scattering_velocity * sin(theta)*sin(phi);
      electron_velocity[array_index+2] = scattering_velocity * cos(phi); 
    
     // electron_ee_scattering_list[e][0] = 1;

      #pragma omp critical(eascattering)
      {
      if (e_energy < transport_cutoff) ea_core_scattering_count++;
      else ea_transport_scattering_count++;
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
  const static double q_sq = 0.25;
  const static double DoS_width = 4;//AJ

  #pragma omp parallel reduction(+:e_e_scattering_count, ee_core_scattering_count, ee_transport_scattering_count) 
  {
    
    const int thread = omp_get_thread_num();

  for(int l = 0; l < cells_per_thread; l++) {
    const unsigned int cell = lattice_cells_per_omp[thread][l];
    const unsigned int size = cell_integration_lists[cell][0];
    
    #pragma omp barrier 

    for(int e = 1; e < size; e++) {
      const unsigned int electron = cell_integration_lists[cell][e];
        if (electron_ee_scattering_list[electron][0] == 1) continue;
      const int scattering_size = electron_ee_scattering_list[electron][1];

      for(int a = 2; a < scattering_size; a++) {
        
        int electron_collision = electron_ee_scattering_list[electron][a];
        if(electron_ee_scattering_list[electron_collision][0] == 1) continue;

        const double e_energy = electron_potential[electron];
        const double d_e_energy = electron_potential[electron_collision];
       if(equilibrium_step) {  
          int blocking_lr;
          int blocking_hr;
          double DoS_lhs;
          double DoS_rhs;
           
          double deltaE = omp_uniform_random[thread]()*(e_energy - d_e_energy);
          //if(test_uniform(gen) < 0.5*exp(abs(deltaE)/(-0.5))) deltaE *= -1.0;
          if(e_energy - deltaE < core_cutoff || d_e_energy + deltaE < core_cutoff) continue;

          if( (e_energy - deltaE) > (E_f_A+37.0) || (d_e_energy + deltaE > (E_f_A+37.0))) continue;
         // const double FD_width = (electron_nearest_electron_list[electron][0]-1.0)/(3.0*constants::kB_r*300.0+E_f_A - core_cutoff);
          //if(electron_nearest_electron_list[electron][0]-1.0 < ee_density/3.0) std::cout << electron_nearest_electron_list[electron][0]-1.0 << ", " << ee_density/3.0 << ", " << FD_width << std::endl;
          double e_occupation;
          // double e_dos;
          double d_e_occupation;
          // double d_e_dos;
          // int hist;
          // if(e_energy - deltaE < transport_cutoff) {
          //   hist = int(std::max(0.0, floor((e_energy - deltaE - core_cutoff)/4.0))); 
          //   e_dos = DoS_width*FD_width;
          // } else {
          //   hist = int(std::min(ee_dos_hist[0].size() - 1.0, ((transport_cutoff-core_cutoff)/4.0) + floor((e_energy - deltaE - transport_cutoff)/1.0))); 
          //   e_dos = DoS_width*FD_width/4.0;
          // }
          // e_occupation =  1.0 - ee_dos_hist[electron].at(hist)/e_dos;  
          // if(d_e_energy + deltaE < transport_cutoff) {
          //   hist = int(std::max(0.0, floor((d_e_energy + deltaE - core_cutoff)/4.0))); 
          //   d_e_dos = DoS_width*FD_width; 
          // } else {
          //   hist = int(std::min(ee_dos_hist[0].size() - 1.0, ((transport_cutoff-core_cutoff)/4.0) + floor((d_e_energy + deltaE - transport_cutoff)/1.0))); 
          //   d_e_dos = DoS_width*FD_width/4.0;
          // }
          e_occupation = 1.0 - return_phonon_distribution((e_energy-deltaE-E_f_A)/E_f_A, constants::kB_r*Te/E_f_A);
          d_e_occupation = 1.0 - return_phonon_distribution((d_e_energy+deltaE-E_f_A)/E_f_A, constants::kB_r*Te/E_f_A);
          
          // d_e_occupation = 1.0 - ee_dos_hist[electron_collision].at(hist)/d_e_dos;
          // //if( (e_occupation < 0.0 || d_e_occupation < 0.0) && (d_e_energy + deltaE > transport_cutoff || e_energy - deltaE > transport_cutoff)) {
          //   // #pragma omp critical
          //  // std::cout << e_occupation << ", " << e_dos << ", " << d_e_occupation << ", " << d_e_dos << ", " << FD_width << std::endl;
          // //}
          // //if(hist == 11 || hist == 10) std::cout << "ee d_occupation problem: " << d_e_occupation << ", " << d_e_dos << ", " << hist << ", " << d_e_energy << std::endl; 
          // if(e_occupation < (15.0/195.0) || d_e_occupation < (15.0/195.0)) continue; 

          const int array_index = 3*electron;
          const int array_index_i = 3*electron_collision;
          double theta = atan2(electron_velocity[array_index+1], electron_velocity[array_index]);
              if(theta != theta) theta = 0.0;
          double phi = acos(electron_velocity[array_index+2] / sqrt(2.0*constants::m_e_r_i*e_energy));
          const double d_theta = 4.0*M_PI*ee_scattering_angle*(omp_uniform_random[thread]() - 0.5);
          const double d_phi = 2.0*M_PI*ee_scattering_angle*(omp_uniform_random[thread]() - 0.5);
          const double v_1 = sqrt(2.0*constants::m_e_r_i*(e_energy-deltaE));

          const double k_1_x = electron_velocity[array_index]  *constants::m_over_hbar_sqrt;
          const double k_1_y = electron_velocity[array_index+1]*constants::m_over_hbar_sqrt;
          const double k_1_z = electron_velocity[array_index+2]*constants::m_over_hbar_sqrt;
          const double k_2_x = v_1*cos(theta+d_theta)*sin(phi+d_phi)*constants::m_over_hbar_sqrt;
          const double k_2_y = v_1*sin(theta+d_theta)*sin(phi+d_phi)*constants::m_over_hbar_sqrt;
          const double k_2_z = v_1*cos(phi+d_phi)*constants::m_over_hbar_sqrt;

          // const double deltaK = (electron_velocity[array_index]  -electron_velocity[array_index_i])  *(electron_velocity[array_index]  -electron_velocity[array_index_i])\
                              + (electron_velocity[array_index+1]-electron_velocity[array_index_i+1])*(electron_velocity[array_index+1]-electron_velocity[array_index_i+1])\
                              + (electron_velocity[array_index+2]-electron_velocity[array_index_i+2])*(electron_velocity[array_index+2]-electron_velocity[array_index_i+2]);
          const double deltaK = (k_1_x-k_2_x)*(k_1_x-k_2_x) + (k_1_y-k_2_y)*(k_1_y-k_2_y) + (k_1_z-k_2_z)*(k_1_z-k_2_z);
        if(omp_uniform_random[thread]() > exp(-1.0*ee_rate*e_occupation*d_e_occupation/((q_sq+(deltaK))*(q_sq+(deltaK))))) {
        //      #pragma omp critical
        //  if(e_energy-deltaE < transport_cutoff) std::cout << e_energy << ", " << e_energy-deltaE << ", " << e_occupation << ", " << 1.0 - return_phonon_distribution((e_energy-deltaE-E_f_A)/E_f_A, constants::kB_r*Te/E_f_A) << ", " << deltaK << ", " <<  exp(-1.0*ee_rate*e_occupation*d_e_occupation/((q_sq+(deltaK))*(q_sq+(deltaK)))) << std::endl;
        //   #pragma omp critical
        //   if(d_e_energy+deltaE < transport_cutoff) std::cout << d_e_energy <<  ", " << deltaE+d_e_energy << ", " << d_e_occupation << ", " << 1.0 - return_phonon_distribution((d_e_energy+deltaE-E_f_A)/E_f_A, constants::kB_r*Te/E_f_A) << ", " << deltaK << ", " <<  exp(-1.0*ee_rate*e_occupation*d_e_occupation/((q_sq+(deltaK))*(q_sq+(deltaK)))) << std::endl;
          
            if (e_energy < transport_cutoff) ee_core_scattering_count++;
            else ee_transport_scattering_count++;
            if (d_e_energy < transport_cutoff) ee_core_scattering_count++;
            else ee_transport_scattering_count++;

            electron_potential[electron] -= deltaE;
            electron_potential[electron_collision]   += deltaE;

            electron_transport_list[electron] = true;
            electron_transport_list[electron_collision] = true;
            if (electron_potential[electron] < transport_cutoff)  electron_transport_list[electron] = false;
            if (electron_potential[electron_collision] < transport_cutoff) electron_transport_list[electron_collision] = false;
            
            const double v_2 = sqrt(2.0*constants::m_e_r_i*electron_potential[electron_collision]);

            electron_velocity[array_index]   = v_1*cos(theta+d_theta)*sin(phi+d_phi);// k_2_x*constants::hbar_over_me_sqrt;
            electron_velocity[array_index+1] = v_1*sin(theta+d_theta)*sin(phi+d_phi);// k_2_y*constants::hbar_over_me_sqrt;
            electron_velocity[array_index+2] = v_1*cos(phi+d_phi);// k_2_z*constants::hbar_over_me_sqrt;

            theta = atan2(electron_velocity[array_index_i+1], electron_velocity[array_index_i]);
            phi = acos(electron_velocity[array_index_i+2] / sqrt(2.0*constants::m_e_r_i*d_e_energy));

            electron_velocity[array_index_i]   = v_2*cos(theta-d_theta)*sin(phi-d_phi);
            electron_velocity[array_index_i+1] = v_2*sin(theta-d_theta)*sin(phi-d_phi);
            electron_velocity[array_index_i+2] = v_2*cos(phi-d_phi);

            e_e_scattering_count += 2;

            electron_ee_scattering_list[electron][0] = 1;
            electron_ee_scattering_list[electron_collision][0] = 1;
            break;
          }
        } else if (int(ee_elastic(electron, electron_collision, e_energy, d_e_energy, omp_uniform_random[thread]() )) == 2) {
          e_e_scattering_count += 2;
          break;
       }
      }
      }
    } 
  }   if (err::check) std::cout << "ee_scattering done." << std::endl;
  
}

int ee_energy_conserved(const int electron, const int electron_collision, const double deltaE){

        std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> test_int(0, conduction_electrons -1);
    std::uniform_real_distribution<double> test_uniform;

      const unsigned int array_index = 3*electron;
      const unsigned int array_index_i = 3*electron_collision;
       double theta = 2.0*M_PI*test_uniform(gen);
       double phi = M_PI*test_uniform(gen);

      const double v_1 = sqrt(2.0*constants::m_e_r_i*electron_potential[electron]);
      const double v_2 = sqrt(2.0*constants::m_e_r_i*electron_potential[electron_collision]);

      electron_velocity[array_index]   = v_1*cos(theta)*sin(phi);
      electron_velocity[array_index+1] = v_1*sin(theta)*sin(phi);
      electron_velocity[array_index+2] = v_1*cos(phi);

       theta = 2.0*M_PI*test_uniform(gen);
        phi = M_PI*test_uniform(gen);
      electron_velocity[array_index_i]   = v_2*cos(theta)*sin(phi);
      electron_velocity[array_index_i+1] = v_2*sin(theta)*sin(phi);
      electron_velocity[array_index_i+2] = v_2*cos(phi);

    return 2;
}

int ee_final_momentum_conserved(const int electron, const int electron_collision, const double deltaE, const double e_energy, const double d_e_energy){

        std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> test_int(0, conduction_electrons -1);
    std::uniform_real_distribution<double> test_uniform;

      const unsigned int array_index = 3*electron;
      const unsigned int array_index_i = 3*electron_collision;
      double theta = atan2(electron_velocity[array_index+1], electron_velocity[array_index]);
      double phi = acos(electron_velocity[array_index+2] / sqrt(2.0*constants::m_e_r_i*e_energy));
      const double d_theta = 4.0*M_PI*ee_scattering_angle*(test_uniform(gen) - 0.5);
      const double d_phi = 2.0*M_PI*ee_scattering_angle*(test_uniform(gen) - 0.5);

      const double v_1 = sqrt(2.0*constants::m_e_r_i*electron_potential[electron]);
      const double v_2 = sqrt(2.0*constants::m_e_r_i*electron_potential[electron_collision]);

      electron_velocity[array_index]   = v_1*cos(theta+d_theta)*sin(phi+d_phi);
      electron_velocity[array_index+1] = v_1*sin(theta+d_theta)*sin(phi+d_phi);
      electron_velocity[array_index+2] = v_1*cos(phi+d_phi);

       theta = atan2(electron_velocity[array_index_i+1], electron_velocity[array_index_i]);
       phi = acos(electron_velocity[array_index_i+2] / sqrt(2.0*constants::m_e_r_i*d_e_energy));
      electron_velocity[array_index_i]   = v_2*cos(theta-d_theta)*sin(phi-d_phi);
      electron_velocity[array_index_i+1] = v_2*sin(theta-d_theta)*sin(phi-d_phi);
      electron_velocity[array_index_i+2] = v_2*cos(phi-d_phi);

    return 2;
}

int ee_elastic(const int electron, const int electron_collision, const double e_energy, const double d_e_energy, const double probability) {
   const unsigned int array_index = 3*electron;
   const unsigned int array_index_i = 3*electron_collision;

    double x_distance = electron_position[array_index]  - electron_position[array_index_i];
    double y_distance = electron_position[array_index+1]- electron_position[array_index_i+1];
    double z_distance = electron_position[array_index+2]- electron_position[array_index_i+2]; 

    if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance += lattice_width;
    else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance -= lattice_width;

    if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance += lattice_depth;
    else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance -= lattice_depth;
    
    if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance += lattice_height;
    else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance -= lattice_height;
              
    const double length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
      if(length == 0.0) return 0;
    const double v_x = electron_velocity[array_index];
    const double v_y = electron_velocity[array_index+1];
    const double v_z = electron_velocity[array_index+2];
    const double v_x_2 = electron_velocity[array_index_i];
    const double v_y_2 = electron_velocity[array_index_i+1];
    const double v_z_2 = electron_velocity[array_index_i+2];

    const double v_x_dot_product = (( v_x - v_x_2)*x_distance\
                                  + (v_y - v_y_2)*y_distance\
                                  + (v_z - v_z_2)*z_distance) /length;

      if(v_x_dot_product == 0.0) return 0;
    const double normalised_dot_product = v_x_dot_product;

    const double d_v_x_1 = v_x_2 + x_distance*normalised_dot_product;
    const double d_v_y_1 = v_y_2 + y_distance*normalised_dot_product;
    const double d_v_z_1 = v_z_2 + z_distance*normalised_dot_product;
    const double d_v_x_2 = v_x   - x_distance*normalised_dot_product;
    const double d_v_y_2 = v_y   - y_distance*normalised_dot_product;
    const double d_v_z_2 = v_z   - z_distance*normalised_dot_product;

    const double deltaE = -1.0*(0.5*constants::m_e_r*(d_v_x_1*d_v_x_1 + d_v_y_1*d_v_y_1 + d_v_z_1*d_v_z_1) -\
                           0.5*constants::m_e_r*(d_v_x_2*d_v_x_2 + d_v_y_2*d_v_y_2 + d_v_z_2*d_v_z_2));
    // #pragma omp critical
    if(deltaE != deltaE) std::cout << "deltaE " << deltaE << ", " << normalised_dot_product << ", " << length << std::endl;

      //if(e_energy - deltaE < transport_cutoff || d_e_energy + deltaE < transport_cutoff) return 0;
      if((e_energy - deltaE < core_cutoff) || (d_e_energy + deltaE < core_cutoff)) return 0;
      else {
        int blocking_lr;
        int blocking_hr;
        double DoS_lhs;
        double DoS_rhs;
        int DoS1_n;
        int DoS2_n;
        int DoS1; 
        int DoS2;
        // const double DoS_width = 4;
        if( (e_energy - deltaE) > (E_f_A+37.0) || (d_e_energy + deltaE > (E_f_A+37.0))) return 0;
        //   const double FD_width = (electron_nearest_electron_list[electron][0]-1.0)/(3.0*constants::kB_r*300.0+E_f_A - core_cutoff);
        //   //if(electron_nearest_electron_list[electron][0]-1.0 < ee_density/3.0) std::cout << electron_nearest_electron_list[electron][0]-1.0 << ", " << ee_density/3.0 << ", " << FD_width << std::endl;
          double e_occupation;
        //   double e_dos;
          double d_e_occupation;
        //   double d_e_dos;
        //   int hist;
        //   if(e_energy - deltaE < transport_cutoff) {
        //     hist = int(std::max(0.0, floor((e_energy - deltaE - core_cutoff)/4.0))); 
        //     e_dos = DoS_width*FD_width;
        //   } else {
        //     hist = int(std::min(ee_dos_hist[0].size() - 1.0, ((transport_cutoff-core_cutoff)/4.0) + floor((e_energy - deltaE - transport_cutoff)/1.0))); 
        //     e_dos = DoS_width*FD_width/4.0;
        //   }

        //   e_occupation =  1.0 - ee_dos_hist[electron].at(hist)/e_dos;  

        //   if(d_e_energy + deltaE < transport_cutoff) {
        //     hist = int(std::max(0.0, floor((d_e_energy + deltaE - core_cutoff)/4.0))); 
        //     d_e_dos = DoS_width*FD_width; 
        //   } else {
        //     hist = int(std::min(ee_dos_hist[0].size() - 1.0, ((transport_cutoff-core_cutoff)/4.0) + floor((d_e_energy + deltaE - transport_cutoff)/1.0))); 
        //     d_e_dos = DoS_width*FD_width/4.0;
        //   }

        //   d_e_occupation = 1.0 - ee_dos_hist[electron_collision].at(hist)/d_e_dos;

        // if(e_energy - deltaE < transport_cutoff) {
        //       DoS_rhs = 0.5*DoS_width + std::min(0.0, transport_cutoff - (e_energy - deltaE) - 0.5*DoS_width)*0.75; //AJ
        //       DoS_lhs = 0.5*DoS_width; //AJ
        //       DoS1 = int(round((DoS_lhs+DoS_rhs)*ee_density/3.0/(3*constants::kB_r*300.0+E_f_A - core_cutoff))); //AJ
        //       DoS1_n = DoS1;
        //       blocking_lr = int((15.0/195.0)*DoS1);
        //       for (int e = 1; e < electron_nearest_electron_list[electron][0]; ++e) {
        //        if(  ((electron_potential[electron_nearest_electron_list[electron][e]] < (e_energy - deltaE + DoS_rhs)) \
        //        && ( (electron_potential[electron_nearest_electron_list[electron][e]] > (e_energy - deltaE -  DoS_lhs) ) ) ) ){
        //           DoS1 -= 1;
        //           if(DoS1 < blocking_lr) {DoS1 = 0; break;}
        //         }
        //       }
        //    } else {
        //       DoS_lhs = 0.25*0.5*DoS_width - std::min(0.0, e_energy - transport_cutoff - deltaE - 0.25*0.5*DoS_width)*DoS_width*0.75;
        //       DoS_rhs = 0.25*0.5*DoS_width;
        //       DoS1 = int(round((DoS_lhs+DoS_rhs)*ee_density/3.0/(3*constants::kB_r*300.0+E_f_A - core_cutoff)));
        //       DoS1_n = DoS1;
        //       blocking_hr = int((15.0/195.0)*DoS1);
        //       for (int e = 1; e < electron_nearest_electron_list[electron][0]; ++e) {
        //         if( electron_potential[electron_nearest_electron_list[electron][e]] < (e_energy - deltaE + DoS_rhs)  \
        //         &&  (electron_potential[electron_nearest_electron_list[electron][e]] > (e_energy - deltaE - DoS_lhs) )) {
        //           DoS1 -= 1;
        //           if(DoS1 < blocking_hr) { DoS1 = 0; break;}
        //         }
        //       }
        //      //  std::cout << DoS1_n << ", " << DoS1 << ", " << electron << ", " << deltaE << std::endl;
        //     } 
        //     if(d_e_energy + deltaE < transport_cutoff) {
        //       DoS_rhs = 0.5*DoS_width + std::min(0.0, transport_cutoff - (d_e_energy + deltaE) - 0.5*DoS_width)*0.75;
        //       DoS_lhs = 0.5*DoS_width;
        //       DoS2 = int(round((DoS_lhs+DoS_rhs)*ee_density/3.0/(3*constants::kB_r*300.0+E_f_A - core_cutoff)));
        //       DoS2_n = DoS2;
        //       blocking_lr = int((15.0/195.0)*DoS2);
        //       for (int i = 1; i < electron_nearest_electron_list[electron_collision][0]; ++i) {
        //         if( ((electron_potential[electron_nearest_electron_list[electron_collision][i]] < (d_e_energy + deltaE + DoS_rhs)) \
        //         &&  ( (electron_potential[electron_nearest_electron_list[electron_collision][i]] > (d_e_energy + deltaE -  DoS_lhs) ) ) ) ) {
        //           DoS2 -= 1;
        //           if(DoS2 < blocking_lr) {DoS2 = 0; break;}
        //         }
        //       }
        //     } else {
        //       DoS_lhs = 0.25*0.5*DoS_width - std::min(0.0, d_e_energy - transport_cutoff + deltaE - 0.25*0.5*DoS_width)*DoS_width*0.75;
        //       DoS_rhs = 0.25*0.5*DoS_width;
        //       DoS2 = int(round((DoS_lhs+DoS_rhs)*ee_density/3.0/(3*constants::kB_r*300.0+E_f_A - core_cutoff)));
        //       DoS2_n = DoS2;
        //       blocking_hr = int((15.0/195.0)*DoS2);
        //       for (int i = 1; i < electron_nearest_electron_list[electron_collision][0]; ++i) {
        //         if((electron_potential[electron_nearest_electron_list[electron_collision][i]] < (d_e_energy + deltaE + DoS_rhs)) \
        //         && ( (electron_potential[electron_nearest_electron_list[electron_collision][i]] > (d_e_energy + deltaE - DoS_lhs) ) )) {
        //           DoS2 -= 1;
        //        ///   if(DoS2 < -10) std::cout << "need tighter bounds hr ee2 " << DoS2 << ", " << electron_nearest_electron_list[electron_collision][0] << ", " << d_e_energy << ", " << d_e_energy + deltaE - 0.004*E_f_A << ", " << d_e_energy + deltaE + 0.004*E_f_A << std::endl;
        //           if(DoS2 < blocking_hr) { DoS2 = 0; break;}
        //         }
        //       }
        //    //   std::cout << DoS2_n << ", " << DoS2 << ", " << electron_collision << ", " << deltaE << std::endl;
        //     } 
          e_occupation = 1.0 - return_phonon_distribution((e_energy-deltaE-E_f_A)/E_f_A, constants::kB_r*Te/E_f_A);
          d_e_occupation = 1.0 - return_phonon_distribution((d_e_energy+deltaE-E_f_A)/E_f_A, constants::kB_r*Te/E_f_A);
        //   if(e_occupation < (15.0/195.0) || d_e_occupation < (15.0/195.0)) return 0; 
          const double k_1_x = v_x*constants::m_over_hbar_sqrt;
          const double k_1_y = v_y*constants::m_over_hbar_sqrt;
          const double k_1_z = v_z*constants::m_over_hbar_sqrt;
          const double k_2_x = v_x - x_distance*normalised_dot_product;
          const double k_2_y = v_y - y_distance*normalised_dot_product;
          const double k_2_z = v_z - z_distance*normalised_dot_product;  

         const double deltaK = (k_1_x-k_2_x)*(k_1_x-k_2_x) + (k_1_y-k_2_y)*(k_1_y-k_2_y) + (k_1_z-k_2_z)*(k_1_z-k_2_z);

      if(probability > exp(-1.0*ee_rate*e_occupation*d_e_occupation/((0.25+(deltaK))*(0.25+(deltaK))))) {
       
        electron_potential[electron] -= deltaE;
        electron_potential[electron_collision]   += deltaE;

        electron_ee_scattering_list[electron][0] == 1;
        electron_ee_scattering_list[electron_collision][0] == 1;

        electron_transport_list[electron] = true;
        electron_transport_list[electron_collision] = true;
        if (electron_potential[electron] < transport_cutoff)  electron_transport_list[electron] = false;
        if (electron_potential[electron_collision] < transport_cutoff) electron_transport_list[electron_collision] = false;
    
        electron_velocity[array_index_i]  += x_distance*normalised_dot_product;
        electron_velocity[array_index_i+1]+= y_distance*normalised_dot_product; 
        electron_velocity[array_index_i+2]+= z_distance*normalised_dot_product;
    
        electron_velocity[array_index]  -= x_distance*normalised_dot_product;
        electron_velocity[array_index+1]-= y_distance*normalised_dot_product;
        electron_velocity[array_index+2]-= z_distance*normalised_dot_product; 

      #pragma omp critical 
      {
      if (electron_potential[electron] < transport_cutoff) ee_core_scattering_count++;
      else ee_transport_scattering_count++;
      if (electron_potential[electron_collision] < transport_cutoff) ee_core_scattering_count++;
      else ee_transport_scattering_count++;
      }
        return 2;

      }
     return 0;
    }
}


} //end CASTLE namespace

