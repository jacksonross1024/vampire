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
            if (err::check) std::cout << "Updating new electron position." << std::endl;
    update_position();

            if (err::check) std::cout << "Forces, spins, and velocities update..." << std::endl;
    update_dynamics();
          if (err::check) std::cout << "Output mean data" << std::endl;
    
    if (current_time_step % CASTLE_output_rate == 0)   output_data(); //std::cout << "x_flux: " << x_flux / CASTLE_output_rate << "\n"; x_flux = 0;
    
    //reset integration
    CASTLE_real_time += dt;
    current_time_step++;

   // electron_position.swap(electron_position);
   
    return EXIT_SUCCESS;
}

void setup_output() {

    CASTLE_output_data = true;
}


void update_position(){

    int array_index,array_index_y,array_index_z;
    double x_pos,y_pos,z_pos;

    #pragma omp parallel for private(array_index,array_index_y,array_index_z, x_pos,y_pos,z_pos) schedule(guided) reduction(+:x_flux,y_flux,z_flux)
    for (int e = 0; e < conduction_electrons; e++){ 
        if(!electron_transport_list[e]) continue;

        array_index = 3*e;
        array_index_y = array_index + 1;
        array_index_z = array_index + 2;
    
        x_pos = electron_position[array_index]   + (electron_velocity[array_index]   * dt);// + (electron_force[array_index]   * dt * dt * constants::K_A / 2); // x superarray component
        y_pos = electron_position[array_index_y] + (electron_velocity[array_index_y] * dt);// + (electron_force[array_index_y] * dt * dt * constants::K_A / 2); // y superarray component
        z_pos = electron_position[array_index_z] + (electron_velocity[array_index_z] * dt);// + (electron_force[array_index_z] * dt * dt * constants::K_A / 2); // z superarray component
        
     //   std::cout << electron_position[array_index] << ", " << electron_position[array_index_y] << ", " << electron_position[array_index_z] << ", " << electron_velocity[array_index] << ", " << electron_velocity[array_index_y] << ", " << electron_velocity[array_index_z] << std::endl;
        
        if (x_pos < 0.0) {x_pos += lattice_width; x_flux--;}
        else if (x_pos > lattice_width) {x_pos -= lattice_width; x_flux++;}

        if (y_pos < 0.0) {y_pos += lattice_depth; y_flux--;}
        else if (y_pos > lattice_depth) {y_pos -= lattice_depth; y_flux++;}

        if (z_pos < 0.0) {z_pos += lattice_height; z_flux--;}
        else if (z_pos > lattice_height) {z_pos -= lattice_height; z_flux++;}

        electron_position[array_index]   = x_pos;
        electron_position[array_index_y] = y_pos;
        electron_position[array_index_z] = z_pos;

        // if(!equilibrium_step) {
        //   // if(heat_pulse_sim) {
        //   //   if(x_pos < 22.0 && x_pos > 14.0 && y_pos > 14.0 && y_pos < 22.0 && z_pos > 14.0 && z_pos < 22.0 ) {
        //   //     external_interaction_list[e] = true;
        //   //     external_interaction_list_count++;
        //   //   }
        //     // x_pos = atom_position[array_index];
        //     // y_pos = atom_position[array_index_y];
        //     // z_pos = atom_position[array_index_z];
        //    // if(x_pos < 22.0 && x_pos > 14.0 && y_pos > 14.0 && y_pos < 22.0 && z_pos > 14.0 && z_pos < 22.0 ) {
        //       // const static double sigma = 0.001;
        //       // pump = heat_pulse * sigma * 7.5e6 * exp(-0.5*sigma*sigma*double((current_time_step - 4000)*(current_time_step - 4000))); // AJ/fs
        //       // particle_heat = pump * dt / double(conduction_electrons); // AJ / particle
        //       // atom_potential[e] += particle_heat;
        //   if(applied_voltage_sim) electron_applied_voltage(e, array_index);
        // }
    }
}

void update_dynamics() {
  
    int array_index;
    const static double photon_energy = sim::applied_voltage*constants::eV_to_AJ; //AJ/hv
    //AJ/fs/nm**2 -> [1e-3volume(A**3)/(AJ/hv)] hv/fs 
    const static double photon_rate = 1e-2*power_density*lattice_width*lattice_depth/photon_energy; //hv/fs
    double photons_at_dt = 0.0; //hv*dt
    double pump = 0.0; // AJ / fs
    double external_potential = 0.0; //AJ/e-
    const static double sigma = dt * 0.1;
    int count = 0;
   
    if(!equilibrium_step) {
      if(heat_pulse_sim) {
        //hv(dt)/fs
        photons_at_dt = int(round(photon_rate*dt*exp(-0.5*sigma*sigma*((double(current_time_step) - ((40.0/dt)+sim::equilibration_time))*(double(current_time_step) - ((40.0 / dt)+sim::equilibration_time)))))); // AJ/fs/nm**3
        pump = 1e3*photons_at_dt*photon_energy/(dt*lattice_depth*lattice_height*lattice_width); //AJ/fs/nm**3
        external_potential = photon_energy; //AJ/hv     ;//1e27*pump*dt/ n_f; // AJ / particle
     
        TTMe = d_TTMe;
        TTMp = d_TTMp;
          //AJ/fs/K -> g(T-T)=C/t
        d_TTMe = ((G*(TTMp - TTMe)+pump)*dt*e_heat_capacity_i*300.0/TTMe) + TTMe;
       
        d_TTMp = ( G*(TTMe - TTMp)      *dt*a_heat_capacity_i)      + TTMp;
      }
    }

    #pragma omp parallel for private(array_index) schedule(dynamic, 10) 
    for (int e = 0; e < conduction_electrons; e++) {
      array_index = 3*e;        

      if(current_time_step % full_int_var == 0) e_e_coulomb(e, array_index); 
      else if (current_time_step % half_int_var == 0 && electron_transport_list[e]) e_e_coulomb(e, array_index); 
      else if(electron_transport_list[e]) neighbor_e_e_coulomb(e, array_index);
      
      if((photons_at_dt > count) && (omp_int_random[omp_get_thread_num()]() % conduction_electrons < photons_at_dt) \
      && (return_phonon_distribution((electron_potential[e] + 25.0 - E_f_A)/E_f_A, constants::kB_r*Te/E_f_A) < 1.0) \
      && electron_ea_scattering_list[e][0] != 1) {
        
      electron_thermal_field(e, array_index, external_potential);
        
        #pragma omp atomic
        count++;
      }

      //  if(!equilibrium_step) electron_applied_voltage(e, array_index, external_potential);
      if(!equilibrium_step && electron_transport_list[e]) ea_scattering(e, array_index);
   } 
  //  TEKE += external_potential;
      if(err::check) std::cout << " nearest neighbor complete. ee prep" << std::endl;
    ee_scattering();

  Tp = Tp +  a_heat_capacity_i*1e-27*TLE *n_f/lattice_atoms;
  if(Te > 1.0) Te = Te + (e_heat_capacity_i*1e-27*TEKE*n_f*300.0/conduction_electrons/Te) + (e_heat_capacity_i*pump*dt*300.0/Te);
  else         Te = Te + (e_heat_capacity_i*1e-27*TEKE*n_f*300.0/conduction_electrons)    + (e_heat_capacity_i*pump*dt*300.0);
        if (err::check) std::cout << "reset scattering." << std::endl;
    
  
    for(int e = 0; e < conduction_electrons; e++) {
      electron_ee_scattering_list[e][0] = 0;
    }   
}

void electron_thermal_field(const int e, const int array_index, const double EKE) {
        
        //  old_vel += EKE;
   // external_interaction_list[e] = false;
    // int array_index_y = array_index + 1;
    // int array_index_z = array_index + 2;
   // double x_vel = electron_velocity[array_index];//   + ((electron_force[array_index]   + new_electron_force[array_index])   * dt  * constants::K_A / 2); 
    //double y_vel = electron_velocity[array_index+1];// + ((electron_force[array_index_y] + new_electron_force[array_index_y]) * dt  * constants::K_A / 2);
    //double z_vel = electron_velocity[array_index+2];// + ((electron_force[array_index_z] + new_electron_force[array_index_z]) * dt  * constants::K_A / 2);
 //   double vel = sqrt((x_vel*x_vel)+(y_vel*y_vel)+(z_vel*z_vel));
  //  double theta = atan(y_vel / x_vel);
   // double phi = acos(z_vel / vel);
    //if(x_vel < 0.0) theta += M_PI;
    double theta = omp_uniform_random[omp_get_thread_num()]() * 2.0 * M_PI;
    double phi   = omp_uniform_random[omp_get_thread_num()]() * M_PI; 

    electron_potential[e] += EKE;
    
    double vel = sqrt(2.0*electron_potential[e]*constants::m_e_r_i);
      if(electron_potential[e] < 0.9817*E_f_A) std::cout << electron_potential[e] << ", " << EKE << ", " << vel << std::endl;
    electron_velocity[array_index]   = vel*cos(theta)*sin(phi);
    electron_velocity[array_index+1] = vel*sin(theta)*sin(phi);
    electron_velocity[array_index+2] = vel*cos(phi);

  //  electron_ee_scattering_list[e][0] = 1;
    electron_ea_scattering_list[e][0] = 1;
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
//         x_distance = electron_position[array_index]     - atom_position[array_index_a];
//         y_distance = electron_position[array_index + 1] - atom_position[array_index_a + 1];
//         z_distance = electron_position[array_index + 2] - atom_position[array_index_a + 2]; 
       
//         if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
//         else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

//         if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
//         else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

//         if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
//         else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

//         length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
        
//         if(length > e_a_neighbor_cutoff) continue;

//         atomic_nearest_electron_list[e][nearest_electron_count] = array_index_a;
//         nearest_electron_count++;

//         if (length > e_a_coulomb_cutoff) continue;
     
//         if(mean_radius[2*e] > length) mean_radius[2*e] = length;

//         if(ea_coupling) {
//             electron_ea_scattering_list[e][count] = a;
//             count++;
//         }
//   }

//  atomic_nearest_electron_list[e][0] = nearest_electron_count;
//  electron_ea_scattering_list[e][1] = count;
 
// }

// void neighbor_e_a_coulomb(const int& e, const int& array_index){
                     
//     double x_distance;
//     double y_distance;
//     double z_distance;
    
//     double length;
//     int size = atomic_nearest_electron_list[e][0];
//     int count = 2;
  
   
//     for (int a = 1; a < size; a++) {
      
//         int array_index_a = atomic_nearest_electron_list[e][a];
        
//         x_distance = electron_position[array_index]     - atom_position[array_index_a];
//         y_distance = electron_position[array_index + 1] - atom_position[array_index_a + 1];
//         z_distance = electron_position[array_index + 2] - atom_position[array_index_a + 2]; 
      
//         if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
//         else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

//         if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
//         else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

//         if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
//         else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

//       length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
   
//       if (length > e_a_coulomb_cutoff) continue;
     
//       if(mean_radius[2*e] > length) mean_radius[2*e] = length;
      
//       if(ea_coupling) {
//         electron_ea_scattering_list[e][count] = array_index_a / 3;
//         count++;
//       }
//     }
//    electron_ea_scattering_list[e][1] = count;
// }
*/
void e_e_coulomb(const int e, const int array_index) {
    
   // int array_index_i;
   // double x_distance,y_distance,z_distance,length;
    int integration_count = 1;
    int neighbor_count = 1;
    int scattering_count = 2;
  //if(omp_get_thread_num() == 1) std::cout << e << std::endl;
        for (int i = 0; i < conduction_electrons; i++) {
            if (i == e) continue; //no self repulsion
           
            int array_index_i = 3*i;

            double x_distance = electron_position[array_index]   - electron_position[array_index_i];
            double y_distance = electron_position[array_index+1] - electron_position[array_index_i + 1];
            double z_distance = electron_position[array_index+2] - electron_position[array_index_i + 2]; 

            if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance = x_distance + lattice_width;
            else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance = x_distance - lattice_width;

            if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance = y_distance + lattice_depth;
            else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance = y_distance - lattice_depth;

            if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance = z_distance + lattice_height;
            else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance = z_distance - lattice_height;

            double length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
            
            if(length > e_e_integration_cutoff) continue;
            electron_integration_list[e][integration_count] = array_index_i;
            integration_count++;
              if(integration_count >= electron_integration_list[e].size() - 2) {std::cout << e << ", " << integration_count << " > " << electron_integration_list[e].size() << ", " << length << std::endl;
              break; }

            if (length > e_e_neighbor_cutoff) continue;
            electron_nearest_electron_list[e][neighbor_count] = array_index_i/3;
            neighbor_count++;
          //   std::cout << e << ", " << i << ", " << length << std::endl;
    
              if(neighbor_count >= electron_nearest_electron_list[e].size() - 2) {std::cout << e << ", " << neighbor_count << " > " << electron_nearest_electron_list[e].size() << ", " << length << std::endl;
              break; }

            if(length > e_e_coulomb_cutoff) continue;
           //  std::cout << e << ", " << i << ", " << length << std::endl;
            electron_ee_scattering_list[e][scattering_count] = array_index_i/3;
            scattering_count++;
      
        }
    electron_integration_list[e][0] = integration_count;
    electron_nearest_electron_list[e][0] = neighbor_count;
    electron_ee_scattering_list[e][1] = scattering_count;
 
}

void neighbor_e_e_coulomb(const int e, const int array_index) {
    
    double x_distance,y_distance,z_distance, length;
    int size = electron_integration_list[e][0];
    int array_index_i;
    int neighbor_count = 1;
    int scattering_count = 2;
 //if(omp_get_thread_num() == 1) std::cout << e << ", neighbor" << std::endl;
    for (int i = 1; i < size; i++) {
        
        array_index_i = electron_integration_list[e][i];

        x_distance = electron_position[array_index]   - electron_position[array_index_i];
        y_distance = electron_position[array_index+1] - electron_position[array_index_i + 1];
        z_distance = electron_position[array_index+2] - electron_position[array_index_i + 2]; 
    
        if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance = x_distance + lattice_width;
        else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance = x_distance - lattice_width;

        if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance = y_distance + lattice_depth;
        else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance = y_distance - lattice_depth;

        if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance = z_distance + lattice_height;
        else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance = z_distance - lattice_height;
        
        length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
        
        if(length > e_e_neighbor_cutoff) continue;
        electron_nearest_electron_list[e][neighbor_count] = array_index_i/3;
        neighbor_count++;
       //  std::cout << e << ", " << i << ", " << length << std::endl;
    
          if(neighbor_count >= electron_nearest_electron_list[e].size() - 2) {std::cout << e << ", " << neighbor_count << " > " << electron_nearest_electron_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
            break; }

        if (length > e_e_coulomb_cutoff) continue; 
        electron_ee_scattering_list[e][scattering_count] = array_index_i/3;
        scattering_count++;
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
//             y_distance -= atom_position[array_index + 1];
//             z_distance -= atom_position[array_index + 2];

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
//             atom_position[array_index + 1] += force * sin(theta)*sin(phi) / atomic_mass;
//             atom_position[array_index + 2] += force * cos(phi) / atomic_mass;
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
           
//             array_index_i = atomic_nearest_atom_list[a][i];

//             x_distance -= atom_position[array_index_i];
//             y_distance -= atom_position[array_index_i + 1];
//             z_distance -= atom_position[array_index_i + 2]; 

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

//       if(omp_uniform_random[omp_get_thread_num()]() > exp(aa_rate*scattering_prob)) {
        
//         deltaE *= 0.5-0.05*exp(-2.0*deltaE/(E_f_A-3.0*constants::kB_r*Tp));
//         int size = atomic_nearest_atom_list[e][1] - 2;
//         size = 2+ (omp_int_random[omp_get_thread_num()]() % size);
//       //  if(size >= atomic_nearest_atom_list[e].size() -1) std::cout << "bus error: " << size << " < " << atomic_nearest_atom_list[e].size() << std::endl; 
//         int a = atomic_nearest_atom_list[e][size];
      
//        // if(a >= conduction_electrons -1) std::cout << "bus error: " << a << " < " << conduction_electrons-1 << std::endl; 
//         if(atomic_nearest_atom_list[a][0]) continue;

//         deltaE *= 0.5;
//         atom_potential[e] -= deltaE;
//         atom_potential[a] += deltaE;  

//         atomic_nearest_atom_list[e][0] = 1;
//         atomic_nearest_atom_list[a][0] = 1;
        
//         a_a_scattering_count++;
//       }
//     }
//     if (err::check) std::cout << "aa_scattering done." << std::endl;
// }
*/
void ea_scattering(const int e, const int array_index) {

    double scattering_velocity = electron_potential[e];
    double deltaE = sqrt(phonon_energy*E_f_A);
  
    if(omp_uniform_random[omp_get_thread_num()]() > exp(ea_rate/scattering_velocity)) {
     
      int size = electron_nearest_electron_list[e][0];
      int DoS1 = size - 1;
      int DoS2 = size - 1;
      for (int i = 1; i < size ; i++) {
        if((electron_potential[electron_nearest_electron_list[e][i]] < scattering_velocity + 1.2*deltaE) \
        && (electron_potential[electron_nearest_electron_list[e][i]] > scattering_velocity + 0.8*deltaE)) {
          DoS1 -= 1.0;
        }
        else if ((electron_potential[electron_nearest_electron_list[e][i]] < 1.2*scattering_velocity) \
        && (electron_potential[electron_nearest_electron_list[e][i]] > 0.8*scattering_velocity)) {
          DoS2 -= 1.0;
        }
      }
      if(DoS1 < DoS2) deltaE *= -1.0;
      if(Tp < Te) deltaE *= -1.0;
      if((scattering_velocity + deltaE) < 0.9817*E_f_A) deltaE = 0.0;

      electron_potential[e] += deltaE;

      if(electron_potential[e] > 0.99*E_f_A) electron_transport_list[e] = true;
      else electron_transport_list[e] = false;

      double unit_vec = sqrt(2.0*scattering_velocity*constants::m_e_r_i);
      double theta = atan(electron_velocity[array_index+1] / electron_velocity[array_index]) + (0.10*M_PI*omp_uniform_random[omp_get_thread_num()]() - 0.05*M_PI);
      double phi   = acos(electron_velocity[array_index+2] / unit_vec) + (0.10*M_PI*omp_uniform_random[omp_get_thread_num()]() - 0.05*M_PI);    
      if(electron_velocity[array_index+2] < 0.0) theta += M_PI;

      scattering_velocity = sqrt(2.0*electron_potential[e]*constants::m_e_r_i);
      electron_velocity[array_index]   = scattering_velocity * cos(theta)*sin(phi);
      electron_velocity[array_index+1] = scattering_velocity * sin(theta)*sin(phi);
      electron_velocity[array_index+2] = scattering_velocity * cos(phi); 
    
     // electron_ee_scattering_list[e][0] = 1;

      #pragma omp critical 
      {
      TEKE += deltaE;
      TLE -= deltaE;
      e_a_scattering_count++;
      }
    }  
}

void ee_scattering() {
         if (err::check) std::cout << "ee_scattering." << std::endl;


  for(int e = 0; e < conduction_electrons; e++) {
    
   // if(electron_transport_list[e]) std::cout << e << ", transport: " << electron_transport_list[e];
    if(!electron_transport_list[e]) continue;
  //  if(!electron_ee_scattering_list[e][0]) std::cout << ", scattered: " << electron_ee_scattering_list[e][0];
    if(electron_ee_scattering_list[e][0]) {
    //  std::cout << e << ", scattered: " << electron_ee_scattering_list[e][0] << std::endl;
      continue;
    }
    int size = electron_ee_scattering_list[e][1] - 2;
    if(size < 3) {
     // std::cout << e << ", " << size << std::endl;
      continue;
    }
    int electron_collision = 2 + (int_random() % size);
    electron_collision = electron_ee_scattering_list[e][electron_collision];
  
    if(electron_ee_scattering_list[electron_collision][0])  {
     // std::cout << e << ", scattered pair: " << electron_ee_scattering_list[electron_collision][0] << std::endl;
      continue;
    }
    electron_ee_scattering_list[e][0] = 1;
    electron_ee_scattering_list[electron_collision][0] = 1;

    double e_energy = electron_potential[e];
    double d_e_energy = electron_potential[electron_collision];
    double deltaE = uniform_random()*(e_energy - d_e_energy);
  
    if(uniform_random() < 0.5*exp(abs(deltaE)/(-25.0))) deltaE *= -1.0;
 
    double DoS1 = electron_nearest_electron_list[e][0] - 1;
    double DoS2 = electron_nearest_electron_list[electron_collision][0] - 1;
 
    if((e_energy - deltaE < 0.9817*E_f_A) || (d_e_energy + deltaE < 0.9817*E_f_A)) DoS2 = 0.0;
    else {
    #pragma omp sections 
    {
      #pragma omp section 
      {
      for (int i = 1; i < electron_nearest_electron_list[e][0]; i++) {
        if((electron_potential[electron_nearest_electron_list[e][i]] < (e_energy - fmin(0.8*deltaE, 1.2*deltaE))) \
        && (electron_potential[electron_nearest_electron_list[e][i]] > (e_energy - fmax(0.8*deltaE, 1.2*deltaE)))) {
          DoS1 -= 1.0;
        }
      }
      }
      #pragma omp section 
      {
      for (int i = 1; i < electron_nearest_electron_list[electron_collision][0]; i++) {
        if((electron_potential[electron_nearest_electron_list[electron_collision][i]] < (d_e_energy + fmax(0.8*deltaE, 1.2*deltaE))) \
        && (electron_potential[electron_nearest_electron_list[electron_collision][i]] > (d_e_energy + fmin(0.8*deltaE, 1.2*deltaE)))) {
          DoS2 -= 1.0;
        }
      }
      }
    }
    }
  //std::cout  << ", RTA: " << exp(ee_rate*DoS1*DoS2/(0.6666+(deltaE*deltaE))) << std::endl;

   if(uniform_random() > exp(ee_rate*DoS1*DoS2/(0.6666+(deltaE*deltaE)))) {
      int array_index = 3*e;

      electron_potential[e] -= deltaE;
      if(electron_potential[e] > 0.99*E_f_A) electron_transport_list[e] = true;
      else electron_transport_list[e] = false;

      double unit_vec = sqrt(2.0*e_energy*constants::m_e_r_i);
      double theta = atan(electron_velocity[array_index+1] / electron_velocity[array_index]) + (0.1*M_PI*uniform_random() - 0.05*M_PI);
      double phi   = acos(electron_velocity[array_index+2] / unit_vec) + (0.1*M_PI*uniform_random() - 0.05*M_PI);
        if(electron_velocity[array_index+2] < 0.0) theta += M_PI;
        
      double scattering_velocity = sqrt(2.0*electron_potential[e]*constants::m_e_r_i);
          if (scattering_velocity < 0.0) std::cout << scattering_velocity << ", " << deltaE << ", " << e_energy << std::endl;
      electron_velocity[array_index]   = scattering_velocity * cos(theta)*sin(phi);
      electron_velocity[array_index+1] = scattering_velocity * sin(theta)*sin(phi);
      electron_velocity[array_index+2] = scattering_velocity * cos(phi);
      
      electron_potential[electron_collision]   += deltaE;
      if(electron_potential[electron_collision] > 0.99*E_f_A) electron_transport_list[e] = true;
      else electron_transport_list[e] = false;

      unit_vec = sqrt(2.0*d_e_energy*constants::m_e_r_i);
      theta = atan(electron_velocity[3*electron_collision+1] / electron_velocity[3*electron_collision]) + (0.1*M_PI*uniform_random() - 0.05*M_PI);
      phi   = acos(electron_velocity[3*electron_collision+2] / unit_vec) + (0.1*M_PI*uniform_random() - 0.05*M_PI);
        if(electron_velocity[3*electron_collision+2] < 0.0) theta += M_PI;
      
      scattering_velocity = sqrt(2.0*electron_potential[electron_collision]*constants::m_e_r_i);
          if (scattering_velocity < 0.0) std::cout << scattering_velocity << ", " << deltaE << ", " << d_e_energy << std::endl;
      electron_velocity[3*electron_collision]   = scattering_velocity * cos(theta)*sin(phi);
      electron_velocity[3*electron_collision+1] = scattering_velocity * sin(theta)*sin(phi);
      electron_velocity[3*electron_collision+2] = scattering_velocity * cos(phi);
      
      e_e_scattering_count++;
    }
  }
  
  if (err::check) std::cout << "ee_scattering done." << std::endl;
} 

} //end CASTLE namespace

