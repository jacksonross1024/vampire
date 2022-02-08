// ==============================================
//   Coupled Atomistic Spin Thermalized Lattice Environment
//
//  =========     ========      ========   ============   ||           =========
// ||           ||        ||   ||               ||        ||          ||
// ||           ||        ||   ||               ||        ||          ||
// ||           ||        ||   ||               ||        ||          ||
// ||           || ====== ||   =========        ||        ||           =========
// ||           ||        ||           ||       ||        ||          ||
// ||           ||        ||           ||       ||        ||          ||
// ||           ||        ||           ||       ||        ||          ||
//  =========                   ========                   =========   =========
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


int velocity_verlet_step(double time_step) {
    
    TEKE = 0.0;
    TLE = 0.0;
            if (err::check) std::cout << "Updating new electron position." << std::endl;
    update_position();

            if (err::check) std::cout << "Forces, spins, and velocities update..." << std::endl;
    update_dynamics();
          if (err::check) std::cout << "Output mean data" << std::endl;
    
    if (current_time_step % CASTLE_output_rate == 0)   output_data(); //std::cout << "x_flux: " << x_flux / CASTLE_output_rate << "\n"; x_flux = 0;
    

  
    //reset integration
    CASTLE_real_time += dt;
    current_time_step++;

    electron_position.swap(new_electron_position);
  
   
    return EXIT_SUCCESS;
}

void setup_output() {

    CASTLE_output_data = true;
}


void update_position(){

    int array_index,array_index_y,array_index_z;
    double x_pos,y_pos,z_pos;
    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> phonon_transfer_chance(0,1);
  //  int count = 0;
    double excitation_constant;
    double pump, particle_heat;
   // external_interaction_list_count = 0;
    #pragma omp parallel for private(array_index,array_index_y,array_index_z, x_pos,y_pos,z_pos) \
    schedule(static) reduction(+:x_flux,y_flux,z_flux)
    for (int e = 0; e < conduction_electrons; e++){ 

        array_index = 3*e;
        array_index_y = array_index + 1;
        array_index_z = array_index + 2;
    
        x_pos = electron_position[array_index]   + (electron_velocity[array_index]   * dt);// + (electron_force[array_index]   * dt * dt * constants::K_A / 2); // x superarray component
        y_pos = electron_position[array_index_y] + (electron_velocity[array_index_y] * dt);// + (electron_force[array_index_y] * dt * dt * constants::K_A / 2); // y superarray component
        z_pos = electron_position[array_index_z] + (electron_velocity[array_index_z] * dt);// + (electron_force[array_index_z] * dt * dt * constants::K_A / 2); // z superarray component
        
        if (x_pos < 0.0) {x_pos += lattice_width; x_flux--;}
        else if (x_pos > lattice_width) {x_pos -= lattice_width; x_flux++;}

        if (y_pos < 0.0) {y_pos += lattice_depth; y_flux--;}
        else if (y_pos > lattice_depth) {y_pos -= lattice_depth; y_flux++;}

        if (z_pos < 0.0) {z_pos += lattice_height; z_flux--;}
        else if (z_pos > lattice_height) {z_pos -= lattice_height; z_flux++;}

        new_electron_position[array_index]   = x_pos;
        new_electron_position[array_index_y] = y_pos;
        new_electron_position[array_index_z] = z_pos;

        if(!equilibrium_step) {
          // if(heat_pulse_sim) {
          //   if(x_pos < 22.0 && x_pos > 14.0 && y_pos > 14.0 && y_pos < 22.0 && z_pos > 14.0 && z_pos < 22.0 ) {
          //     external_interaction_list[e] = true;
          //     external_interaction_list_count++;
          //   }
            // x_pos = atom_position[array_index];
            // y_pos = atom_position[array_index_y];
            // z_pos = atom_position[array_index_z];

           // if(x_pos < 22.0 && x_pos > 14.0 && y_pos > 14.0 && y_pos < 22.0 && z_pos > 14.0 && z_pos < 22.0 ) {
              const static double sigma = 0.001;
              pump = heat_pulse * sigma * 7.5e6 * exp(-0.5*sigma*sigma*double((current_time_step - 4000)*(current_time_step - 4000))); // AJ/fs
              particle_heat = pump * dt / double(conduction_electrons); // AJ / particle
         
              atom_potential[e] += particle_heat;
             
            
          
          if(applied_voltage_sim) electron_applied_voltage(e, array_index);
        }
    }
   // std::cout << count << std::endl;
    
    // std::forward_list<int> scattering_reset_list;
    
  //  #pragma omp parallel for schedule(static) reduction(+:a_a_scattering_count)
    for(int e = 0; e < conduction_electrons; e++) {
      
      const static double aa_rate = -1.0*dt/1800.0;
      double scattering_prob;
      if(atomic_nearest_atom_list[e][0]) continue;

      int size = atomic_nearest_atom_list[e][1];
      int a = atomic_nearest_atom_list[e][2];
      double max_dif = 0.0;
      // for(int i = 2; i < size; i++) {
      //   if(atomic_nearest_atom_list[i][0]) continue;         
      //   excitation_constant = atom_potential[e] - atom_potential[atomic_nearest_atom_list[e][i]];
      //   if(excitation_constant > max_dif) {
      //     max_dif = excitation_constant;
      //     a = atomic_nearest_atom_list[e][i];
      //   }
      // }

      excitation_constant = atom_potential[e] - E_f_A;
      if(excitation_constant < 8.0) scattering_prob = 8.0;
      else scattering_prob = excitation_constant*excitation_constant;

      if(phonon_transfer_chance(gen) >  exp(aa_rate*scattering_prob)) {

        std::uniform_int_distribution<> atom_collision_vector(2,size);
        a = atomic_nearest_atom_list[e][atom_collision_vector(gen)];
        if(atomic_nearest_atom_list[a][0]) continue;

        //if(max_dif > E_f_A) max_dif = E_f_A;
        excitation_constant *= 0.5;

      //  #pragma omp critical 
        {
        atom_potential[e] -= excitation_constant;
        atom_potential[a] += excitation_constant;  

        atomic_nearest_atom_list[e][0] = 1;
        atomic_nearest_atom_list[a][0] = 1;
        }
        a_a_scattering_count++;
        
      }
    }

    // #pragma omp parallel for schedule(dynamic) 
    // for(int e = 0; e< conduction_electrons; e++) {
    //   //electron_ee_scattering_list[e][0] = 0;
    //   //electron_ea_scattering_list[e][0] = 0;
    //   atomic_nearest_atom_list[e][0] = 0;
    // }
}

void update_dynamics() {
  
    int array_index;
    double pump = 0.0; // AJ / fs
    double particle_heat = 0.0;
    if(!equilibrium_step && heat_pulse_sim) {
     
      const static double sigma = 0.001;
      pump = heat_pulse * sigma * 7.5e6 * exp(-0.5*sigma*sigma*double((current_time_step - 4000)*(current_time_step - 4000))); // AJ/fs
      particle_heat = pump * dt / double(conduction_electrons); // AJ / particle
          
      TTMe = d_TTMe;
      TTMp = d_TTMp;

      if(TTMe > 1.0)  d_TTMe = ((2.15e-1*(TTMp - TTMe)+pump)*dt*e_heat_capacity)      + TTMe;
      else            d_TTMe = ((2.15e-1*(TTMp - TTMe)+pump)*dt*e_heat_capacity)      + TTMe;
                      d_TTMp = ( 2.15e-1*(TTMe - TTMp)      *dt*a_heat_capacity)      + TTMp;
    }

    #pragma omp parallel for private(array_index) schedule(static)
    for (int e = 0; e < conduction_electrons; e++) {
        array_index = 3*e;        

        if(current_time_step % 15 == 0) {
          e_a_coulomb(e, array_index);
          e_e_coulomb(e, array_index); 
        } else {
          neighbor_e_a_coulomb(e, array_index);
          neighbor_e_e_coulomb(e, array_index);   
        }
    }    
 
    #pragma omp parallel for private(array_index) schedule(static)
    for(int e = 0; e < conduction_electrons; e++) {
      array_index = 3*e;
      update_velocity(e, array_index, particle_heat);
    }

    if(ea_coupling) ea_scattering();
    if(ee_coupling) ee_scattering();

    #pragma omp parallel for reduction(+:TEKE,TLE) schedule(static)
    for(int e = 0; e < conduction_electrons; e++) {
      TEKE += electron_potential[e];
      TLE  += atom_potential[e];
    }
}

void update_velocity(const int& e, const int& array_index, const double& EKE) {
        
        //  old_vel += EKE;
   // external_interaction_list[e] = false;
    int array_index_y = array_index + 1;
    int array_index_z = array_index + 2;
      
    double x_vel = electron_velocity[array_index];//   + ((electron_force[array_index]   + new_electron_force[array_index])   * dt  * constants::K_A / 2); 
    double y_vel = electron_velocity[array_index_y];// + ((electron_force[array_index_y] + new_electron_force[array_index_y]) * dt  * constants::K_A / 2);
    double z_vel = electron_velocity[array_index_z];// + ((electron_force[array_index_z] + new_electron_force[array_index_z]) * dt  * constants::K_A / 2);

    double vel = sqrt((x_vel*x_vel)+(y_vel*y_vel)+(z_vel*z_vel));
    double theta = atan(y_vel / x_vel);
    double phi = acos(z_vel / vel);
    if(x_vel < 0.0) theta += M_PI;
    
    electron_potential[e] += EKE;

    vel = sqrt(2.0*electron_potential[e]*constants::m_e_r_i);

    electron_velocity[array_index]   = vel*cos(theta)*sin(phi);
    electron_velocity[array_index_y] = vel*sin(theta)*sin(phi);
    electron_velocity[array_index_z] = vel*cos(phi);
  
}

void e_a_coulomb(const int& e, const int& array_index){

    double x_distance;
    double y_distance;
    double z_distance;
   
    double length;
 
    int array_index_a, nearest_electron_count = 1;
    int count = 2;

  for (int a = 0; a < lattice_atoms; a++) {

        array_index_a = 3*a;
        x_distance = new_electron_position[array_index]     - atom_position[array_index_a];
        y_distance = new_electron_position[array_index + 1] - atom_position[array_index_a + 1];
        z_distance = new_electron_position[array_index + 2] - atom_position[array_index_a + 2]; 
       
        if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
        else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

        if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
        else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

        if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
        else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

      length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
        
        if(length > e_a_neighbor_cutoff) continue;

        atomic_nearest_electron_list[e][nearest_electron_count] = array_index_a;
        nearest_electron_count++;

        if (length > e_a_coulomb_cutoff) continue;
     
        if(mean_radius[2*e] > length) mean_radius[2*e] = length;

        if(ea_coupling) {
            electron_ea_scattering_list[e][count] = a;
            count++;
        }
  }

  atomic_nearest_electron_list[e][0] = nearest_electron_count;
  electron_ea_scattering_list[e][1] = count;
 
}

void neighbor_e_a_coulomb(const int& e, const int& array_index){
                     

    double x_distance;
    double y_distance;
    double z_distance;
    
    double length;
  //  double force,  phi,theta, PE = 0;
    //int array_index_a;
  //  bool collision = false;
    int size = atomic_nearest_electron_list[e][0];
    int count = 2;
  
   
    for (int a = 1; a < size; a++) {
      
        int array_index_a = atomic_nearest_electron_list[e][a];
        
        x_distance = new_electron_position[array_index]     - atom_position[array_index_a];
        y_distance = new_electron_position[array_index + 1] - atom_position[array_index_a + 1];
        z_distance = new_electron_position[array_index + 2] - atom_position[array_index_a + 2]; 
      
        if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
        else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

        if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
        else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

        if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
        else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

      length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
   
      if (length > e_a_coulomb_cutoff) continue;
     
      if(mean_radius[2*e] > length) mean_radius[2*e] = length;
      
      if(ea_coupling) {
        electron_ea_scattering_list[e][count] = array_index_a / 3;
        count++;
      }
    }
    electron_ea_scattering_list[e][1] = count;
}

void e_e_coulomb(const int& e, const int& array_index) {
    
    int array_index_i;
    double x_distance,y_distance,z_distance,length;
    int neighbor_count = 1;
    int count = 2;

        for (int i = 0; i < conduction_electrons; i++) {
            if (i == e) continue; //no self repulsion
           
            array_index_i = 3*i;

            x_distance = new_electron_position[array_index]   - new_electron_position[array_index_i];
            y_distance = new_electron_position[array_index+1] - new_electron_position[array_index_i + 1];
            z_distance = new_electron_position[array_index+2] - new_electron_position[array_index_i + 2]; 

            if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
            else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

            if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
            else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

            if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
            else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

            length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
            
            if (length > e_e_neighbor_cutoff) continue;

            electron_nearest_electron_list[e][neighbor_count] = array_index_i;
            neighbor_count++;
            
            if(length > e_e_coulomb_cutoff) continue;

            if(mean_radius[2*e + 1] > length) mean_radius[2*e + 1] = length;
            if(ee_coupling) {
                electron_ee_scattering_list[e][count] = i;
                count++;
            }
        }

    electron_nearest_electron_list[e][0] = neighbor_count;
    electron_ee_scattering_list[e][1] = count;
 
}

void neighbor_e_e_coulomb(const int& e, const int& array_index) {
    
    double x_distance,y_distance,z_distance, length;
    int size = electron_nearest_electron_list[e][0]; //.size();
    int array_index_i;
    int count = 2;

    for (int i = 1; i < size; i++) {
        
        array_index_i = electron_nearest_electron_list[e][i];

        x_distance = new_electron_position[array_index]   - new_electron_position[array_index_i];
        y_distance = new_electron_position[array_index+1] - new_electron_position[array_index_i + 1];
        z_distance = new_electron_position[array_index+2] - new_electron_position[array_index_i + 2]; 
    
        if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
        else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

        if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
        else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

        if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
        else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;
        
        length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
        if (length > e_e_coulomb_cutoff) continue; 
        
        if(mean_radius[2*e + 1] > length) mean_radius[2*e + 1] = length;
  
        if(ee_coupling) {
          electron_ee_scattering_list[e][count] = array_index_i/3;
          count++;
        }
  
    }
    electron_ee_scattering_list[e][1] = count;
   
}

void a_a_coulomb(const int a, const int array_index, \
                        double& a_x_force, double& a_y_force, double& a_z_force, double& LPE) {
   
    double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
 //   int neighbor_count = 1;

        x_distance = atom_position[array_index];
        y_distance = atom_position[array_index+1];
        z_distance = atom_position[array_index+2];
      
               // if (err::check)  std::cout << "Calculating atomic interactions" << std::endl;

       // neighbor_count = 1;

        for (int i = 0; i < lattice_atoms; i++) {
            if (i == a) continue; //no self repulsion

            x_distance -= atom_position[array_index];
            y_distance -= atom_position[array_index + 1];
            z_distance -= atom_position[array_index + 2];

            length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
           
            if(length > 6) continue;

            //length = sqrt(length);
         //   force = -2000*(length - 2);
           // PE = 1000*(length - 2)*(length - 2);
        //   if(e == 0) std::cout << force << ", " << length <<  std::endl;
 
            //phi   = acos(z_distance / length);
            //theta = atan(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            atom_position[array_index]     += force * cos(theta)*sin(phi) / atomic_mass;
            atom_position[array_index + 1] += force * sin(theta)*sin(phi) / atomic_mass;
            atom_position[array_index + 2] += force * cos(phi) / atomic_mass;
        }
   // atomic_nearest_atom_list[a][0] = neighbor_count;
    LPE += PE/2;
   // if(a == 100) std::cout << neighbor_count << std::endl;
}

void neighbor_a_a_coulomb(const int a, const int array_index, \
                        double& a_x_force, double& a_y_force, double& a_z_force, double& LPE) {
    
    int array_index_i;
    double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
    double size = atomic_nearest_atom_list[a][0];

        x_distance = atom_position[array_index];
        y_distance = atom_position[array_index+1];
        z_distance = atom_position[array_index+2];

        for (int i = 1; i < size; i++) {
           // if (i == a) continue; //no self repulsion
           
            array_index_i = atomic_nearest_atom_list[a][i];

            x_distance -= atom_position[array_index_i];
            y_distance -= atom_position[array_index_i + 1];
            z_distance -= atom_position[array_index_i + 2]; 

            if (x_distance < -30)     x_distance = x_distance + 40;
            else if (x_distance > 30) x_distance = x_distance - 40;
            if (y_distance < -30)     y_distance = y_distance + 40;
            else if (y_distance > 30) y_distance = y_distance - 40;
            if (z_distance <  -30)    z_distance = z_distance + 40;
            else if (z_distance > 30) z_distance = z_distance - 40;

            length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
             
            // if (length > a_a_coulomb_cutoff) continue; 

            length = sqrt(length);
            force = 4*(-2*exp(4 - 2*length) + 2*exp(length + 2));
            PE = 4*(exp(4 - 2*length) - (2*exp(2 -length)) + 1);
            
            phi   = acos(z_distance / length);
            theta = atan(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            a_x_force += force * cos(theta)*sin(phi) / atomic_mass;
            a_y_force += force * sin(theta)*sin(phi) / atomic_mass;
            a_z_force += force * cos(phi) / atomic_mass;

        }
    LPE += PE/2;
   // if(a == 100) std::cout << size << std::endl;
}

void electron_applied_voltage(const int& e, const int& array_index) {
    
  double vel = applied_voltage*1.0e-1*dt*constants::e_A/constants::m_e_r;
  electron_potential[e] += vel*vel*0.5*constants::m_e_r;
  electron_velocity[array_index] += vel;

}

void ea_scattering() {
            std::srand(std::time(nullptr));
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_real_distribution<double> scattering_chance(0,1);
            
            std::uniform_real_distribution<double> Theta_pos_distrib(0.0,2.0*M_PI);
            std::uniform_real_distribution<double> Phi_pos_distrib(0.0,M_PI);
 
  

 // #pragma omp parallel for schedule(guided) reduction(+:e_a_scattering_count)
  for(int e = 0; e < conduction_electrons; e++) {

    int array_index;
    int atom_array;
    int size;
    double scattering_velocity;
    double lattice_energy;
    size = electron_ea_scattering_list[e][1];

    std::uniform_int_distribution<> phonon_scattering_vector(2,size);

    atom_array = electron_ea_scattering_list[e][phonon_scattering_vector(gen)];
    if(electron_ea_scattering_list[atom_array][0]) continue;
    
    lattice_energy = atom_potential[atom_array];
    scattering_velocity = electron_potential[e];

    if(scattering_chance(gen) > exp(ea_rate/sqrt(scattering_velocity))) {
      
      array_index = 3*e;
      double deltaE = 0.5*(scattering_velocity - E_f_A);
      
      if(scattering_velocity < lattice_energy) deltaE = 0.5*(E_f_A - lattice_energy);

      double theta = Theta_pos_distrib(gen);
      double phi   = Phi_pos_distrib(gen);
      scattering_velocity = sqrt(2.0*(scattering_velocity - deltaE)*constants::m_e_r_i);

      electron_potential[e] -= deltaE;
      electron_velocity[array_index]   = scattering_velocity * cos(theta)*sin(phi);
      electron_velocity[array_index+1] = scattering_velocity * sin(theta)*sin(phi);
      electron_velocity[array_index+2] = scattering_velocity * cos(phi);     
      
    //  #pragma omp atomic
      atom_potential[atom_array] += deltaE;

    //  #pragma omp atomic write
      electron_ea_scattering_list[atom_array][0] = 1;

      electron_ee_scattering_list[e][0] = 1;

      e_a_scattering_count++;
    }
  }
}

void ee_scattering() {

    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> scattering_chance(0,1);
   
    std::uniform_real_distribution<double> theta_distrib(0.0,2.0*M_PI);
    std::uniform_real_distribution<double> phi_distrib(0.0,M_PI);

 // #pragma omp parallel for reduction(+:e_e_scattering_count) schedule(guided)
  for(int e = 0; e < conduction_electrons; e++) {
    
    int array_index;
    int size;
    int electron_collision;
    double scattering_prob;
    double e_energy;
    double deltaE;

    if(electron_ee_scattering_list[e][0]) continue;
   
    e_energy = electron_potential[e];
    deltaE = e_energy - E_f_A;
    
    if(deltaE < 8.0) scattering_prob = 8.0;
    else scattering_prob = deltaE*deltaE;

    if(scattering_chance(gen) > exp(ee_rate*scattering_prob)) {
   
      size = electron_ee_scattering_list[e][1];
      std::uniform_int_distribution<> electron_collision_vector(2,size);
      electron_collision = electron_ee_scattering_list[e][electron_collision_vector(gen)];

      if(electron_ee_scattering_list[electron_collision][0])  continue;
     
        double d_e_energy = electron_potential[electron_collision];

        array_index = 3*e;
        deltaE *= 0.5;

        double theta = theta_distrib(gen); 
        double phi = phi_distrib(gen); 
        double scattering_velocity = sqrt(2.0*(e_energy - deltaE)*constants::m_e_r_i);
        double x_vec = cos(theta)*sin(phi);
        double z_vec = cos(phi);

        double y_vec = sin(theta)*sin(phi);
     //   #pragma omp critical
        {
        electron_potential[e] -= deltaE;
        electron_velocity[array_index]   = scattering_velocity * x_vec;
        electron_velocity[array_index+1] = scattering_velocity * y_vec;
        electron_velocity[array_index+2] = scattering_velocity * z_vec;
        
        scattering_velocity = -1.0*sqrt(2.0*(d_e_energy + deltaE)*constants::m_e_r_i);
          
        
        electron_potential[electron_collision]  += deltaE;
        electron_velocity[3*electron_collision]   = scattering_velocity * x_vec;
        electron_velocity[3*electron_collision+1] = scattering_velocity * y_vec;
        electron_velocity[3*electron_collision+2] = scattering_velocity * z_vec;
  
        e_e_scattering_count++;

        electron_ee_scattering_list[e][0] = 1;
        electron_ee_scattering_list[electron_collision][0] = 1;
        }
    }
  }

  #pragma omp parallel for schedule(dynamic) 
  for(int e = 0; e< conduction_electrons; e++) {
    electron_ee_scattering_list[e][0] = 0;
    electron_ea_scattering_list[e][0] = 0;
    atomic_nearest_atom_list[e][0] = 0;
  }
} 

} //end CASTLE namespace

