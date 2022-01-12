// ==============================================
//   Coupled Atomistic and Spintronic Thermal Lattice Ensemble
//
//  =========       ========      ========   ============   ||           =========
// ||             ||        ||   ||               ||        ||          ||
// ||             ||        ||   ||               ||        ||          ||
// ||             ||        ||   ||               ||        ||          ||
// ||             || ====== ||   =========        ||        ||           =========
// ||             ||        ||           ||       ||        ||          ||
// ||             ||        ||           ||       ||        ||          ||
// ||             ||        ||           ||       ||        ||          ||
//  =========                     ========                   =========   =========
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
    
    TEPE = 0;
    TEKE = 0;
    TLE  = 0;
    

   // chosen_electron = 0;

            if (err::check) std::cout << "Calculating verlet integration...Time step: " << dt << std::endl;

            if(err::check) std::cout << "Initializing output files..." << std::endl;

            if (err::check) std::cout << "Updating new electron position." << std::endl;
    update_position();

            if (err::check) std::cout << "Forces, spins, and velocities update..." << std::endl;
    update_dynamics();
   // MLE += reinitialize_electron_conserve_momentum(captured_electron_list);
  //  std::cout << MPE/CASTLE_output_rate << ", " << TPE << ", " << MKE/CASTLE_output_rate << ", " << TKE << ", " << (MPE+MKE)/CASTLE_output_rate << ", " << TPE +TKE << ", " << std::endl;
            if (err::check) std::cout << "Output mean data" << std::endl;
     
    if (current_time_step % CASTLE_output_rate == 0)   output_data(); //std::cout << "x_flux: " << x_flux / CASTLE_output_rate << "\n"; x_flux = 0;
    

    //if(current_time_step > 10000) dt = 1e-3;
    //reset integration
    CASTLE_real_time += dt;
    current_time_step++;

   
    electron_position.swap(new_electron_position);
    electron_force.swap(new_electron_force);
    electron_velocity.swap(new_electron_velocity);
    electron_potential.swap(new_electron_potential);
  //  atom_potential.swap(new_atom_potential);
   // atom_position.swap(new_atom_position);
    //atom_velocity.swap(new_atom_velocity);
    //atom_force.swap(new_atom_force);
   
    return EXIT_SUCCESS;
}

void setup_output() {

    CASTLE_output_data = true;
}


void update_position(){

    int array_index,array_index_y,array_index_z, i;
    external_interaction_list_count = 0;
    double x_pos,y_pos,z_pos;
    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> phonon_transfer_chance(0,1);
    std::uniform_int_distribution<> phonon_transfer_vector(1,6);
    double excitation_constant;
    //double EKE = 0.0;
    #pragma omp parallel for private(i,array_index,array_index_y,array_index_z, x_pos,y_pos,z_pos, excitation_constant) \
    schedule(static) reduction(+:x_flux,y_flux,z_flux, external_interaction_list_count, TEKE)
    for (int e = 0; e < conduction_electrons; e++) { 

        array_index = 3*e;
        array_index_y = array_index + 1;
        array_index_z = array_index + 2;

      // std::cout << "why are you dying...." << e << "\n";
    
        x_pos = electron_position[array_index]   + (electron_velocity[array_index]   * dt);// + (electron_force[array_index]   * dt * dt * constants::K_A / 2); // x superarray component
        y_pos = electron_position[array_index_y] + (electron_velocity[array_index_y] * dt);// + (electron_force[array_index_y] * dt * dt * constants::K_A / 2); // y superarray component
        z_pos = electron_position[array_index_z] + (electron_velocity[array_index_z] * dt);// + (electron_force[array_index_z] * dt * dt * constants::K_A / 2); // z superarray component
        

      //  if(e==0) std::cout << sqrt((electron_velocity[array_index]* electron_velocity[array_index]) + (electron_velocity[array_index_y]*electron_velocity[array_index_y]) + (electron_velocity[array_index_z]*electron_velocity[array_index_z])) << ", " << \
        sqrt(((electron_velocity[array_index]   + (electron_force[array_index]  * dt * constants::K_A))*(electron_velocity[array_index]   + (electron_force[array_index]  * dt * constants::K_A))) + \
             ((electron_velocity[array_index_y] + (electron_force[array_index_y]* dt * constants::K_A))*(electron_velocity[array_index_y] + (electron_force[array_index_y]* dt * constants::K_A))) + \
             ((electron_velocity[array_index_z] + (electron_force[array_index_z]* dt * constants::K_A))*(electron_velocity[array_index_z] + (electron_force[array_index_z]* dt * constants::K_A)))) << std::endl; // x superarray component
     //  if(e==100) std::cout << e << ", " << electron_velocity[array_index] << " , " << electron_velocity[array_index + 1] << " , " << electron_velocity[array_index + 2] << std::endl;
        if (x_pos < 0.0) {
            x_pos += 40.0;
            x_flux--;
        }
        else if (x_pos > 40.0) {
            x_pos -= 40.0;
            x_flux++;
        }
      

	    if (y_pos < 0.0) {
            y_pos += 40.0;
            y_flux--;
        }
        else if (y_pos > 40.0) {
            y_pos -= 40.0;
            y_flux++;
        }

	    if (z_pos < 0.0) {
            z_pos += 40.0;
            z_flux--;
        }
        else if (z_pos > 40.0) {
            z_pos -= 40.0;
            z_flux++;
        }

        new_electron_position[array_index]   = x_pos;
        new_electron_position[array_index_y] = y_pos;
        new_electron_position[array_index_z] = z_pos;

      if(!equilibrium_step) {
          if(heat_pulse_sim) {
            if(x_pos < 22.0 && x_pos > 14.0 && y_pos > 14.0 && y_pos < 22.0 && z_pos > 14.0 && z_pos < 22.0 ) {
              external_interaction_list[e] = true;
              external_interaction_list_count++;
            }
          }
          if(applied_voltage_sim) TEKE += electron_applied_voltage(e, array_index);
      }
      ///  new_atom_position[array_index]   = atom_position[array_index]   + (atom_velocity[array_index]   * dt) + (atom_force[array_index]   * dt * dt * constants::K_A / 2); // x superarray component
       // new_atom_position[array_index_y] = atom_position[array_index_y] + (atom_velocity[array_index_y] * dt) + (atom_force[array_index_y] * dt * dt * constants::K_A / 2); // y superarray component
        //new_atom_position[array_index_z] = atom_position[array_index_z] + (atom_velocity[array_index_z] * dt) + (atom_force[array_index_z] * dt * dt * constants::K_A / 2); // z superarray component
        
      //   std::cout << e << ", " << atom_force[array_index] << ", " << atom_force[array_index_y] << ", " << atom_force[array_index_z]  << std::endl; // x superarray component
        i = atomic_nearest_atom_list[e][phonon_transfer_vector(gen)];
        excitation_constant = atom_potential[e] - atom_potential[i];
        if(excitation_constant < 0.0) continue;
        
        if(phonon_transfer_chance(gen) >  exp(-1.0*dt*excitation_constant / mu_f)) {
            if(excitation_constant > E_f_A) excitation_constant = E_f_A;
            #pragma omp critical(aa_update)
            {
            atom_potential[e] -= excitation_constant;
            atom_potential[i] += excitation_constant;
            }
        }
    }
}

void update_dynamics() {
   // double applied_electronic_field = {0.0, 0.0, 1.0, 1.0}; //x, y, z, strength
    int array_index;
    double e_x_force,e_y_force,e_z_force,EPE, EKE;
    double AKE = 0.0;
    double TEPE = 0.0;
    // double TEKE = 0.0;
   // TLE = 0.0;

    if(!equilibrium_step && heat_pulse_sim){
        
          const static double sigma = 0.001;
          double en_scale = heat_pulse * sigma * sqrt(5e7 / (constants::m_e_r * M_PI)) / double(external_interaction_list_count);
          AKE = en_scale* exp(-0.5*sigma*sigma*(current_time_step - 4000)*(current_time_step - 4000));
      //  std::cout << AKE << std::endl;
     //   std::cout << sigma << ", " << en_scale << ", " << AKE << ", " << -0.5*sigma*sigma*(current_time_step - 4000)*(current_time_step - 4000) << std::endl;
    }   
    #pragma omp parallel for private(array_index, e_x_force, e_y_force, e_z_force, EPE)\
     schedule(static) reduction(+:TEPE)
    for (int e = 0; e < conduction_electrons; e++) {
        array_index = 3*e;
        e_x_force = 0.0;
        e_y_force = 0.0;
        e_z_force = 0.0;
        EPE = 0.0;
       
        if(current_time_step % 15 == 0) {
            e_a_coulomb(e, array_index, e_x_force, e_y_force, e_z_force, EPE);
            e_e_coulomb(e, array_index, e_x_force, e_y_force, e_z_force, EPE);
          //  a_a_coulomb(e, array_index, a_x_force,a_y_force,a_z_force, LPE);
        
        } else {
         
            neighbor_e_a_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, EPE);
            neighbor_e_e_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, EPE);
            
        }
        
        new_electron_force[array_index]     = e_x_force;
        new_electron_force[array_index + 1] = e_y_force;
        new_electron_force[array_index + 2] = e_z_force;
        
        TEPE += EPE;  
    }

    #pragma omp parallel for private(array_index, EKE) schedule(static) reduction(+:TEKE)
    for(int e = 0; e < conduction_electrons; e++) {
      array_index = e*3;
      EKE = 0.0;
      update_velocity(e, array_index, EKE, AKE);
      TEKE += EKE;
    }

    ea_scattering();
    ee_scattering();

    MLE += std::accumulate(atom_potential.begin(),atom_potential.end(), 0.0);
    MEPE += TEPE;
    MEKE += TEKE;
  //  MLE  += TLE;
    //MLPE += TLPE;
   // MLKE += TLKE;
}

void update_velocity(const int& e, const int& array_index, double& EKE, const double& AKE) {
    
        int array_index_y = array_index + 1;
        int array_index_z = array_index + 2;
        double old_vel = sqrt((electron_velocity[array_index]*electron_velocity[array_index])+(electron_velocity[array_index_y]*electron_velocity[array_index_y])+(electron_velocity[array_index_z]*electron_velocity[array_index_z]));
     //   if (e == 0) std::cout << new_electron_velocity[array_index] << " " << electron_velocity[array_index] << "    " <<  electron_force[array_index]  << "    " << new_force_array[array_index]  << "    " <<  dt * 0.5  * constants::K / constants::m_e << "\n"; 
        double x_vel = electron_velocity[array_index]   + ((electron_force[array_index]   + new_electron_force[array_index])   * dt  * constants::K_A / 2); 
        double y_vel = electron_velocity[array_index_y] + ((electron_force[array_index_y] + new_electron_force[array_index_y]) * dt  * constants::K_A / 2);
        double z_vel = electron_velocity[array_index_z] + ((electron_force[array_index_z] + new_electron_force[array_index_z]) * dt  * constants::K_A / 2);
        double vel   = sqrt((x_vel*x_vel) + (y_vel*y_vel) + (z_vel*z_vel));
   
        double theta = atan(y_vel / x_vel);
        double phi = acos(z_vel / vel);
        if(x_vel < 0.0) theta += M_PI;
      
        if(!equilibrium_step && external_interaction_list[e]) {
          old_vel += AKE;
          external_interaction_list[e] = false;
        }

    new_electron_velocity[array_index]   = old_vel * cos(theta)*sin(phi);
    new_electron_velocity[array_index_y] = old_vel * sin(theta)*sin(phi);
    new_electron_velocity[array_index_z] = old_vel * cos(phi);
   
    EKE += old_vel*old_vel;
   
    
}

void e_a_coulomb(const int& e, const int& array_index, double& e_x_force, double& e_y_force, double& e_z_force, double& EPE){
                      //  double& a_x_force, double& a_y_force, double& a_z_force, double& EPE, double& LPE) {

    double x_distance;
    double y_distance;
    double z_distance;
   
    double length;
    double  force,  phi,theta, PE = 0;
    int array_index_a, nearest_electron_count = 1;
    int count = 2;
 //   int nearest_atom_count = 1;
    //bool collision = false;
  /*  std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> scattering_chance(0,1); */
  // std::uniform_int_distribution<> phonon_scattering_vector(0,phonon_size);
  //  int count = 0;
    for (int a = 0; a < lattice_atoms; a++) {

        array_index_a = 3*a;
        x_distance = new_electron_position[array_index]     - atom_position[array_index_a];
        y_distance = new_electron_position[array_index + 1] - atom_position[array_index_a + 1];
        z_distance = new_electron_position[array_index + 2] - atom_position[array_index_a + 2]; 
       
    //   if(a==100) std::cout << array_index_e / 3 << ", " << x_distance << ", " << y_distance << ", " << z_distance <<  std::endl;
        if (x_distance < -30.0)     x_distance = x_distance + 40.0;
        else if (x_distance > 30.0) x_distance = x_distance - 40.0;
        if (y_distance < -30.0)     y_distance = y_distance + 40.0;
        else if (y_distance > 30.0) y_distance = y_distance - 40.0;
        if (z_distance <  -30.0)    z_distance = z_distance + 40.0;
        else if (z_distance > 30.0) z_distance = z_distance - 40.0;  

        length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
        
        if(length > e_a_neighbor_cutoff) continue;

        atomic_nearest_electron_list[e][nearest_electron_count] = array_index_a;
        nearest_electron_count++;

       if (length > e_a_coulomb_cutoff) continue;
        length = sqrt(length);
     

        if(mean_radius[e] > length) {
            mean_radius[e] = length;
        }

        force = -1*(1/(length * length) - 8*150*exp(-8*length));           
        PE += 150*exp(-8*length) - (1 / length);
        
        phi   = acos(z_distance / length);
        theta = atan(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        e_x_force += mu_r * force * cos(theta)*sin(phi); 
        e_y_force += mu_r * force * sin(theta)*sin(phi);
        e_z_force += mu_r * force * cos(phi); 

        if(ea_coupling) {
        if(length < 3.0) {
            electron_ea_scattering_list[e][count] = array_index_a / 3;
            count++;
        }
      }
    }
  
    EPE += PE;
    atomic_nearest_electron_list[e][0] = nearest_electron_count;
    new_electron_potential[e] = PE;
     electron_ea_scattering_list[e][1] = count;
}

void neighbor_e_a_coulomb(const int& e, const int& array_index, double& e_x_force, double& e_y_force, double& e_z_force, double& EPE){

    double x_distance;
    double y_distance;
    double z_distance;
    
    double length;
    double force,  phi,theta, PE = 0;
    int array_index_a;
  
    int size = atomic_nearest_electron_list[e][0];
    int count = 2;
  
    for (int a = 1; a < size; a++) {
      
        array_index_a = atomic_nearest_electron_list[e][a];
        
        x_distance = new_electron_position[array_index]     - atom_position[array_index_a];
        y_distance = new_electron_position[array_index + 1] - atom_position[array_index_a + 1];
        z_distance = new_electron_position[array_index + 2] - atom_position[array_index_a + 2]; 
      
         if (x_distance < -30.0)    x_distance += 40.0;
        else if (x_distance > 30.0) x_distance -= 40.0;
        if (y_distance < -30.0)     y_distance += 40.0;
        else if (y_distance > 30.0) y_distance -= 40.0;
        if (z_distance <  -30.0)    z_distance += 40.0;
        else if (z_distance > 30.0) z_distance -= 40.0;  

        length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
   
        if (length > e_a_coulomb_cutoff) continue;
  
        length = sqrt(length);
      
       if(mean_radius[e] > length) {
            mean_radius[e] = length;
        }
        force = -1*(1/(length * length) - 8*150*exp(-8*length));
           
        PE += 150*exp(-8*length) - (1 / length);
        
        phi   = acos(z_distance / length);
        theta = atan(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        e_x_force += mu_r * force * cos(theta)*sin(phi); 
        e_y_force += mu_r * force * sin(theta)*sin(phi);
        e_z_force += mu_r * force * cos(phi); 
        
      if(ea_coupling) {
        if(length < 3.0) {
            electron_ea_scattering_list[e][count] = array_index_a / 3;
            count++;
        }
      }
    }
    EPE += PE;
    new_electron_potential[e] = PE;
    electron_ea_scattering_list[e][1] = count;
}

void e_e_coulomb(const int& e, const int& array_index, double& e_x_force, double& e_y_force, double& e_z_force, \
                double& EPE) {
    
    int array_index_i;
    double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
    int neighbor_count = 1;
    int count = 2;
        for (int i = 0; i < conduction_electrons; i++) {
            if (i == e) continue; //no self repulsion
           
            array_index_i = 3*i;

            x_distance = new_electron_position[array_index]   - new_electron_position[array_index_i];
            y_distance = new_electron_position[array_index+1] - new_electron_position[array_index_i + 1];
            z_distance = new_electron_position[array_index+2] - new_electron_position[array_index_i + 2]; 

            if (x_distance < -30)     x_distance = x_distance + 40;
            else if (x_distance > 30) x_distance = x_distance - 40;
            if (y_distance < -30)     y_distance = y_distance + 40;
            else if (y_distance > 30) y_distance = y_distance - 40;
            if (z_distance <  -30)    z_distance = z_distance + 40;
            else if (z_distance > 30) z_distance = z_distance - 40;

            length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
            
            if (length > e_e_neighbor_cutoff) continue;

            electron_nearest_electron_list[e][neighbor_count] = array_index_i;
            neighbor_count++;
            
            if (length > e_e_coulomb_cutoff) continue; 

            length = sqrt(length);
       // if(length < 0.11) length = 0.11;
        force = 1 / (length*length);
        PE += force * length;
            
            phi   = acos(z_distance / length);
            theta = atan(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            e_x_force += force * cos(theta)*sin(phi) * constants::m_e_r_i;
            e_y_force += force * sin(theta)*sin(phi) * constants::m_e_r_i;
            e_z_force += force * cos(phi) * constants::m_e_r_i;

        if(ee_coupling) {
          if(length < 3.0) {
            electron_ee_scattering_list[e][count] = array_index_i;
            count++;
          }
        }
      }
    electron_nearest_electron_list[e][0] = neighbor_count;
    EPE += PE/2;
    new_electron_potential[e] += PE;
    electron_ee_scattering_list[e][1] = count;
   // if(e == 100) std::cout << neighbor_count << std::endl;
}

void neighbor_e_e_coulomb(const int& e, const int& array_index, double& e_x_force, double& e_y_force, double& e_z_force, \
                                double& EPE) {
    
  double x_distance,y_distance,z_distance, length, force, theta,phi, PE = 0;
  int size = electron_nearest_electron_list[e][0]; //.size();
  int array_index_i;
  int count = 2;

  for (int i = 1; i < size; i++) {
        
    array_index_i = electron_nearest_electron_list[e][i];

    x_distance = new_electron_position[array_index]   - new_electron_position[array_index_i];
    y_distance = new_electron_position[array_index+1] - new_electron_position[array_index_i + 1];
    z_distance = new_electron_position[array_index+2] - new_electron_position[array_index_i + 2]; 
    
    if (x_distance < -30)     x_distance = x_distance + 40;
    else if (x_distance > 30) x_distance = x_distance - 40;
    if (y_distance < -30)     y_distance = y_distance + 40;
    else if (y_distance > 30) y_distance = y_distance - 40;
    if (z_distance <  -30)    z_distance = z_distance + 40;
    else if (z_distance > 30) z_distance = z_distance - 40;
        
    length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
    if (length > e_e_coulomb_cutoff) continue; 

    length = sqrt(length);
       
    force = 1 / (length*length);
    PE += force * length;

    phi   = acos(z_distance / length);
    theta = atan(y_distance / x_distance);
    if (x_distance < 0) theta += M_PI;

    e_x_force += force * cos(theta)*sin(phi) * constants::m_e_r_i;
    e_y_force += force * sin(theta)*sin(phi)* constants::m_e_r_i;
    e_z_force += force * cos(phi) * constants::m_e_r_i;

    if(ee_coupling) {
      if(length < 3.0) {
        electron_ee_scattering_list[e][count] = array_index_i;
        count++;
      }
    }
  }
  EPE += PE/2;
  new_electron_potential[e] += PE;
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

double electron_applied_voltage(const int& e, const int& array_index) {
    
  double vel = applied_voltage *  1.0e-1*dt*constants::e_A/constants::m_e_r;
  electron_velocity[array_index] += vel;

  return vel*vel*0.5*constants::m_e_r;
}

void ea_scattering() {
            std::srand(std::time(nullptr));
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_real_distribution<double> scattering_chance(0,1);
            
            std::uniform_real_distribution<double> Theta_pos_distrib(0.0,2.0*M_PI);
            std::uniform_real_distribution<double> Phi_pos_distrib(0.0,M_PI);
    
  int phonon_collision;
  int array_index;
  int atom_array;
  int size;
  double scattering_velocity;
  std::forward_list<int> scattering_reset_list;

  for(int e = 0; e < conduction_electrons; e++) {

    size = electron_ea_scattering_list[e][1];

    std::uniform_int_distribution<> phonon_scattering_vector(2,size);
    phonon_collision = phonon_scattering_vector(gen);

    atom_array = electron_ea_scattering_list[e][phonon_collision];
    if(electron_ea_scattering_list[atom_array][0]) continue;

    array_index = 3*e;
    scattering_velocity = 0.5*constants::m_e_r*((new_electron_velocity[array_index]  *new_electron_velocity[array_index])   \
                                              + (new_electron_velocity[array_index+1]*new_electron_velocity[array_index+1]) \
                                              + (new_electron_velocity[array_index+2]*new_electron_velocity[array_index+2]));
    
  //  std::cout << exp(ea_rate/sqrt(scattering_velocity)) << ", " << ea_rate << std::endl;

    if(scattering_chance(gen) > exp(ea_rate/sqrt(scattering_velocity))) {
      
      double deltaE = scattering_velocity - atom_potential[atom_array];
      if (deltaE > E_f_A) deltaE = E_f_A;
      else if(deltaE < 0.0) deltaE = fmax(E_f_A - atom_potential[atom_array], -1.0 * E_f_A);        

      double theta = Theta_pos_distrib(gen);
      double phi   = Phi_pos_distrib(gen);
      scattering_velocity = sqrt(2.0*(scattering_velocity - deltaE)*constants::m_e_r_i);
                
      double x_vec = scattering_velocity * cos(theta)*sin(phi);
      double y_vec = scattering_velocity * sin(theta)*sin(phi);
      double z_vec = scattering_velocity * cos(phi);

      new_electron_velocity[array_index]   = x_vec;
      new_electron_velocity[array_index+1] = y_vec;
      new_electron_velocity[array_index+2] = z_vec;
                       // std::cout << scattering_velocity << std::endl;
   
      
      atom_potential[atom_array] += deltaE;
      electron_ea_scattering_list[atom_array][0] = 1;
      e_a_scattering++;
      scattering_reset_list.push_front(atom_array);
    }
  }
  while(!scattering_reset_list.empty()) {
    electron_ee_scattering_list[scattering_reset_list.front()][0] = 0;
    scattering_reset_list.pop_front();
  }

}

void ee_scattering() {

    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> scattering_chance(0,1);
   
    std::uniform_real_distribution<double> theta_distrib(0.0,2.0*M_PI);
    std::uniform_real_distribution<double> phi_distrib(0.0,M_PI);

  std::forward_list<int> scattering_reset_list;
  int electron_collision;
  int array_index;
  int array_index_i;
  int size;
  for(int e = 0; e < conduction_electrons; e++) {
    
    if(electron_ee_scattering_list[e][0]) continue;

    size = electron_ee_scattering_list[e][1];
    std::uniform_int_distribution<> electron_scattering_vector(2,size);
    electron_collision = electron_scattering_vector(gen);//atomic_nearest_electron_list[e][phonon_scattering_vector(gen)];

    array_index_i = electron_ee_scattering_list[e][electron_collision];

    if(electron_ee_scattering_list[array_index_i/3][0]) continue;
    array_index = 3*e;

    double e_energy = 0.5*constants::m_e_r*((new_electron_velocity[array_index]  *new_electron_velocity[array_index])   \
                                          + (new_electron_velocity[array_index+1]*new_electron_velocity[array_index+1]) \
                                          + (new_electron_velocity[array_index+2]*new_electron_velocity[array_index+2]));
    double deltaE = e_energy - E_f_A;
      
   // if(e == 3792) std::cout << exp(ee_rate*deltaE*deltaE) << ", " << deltaE << ", " << deltaE*deltaE << ", " << ee_rate*deltaE*deltaE << ", " << e_energy << std::endl;

    if(scattering_chance(gen) > exp(ee_rate*deltaE*deltaE)) {
   
      double d_e_energy = 0.5*constants::m_e_r*((new_electron_velocity[array_index_i]  *new_electron_velocity[array_index_i])   \
                                              + (new_electron_velocity[array_index_i+1]*new_electron_velocity[array_index_i+1]) \
                                              + (new_electron_velocity[array_index_i+2]*new_electron_velocity[array_index_i+2]));
      deltaE = e_energy - d_e_energy;
               
      if(deltaE > E_f_A) deltaE = E_f_A;
      else if (deltaE < 0.0)  deltaE = fmax(E_f_A - d_e_energy, -1.0*E_f_A);         
     // if(e == 3792) std::cout << deltaE << ", " << d_e_energy << std::endl;


      double theta = theta_distrib(gen); 
      double phi = phi_distrib(gen); 
      double scattering_velocity = sqrt(2.0*(e_energy - deltaE)*constants::m_e_r_i);
      double x_vec = cos(theta)*sin(phi);
      double y_vec = sin(theta)*sin(phi);
      double z_vec = cos(phi);

      new_electron_velocity[array_index]   = scattering_velocity * x_vec;
      new_electron_velocity[array_index+1] = scattering_velocity * y_vec;
      new_electron_velocity[array_index+2] = scattering_velocity * z_vec;
        
      scattering_velocity = -1.0*sqrt(2.0*(d_e_energy + deltaE)*constants::m_e_r_i);
          
      new_electron_velocity[array_index_i]   = scattering_velocity * x_vec;
      new_electron_velocity[array_index_i+1] = scattering_velocity * y_vec;
      new_electron_velocity[array_index_i+2] = scattering_velocity * z_vec;
  
      e_e_scattering++;

      electron_ee_scattering_list[e][0] = 1;
      electron_ee_scattering_list[array_index_i/3][0] = 1;
      scattering_reset_list.push_front(e);
      scattering_reset_list.push_front(array_index_i/3);
    }
  }

  while(!scattering_reset_list.empty()) {
    electron_ee_scattering_list[scattering_reset_list.front()][0] = 0;
    scattering_reset_list.pop_front();
  }
}

} //end CASTLE namespace

