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

    #pragma omp parallel for private(array_index,array_index_y,array_index_z, x_pos,y_pos,z_pos) schedule(static) reduction(+:x_flux,y_flux,z_flux)
    for (int e = 0; e < conduction_electrons; e++){ 

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
    
   // std::cout << "power density: " << power_density << ", "<< photon_energy << ", " << photon_rate << ", ";
    if(!equilibrium_step) {
      if(heat_pulse_sim) {
        //hv(dt)/fs
        photons_at_dt = int(round(photon_rate*dt*exp(-0.5*sigma*sigma*((double(current_time_step) - (40.0 / dt))*(double(current_time_step) - (40.0 / dt)))))); // AJ/fs/nm**3
        pump = 1e3*photons_at_dt*photon_energy/(dt*lattice_depth*lattice_height*lattice_width); //AJ/fs/nm**3
        external_potential = photon_energy; //AJ/hv     ;//1e27*pump*dt/ n_f; // AJ / particle
      //  std::cout << exp(-0.5*sigma*sigma*((double(current_time_step) - (40.0 / dt))*(double(current_time_step) - (40.0 / dt)))) << ", " << photons_at_dt << ", " << pump << std::endl;
        TTMe = d_TTMe;
        TTMp = d_TTMp;
          //AJ/fs/K -> g(T-T)=C/t
        d_TTMe = ((G*(TTMp - TTMe)+pump)*dt*e_heat_capacity_i*300.0/TTMe) + TTMe;
        // else            d_TTMe = ((G*(TTMp - TTMe)+pump)*dt*e_heat_capacity_i)      + TTMe;
        d_TTMp = ( G*(TTMe - TTMp)      *dt*a_heat_capacity_i)      + TTMp;
      }
    }

    #pragma omp parallel for private(array_index) schedule(guided) 
    for (int e = 0; e < conduction_electrons; e++) {
        array_index = 3*e;        

        if(current_time_step % 15 == 0) {
          //e_a_coulomb(e, array_index);
          e_e_coulomb(e, array_index); 
        } else {
          //neighbor_e_a_coulomb(e, array_index);
          neighbor_e_e_coulomb(e, array_index);   
        }
      if((int_random() % conduction_electrons) < photons_at_dt && (return_phonon_distribution(E_f_A + 25.0, constants::kB_r*Te) < 1.0)) electron_thermal_field(e, array_index, external_potential);
        // if(int_random() % conduction_electrons < 320) 
        //if(!equilibrium_step) electron_applied_voltage(e, array_index, external_potential);
        // if(ea_coupling) 
      //  ea_scattering(e, array_index);
    } 
    TEKE += external_potential;
   
 //  if(ee_coupling)
  //#pragma omp parallel sections 
 // {
   // #pragma omp section
   // create_phonon_distribution(atom_potential, constants::kB_r*Tp/E_f_A);
     // 
    
  //  #pragma omp section
    ee_scattering();
      //aa_scattering();
 // }
  Tp = Tp +  a_heat_capacity_i*1e-27*TLE *n_f/lattice_atoms;
  Te = Te + (e_heat_capacity_i*1e-27*TEKE*n_f*300.0/conduction_electrons/Te) ;//+ (e_heat_capacity_i*pump*dt*300.0/Te);
  
        if (err::check) std::cout << "reset scattering." << std::endl;
    for(int e = 0; e < conduction_electrons; e++) {
      electron_ee_scattering_list[e][0] = 0;
     // atomic_nearest_atom_list[e][0] = 0;
    //  atomic_nearest_atom_list[e][0] = 0;
    }   
}

void electron_thermal_field(const int& e, const int& array_index, const double& EKE) {
        
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

    electron_velocity[array_index]   = vel*cos(theta)*sin(phi);
    electron_velocity[array_index+1] = vel*sin(theta)*sin(phi);
    electron_velocity[array_index+2] = vel*cos(phi);

    electron_ee_scattering_list[e][0] = 1;
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
        x_distance = electron_position[array_index]     - atom_position[array_index_a];
        y_distance = electron_position[array_index + 1] - atom_position[array_index_a + 1];
        z_distance = electron_position[array_index + 2] - atom_position[array_index_a + 2]; 
       
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
    int size = atomic_nearest_electron_list[e][0];
    int count = 2;
  
   
    for (int a = 1; a < size; a++) {
      
        int array_index_a = atomic_nearest_electron_list[e][a];
        
        x_distance = electron_position[array_index]     - atom_position[array_index_a];
        y_distance = electron_position[array_index + 1] - atom_position[array_index_a + 1];
        z_distance = electron_position[array_index + 2] - atom_position[array_index_a + 2]; 
      
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

            x_distance = electron_position[array_index]   - electron_position[array_index_i];
            y_distance = electron_position[array_index+1] - electron_position[array_index_i + 1];
            z_distance = electron_position[array_index+2] - electron_position[array_index_i + 2]; 

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
            
            // if(neighbor_count > 1361) {
            //   std::cout << neighbor_count << ", " << e << std::endl; 
            //   break;
            // }
            if(length > e_e_coulomb_cutoff) continue;

            electron_ee_scattering_list[e][count] = array_index_i/3;
            count++;
            
            // if(count > 1361) {
            //   std::cout << count << " > " << neighbor_count << ', ' << e << std::endl;
            //   break;
            // }
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

        x_distance = electron_position[array_index]   - electron_position[array_index_i];
        y_distance = electron_position[array_index+1] - electron_position[array_index_i + 1];
        z_distance = electron_position[array_index+2] - electron_position[array_index_i + 2]; 
    
        if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
        else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

        if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
        else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

        if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
        else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;
        
        length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
        if (length > e_e_coulomb_cutoff) continue; 
        
        electron_ee_scattering_list[e][count] = array_index_i/3;
        count++;

        // if(count > 1361) {
        //       std::cout << count << " > " << 502 << ", " << size << ", " << e << ", " << length << std::endl;
        //       break;
        // }
    }
    electron_ee_scattering_list[e][1] = count;
  // if(count >= electron_ee_scattering_list[e].size() - 1)  electron_ee_scattering_list[e][1] = 0;
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

double electron_applied_voltage(const int& e, const int& array_index, double& external_potential) {
  
  //const static double k = 1e-10*2.0*M_PI*sim::applied_voltage*1.60218e-19/(constants::h*3.08e10);
  const static double w = 1e-15*2.0*M_PI*sim::applied_voltage*1.60218e-19/constants::h;
  //double E_field = applied_voltage*sin(w*current_time_step)*exp(-0.5*dt * 0.1*dt * 0.1*((double(current_time_step) - (40.0 / dt))*(double(current_time_step) - (40.0 / dt))));
  double E_field = applied_voltage;
  double vel = 5e-11*dt*E_field*constants::e_A*constants::m_e_r_i; //V/m * q = F_x/kg = 0.5*a_x*dt = d_v_x
  
 // if((vel * electron_velocity[array_index]) > 0.0) {
    electron_velocity[array_index] += vel;
    double val = 0.5*constants::m_e_r*((electron_velocity[array_index]*electron_velocity[array_index])+(electron_velocity[array_index+1]*electron_velocity[array_index+1])+(electron_velocity[array_index+2]*electron_velocity[array_index+2])) - electron_potential[e];
    external_potential += val;
    electron_potential[e] = 0.5*constants::m_e_r*((electron_velocity[array_index]*electron_velocity[array_index])+(electron_velocity[array_index+1]*electron_velocity[array_index+1])+(electron_velocity[array_index+2]*electron_velocity[array_index+2]));
 // }
  return 0.0;
}

void aa_scattering() {
        if (err::check) std::cout << "aa_scattering." << std::endl;
  double scattering_prob;
  
  const static double aa_rate = -10.0*dt/sim::ee_coupling_strength;
    for(int e = 0; e < conduction_electrons; e++) {
  
      if(atomic_nearest_atom_list[e][0]) continue;
      
      double deltaE = atom_potential[e] - E_f_A; - 3.0*constants::kB_r*Tp;

     // if(deltaE < 8.0) scattering_prob = 8.0;
       scattering_prob = deltaE*deltaE;

      if(omp_uniform_random[omp_get_thread_num()]() > exp(aa_rate*scattering_prob)) {
        
        deltaE *= 0.5-0.05*exp(-2.0*deltaE/(E_f_A-3.0*constants::kB_r*Tp));
        int size = atomic_nearest_atom_list[e][1] - 2;
        size = 2+ (omp_int_random[omp_get_thread_num()]() % size);
      //  if(size >= atomic_nearest_atom_list[e].size() -1) std::cout << "bus error: " << size << " < " << atomic_nearest_atom_list[e].size() << std::endl; 
        int a = atomic_nearest_atom_list[e][size];
      
       // if(a >= conduction_electrons -1) std::cout << "bus error: " << a << " < " << conduction_electrons-1 << std::endl; 
        if(atomic_nearest_atom_list[a][0]) continue;

        deltaE *= 0.5;
        atom_potential[e] -= deltaE;
        atom_potential[a] += deltaE;  

        atomic_nearest_atom_list[e][0] = 1;
        atomic_nearest_atom_list[a][0] = 1;
        
        a_a_scattering_count++;
      }
    }
    if (err::check) std::cout << "aa_scattering done." << std::endl;
}

void ea_scattering(const int& e, const int& array_index) {

    double scattering_velocity = electron_potential[e];
    double scattering_const = M_B_distrib((scattering_velocity - E_f_A - 1.0*constants::kB_r*Te)/E_f_A - 0.1, constants::kB_r*Te/E_f_A);
    double deltaE = sqrt(phonon_energy*(E_f_A+ constants::kB_r*Te));
 
    if(omp_uniform_random[omp_get_thread_num()]() > exp(ea_rate*scattering_const)) {
     
      //deltaE = fmin(abs(scattering_velocity - E_f_A + 0.0*constants::kB_r*Te), deltaE);
      if(scattering_velocity < (E_f_A + constants::kB_r*Te)) deltaE *= -1.0;
      if(Tp > Te) deltaE *= -1.0;

      double theta = omp_uniform_random[omp_get_thread_num()]() * 2.0 * M_PI;
      double phi   = omp_uniform_random[omp_get_thread_num()]() * M_PI; 
      scattering_velocity = sqrt(2.0*(scattering_velocity - deltaE)*constants::m_e_r_i);
     
      electron_potential[e] -= deltaE;
      
      electron_velocity[array_index]   = scattering_velocity * cos(theta)*sin(phi);
      electron_velocity[array_index+1] = scattering_velocity * sin(theta)*sin(phi);
      electron_velocity[array_index+2] = scattering_velocity * cos(phi);     
      
      electron_ee_scattering_list[e][0] = 1;

    //  std::cout << TEKE << ", " << TLE << ", " << e_a_scattering_count << std::endl;
      #pragma omp critical 
      {
      TEKE -= deltaE;
      TLE += deltaE;
     // atom_potential[atom_collision] += deltaE;
     // atomic_nearest_atom_list[atom_collision][0] = 1;
      e_a_scattering_count++;
      }
  }  
}

void ee_scattering() {
         if (err::check) std::cout << "ee_scattering." << std::endl;
  for(int e = 0; e < conduction_electrons; e++) {
    
    if(electron_ee_scattering_list[e][0]) continue;
    int size = electron_ee_scattering_list[e][1] - 2;
    
    if(size < 3) continue;
    int electron_collision = 2 + (omp_int_random[omp_get_thread_num()]() % size);
    //if(electron_collision > electron_ee_scattering_list[e].size()) std::cout << "is this the bus error?" << electron_collision << " < " << size << std::endl; 
    electron_collision = electron_ee_scattering_list[e][electron_collision];
   // if(electron_collision > electron_ee_scattering_list.size()) std::cout << "is this the bus error?" << electron_collision << " < " << size << " < " << conduction_electrons - 1 << std::endl; 
    if(electron_ee_scattering_list[electron_collision][0])  continue;
    
    double scattering_prob;
    double e_energy = electron_potential[e];
    double d_e_energy = electron_potential[electron_collision];
    double deltaE = omp_uniform_random[omp_get_thread_num()]()*(e_energy - d_e_energy);
    double thermal_factor = return_phonon_distribution((e_energy + deltaE)/E_f_A, constants::kB_r*Te/E_f_A) - return_phonon_distribution((d_e_energy - deltaE)/E_f_A, constants::kB_r*Te/E_f_A);
   // if(e_energy < E_f_A - (3.0*constants::kB_r*Te)) std::cout << e_energy << ", " << E_f_A - 3.0*constants::kB_r*Te << "< " << e << std::endl;
    if(abs(thermal_factor) > 0.001) {
      std::cout << e_energy << ", " << d_e_energy << ", " << deltaE << ", " << abs(thermal_factor) << ", " << return_phonon_distribution(e_energy + deltaE, constants::kB_r*Te) << ", " << return_phonon_distribution(d_e_energy - deltaE, constants::kB_r*Te) << std::endl;
    //if(omp_uniform_random[omp_get_thread_num()]() > exp(ee_rate*deltaE*deltaE)) {

    
      //if(deltaE < 0.0) std::cout << "ee DeltaE: " << deltaE << " > " << 0.5*(e_energy - E_f_A + 3.0*constants::kB_r*Te) << "; " << e_energy - deltaE << " < " << E_f_A - 3.0*constants::kB_r*Te << std::endl;
      
   //   deltaE *= 0.5 - (0.5*exp(-8.0*constants::kB_r*Te/abs(deltaE)));
      
        
        int array_index = 3*e;
        double theta = omp_uniform_random[omp_get_thread_num()]()*2.0*M_PI; //theta_distrib(gen); 
        double phi   = omp_uniform_random[omp_get_thread_num()]()*M_PI; //phi_distrib(gen); 

        double scattering_velocity = sqrt(2.0*(e_energy - deltaE)*constants::m_e_r_i);
        double x_vec = cos(theta)*sin(phi);
        double y_vec = sin(theta)*sin(phi);
        double z_vec = cos(phi);

        electron_potential[e] -= deltaE;
        electron_velocity[array_index]   = scattering_velocity * x_vec;
        electron_velocity[array_index+1] = scattering_velocity * y_vec;
        electron_velocity[array_index+2] = scattering_velocity * z_vec;
        
        scattering_velocity = -1.0*sqrt(2.0*(d_e_energy + deltaE)*constants::m_e_r_i);
          
        electron_potential[electron_collision]   += deltaE;
        electron_velocity[3*electron_collision]   = scattering_velocity * x_vec;
        electron_velocity[3*electron_collision+1] = scattering_velocity * y_vec;
        electron_velocity[3*electron_collision+2] = scattering_velocity * z_vec;
  
        e_e_scattering_count++;

        electron_ee_scattering_list[e][0] = 1;
        electron_ee_scattering_list[electron_collision][0] = 1;
    }
  }
  
  if (err::check) std::cout << "ee_scattering done." << std::endl;
} 

} //end CASTLE namespace

