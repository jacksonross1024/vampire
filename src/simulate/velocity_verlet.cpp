// ==============================================
//   Coupled Atomistic and Spintronic Thermal Lattice Ensemble
//
//         =========       ========     =========   ============   ||           =========
//        ||             ||        ||   ||               ||        ||          ||
//        ||             ||        ||   ||               ||        ||          ||
//        ||             ||        ||   ||               ||        ||          ||
//        ||             || ====== ||   =========        ||        ||           =========
//        ||             ||        ||           ||       ||        ||          ||
//        ||             ||        ||           ||       ||        ||          ||
//        ||             ||        ||           ||       ||        ||          ||
//         =========                    =========                   =========   =========
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
    TLPE = 0;
    TLKE = 0;

    chosen_electron = 0;

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
    

    if(current_time_step > 4000) dt = 1e-3;
    //reset integration
    CASTLE_real_time += dt;
    current_time_step++;

    electron_position.swap(new_electron_position);
    electron_force.swap(new_electron_force);
    electron_velocity.swap(new_electron_velocity);
    atom_position.swap(new_atom_position);
    atom_velocity.swap(new_atom_velocity);
    atom_force.swap(new_atom_force);
   
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

      // std::cout << "why are you dying...." << e << "\n";
    
        x_pos = electron_position[array_index]   + (electron_velocity[array_index]   * dt) + (electron_force[array_index]   * dt * dt * constants::K_A / 2); // x superarray component
        y_pos = electron_position[array_index_y] + (electron_velocity[array_index_y] * dt) + (electron_force[array_index_y] * dt * dt * constants::K_A / 2); // y superarray component
        z_pos = electron_position[array_index_z] + (electron_velocity[array_index_z] * dt) + (electron_force[array_index_z] * dt * dt * constants::K_A / 2); // z superarray component
        
      //  if(e==0) std::cout << x_pos << ", " << electron_position[array_index] << ", " << (electron_velocity[array_index]   * dt) << ", " << (electron_force[array_index]   * dt * dt * constants::K_A / 2) << std::endl; // x superarray component
     //  if(e==100) std::cout << e << ", " << electron_velocity[array_index] << " , " << electron_velocity[array_index + 1] << " , " << electron_velocity[array_index + 2] << std::endl;
        if (x_pos < 0) {
            x_pos += 40;
            x_flux--;
        }
        else if (x_pos > 40) {
            x_pos -= 40;
            x_flux++;
        }

	    if (y_pos < 0) {
            y_pos += 40;
            y_flux--;
        }
        else if (y_pos > 40) {
            y_pos -= 40;
            y_flux++;
        }

	    if (z_pos < 0) {
            z_pos += 40;
            z_flux--;
        }
        else if (z_pos > 40) {
            z_pos -= 40;
            z_flux++;
        }

        new_electron_position[array_index]   = x_pos;
        new_electron_position[array_index_y] = y_pos;
        new_electron_position[array_index_z] = z_pos;

        new_atom_position[array_index]   = atom_position[array_index]   + (atom_velocity[array_index]   * dt) + (atom_force[array_index]   * dt * dt * constants::K_A / 2); // x superarray component
        new_atom_position[array_index_y] = atom_position[array_index_y] + (atom_velocity[array_index_y] * dt) + (atom_force[array_index_y] * dt * dt * constants::K_A / 2); // y superarray component
        new_atom_position[array_index_z] = atom_position[array_index_z] + (atom_velocity[array_index_z] * dt) + (atom_force[array_index_z] * dt * dt * constants::K_A / 2); // z superarray component
        
      //   std::cout << e << ", " << atom_force[array_index] << ", " << atom_force[array_index_y] << ", " << atom_force[array_index_z]  << std::endl; // x superarray component
        
        
        
    } 
}

void update_dynamics() {
   // double applied_electronic_field = {0.0, 0.0, 1.0, 1.0}; //x, y, z, strength
    int array_index;
    double e_x_force,e_y_force,e_z_force, a_x_force,a_y_force,a_z_force, EPE, LPE, EKE, LKE ,x,y,z;
    double TEPE = 0;
    double TLPE = 0;
    double TEKE = 0;
    double TLKE = 0;

  //  #pragma omp parallel for private(array_index, e_x_force,e_y_force,e_z_force, a_x_force,a_y_force,a_z_force, EPE, LPE, EKE, LKE)\
     schedule(static) reduction(+:TEPE,TEKE,TLPE,TLKE)
    for (int e = 0; e < conduction_electrons; e++) {
        array_index = 3*e;
        e_x_force = 0;
        e_y_force = 0;
        e_z_force = 0;

        a_x_force = 0;
        a_y_force = 0;
        a_z_force = 0;

        EPE = 0;
        LPE = 0;
        EKE = 0;
        LKE = 0;

        x = new_atom_position[array_index];
        x = new_atom_position[array_index + 1];
        z = new_atom_position[array_index + 2];

        if(current_time_step % 2 == 0) {
            e_e_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, EPE);
            e_a_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, a_x_force,a_y_force,a_z_force, EPE, LPE);
            a_a_coulomb(e, array_index, a_x_force,a_y_force,a_z_force, LPE);
        
        } else {
            e_e_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, EPE);
           // neighbor_e_e_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, EPE);
            //e_a_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, a_x_force,a_y_force,a_z_force, EPE, LPE);
            neighbor_e_a_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, a_x_force,a_y_force,a_z_force, EPE, LPE);
            a_a_coulomb(e, array_index, a_x_force,a_y_force,a_z_force, LPE);
           // neighbor_a_a_coulomb(e, array_index, a_x_force,a_y_force,a_z_force, LPE);
        }
        
        new_electron_force[array_index]     = e_x_force;
        new_electron_force[array_index + 1] = e_y_force;
        new_electron_force[array_index + 2] = e_z_force;

        new_atom_force[array_index]     = a_x_force;
        new_atom_force[array_index + 1] = a_y_force;
        new_atom_force[array_index + 2] = a_z_force;
        
        update_velocity(array_index, EKE, LKE);
        
        TEPE += EPE;
        TEKE += EKE;
        TLPE += LPE;
        TLKE += LKE;
    }

    MEPE += TEPE;
    MEKE += TEKE;
    MLPE += TLPE;
    MLKE += TLKE;
}

void update_velocity(int array_index, double& EKE, double& LKE) {
        int array_index_y = array_index + 1;
        int array_index_z = array_index + 2;

     //   if (e == 0) std::cout << new_electron_velocity[array_index] << " " << electron_velocity[array_index] << "    " <<  electron_force[array_index]  << "    " << new_force_array[array_index]  << "    " <<  dt * 0.5  * constants::K / constants::m_e << "\n"; 
        double x_vel = new_electron_velocity[array_index]   = electron_velocity[array_index]   + ((electron_force[array_index]   + new_electron_force[array_index])   * dt  * constants::K_A / 2); 
        double y_vel = new_electron_velocity[array_index_y] = electron_velocity[array_index_y] + ((electron_force[array_index_y] + new_electron_force[array_index_y]) * dt  * constants::K_A / 2);
        double z_vel = new_electron_velocity[array_index_z] = electron_velocity[array_index_z] + ((electron_force[array_index_z] + new_electron_force[array_index_z]) * dt  * constants::K_A / 2);
        double vel = (x_vel*x_vel) + (y_vel*y_vel) + (z_vel*z_vel);
        EKE += vel;

        x_vel = new_atom_velocity[array_index]   = atom_velocity[array_index]   + ((atom_force[array_index]   + new_atom_force[array_index])   * dt  * constants::K_A / 2); 
        y_vel = new_atom_velocity[array_index_y] = atom_velocity[array_index_y] + ((atom_force[array_index_y] + new_atom_force[array_index_y]) * dt  * constants::K_A / 2);
        z_vel = new_atom_velocity[array_index_z] = atom_velocity[array_index_z] + ((atom_force[array_index_z] + new_atom_force[array_index_z]) * dt  * constants::K_A / 2);
        vel = (x_vel*x_vel) + (y_vel*y_vel) + (z_vel*z_vel);
        LKE += vel;

    if(current_time_step > 4000) {
        for(int a = 0; a < conduction_electrons; a++) {
            x_vel = electron_position[3*a];
            y_vel = electron_position[3*a+1];
            z_vel = electron_position[3*a+2];
            if((x_vel < 22 && x_vel > 14) && (y_vel > 14 && y_vel < 22) && (z_vel > 14 && z_vel < 22)) {
                new_electron_velocity[3*a]   *= 100 * 0.333333333333*sqrtl(2*E_f_A/constants::m_e_r) * exp(-1*pow(current_time_step - 40000, 2));
                new_electron_velocity[3*a+1] *= 100 * 0.333333333333*sqrtl(2*E_f_A/constants::m_e_r) * exp(-1*pow(current_time_step - 40000, 2));
                new_electron_velocity[3*a+2] *= 100 * 0.333333333333*sqrtl(2*E_f_A/constants::m_e_r) * exp(-1*pow(current_time_step - 40000, 2));
            }
        }
    }
}

void e_a_coulomb(const int a, const int& array_index, double& e_x_force, double& e_y_force, double& e_z_force, \
                        double& a_x_force, double& a_y_force, double& a_z_force, double& EPE, double& LPE) {

    double x_distance;
    double y_distance;
    double z_distance;
    double x;
    double y;
    double z;
    double length;
    double x_force, y_force, z_force, force,  phi,theta, PE = 0;
    int array_index_e, nearest_electron_count = 1;
    
    x = new_atom_position[array_index];
    x = new_atom_position[array_index + 1];
    z = new_atom_position[array_index + 2];
    int count = 0;
    for (int e = 0; e < conduction_electrons; e++) {

        array_index_e = 3*e;
        x_distance = new_atom_position[array_index] - new_electron_position[array_index_e];
        y_distance = new_atom_position[array_index + 1] - new_electron_position[array_index_e + 1];
        z_distance = new_atom_position[array_index + 2] - new_electron_position[array_index_e + 2]; 
       
       if(a==100) std::cout << array_index_e / 3 << ", " << x_distance << ", " << y_distance << ", " << z_distance <<  std::endl;
        if (x_distance < (-30.0))     x_distance = x_distance + 40.0;
        else if (x_distance > 30.0) x_distance = x_distance - 40.0;
        if (y_distance < (-30.0))     y_distance = y_distance + 40.0;
        else if (y_distance > 30.0) y_distance = y_distance - 40.0;
        if (z_distance <  (-30.0))    z_distance = z_distance + 40.0;
        else if (z_distance > 30.0) z_distance = z_distance - 40.0;  

        length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
        
        if(length > e_a_neighbor_cutoff) continue;

         if(a == 100 ) std::cout << array_index_e / 3 << ", " << new_atom_position[array_index] << ", " << new_atom_position[array_index+1] << ", " << new_atom_position[array_index+2] << ", " << sqrtl(length) << std::endl;
        if(a == 100 ) std::cout << array_index_e / 3 << ", " << new_electron_position[array_index_e] << ", " << new_electron_position[array_index_e+1] << ", " << new_electron_position[array_index_e+2] << ", " << sqrtl(length) << std::endl;
         if(a == 100 ) std::cout << array_index_e / 3 << ", " << new_atom_position[array_index] - new_electron_position[array_index_e] << ", " << new_atom_position[array_index+1] - new_electron_position[array_index_e+1]<< ", " << new_atom_position[array_index+2] - new_electron_position[array_index_e+2] << std::endl;
        // if(a == 100 ) std::cout << array_index_e / 3 << ", " << (new_atom_position[array_index] - new_electron_position[array_index_e])*(new_atom_position[array_index] - new_electron_position[array_index_e]) << ", " << (new_atom_position[array_index+1] - new_electron_position[array_index_e+1])*(new_atom_position[array_index+1] - new_electron_position[array_index_e+1]) << ", " << (new_atom_position[array_index+2] - new_electron_position[array_index_e+2])*(new_atom_position[array_index+2] - new_electron_position[array_index_e+2]) << ", " << sqrtl(length) << std::endl;
       // if(a == 100 ) std::cout << array_index_e / 3 << ", " << (new_atom_position[array_index] - new_electron_position[array_index_e])*(new_atom_position[array_index] - new_electron_position[array_index_e]) + (new_atom_position[array_index+1] - new_electron_position[array_index_e+1])*(new_atom_position[array_index+1] - new_electron_position[array_index_e+1]) + (new_atom_position[array_index+2] - new_electron_position[array_index_e+2])*(new_atom_position[array_index+2] - new_electron_position[array_index_e+2]) << ", " << sqrtl(length) << std::endl;
        //if(a == 100 ) std::cout << array_index_e / 3 << ", " << sqrt(((new_atom_position[array_index] - new_electron_position[array_index_e])*(new_atom_position[array_index] - new_electron_position[array_index_e])) + ((new_atom_position[array_index+1] - new_electron_position[array_index_e+1])*(new_atom_position[array_index+1] - new_electron_position[array_index_e+1]))+ ((new_atom_position[array_index+2] - new_electron_position[array_index_e+2])*(new_atom_position[array_index+2] - new_electron_position[array_index_e+2]))) << ", " << sqrtl(length) << std::endl;
        if(a==100) std::cout << array_index_e / 3 << ", " << x_distance << ", " << y_distance << ", " << z_distance << "\n" <<  std::endl;
        atomic_nearest_electron_list[a][nearest_electron_count] = e;
        nearest_electron_count++;

        if (length > e_a_coulomb_cutoff) continue;
        count++;

       // if(a == 100) std::cout << e << ", " << sqrtl(length) << std::endl;
        length = sqrtl(length);
       // if(a==100) std::cout << e << ", " << length << std::endl;
        if(length < 0.11) length = 0.11;

        if(mean_radius[e] > length) {
            #pragma omp critical
            mean_radius[e] = length;
        }

        force = -1.06*(1 - (27*exp(-5*length)*(5*length + 1))) / (length * length);
                        //q*k*k * exp(-15(A**-1) * length (A));
        PE += 1.06*((27*expl(-5*length) / length) - (1 / length));
        

        phi   = acosl(z_distance / length);
        theta = atanl(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        x_force = force * cosl(theta)*sinl(phi); 
        y_force = force * sinl(theta)*sinl(phi);
        z_force = force * cosl(phi);

        a_x_force += x_force * combined_mass;
        a_y_force += y_force * combined_mass;
        a_z_force += z_force * combined_mass;

        e_x_force -= x_force * mu_r;
        e_y_force -= y_force * mu_r;
        e_z_force -= z_force * mu_r;
       /* a_x_force += force * cosl(theta)*sinl(phi);
        a_y_force += force * sinl(theta)*sinl(phi);
        a_z_force += force * cosl(phi);

        e_x_force += -1*force * cosl(theta)*sinl(phi);
        e_y_force += -1*force * sinl(theta)*sinl(phi);
        e_z_force += -1*force * cosl(phi); */
     /*   if(length < 0.7 && current_time_step > 500) {
            if(electron_nearest_atom_list[e][0] == a && electron_nearest_atom_list[e][1]) continue;
            else if (electron_nearest_atom_list[e][0] != a) {
                electron_nearest_atom_list[e][0] = a;
                electron_nearest_atom_list[e][1] = false;
            }
            double excitation_constant = 1;
            if(current_time_step > 4000) excitation_constant = 0.1;
            double scattering_velocity = sqrtl((electron_velocity[array_index_e]*electron_velocity[array_index_e]) + (electron_velocity[array_index_e+1]*electron_velocity[array_index_e+1]) + (electron_velocity[array_index_e+2]*electron_velocity[array_index_e+2])) - sqrtl(2*E_f_A/constants::m_e_r);
            double excitation_energy = electron_potential[e]*constants::K_A + constants::m_e_r*0.5*((electron_velocity[array_index_e]*electron_velocity[array_index_e]) + (electron_velocity[array_index_e+1]*electron_velocity[array_index_e+1]) + (electron_velocity[array_index_e+2]*electron_velocity[array_index_e+2]));
            excitation_energy -= E_f_A;

            if (excitation_energy > 0 && scattering_velocity > 0) {
                if(scattering_chance(gen) < (1 - exp(-1*excitation_energy*excitation_constant))) {
                double vel = sqrtl((electron_velocity[array_index_e]*electron_velocity[array_index_e]) + (electron_velocity[array_index_e+1]*electron_velocity[array_index_e+1]) + (electron_velocity[array_index_e+2]*electron_velocity[array_index_e+2]));
                double theta = atanl(electron_velocity[array_index_e+1] / electron_velocity[array_index_e]);
                double phi = acosl(electron_velocity[array_index_e+2] / vel);
                if (electron_velocity[array_index_e] < 0) theta += M_PI;
                #pragma omp critical
                {
                electron_velocity[array_index_e]   = scattering_velocity * cosl(theta)*sinl(phi);
                electron_velocity[array_index_e+1] = scattering_velocity * sinl(theta)*sinl(phi);
                electron_velocity[array_index_e+2] = scattering_velocity * cosl(phi);
                
                }

                vel = sqrtl((atom_velocity[array_index_e]*atom_velocity[array_index_e]) + (atom_velocity[array_index_e+1]*atom_velocity[array_index_e+1]) + (atom_velocity[array_index_e+2]*atom_velocity[array_index_e+2]));
                theta = atanl(atom_velocity[array_index_e+1] / atom_velocity[array_index_e]);
                phi = acosl(atom_velocity[array_index_e+2] / vel);
                if(atom_velocity[array_index_e] < 0) theta += M_PI;
                scattering_velocity += sqrtl(2*E_f_A/atomic_mass);
                #pragma omp critical
                {
                chosen_electron++;
                atom_velocity[array_index_e]   = scattering_velocity * cosl(theta)*sinl(phi);
                atom_velocity[array_index_e+1] = scattering_velocity * sinl(theta)*sinl(phi);
                atom_velocity[array_index_e+2] = scattering_velocity * cosl(phi);
                }
                electron_nearest_atom_list[e][1] = false;
                }
            }
            
        } */
    }
    LPE += PE /2;
    EPE += PE /2;
    atomic_nearest_electron_list[a][0] = nearest_electron_count;
    if(a == 100) std::cout << count << std::endl;
    if(a == 100) std::cout << atomic_nearest_electron_list[a][0] << std::endl;
}

void neighbor_e_a_coulomb(const int a, const int& array_index, double& e_x_force, double& e_y_force, double& e_z_force, \
                        double& a_x_force, double& a_y_force, double& a_z_force, double& EPE, double& LPE) {

    double x_distance;
    double y_distance;
    double z_distance;
    
    double length;
    double x_force, y_force, z_force, force,  phi,theta, PE = 0;
    int array_index_e;

    int size = atomic_nearest_electron_list[a][0];
 //   if(a == 100) std::cout << atomic_nearest_electron_list[a][0] << std::endl;

    int count = 0;
    for (int e = 1; e < size; e++) {
       // if(atomic_nearest_electron_list[a][e] < 0) std::cout << a << ", " << e << std::endl;
        array_index_e = 3*atomic_nearest_electron_list[a][e];
        
        x_distance = new_atom_position[array_index] - new_electron_position[array_index_e];
        y_distance = new_atom_position[array_index + 1] - new_electron_position[array_index_e + 1];
        z_distance = new_atom_position[array_index + 2] - new_electron_position[array_index_e + 2]; 
      
         if (x_distance < (-30.0))     x_distance += 40.0;
        else if (x_distance > 30.0) x_distance -= 40.0;
        if (y_distance < (-30.0))     y_distance += 40.0;
        else if (y_distance > 30.0) y_distance -= 40.0;
        if (z_distance <  (-30.0))    z_distance += 40.0;
        else if (z_distance > 30.0) z_distance -= 40.0;  

        length = ((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance)); //Angstroms
          if(a == 100 ) std::cout << array_index_e / 3 << ", " << new_atom_position[array_index] << ", " << new_atom_position[array_index+1] << ", " << new_atom_position[array_index+2] << ", " << sqrtl(length) << std::endl;
        if(a == 100 ) std::cout << array_index_e / 3 << ", " << new_electron_position[array_index_e] << ", " << new_electron_position[array_index_e+1] << ", " << new_electron_position[array_index_e+2] << ", " << sqrtl(length) << std::endl;
         if(a == 100 ) std::cout << array_index_e / 3 << ", " << new_atom_position[array_index] - new_electron_position[array_index_e] << ", " << new_atom_position[array_index+1] - new_electron_position[array_index_e+1]<< ", " << new_atom_position[array_index+2] - new_electron_position[array_index_e+2] << ", " << sqrtl(length) << std::endl;
         //if(a == 100 ) std::cout << array_index_e / 3 << ", " << (new_atom_position[array_index] - new_electron_position[array_index_e])*(new_atom_position[array_index] - new_electron_position[array_index_e]) << ", " << (new_atom_position[array_index+1] - new_electron_position[array_index_e+1])*(new_atom_position[array_index+1] - new_electron_position[array_index_e+1]) << ", " << (new_atom_position[array_index+2] - new_electron_position[array_index_e+2])*(new_atom_position[array_index+2] - new_electron_position[array_index_e+2]) << ", " << sqrtl(length) << std::endl;
       // if(a == 100 ) std::cout << array_index_e / 3 << ", " << (new_atom_position[array_index] - new_electron_position[array_index_e])*(new_atom_position[array_index] - new_electron_position[array_index_e]) + (new_atom_position[array_index+1] - new_electron_position[array_index_e+1])*(new_atom_position[array_index+1] - new_electron_position[array_index_e+1]) + (new_atom_position[array_index+2] - new_electron_position[array_index_e+2])*(new_atom_position[array_index+2] - new_electron_position[array_index_e+2]) << ", " << sqrtl(length) << std::endl;
       // if(a == 100 ) std::cout << array_index_e / 3 << ", " << sqrtl(((new_atom_position[array_index] - new_electron_position[array_index_e])*(new_atom_position[array_index] - new_electron_position[array_index_e])) + ((new_atom_position[array_index+1] - new_electron_position[array_index_e+1])*(new_atom_position[array_index+1] - new_electron_position[array_index_e+1]))+ ((new_atom_position[array_index+2] - new_electron_position[array_index_e+2])*(new_atom_position[array_index+2] - new_electron_position[array_index_e+2]))) << ", " << sqrtl(length) << std::endl;
        if(a==100) std::cout << array_index_e / 3 << ", " << x_distance << ", " << y_distance << ", " << z_distance << "\n" << std::endl;
        if (length > 144) continue;
        count++;
        
        length = sqrtl(length);
       // if(a==100) std::cout << array_index_e / 3 << ", " << length << std::endl;

        if(length < 0.11) length = 0.11;
        if(mean_radius[e] > length) {
            #pragma omp critical
            mean_radius[e] = length;
        } /*
        if(length < 0.7) {
            double scattering_velocity = sqrtl((electron_velocity[array_index_e]*electron_velocity[array_index_e]) + (electron_velocity[array_index_e+1]*electron_velocity[array_index_e+1]) + (electron_velocity[array_index_e+2]*electron_velocity[array_index_e+2])) - sqrtl(2*E_f_A/constants::m_e_r);
            double excitation_energy = electron_potential[e]*constants::K_A + constants::m_e_r*0.5*((electron_velocity[array_index_e]*electron_velocity[array_index_e]) + (electron_velocity[array_index_e+1]*electron_velocity[array_index_e+1]) + (electron_velocity[array_index_e+2]*electron_velocity[array_index_e+2]));
            excitation_energy -= E_f_A;

            if (excitation_energy > 0 && scattering_velocity > 0) {

                double vel = sqrtl((electron_velocity[array_index_e]*electron_velocity[array_index_e]) + (electron_velocity[array_index_e+1]*electron_velocity[array_index_e+1]) + (electron_velocity[array_index_e+2]*electron_velocity[array_index_e+2]));
                double theta = atanl(electron_velocity[array_index_e+1] / electron_velocity[array_index_e]);
                double phi = acosl(electron_velocity[array_index_e+2] / vel);
                if (electron_velocity[array_index_e] < 0) theta += M_PI;
                #pragma omp critical
                {
                electron_velocity[array_index_e]   = scattering_velocity * cosl(theta)*sinl(phi);
                electron_velocity[array_index_e+1] = scattering_velocity * sinl(theta)*sinl(phi);
                electron_velocity[array_index_e+2] = scattering_velocity * cosl(phi);
                }

                vel = sqrtl((atom_velocity[array_index_e]*atom_velocity[array_index_e]) + (atom_velocity[array_index_e+1]*atom_velocity[array_index_e+1]) + (atom_velocity[array_index_e+2]*atom_velocity[array_index_e+2]));
                theta = atanl(atom_velocity[array_index_e+1] / atom_velocity[array_index_e]);
                phi = acosl(atom_velocity[array_index_e+2] / vel);
                if(atom_velocity[array_index_e] < 0) theta += M_PI;
                scattering_velocity += sqrtl(2*E_f_A/constants::m_e_r);
                #pragma omp critical
                {
                chosen_electron++;
                atom_velocity[array_index_e]   = scattering_velocity * cosl(theta)*sinl(phi);
                atom_velocity[array_index_e+1] = scattering_velocity * sinl(theta)*sinl(phi);
                atom_velocity[array_index_e+2] = scattering_velocity * cosl(phi);
                }
            }
             
        } */

        force = -1.06*(1 - (27*exp(-5*length)*(5*length + 1))) / (length * length);
                        //q*k*k * exp(-15(A**-1) * length (A));
        PE += 1.06*((27*expl(-5*length) / length) - (1 / length));
        
        phi   = acosl(z_distance / length);
        theta = atanl(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        x_force = force * cosl(theta)*sinl(phi); 
        y_force = force * sinl(theta)*sinl(phi);
        z_force = force * cosl(phi);

        a_x_force += x_force * combined_mass;
        a_y_force += y_force * combined_mass;
        a_z_force += z_force * combined_mass;

        e_x_force -= x_force * mu_r;
        e_y_force -= y_force * mu_r;
        e_z_force -= z_force * mu_r;
       /* a_x_force += force * cosl(theta)*sinl(phi);
        a_y_force += force * sinl(theta)*sinl(phi);
        a_z_force += force * cosl(phi);

        e_x_force += -1*force * cosl(theta)*sinl(phi);
        e_y_force += -1*force * sinl(theta)*sinl(phi);
        e_z_force += -1*force * cosl(phi); */
    }
    EPE += PE/2;
    LPE += PE/2;
    if(a == 100) std::cout << count << std::endl;
}

void e_e_coulomb(const int e, const int array_index, double& e_x_force, double& e_y_force, double& e_z_force, \
                double& EPE) {
    
    int array_index_i;
    double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
    int neighbor_count = 1;

        x_distance = new_electron_position[array_index];
        y_distance = new_electron_position[array_index+1];
        z_distance = new_electron_position[array_index+2];

        for (int i = 0; i < conduction_electrons; i++) {
            if (i == e) continue; //no self repulsion
           
            array_index_i = 3*i;

            x_distance -= new_electron_position[array_index_i];
            y_distance -= new_electron_position[array_index_i + 1];
            z_distance -= new_electron_position[array_index_i + 2]; 

            if (x_distance < -30)     x_distance = x_distance + 40;
            else if (x_distance > 30) x_distance = x_distance - 40;
            if (y_distance < -30)     y_distance = y_distance + 40;
            else if (y_distance > 30) y_distance = y_distance - 40;
            if (z_distance <  -30)    z_distance = z_distance + 40;
            else if (z_distance > 30) z_distance = z_distance - 40;

            length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
            
          //  if (length > e_e_neighbor_cutoff) continue;

        //    electron_nearest_electron_list[e][neighbor_count] = array_index_i;
      //      neighbor_count++;
            
            if (length > e_e_coulomb_cutoff) continue; 

            force = 1 / length;
            length = sqrtl(length);
            if(length < 0.11) length = 0.11;
            PE += force * length;
            
            phi   = acosl(z_distance / length);
            theta = atanl(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            e_x_force += force * cosl(theta)*sinl(phi) / constants::m_e_r;
            e_y_force += force * sinl(theta)*sinl(phi)/ constants::m_e_r;
            e_z_force += force * cosl(phi) / constants::m_e_r;

        }
  //  electron_nearest_electron_list[e][0] = neighbor_count;
    EPE += PE/2;
   // if(e == 100) std::cout << neighbor_count << std::endl;
}

void neighbor_e_e_coulomb(const int e, const int array_index, double& e_x_force, double& e_y_force, double& e_z_force, \
                                double& EPE) {
    
    double x_distance,y_distance,z_distance, length, force, theta,phi, PE = 0;
    int size = electron_nearest_electron_list[e][0]; //.size();
    int array_index_i;

    x_distance = new_electron_position[array_index];
    y_distance = new_electron_position[array_index+1];
    z_distance = new_electron_position[array_index+2];

    for (int i=1; i < size; i++) {
        
        array_index_i = electron_nearest_electron_list[e][i];

        x_distance -= new_electron_position[array_index_i]; //don't get caching but cut down iterations by large margin
        y_distance -= new_electron_position[array_index_i + 1];
        z_distance -= new_electron_position[array_index_i + 2];
    
        if (x_distance < -30)     x_distance = x_distance + 40;
        else if (x_distance > 30) x_distance = x_distance - 40;
        if (y_distance < -30)     y_distance = y_distance + 40;
        else if (y_distance > 30) y_distance = y_distance - 40;
        if (z_distance <  -30)    z_distance = z_distance + 40;
        else if (z_distance > 30) z_distance = z_distance - 40;
        
        length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
        if (length > e_e_coulomb_cutoff) continue; 

        force = 1 / length;
        length = sqrtl(length);
        if(length < 0.11) length = 0.11;
        PE += force * length;

        phi   = acosl(z_distance / length);
        theta = atanl(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        e_x_force += force * cosl(theta)*sinl(phi) / constants::m_e_r;
            e_y_force += force * sinl(theta)*sinl(phi)/ constants::m_e_r;
            e_z_force += force * cosl(phi) / constants::m_e_r;

    }
    EPE += PE/2;
  //  if(e == 100) std::cout << size << std::endl;
}

void a_a_coulomb(const int a, const int array_index, \
                        double& a_x_force, double& a_y_force, double& a_z_force, double& LPE) {
   
    double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
 //   int neighbor_count = 1;

        x_distance = new_atom_position[array_index];
        y_distance = new_atom_position[array_index+1];
        z_distance = new_atom_position[array_index+2];
      
               // if (err::check)  std::cout << "Calculating atomic interactions" << std::endl;

       // neighbor_count = 1;

        for (int i = 0; i < lattice_atoms; i++) {
            if (i == a) continue; //no self repulsion

            x_distance -= new_atom_position[array_index];
            y_distance -= new_atom_position[array_index + 1];
            z_distance -= new_atom_position[array_index + 2];

            length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
           
            if(length > 6) continue;

            length = sqrt(length);
            force = -2000*(length - 2);
            PE = 1000*(length - 2)*(length - 2);
        //   if(e == 0) std::cout << force << ", " << length <<  std::endl;
 
            phi   = acosl(z_distance / length);
            theta = atanl(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            new_atom_force[array_index]     += force * cos(theta)*sin(phi) / atomic_mass;
            new_atom_force[array_index + 1] += force * sin(theta)*sin(phi) / atomic_mass;
            new_atom_force[array_index + 2] += force * cos(phi) / atomic_mass;
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

        x_distance = new_atom_position[array_index];
        y_distance = new_atom_position[array_index+1];
        z_distance = new_atom_position[array_index+2];

        for (int i = 1; i < size; i++) {
           // if (i == a) continue; //no self repulsion
           
            array_index_i = atomic_nearest_atom_list[a][i];

            x_distance -= new_atom_position[array_index_i];
            y_distance -= new_atom_position[array_index_i + 1];
            z_distance -= new_atom_position[array_index_i + 2]; 

            if (x_distance < -30)     x_distance = x_distance + 40;
            else if (x_distance > 30) x_distance = x_distance - 40;
            if (y_distance < -30)     y_distance = y_distance + 40;
            else if (y_distance > 30) y_distance = y_distance - 40;
            if (z_distance <  -30)    z_distance = z_distance + 40;
            else if (z_distance > 30) z_distance = z_distance - 40;

            length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
             
            if (length > a_a_coulomb_cutoff) continue; 

            length = sqrtl(length);
            force = 4*(-2*exp(4 - 2*length) + 2*exp(length + 2));
            PE = 4*(exp(4 - 2*length) - (2*exp(2 -length)) + 1);
            
            phi   = acosl(z_distance / length);
            theta = atanl(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            a_x_force += force * cosl(theta)*sinl(phi) / atomic_mass;
            a_y_force += force * sinl(theta)*sinl(phi) / atomic_mass;
            a_z_force += force * cosl(phi) / atomic_mass;

        }
    LPE += PE/2;
   // if(a == 100) std::cout << size << std::endl;
}

double electron_applied_voltage(int array_index, double& x_force, double& y_force, double& z_force) {
    
  //  x_force -= 4e-9;
  //  new_electron_potential[array_index/3] += -4e-9;
    return 0;//-4e-9; //-10.0 * 1e-10 / constants::e ;
}


} //end CASTLE namespace

