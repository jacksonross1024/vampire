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
    

    if(current_time_step > 1000) dt = 2e-3;
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
    long double x_pos,y_pos,z_pos;
    #pragma omp parallel for private(array_index,array_index_y,array_index_z, x_pos,y_pos,z_pos) schedule(static) reduction(+:x_flux,y_flux,z_flux)
    for (int e = 0; e < conduction_electrons; e++){ 


        array_index = 3*e;
        array_index_y = array_index + 1;
        array_index_z = array_index + 2;

      // std::cout << "why are you dying...." << e << "\n";
    
        x_pos = electron_position[array_index]   + (electron_velocity[array_index]   * dt) + (electron_force[array_index]   * dt * dt * constants::K_A / (2*constants::m_e_r)); // x superarray component
        y_pos = electron_position[array_index_y] + (electron_velocity[array_index_y] * dt) + (electron_force[array_index_y] * dt * dt * constants::K_A / (2*constants::m_e_r)); // y superarray component
        z_pos = electron_position[array_index_z] + (electron_velocity[array_index_z] * dt) + (electron_force[array_index_z] * dt * dt * constants::K_A / (2*constants::m_e_r)); // z superarray component
        
       // if(e==0) std::cout << x_pos << ", " << electron_position[array_index] << ", " << (electron_velocity[array_index]   * dt) << ", " << (electron_force[array_index]   * dt * dt * constants::K_A / (2*constants::m_e_r)) << std::endl; // x superarray component
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

        x_pos = atom_position[array_index]   + (atom_velocity[array_index]   * dt) + (atom_force[array_index]   * dt * dt * constants::K_A / (2*atomic_mass)); // x superarray component
        y_pos = atom_position[array_index_y] + (atom_velocity[array_index_y] * dt) + (atom_force[array_index_y] * dt * dt * constants::K_A / (2*atomic_mass)); // y superarray component
        z_pos = atom_position[array_index_z] + (atom_velocity[array_index_z] * dt) + (atom_force[array_index_z] * dt * dt * constants::K_A / (2*atomic_mass)); // z superarray component
        
      //  if(e==0) std::cout << x_pos << ", " << atom_position[array_index] << ", " << (atom_velocity[array_index]   * dt) << ", " << (atom_force[array_index]   * dt * dt * constants::K_A / (2*atomic_mass)) << std::endl; // x superarray component
        new_atom_position[array_index]   = x_pos;
        new_atom_position[array_index_y] = y_pos;
        new_atom_position[array_index_z] = z_pos;
    } 
}

void update_dynamics() {
   // double applied_electronic_field = {0.0, 0.0, 1.0, 1.0}; //x, y, z, strength
    int array_index;
    long double e_x_force,e_y_force,e_z_force, a_x_force,a_y_force,a_z_force, EPE, LPE, EKE, LKE;
    long double TEPE = 0;
    long double TLPE = 0;
    long double TEKE = 0;
    long double TLKE = 0;

    #pragma omp parallel for private(array_index, e_x_force,e_y_force,e_z_force, a_x_force,a_y_force,a_z_force, EPE, LPE, EKE, LKE)\
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

        if(current_time_step % 5 == 0) {
            e_e_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, EPE);
            e_a_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, a_x_force,a_y_force,a_z_force, EPE, LPE);
            a_a_coulomb(e, array_index, a_x_force,a_y_force,a_z_force, LPE);
        
        } else {
            neighbor_e_e_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, EPE);
            neighbor_e_a_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, a_x_force,a_y_force,a_z_force, EPE, LPE);
            neighbor_a_a_coulomb(e, array_index, a_x_force,a_y_force,a_z_force, LPE);
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

void update_velocity(int array_index, long double& EKE, long double& LKE) {
        int array_index_y = array_index + 1;
        int array_index_z = array_index + 2;

     //   if (e == 0) std::cout << new_electron_velocity[array_index] << " " << electron_velocity[array_index] << "    " <<  electron_force[array_index]  << "    " << new_force_array[array_index]  << "    " <<  dt * 0.5  * constants::K / constants::m_e << "\n"; 
        long double x_vel = new_electron_velocity[array_index]   = electron_velocity[array_index]   + ((electron_force[array_index]   + new_electron_force[array_index])   * dt  * constants::K_A / (2*constants::m_e_r)); 
        long double y_vel = new_electron_velocity[array_index_y] = electron_velocity[array_index_y] + ((electron_force[array_index_y] + new_electron_force[array_index_y]) * dt  * constants::K_A / (2*constants::m_e_r));
        long double z_vel = new_electron_velocity[array_index_z] = electron_velocity[array_index_z] + ((electron_force[array_index_z] + new_electron_force[array_index_z]) * dt  * constants::K_A / (2*constants::m_e_r));
        long double vel = (x_vel*x_vel) + (y_vel*y_vel) + (z_vel*z_vel);
        EKE += vel;

        x_vel = new_atom_velocity[array_index]   = atom_velocity[array_index]   + ((atom_force[array_index]   + new_atom_force[array_index])   * dt  * constants::K_A / (2*atomic_mass)); 
        y_vel = new_atom_velocity[array_index_y] = atom_velocity[array_index_y] + ((atom_force[array_index_y] + new_atom_force[array_index_y]) * dt  * constants::K_A / (2*atomic_mass));
        z_vel = new_atom_velocity[array_index_z] = atom_velocity[array_index_z] + ((atom_force[array_index_z] + new_atom_force[array_index_z]) * dt  * constants::K_A / (2*atomic_mass));
        vel = (x_vel*x_vel) + (y_vel*y_vel) + (z_vel*z_vel);
        LKE += vel;
}

void e_a_coulomb(const int a, const int array_index, long double& e_x_force, long double& e_y_force, long double& e_z_force, \
                          long double& a_x_force, long double& a_y_force, long double& a_z_force, long double& EPE, long double& LPE) {

    long double  x_distance,y_distance,z_distance, length, x_force, y_force, z_force, force, phi,theta, PE = 0;
    int array_index_e, nearest_electron_count = 1;

    x_distance = new_atom_position[array_index];
    y_distance = new_atom_position[array_index+1];
    z_distance = new_atom_position[array_index+2];

    for (int e = 0; e < conduction_electrons; e++) {

        array_index_e = 3*e;
        x_distance -= new_electron_position[array_index_e];
        y_distance -= new_electron_position[array_index_e +1];
        z_distance -= new_electron_position[array_index_e +2]; 

        if (x_distance < -30)     x_distance = x_distance + 40;
        else if (x_distance > 30) x_distance = x_distance - 40;
        if (y_distance < -30)     y_distance = y_distance + 40;
        else if (y_distance > 30) y_distance = y_distance - 40;
        if (z_distance <  -30)    z_distance = z_distance + 40;
        else if (z_distance > 30) z_distance = z_distance - 40; 

        length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms

        if(length > e_a_neighbor_cutoff) continue;

        atomic_nearest_electron_list[a][nearest_electron_count] = array_index_e;
        nearest_electron_count++;

        if (length > e_a_coulomb_cutoff) continue;

        length = sqrtl(length);
  
        force = -24*((3.3*3.3* expl(-3.3 * length)) - (expl(-1 * length))) + (1 / (length*length));
        PE += 24*((3.3*expl(-3.3*length)) - (expl(-1*length))) - (1 / length);
        

        phi   = acosl(z_distance / length);
        theta = atanl(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        x_force = force * cosl(theta)*sinl(phi); 
        y_force = force * sinl(theta)*sinl(phi);
        z_force = force * cosl(phi);

        a_x_force += x_force;
        a_y_force += y_force;
        a_z_force += z_force;

        e_x_force -= x_force;
        e_y_force -= y_force;
        e_z_force -= z_force;
       /* a_x_force += force * cosl(theta)*sinl(phi);
        a_y_force += force * sinl(theta)*sinl(phi);
        a_z_force += force * cosl(phi);

        e_x_force += -1*force * cosl(theta)*sinl(phi);
        e_y_force += -1*force * sinl(theta)*sinl(phi);
        e_z_force += -1*force * cosl(phi); */

    }
    LPE += PE /2;
    EPE += PE /2;
    atomic_nearest_electron_list[a][0] = nearest_electron_count;
}

void neighbor_e_a_coulomb(const int a, const int array_index, long double& e_x_force, long double& e_y_force, long double& e_z_force, \
                          long double& a_x_force, long double& a_y_force, long double& a_z_force, long double& EPE, long double& LPE) {

    long double  x_distance,y_distance,z_distance, length, x_force, y_force, z_force, force,  phi,theta, PE = 0;
    int array_index_e;

    int size = 1 + atomic_nearest_electron_list[a][0];
    
    x_distance = new_atom_position[array_index];
    y_distance = new_atom_position[array_index + 1];
    z_distance = new_atom_position[array_index + 2];

    for (int e = 1; e < size; e++) {
       // if(nearest_atom_list[e][a] < 0) break;
        array_index_e = atomic_nearest_electron_list[a][e];
        
        x_distance -= new_electron_position[array_index_e];
        y_distance -= new_electron_position[array_index_e + 1];
        z_distance -= new_electron_position[array_index_e + 2]; 

        if (x_distance < -30)     x_distance = x_distance + 40;
        else if (x_distance > 30) x_distance = x_distance - 40;
        if (y_distance < -30)     y_distance = y_distance + 40;
        else if (y_distance > 30) y_distance = y_distance - 40;
        if (z_distance <  -30)    z_distance = z_distance + 40;
        else if (z_distance > 30) z_distance = z_distance - 40; 

        length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms

        if (length > e_a_coulomb_cutoff) continue;

        length = sqrtl(length);
  
        force = -24*((3.3*3.3* expl(-3.3 * length)) - (expl(-1 * length))) + (1 / (length*length));
        PE += 24*((3.3*expl(-3.3*length)) - (expl(-1*length))) - (1 / length);
        
        phi   = acosl(z_distance / length);
        theta = atanl(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        x_force = force * cosl(theta)*sinl(phi); 
        y_force = force * sinl(theta)*sinl(phi);
        z_force = force * cosl(phi);

        a_x_force += x_force;
        a_y_force += y_force;
        a_z_force += z_force;

        e_x_force -= x_force;
        e_y_force -= y_force;
        e_z_force -= z_force;
       /* a_x_force += force * cosl(theta)*sinl(phi);
        a_y_force += force * sinl(theta)*sinl(phi);
        a_z_force += force * cosl(phi);

        e_x_force += -1*force * cosl(theta)*sinl(phi);
        e_y_force += -1*force * sinl(theta)*sinl(phi);
        e_z_force += -1*force * cosl(phi); */
    }
    EPE += PE/2;
    LPE += PE/2;
}

void e_e_coulomb(const int e, const int array_index, long double& e_x_force, long double& e_y_force, long double& e_z_force, \
                long double& EPE) {
    
    int array_index_i;
    long double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
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
            
            if (length > e_e_neighbor_cutoff) continue;

            electron_nearest_electron_list[e][neighbor_count] = array_index_i;
            neighbor_count++;
            
            if (length > e_e_coulomb_cutoff) continue; 

            force = 1 / length;
            length = sqrtl(length);
            PE += force * length;
            
            phi   = acosl(z_distance / length);
            theta = atanl(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            e_x_force += force * cosl(theta)*sinl(phi);
            e_y_force += force * sinl(theta)*sinl(phi);
            e_z_force += force * cosl(phi);

        }
    electron_nearest_electron_list[e][0] = neighbor_count;
    EPE += PE/2;
}

void neighbor_e_e_coulomb(const int e, const int array_index, long double& e_x_force, long double& e_y_force, long double& e_z_force, \
                                long double& EPE) {
    
    long double x_distance,y_distance,z_distance, length, force, theta,phi, PE = 0;
    int size = 1+electron_nearest_electron_list[e][0]; //.size();
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
        PE += force * length;

        phi   = acosl(z_distance / length);
        theta = atanl(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        e_x_force += force * cosl(theta)*sinl(phi);
        e_y_force += force * sinl(theta)*sinl(phi);
        e_z_force += force * cosl(phi);

    }
    EPE += PE/2;
}

void a_a_coulomb(const int a, const int array_index, \
                          long double& a_x_force, long double& a_y_force, long double& a_z_force, long double& LPE) {
    
    int array_index_i;
    long double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
    int neighbor_count = 1;

        x_distance = new_atom_position[array_index];
        y_distance = new_atom_position[array_index+1];
        z_distance = new_atom_position[array_index+2];

        for (int i = 0; i < lattice_atoms; i++) {
            if (i == a) continue; //no self repulsion
           
            array_index_i = 3*i;

            x_distance -= new_atom_position[array_index_i];
            y_distance -= new_atom_position[array_index_i + 1];
            z_distance -= new_atom_position[array_index_i + 2]; 

            if (x_distance < -30)     x_distance = x_distance + 40;
            else if (x_distance > 30) x_distance = x_distance - 40;
            if (y_distance < -30)     y_distance = y_distance + 40;
            else if (y_distance > 30) y_distance = y_distance - 40;
            if (z_distance <  -30)    z_distance = z_distance + 40;
            else if (z_distance > 30) z_distance = z_distance - 40;

            length = ((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance));
            
            if (length > a_a_neighbor_cutoff) continue;

            atomic_nearest_atom_list[a][neighbor_count] = array_index_i;
            neighbor_count++;
             
            if (length > a_a_coulomb_cutoff) continue; 

            length = sqrtl(length);
            force = 2.1*3.3*27*exp(-3.3*length) - (1 / (length*length));
            PE = 2.1*27*exp(-3.3*length) - (1 / length);
            
            phi   = acosl(z_distance / length);
            theta = atanl(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            a_x_force += force * cosl(theta)*sinl(phi);
            a_y_force += force * sinl(theta)*sinl(phi);
            a_z_force += force * cosl(phi);

        }
    atomic_nearest_atom_list[a][0] = neighbor_count;
    LPE += PE/2;
}

void neighbor_a_a_coulomb(const int a, const int array_index, \
                          long double& a_x_force, long double& a_y_force, long double& a_z_force, long double& LPE) {
    
    int array_index_i;
    long double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
    long double size = 1+atomic_nearest_atom_list[a][0];

        x_distance = new_atom_position[array_index];
        y_distance = new_atom_position[array_index+1];
        z_distance = new_atom_position[array_index+2];

        for (int i = 1; i < size; i++) {
            if (i == a) continue; //no self repulsion
           
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

            length = ((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance));
             
            if (length > a_a_coulomb_cutoff) continue; 

            length = sqrtl(length);
            force = 2.1*3.3*27*exp(-3.3*length) - (1 / (length*length));
            PE = 2.1*27*exp(-3.3*length) - (1 / length);
            
            phi   = acosl(z_distance / length);
            theta = atanl(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            a_x_force += force * cosl(theta)*sinl(phi);
            a_y_force += force * sinl(theta)*sinl(phi);
            a_z_force += force * cosl(phi);

        }
    LPE += PE/2;
}

long double electron_applied_voltage(int array_index, long double& x_force, long double& y_force, long double& z_force) {
    
  //  x_force -= 4e-9;
  //  new_electron_potential[array_index/3] += -4e-9;
    return 0;//-4e-9; //-10.0 * 1e-10 / constants::e ;
}


} //end CASTLE namespace

