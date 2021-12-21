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
    //atom_potential.swap(new_atom_potential);
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
    double x_pos,y_pos,z_pos;
    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> phonon_transfer_chance(0,1);
    std::uniform_int_distribution<> phonon_transfer_vector(1,6);
    double excitation_constant;
    double excitation_energy = mu_f;
    #pragma omp parallel for private(i, array_index,array_index_y,array_index_z, x_pos,y_pos,z_pos, excitation_constant) schedule(static) reduction(+:x_flux,y_flux,z_flux)
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

      ///  new_atom_position[array_index]   = atom_position[array_index]   + (atom_velocity[array_index]   * dt) + (atom_force[array_index]   * dt * dt * constants::K_A / 2); // x superarray component
       // new_atom_position[array_index_y] = atom_position[array_index_y] + (atom_velocity[array_index_y] * dt) + (atom_force[array_index_y] * dt * dt * constants::K_A / 2); // y superarray component
        //new_atom_position[array_index_z] = atom_position[array_index_z] + (atom_velocity[array_index_z] * dt) + (atom_force[array_index_z] * dt * dt * constants::K_A / 2); // z superarray component
        
      //   std::cout << e << ", " << atom_force[array_index] << ", " << atom_force[array_index_y] << ", " << atom_force[array_index_z]  << std::endl; // x superarray component
        i = atomic_nearest_atom_list[e][phonon_transfer_vector(gen)];
        excitation_constant = atom_potential[e] - atom_potential[i];
        if(excitation_constant < 0.0) continue;
        if(phonon_transfer_chance(gen) >  exp(-1.0*dt*excitation_constant / mu_f)) {
            if(excitation_constant > E_f_A) excitation_constant = E_f_A;
            #pragma omp critical
            {
           //   std::cout << excitation_constant << ", " << excitation_energy << std::endl;
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
    double TEPE = 0;
    double TEKE = 0;
    TLE = 0;

    #pragma omp parallel for private(array_index, e_x_force,e_y_force,e_z_force, EPE, EKE)\
     schedule(static) reduction(+:TEPE,TEKE, TLE)
    for (int e = 0; e < conduction_electrons; e++) {
        array_index = 3*e;
        e_x_force = 0.0;
        e_y_force = 0.0;
        e_z_force = 0.0;

        EPE = 0.0;
        EKE = 0.0;
        

        if(current_time_step % 20 == 0) {
            e_a_coulomb(e, array_index, e_x_force,e_y_force,e_z_force,EPE);
            e_e_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, EPE);
            
          //  a_a_coulomb(e, array_index, a_x_force,a_y_force,a_z_force, LPE);
        
        } else {
          //  std::cout << e << "\n";
           // e_e_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, EPE);
            neighbor_e_a_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, EPE);
           neighbor_e_e_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, EPE);
            //e_a_coulomb(e, array_index, e_x_force,e_y_force,e_z_force, a_x_force,a_y_force,a_z_force, EPE, LPE);
           
           // a_a_coulomb(e, array_index, a_x_force,a_y_force,a_z_force, LPE);
           // neighbor_a_a_coulomb(e, array_index, a_x_force,a_y_force,a_z_force, LPE);
        }
        
       // atom_potential[e] += new_atom_potential[e];
        TLE += atom_potential[e]; 
        //new_atom_potential[e] = 0;

        new_electron_force[array_index]     = e_x_force;
        new_electron_force[array_index + 1] = e_y_force;
        new_electron_force[array_index + 2] = e_z_force;

     //   new_atom_force[array_index]     = a_x_force;
       // new_atom_force[array_index + 1] = a_y_force;
        //new_atom_force[array_index + 2] = a_z_force;
        
        update_velocity(array_index, EKE);
        
        TEPE += EPE;
        TEKE += EKE;
     //   TLPE += LPE;
       // TLKE += LKE;
    }

    MEPE += TEPE;
    MEKE += TEKE;
    MLE += TLE;
    //MLPE += TLPE;
   // MLKE += TLKE;
}

void update_velocity(int array_index, double& EKE) {
        int array_index_y = array_index + 1;
        int array_index_z = array_index + 2;
        
     //   if (e == 0) std::cout << new_electron_velocity[array_index] << " " << electron_velocity[array_index] << "    " <<  electron_force[array_index]  << "    " << new_force_array[array_index]  << "    " <<  dt * 0.5  * constants::K / constants::m_e << "\n"; 
        double x_vel = electron_velocity[array_index]   + ((electron_force[array_index]   + new_electron_force[array_index])   * dt  * constants::K_A / 2); 
        double y_vel = electron_velocity[array_index_y] + ((electron_force[array_index_y] + new_electron_force[array_index_y]) * dt  * constants::K_A / 2);
        double z_vel = electron_velocity[array_index_z] + ((electron_force[array_index_z] + new_electron_force[array_index_z]) * dt  * constants::K_A / 2);
        double vel = sqrt((x_vel*x_vel)+(y_vel*y_vel)+(z_vel*z_vel));
        
    /*    if(vel*dt > 0.1) {
            
            double theta = atan(y_vel / x_vel);
            double phi = acos(z_vel / vel);
            if(x_vel < 0) theta += M_PI;
           // vel -= sqrt(2*mu_f*1e20*(1 - 0.07 / (dt*vel)) / constants::m_e_r);
            #pragma omp critical
           // new_atom_potential[electron_nearest_atom_list[array_index/3][(array_index / 3) % 5 + 10]] += sqrt(2*(vel - (0.1 / dt))/constants::m_e_r);
            //if(vel * dt > 0.1) vel = 1e-5;
            x_vel = vel*cos(theta)*sin(phi);
            y_vel = vel*sin(theta)*sin(phi);
            z_vel = vel*cos(phi);
        } */
    if(!equilibrium_step) {
        double x = new_electron_position[array_index];
        double y = new_electron_position[array_index_y];
        double z = new_electron_position[array_index_z];
        if(x < 22.0 && x > 14.0 && y > 14.0 && y < 22.0 && z > 14.0 && z < 22.0 ) {
            
            double theta = atan(y_vel / x_vel);
            double phi = acos(z_vel / vel);
            if(x_vel < 0.0) theta += M_PI;

          const static double sigma = 1 / 1e3;
          double en_scale = sigma * sqrt(2.0*5e7 / constants::m_e_r)/sqrt(2.0 * M_PI);
          vel += en_scale* exp(-0.5*sigma*sigma*double(current_time_step - 40000)*double(current_time_step - 40000));
         //   std::cout << current_time_step / (sim::equilibration_time+40000.0) << std::endl;

        x_vel = vel * cos(theta)*sin(phi);
        y_vel = vel * sin(theta)*sin(phi);
        z_vel = vel * cos(phi);
        }
    }
        
        new_electron_velocity[array_index]   = x_vel;
        new_electron_velocity[array_index_y] = y_vel;
        new_electron_velocity[array_index_z] = z_vel;

        EKE = vel*vel;
        
     //   x_vel = new_atom_velocity[array_index]   = atom_velocity[array_index]   + ((atom_force[array_index]   + new_atom_force[array_index])   * dt  * constants::K_A / 2); 
       // y_vel = new_atom_velocity[array_index_y] = atom_velocity[array_index_y] + ((atom_force[array_index_y] + new_atom_force[array_index_y]) * dt  * constants::K_A / 2);
        //z_vel = new_atom_velocity[array_index_z] = atom_velocity[array_index_z] + ((atom_force[array_index_z] + new_atom_force[array_index_z]) * dt  * constants::K_A / 2);
        //vel = (x_vel*x_vel) + (y_vel*y_vel) + (z_vel*z_vel);
       // LKE += vel;

}

void e_a_coulomb(const int e, const int& array_index, double& e_x_force, double& e_y_force, double& e_z_force, double& EPE){
                      //  double& a_x_force, double& a_y_force, double& a_z_force, double& EPE, double& LPE) {

    double x_distance;
    double y_distance;
    double z_distance;
   
    double length;
    double  force,  phi,theta, PE = 0;
    int array_index_a, nearest_electron_count = 1;
    int nearest_atom_count = 1;
   // bool collision = false;
    // std::srand(std::time(nullptr));
    // std::random_device rd;  //Will be used to obtain a seed for the random number engine
    // std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    // std::uniform_real_distribution<double> scattering_chance(0,1);
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

    //      if(a == 100 ) std::cout << array_index_e / 3 << ", " << atom_position[array_index] << ", " << atom_position[array_index+1] << ", " << atom_position[array_index+2] << ", " << sqrt(length) << std::endl;
    //     if(a == 100 ) std::cout << array_index_e / 3 << ", " << new_electron_position[array_index_e] << ", " << new_electron_position[array_index_e+1] << ", " << new_electron_position[array_index_e+2] << ", " << sqrt(length) << std::endl;
    //     if(a == 100 ) std::cout << array_index_e / 3 << ", " << atom_position[array_index] - new_electron_position[array_index_e] << ", " << atom_position[array_index+1] - new_electron_position[array_index_e+1]<< ", " << atom_position[array_index+2] - new_electron_position[array_index_e+2] << std::endl;
    //     // if(a == 100 ) std::cout << array_index_e / 3 << ", " << (new_atom_position[array_index] - new_electron_position[array_index_e])*(new_atom_position[array_index] - new_electron_position[array_index_e]) << ", " << (new_atom_position[array_index+1] - new_electron_position[array_index_e+1])*(new_atom_position[array_index+1] - new_electron_position[array_index_e+1]) << ", " << (new_atom_position[array_index+2] - new_electron_position[array_index_e+2])*(new_atom_position[array_index+2] - new_electron_position[array_index_e+2]) << ", " << sqrt(length) << std::endl;
    //    // if(a == 100 ) std::cout << array_index_e / 3 << ", " << (new_atom_position[array_index] - new_electron_position[array_index_e])*(new_atom_position[array_index] - new_electron_position[array_index_e]) + (new_atom_position[array_index+1] - new_electron_position[array_index_e+1])*(new_atom_position[array_index+1] - new_electron_position[array_index_e+1]) + (new_atom_position[array_index+2] - new_electron_position[array_index_e+2])*(new_atom_position[array_index+2] - new_electron_position[array_index_e+2]) << ", " << sqrt(length) << std::endl;
    //     //if(a == 100 ) std::cout << array_index_e / 3 << ", " << sqrt(((new_atom_position[array_index] - new_electron_position[array_index_e])*(new_atom_position[array_index] - new_electron_position[array_index_e])) + ((new_atom_position[array_index+1] - new_electron_position[array_index_e+1])*(new_atom_position[array_index+1] - new_electron_position[array_index_e+1]))+ ((new_atom_position[array_index+2] - new_electron_position[array_index_e+2])*(new_atom_position[array_index+2] - new_electron_position[array_index_e+2]))) << ", " << sqrt(length) << std::endl;
     //    if(a==100) std::cout << array_index_e / 3 << ", " << x_distance << ", " << y_distance << ", " << z_distance << "\n" <<  std::endl;
        atomic_nearest_electron_list[e][nearest_electron_count] = array_index_a;
        nearest_electron_count++;

        if (length > e_a_coulomb_cutoff) continue;
       // count++;

       // if(a == 100) std::cout << e << ", " << sqrt(length) << std::endl;
        length = sqrt(length);
        if(length < 0.1) length = 0.1;
       // if(a==100) std::cout << e << ", " << length << std::endl;
        // if(length < 0.11) length = 0.11;

        if(mean_radius[e] > length) {
            #pragma omp critical
            mean_radius[e] = length;
        }

        force = -1*(1/(length * length) - 8*150*exp(-8*length));
                        //q*k*k * exp(-15(A**-1) * length (A));
          //  std::cout << force << std::endl;
        PE += 150*exp(-8*length) - (1 / length);
        

        phi   = acos(z_distance / length);
        theta = atan(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        e_x_force += mu_r * force * cos(theta)*sin(phi); 
        e_y_force += mu_r * force * sin(theta)*sin(phi);
        e_z_force += mu_r * force * cos(phi);

    /* //   a_x_force += x_force * combined_mass;
       // a_y_force += y_force * combined_mass;
        //a_z_force += z_force * combined_mass;

        // e_x_force -= x_force * mu_r;
        // e_y_force -= y_force * mu_r;
        // e_z_force -= z_force * mu_r;
        a_x_force += force * cos(theta)*sin(phi);
        a_y_force += force * sin(theta)*sin(phi);
        a_z_force += force * cos(phi);

        e_x_force += -1*force * cos(theta)*sin(phi);
        e_y_force += -1*force * sin(theta)*sin(phi);
        e_z_force += -1*force * cos(phi); */
        
    /*    if(a == phonon_scattering_vector(gen)) {
            
            // if(electron_nearest_atom_list[e][nearest_electron_count*2] == a && electron_nearest_atom_list[e][nearest_electron_count*2 + 1]) continue;
            // else if (electron_nearest_atom_list[e][nearest_electron_count*2] != a) {
            //     electron_nearest_atom_list[e][nearest_electron_count*2] = a;
            //     electron_nearest_atom_list[e][nearest_electron_count*2 + 1] = false;
            // }
            //double excitation_constant = -1.0 * dt / E_f_A;
            //if(current_time_step > 4000) excitation_constant = 0.1;
                double scattering_velocity = (electron_velocity[array_index]*electron_velocity[array_index]) + (electron_velocity[array_index+1]*electron_velocity[array_index+1]) + (electron_velocity[array_index+2]*electron_velocity[array_index+2]);
                double excitation_constant = electron_potential[e]*constants::K_A + constants::m_e_r*0.5*scattering_velocity - atom_potential[a];
                double excitation_energy   = excitation_constant;
                if(abs(excitation_energy)) {
                   //  std::cout << excitation_energy << std::endl;
                    if(abs(excitation_constant/(mu_f)) > 1) excitation_energy = mu_f; 
                    scattering_velocity = sqrt(scattering_velocity) - sqrt(2*excitation_energy/constants::m_e_r);
                    if(excitation_constant < 0) {
                        scattering_velocity += 2*sqrt(excitation_energy*2/constants::m_e_r);
                    }
                    
                    if(scattering_velocity > 0) {
                    if(scattering_chance(gen) > exp(-1*dt*abs(excitation_constant) / E_f_A)) {
                    
                    // electron_nearest_atom_list[e][2*a+1] = true;
                        collision = true;
                        double vel = sqrt((electron_velocity[array_index]*electron_velocity[array_index]) + (electron_velocity[array_index+1]*electron_velocity[array_index+1]) + (electron_velocity[array_index+2]*electron_velocity[array_index+2]));
                        double theta = atan(electron_velocity[array_index+1] / electron_velocity[array_index]);
                        double phi = acos(electron_velocity[array_index+2] / vel);
                        if (electron_velocity[array_index] < 0) theta += M_PI;
                        electron_velocity[array_index]   = scattering_velocity * cos(theta)*sin(phi);
                        electron_velocity[array_index+1] = scattering_velocity * sin(theta)*sin(phi);
                        electron_velocity[array_index+2] = scattering_velocity * cos(phi);
                       // std::cout << scattering_velocity << std::endl;
                        #pragma omp critical
                        {
                        chosen_electron++;
                        new_atom_potential[a] += excitation_energy;
                        }
                    }
                }
            }
        } */
    }
  //  LPE += PE /2;
    EPE += PE;
    atomic_nearest_electron_list[e][0] = nearest_electron_count;
    new_electron_potential[e] = PE;
   // if(a == 100) std::cout << count << std::endl;
    //if(a == 100) std::cout << atomic_nearest_electron_list[a][0] << std::endl;
}

void neighbor_e_a_coulomb(const int e, const int& array_index, double& e_x_force, double& e_y_force, double& e_z_force, double& EPE){
                      //  double& a_x_force, double& a_y_force, double& a_z_force, double& EPE, double& LPE) {

    double x_distance;
    double y_distance;
    double z_distance;
    
    double length;
    double force,  phi,theta, PE = 0;
    int array_index_a;
  //  bool collision = false;
    int size = atomic_nearest_electron_list[e][0];
    int nearest_atom_count = 0;
std::srand(std::time(nullptr));
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_real_distribution<double> scattering_chance(0,1);
            std::uniform_int_distribution<> phonon_scattering_vector(1,27);
    int phonon_collision = phonon_scattering_vector(gen);// atomic_nearest_electron_list[e][phonon_scattering_vector(gen)];
 //   if(a == 100) std::cout << atomic_nearest_electron_list[a][0] << std::endl;

     int count = 0;
    for (int a = 1; a < size; a++) {
       // if(atomic_nearest_electron_list[a][e] < 0) std::cout << a << ", " << e << std::endl;
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
    //       if(a == 100 ) std::cout << array_index_e / 3 << ", " << atom_position[array_index] << ", " << atom_position[array_index+1] << ", " << atom_position[array_index+2] << ", " << sqrt(length) << std::endl;
    //     if(a == 100 ) std::cout << array_index_e / 3 << ", " << new_electron_position[array_index_e] << ", " << new_electron_position[array_index_e+1] << ", " << new_electron_position[array_index_e+2] << ", " << sqrt(length) << std::endl;
      //    if(a == 100 ) std::cout << array_index_e / 3 << ", " << atom_position[array_index] - new_electron_position[array_index_e] << ", " << atom_position[array_index+1] - new_electron_position[array_index_e+1]<< ", " << atom_position[array_index+2] - new_electron_position[array_index_e+2] << ", " << sqrt(length) << std::endl;
    //     //if(a == 100 ) std::cout << array_index_e / 3 << ", " << (new_atom_position[array_index] - new_electron_position[array_index_e])*(new_atom_position[array_index] - new_electron_position[array_index_e]) << ", " << (new_atom_position[array_index+1] - new_electron_position[array_index_e+1])*(new_atom_position[array_index+1] - new_electron_position[array_index_e+1]) << ", " << (new_atom_position[array_index+2] - new_electron_position[array_index_e+2])*(new_atom_position[array_index+2] - new_electron_position[array_index_e+2]) << ", " << sqrt(length) << std::endl;
    //    // if(a == 100 ) std::cout << array_index_e / 3 << ", " << (new_atom_position[array_index] - new_electron_position[array_index_e])*(new_atom_position[array_index] - new_electron_position[array_index_e]) + (new_atom_position[array_index+1] - new_electron_position[array_index_e+1])*(new_atom_position[array_index+1] - new_electron_position[array_index_e+1]) + (new_atom_position[array_index+2] - new_electron_position[array_index_e+2])*(new_atom_position[array_index+2] - new_electron_position[array_index_e+2]) << ", " << sqrt(length) << std::endl;
    //    // if(a == 100 ) std::cout << array_index_e / 3 << ", " << sqrt(((new_atom_position[array_index] - new_electron_position[array_index_e])*(new_atom_position[array_index] - new_electron_position[array_index_e])) + ((new_atom_position[array_index+1] - new_electron_position[array_index_e+1])*(new_atom_position[array_index+1] - new_electron_position[array_index_e+1]))+ ((new_atom_position[array_index+2] - new_electron_position[array_index_e+2])*(new_atom_position[array_index+2] - new_electron_position[array_index_e+2]))) << ", " << sqrt(length) << std::endl;
     //    if(a==100) std::cout << array_index_e / 3 << ", " << x_distance << ", " << y_distance << ", " << z_distance << "\n" << std::endl;
        if (length > e_a_coulomb_cutoff) continue;
        // count++;
        
        length = sqrt(length);
        if(length < 0.1) length = 0.1;
       // if(a==100) std::cout << array_index_e / 3 << ", " << length << std::endl;
       if(mean_radius[e] > length) {
            #pragma omp critical
            mean_radius[e] = length;
        }
        force = -1*(1/(length * length) - 8*150*exp(-8*length));
                        //q*k*k * exp(-15(A**-1) * length (A));
          //  std::cout << force << std::endl;
            PE += 150*exp(-8*length) - (1 / length);
        
        phi   = acos(z_distance / length);
        theta = atan(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        e_x_force += mu_r * force * cos(theta)*sin(phi); 
        e_y_force += mu_r * force * sin(theta)*sin(phi);
        e_z_force += mu_r * force * cos(phi);

 /*  //     a_x_force += x_force * combined_mass;
     //   a_y_force += y_force * combined_mass;
       // a_z_force += z_force * combined_mass;

        // e_x_force -= x_force * mu_r;
        // e_y_force -= y_force * mu_r;
        // e_z_force -= z_force * mu_r;
        a_x_force += force * cos(theta)*sin(phi);
        a_y_force += force * sin(theta)*sin(phi);
        a_z_force += force * cos(phi);

        e_x_force += -1*force * cos(theta)*sin(phi);
        e_y_force += -1*force * sin(theta)*sin(phi);
        e_z_force += -1*force * cos(phi); */
        if(length < 3.0) {
            count++;
            if(count == phonon_collision) {
            
            // if(electron_nearest_atom_list[e][nearest_electron_count*2] == a && electron_nearest_atom_list[e][nearest_electron_count*2 + 1]) continue;
            // else if (electron_nearest_atom_list[e][nearest_electron_count*2] != a) {
            //     electron_nearest_atom_list[e][nearest_electron_count*2] = a;
            //     electron_nearest_atom_list[e][nearest_electron_count*2 + 1] = false;
            // }
            //double excitation_constant = -1.0 * dt / E_f_A;
            //if(current_time_step > 4000) excitation_constant = 0.1;
               double scattering_velocity = 0.5*constants::m_e_r*((electron_velocity[array_index]*electron_velocity[array_index]) + (electron_velocity[array_index+1]*electron_velocity[array_index+1]) + (electron_velocity[array_index+2]*electron_velocity[array_index+2]));
               // double EE = electron_potential[e]*constants::K_A + constants::m_e_r*0.5*scattering_velocity;
                // - atom_potential[array_index_a/3];
               // double excitation_energy   = deltaE;
                  //  excitation_energy = mu_f + mu_f * floor(abs(excitation_constant / E_f_A))*floor(abs(excitation_constant / E_f_A)); 
                //if(deltaE > 0) scattering_velocity = sqrt(scattering_velocity) - (sqrt(2.0*deltaE/constants::m_e_r));
                //else scattering_velocity = sqrt(scattering_velocity) + (sqrt(2.0*abs(deltaE)/constants::m_e_r));  
                //if(scattering_velocity < 0) continue;
                //scattering_velocity = constants::m_e_r*0.5*scattering_velocity;
                if(scattering_chance(gen) > exp(-1.0*dt*sqrt(scattering_velocity / atom_potential[array_index_a/3]) / 27.7)) {
                   // std::cout << -1.0*dt*sqrt(atom_potential[array_index_a/3] / abs(electron_potential[e]*constants::K_A + constants::m_e_r*0.5*scattering_velocity)) / 27.7 << std::endl;
                    double deltaE = scattering_velocity - atom_potential[array_index_a/3];
                    if (deltaE > E_f_A) deltaE = E_f_A;
                    else if(deltaE < 0.0) deltaE = fmax(E_f_A - atom_potential[array_index_a/3], -1.0 * E_f_A);
                    
                    std::uniform_real_distribution<double> Theta_pos_distrib(0.0,2.0*M_PI);
                    std::uniform_real_distribution<double> Phi_pos_distrib(0.0,M_PI);
                    double theta = Theta_pos_distrib(gen);//*M_PI;
                    double phi   = Phi_pos_distrib(gen);// * M_PI;
                    
                    scattering_velocity = sqrt(2.0*(scattering_velocity - deltaE)/constants::m_e_r);
                    
                    electron_velocity[array_index]   = scattering_velocity * cos(theta)*sin(phi);
                    electron_velocity[array_index+1] = scattering_velocity * sin(theta)*sin(phi);
                    electron_velocity[array_index+2] = scattering_velocity * cos(phi);
                       // std::cout << scattering_velocity << std::endl;
                    #pragma omp critical
                    {
                    chosen_electron++;
                    atom_potential[array_index_a/3] += deltaE;
                    }
                }
            }
        }
    }
    EPE += PE;
    new_electron_potential[e] = PE;
  //  LPE += PE/2;
   // if(a == 100) std::cout << count << std::endl;
}

void e_e_coulomb(const int e, const int array_index, double& e_x_force, double& e_y_force, double& e_z_force, \
                double& EPE) {
    
    int array_index_i;
    double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
    int neighbor_count = 1;

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
            if(length < 0.1) length = 0.1;
       // if(length < 0.11) length = 0.11;
        force = 1 / (length*length);
        //     if(e == 100 ) std::cout << array_index_i / 3 << ", " << new_electron_position[array_index] - new_electron_position[array_index_i] << ", " << new_electron_position[array_index+1] - new_electron_position[array_index_i+1]<< ", " << new_electron_position[array_index+2] - new_electron_position[array_index_i+2] << ", " << sqrt(length) << std::endl;
        //    if(e==100) std::cout << array_index_i / 3 << ", " << x_distance << ", " << y_distance << ", " << z_distance << "\n" << std::endl;
        PE += force * length;
            
            phi   = acos(z_distance / length);
            theta = atan(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            e_x_force += force * cos(theta)*sin(phi) / constants::m_e_r;
            e_y_force += force * sin(theta)*sin(phi)/ constants::m_e_r;
            e_z_force += force * cos(phi) / constants::m_e_r;

        }
    electron_nearest_electron_list[e][0] = neighbor_count;
    EPE += PE/2;
    new_electron_potential[e] += PE;
   // if(e == 100) std::cout << neighbor_count << std::endl;
}

void neighbor_e_e_coulomb(const int e, const int array_index, double& e_x_force, double& e_y_force, double& e_z_force, \
                                double& EPE) {
    
    double x_distance,y_distance,z_distance, length, force, theta,phi, PE = 0;
    int size = electron_nearest_electron_list[e][0]; //.size();
    int array_index_i;

    //std::cout << size << std::endl;

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
        if(length < 0.1) length = 0.1;
       // if(length < 0.11) length = 0.11;
        force = 1 / (length*length);

        // if(e == 100 ) std::cout << array_index_i / 3 << ", " << new_electron_position[array_index] - new_electron_position[array_index_i] << ", " << new_electron_position[array_index+1] - new_electron_position[array_index_i+1]<< ", " << new_electron_position[array_index+2] - new_electron_position[array_index_i+2] << ", " << sqrt(length) << std::endl;
        // if(e==100) std::cout << array_index_i / 3 << ", " << x_distance << ", " << y_distance << ", " << z_distance << "\n" << std::endl;
       
        
        PE += force * length;

        phi   = acos(z_distance / length);
        theta = atan(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        e_x_force += force * cos(theta)*sin(phi) / constants::m_e_r;
            e_y_force += force * sin(theta)*sin(phi)/ constants::m_e_r;
            e_z_force += force * cos(phi) / constants::m_e_r;

    }
    EPE += PE/2;
    new_electron_potential[e] += PE;
  //  if(e == 100) std::cout << size << std::endl;
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

double electron_applied_voltage(int array_index, double& x_force, double& y_force, double& z_force) {
    
  //  x_force -= 4e-9;
  //  new_electron_potential[array_index/3] += -4e-9;
    return 0;//-4e-9; //-10.0 * 1e-10 / constants::e ;
}


} //end CASTLE namespace

