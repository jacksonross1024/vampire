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
    
    
    TPE = 0;
    TKE = 0;
    TLE = 0;

            if (err::check) std::cout << "Calculating verlet integration...Time step: " << dt << std::endl;

   // std::ofstream electron_position_output_up;
    //std::ofstream electron_position_output_down;
    //std::ofstream electron_velocity_output;
    //std::ofstream electron_spin_output;

            if(err::check) std::cout << "Initializing output files..." << std::endl;
    if (current_time_step % CASTLE_output_rate == 0) setup_output();

            if (err::check) std::cout << "Updating new electron position." << std::endl;
    update_position();

            if (err::check) std::cout << "Forces, spins, and velocities update..." << std::endl;
    update_dynamics();
   // MLE += reinitialize_electron_conserve_momentum(captured_electron_list);
  //  std::cout << MPE/CASTLE_output_rate << ", " << TPE << ", " << MKE/CASTLE_output_rate << ", " << TKE << ", " << (MPE+MKE)/CASTLE_output_rate << ", " << TPE +TKE << ", " << std::endl;
            if (err::check) std::cout << "Output mean data" << std::endl;
    if (CASTLE_output_data)   output_data(); //std::cout << "x_flux: " << x_flux / CASTLE_output_rate << "\n"; x_flux = 0;
    


    //reset integration
    CASTLE_real_time += dt;
    current_time_step++;

    electron_position.swap(new_electron_position);
    electron_force.swap(new_force_array);
    electron_velocity.swap(new_electron_velocity);
    electron_potential.swap(new_electron_potential);

   // std::fill(electron_potential.begin(), electron_potential.end(), 0);
  //  std::fill(electron_capture.begin(), electron_capture.end(), false);
   
    return EXIT_SUCCESS;
 
        
      /*  for (int i = 0; i < lattice_electrons.size(); i++) {
            array_index_i = 3*i;
       //     electron_spin_two = lattice_electron_spin[i];
            x_distance = x - lattice_electrons[array_index_i];
            y_distance = y - lattice_electrons[array_index_i_y];
            z_distance = z - lattice_electrons[array_index_i_z];

            length = sqrt((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance));
            if (length < 0.00001) length = 0.00001;
            
            if (electron_spin == electron_spin_two) {
                force = 1e10 * ((constants::K / (length * length * constants::m_e))- (635 * constants::kB / (length * constants::m_e)));
                
            }
            else { 
                force = 1/ (length * length);

            
            TPE += -2 * force * length;

            new_force_array[array_index]   += x_distance * force / length; //Angstroms
            new_force_array[array_index_y] += y_distance * force / length;
            new_force_array[array_index_z] += z_distance * force / length;
        */
}

void setup_output() {

    CASTLE_output_data = true;
    time_stamp = std::to_string(current_time_step);

      //  electron_position_output_up.open("CASTLE/Electron_Position_Up/" + time_stamp + "U.xyz");
        //    electron_position_output_up << total_spin_up << "\n";
          //  electron_position_output_up << time_stamp << "\n";

    electron_position_output_down.open("CASTLE/Electron_Position/" + time_stamp + ".xyz");
    electron_position_output_down << conduction_electrons << "\n";
    electron_position_output_down << time_stamp << "\n";

    electron_velocity_output.open("CASTLE/Electron_Velocity/" + time_stamp + ".csv");
    electron_velocity_output << "Electron number,    x-component,     y-component,    z-component,     length" << "\n";

    
      //  electron_spin_output.open("CASTLE/Electron_Spin/" + time_stamp + ".txt");
       //     electron_spin_output << conduction_electrons << "\n";
       //     electron_spin_output << "       1 is up, 0 is down" << "\n";
        // acceleration_output.open("CASTLE/test" + time_stamp);
 
}


void update_position(){

    int array_index,array_index_y,array_index_z = 0;
    long double x_pos,y_pos,z_pos = 0;
    #pragma omp parallel for private(array_index,array_index_y,array_index_z, x_pos,y_pos,z_pos) schedule(static) reduction(+:x_flux,y_flux,z_flux)
    for (int e = 0; e < conduction_electrons; e++){ 

      //  TLE += atomic_phonon_energy[2*e];

        array_index = 3*e;
        array_index_y = array_index + 1;
        array_index_z = array_index + 2;

      // std::cout << "why are you dying...." << e << "\n";
    
        x_pos = electron_position[array_index]   + (electron_velocity[array_index]   * dt) + (electron_force[array_index]   * dt * dt * constants::K_A / (2*constants::m_e_r)); // x superarray component
        y_pos = electron_position[array_index_y] + (electron_velocity[array_index_y] * dt) + (electron_force[array_index_y] * dt * dt * constants::K_A / (2*constants::m_e_r)); // y superarray component
        z_pos = electron_position[array_index_z] + (electron_velocity[array_index_z] * dt) + (electron_force[array_index_z] * dt * dt * constants::K_A / (2*constants::m_e_r)); // z superarray component
        
      //  if(e==100) std::cout << x_pos << ", " << electron_position[array_index] << ", " << (electron_velocity[array_index]   * dt) << ", " << (electron_force[array_index]   * dt * dt * constants::K_A / (2*constants::m_e_r)) << std::endl; // x superarray component
     //  if(e==100) std::cout << e << ", " << electron_velocity[array_index] << " , " << electron_velocity[array_index + 1] << " , " << electron_velocity[array_index + 2] << std::endl;
        if (x_pos < 0) {
            x_pos += 40;
            x_flux--;
        }
        else if (x_pos > 40) {
            x_pos -= 40;
            x_flux++;
        }
	  /*  if(abs(x_pos) > 40.0) {
		    terminaltextcolor(RED);
		    std::cout << "Latice bounds exceeded for " << e << " on " << current_time_step << ". " << x_pos << std::endl;
        } */

	    if (y_pos < 0) {
            y_pos += 40;
            y_flux--;
        }
        else if (y_pos > 40) {
            y_pos -= 40;
            y_flux++;
        }
	/*    if(abs(y_pos) > 40.0) {
		    terminaltextcolor(RED);
		    std::cout << "Latice bounds exceeded for " << e << " on " << current_time_step << ". " << y_pos << std::endl;
        } */

	    if (z_pos < 0) {
            z_pos += 40;
            z_flux--;
        }
        else if (z_pos > 40) {
            z_pos -= 40;
            z_flux++;
        }
	 /*   if(abs(z_pos) > 40.0) {
		    terminaltextcolor(RED);
		    std::cout << "Latice bounds exceeded for " << e << " on " << current_time_step << ". " << z_pos << std::endl;
	    } */

        new_electron_position[array_index]   = x_pos;
        new_electron_position[array_index_y] = y_pos;
        new_electron_position[array_index_z] = z_pos;

       // if (e == 0) std::cout << new_electron_position[array_index] << "   " << electron_position[array_index]  << "   " << (electron_velocity[array_index] * dt) << "   " << (electron_force[array_index] * dt * dt * 0.5  * 1e30 * constants::K / constants::m_e) << "\n";
        //symmetry_list[e].resize(conduction_electrons, false);
    } 

    int b;
    int p_p_coupling;
    int size;
    long double d_x,d_y,d_z, excitation_energy;
    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
   std::uniform_int_distribution<> phonon_random_walk(1,18);
        std::uniform_real_distribution<long double> phonon_jump_distrib(0,1);
   // #pragma omp parallel for private(b,array_index,excitation_energy, size) schedule(dynamic) reduction(+:TLE)
    for(int a = 0; a < lattice_atoms; a++) {
        
        b = 0;
        TLE += atomic_phonon_energy[a];
       // std::cout << atomic_phonon_energy[2*a] << std::endl;
        size = atomic_nearest_particle_list[a][0];
      //  std::cout << size << std::endl;
     //   x = atom_position[3*a];
       // y = atom_position[3*a + 1];
        //z = atom_position[3*a + 2]; 
        /*
        while(b != size) {
            int choice = phonon_random_walk(gen);
            
            array_index = atomic_nearest_particle_list[a][choice];
            excitation_energy = 0.25*(atomic_phonon_energy[a] - atomic_phonon_energy[array_index]);
           // d_x = x - atom_position[array_index];
            //d_y = y - atom_position[array_index+1];
            //d_z = z - atom_position[array_index+2]; 

          //  length = sqrtl((d_x*d_x) + (d_y*d_y) + (d_z*d_z));
            if(excitation_energy > 0) {
                if(phonon_jump_distrib(gen) < (0.2 - 0.2*exp(-1*excitation_energy))) { //*atomic_phonon_energy[a*2+1])) {
                   // #pragma omp critical
                //    std::cout << atomic_phonon_energy[a] << ", " << atomic_phonon_energy[array_index] << ", " << array_index << ", " << choice << std::endl;
                    atomic_phonon_energy[array_index] += excitation_energy;
                    atomic_phonon_energy[a] -= excitation_energy;
                    
                  //  atomic_phonon_energy[2*a + 1] = 0;
                  //  TLE += excitation_energy;
                }
            }
            b++; */
        
    }// std::cout << TLE;

}

long double update_dynamics() {
   // double applied_electronic_field = {0.0, 0.0, 1.0, 1.0}; //x, y, z, strength
    int array_index;
    long double x_force,y_force,z_force, x,y,z;
   
 /*   std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> random_electron(0,8000);
    chosen_electron = random_electron(gen); */
    #pragma omp parallel for private(array_index, x_force,y_force,z_force, x,y,z)  schedule(dynamic) reduction(+:TPE,TKE)
    for (int e = 0; e < conduction_electrons; e++) {
        array_index = 3*e;
        x_force = 0;
        y_force = 0;
        z_force = 0;
        x = new_electron_position[array_index];
        y = new_electron_position[array_index + 1];
        z = new_electron_position[array_index + 2];

        if(current_time_step % 10 == 0) {
            TPE += electron_e_a_coulomb(e, array_index, x_force,y_force,z_force, x,y,z);
        //if (e == 1050) std::cout << x << "  " << y << " " << z <<   "    " << std::endl;
        } else {
            TPE += neighbor_e_a_coulomb(e, array_index, x_force,y_force,z_force, x,y,z);
        } 
        if(current_time_step % 10 == 0) {
            TPE += electron_e_e_coulomb(e, array_index, x_force,y_force,z_force, x,y,z);
           // std::cout << nearest_neighbor_list[0].size() << std::endl;
        } else {
            TPE += neighbor_e_e_coulomb(e, array_index, x_force,y_force,z_force, x,y,z);
        }
      //  TLE += e_e_scattering(scattering_list);
        if(!equilibrium_step) TPE += electron_applied_voltage(array_index, x_force,y_force,z_force);
       // if (e ==1050) std::cout << x_force << "    " << y_force << "   " << z_force << " normalization: " << constants::K / constants::m_e <<  std::endl;

       // #pragma omp critical 
        
        new_force_array[array_index]     = x_force;
        new_force_array[array_index + 1] = y_force;
        new_force_array[array_index + 2] = z_force;
        
      //  if (e ==1050) std::cout << new_force_array[array_index] << "    " << new_force_array[array_index+1] << "   " << new_force_array[array_index+2] << " normalization: " << constants::K / constants::m_e <<  std::endl;
        TKE += update_velocity(array_index);
        
    }
    MPE += TPE;
    MKE += TKE;
    MLE += TLE;
    //std::cout << ", " << TLE << std::endl;
   // std::cout << "Count " << current_time_step << ", TPE: " << TPE * 1e10 * constants::K  << ", TKE: " << TKE * 1e10 * constants::m_e  << ", mean-MPE: " << MPE * 1e10 * constants::K / CASTLE_output_rate << ", mean-MKE: " << MKE * 1e10 * constants::m_e / CASTLE_output_rate << std::endl;

        //spontaneous spin flip

       
        
    //    double spin_chance = Spin_distrib(gen) * 0.0002;
     //   symmetry_list[e].resize(conduction_electrons, false);
        
      // if (e == 0) std::cout << "deltaV " << deltaV << " eps " << (0.5 * constants::m_e * velocity_length * velocity_length) / (constants::kB * temperature) << " E_f " <<  E_f / (constants::kB * temperature) << "\n";
 /*       double flip_chance = 1.0;
        if (((0.5 * constants::m_e * velocity_length * velocity_length) - E_f) <  0) flip_chance = exp(((0.5 * constants::m_e * velocity_length * velocity_length) - E_f) / (constants::kB * temperature));
       // std::cout << "flip chance " << flip_chance << " spin chance " << spin_chance << "\n";
        if (spin_chance > flip_chance) {        
            double deltaV = sqrt(2 * mu_f / constants::m_e) * 1e10;
            bool new_spin = conduction_electron_spin[e] = !conduction_electron_spin[e];
            double scaling_factor, new_velocity_length = 0;
            if (new_spin) {
                total_spin_up += 1;
                total_spin_down -= 1;
                if (total_spin_up > total_spin_down ) new_velocity_length = (1e10 * velocity_length) - deltaV;
                else new_velocity_length = (1e10 * velocity_length) + deltaV;          
            } else {
                total_spin_down += 1;
                total_spin_up -= 1;
                if (total_spin_up > total_spin_down) new_velocity_length = (1e10 * velocity_length) + deltaV;
                else new_velocity_length = (1e10 * velocity_length) - deltaV;
            }
            scaling_factor = new_velocity_length / (1e10 * velocity_length);
            electron_velocity[array_index] *= scaling_factor;
            electron_velocity[array_index + 1] *= scaling_factor;
            electron_velocity[array_index + 2] *= scaling_factor;
            velocity_length = 1e-10 * new_velocity_length;
        } */
       // TKE += 0.5 * constants::m_e * velocity_length * velocity_length * 1e20; //energy in Angstroms

          //  if (err::check) std::cout << "Forces advancing..." << "\n";
     //   electron_spin = conduction_electron_spin[e];
      //  if (CASTLE_output_data) electron_spin_output << e << "  " << electron_spin << "\n";
    return 0;///TLE;
}

long double update_velocity(int array_index) {
        int array_index_y = array_index + 1;
        int array_index_z = array_index + 2;

     //   if (e == 0) std::cout << new_electron_velocity[array_index] << " " << electron_velocity[array_index] << "    " <<  electron_force[array_index]  << "    " << new_force_array[array_index]  << "    " <<  dt * 0.5  * constants::K / constants::m_e << "\n"; 
        long double x_vel = new_electron_velocity[array_index]   = electron_velocity[array_index]   + ((electron_force[array_index]   + new_force_array[array_index])   * dt  * constants::K_A / (2*constants::m_e_r)); 
        long double y_vel = new_electron_velocity[array_index_y] = electron_velocity[array_index_y] + ((electron_force[array_index_y] + new_force_array[array_index_y]) * dt  * constants::K_A / (2*constants::m_e_r));
        long double z_vel = new_electron_velocity[array_index_z] = electron_velocity[array_index_z] + ((electron_force[array_index_z] + new_force_array[array_index_z]) * dt  * constants::K_A / (2*constants::m_e_r));
        
    //    if(array_index/3 ==100) std::cout << x_vel << ", " << electron_velocity[array_index] << ", " << y_vel << ", " << electron_velocity[array_index_y] << ", " << z_vel << ", " << electron_velocity[array_index_z] << std::endl;
       // new_electron_velocity[array_index]   = x_vel;
       // new_electron_velocity[array_index_y] = y_vel;
       // new_electron_velocity[array_index_z] = z_vel;
        /*
        velocity_length_hist[array_index]   += abs(x_vel);
        velocity_length_hist[array_index_y] += abs(y_vel);
        velocity_length_hist[array_index_z] += abs(z_vel); 
       */// long double velocity_length = (x_vel*x_vel) + (y_vel*y_vel) + (z_vel*z_vel);
        //velocity_length_hist[array_index + 3] += sqrt(velocity_length); 
        //if(array_index / 3 ==0) std::cout << velocity_length_hist[array_index/3] << "\n";
    return ((x_vel*x_vel) + (y_vel*y_vel) + (z_vel*z_vel));
}

long double electron_e_a_coulomb(int e, int array_index, long double& x_force, long double& y_force, long double& z_force, const long double& x, const long double& y, const long double& z) {
  //set e-a attraction
        //calculate nearest neighbor;
    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<long double> range_distrib(-1,1);
    std::uniform_real_distribution<long double> scattering_prob_distrib(0,1);
    std::normal_distribution<long double> velocity_gaussian_distrib(0.1,0.1);

    long double  x_distance,y_distance,z_distance,d_x,d_y,d_z, l_x,l_y,l_z, x_mod,y_mod,z_mod, length, force;
    long double phi,theta, PE = 0;
    int atom_array, neighbor_count = 1;
    bool collision = false;

    long double excitation_constant;// = velocity_gaussian_distrib(gen);
    long double Px;// = electron_velocity[array_index];
    long double Py;// = electron_velocity[array_index+1];
    long double Pz;// = electron_velocity[array_index+2];
    long double P;// = sqrtl(Px*Px+Py*Py+Pz*Pz);
    long double P_p;// = sqrtl(Px*Px+Py*Py);
    long double P_i;// = sqrtl(Pz*Pz);
    long double scattering_velocity;// = P;
    long double electron_TE;// = 1e10*(electron_potential[e]*constants::K + P*P*constants::m_e*0.5);
   
   

    long double d_p;// = sqrtl(l_x*l_x+l_y*l_y);
    long double d_i;// = sqrtl(l_z*l_z);
    long double d_r;// = sqrtl((l_x*l_x)+(l_y*l_y)+(l_z*l_z));
    long double polar_value;// = M_PI * range_distrib(gen);
    long double incline_value;// = M_PI * range_distrib(gen);
    long double normal_polar_angle;// = acosl((l_x*Px+l_y*Py)/(P_p*d_p));
    long double normal_incline_angle ;//= acosl((l_z*Pz) / (P_i*d_i));
    long double polar_scattering_angle;
    long double excitation_energy; //Angstroms
    // = abs(excitation_constant)*(electron_TE - atomic_phonon_energy[2*a]);
    long double incline_scattering_angle;
    long double polar_prob;
    long double incline_prob;
    /*std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<long double> capture_chance_distrib(0,1); */

  //  d_x = x - ((atomic_size * round((x-1) / atomic_size)) + 1); //closest x atom index
    //d_y = y - ((atomic_size * round((y-1) / atomic_size)) + 1); //closest y atom index
    //d_z = z - ((atomic_size * round((z-1) / atomic_size)) + 1); //closest z atom index

    for (int a = 0; a < lattice_atoms; a++) {
        atom_array = 3*a;
        x_distance = atom_position[atom_array]   - x;
        y_distance = atom_position[atom_array+1] - y;
        z_distance = atom_position[atom_array+2] - z; 

    //mean_radius[(array_index/3)*2] = radius;
    //mean_radius[(array_index/3)*2 + 1] =  28*((3.3*expl(-3.3*radius)) - (expl(-1*radius)));
    //if (array_index / 3 == 1050) std::cout << "electron" << x << "    " << y <<  "    " << z << " " << std::endl;
    //if (array_index / 3 == 1050) std::cout << "distance from atom" << d_x << "    " << d_y <<  "    " << d_z << " " << std::endl;
         //atoms go ±1 from there
    /*    for (int a = 0; a < 1331; a++) {
            x_mod = (atomic_size * (a % 11)) - (5*atomic_size);
            y_mod = (atomic_size * ((int(floor(a/11))) % 11)) - (5*atomic_size);
            z_mod = (atomic_size * floor(a / 121)) - (5*atomic_size);
            x_distance = x_mod - d_x;
            y_distance = y_mod - d_y;
            z_distance = z_mod - d_z;  */

        if (x_distance < -30)     x_distance = x_distance + 40;
        else if (x_distance > 30) x_distance = x_distance - 40;
        if (y_distance < -30)     y_distance = y_distance + 40;
        else if (y_distance > 30) y_distance = y_distance - 40;
        if (z_distance <  -30)    z_distance = z_distance + 40;
        else if (z_distance > 30) z_distance = z_distance - 40; 

        length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
        if(length > 250) continue;
        nearest_atom_list[e][neighbor_count] = a;
        neighbor_count++;
        if (length > 100) continue;

        length = sqrtl(length);
           // #pragma omp critical
         //    std::cout << x_distance << ", " << y_distance << ", " << z_distance << ", " << length << std::endl;
            //if (array_index / 3 == 1050) std::cout << "distance from lattice atom" << x_distance << "    " << y_distance <<  "    " << z_distance << " " << length << std::endl;
           // if (length < 0.00001) length = 0.00001;
          //  if(length < 0.2) chosen_electron++;
        
       // std::cout << length << std::endl;
               // terminaltextcolor(RED);
                //std::cout << "Scattering event at electron " << array_index/3 << ". Length: " << length << std::endl;
                //terminaltextcolor(WHITE);
           // }
        force = -28*((3.3* 3.3* expl(-3.3 * length)) - (expl(-1 * length)));
                        //q*k*k * exp(-15(A**-1) * length (A));
    
        PE += 28*((3.3*expl(-3.3*length)) - (expl(-1*length)));

           
       /*     if (array_index/3 == chosen_electron){
                if (x_distance > 0) {
                    int bin = int(floor(x_distance *10));
                    if(bin < 0) std::cout << "e-a, e " << array_index/3 << " bin " << bin <<  " length " << length << " x distance " << x_distance << std::endl;
                 //   charge_distrib[bin] += -1*(2*((420*exp(-15*length)) - (exp(-1*length))))*x_distance / length;
                }
            } */
         /*   x_distance *= -1;
            y_distance *= -1;
            z_distance *= -1; */

            //if (array_index / 3 == 1050) std::cout << "old force" << x_force << "  " << y_force << "   " << z_force << "   " << force << std::endl;
        phi   = acosl(z_distance / length);
        theta = atanl(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        x_force += force * cosl(theta)*sinl(phi);
        y_force += force * sinl(theta)*sinl(phi);
        z_force += force * cosl(phi);
            //if (array_index / 3 == 1050) std::cout << "new force" << x_force << "  " << y_force << "   " << z_force << "   " << force << std::endl;

        if( length < 2 && !collision && current_time_step > 10) {
            
            l_x = x_distance;
            l_y = y_distance;
            l_z = z_distance;

            excitation_constant = abs(velocity_gaussian_distrib(gen));
            Px = electron_velocity[array_index];
            Py = electron_velocity[array_index+1];
            Pz = electron_velocity[array_index+2];
            P = sqrtl((Px*Px)+(Py*Py)+(Pz*Pz));
            P_p = sqrtl(Px*Px+Py*Py);
            P_i = sqrtl(Pz*Pz);
            scattering_velocity = P;
            electron_TE = electron_potential[e]*constants::K_A + P*P*constants::m_e_r/2;
            
            d_p = sqrtl(l_x*l_x+l_y*l_y);
            d_i = sqrtl(l_z*l_z);
            d_r = sqrtl((l_x*l_x)+(l_y*l_y)+(l_z*l_z));
            polar_value = M_PI * range_distrib(gen);
            incline_value = M_PI * range_distrib(gen);
            normal_polar_angle = acosl((l_x*Px+l_y*Py)/(P_p*d_p));
            normal_incline_angle = acosl((l_z*Pz) / (P_i*d_i));

            excitation_energy = excitation_constant*(electron_TE - E_f_A);
            scattering_velocity = P - sqrtl(2*excitation_energy/constants::m_e_r);
          //  std::cout << electron_TE << ", " << E_f*1e20 << std::endl;
            if ( excitation_energy > 0 && scattering_velocity > 0) {
      //  std::cout << electron_TE << ", " << atomic_phonon_energy[2*a] << ", " << excitation_energy << ", " << electron_TE - excitation_energy << std::endl;
              //  excitation_energy = electron_TE - excitation_energy;
        
               
       // std::cout << P - scattering_velocity << std::endl;
      // if(scattering_velocity < 0) std::cout << P << ", " << scattering_velocity << ", " << P - scattering_velocity <<  std::endl;
    
                if(normal_polar_angle > M_PI) normal_polar_angle = -1*(normal_polar_angle-M_PI);
                else if (normal_polar_angle < -1*M_PI) normal_polar_angle = -1*(normal_polar_angle+M_PI);

                if(normal_incline_angle > M_PI) normal_incline_angle = -1*(normal_incline_angle-M_PI);
                else if (normal_incline_angle < -1*M_PI) normal_incline_angle = -1*(normal_incline_angle+M_PI);

                polar_scattering_angle = atanl(Py/Px);
                if(polar_scattering_angle < 0) polar_scattering_angle += M_PI;
                incline_scattering_angle = acosl(Pz/P);

                polar_prob = excitation_constant*polar_value*polar_value* ( ( ((M_PI/2)+normal_polar_angle) * exp( (polar_value-1)*(polar_value-1)/(-8)) )+( ((M_PI/2)-normal_polar_angle) * exp( (polar_value+1)*(polar_value+1)/(-8)) ) )/(d_r*2*M_PI*sqrt(2*M_PI));
                incline_prob = excitation_constant*incline_value*incline_value* ( ( ((M_PI/2)+normal_incline_angle)  * exp( (incline_value-1)*(incline_value-1)/(-8)) ) + ( ((M_PI/2)-normal_incline_angle) * exp( (incline_value-1)*(incline_value-1)/(-8)) ) )/(d_r*2*M_PI*sqrt(2*M_PI));
             //   std::cout << polar_prob << ", " << incline_prob << std::endl;
                if(scattering_prob_distrib(gen) < polar_prob) {
                    polar_scattering_angle = polar_value - normal_polar_angle;
      //  scattering_velocity *= velocity_gaussian_distrib(gen);
                    collision = true;
                }
                if(scattering_prob_distrib(gen) < incline_prob) {
                    incline_scattering_angle = incline_value - normal_incline_angle;
        //scattering_velocity *= ;
                    collision = true;
                }
   // std::cout << scattering_velocity << std::endl;
                if(collision) {
                   // if(scattering_velocity > P) std::cout << scattering_velocity << ", " << P << std::endl;
                   
                    #pragma omp critical
                    {
                    chosen_electron++;
                    atomic_phonon_energy[a] += excitation_energy;
                    TLE += excitation_energy;
                   // std::cout << excitation_energy << std::endl;
                   // std::cout << excitation_energy*1e-20 << " , " << sqrtl(1e10*2*excitation_energy/constants::m_e) << ", " << 1e5*(P - scattering_velocity) << ", " << (P*P - scattering_velocity*scattering_velocity)*constants::m_e*1e10/2 << std::endl;
                    }
                    P = scattering_velocity;
                    electron_velocity[array_index]   = P * cos(polar_scattering_angle)*sin(incline_scattering_angle);
                    electron_velocity[array_index+1] = P * sin(polar_scattering_angle)*sin(incline_scattering_angle);
                    electron_velocity[array_index+2] = P * cos(incline_scattering_angle); 
                } 
        //#pragma omp critical
     //   std::cout << "Scattering Velocity: " << scattering_velocity << ", incoming_velocity" << P << std::endl; //", polar_probability " << polar_prob << ", incline_prob " << incline_prob << ", normal_polar_angle " << normal_polar_angle << ", " << ", polar_scattering_angle " << polar_scattering_angle << ", incline_scattering_angle" << incline_scattering_angle << std::endl;
        //#pragma omp critical
      //  TLE += P*P - scattering_velocity*scattering_velocity;
       // #pragma omp critical
        //if(TLE != 0) std::cout << TLE << std::endl;
     //   std::cout << scattering_velocity*1e5 << std::endl;
              
            }
            //    e_a_scattering(e, a, x_distance,y_distance,z_distance);
               // e_p_scattering(e, a, x_distance, y_distance, z_distance);
              // e_p_scattering(e, a, x_distance,y_distance,z_distance);
        }
          /*  if(length < 0.9) {
               // e_e_scattering(array_index/3, x_distance,y_distance,z_distance);
                long double capture = capture_chance_distrib(gen);
                if( capture < (exp(-1*pow(length / 0.3, 2)))) {
                    if (!electron_capture[array_index/3]) {
                        #pragma omp critical 
                        {
                       //std::cout << length << ", " << exp(-1*pow(length / 0.3, 2)) << ", " << capture << std::endl;
                        captured_electron_list.push_back(array_index/3);
                        captured_electron_list.push_back(length);

                        }
                       // std::cout << "Electron " << array_index/3 << " captured." << std::endl;
                       electron_capture[array_index/3] = true;
                    }
                }
            } */
    }
    new_electron_potential[e] = PE;
    nearest_atom_list[e][0] = neighbor_count;
  //std::cout << electron_potential[e] << std::endl;
    return PE;
}

long double neighbor_e_a_coulomb(int e, int array_index, long double& x_force, long double& y_force, long double& z_force, const long double& x, const long double& y, const long double& z) {
  //set e-a attraction
        //calculate nearest neighbor;
     std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<long double> range_distrib(-1,1);
    std::uniform_real_distribution<long double> scattering_prob_distrib(0,1);
    std::normal_distribution<long double> velocity_gaussian_distrib(0.1,0.1);

    long double  x_distance,y_distance,z_distance,d_x,d_y,d_z, l_x,l_y,l_z, x_mod,y_mod,z_mod, length, force;
    long double phi,theta, PE = 0;
    int atom_array, neighbor_count = 1;
    bool collision= false;

    long double excitation_constant;// = velocity_gaussian_distrib(gen);
    long double Px;// = electron_velocity[array_index];
    long double Py;// = electron_velocity[array_index+1];
    long double Pz;// = electron_velocity[array_index+2];
    long double P;// = sqrtl(Px*Px+Py*Py+Pz*Pz);
    long double P_p;// = sqrtl(Px*Px+Py*Py);
    long double P_i;// = sqrtl(Pz*Pz);
    long double scattering_velocity;// = P;
    long double electron_TE;// = 1e10*(electron_potential[e]*constants::K + P*P*constants::m_e*0.5);
   
    long double d_p;// = sqrtl(l_x*l_x+l_y*l_y);
    long double d_i;// = sqrtl(l_z*l_z);
    long double d_r;// = sqrtl((l_x*l_x)+(l_y*l_y)+(l_z*l_z));
    long double polar_value;// = M_PI * range_distrib(gen);
    long double incline_value;// = M_PI * range_distrib(gen);
    long double normal_polar_angle;// = acosl((l_x*Px+l_y*Py)/(P_p*d_p));
    long double normal_incline_angle ;//= acosl((l_z*Pz) / (P_i*d_i));
    long double incline_prob;
    long double polar_prob;
    long double incline_scattering_angle;
    long double polar_scattering_angle;
    long double excitation_energy; //Angstroms
    // = abs(excitation_constant)*(electron_TE - atomic_phonon_energy[2*a]);

    int size = 1+nearest_atom_list[e][0];
    /*std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<long double> capture_chance_distrib(0,1); */

/*    d_x = x - ((atomic_size * round((x-1) / atomic_size)) + 1); //closest x atom index
    d_y = y - ((atomic_size * round((y-1) / atomic_size)) + 1); //closest y atom index
    d_z = z - ((atomic_size * round((z-1) / atomic_size)) + 1); //closest z atom index
*/
    for (int a = 1; a < size; a++) {
       // if(nearest_atom_list[e][a] < 0) break;
        atom_array = nearest_atom_list[e][a]*3;
        x_distance = atom_position[atom_array]   - x;
        y_distance = atom_position[atom_array+1] - y;
        z_distance = atom_position[atom_array+2] - z; 

    //mean_radius[(array_index/3)*2] = radius;
    //mean_radius[(array_index/3)*2 + 1] =  28*((3.3*expl(-3.3*radius)) - (expl(-1*radius)));
    //if (array_index / 3 == 1050) std::cout << "electron" << x << "    " << y <<  "    " << z << " " << std::endl;
    //if (array_index / 3 == 1050) std::cout << "distance from atom" << d_x << "    " << d_y <<  "    " << d_z << " " << std::endl;
         //atoms go ±1 from there
     /*   for (int a = 0; a < 1331; a++) {
            x_mod = (atomic_size * (a % 11)) - (5*atomic_size);
            y_mod = (atomic_size * ((int(floor(a/11))) % 11)) - (5*atomic_size);
            z_mod = (atomic_size * floor(a / 121)) - (5*atomic_size);
            x_distance = x_mod - d_x;
            y_distance = y_mod - d_y;
            z_distance = z_mod - d_z; */

        if (x_distance < -30)     x_distance = x_distance + 40;
        else if (x_distance > 30) x_distance = x_distance - 40;
        if (y_distance < -30)     y_distance = y_distance + 40;
        else if (y_distance > 30) y_distance = y_distance - 40;
        if (z_distance <  -30)    z_distance = z_distance + 40;
        else if (z_distance > 30) z_distance = z_distance - 40;

        length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
        if (length > 100) continue;

        length = sqrtl(length);
           // #pragma omp critical
         //    std::cout << x_distance << ", " << y_distance << ", " << z_distance << ", " << length << std::endl;
            //if (array_index / 3 == 1050) std::cout << "distance from lattice atom" << x_distance << "    " << y_distance <<  "    " << z_distance << " " << length << std::endl;
           // if (length < 0.00001) length = 0.00001;
          //  if(length < 0.2) chosen_electron++;
            
               // terminaltextcolor(RED);
                //std::cout << "Scattering event at electron " << array_index/3 << ". Length: " << length << std::endl;
                //terminaltextcolor(WHITE);
           // }
        force = -28*((3.3* 3.3* expl(-3.3 * length)) - (expl(-1 * length)));
                        //q*k*k * exp(-15(A**-1) * length (A));
    
        PE += 28*((3.3*expl(-3.3*length)) - (expl(-1*length)));
          //   std::cout << 28*((3.3*expl(-3.3*length)) - (expl(-1*length))) << std::endl;
       /*     if (array_index/3 == chosen_electron){
                if (x_distance > 0) {
                    int bin = int(floor(x_distance *10));
                    if(bin < 0) std::cout << "e-a, e " << array_index/3 << " bin " << bin <<  " length " << length << " x distance " << x_distance << std::endl;
                 //   charge_distrib[bin] += -1*(2*((420*exp(-15*length)) - (exp(-1*length))))*x_distance / length;
                }
            } */
         /*   x_distance *= -1;
            y_distance *= -1;
            z_distance *= -1; */

            //if (array_index / 3 == 1050) std::cout << "old force" << x_force << "  " << y_force << "   " << z_force << "   " << force << std::endl;
        phi   = acosl(z_distance / length);
        theta = atanl(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        x_force += force * cosl(theta)*sinl(phi);
        y_force += force * sinl(theta)*sinl(phi);
        z_force += force * cosl(phi);
            //if (array_index / 3 == 1050) std::cout << "new force" << x_force << "  " << y_force << "   " << z_force << "   " << force << std::endl;

        if(length < 2 && !collision && current_time_step > 10 ) {
               // e_a_scattering(e, a, x_distance,y_distance,z_distance);
           
            l_x = x_distance;
            l_y = y_distance;
            l_z = z_distance;

            excitation_constant = abs(velocity_gaussian_distrib(gen));
            Px = electron_velocity[array_index];
            Py = electron_velocity[array_index+1];
            Pz = electron_velocity[array_index+2];
            P = sqrtl(Px*Px+Py*Py+Pz*Pz);
            P_p = sqrtl((Px*Px)+(Py*Py));
            P_i = sqrtl(Pz*Pz);
            scattering_velocity = P;
            electron_TE = electron_potential[e]*constants::K_A + P*P*constants::m_e_r/2;
            
            d_p = sqrtl((l_x*l_x)+(l_y*l_y));
            d_i = sqrtl(l_z*l_z);
            d_r = sqrtl((l_x*l_x)+(l_y*l_y)+(l_z*l_z));
            polar_value = M_PI * range_distrib(gen);
            incline_value = M_PI * range_distrib(gen);
            normal_polar_angle = acosl(((l_x*Px)+(l_y*Py))/(P_p*d_p));
            normal_incline_angle = acosl((l_z*Pz) / (P_i*d_i));

            excitation_energy = excitation_constant*(electron_TE - E_f_A);
             scattering_velocity = P - sqrtl(2*excitation_energy/constants::m_e_r);
          //  std::cout << electron_TE << ", " << E_f*1e20 << std::endl;
            if ( excitation_energy > 0 && scattering_velocity > 0)  {
      //  std::cout << electron_TE << ", " << atomic_phonon_energy[2*a] << ", " << excitation_energy << ", " << electron_TE - excitation_energy << std::endl;
              //  excitation_energy = electron_TE - excitation_energy;
        
            
                    
       // std::cout << P - scattering_velocity << std::endl;
      // if(scattering_velocity < 0) std::cout << P << ", " << scattering_velocity << ", " << P - scattering_velocity <<  std::endl;
    

                if(normal_polar_angle > M_PI) normal_polar_angle = -1*(normal_polar_angle-M_PI);
                else if (normal_polar_angle < -1*M_PI) normal_polar_angle = -1*(normal_polar_angle+M_PI);

                if(normal_incline_angle > M_PI) normal_incline_angle = -1*(normal_incline_angle-M_PI);
                else if (normal_incline_angle < -1*M_PI) normal_incline_angle = -1*(normal_incline_angle+M_PI);

    
                polar_scattering_angle = atanl(Py/Px);
                if(polar_scattering_angle < 0) polar_scattering_angle += M_PI;
                incline_scattering_angle = acosl(Pz/P);

                polar_prob = excitation_constant*polar_value*polar_value* ( ( ((M_PI/2)+normal_polar_angle) * exp( (polar_value-1)*(polar_value-1)/(-8)) )+( ((M_PI/2)-normal_polar_angle) * exp( (polar_value+1)*(polar_value+1)/(-8)) ) )/(d_r*2*M_PI*sqrt(2*M_PI));
                incline_prob = excitation_constant*incline_value*incline_value* ( ( ((M_PI/2)+normal_incline_angle)  * exp( (incline_value-1)*(incline_value-1)/(-8)) ) + ( ((M_PI/2)-normal_incline_angle) * exp( (incline_value-1)*(incline_value-1)/(-8)) ) )/(d_r*2*M_PI*sqrt(2*M_PI));
             //    std::cout << d_p << ", " << P_p << ", " << normal_polar_angle << ", " << normal_incline_angle << std::endl;
                if(scattering_prob_distrib(gen) < polar_prob) {
                    polar_scattering_angle = polar_value - normal_polar_angle;
      //  scattering_velocity *= velocity_gaussian_distrib(gen);
                    collision = true;
                }
                if(scattering_prob_distrib(gen) < incline_prob) {
                    incline_scattering_angle = incline_value - normal_incline_angle;
        //scattering_velocity *= ;
                    collision = true;
                }
   // std::cout << scattering_velocity << std::endl;
                if(collision) {
                   
                    #pragma omp critical
                    {
                    chosen_electron++;
                    atomic_phonon_energy[a] += excitation_energy;
                    TLE += excitation_energy;
                   // std::cout << excitation_energy << std::endl;
               //    std::cout << excitation_energy*1e-20 << " , " << sqrtl(1e-20*2*excitation_energy/constants::m_e) << ", " << 1e5*(P - scattering_velocity) << ", " << (P-scattering_velocity)*(P-scattering_velocity)*constants::m_e*1e10/2 << std::endl;
                    }
                    // P = scattering_velocity;
                      electron_velocity[array_index]   = scattering_velocity * cos(polar_scattering_angle)*sin(incline_scattering_angle);
                electron_velocity[array_index+1] = scattering_velocity* sin(polar_scattering_angle)*sin(incline_scattering_angle);
                electron_velocity[array_index+2] = scattering_velocity * cos(incline_scattering_angle); 

                //#pragma omp critical
              //  std::cout << sqrtl((electron_velocity[array_index]*electron_velocity[array_index]) + (electron_velocity[array_index+1]*electron_velocity[array_index+1]) + (electron_velocity[array_index+2]*electron_velocity[array_index+2])) << ", " << P-scattering_velocity << ", " << sqrtl(2*excitation_energy/constants::m_e_r) << std::endl;
                }
        //#pragma omp critical
     //   std::cout << "Scattering Velocity: " << scattering_velocity << ", incoming_velocity" << P << std::endl; //", polar_probability " << polar_prob << ", incline_prob " << incline_prob << ", normal_polar_angle " << normal_polar_angle << ", " << ", polar_scattering_angle " << polar_scattering_angle << ", incline_scattering_angle" << incline_scattering_angle << std::endl;
        //#pragma omp critical
      //  TLE += P*P - scattering_velocity*scattering_velocity;
       // #pragma omp critical
        //if(TLE != 0) std::cout << TLE << std::endl;
     //   std::cout << scattering_velocity*1e5 << std::endl;
              
            }
                //e_p_scattering(e, a, x_distance, y_distance, z_distance);
        }
          /*  if(length < 0.9) {
               // e_e_scattering(array_index/3, x_distance,y_distance,z_distance);
                long double capture = capture_chance_distrib(gen);
                if( capture < (exp(-1*pow(length / 0.3, 2)))) {
                    if (!electron_capture[array_index/3]) {
                        #pragma omp critical 
                        {
                       //std::cout << length << ", " << exp(-1*pow(length / 0.3, 2)) << ", " << capture << std::endl;
                        captured_electron_list.push_back(array_index/3);
                        captured_electron_list.push_back(length);

                        }
                       // std::cout << "Electron " << array_index/3 << " captured." << std::endl;
                       electron_capture[array_index/3] = true;
                    }
                }
            } */
    }
    new_electron_potential[e] = PE;
  //  std::cout << electron_potential[e] << std::endl;
    return PE;
}

long double electron_e_e_coulomb(int e, int array_index, long double& x_force, long double& y_force, long double& z_force, const long double& x, const long double& y, const long double& z) {
    
    int array_index_i;
    long double x_distance,y_distance,z_distance,length;
    long double force, theta,phi, PE = 0;
    int neighbor_count = 1;

  //set e-e repulsion
      //  #pragma omp parallel for
   //     std::ofstream nearest_neighbors;
     //   int count = 0;
    //    if(e==100) nearest_neighbors.open("CASTLE/Electron_Neighbors.txt");
    
        for (int i = 0; i < conduction_electrons; i++) {
            if (i == e) continue; //no self repulsion
            //  if (symmetry_list[e][i]) continue;  //make use of symmetry
            array_index_i = 3*i;
            //   electron_spin_two = conduction_electron_spin[i];
            //   array_index_i = 3*i;
            //   array_index_i_y = array_index_i + 1;
            //   array_index_i_z = array_index_i + 2;

            x_distance = x - new_electron_position[array_index_i];
            y_distance = y - new_electron_position[array_index_i + 1];
            z_distance = z - new_electron_position[array_index_i + 2]; 

            if (x_distance < -30)     x_distance = x_distance + 40;
            else if (x_distance > 30) x_distance = x_distance - 40;
            if (y_distance < -30)     y_distance = y_distance + 40;
            else if (y_distance > 30) y_distance = y_distance - 40;
            if (z_distance <  -30)    z_distance = z_distance + 40;
            else if (z_distance > 30) z_distance = z_distance - 40;

            length = ((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance));
            
            if (length > 250) continue;
            nearest_neighbor_list[e][neighbor_count] = i;
            neighbor_count++;
             
          
            if (length > 100) continue; 
           // count++;
        //   if (e==100) nearest_neighbors << length << "\n";
            length = sqrtl(length);
            /*   if (electron_spin == electron_spin_two) {
                force = 1e10 * ((constants::K / (length * length * constants::m_e))- (635 * constants::kB / (length * constants::m_e)));
                
            }
            else { */
            force = 1 / (length * length);
            // if(array_index / 3 == 0) std::cout << "force " << force << std::endl;
            PE += force * length;
            
            
            
       /*     if (e == chosen_electron){
                if(x_distance < 0) {
                    int bin = int(floor(abs(x_distance)*10));
                    if(bin<0) std::cout << "e-e " << e << " bin " << bin <<  " length " << length << " x distance " << x_distance << std::endl;
                   // charge_distrib[bin] += x_distance * length;
                }
            } */
            phi   = acosl(z_distance / length);
            theta = atanl(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            x_force += force * cosl(theta)*sinl(phi);
            y_force += force * sinl(theta)*sinl(phi);
            z_force += force * cosl(phi);

            /*    //make use of symmetry
            new_force_array[array_index_i]   += -1 * x_unit;
            new_force_array[array_index_i_y] += -1 * y_unit;
            new_force_array[array_index_i_z] += -1 * z_unit;
    
            symmetry_list[i][e] = true; //set symmetry flag */
        }
    nearest_neighbor_list[e][0] = neighbor_count;
  //  nearest_neighbors << neighbor_count << std::endl;
    //nearest_neighbors.close();
  //  if(e==100) std::cout << neighbor_count << std::endl;
    new_electron_potential[e] += PE;
    return PE / 2;
}

long double neighbor_e_e_coulomb(int e, int array_index, long double &x_force, long double &y_force, long double &z_force, const long double& x, const long double& y, const long double& z) {
    long double x_distance,y_distance,z_distance, length;
    long double force, theta,phi, PE = 0;
    int size = 1+nearest_neighbor_list[e][0]; //.size();
    int array_index_i;

  //  std::ofstream nearest_electrons;
    //int count = 0;
  //  if(e==100) nearest_electrons.open("CASTLE/Nearest_Electron.txt");
        
   // if(e==100) std::cout << size << std::endl;
    for (int i=1; i < size; i++) { //a better way exists to iterate this, just can't remember
     //   if(nearest_neighbor_list[e][i] < 0) { //just in case
           // if(e==100) std::cout << i << std::endl;
       //     break; //terminate when reach end of nearest electron
        //}
        
        array_index_i = nearest_neighbor_list[e][i] * 3;

        x_distance = x - new_electron_position[array_index_i]; //don't get caching but cut down iterations by large margin
        y_distance = y - new_electron_position[array_index_i + 1];
        z_distance = z - new_electron_position[array_index_i + 2];
    
        if (x_distance < -30)     x_distance = x_distance + 40;
        else if (x_distance > 30) x_distance = x_distance - 40;
        if (y_distance < -30)     y_distance = y_distance + 40;
        else if (y_distance > 30) y_distance = y_distance - 40;
        if (z_distance <  -30)    z_distance = z_distance + 40;
        else if (z_distance > 30) z_distance = z_distance - 40;
        
        length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
        if (length > 100) continue; //exclude zone 2
       // if(e==100) nearest_electrons << length << "\n";
        //    count++;
        
        length = sqrtl(length);
            /*   if (electron_spin == electron_spin_two) {
                force = 1e10 * ((constants::K / (length * length * constants::m_e))- (635 * constants::kB / (length * constants::m_e)));
                
            }
            else { */
        force = 1 / (length * length);
            // if(array_index / 3 == 0) std::cout << "force " << force << std::endl;
        PE += force * length;
        
        
        
     /*   if (e == chosen_electron){
            if(x_distance < 0) {
                int bin = int(floor(abs(x_distance)*10));
                 if(bin<0) std::cout << "e-e " << e << " bin " << bin <<  " length " << length << " x distance " << x_distance << std::endl;
             //   charge_distrib[bin] += x_distance *force;
            }
        } */

        phi   = acosl(z_distance / length);
        theta = atanl(y_distance / x_distance);
        if (x_distance < 0) theta += M_PI;

        x_force += force * cosl(theta)*sinl(phi);
        y_force += force * sinl(theta)*sinl(phi);
        z_force += force * cosl(phi);

    }
    
   // nearest_electrons << size << std::endl;
    //nearest_electrons.close();
  //  if(e==100) std::cout << count << std::endl;
    new_electron_potential[e] += PE;
   
    return PE / 2;
}

long double electron_applied_voltage(int array_index, long double& x_force, long double& y_force, long double& z_force) {
    
    x_force -= 2;

    return 2; //-10.0 * 1e-10 / constants::e ;
}


} //end CASTLE namespace

