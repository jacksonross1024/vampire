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
    
    dt = 1e-19;
    TPE = 0;
    TKE = 0;

            if (err::check) std::cout << "Calculating verlet integration...Time step: " << dt << std::endl;

   // std::ofstream electron_position_output_up;
    //std::ofstream electron_position_output_down;
    //std::ofstream electron_velocity_output;
    //std::ofstream electron_spin_output;

            if(err::check) std::cout << "Initializing output files..." << std::endl;
    if (equilibrium_step && current_time_step % CASTLE_output_rate == 0) setup_output();

            if (err::check) std::cout << "Updating new electron position." << std::endl;
    update_position();

            if (err::check) std::cout << "Forces, spins, and velocities update..." << std::endl;
    update_dynamics();

            if (err::check) std::cout << "Output mean data" << std::endl;
    if (CASTLE_output_data)   output_data();



    //reset integration
    current_time_step += 1;

    electron_position = new_electron_position;
    electron_force = new_force_array;
    electron_velocity = new_electron_velocity;

   
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
    std::string time_stamp = std::to_string(current_time_step);

      //  electron_position_output_up.open("CASTLE/Electron_Position_Up/" + time_stamp + "U.xyz");
        //    electron_position_output_up << total_spin_up << "\n";
          //  electron_position_output_up << time_stamp << "\n";

    electron_position_output_down.open("CASTLE/Electron_Position/" + time_stamp + ".xyz");
    electron_position_output_down << conduction_electrons << "\n";
    electron_position_output_down << time_stamp << "\n";

    electron_velocity_output.open("CASTLE/Electron_Velocity/" + time_stamp + ".txt");
    electron_velocity_output << "Electron number    x-component     y-component    z-component     length" << "\n";

    
      //  electron_spin_output.open("CASTLE/Electron_Spin/" + time_stamp + ".txt");
       //     electron_spin_output << conduction_electrons << "\n";
       //     electron_spin_output << "       1 is up, 0 is down" << "\n";
        // acceleration_output.open("CASTLE/test" + time_stamp);
 
}


void update_position(){

    int array_index,array_index_y,array_index_z = 0;
    double x_pos,y_pos,z_pos = 0.0;
    #pragma omp parallel for private(array_index,array_index_y,array_index_z, x_pos,y_pos,z_pos) schedule(static)
    for (int e = 0; e < conduction_electrons; e++){ 
        array_index = 3*e;
        array_index_y = array_index + 1;
        array_index_z = array_index + 2;

      // std::cout << "why are you dying...." << e << "\n";
    
        x_pos = electron_position[array_index]   + (electron_velocity[array_index]   * dt) + (electron_force[array_index]   * dt * dt * 0.5  * constants::K / constants::m_e); // x superarray component
        y_pos = electron_position[array_index_y] + (electron_velocity[array_index_y] * dt) + (electron_force[array_index_y] * dt * dt * 0.5  * constants::K / constants::m_e); // y superarray component
        z_pos = electron_position[array_index_z] + (electron_velocity[array_index_z] * dt) + (electron_force[array_index_z] * dt * dt * 0.5  * constants::K / constants::m_e); // z superarray component

        if (x_pos < 0.0) x_pos += 40.0;
        else if (x_pos > 40.0) x_pos -= 40.0;

        if (y_pos < 0.0) y_pos += 40.0;
        else if (y_pos > 40.0) y_pos -= 40.0;

        if (z_pos < 0.0) z_pos += 40.0;
        else if (z_pos > 40.0) z_pos -= 40.0;

        new_electron_position[array_index]   = x_pos;
        new_electron_position[array_index_y] = y_pos;
        new_electron_position[array_index_z] = z_pos;

       // if (e == 0) std::cout << new_electron_position[array_index] << "   " << electron_position[array_index]  << "   " << (electron_velocity[array_index]   * dt) << "   " << (electron_force[array_index]   * dt * dt * 0.5  * constants::K / constants::m_e) << "\n";
        //symmetry_list[e].resize(conduction_electrons, false);
    }
}

void update_dynamics() {
   // double applied_electronic_field = {0.0, 0.0, 1.0, 1.0}; //x, y, z, strength
    int array_index;
    double x_force,y_force,z_force, x,y,z;
    #pragma omp parallel for private(array_index, x_force,y_force,z_force, x,y,z)  schedule(static) reduction(+:TPE,TKE)
    for (int e = 0; e < conduction_electrons; e++) {
        array_index = 3*e;
        x_force = y_force = z_force = 0.0;
        x = new_electron_position[array_index];
        y = new_electron_position[array_index + 1];
        z = new_electron_position[array_index + 2];

        //if (e == 1050) std::cout << x << "  " << y << " " << z <<   "    " << std::endl;
        TPE += electron_e_a_coulomb(array_index, x_force,y_force,z_force, x,y,z);
        TPE += electron_e_e_coulomb(e, array_index, x_force,y_force,z_force, x,y,z);
        TPE += electron_applied_voltage(array_index, x_force,y_force,z_force);
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

}

long double update_velocity(int array_index) {
        int array_index_y = array_index + 1;
        int array_index_z = array_index + 2;

     //   if (e == 0) std::cout << new_electron_velocity[array_index] << " " << electron_velocity[array_index] << "    " <<  electron_force[array_index]  << "    " << new_force_array[array_index]  << "    " <<  dt * 0.5  * constants::K / constants::m_e << "\n"; 
        new_electron_velocity[array_index]   = electron_velocity[array_index]   + ((electron_force[array_index]   + new_force_array[array_index])   * dt * 0.5  * constants::K / constants::m_e); 
        new_electron_velocity[array_index_y] = electron_velocity[array_index_y] + ((electron_force[array_index_y] + new_force_array[array_index_y]) * dt * 0.5  * constants::K / constants::m_e);
        new_electron_velocity[array_index_z] = electron_velocity[array_index_z] + ((electron_force[array_index_z] + new_force_array[array_index_z]) * dt * 0.5  * constants::K / constants::m_e);
    
    long double velocity_length = (new_electron_velocity[array_index]*new_electron_velocity[array_index]) + (new_electron_velocity[array_index_y]*new_electron_velocity[array_index_y]) + (new_electron_velocity[array_index_z]*new_electron_velocity[array_index_z]); //Angstroms
    return (0.5 * velocity_length);
}

long double electron_e_a_coulomb(int array_index, double& x_force, double& y_force, double& z_force, const double& x, const double& y, const double& z) {
  //set e-a attraction
        //calculate nearest neighbor;
    double d_x,d_y,d_z, x_mod,y_mod,z_mod, x_distance,y_distance,z_distance, length  = 0.0;
    long double PE, force = 0.0;
    d_x = x - ((atomic_size * round(x / atomic_size)) + 1); //closest x atom index
    d_y = y - ((atomic_size * round(y / atomic_size)) + 1); //closest y atom index
    d_z = z - ((atomic_size * round(z / atomic_size)) + 1); //closest z atom index

    //if (array_index / 3 == 1050) std::cout << "electron" << x << "    " << y <<  "    " << z << " " << std::endl;
    //if (array_index / 3 == 1050) std::cout << "distance from atom" << d_x << "    " << d_y <<  "    " << d_z << " " << std::endl;
         //atoms go ±1 from there
        for (int a = 0; a < 27; a++) {
            x_mod = (atomic_size * (a % 3)) - atomic_size;
            y_mod = (atomic_size * ((int(floor(a/3))) % 3)) - atomic_size;
            z_mod = (atomic_size * floor(a / 9)) - atomic_size;
            x_distance = x_mod - d_x;
            y_distance = y_mod - d_y;
            z_distance = z_mod - d_z;


            length = sqrt((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance)); //Angstroms
            //if (array_index / 3 == 1050) std::cout << "distance from lattice atom" << x_distance << "    " << y_distance <<  "    " << z_distance << " " << length << std::endl;
            if (length < 0.0001) length = 0.0001;
            
            force = (28 * exp(-20 * length)) - exp(-length);
                /*
                velocity_length = sqrt( (electron_velocity[array_index]*electron_velocity[array_index]) + (electron_velocity[array_index_y]*electron_velocity[array_index_y]) + (electron_velocity[array_index_z]*electron_velocity[array_index_z]) );
                int angle = Vel_distrib(gen);
                double theta = M_PI * angle / 180;
                double phi = M_PI * angle / 360;
                electron_velocity[array_index]   = cos(theta)*sin(phi) * velocity_length;
                electron_velocity[array_index_y] = sin(theta)*sin(phi) * velocity_length;
                electron_velocity[array_index_z] = cos(phi)            * velocity_length;

                new_electron_position[array_index]   = a_x + (cos(theta)*sin(phi) * screening_depth);
                new_electron_position[array_index_y] = a_y + (sin(theta)*sin(phi) * screening_depth);
                new_electron_position[array_index_z] = a_z + (cos(phi)            * screening_depth);

                length = screening_depth;

                force = -1e10 / (length * length * constants::m_e);
                TPE += -2 * force * length;
            } */
            PE += exp(-1) - ((7 / 5) * exp(-20*length));
            }
            //if (array_index / 3 == 1050) std::cout << "old force" << x_force << "  " << y_force << "   " << z_force << "   " << force << std::endl;
            x_force += force * x_distance / length;
            y_force += force * y_distance / length;
            z_force += force * z_distance / length;
            //if (array_index / 3 == 1050) std::cout << "new force" << x_force << "  " << y_force << "   " << z_force << "   " << force << std::endl;
        
    return PE * 1e10;
}

long double electron_e_e_coulomb(int e, int array_index, double& x_force, double& y_force, double& z_force, const double& x, const double& y, const double& z) {
    
    int array_index_i;
    double d_x,d_y,d_z, x_mod,y_mod,z_mod, x_distance,y_distance,z_distance, length = 0.0;
    long double force, PE = 0.0;
  //set e-e repulsion
      //  #pragma omp parallel for
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
            if (length > 400) continue;

            length = sqrt(length);
            /*   if (electron_spin == electron_spin_two) {
                force = 1e10 * ((constants::K / (length * length * constants::m_e))- (635 * constants::kB / (length * constants::m_e)));
                
            }
            else { */
            force = 1 / (length * length);
            
            PE += force * length;

            x_force += x_distance * force / length;
            y_force += y_distance * force / length;
            z_force += z_distance * force / length;

            /*    //make use of symmetry
            new_force_array[array_index_i]   += -1 * x_unit;
            new_force_array[array_index_i_y] += -1 * y_unit;
            new_force_array[array_index_i_z] += -1 * z_unit;
    
            symmetry_list[i][e] = true; //set symmetry flag */
        }
    return (-0.5 * PE * 1e10);
}

long double electron_applied_voltage(int array_index, double& x_force, double& y_force, double& z_force) {
    
    z_force -= 1e-11;

    return -1.0;
}


} //end CASTLE namespace

