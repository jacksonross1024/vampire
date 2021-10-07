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


int velocity_verlet_step(double dt) {
    
  
dt = 1e-20;
TPE = 0;
TKE = 0;



        if (err::check) std::cout << "Calculating verlet integration...Time step: " << dt << "\n";
    // std::ofstream acceleration_output;
    
       
    std::ofstream electron_position_output_up;
    std::ofstream electron_position_output_down;
    std::ofstream electron_velocity_output;
    std::ofstream electron_spin_output;

    setup_output();
    update_position();
    update_forces();
    update_velocity();
    output_data();
        if(err::check) std::cout << "Initializing output files..." << "\n";
    if (current_time_step % CASTLE_output_rate == 0) {
       
       
        CASTLE_output_data = true;
        std::string time_stamp = std::to_string(current_time_step);

        electron_position_output_up.open("CASTLE/Electron_Position_Up/" + time_stamp + "U.xyz");
            electron_position_output_up << total_spin_up << "\n";
            electron_position_output_up << time_stamp << "\n";

        electron_position_output_down.open("CASTLE/Electron_Position_Down/" + time_stamp + "D.xyz");
            electron_position_output_down << total_spin_down << "\n";
            electron_position_output_down << time_stamp << "\n";

        electron_velocity_output.open("CASTLE/CASTLE_Electron_velocities/" + time_stamp + ".txt");
            electron_velocity_output << "Electron number    x-component     y-component    z-component     length" << "\n";

        electron_spin_output.open("CASTLE/Electron_Spin/" + time_stamp + ".txt");
            electron_spin_output << conduction_electrons << "\n";
            electron_spin_output << "       1 is up, 0 is down" << "\n";
   // acceleration_output.open("CASTLE/test" + time_stamp);

        
    }

     int array_index_z, array_index_y, array_index = 0; //equivalent to current electron x component vector

    
   // #pragma omp parallel for
            if (err::check) std::cout << "Updating new electron position." << "\n";
       // std::cout << "conduction electrons" << conduction_electrons << "\n";

    int array_index_i, array_index_i_y, array_index_i_z = 0;
    double x,y,z, a_x, a_y, a_z, d_x, d_y, d_z, x_distance, y_distance, z_distance = 0.0;
    double velocity_length = 0.0;
    double x_unit, y_unit, z_unit, length, force, modifier = 0.0;
    bool electron_spin = false;
    bool electron_spin_two = false;
    for (int e = 0; e < conduction_electrons; e++){ 
        array_index = 3*e;
        array_index_y = array_index + 1;
        array_index_z = array_index + 2;

      // std::cout << "why are you dying...." << e << "\n";
    
        new_electron_position[array_index]   = electron_position[array_index]   + (electron_velocity[array_index]   * dt) + (electron_force[array_index]   * dt * dt * 0.5); // x superarray component
        new_electron_position[array_index_y] = electron_position[array_index_y] + (electron_velocity[array_index_y] * dt) + (electron_force[array_index_y] * dt * dt * 0.5); // y superarray component
        new_electron_position[array_index_z] = electron_position[array_index_z] + (electron_velocity[array_index_z] * dt) + (electron_force[array_index_z] * dt * dt * 0.5); // z superarray component
    

       //  if (e == 0) std::cout << new_electron_position[array_index] << "   " << electron_position[array_index]  << "   " << (electron_velocity[array_index]   * dt) << "   " << (electron_force[array_index]   * dt * dt * 0.5) << "\n";
        symmetry_list[e].resize(conduction_electrons, false);

        if (CASTLE_output_data)   {
            if (!conduction_electron_spin[e]) electron_position_output_up << "H" << "    " << new_electron_position[array_index] << "    " << new_electron_position[array_index_y] << "    " << new_electron_position[array_index_z] << "\n";
          //  else                              electron_position_output_down << "H" << "    " << new_electron_position[array_index] << "    " << new_electron_position[array_index_y] << "    " << new_electron_position[array_index_z] << "\n";
        } 
  //  #pragma omp parallel for
    }
    
            if (err::check) std::cout << "Positions updated... next step and array index: " << "\n";


    //set forces local variables


    //forces and velocity integration
   // #pragma omp parallel for
            if (err::check) std::cout << "Forces, spins, and velocities update..." << "\n";
    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
 //   std::uniform_int_distribution<> Spin_distrib(1, 10000);
    std::uniform_int_distribution<> Vel_distrib(0, 360);
    for (int e = 0; e < conduction_electrons; e++) {
        array_index = 3*e;
        array_index_y = array_index + 1;
        array_index_z = array_index + 2;

        //spontaneous spin flip
        
        velocity_length = sqrt( (electron_velocity[array_index]*electron_velocity[array_index]) + (electron_velocity[array_index_y]*electron_velocity[array_index_y]) + (electron_velocity[array_index_z]*electron_velocity[array_index_z]) ); //meters
        
    //    double spin_chance = Spin_distrib(gen) * 0.0002;
        symmetry_list[e].resize(conduction_electrons, false);
        
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

       
        
    //set e-a attraction
        //calculate nearest neighbor;
    
        x = new_electron_position[array_index];
        y = new_electron_position[array_index_y];
        z = new_electron_position[array_index_z];

        a_x = atomic_size * round(x / atomic_size); //closest x atom index
        a_y = atomic_size * round(y / atomic_size); //closest y atom index
        a_z = atomic_size * round(z / atomic_size); //closest z atom index

        d_x = x - a_x; //beware the sign change in the next step
        d_y = y - a_y; //beware the sign change in the next step
        d_z = z - a_z; //beware the sign change in the next step
    
  
         //atoms go ±1 from there
        for (int a = 0; a < 9; a++) {
            modifier = (atomic_size * (a % 3)) - atomic_size;
            x_distance = modifier - d_x;
            y_distance = modifier - d_y;
            z_distance = modifier - d_z;
            
            length = sqrt((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance)); //meters
            if (length < 0.00001) length = 0.00001;
            
            if (length > screening_depth) {
                force = 1e-10 / (length * length * length * length * constants::m_e); //Angstroms
                TPE += -2 * force * length;
            } else {
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
            }

            new_force_array[array_index]   += force * x_distance / length;
            new_force_array[array_index_y] += force * y_distance / length;
            new_force_array[array_index_z] += force * z_distance / length;
        }
  
    //set e-e repulsion
      //  #pragma omp parallel for
        for (int i = 0; i < conduction_electrons; i++) {
            if (i == e) continue; //no self repulsion
            if (symmetry_list[e][i]) continue;  //make use of symmetry

         //   electron_spin_two = conduction_electron_spin[i];
            array_index_i = 3*i;
            array_index_i_y = array_index_i + 1;
            array_index_i_z = array_index_i + 2;

            x_distance = x - new_electron_position[array_index_i];
            y_distance = y - new_electron_position[array_index_i_y];
            z_distance = z - new_electron_position[array_index_i_z]; 


            length = sqrt((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance));
            if (length < 0.00001) length = 0.00001;
            
         /*   if (electron_spin == electron_spin_two) {
                force = 1e10 * ((constants::K / (length * length * constants::m_e))- (635 * constants::kB / (length * constants::m_e)));
                
            }
            else { */
                force = 1 / (length * length);

            
            TPE += -2 * force * length;

            x_unit = x_distance * force / length;
            y_unit = y_distance * force / length;
            z_unit = z_distance * force / length;
            

            new_force_array[array_index]     += x_unit;
            new_force_array[array_index + 1] += y_unit;
            new_force_array[array_index + 2] += z_unit;

            //make use of symmetry
            new_force_array[array_index_i]   += -1 * x_unit;
            new_force_array[array_index_i_y] += -1 * y_unit;
            new_force_array[array_index_i_z] += -1 * z_unit;
    
            symmetry_list[i][e] = true; //set symmetry flag
        }
        for (int i = 0; i < lattice_electrons.size(); i++) {
            array_index_i = 3*i;
       //     electron_spin_two = lattice_electron_spin[i];
            x_distance = x - lattice_electrons[array_index_i];
            y_distance = y - lattice_electrons[array_index_i_y];
            z_distance = z - lattice_electrons[array_index_i_z];

            length = sqrt((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance));
            if (length < 0.00001) length = 0.00001;
            
         /*   if (electron_spin == electron_spin_two) {
                force = 1e10 * ((constants::K / (length * length * constants::m_e))- (635 * constants::kB / (length * constants::m_e)));
                
            }
            else { */
                force = 1/ (length * length);

            
            TPE += -2 * force * length;

            new_force_array[array_index]   += x_distance * force / length; //Angstroms
            new_force_array[array_index_y] += y_distance * force / length;
            new_force_array[array_index_z] += z_distance * force / length;
        }
  
    

        new_electron_velocity[array_index]   = electron_velocity[array_index]   + (electron_force[array_index]   + new_force_array[array_index])   * dt * 0.5; 
        new_electron_velocity[array_index_y] = electron_velocity[array_index_y] + (electron_force[array_index_y] + new_force_array[array_index_y]) * dt * 0.5;
        new_electron_velocity[array_index_z] = electron_velocity[array_index_z] + (electron_force[array_index_z] + new_force_array[array_index_z]) * dt * 0.5;

        velocity_length = sqrt( (new_electron_velocity[array_index]*new_electron_velocity[array_index]) + (new_electron_velocity[array_index_y]*new_electron_velocity[array_index_y]) + (new_electron_velocity[array_index_z]*new_electron_velocity[array_index_z]) );
        TKE += 0.5 * velocity_length * velocity_length;

        //lattice boundary conserves momentum and is one atom less than actual lattice to preserve boundary effect
        if (new_electron_position[array_index]   < (0.00001 + 3*atomic_size) || new_electron_position[array_index]   > (lattice_width - 3*atomic_size - 0.00001) ) new_electron_velocity[array_index]   *= -1.0;
        if (new_electron_position[array_index_y] < (0.00001 + 3*atomic_size) || new_electron_position[array_index_y] > (lattice_depth - 3*atomic_size - 0.00001) ) new_electron_velocity[array_index_y] *= -1.0;
        if (new_electron_position[array_index_z] < (0.00001 + 3*atomic_size) || new_electron_position[array_index_z] > (lattice_height - 3*atomic_size - 0.00001)) new_electron_velocity[array_index_z] *= -1.0;

        if (CASTLE_output_data) {
            electron_velocity_output << e << "  " << new_electron_velocity[array_index] << "    " << new_electron_velocity[array_index_y] << "  " << new_electron_velocity[array_index_z] << "  " << velocity_length << "\n";
         //   acceleration_output << e << "   " <<(electron_force[array_index] << " " <<(electron_force[array_index_y] <<   "   " <<(electron_force[array_index_z] << "\n";  

        }  
    }
         //   if (err::check) std::cout << "Velocity increased..." << "\n";
    
        current_time_step += 1;
        electron_position = new_electron_position;
        electron_force = new_force_array;
        electron_velocity = new_electron_velocity;

        electron_position_output_up.close();
        electron_position_output_down.close();
        electron_velocity_output.close();
        electron_spin_output.close();
       
        fill(new_force_array.begin(), new_force_array.end(), 0.0);

        if (equilibrium_step) {
        mean_data_array[current_time_step*5]     = TKE * constants::K * 1e-20;
        mean_data_array[current_time_step*5 + 1] = TPE * constants::K * 1e-10;
        mean_data_array[current_time_step*5 + 2] = TPE * 1e-10/ TKE;
        mean_data_array[current_time_step*5 + 3] = sqrt(2 * TPE * 1e-10 / constants::m_e);
     //   mean_data_array[current_time_step*5 + 4] = (total_spin_up - total_spin_down);
        }
        else {
        mean_data_array[0]      += TKE * constants::K * 1e-20 / conduction_electrons;
        mean_data_array[1] += TPE * constants::K * 1e-10 / conduction_electrons;
        mean_data_array[2] += TPE * 1e-10/ TKE;
        mean_data_array[3] += sqrt(2 * TPE * 1e-10 / constants::m_e);
     //   mean_data_array[4] += (total_spin_up - total_spin_down);
       
        }
        //    if (err::check) std::cout << "Acceleration arrays reset" << "\n";
        if (CASTLE_output_data) std::cout << (current_time_step / total_time_steps) * 100 << "%" << "\n";
        CASTLE_output_data = false;
    return EXIT_SUCCESS;
}
void setup_output() {

}

void update_position(){

}

void update_forces() {

}

void update_cells() {

}

void update_velocity() {

}

void electron_e_a_coulomb() {

}

void electron_e_e_coulomb() {

}


} //end CASTLE namespace

