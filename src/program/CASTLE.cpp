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


#include "CASTLE.hpp"

namespace CASTLE {


void create() {
        if (err::check) std::cout << "Creating CASTLE..." << std::endl; 
        if (err::check) std::cout << "  ";
  //  atoms::num_atoms = 15 * 15 * 15;
    conduction_electrons = 1728; //15 * 15 * 15;

   // sim::temperature = sim::Teq;
  //  sim::temperature 
  
    temperature = 300;
    equilibrium_step = true;
        if (err::check) std::cout << "Prepare to initialize..." << std::endl;
        
    total_time_steps = sim::equilibration_time; //100
    loop_time = sim::loop_time;
    lattice_atoms = 20 * 20 * 20;

    std::ofstream lattice_output;
    std::ofstream electron_position_output_up;
    std::ofstream electron_position_output_down;
    std::ofstream electron_velocity_output;
    std::ofstream mean_data;
    std::ofstream electron_spin_output;

    std::cout << "Building CASTLE..." <<std::endl;
    stopwatch_t castle_watch;
    castle_watch.start();
    
    mean_data.open("CASTLE/mean_data");
    mean_data << "step  mean-KE     mean-PE     total-KE    total-PE    mean-spin" << std::endl;
    initialize(atoms::num_atoms, conduction_electrons);
    std::cout << "CASTLE build time[s]: " << castle_watch.elapsed_seconds() << std::endl;
    std::cout << "Storming CASTLE..." << std::endl;
   
    sim::integrate(total_time_steps);
     std::cout << "Average time step[s]:  " << (castle_watch.elapsed_seconds()) / total_time_steps << std::endl;
    //create_gif();
    for (int t = 0; t < total_time_steps; t++) {
        mean_data << t << "     " << mean_data_array[5*t] << "  " << mean_data_array[5*t + 1] << "  " << mean_data_array[5*t + 2] << "  " << mean_data_array[5*t + 3] << "  " << mean_data_array[5*t + 4] << std::endl;  
    }
    mean_data.close();

    std::cout << "Equilibrium step complete. Averaging CASTLE..." << std::endl;
     double current_time = castle_watch.elapsed_seconds();
    equilibrium_step = false;
    current_time_step = 0;
    total_time_steps = sim::loop_time;

    mean_data.open("CASTLE/loop_data");
    mean_data << "mean-KE     mean-PE     total-KE    total-PE    mean-spin" << std::endl;
    mean_data_array.resize(5, 0.0);
    sim::integrate(loop_time);
    mean_data << mean_data_array[0] / loop_time  << "  " << mean_data_array[1] / loop_time  << "    " << mean_data_array[2] / loop_time << "  " << mean_data_array[3] / loop_time << "  " << mean_data_array[4] / loop_time << std::endl;  
    mean_data.close();

    std::cout << "Averaging complete. " << current_time - castle_watch.elapsed_seconds() << " [s] elapsed." << std::endl;
}

void initialize (const double num_atoms, const double num_electrons) {
        // string lattice_type) { awaiting future development

            if (err::check) std::cout << "Initializing CASTLE..."  << std::endl;

    //initialize all arrays and variables
    
    lattice_output.open("CASTLE/CASTLE_Lattice.xyz");
     // output lattice atoms and locations
    lattice_output << lattice_atoms << std::endl;
    lattice_output << std::endl;

    
    electron_position_output_up.open("CASTLE/Electron_position_Up/initU.xyz");
    electron_position_output_up << conduction_electrons << std::endl;  
    electron_position_output_up << "Initial positions for spin up electrons" << std::endl;  

    electron_position_output_down.open("CASTLE/Electron_Position_Down/initD.xyz");
    electron_position_output_down << conduction_electrons << std::endl;  
    electron_position_output_down << "Initial positions for spin down electrons" << std::endl;  

    
    electron_velocity_output.open("CASTLE/Electron_Velocities/init.txt");
    electron_velocity_output << "electron number    x-component     y-component     z-component     length" << std::endl;  

   
    electron_spin_output.open("CASTLE/Electron_Spin/init.txt");
   
    double x, y, z, a_x, a_y, a_z, d_x, d_y, d_z, x_distance, y_distance, z_distance, e_distance_x, e_distance_y, e_distance_z;
    double x_unit, y_unit, z_unit, length, force, mean_velocity = 0;
    int array_index_i;

            if (err::check) std::cout << "Lattice output file and electron position file opened..." << std::endl;

    current_time_step = 0;
    lattice_height = 38.0; //A
    lattice_width = 38.0; // A
    lattice_depth = 38.0; // A

    CASTLE_output_rate = 1;

    lattice_atoms = 20 * 20 * 20;
    double electrons = 20 * 20 * 20; //15 * 15 * 15;
    conduction_electrons = 1728.0;
    atomic_size = 2.0; //Angst diameter. 
    screening_depth = atomic_size * 0.25 * 0.5 * 1e-10; //meters

    atom_position.resize(lattice_atoms * 3, 0.0);
    double n_f = 1e30 * electrons / (lattice_width * lattice_height * lattice_depth);
    E_f = 3 * constants::h * constants::h * electrons * pow((3 * n_f / (8 * M_PI)), 0.666666667) / (10 * constants::m_e); //Fermi-energy
    mu_f = 5 * E_f / (3 * electrons);//Fermi-level
    v_f = sqrt(2 * E_f / constants::m_e); //meters
    TKE = 0;
    TPE = 0;
  //  total_spin_up = 0;
    total_spin_down = 0;

   // std::cout << "n_f" << n_f << "mu_f: " << mu_f << " E_f: " << E_f << " v_f: " << v_f <<  std::endl;

    electron_position.resize(conduction_electrons * 3, 0.0); // ""'Memory is cheap. Time is expensive' -Steve Jobs; probably" -Michael Scott." -Headcannon.
    new_electron_position.resize(conduction_electrons * 3, 0.0);
    lattice_electrons.resize((electrons - conduction_electrons), 0.0);
    electron_velocity.resize(conduction_electrons * 3, 0.0); //Angstroms
    new_electron_velocity.resize(conduction_electrons * 3, 0.0); //Angstroms
    electron_acc.resize(conduction_electrons * 3, 0.0); //current and future arrays
    new_acc_array.resize(conduction_electrons * 3, 0.0);
    mean_data_array.resize(total_time_steps*5 + 5, 0.0);
    conduction_electron_spin.resize(conduction_electrons, false);
    lattice_electron_spin.resize((electrons - conduction_electrons), false);

    int array_index = 0;
    symmetry_list.resize(conduction_electrons);
     //   std::cout << "conduction electrons" << conduction_electrons << std::endl;
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> Pos_distrib(1, 360);

    for (int a = 0; a < lattice_atoms; ++a) {
        array_index = 3*a;
        atom_position[array_index]     = atomic_size * (a % 20); //lattice creation awaiting future development
        atom_position[array_index + 1] = atomic_size * ((int(floor(a / 20))) % 20);
        atom_position[array_index + 2] = atomic_size * floor(a / 400);
        lattice_output << "Ni" << "     " << atom_position[array_index] << "     " << atom_position[array_index + 1] << "   " << atom_position[array_index + 2] << std::endl;  
    }
    lattice_output.close();

        if (err::check) std::cout << "Lattice built " << std::endl;

    double velocity_length = 0.0;
    double phi = 0.0;
    double theta = 0.0;
    double  modifier, x_pos, y_pos, z_pos = 0.0;
    std::srand(std::time(nullptr));
    //super loop for each conducting electron
    int array_index_l = 0;
    conduction_electrons = 0;
    for (int e = 0; e < electrons; e++) {
        //    if (err::check) std::cout << "Raising structs..." << std::endl;
        
        //random program for velocity initialization
        

        phi = M_PI * Pos_distrib(gen) / 180;
        theta = M_PI * Pos_distrib(gen) / 180.0;
        //    if (err::check) std::cout << "Prepare to set positions: " << conduction_electrons * 3 * total_time_steps << std::endl;

        //initialize and output electron posititons
        x_pos =  atomic_size * (e % 20); 
        y_pos =  atomic_size * ((int(floor(e / 20))) % 20);
        z_pos =  atomic_size * floor(e/ 400);
     
        if (x_pos < 6.001 || x_pos > 31.999) {
            array_index_l += 1;
            lattice_electrons[array_index_l]     = x_pos;
            lattice_electrons[array_index_l + 1] = y_pos;
            lattice_electrons[array_index_l + 2] = z_pos;
           // if (Pos_distrib(gen) > 180) lattice_electron_spin[array_index_l] = true;
        }
        else if (y_pos < 6.001 || y_pos > 31.999) {
            array_index_l += 1;
            lattice_electrons[array_index_l]     = x_pos;
            lattice_electrons[array_index_l + 1] = y_pos;
            lattice_electrons[array_index_l + 2] = z_pos;
          //  if (Pos_distrib(gen) > 180) lattice_electron_spin[array_index_l] = true;
        }
        else if (z_pos < 6.001 || z_pos > 31.999) {
            array_index_l += 1;
            lattice_electrons[array_index_l]     = x_pos;
            lattice_electrons[array_index_l + 1] = y_pos;
            lattice_electrons[array_index_l + 2] = z_pos;
          //  if (Pos_distrib(gen) > 180) lattice_electron_spin[array_index_l] = true;
        } else {
            
            array_index = conduction_electrons*3;
            electron_position[array_index]     = x_pos + sin(phi)/10;
            electron_position[array_index + 1] = y_pos + cos(theta)/10;
            electron_position[array_index + 2] = z_pos + sin(theta)/10;

         /*   if (Pos_distrib(gen) > 180) { 
                conduction_electron_spin[conduction_electrons] = true;
                total_spin_up += 1;
                electron_position_output_up << "H" << "    " << electron_position[array_index] << "    " << electron_position[array_index + 1] << "    " << electron_position[array_index + 2] << std::endl;  
            } else { */
                total_spin_down += 1;
                electron_position_output_down << "H" << "    " << electron_position[array_index] << "    " << electron_position[array_index + 1] << "    " << electron_position[array_index + 2] << std::endl;  
            
            conduction_electrons +=1;
          
        
           
            //std::random_device rd;  //Will be used to obtain a seed for the random number engine
            //std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_int_distribution<> Vel_distrib(0, 360);
            int random_angle = Vel_distrib(gen);
            phi = M_PI * random_angle / 360.0;
            theta = M_PI * random_angle / 180.0;

          //  if (err::check) std::cout << "Random angles chosen. v_f: %d. cos(theta): %d. sin(theta): %d. sin(0.5 - phi): %d" << v_f << cos(theta) << sin(theta) << sin(0.5 - phi) << std::endl;

            electron_velocity[array_index]     = cos(theta)*sin(phi) * v_f * 1e10; //gotta get back to Angstrom's
            electron_velocity[array_index + 1] = sin(theta)*sin(phi) * v_f * 1e10;
            electron_velocity[array_index + 2] = cos(phi)            * v_f * 1e10;
            velocity_length = 1e-10 * sqrt( (electron_velocity[array_index]*electron_velocity[array_index]) + (electron_velocity[array_index + 1]*electron_velocity[array_index + 1]) + (electron_velocity[array_index + 2]*electron_velocity[array_index + 2]) ); //meters
        //    std::cout << "v_f" << v_f << " velocity " << velocity_length << std::endl;
            electron_velocity_output <<  conduction_electrons << "      " << electron_velocity[array_index] << "    " << electron_velocity[array_index + 1] << "    " << electron_velocity[array_index + 2] << "    " << velocity_length*1e10 << std::endl;
          // if (err::check) std::cout << "Velocities randomized..." << std::endl;
         //   electron_spin_output << conduction_electrons << "  " << conduction_electron_spin[conduction_electrons - 1] << std::endl;
        }
    }
                
        //    std::uniform_int_distribution<> Spin_distrib(0, 10000);
    for (int e = 0; e < conduction_electrons; e++) {
        //spontaneous spin flip

         //   double spin_chance = Spin_distrib(gen) * 0.0002;
            symmetry_list[e].resize(conduction_electrons, false);
        
      // if (e == 0) std::cout << "deltaV " << deltaV << " eps " << (0.5 * constants::m_e * velocity_length * velocity_length) / (constants::kB * temperature) << " E_f " <<  E_f / (constants::kB * temperature) << std::endl;
    /*    double flip_chance = 1.0;
        if (((0.5 * constants::m_e * velocity_length * velocity_length) - E_f) <  0) flip_chance = exp(((0.5 * constants::m_e * velocity_length * velocity_length) - E_f) / (constants::kB * temperature));
       // std::cout << "flip chance " << flip_chance << " spin chance " << spin_chance << std::endl;
        if (spin_chance > flip_chance) {
            double deltaV = sqrt(2 * mu_f / constants::m_e) * 1e10;
            bool old_spin = conduction_electron_spin[e];
            bool new_spin = conduction_electron_spin[e] = !old_spin;
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
        TKE += 0.5 * constants::m_e * velocity_length * velocity_length * 1e20; //energy in Angstroms

            //calculate e-a force
        x = electron_position[array_index];
        y = electron_position[array_index + 1];
        z = electron_position[array_index + 2];

        a_x = atomic_size * round(x / atomic_size); //closest x atom index
        a_y = atomic_size * round(y / atomic_size); //closest y atom index
        a_z = atomic_size * round(z / atomic_size); //closest z atom index

        d_x = x - a_x;
        d_y = y - a_y;
        d_z = z - a_z;
      
            //cube around electron
        for (int a = 0; a < 9; a++) {
            modifier = (atomic_size * (a % 3)) - atomic_size;
            x_distance = modifier - d_x;
            y_distance = modifier - d_y;
            z_distance = modifier - d_z;

            length = 1e-10 * sqrt((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance));
             if (length < 1e-14) length = 1e-14;
           
            if (length > screening_depth) {
                force = 1e30 * constants::K / (length * length * length * length * constants::m_e);
                TPE += -2 * force * length * 1e10; //Angstroms
            } else {
                force = -1e10 * constants::K / (length * length * constants::m_e);
                TPE += -2 * force * length * 1e10;
            }

            electron_acc[array_index]     += 1e-10 * force * x_distance / length;
            electron_acc[array_index + 1] += 1e-10 * force * y_distance / length;
            electron_acc[array_index + 2] += 1e-10 * force * z_distance / length;

        }
}
            
             if (err::check) std::cout << "Nearest atoms set..." << std::endl;

    electron_position_output_up.close();
    electron_position_output_down.close();
    electron_velocity_output.close();
    electron_spin_output.close();
    //set e-e repulsion local variables
            if (err::check) std::cout << "Initializing e-e repulsion symmetry." << std::endl;
       //     std::cout << "conduction electrons" << conduction_electrons << std::endl;
  
   // std::ofstream acceleration_output;
   // acceleration_output.open("CASTLE/test");
  //  bool electron_spin;
  //  bool electron_spin_two;
    for (int e = 0; e < conduction_electrons; e++) {
        array_index = 3*e;
    
       //set e-e repulsion
        e_distance_x = electron_position[array_index];
        e_distance_y = electron_position[array_index + 1];
        e_distance_z = electron_position[array_index + 2];

      //  electron_spin = conduction_electron_spin[e];
            if (err::check) if(e ==0) std::cout << "Calculating conduction electron repulsion" << std::endl;
        
        for (int i = 0; i < conduction_electrons; i++) {
            if (i == e) continue; //no self repulsion
            if (symmetry_list[e][i]) continue;  //make use of symmetry

         //   electron_spin_two = conduction_electron_spin[i];

            array_index_i = 3*i;
            x_distance = e_distance_x - electron_position[array_index_i];
            y_distance = e_distance_y - electron_position[array_index_i + 1];
            z_distance = e_distance_z - electron_position[array_index_i + 2];

            length = 1e-10 * sqrt((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance));
            //if (length < 0.000001) length = 0.000001;

          /*  if (electron_spin == electron_spin_two) {
                force = 1e10 * ((constants::K / (length * length * constants::m_e))- (635 * constants::kB / (length * constants::m_e)));
               // std::cout << " 1, 0 " << 1e10 * (constants::K / (length * length * constants::m_e)) << " delta1,1 " << (635 * constants::kB / (length * constants::m_e)) << std::endl;
            }
            else { */
            force = 1e10 * constants::K / (length * length * constants::m_e);
            
            TPE += -2 * force * length * 1e10;

            x_unit = 1e-10 * x_distance * force / length;
            y_unit = 1e-10 * y_distance * force / length;
            z_unit = 1e-10 * z_distance * force / length;
            //    if (err::check) std::cout << "Distances grabbed" << std::endl; 

            electron_acc[array_index]     += x_unit;
            electron_acc[array_index + 1] += y_unit;
            electron_acc[array_index + 2] += z_unit;

            //make use of symmetry
            electron_acc[array_index_i]     += -1 * x_unit;
            electron_acc[array_index_i + 1] += -1 * y_unit;
            electron_acc[array_index_i + 2] += -1 * z_unit;
            symmetry_list[i][e] = true; //apply symmetry flag


            //  if (e ==0)   std::cout << electron_acc[array_index] << std::endl;
        } 
                if(err::check) if(e ==0) std::cout << "Calculating conduction-lattice repulsion" << std::endl;
        for (int i = 0; i < lattice_electrons.size(); i++) {

        //    electron_spin_two = lattice_electron_spin[i];

            array_index_i = 3*i;
            x_distance = e_distance_x - lattice_electrons[array_index_i];
            y_distance = e_distance_y - lattice_electrons[array_index_i + 1];
            z_distance = e_distance_z - lattice_electrons[array_index_i + 2];

            length = 1e-10 * sqrt((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance));
            //if (length < 0.000001) length = 0.000001;

         /*   if (electron_spin == electron_spin_two) {
                force = 1e10 * ((constants::K / (length * length * constants::m_e))- (635 * constants::kB / (length * constants::m_e)));
                
            }
            else { */
                force = 1e10 * constants::K / (length * length);

            
            TPE += -2 * force * length * 1e10;

            electron_acc[array_index]     += 1e-10 * x_distance * force / length;
            electron_acc[array_index + 1] += 1e-10 * y_distance * force / length;
            electron_acc[array_index + 2] += 1e-10 * z_distance * force / length;
        }
        //electron potential walls of lattice
       /*
            electron_acc[array_index]     += 2 / ((e_distance_x - 1.99) * (e_distance_x - 1.99));
            electron_acc[array_index]     += -2 / ((28.01 - e_distance_x) * (28.01 - e_distance_x));
            electron_acc[array_index + 1] += 2 / ((e_distance_y - 1.99) * (e_distance_y - 1.99));
            electron_acc[array_index + 1] += -2 / ((28.01 - e_distance_y) * (28.01 - e_distance_y));
            electron_acc[array_index + 2] += 2 / ((e_distance_z - 1.99) * (e_distance_z - 1.99));
            electron_acc[array_index + 2] += -2 / ((28.01 - e_distance_z) * (28.01 - e_distance_z));
        */
    //acceleration_output << e << "   " << electron_acc[array_index] << " " << electron_acc[array_index + 1] <<   "   " << electron_acc[array_index + 2] << std::endl;  
    }
            if (err::check) std::cout << "Forces ready..." << std::endl;
    
            if (err::check) std::cout << "Initialization complete." << std::endl;
     mean_data_array[0] = TKE ;
     mean_data_array[1] = TPE;
     mean_data_array[2] = TPE / TKE;
     mean_data_array[3] = 1e10 * sqrt(2 * TPE / constants::m_e);
  //   mean_data_array[4] = (total_spin_up - total_spin_down);
     
}
 


void create_gif() {

}

} //end of CASTLE namespace



