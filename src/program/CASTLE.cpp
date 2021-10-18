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
  

        std::cout << "Building CASTLE..." <<std::endl;

    stopwatch_t castle_watch;
    castle_watch.start();
    
    equilibrium_step = true;
    //========
    // Initialize the lattice, electrons, and starting forces
    //========
            if (err::check) std::cout << "Prepare to initialize..." << std::endl;

    initialize();
        //omp_set_dynamic(0);
        //omp_set_num_threads(8);
        std::cout << "CASTLE build time[s]: " << castle_watch.elapsed_seconds() << std::endl;
        #pragma omp parallel 
        {
          //  std::cout << "OpenMP capability detected. Parallelizing integration. Thread " << omp_get_thread_num() <<  " of threads: " << omp_get_num_threads() << std::endl;
        }
        std::cout << "Storming CASTLE..." << std::endl;
   
    
    //========
    // Integrate total time steps
    //========
    castle_watch.start();
    sim::integrate(total_time_steps);
    
        std::cout << "Average time step[s]:  " << (castle_watch.elapsed_seconds()) / total_time_steps << std::endl;

        std::cout << "Equilibrium step complete. Averaging CASTLE..." << std::endl;
    
    //=========
    // Set up CASTLE for loop averaging step
    //=========
    castle_watch.start();
    equilibrium_step = false;
    CASTLE_output_rate = 1;
    total_time_steps = sim::loop_time;
    //dt = 1e-20;

    //========
    // Run averaging step
    //========
    sim::integrate(total_time_steps);

    //========
    // Output data
    //========

    mean_data.close();
        std::cout << "Averaging complete. " << castle_watch.elapsed_seconds() << " [s] elapsed." << std::endl;
}

//====================================
// Initialize function to set up CASTLE
//====================================
void initialize () {
       
            if (err::check) std::cout << "Initializing CASTLE..."  << std::endl;

    
    //========
    // Initialize lattice
    //========
    initialize_lattice();
            if (err::check) std::cout << "Lattice built " << std::endl;
    
    //========
    // Initialzie variables used by all functions
    //========

    //=========
    // Grab simulation variables from VAMPIRE
    //=========
    conduction_electrons = 20*20*20;  //sim::conduction_electrons;
    CASTLE_output_rate = 10; //sim::CASTLE_output_rate;
    dt = 1e-18;
    temperature = 300; //sim::temperature;
    total_time_steps = sim::equilibration_time; //100
    x_flux = 0.0;
    current_time_step = 0;
    CASTLE_real_time = 0.0;
    
    //========
    // initialize electrons: lattice and conduction bands, velocity, spin, etc.
    //=======
    initialize_electrons();
            if (err::check) std::cout << "Electrons ready..." << std::endl;

    //=======
    // Calls forces set up
    //=======
    initialize_forces();
            if (err::check) std::cout << "Forces ready..." << std::endl;


    mean_data.open("CASTLE/mean_data.csv");
    mean_data << "time, step, mean-KE, mean-PE, mean-TE, -lambda, +lambda, mean-lambda, mean-current_density, mean-drift_velocity, current, resistance" << "\n";

}

//====================================
// Creates and outputs atomic lattice
//      Currently static lattice
//====================================
void initialize_lattice() {
    
    atomic_size = 2.0; // sim::atomic_size; //Angst diameter. 
    screening_depth = 0.36; //sim::screening_depth; //Angstroms 
    lattice_atoms = 20 * 20 * 20; //Better lattice creation will come from VAMPIRE in future

    lattice_height = 40.0; //A
    lattice_width  = 40.0; // A
    lattice_depth  = 40.0; // A

    lattice_output.open("CASTLE/CASTLE_Lattice.xyz");

     // output lattice atoms and locations
    lattice_output << lattice_atoms << "\n"; //xyz file requires first-line for number of elements
    lattice_output << " static lattice" "\n"; //comment line
    
    atom_position.resize(lattice_atoms * 3, 0.0);

    int array_index = 0; //local loop index variable
    for (int a = 0; a < lattice_atoms; ++a) {  
        array_index = 3*a;
        atom_position[array_index]     = 1+ atomic_size * (a % 20); //lattice creation awaiting future development
        atom_position[array_index + 1] = 1+ atomic_size * ((int(floor(a / 20))) % 20);
        atom_position[array_index + 2] = 1+ atomic_size * floor(a / 400);
        lattice_output << "Ni" << "     " << atom_position[array_index] << "     " << atom_position[array_index + 1] << "   " << atom_position[array_index + 2] << "\n";  
    }
    lattice_output.close();

}

//====================================
// Function call to initialize basic electron structures
//      Will set up lattice and conduction electrons and velocites
//      Any additional features, e.g. spin, will need to be added here through function calls
//====================================
void initialize_electrons() {

    // Open output isntances
    /*    electron_position_output_up.open("CASTLE/Electron_position_Up/initU.xyz");
        electron_position_output_up << conduction_electrons << std::endl;  
        electron_position_output_up << "Initial positions for spin up electrons" << std::endl;  
    */

    electron_position_output_down.open("CASTLE/Electron_Position/init.xyz");
    electron_position_output_down << conduction_electrons << "\n";  
    electron_position_output_down << "Initial positions for electrons" << "\n";  

    
    electron_velocity_output.open("CASTLE/Electron_Velocity/init.txt");
    electron_velocity_output << "electron number    x-component     y-component     z-component     length" << std::endl;  
    
            if (err::check) std::cout << "Lattice output file and electron position file opened..." << std::endl;

    //=========
    // Initialize arrays for holding electron variables
    //      Arrays in super array format to take advantage of caching
    //========
    electron_position.resize(conduction_electrons * 3, 0.0); // ""'Memory is cheap. Time is expensive' -Steve Jobs; probably" -Michael Scott." -Headcannon.
    new_electron_position.resize(conduction_electrons * 3, 0.0);
    // lattice_electrons.resize((total_electrons - conduction_electrons), 0.0);
    electron_velocity.resize(conduction_electrons * 3, 0.0); //Angstroms
    new_electron_velocity.resize(conduction_electrons * 3, 0.0); //Angstroms
    electron_force.resize(conduction_electrons * 3, 0.0); //current and future arrays
    new_force_array.resize(conduction_electrons * 3, 0.0);
    //mean_data_array.resize(total_time_steps*5 + 5, 0.0);
    // conduction_electron_spin.resize(conduction_electrons, false);
    // lattice_electron_spin.resize((total_electrons - conduction_electrons), false);
    // symmetry_list.resize(conduction_electrons);
    velocity_length_hist.resize(conduction_electrons*4, 0.0);
     //   std::cout << "conduction electrons" << conduction_electrons << std::endl;
    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> Pos_distrib(1, 359);

    n_f = 1e30 * conduction_electrons / (lattice_width * lattice_height * lattice_depth); // e- / m**3
    E_f = constants::h * constants::h * pow(3 * M_PI * M_PI * n_f, 0.666666666666667) / (8 * M_PI * M_PI * constants::m_e); //Fermi-energy // meters
    mu_f = 5 * E_f / (3 * conduction_electrons);//Fermi-level //meters
    v_f = sqrt(2 * E_f / constants::m_e); //meters

    std::cout << "n_f: " << n_f <<  ", E_f: " << E_f << ", E_f(eV): " << E_f*6.242e18 << ", v_f " << v_f << std::endl;
    //TKE = 0; //total Kinetic energy, meters
    //  total_spin_up = 0;
   // total_spin_down = 0;
   

    long double phi,theta, x_pos,y_pos,z_pos, velocity_length = 0.0;
    //super loop for each conducting electron
    int array_index = 0;
    nearest_neighbor_list.resize(conduction_electrons);
    for (int e = 0; e < conduction_electrons; e++) {
        nearest_neighbor_list[e].resize(conduction_electrons * 0.5, -1);
        array_index = 3*e;
        //random program for velocity initialization
        phi   = M_PI * Pos_distrib(gen) / 180;
        theta = M_PI * Pos_distrib(gen) / 180.0;

                if (err::check) std::cout << "Prepare to set positions: " << e << std::endl;

        //initialize and output electron posititons
        x_pos = (atomic_size * (e % 20)) + (cos(theta)*sin(phi) * screening_depth); //Angstroms
        y_pos = (atomic_size * ((int(floor(e / 20))) % 20)) + (sin(theta)*sin(phi) * screening_depth); //Sets on radius of screening depth from nucleus
        z_pos = (atomic_size * floor(e/ 400)) + (cos(phi) * screening_depth);
     
        if (x_pos < 0.0) x_pos += 40.0;
        else if (x_pos > 40.0) x_pos -= 40.0;

        if (y_pos < 0.0) y_pos += 40.0;
        else if (y_pos > 40.0) y_pos -= 40.0;

        if (z_pos < 0.0) z_pos += 40.0;
        else if (z_pos > 40.0) z_pos -= 40.0;

       /* if (x_pos < 6.001 || x_pos > 31.999) {
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
        } else { */
            
            
        electron_position[array_index]     = x_pos;
        electron_position[array_index + 1] = y_pos;
        electron_position[array_index + 2] = z_pos;

         /*   if (Pos_distrib(gen) > 180) { 
                conduction_electron_spin[conduction_electrons] = true;
                total_spin_up += 1;
                electron_position_output_up << "H" << "    " << electron_position[array_index] << "    " << electron_position[array_index + 1] << "    " << electron_position[array_index + 2] << std::endl;  
            } else { */
           
        electron_position_output_down << "H" << "    " << electron_position[array_index] << "    " << electron_position[array_index + 1] << "    " << electron_position[array_index + 2] << "\n";    
            
          
            //std::random_device rd;  //Will be used to obtain a seed for the random number engine
            //std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> Vel_distrib(0, 360);
        int random_angle = Vel_distrib(gen);
        phi   = M_PI * random_angle / 360.0;
        theta = M_PI * random_angle / 180.0;

                    //  if (err::check) std::cout << "Random angles chosen. v_f: %d. cos(theta): %d. sin(theta): %d. sin(0.5 - phi): %d" << v_f << cos(theta) << sin(theta) << sin(0.5 - phi) << std::endl;

        electron_velocity[array_index]     = cos(theta)*sin(phi) * v_f * 1e10; //gotta get back to Angstrom's
        electron_velocity[array_index + 1] = sin(theta)*sin(phi) * v_f * 1e10;
        electron_velocity[array_index + 2] = cos(phi)            * v_f * 1e10;
        velocity_length = (electron_velocity[array_index]*electron_velocity[array_index]) + (electron_velocity[array_index + 1]*electron_velocity[array_index + 1]) + (electron_velocity[array_index + 2]*electron_velocity[array_index + 2]); //Angstroms
                  // std::cout << "v_f" << v_f << " velocity " << velocity_length * 1e-10 << std::endl;
        electron_velocity_output << e << "      " << 1e-10*electron_velocity[array_index] << "    " << 1e-10*electron_velocity[array_index + 1] << "    " << 1e-10*electron_velocity[array_index + 2] << "    " << 1e-20*velocity_length << std::endl;
                    // if (err::check) std::cout << "Velocities randomized..." << std::endl;
            //   electron_spin_output << conduction_electrons << "  " << conduction_electron_spin[conduction_electrons - 1] << std::endl;
            
        }
    
    //electron_position_output_up.close();
    electron_position_output_down.close();
    electron_velocity_output.close();
    //electron_spin_output.close();
                
}
/*
void initialize_cells(int electron) {
    num_cells = 64;

    electrons_in_cell.resize(num_cells);
    electron_cell.resize(total_electrons, 0);
    std::vector<std::vector<std::vector<int> > > cell_list;
    cell_list.resize(4);
    int count = 0;
    for (int j = 0; j < 4; j++) {
        cell_list[j].resize(4);
        for (int k = 0; k < 4; k++) {
            cell_list[j][k].resize(4, 0.0);
            for (int i = 0; i < 4; i ++) {
                cell_list[j][k][i] = count;
                count++;
            }
        }
    }
    x = electron_position[electron*3];
    y = electron_position[electron*3 + 1];
    z = electron_position[electron*3 + 2];
    cell_x = x / (atomic_size * 4); 
    cell_y = y / (atomic_size * 4);
    cell_z = z / (atomic_size * 4);
    x_cell = floor(cell_x);
    y_cell = floor(cell_y);
    z_cell = floor(cell_z);

    current_cell = electron_cell[electron] = cell_list[x_cell][y_cell][z_cell];
    electrons_per_cell[current_cell] += 1;
    for (int cell = 0; cell < 27; cell++) {
        current = 
        electrons_per_cell[current]
        for (int e = 0; e < )
        electron_neightbors[electron] = 
    }
    electron_neighbors[electron] 
    
} */


void initialize_forces() {

    int array_index, array_index_i = 0;
    
    double x,y,z, d_x,d_y,d_z, x_mod,y_mod,z_mod, x_distance,y_distance,z_distance, length, force = 0.0;
    double x_unit,y_unit,z_unit = 0.0; 
    for (int e = 0; e < conduction_electrons; e++) {
        array_index = 3*e;
        //spontaneous spin flip

         //   double spin_chance = Spin_distrib(gen) * 0.0002;
         //symmetry_list[e].resize(conduction_electrons, false);
        
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
        
        

            //calculate e-a force
        x = electron_position[array_index];
        y = electron_position[array_index + 1];
        z = electron_position[array_index + 2];

        d_x = x - ((atomic_size * round(x / atomic_size)) + 1); //closest x atom index
        d_y = y - ((atomic_size * round(y / atomic_size)) + 1); //closest y atom index
        d_z = z - ((atomic_size * round(z / atomic_size)) + 1); //closest z atom index

      
            //cube around electron
        for (int a = 0; a < 1331; a++) {
            x_mod = (atomic_size * (a % 11)) - 5*atomic_size;
            y_mod = (atomic_size * ((int(floor(a/11))) % 11)) - 5*atomic_size;
            z_mod = (atomic_size * floor(a / 121)) - 5*atomic_size;
            x_distance = x_mod - d_x;
            y_distance = y_mod - d_y;
            z_distance = z_mod - d_z;

            length = sqrt((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance));
            if (length < 0.00001) length = 0.00001;
        
            // force
            force = -1*((1835.08 * exp(-25.6 * length)) - (2.4*exp(-length)));
         //   if(e == 0) std::cout << "force " << force << std::endl;

            electron_force[array_index]     += force * x_distance / length;
            electron_force[array_index + 1] += force * y_distance / length;
            electron_force[array_index + 2] += force * z_distance / length;

        }

             if (err::check) std::cout << "Nearest atoms set..." << std::endl;


    //set e-e repulsion local variables
            if (err::check) std::cout << "Initializing e-e repulsion symmetry." << std::endl;
       //     std::cout << "conduction electrons" << conduction_electrons << std::endl;

    
    
       //set e-e repulsion
    

      //  electron_spin = conduction_electron_spin[e];
                if (err::check) if(e ==0) std::cout << "Calculating conduction electron repulsion" << std::endl;
        int neighbor_count = 0;
        for (int i = 0; i < conduction_electrons; i++) {
            if (i == e) continue; //no self repulsion
            //if (symmetry_list[e][i]) continue;  //make use of symmetry

            //   electron_spin_two = conduction_electron_spin[i];

            array_index_i = 3*i;
            x_distance = x - electron_position[array_index_i];
            y_distance = y - electron_position[array_index_i + 1];
            z_distance = z - electron_position[array_index_i + 2];

            if (x_distance < -30)     x_distance = x_distance + 40;
            else if (x_distance > 30) x_distance = x_distance - 40;

            if (y_distance < -30)     y_distance = y_distance + 40;
            else if (y_distance > 30) y_distance = y_distance - 40;

            if (z_distance <  -30)    z_distance = z_distance + 40;
            else if (z_distance > 30) z_distance = z_distance - 40;

            length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);

            if (length > 200) continue;
            
            nearest_neighbor_list[e][neighbor_count] = i;
            neighbor_count++;
            if (length > 100) continue;
            length = sqrt(length);
            //if (length < 0.000001) length = 0.000001;

          /*  if (electron_spin == electron_spin_two) {
                force = 1e10 * ((constants::K / (length * length * constants::m_e))- (635 * constants::kB / (length * constants::m_e)));
               // std::cout << " 1, 0 " << 1e10 * (constants::K / (length * length * constants::m_e)) << " delta1,1 " << (635 * constants::kB / (length * constants::m_e)) << std::endl;
            }
            else { */
            force = 1 / (length * length);
           // if(e == 0) std::cout << "force " << force << std::endl;
            x_unit = x_distance * force / length;
            y_unit = y_distance * force / length;
            z_unit = z_distance * force / length;
            //    if (err::check) std::cout << "Distances grabbed" << std::endl; 

            electron_force[array_index]     += x_unit;
            electron_force[array_index + 1] += y_unit;
            electron_force[array_index + 2] += z_unit;

         /*   //make use of symmetry
            electron_force[array_index_i]     += -1 * x_unit;
            electron_force[array_index_i + 1] += -1 * y_unit;
            electron_force[array_index_i + 2] += -1 * z_unit;
            symmetry_list[i][e] = true; //apply symmetry flag
        */

            //  if (e ==0)   std::cout << electron_force[array_index] << std::endl;
        } 
                if(err::check) if(e ==0) std::cout << "Calculating conduction-lattice repulsion" << std::endl;
     /*   for (int i = 0; i < lattice_electrons.size(); i++) {

        //    electron_spin_two = lattice_electron_spin[i];

            array_index_i = 3*i;
            x_distance = e_distance_x - lattice_electrons[array_index_i];
            y_distance = e_distance_y - lattice_electrons[array_index_i + 1];
            z_distance = e_distance_z - lattice_electrons[array_index_i + 2];

            length = sqrt((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance));
            //if (length < 0.000001) length = 0.000001;

            if (electron_spin == electron_spin_two) {
                force = 1e10 * ((constants::K / (length * length * constants::m_e))- (635 * constants::kB / (length * constants::m_e)));
                
            }
            else {
                force = 1e-10 / (length * length * constants::m_e);

            
            TPE += -2 * force * length;

            electron_force[array_index]     += x_distance * force / length;
            electron_force[array_index + 1] += y_distance * force / length;
            electron_force[array_index + 2] += z_distance * force / length;
        } 
        */
        //electron potential walls of lattice
       /*
            electron_force[array_index]     += 2 / ((e_distance_x - 1.99) * (e_distance_x - 1.99));
            electron_force[array_index]     += -2 / ((28.01 - e_distance_x) * (28.01 - e_distance_x));
            electron_force[array_index + 1] += 2 / ((e_distance_y - 1.99) * (e_distance_y - 1.99));
            electron_force[array_index + 1] += -2 / ((28.01 - e_distance_y) * (28.01 - e_distance_y));
            electron_force[array_index + 2] += 2 / ((e_distance_z - 1.99) * (e_distance_z - 1.99));
            electron_force[array_index + 2] += -2 / ((28.01 - e_distance_z) * (28.01 - e_distance_z));
        */
    //acceleration_output << e << "   " << electron_force[array_index] << " " << electron_force[array_index + 1] <<   "   " << electron_force[array_index + 2] << std::endl;  
    }
}

void output_data() {
        
    //=========
    // Output equilibration step data
    //=========
    mean_data.precision(10);
    electron_position_output_down.precision(10);
    electron_velocity_output.precision(10);

    mean_data << std::scientific;
    electron_position_output_down << std::fixed;
    electron_velocity_output << std::scientific;
    
    int array_index, array_index_y, array_index_z = 0;
    long double x_pos, y_pos, z_pos = 0.0;
    long double x_vel, y_vel ,z_vel, velocity_length, Hfac, lambda = 0.0;
    int num_bins = 1e6;
    int max = 1e7;
    int min = 0;
    int lattice_constant = 2;
    double width = (max - min) / num_bins;
    std::vector <int> histogram;
    histogram.resize(num_bins*4, 0);
    int bin = 0;
    long double x_Hfac = 0.0;
    long double y_Hfac = 0.0;
    long double z_Hfac = 0.0;
    long double x_lambda = 0.0;
    long double y_lambda = 0.0;
    long double z_lambda = 0.0;
    long double calc_lambda = 1/ sqrt(conduction_electrons);
    #pragma omp parallel for private(array_index,array_index_y,array_index_z, x_pos,y_pos,z_pos, x_vel,y_vel,z_vel \
    ) reduction(+:x_lambda,y_lambda,z_lambda) schedule(dynamic) num_threads(3)
    for (int e = 0; e < conduction_electrons; e++) {
        array_index = 3*e;
        array_index_y = array_index + 1;
        array_index_z = array_index + 2;

        x_pos = new_electron_position[array_index];
        y_pos = new_electron_position[array_index_y]; 
        z_pos = new_electron_position[array_index_z];

        x_lambda += cos(4*M_PI * x_pos / lattice_constant);
        y_lambda += cos(4*M_PI * y_pos / lattice_constant);
        z_lambda += cos(4*M_PI * z_pos / lattice_constant);
        

        x_vel = 1e-10*new_electron_velocity[array_index];
        y_vel = 1e-10*new_electron_velocity[array_index_y];
        z_vel = 1e-10*new_electron_velocity[array_index_z];
        velocity_length = sqrt((x_vel*x_vel) + (y_vel*y_vel) + (z_vel*z_vel));

        #pragma omp critical
        {
            electron_position_output_down << "H" << ", " << x_pos << ", " << y_pos << ", " << z_pos << "\n";
            electron_velocity_output      << e   << ", " << x_vel << ", " << y_vel << ", " << z_vel << ", " << velocity_length << "\n";
        }
/*
        bin = int(floor(x_vel / width));
        //if(bin > histogram.size()/4) bin = histogram.size()/4;
        histogram[bin]++;
        if (bin == 0) std::cout << velocity_length_hist[array_index] << std::endl;

        bin = int(floor(y_vel / width));
        //if(bin > histogram.size()/4) bin = histogram.size()/4 + 1;
	    histogram[bin + 1]++;
        if (bin == 0) std::cout << velocity_length_hist[array_index_y] << std::endl;

        bin = int(floor(z_vel / width));
        //if(bin > histogram.size()/4) bin = histogram.size()/4 + 2;
	    histogram[bin + 2]++;
        if (bin == 0) std::cout << velocity_length_hist[array_index_z] << std::endl;

        bin = int(floor( velocity_length / width));
        //if(bin > histogram.size()/4) bin = histogram.size()/4 + 3;
	    histogram[bin + 3]++;
        if (bin == 0) std::cout << velocity_length_hist[array_index+3] << std::endl;

        velocity_length_hist[array_index]   = 0.0;
        velocity_length_hist[array_index_y] = 0.0;
        velocity_length_hist[array_index_z] = 0.0;
        velocity_length_hist[array_index+3] = 0.0;
       // if (e==0) std::cout << "v_f" << v_f << " velocity " << velocity_length << std::endl; */
    }
    lambda = 0.333333333333 * (x_lambda + y_lambda + z_lambda) / (CASTLE_output_rate * conduction_electrons);
    std::cout << "  " << current_time_step / total_time_steps * 100 << "%" << "\n";
    //  electron_position_output_up.close();
    electron_position_output_down.close();
    electron_velocity_output.close();
  /*  std::string time_stamp = std::to_string(current_time_step);
    electron_velocity_output.open("CASTLE/Electron_Velocity/" + time_stamp + "_Hfactor.csv");
    int x_binval,y_binval,z_binval = 0;
    for (int n=0; n < num_bins; n++) {
        x_binval = histogram[n*4];
        y_binval = histogram[n*4 + 1];
        z_binval = histogram[n*4 + 2];
        if (x_binval > 0) x_Hfac += x_binval * log(x_binval / conduction_electrons);
        if (y_binval > 0) y_Hfac += y_binval * log(y_binval / conduction_electrons);
        if (z_binval > 0) z_Hfac += z_binval * log(z_binval / conduction_electrons);
           // std::cout << x_binval << y_binval << z_binval << "\n";
        histogram[n*4] = 0;
        histogram[n*4 + 1] = 0;
        histogram[4*n + 2] = 0;
        histogram[4*n + 3] = 0;
        
        
           // std::cout << x_Hfac << y_Hfac << z_Hfac << "\n";
        electron_velocity_output << min + width*n << ",  " << x_binval  << ", " << y_binval  << ", " << z_binval << ", " << histogram[n*4 + 3]  << "\n";
    }
    electron_velocity_output.close();
    Hfac = 0.333333333 * (x_Hfac + y_Hfac + z_Hfac) * width / (conduction_electrons * CASTLE_output_rate); */
    double j = x_flux * constants::e * 1e20 / (1600 * CASTLE_output_rate * dt); //current density
    double nu = j / (n_f * constants::e); //drift velocity
    double I = n_f * 1600 * 1e-20 * nu * constants::e; //current
    
    if (current_time_step > 0) {
    mean_data << CASTLE_real_time << ", " << current_time_step << ", " \
        << MKE * 1e-20 * constants::m_e / CASTLE_output_rate << ", " \
        << MPE * 1e10 * constants::K / CASTLE_output_rate << ", "    \
        << ((MPE*constants::K*1e10) + (MKE*1e-20 * constants::m_e)) / CASTLE_output_rate << ", " \
        << -1* calc_lambda << ", " << calc_lambda << ", " << lambda << ", " \
        << j <<  ", " << nu << ", " << I << ", " << 4000 *1e-10/ I \
        << std::endl;
       // << Hfac << ", ?" << 
    }   
    long double mean_vel = sqrt(2*MKE) / conduction_electrons; //Still Angstroms
   // std::cout << mean_vel*dt << std::endl;
    if (mean_vel * dt > 0.0005) {
        terminaltextcolor(RED); 
        std::cout << "0.005 " << mean_vel << ", " << dt << ", " << dt*1e-1 << std::endl;
        dt = dt * 1e-1;
        terminaltextcolor(WHITE);
    }
    else if (mean_vel * dt > 0.00001) {
        terminaltextcolor(YELLOW);
        std::cout << "0.0001 " << mean_vel << ", " << dt << ", " << dt*1e-1 << std::endl;
        terminaltextcolor(WHITE);
    }
   /* bool one_time_slowdown = false;
    if (sqrt(MKE*2) > (10/dt)) {
        if (!one_time_slowdown) {
            dt = dt * 1e-2;
            one_time_slowdown = true;
            std::cout << sqrt(MKE*2) << " KE exceeding time-step of " << 1e-1 * dt << " Adjusting parameter to: " << dt << std::endl;
        }
    } */
  //  else if (sqrt(MKE*2) > 5/dt) std::cout << "KE begining to exceed time-step..." << 5/dt << std::endl;

    MKE = MPE = 0.0;
    
     //   electron_spin_output.close();
    CASTLE_output_data = false;
}




} //end of CASTLE namespace



