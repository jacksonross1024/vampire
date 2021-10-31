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
        omp_set_dynamic(0);
        omp_set_num_threads(8);
        std::cout << "CASTLE build time[s]: " << castle_watch.elapsed_seconds() << std::endl;
        #pragma omp parallel 
        {
            #pragma omp critical
            std::cout << "OpenMP capability detected. Parallelizing integration. Thread " << omp_get_thread_num() <<  " of threads: " << omp_get_num_threads() << std::endl;
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
    CASTLE_output_rate = 50;
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
  //  namespace fs = std::filesystem; {
      
     // filesystem::remove_all("CASTLE/Electron_Position");
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
    CASTLE_output_rate = 50; //sim::CASTLE_output_rate;
    dt = 1e-4; //reducd seconds (e10 scale factor), femptoSeconds
    temperature = 300; //sim::temperature;
    total_time_steps = sim::equilibration_time; //100
    x_flux = 0;
    y_flux = 0;
    z_flux = 0;
    current_time_step = 0;
    CASTLE_real_time = 0;
    
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

    initialize_velocity();

    std::cout << "E_f(AJ): " << E_f << ", TE(AJ): " << E_f*conduction_electrons << ", KE(AJ):" << TKE*1e10*constants::m_e/2 << ", PE " << TPE*1e10*constants::K/2 << ", PE+KE " <<  ((TKE*constants::m_e) + (TPE*constants::K))*1e10/2 << std::endl;
    mean_data.open("CASTLE/mean_data.csv");
    mean_data << "time, step, mean-KE, KE, mean-PE, PE, mean-TE, TE, -lambda, +lambda, mean-lambda, mean-x_flux, mean-y_flux, mean-z_flux, current, resistance" << "\n";

}

//====================================
// Creates and outputs atomic lattice
//      Currently static lattice
//====================================
void initialize_lattice() {
    
    atomic_size = 2; // sim::atomic_size; //Angst diameter. 
    screening_depth = 1; //sim::screening_depth; //Angstroms 
    lattice_atoms = 20 * 20 * 20; //Better lattice creation will come from VAMPIRE in future
    long double omega = 3000 * pow(6*M_PI*M_PI*n_f, 0.333333333333333333333);
    long double harmonic_constant = constants::h * omega;
    lattice_height = 40; //A
    lattice_width  = 40; // A
    lattice_depth  = 40; // A

   
    lattice_output.open("CASTLE/CASTLE_Lattice.xyz");

     // output lattice atoms and locations
    lattice_output << lattice_atoms << "\n"; //xyz file requires first-line for number of elements
    lattice_output << " static lattice" "\n"; //comment line
    
    atom_position.resize(lattice_atoms * 3, 0.0);
    atomic_phonon_energy.resize(lattice_atoms*3, 0);
    int array_index = 0; //local loop index variable
    for (int a = 0; a < lattice_atoms; ++a) {  
        array_index = 3*a;
        atomic_phonon_energy[array_index]   = harmonic_constant;
        atomic_phonon_energy[array_index+1] = harmonic_constant;
        atomic_phonon_energy[array_index+2] = harmonic_constant;

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

    
     
            if (err::check) std::cout << "Lattice output file and electron position file opened..." << std::endl;

    //=========
    // Initialize arrays for holding electron variables
    //      Arrays in super array format to take advantage of caching
    //========
    electron_position.resize(conduction_electrons * 3, 0); // ""'Memory is cheap. Time is expensive' -Steve Jobs; probably" -Michael Scott." -Headcannon.
    new_electron_position.resize(conduction_electrons * 3);
    // lattice_electrons.resize((total_electrons - conduction_electrons), 0.0);
    electron_velocity.resize(conduction_electrons * 3, 0); //Angstroms
    new_electron_velocity.resize(conduction_electrons * 3); //Angstroms
    electron_force.resize(conduction_electrons * 3, 0); //current and future arrays
    new_force_array.resize(conduction_electrons * 3);
    electron_potential.resize(conduction_electrons, 0);
    charge_distrib.resize(101, 0);
    mean_radius.resize(conduction_electrons*2);
    //mean_data_array.resize(total_time_steps*5 + 5, 0.0);
    // conduction_electron_spin.resize(conduction_electrons, false);
    // lattice_electron_spin.resize((total_electrons - conduction_electrons), false);
    // symmetry_list.resize(conduction_electrons);
   // velocity_length_hist.resize(conduction_electrons*4, 0.0);
     //   std::cout << "conduction electrons" << conduction_electrons << std::endl;
    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<long double> Theta_pos_distrib(0,2);
    std::uniform_real_distribution<long double> Phi_pos_distrib(0,1);
    std::normal_distribution<long double> radius_mod(1, 0.1);
   

    n_f = 1e10*1e20 * conduction_electrons / (lattice_width * lattice_height * lattice_depth); // e- / m**3
    E_f = constants::h * constants::h * powl(3 * M_PI * M_PI * n_f, 0.66666666666666666666666667) / (8 * M_PI * M_PI * constants::m_e); //Fermi-energy // meters
    mu_f = 5 * E_f / (3 * conduction_electrons);//Fermi-level //meters
    v_f = sqrt(2 * E_f / constants::m_e); //meters

    TPE = 0;
    TKE = 0;

    
    //TKE = 0; //total Kinetic energy, meters
    //  total_spin_up = 0;
   // total_spin_down = 0;
   
    long double phi,theta, x_pos,y_pos,z_pos, velocity_length = 0;
    //super loop for each conducting electron
    int array_index = 0;
    nearest_neighbor_list.resize(conduction_electrons);

    
    
    if (err::check) std::cout << "Prepare to set position: " << std::endl;
    for (int e = 0; e < conduction_electrons; e++) {
        nearest_neighbor_list[e].resize(conduction_electrons * 0.4, -1);
        array_index = 3*e;
        //random program for position initialization

        theta = M_PI*Theta_pos_distrib(gen);
        phi = M_PI*Phi_pos_distrib(gen);


            

        //initialize and output electron posititons
        x_pos = 1+ (atomic_size * (e % 20)) + (cosl(theta)*sinl(phi)*screening_depth*radius_mod(gen)); //Angstroms
        y_pos = 1+ (atomic_size * ((int(floor(e / 20))) % 20)) + (sinl(theta)*sinl(phi)*screening_depth*radius_mod(gen)); //Sets on radius of screening depth from nucleus
        z_pos = 1+ (atomic_size * floor(e/ 400)) + (cosl(phi)*screening_depth)*radius_mod(gen);

        if (x_pos < 0) x_pos += 40;
        else if (x_pos > 40) x_pos -= 40;

        if (y_pos < 0) y_pos += 40;
        else if (y_pos > 40) y_pos -= 40;

        if (z_pos < 0) z_pos += 40;
        else if (z_pos > 40) z_pos -= 40;

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
           
        electron_position_output_down << "H" << "    " << x_pos << "    " << y_pos << "    " << z_pos << "\n";    
            
              //std::random_device rd;  //Will be used to obtain a seed for the random number engine
            //std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd(
                    //  if (err::check) std::cout << "Random angles chosen. v_f: %d. cos(theta): %d. sin(theta): %d. sin(0.5 - phi): %d" << v_f << cos(theta) << sin(theta) << sin(0.5 - phi) << std::endl;
        }
   // std::cout << "X: " << X << "Y: " << Y << " Z: " << Z << " Px: " << Px << " Py: " << Py << " Pz: " << Pz << std::endl;
    
    //electron_position_output_up.close();
    electron_position_output_down.close(); 
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

    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    //std::uniform_int_distribution<> random_electron(0,8000);
   // chosen_electron = random_electron(gen);
    int array_index, array_index_i;
    long double x,y,z, d_x,d_y,d_z, x_mod,y_mod,z_mod, x_distance,y_distance,z_distance, length, force;
    long double theta,phi, PE;
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

        d_x = x - ((atomic_size * round((x-1) / atomic_size)) + 1); //closest x atom index
        d_y = y - ((atomic_size * round((y-1) / atomic_size)) + 1); //closest y atom index
        d_z = z - ((atomic_size * round((z-1) / atomic_size)) + 1); //closest z atom index
        length = sqrt((d_x*d_x)+(d_y*d_y)+(d_z*d_z));
      //  if(e==0) std::cout << 99*((405* exp(-15* length)) - (exp(-1 * length)))/4 << std::endl;
      int count = 0;
        std::cout << "d_r: " << sqrt((d_x*d_x)+(d_y*d_y)+(d_z*d_z)) << std::endl;
            //cube around electron
        long double negative_PE = 0;
        for (int a = 0; a < 1331; a++) {
            
            x_mod = (atomic_size * (a % 11)) - (5*atomic_size);
            y_mod = (atomic_size * ((int(floor(a/11))) % 11)) - (5*atomic_size);
            z_mod = (atomic_size * floor(a / 121)) - (5*atomic_size);
            x_distance = x_mod - d_x;
            y_distance = y_mod - d_y;
            z_distance = z_mod - d_z;

            length = sqrt((x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance));
            if (length < 0.00001) length = 0.00001;
            if(length > 10) continue;
            count++;
            // force
            force = -28*((3.3*3.3 * exp(-3.3 * length)) - (exp(-1 * length)));
           
            PE    =  28*((3.3* exp(-3.3* length)) - (exp(-1 * length)));
            TPE += PE;
            electron_potential[e] += PE;
            if (e == chosen_electron){
                if (x_distance > 0) {
                    int bin = int(floor(x_distance *10));
                  //  charge_distrib[bin] += -1*PE*x_distance / length;
                }
            }
         //   if(e == 0) std::cout << "force " << force << std::endl;
            phi   = acosl(z_distance / length);
            theta = atanl(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

             negative_PE += PE;
            electron_force[array_index]     += force * cos(theta)*sin(phi);
            electron_force[array_index + 1] += force * sin(theta)*sin(phi);
            electron_force[array_index + 2] += force * cos(phi);

        }
       // if (e==0) std::cout << "e-a count: " << count << std::endl;
             if (err::check) if(e==0) std::cout << "Nearest atoms set..." << std::endl;


    //set e-e repulsion local variables
       //     std::cout << "conduction electrons" << conduction_electrons << std::endl;

    
    
       //set e-e repulsion
    
        std::cout << "e-a attraction: " << negative_PE;
      //  electron_spin = conduction_electron_spin[e];
                if (err::check) if(e ==0) std::cout << "Calculating conduction electron repulsion" << std::endl;
        int neighbor_count = 0;
        count = 0;
        long double positive_PE = 0;
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

            if (length > 250) continue;
            
            nearest_neighbor_list[e][neighbor_count] = i;
            neighbor_count++;
            if (length > 100) continue;
            count++;
            length = sqrt(length);
            //if (length < 0.000001) length = 0.000001;

          /*  if (electron_spin == electron_spin_two) {
                force = 1e10 * ((constants::K / (length * length * constants::m_e))- (635 * constants::kB / (length * constants::m_e)));
               // std::cout << " 1, 0 " << 1e10 * (constants::K / (length * length * constants::m_e)) << " delta1,1 " << (635 * constants::kB / (length * constants::m_e)) << std::endl;
            }
            else { */
            force = 1 / (length * length);
            
            PE = force * length;
            electron_potential[e] += PE;
            TPE += PE;
            if (e == chosen_electron){
                if (x_distance < 0) {
                    int bin = int(floorf(abs(x_distance) *10));
                  //  charge_distrib[bin] += PE * x_distance / length;
                }
            }
           // if(e == 0) std::cout << "force " << force << std::endl;
            //    if (err::check) std::cout << "Distances grabbed" << std::endl; 

            phi   = acosl(z_distance / length);
            theta = atanl(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            positive_PE += PE;
            electron_force[array_index]     += force * cos(theta)*sin(phi);
            electron_force[array_index + 1] += force * sin(theta)*sin(phi);
            electron_force[array_index + 2] += force * cos(phi);
         /*   //make use of symmetry
            electron_force[array_index_i]     += -1 * x_unit;
            electron_force[array_index_i + 1] += -1 * y_unit;
            electron_force[array_index_i + 2] += -1 * z_unit;
            symmetry_list[i][e] = true; //apply symmetry flag
        */

            //  if (e ==0)   std::cout << electron_force[array_index] << std::endl;
        } 
       // if(e==0) std::cout << "e-e count: " << count << std::endl;
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
   /*     if (e == chosen_electron) {
            std::ofstream charge_distribution_output;
            charge_distribution_output.open("CASTLE/Electron_Charge/init.csv");
            charge_distribution_output << chosen_electron << "\n";
            charge_distribution_output  << "Initial charge distribution along x-vector." << "\n";

            for (int c = 0; c < 100; c++) {
                charge_distribution_output << c/10 << ", " << charge_distrib[c] << "\n";
                charge_distrib[c] = 0;
            }
            charge_distribution_output.close();
        } */
         std::cout << ". e-e repulsion: " << positive_PE << ". net force: " << positive_PE + negative_PE << std::endl;
    }

}

void initialize_velocity() {
     
    electron_velocity_output.open("CASTLE/Electron_Velocity/init.csv");
    electron_velocity_output << "electron number, x-component, y-component, z-component, length" << std::endl;  

    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<long double> Theta_Vel_distrib(0,2);
    std::uniform_real_distribution<long double> Phi_Vel_distrib(0,1);

    long double phi,theta, vel; //A/fS
    int array_index;
    long double TE,PE;
    for(int e = 0; e < conduction_electrons; e++) {
        array_index = 3*e;
        theta = M_PI * Theta_Vel_distrib(gen);
        phi = M_PI * Phi_Vel_distrib(gen);

      //vel = sqrt(2*KE/m_e) =               TE     -     PE
       // if (e ==0 ) std::cout << electron_potential[e] << std::endl;
        vel = 1e-15*sqrt(abs(2* ((E_f*1e20 - ((electron_potential[e]*constants::K_A)/2))/constants::m_e))); // m/s -> Angstroms / s -> A/fs = 1e-5
      //  if(e == 0) std::cout << "KE: " << 0.5*constants::m_e*v_f*v_f << ", PE: " << electron_potential[e]*1e10*constants::K << std::endl;
        electron_velocity[array_index]     = cosl(theta)*sinl(phi)*vel; 
        electron_velocity[array_index + 1] = sinl(theta)*sinl(phi)*vel;
        electron_velocity[array_index + 2] = cosl(phi)*vel;
        TKE += vel*vel;
                  // std::cout << "v_f" << v_f << " velocity " << velocity_length * 1e-10 << std::endl;
        electron_velocity_output << e << ", " << 1e5*electron_velocity[array_index] << " , " << 1e5*electron_velocity[array_index + 1] << " , " << 1e5*electron_velocity[array_index + 2] << " , " << 1e5*vel << std::endl;
       // if(e==100) std::cout << e << ", " << 1e-10*electron_velocity[array_index] << " , " << 1e-10*electron_velocity[array_index + 1] << " , " << 1e-10*electron_velocity[array_index + 2] << " , " << 1e-10*vel << std::endl;
                    // if (err::check) std::cout << "Velocities randomized..." << std::endl;
            //   electron_spin_output << conduction_electrons << "  " << conduction_electron_spin[conduction_electrons - 1] << std::endl;
    }

    electron_velocity_output.close();

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
    
    int array_index, array_index_y, array_index_z;
    long double x_pos, y_pos, z_pos;
    long double x_vel, y_vel ,z_vel, velocity_length, lambda;
  //  int num_bins = 1e6;
  //  int max = 1e7;
   // int min = 0;
    int lattice_constant = 2;
   // double width = (max - min) / num_bins;
  //  std::vector <int> histogram;
   // histogram.resize(num_bins*4, 0);
    //int bin = 0;
   // long double x_Hfac = 0.0;
   // long double y_Hfac = 0.0;
   // long double z_Hfac = 0.0;
    long double x_lambda = 0.0;
    long double y_lambda = 0.0;
    long double z_lambda = 0.0;
    long double calc_lambda = 1/ sqrt(conduction_electrons);
 /*   long double x_mo = 0.0;
    long double y_mo = 0.0;
    long double z_mo = 0.0; */
  //  #pragma omp parallel for private(velocity_length,array_index,array_index_y,array_index_z, x_pos,y_pos,z_pos, x_vel,y_vel,z_vel \
    ) reduction(+:x_lambda,y_lambda,z_lambda) schedule(dynamic) num_threads(1)
    for (int e = 0; e < conduction_electrons; e++) {
        array_index   = 3*e;
        array_index_y = array_index + 1;
        array_index_z = array_index + 2;

        x_pos = new_electron_position[array_index];
        y_pos = new_electron_position[array_index_y]; 
        z_pos = new_electron_position[array_index_z];

        x_lambda += cos(4*M_PI * x_pos / lattice_constant);
        y_lambda += cos(4*M_PI * y_pos / lattice_constant);
        z_lambda += cos(4*M_PI * z_pos / lattice_constant);
        

        x_vel = 1e5*new_electron_velocity[array_index];
        y_vel = 1e5*new_electron_velocity[array_index_y];
        z_vel = 1e5*new_electron_velocity[array_index_z];
        velocity_length = sqrt((x_vel*x_vel) + (y_vel*y_vel) + (z_vel*z_vel));

      //  x_mo += x_vel;
        //y_mo += y_vel;
        //z_mo += z_vel;
       // #pragma omp critical
        {
            electron_position_output_down << "H" << ", " << x_pos << ", " << y_pos << ", " << z_pos << std::endl; //<< ", " << mean_radius[2*e] << ", " << mean_radius[2*e+1] << "\n";
            electron_velocity_output      << e   << ", " << x_vel << ", " << y_vel << ", " << z_vel << ", " << velocity_length << "\n";
           // if(e==100) std::cout << e   << ", " << x_vel << ", " << y_vel << ", " << z_vel << ", " << velocity_length << std::endl;
        }
  /*      if (e == chosen_electron) {
            //std::cout << omp_get_num_threads() << std::endl;
            std::ofstream charge_distribution_output;
            charge_distribution_output.precision(10);
            charge_distribution_output << std::scientific;
            std::string time_stamp = std::to_string(current_time_step);
            charge_distribution_output.open("CASTLE/Electron_Charge/" + time_stamp + ".csv");
            charge_distribution_output << chosen_electron << "\n";
            charge_distribution_output << time_stamp << ". Charge distribution along x-vector." << "\n";

            for (int c = 0; c < 100; c++) {
                charge_distribution_output << 0.1*c << ", " << charge_distrib[c] << "\n";
                charge_distrib[c] = 0.0;
            }
            charge_distribution_output.close();
        } */ 

/*        bin = int(floor(x_vel / width));
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
    lambda = (x_lambda + y_lambda + z_lambda) / (3*CASTLE_output_rate * conduction_electrons);
    std::cout << "  " << current_time_step / total_time_steps * 100 << "%. " << std::endl; 
   // if (current_time_step % (CASTLE_output_rate*10) == 0) std::cout << "Estimated time remaining: " << (total_time_steps - current_time_step) \n";
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
    long double j = x_flux * constants::e * 1e20 / (1600 * CASTLE_output_rate * dt); //current density
    long double nu = j / (n_f * constants::e); //drift velocity
    long double I = n_f * 1600 * 1e-20 * nu * constants::e; //current
    
    if (current_time_step > 0) {

    mean_data << CASTLE_real_time << ", " << current_time_step << ", " \
        << MKE * 1e10 * constants::m_e / (2*CASTLE_output_rate) << ", " << TKE * 1e10 * constants::m_e/2 << ", " \
        << MPE * 1e10 * constants::K / (2*CASTLE_output_rate)<< ", " << TPE * 1e10 * constants::K/2 << ", " \
        << ((MPE*constants::K*1e10) + (MKE*1e10 * constants::m_e)) / (2*CASTLE_output_rate) << ", " << ((TKE * 1e10 * constants::m_e) + (TPE * 1e10 * constants::K))/2 << ", " \
        << -1* calc_lambda << ", " << calc_lambda << ", " << lambda << ", " \
        << chosen_electron / (CASTLE_output_rate*dt) << ", " << x_flux / CASTLE_output_rate<< ", " << y_flux / CASTLE_output_rate << ", " << z_flux / CASTLE_output_rate << ", " \
        << I << ", " << 4 / I \
        << std::endl;
       // << Hfac << ", ?" << 
    }   
    long double mean_vel = sqrt(MKE) / (CASTLE_output_rate *conduction_electrons); //Still Angstroms
   // std::cout << mean_vel*dt << std::endl;
    if (mean_vel * dt > 0.01) {
        terminaltextcolor(RED); 
        std::cout << "0.01 " << mean_vel << ", " << dt << ", " << dt*1e-1 << std::endl;
        dt = dt * 1e-1;
        terminaltextcolor(WHITE);
    } /*
    else if (mean_vel * dt > 0.0005) {
        terminaltextcolor(YELLOW);
        std::cout << "0.00001 " << mean_vel << ", " << dt << ", " << dt*1e-1 << std::endl;
        terminaltextcolor(WHITE);
    } */
   /* bool one_time_slowdown = false;
    if (sqrt(MKE*2) > (10/dt)) {
        if (!one_time_slowdown) {
            dt = dt * 1e-2;
            one_time_slowdown = true;
            std::cout << sqrt(MKE*2) << " KE exceeding time-step of " << 1e-1 * dt << " Adjusting parameter to: " << dt << std::endl;
        }
    } */
  //  else if (sqrt(MKE*2) > 5/dt) std::cout << "KE begining to exceed time-step..." << 5/dt << std::endl;

    MKE = 0;
    MPE = 0;
    x_flux = 0;
    y_flux = 0;
    z_flux = 0;
    chosen_electron = 0;
     //   electron_spin_output.close();
    CASTLE_output_data = false;
}




} //end of CASTLE namespace



