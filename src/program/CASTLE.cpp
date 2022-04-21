// ==============================================
//   Coupled Atomistic and Spintronic Thermal Lattice Ensemble
//
//  =========       ========      ========   ============   ||           =========
// ||             ||        ||   ||               ||        ||          ||
// ||             ||        ||   ||               ||        ||          ||
// ||             ||        ||   ||               ||        ||          ||
// ||             || ====== ||   =========        ||        ||           =========
// ||             ||        ||           ||       ||        ||          ||
// ||             ||        ||           ||       ||        ||          ||
// ||             ||        ||           ||       ||        ||          ||
//  =========                     ========                   =========   =========
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
       // omp_set_dynamic(1);
        //omp_set_num_threads(8);
        std::cout << "CASTLE build time[s]: " << castle_watch.elapsed_seconds() << std::endl;
        #pragma omp parallel 
            #pragma omp critical
             std::cout << "OpenMP capability detected. Parallelizing integration. Thread " << omp_get_thread_num() <<  " of threads: " << omp_get_num_threads() - 1 << std::endl;
        
        std::cout << "Storming CASTLE..." << std::endl;
   
 /*   double x,y,z,r,r_min, d_x,d_y,d_z, p_x,p_y,p_z, d_p_x,d_p_y,d_p_z, p, force, f_x,f_y,f_z, d_f_x,d_f_y,d_f_z, theta, phi, potential;
    double a_x,a_y, a_d_x,a_d_y,a_d_z, a_p_x,a_p_y, a_f_x,a_f_y, a_d_f_x,a_d_f_y;
    std::ofstream electron_output;
    std::ofstream ballistic_data;
    std::ofstream atom_output;
  
    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> scattering_chance (0,1);

    //dt = 1e-4;

    a_x = 0;
    a_y = 0;

    x = -7;
    y = 0;
    d_z = 0;

    p_x = -1e-5*v_f; 
    p_y = 0;

    a_p_x = 0;
    a_p_y = 0;
    a_d_z = 0;
ballistic_data.open("CASTLE/Ballistic_Data");
for(int a = 0; a < 1500; a++) {
    if(a == 750) continue;
      bool collision = false;
    r = 7;
    theta = M_PI*(1.25 - (a * 3.3333333333e-4));
   
   // ((a/3000) - 0.25)*M_PI;
    x = r*cos(theta);
    y = r*sin(theta);
    int t = 0;
   // std::cout << theta << ", " <<  x << ", " << y << std::endl;

    std::string height  = std::to_string(a);
    electron_output.open("Castle/Ballistic_Electron/" + height + ".xyz");
    electron_output << "Ballistic Electron" << "\n" << "\n";
    atom_output.open("CASTLE/Ballistic_Atom/" + height + ".xyz");
    atom_output << "Ballistic Atom" << "\n" << "\n";

    p_x = 3e-5*v_f; 
    p_y = 0;

    r_min = r = sqrt((x*x)+(y*y));
  //  force = 10*((3.3*3.3*exp(-3.3*r)) - 1 / (r*r));
    force = (1 - (5*exp(-4*r)*(4*r + 1))) / (r * r);
    
    
    theta = atanl(y/x);
    if (x < 0) theta += M_PI;

    ballistic_data << theta << ", " << theta / M_PI << ", ";
    f_x = force * cos(theta);
    f_y = force * sin(theta);

    a_f_x = -1*f_x;
    a_f_y = -1*f_y;

    a_p_x = 0;
    a_p_y = 0;
    a_x = 0;
    a_y = 0;

   // std::cout << force<< ", " << f_x << ", " << f_y << std::endl; 
   // std::cout << "time_step: " << t << ", x_distance: " << x << ", y_distance: " << y << "\n";
    while(r <= 7.01) {
        
       // std::string time = std::to_string(t);
        // new output file
        
        //updte position
        d_x = x + (p_x*dt) + (f_x*dt*dt*constants::K_A / (2*constants::m_e_r));
        d_y = y + (p_y*dt) + (f_y*dt*dt*constants::K_A / (2*constants::m_e_r));
      //  d_z = z + (p_z*dt) + (f_z*dt*dt*constants::K_A / (2*constants::m_e_r));

        a_d_x = a_x + (a_p_x*dt) + (a_f_x*dt*dt*constants::K_A / (2*atomic_mass));
        a_d_y = a_y + (a_p_y*dt) + (a_f_y*dt*dt*constants::K_A / (2*atomic_mass));

        if (t % 10 == 0) electron_output << t << ", " << d_x << ", " << d_y << ", " << d_z << "\n";
        if (t % 10 == 0) atom_output << t << ", " << a_d_x << ", " << a_d_y << ", " << a_d_z << "\n";

        //update forces
        r = sqrt(((d_x - a_d_x)*(d_x - a_d_x)) + ((d_y - a_d_y)*(d_y - a_d_y)));
        if(r < r_min) r_min = r;
        force = -1.06*(1 - (27*exp(-5*r)*(5*r + 1))) / (r * r);
        //potential = 12*((3.3*exp(-1*r)) - exp(-1*r));

     //   phi   = acos(d_z / r);
        theta = atanl(d_y / d_x);
        if(d_x < 0) theta += M_PI;
    

        d_f_x = force * cos(theta);
        d_f_y = force * sin(theta);

        a_d_f_x = -1*d_f_x;
        a_d_f_y = -1*d_f_y;

       // a_d_f_x += -2e7*a_d_x;
        //a_d_f_y += -2e7*a_d_y;

      //  d_f_z = force * cos(phi); 
     /*   if(r < 0.7 && !collision) {
            double excitation_energy = ((p_x*p_x)+(p_y*p_y))*constants::m_e_r*0.5 + 1.06*((27*expl(-5*r) / r) - (1 / r)) - E_f_A;
            double scattering_velocity = sqrt((p_x*p_x)+(p_y*p_y)) - sqrt(2*E_f_A/constants::m_e_r);
            if(excitation_energy > 0 && scattering_velocity > 0) {
                if(scattering_chance(gen) < (1 - exp(-1*excitation_energy))) {
                double vel = sqrt((p_x*p_x)+(p_y*p_y));
                theta = atanl(p_y / p_x);
                if(p_x < 0) theta += M_PI;
                
                
                p_x = scattering_velocity * cos(theta);
                p_y = scattering_velocity * sin(theta);
                

                vel = sqrt((a_p_x*a_p_x)+(a_p_y*a_p_y));
                theta = atanl(a_p_y / a_p_x);
                if(a_p_x < 0) theta += M_PI;
                scattering_velocity += sqrt(2*E_f_A/atomic_mass);
              
                a_p_x  = scattering_velocity * cos(theta);
                a_p_y = scattering_velocity * sin(theta);
                collision = true;
                }
            }
                
        } 
        //update velocity
        p_x += ((d_f_x + f_x)*dt*constants::K_A / (2*constants::m_e_r));
        p_y += ((d_f_y + f_y)*dt*constants::K_A / (2*constants::m_e_r));

        a_p_x += ((a_d_f_x + a_f_x)*dt*constants::K_A / (2*atomic_mass));
        a_p_y += ((a_d_f_y + a_f_y)*dt*constants::K_A / (2*atomic_mass));
       // d_p_z = p_z + ((d_f_z + f_z)*dt*constants::K_A / (2*constants::m_e_r));

        x = d_x;
        y = d_y;
        a_x = a_d_x;
        a_y = a_d_y;

        f_x = d_f_x;
        f_y = d_f_y;
        a_f_x = a_d_f_x;
        a_f_y = a_d_f_y;

        t++;

        phi = atanl(p_y/p_x);
        if(p_x < 0) phi += M_PI;
    }
    electron_output.close();
    atom_output.close();
   // std::cout << "time_step: " << t << ", x_distance: " << d_x << ", y_distance: " << d_y << "\n";
    ballistic_data << theta / M_PI << ", " << r_min << ", " << phi / M_PI << "\n";
} ballistic_data.close();  */
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
   // CASTLE_output_rate = 100;
    total_time_steps = sim::loop_time;
   
    //========
    // Run averaging step
    //========
   sim::integrate(total_time_steps);

    std::cout << "Averaging complete. " << castle_watch.elapsed_seconds() << " [s] elapsed. Step-time: " << castle_watch.elapsed_seconds() / total_time_steps << std::endl;
    
    mean_data.close();
    temp_data.close();
}

//====================================
// Initialize function to set up CASTLE
//====================================
void initialize () {
       
            if (err::check) std::cout << "Initializing CASTLE..."  << std::endl;
 
      
     // filesystem::remove_all("CASTLE/Electron_Position");
    //========
       //=========
    // Grab simulation variables from VAMPIRE
    //=========
    conduction_electrons = atoms::num_atoms;
    lattice_atoms = atoms::num_atoms; //Better lattice creation will come from VAMPIRE in future
    CASTLE_output_rate = sim::partial_time;
    dt = mp::dt_SI * 1e15;//-4; //S -> femptoSeconds
    TTMe = TTMp = d_TTMe = d_TTMp = Tp = Te = sim::temperature;
    //Te = 1300.0;
    total_time_steps = sim::equilibration_time; //100
    CASTLE_MD_rate = sim::CASTLE_MD_rate;
  
    applied_voltage_sim = sim::applied_voltage_sim;
    heat_pulse_sim = sim::heat_pulse_sim;
     //V/m

    ee_coupling = sim::ee_coupling;
    ea_coupling = sim::ea_coupling;

    //========
    // Initialzie variables used by all functions
    //========
    lattice_height = cs::system_dimensions[0] + x_unit_size; //A
    lattice_width  = cs::system_dimensions[1] + y_unit_size; // A
    lattice_depth  = cs::system_dimensions[2] + z_unit_size; // A

    x_flux = 0;
    y_flux = 0;
    z_flux = 0;
    current_time_step = 0;
    CASTLE_real_time = 0;
    uniform_random.seed(2137082040);
    int_random.seed(2137082040);

    mu_r = (atomic_mass + constants::m_e_r) / (atomic_mass * constants::m_e_r );
    combined_mass = 1 / (atomic_mass + constants::m_e_r);
    applied_voltage = sqrt(1.60218e-19*2.0*sim::applied_voltage/constants::eps_0/(1e-30*lattice_width * lattice_height * lattice_depth)); //eV -> V/m
    n_f = 1e30 * conduction_electrons / (lattice_width * lattice_height * lattice_depth); // e- / A**3 -> e-/m**3
    E_f = constants::h * constants::h * pow(3 * M_PI * M_PI * n_f, 0.66666666666666666666666667) / (8 * M_PI * M_PI * constants::m_e); //Fermi-energy // meters
    E_f_A = E_f*1e20; //Angstroms
    mu_f = 1e20*5 * E_f / (3 * conduction_electrons);//Fermi-level //Angstroms
    v_f = sqrt(2 * E_f / constants::m_e); //meters

    a_specific_heat = 25.0 /6.02e3; // //sim::TTCl; //J/K/mol -> [e20/Na] AJ/K/e-   
    a_specific_heat_i = 1.0 / a_specific_heat;
    a_heat_capacity = 1e-27*a_specific_heat * n_f; //AJ/K/particle [] AJ/K/nm**3
    a_heat_capacity_i = 1.0 / a_heat_capacity;

    e_specific_heat = constants::kB_r*3.0/10.0; // gamma; //AJ/K**2/e- 
    e_specific_heat_i = 1.0 / e_specific_heat;
    e_heat_capacity = 1e-27*e_specific_heat * n_f; //AJ/K**2/e- -> [e-/m**3] -> AJ/K**2/nm**3
    e_heat_capacity_i = 1.0 / e_heat_capacity;

    ea_coupling_strength = 1e-6*sim::ea_coupling_strength*constants::eV_to_AJ*constants::eV_to_AJ/(constants::hbar_r*constants::hbar_r); // meV**2 -> 1/fs**2
    phonon_energy = 1e-3*sqrt(sim::ea_coupling_strength/0.084)*constants::eV_to_AJ;

    atomic_mass = 58.69 * 1.6726219e3; // kg * 1e30 for reduced units
    power_density = 1e1*sim::heat_pulse; // mJ/cm**2 -> [e17/e14/e2(fs)] AJ/fs/nm**2
    
    const static double tau = 1.0 /(M_PI*constants::hbar_r*ea_coupling_strength); //AJ**-1 fs
    G = e_heat_capacity/tau; //J/s/K/m**3 [e20/e15/e27] AJ/fs/K/nm**3
    
    //G=Ce/Te-p = pihbar constant (meV**2)Ef*n_f*(1/eps)**2
    ea_rate = -1.0e-1*dt / tau;
    ee_rate = -1.0*dt*sim::ee_coupling_strength/(constants::eV_to_AJ*constants::eV_to_AJ); //fs**-1 AJ**-2

    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> seed_gen(0, int(pow(2,32)));
    for(int i = 0; i < omp_get_num_threads(); i++) {
      omp_uniform_random[i].seed(seed_gen(gen));
      omp_int_random[i].seed(seed_gen(gen));
    }

    // Initialize lattice
    //========
    initialize_positions();
            if (err::check) std::cout << "Lattice built " << std::endl;
    
    //========
    // initialize electrons: lattice and conduction bands, velocity, spin, etc.
    //=======
    initialize_forces();
            if (err::check) std::cout << "Electrons ready..." << std::endl;
    //=======
    // Calls forces set up
    //=======
            if (err::check) std::cout << "Forces ready..." << std::endl;

    initialize_velocities();

             if (err::check) std::cout << "Particles a movin" << std::endl;
  
    std::cout << "E_f(J): " << E_f << ", G(AJ/K/fs/nm**3): " << G << ", Ce@300K (AJ/K/nm**3): " << e_heat_capacity << ", Cl(AJ/K/nm**3): " << a_heat_capacity << ", Tau_ea@300K(fs): " << tau <<  ", Tau_ee: " << ee_rate << ", E_field(V/m): " << applied_voltage << std::endl;
    
    char directory [256];
    if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog." << std::endl;
    }
    
    temp_data.open(string(directory) + "/temp_data.csv");
    mean_data.open(string(directory) + "/mean_data.csv");
    mean_data << "time, step, mean-EKE, mean-LE, mean-Te, mean-Tp,  mean-radius, mean-e-a-collisions, mean-e-e-collisions, mean-x_flux, mean-y_flux, mean-z_flux" << "\n";

}

//====================================
// Creates and outputs atomic lattice
//      Currently static lattice
//====================================
void initialize_positions() {

    initialize_lattice();

     if (err::check)  std::cout << "Lattice" << std::endl;

    initialize_electrons();

     if (err::check)  std::cout << "Electrons" << std::endl;
}


void initialize_lattice() {
    
    atomic_size = 2; // sim::atomic_size; //Angst diameter. 
    screening_depth = 0.875; //sim::screening_depth; //Angstroms 

    e_a_neighbor_cutoff = 100.0;
    e_e_neighbor_cutoff = 100.0;
    e_a_coulomb_cutoff = 9.0;
    e_e_coulomb_cutoff = 9.0;
    a_a_neighbor_cutoff = 36.0;

    char directory [256];
    if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog." << std::endl;
    }
    lattice_output.open(string(directory) + "/CASTLE_Lattice.xyz");

    // output lattice atoms and locations
    lattice_output << lattice_atoms << "\n"; //xyz file requires first-line for number of elements
    lattice_output << " dynamic lattice" "\n"; //comment line
    
    atom_anchor_position.resize(lattice_atoms*3,0);
    atom_position.resize(lattice_atoms * 3, 0);
    // new_atom_position.resize(lattice_atoms * 3, 0);
    //atom_velocity.resize(lattice_atoms * 3,0);
    //new_atom_velocity.resize(lattice_atoms * 3,0);
    // atom_force.resize(lattice_atoms * 3,0);
    // new_atom_force.resize(lattice_atoms,0);
    atom_potential.resize(lattice_atoms,0);
    //new_atom_potential.resize(lattice_atoms,0);
    phonon_distribution.resize(2*lattice_atoms, 0.0);
    
    new_phonon_energy = 0.0;
    a_a_scattering_count = 0;
    
   // create_phonon_distribution();

    int array_index; //local loop index variable
    atomic_nearest_atom_list.resize(lattice_atoms);
    atomic_nearest_electron_list.resize(lattice_atoms);

    int a_density = 200 + int(round(pow(a_a_neighbor_cutoff, 1.5)*1.25*M_PI * n_f * 1e-30));
    int e_density = 200 + int(round(pow(e_a_neighbor_cutoff, 1.5)*1.25*M_PI * n_f * 1e-30));

    #pragma omp parallel for schedule(static) 
    for (int a = 0; a < lattice_atoms; ++a) {  
        atomic_nearest_atom_list[a].resize(a_density);
        atomic_nearest_electron_list[a].resize(e_density);
        array_index = 3*a;

        atom_anchor_position[array_index]     = atoms::x_coord_array[a] + 0.5*x_unit_size;
        atom_anchor_position[array_index + 1] = atoms::y_coord_array[a] + 0.5*y_unit_size;
        atom_anchor_position[array_index + 2] = atoms::z_coord_array[a] + 0.5*z_unit_size;
        
        atom_position[array_index]   = atom_anchor_position[array_index];
        atom_position[array_index+1] = atom_anchor_position[array_index+1];
        atom_position[array_index+2] = atom_anchor_position[array_index+2];

      //  atom_potential[a] = E_f_A;
    }

    for(int a = 0; a < lattice_atoms; a++) {
      array_index = a*3;
      lattice_output << "Ni" << "     " << atom_position[array_index] << "     " << atom_position[array_index + 1] << "   " << atom_position[array_index + 2] << "  " << a << "\n";  
    }
    lattice_output.close(); 

   // Tp = 300.0;
   
}

//====================================
// Function call to initialize basic electron structures
//      Will set up lattice and conduction electrons and velocites
//      Any additional features, e.g. spin, will need to be added here through function calls
//====================================./
void initialize_electrons() {

    char directory [256];
    if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog." << std::endl;
    }

    electron_position_output_down.open(string(directory) + "/Electron_Position_Init.xyz");
    electron_position_output_down << conduction_electrons << "\n";  
    electron_position_output_down << "Initial positions for electrons" << "\n";  

            if (err::check) std::cout << "Lattice output file and electron position file opened..." << std::endl;
    //=========
    // Initialize arrays for holding electron variables
    //      Arrays in super array format to take advantage of caching
    //========
    electron_position.resize(conduction_electrons * 3, 0); // ""'Memory is cheap. Time is expensive' -Steve Jobs; probably" -Michael Scott." -Headcannon.
    new_electron_position.resize(conduction_electrons * 3);
    electron_velocity.resize(conduction_electrons * 3, 0); //Angstroms
    //new_electron_velocity.resize(conduction_electrons * 3); //Angstroms
    //  electron_force.resize(conduction_electrons * 3, 0); //current and future arrays
    // new_electron_force.resize(conduction_electrons * 3);

    electron_potential.resize(conduction_electrons, 0);
    //new_electron_potential.resize(conduction_electrons, 0);
    mean_radius.resize(conduction_electrons*2,1);
  
    e_a_scattering_count = 0;
    e_e_scattering_count = 0;

    e_a_neighbor_cutoff = 100.0;
    e_e_neighbor_cutoff = 169.0;
    e_a_coulomb_cutoff = 9.0;
    e_e_coulomb_cutoff = 9.0;
    
    electron_nearest_electron_list.resize(conduction_electrons);
    electron_nearest_atom_list.resize(conduction_electrons);
    external_interaction_list.resize(conduction_electrons, false);
    electron_ee_scattering_list.resize(conduction_electrons);
    electron_ea_scattering_list.resize(conduction_electrons);
        if (err::check) std::cout << "Prepare to set position: " << std::endl;
    int e_density = 1100 + int(round(pow(e_e_neighbor_cutoff, 1.5)*1.25*M_PI * n_f * 1e-30));
    int ee_density = e_density;
    std::cout << e_density << ", " << (pow(e_e_neighbor_cutoff, 1.5)*1.25*M_PI * n_f * 1e-30) + 1100 << std::endl;
    #pragma omp parallel for schedule(static) 
    for (int e = 0; e < conduction_electrons; e++) {

        electron_nearest_electron_list[e].resize(e_density);
        electron_nearest_atom_list[e].resize(ee_density);
        electron_ee_scattering_list[e].resize(ee_density);
        electron_ea_scattering_list[e].resize(ee_density);
        int array_index = 3*e;

        double theta = 2.0*M_PI*omp_uniform_random[omp_get_thread_num()]();
        double phi = M_PI*omp_uniform_random[omp_get_thread_num()]();

        //initialize and output electron posititons
        double x_pos = atom_anchor_position[3*e]   + cos(theta)*sin(phi)*screening_depth;//*radius_mod(gen)); //Angstroms
        double y_pos = atom_anchor_position[3*e+1] + sin(theta)*sin(phi)*screening_depth;//*radius_mod(gen)); //Sets on radius of screening depth from nucleus
        double z_pos = atom_anchor_position[3*e+2] + cos(phi)*screening_depth;//*radius_mod(gen);

       // electron_potential[e] = E_f_A;
        
        if (x_pos < 0.0) x_pos += lattice_width;
        else if (x_pos > lattice_width) x_pos -= lattice_width;

        if (y_pos < 0.0) y_pos += lattice_depth;
        else if (y_pos > lattice_depth) y_pos -= lattice_depth;

        if (z_pos < 0.0) z_pos += lattice_height;
        else if (z_pos > lattice_height) z_pos -= lattice_height;
            
        electron_position[array_index]     = x_pos;
        electron_position[array_index + 1] = y_pos;
        electron_position[array_index + 2] = z_pos;
    }

    int array_index;
    for(int e = 0; e < conduction_electrons; e++) {
      array_index = 3*e;
      electron_position_output_down << "H" << "    " << electron_position[array_index]  << "    " << electron_position[array_index + 1] << "    " <<  electron_position[array_index + 2]  << "\n";    
    }
    
    electron_position_output_down.close(); 
  
}
/*
double reinitialize_electron_conserve_momentum(std::vector<double>& captured_electron_list) {
    int electron_num;
    double Px,Py,Pz,P, radius,length, d_x,d_y,d_z, x_mod,y_mod,z_mod, array_index_i;

    double theta,phi, force, x_distance,y_distance,z_distance;
    double new_potential = 0;

    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> Theta_Vel_distrib(0,2);
    std::uniform_real_distribution<double> Phi_Vel_distrib(0,1);
    std::uniform_real_distribution<double> Theta_pos_distrib(0,2);
    std::uniform_real_distribution<double> Phi_pos_distrib(0,1);
    std::normal_distribution<double> radius_distrib(1,0.1);

   // double new_potential = 0;
  //  std::cout << captured_electron_list.size() << std::endl;
 // std::cout << captured_electron_list.size()/2 << std::endl;
    for(int e = 0; e < captured_electron_list.size()/2; e++) {
      //  std::cout << electron_num << std::endl;
        electron_num = captured_electron_list[2*e];
     //   std::cout << captured_electron_list[2*e+1] << std::endl;
        radius = captured_electron_list[2*e + 1]*radius_distrib(gen);
        if(radius > 1) radius = 2 - radius;
     //   std::cout << radius << std::endl;
        theta = M_PI * Theta_pos_distrib(gen);
        phi = M_PI * Phi_pos_distrib(gen);

        d_x = cos(theta)*sin(phi)*radius;
        d_y = sin(theta)*sin(phi)*radius;
        d_z = cos(phi)*radius;

        new_electron_position[3*electron_num]     = d_x + ((atomic_size * round((new_electron_position[electron_num*3]-1) / atomic_size)) + 1); //closest x atom index
        new_electron_position[3*electron_num + 1] = d_y + ((atomic_size * round((new_electron_position[electron_num*3+1]-1) / atomic_size)) + 1); //closest x atom index;
        new_electron_position[3*electron_num + 2] = d_z + ((atomic_size * round((new_electron_position[electron_num*3+2]-1) / atomic_size)) + 1); //closest x atom index;
    }

    for(int e = 0; e < captured_electron_list.size()/2; e++) {
        new_potential = 0;
        electron_num = captured_electron_list[2*e];
    
        Px = electron_velocity[electron_num*3];
        Py = electron_velocity[electron_num*3 + 1];
        Pz = electron_velocity[electron_num*3 + 2];
        P = (Px*Px)+(Py*Py)+(Pz*Pz);
        
        d_x = new_electron_position[electron_num*3] - ((atomic_size * round((new_electron_position[electron_num*3]-1) / atomic_size)) + 1); //closest x atom index
        d_y = new_electron_position[electron_num*3+1] - ((atomic_size * round((new_electron_position[electron_num*3+1]-1) / atomic_size)) + 1); //closest y atom index
        d_z = new_electron_position[electron_num*3+2] - ((atomic_size * round((new_electron_position[electron_num*3+2]-1) / atomic_size)) + 1); //closest z atom index
        
      //  if(e==0) std::cout << 99*((405* exp(-15* length)) - (exp(-1 * length)))/4 << std::endl;
        
      // std::cout << "d_r: " << sqrt((d_x*d_x)+(d_y*d_y)+(d_z*d_z)) << std::endl;
            //cube around electron
      //  double negative_PE = 0;
        for (int a = 0; a < 1331; a++) {
            
            x_mod = (atomic_size * (a % 11)) - (5*atomic_size);
            y_mod = (atomic_size * ((int(floor(a/11))) % 11)) - (5*atomic_size);
            z_mod = (atomic_size * floor(a / 121)) - (5*atomic_size);
            x_distance = x_mod - d_x;
            y_distance = y_mod - d_y;
            z_distance = z_mod - d_z;

            length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
          //  if (length < 0.00001) length = 0.00001;
            if(length > 100) continue;
            length = sqrt(length);
         //   std::cout << length << std::endl;
            // force
            force = -28*((3.3*3.3 * expl(-3.3 * length)) - (expl(-1 * length)));
           
             
            //TPE += PE;
            new_potential += 28*((3.3* expl(-3.3* length)) - (expl(-1 * length)));
          //  std::cout << "e-a " << new_potential << ", " << length << std::endl;
         //   if(e == 0) std::cout << "force " << force << std::endl;
            phi   = acos(z_distance / length);
            theta = atanl(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            
            electron_force[electron_num]     += force * cos(theta)*sin(phi);
            electron_force[electron_num + 1] += force * sin(theta)*sin(phi);
            electron_force[electron_num + 2] += force * cos(phi);

        }
        // std::cout << new_potential << std::endl;
        for (int i = 1; i < electron_nearest_electron_list[electron_num][0]+1; i++) {
            if (i == electron_num) continue; //no self repulsion
            //if (symmetry_list[e][i]) continue;  //make use of symmetry

            //   electron_spin_two = conduction_electron_spin[i];
            array_index_i = 3*electron_nearest_electron_list[electron_num][i];
            x_distance = new_electron_position[electron_num*3] - new_electron_position[array_index_i];
            y_distance = new_electron_position[electron_num*3+1] - new_electron_position[array_index_i + 1];
            z_distance = new_electron_position[electron_num*3+2] -new_electron_position[array_index_i + 2];

            if (x_distance < -30)     x_distance = x_distance + 40;
            else if (x_distance > 30) x_distance = x_distance - 40;

            if (y_distance < -30)     y_distance = y_distance + 40;
            else if (y_distance > 30) y_distance = y_distance - 40;

            if (z_distance <  -30)    z_distance = z_distance + 40;
            else if (z_distance > 30) z_distance = z_distance - 40;

            length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);

            if (length > 100) continue;
            
            length = sqrt(length);
           //  std::cout << length << std::endl;
            //if (length < 0.000001) length = 0.000001;

          /*  if (electron_spin == electron_spin_two) {
                force = 1e10 * ((constants::K / (length * length * constants::m_e))- (635 * constants::kB / (length * constants::m_e)));
               // std::cout << " 1, 0 " << 1e10 * (constants::K / (length * length * constants::m_e)) << " delta1,1 " << (635 * constants::kB / (length * constants::m_e)) << std::endl;
            }
            else { 
            force = 1 / (length * length);
            
            
            new_potential += force * length;
            //     std::cout << "e-e " << new_potential << ", " << force << ", " << length <<  std::endl;
           // TPE += PE;
           // if(e == 0) std::cout << "force " << force << std::endl;
            //    if (err::check) std::cout << "Distances grabbed" << std::endl; 

            phi   = acos(z_distance / length);
            theta = atanl(y_distance / x_distance);
            if (x_distance < 0) theta += M_PI;

            
            electron_force[electron_num*3]     += force * cos(theta)*sin(phi);
            electron_force[electron_num*3 + 1] += force * sin(theta)*sin(phi);
            electron_force[electron_num*3 + 2] += force * cos(phi);
        }
       //KE1+PE1 = KE2+PE2
       //v2 = sqrt(2KE2/m_e) = sqrt(2(KE1+PE1-PE2)/m_e)
        P = 1e-6*sqrt(abs(2*((electron_potential[electron_num]*constants::K*1e10) + (P*constants::m_e*1e10/2)-(new_potential*constants::K*1e10))/constants::m_e));
      //  std::cout  << electron_potential[electron_num]*constants::K*1e10 + ((Px*Px)+(Py*Py)+(Pz*Pz))*constants::m_e*1e10/2 - new_potential*constants::K*1e10 - P*P*constants::m_e*1e10/2 <<  std::endl;
        theta = M_PI * Theta_Vel_distrib(gen);
        phi = M_PI * Phi_Vel_distrib(gen);
        
              //vel = sqrt(2*KE/m_e) =               TE     -     PE
       // if (e ==0 ) std::cout << electron_potential[e] << std::endl;
      
      //  if(e == 0) std::cout << "KE: " << 0.5*constants::m_e*v_f*v_f << ", PE: " << electron_potential[e]*1e10*constants::K << std::endl;
        electron_velocity[electron_num]     = cos(theta)*sin(phi)*P; 
        electron_velocity[electron_num + 1] = sin(theta)*sin(phi)*P;
        electron_velocity[electron_num + 2] = cos(phi)*P;
       // std::cout << constants::K*1e10*(electron_potential[electron_num]) << ", " << ((Px*Px)+(Py*Py)+(Pz*Pz))*constants::m_e*1e10/2 << ", " << new_potential*constants::K*1e10 << ", " << P*P*constants::m_e*1e10/2 << std::endl;
        TLE += constants::K*1e10*electron_potential[electron_num]+ ((Px*Px)+(Py*Py)+(Pz*Pz))*constants::m_e*1e10/2 - new_potential*constants::K*1e10 - P*P*constants::m_e*1e10/2;
    }
    captured_electron_list.resize(0);
    std::fill(electron_capture.begin(),electron_capture.end(),false);
  //  std::cout << TLE << std::endl;
    return TLE;
}
*/
void initialize_forces() {

    initialize_electron_interactions();
    if (err::check)  std::cout << "Electron interactions" << std::endl;

    initialize_atomic_interactions();
    if (err::check)  std::cout << "Atomic interactions" << std::endl;

    initialize_electron_atom_interactions();
    if (err::check)  std::cout << "Atomic electron interactions" << std::endl;
}

void initialize_electron_interactions() {

    #pragma omp parallel for schedule(static)
    for (int e = 0; e < conduction_electrons; e++) {
      int array_index, array_index_i, neighbor_count = 1;
      double x_distance,y_distance,z_distance, length;
      array_index = 3*e;
      
                if (err::check) if(e ==0) std::cout << "Calculating conduction electron repulsion" << std::endl;

        neighbor_count = 1;
        for (int i = 0; i < conduction_electrons; i++) {
            if (i == e) continue; //no self repulsion

            array_index_i = 3*i;
            x_distance = electron_position[array_index]     - electron_position[array_index_i];
            y_distance = electron_position[array_index + 1] - electron_position[array_index_i + 1];
            z_distance = electron_position[array_index + 2] - electron_position[array_index_i + 2];

            
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
        }
        electron_nearest_electron_list[e][0] = neighbor_count;
     }
  }

void initialize_atomic_interactions() {

    #pragma omp parallel for schedule(static)
    for (int e = 0; e < lattice_atoms; e++) {

      int array_index_i;
      double x_distance,y_distance,z_distance, length;
      int array_index = 3*e;

      int neighbor_count = 2;

      for (int i = 0; i < lattice_atoms; i++) {
          if (i == e) continue; //no self repulsion

        array_index_i = i*3;

        x_distance = atom_position[array_index]     - atom_position[array_index_i];
        y_distance = atom_position[array_index + 1] - atom_position[array_index_i + 1];
        z_distance = atom_position[array_index + 2] - atom_position[array_index_i + 2];

        if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
        else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

        if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
        else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

        if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
        else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

        length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
           
        if(length > a_a_neighbor_cutoff) continue;
       
        atomic_nearest_atom_list[e][neighbor_count] = i;
        neighbor_count++;      
      }
      atomic_nearest_atom_list[e][1] = neighbor_count; 
    
            if(err::check) if(e ==0) std::cout << "Calculating conduction-lattice repulsion" << std::endl;
    }

}

void initialize_electron_atom_interactions() { //we'll need a more developed algorithmn for #electrons != #atoms
    

    #pragma omp parallel for schedule(static)
    for (int e = 0; e < conduction_electrons; e++) {

      int array_index, array_index_a;
      double x_distance,y_distance,z_distance, length;

      array_index = 3*e;
      int nearest_electron_count = 1;
      
        for (int a = 0; a < lattice_atoms; a++) {
            
            array_index_a = 3*a;
            x_distance = electron_position[array_index]     - atom_position[array_index_a];
            y_distance = electron_position[array_index + 1] - atom_position[array_index_a+1];
            z_distance = electron_position[array_index + 2] - atom_position[array_index_a+2]; 

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
        }
        atomic_nearest_electron_list[e][0] = nearest_electron_count;
    }
}

void initialize_velocities() {
     
    
    char directory [256];
    if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog." << std::endl;
    }
    electron_velocity_output.open(string(directory) + "/Electron_Velocity/init.csv");
    electron_velocity_output << "electron number, x-component, y-component, z-component, length" << std::endl;  

    std::ofstream atom_phonon_output;
    atom_phonon_output.open(string(directory) + "/Atom_Energy/init.csv");
    atom_phonon_output << "atom number, energy" << std::endl;
    

    const std::string n = "Init_E_distrib";
    create_fermi_distribution(n, electron_potential,constants::kB_r*Te/E_f_A);
    const std::string na = "P_distrib";
    create_phonon_distribution(na, atom_potential,constants::kB_r*Te/E_f_A);

    #pragma omp parallel for schedule(static)
    for(int e = 0; e < conduction_electrons; e++) {
      double phi,theta; //A/fS
      
      int array_index = 3*e;
      double vel = electron_potential[e];
       // if(vel < E_f_A) std::cout << vel << std::endl;
      vel = sqrt(2.0*vel*constants::m_e_r_i);

        if(sim::CASTLE_x_vector < 0.0 && sim::CASTLE_y_vector < 0.0 && sim::CASTLE_z_vector < 0.0) {
          theta = 2.0*M_PI*omp_uniform_random[omp_get_thread_num()]();
          phi = M_PI * omp_uniform_random[omp_get_thread_num()]();
        } else {
          double unit = sqrt((sim::CASTLE_x_vector*sim::CASTLE_x_vector)+(sim::CASTLE_y_vector*sim::CASTLE_y_vector)+(sim::CASTLE_z_vector*sim::CASTLE_z_vector));
          theta = atan(sim::CASTLE_y_vector / sim::CASTLE_x_vector);
          phi = acos(sim::CASTLE_z_vector / unit);
          if(sim::CASTLE_z_vector < 0.0) theta += M_PI;
        }
      
            if (err::check) if(e==0) std::cout << "Electron velocity ready..." << std::endl;
        electron_velocity[array_index]     = cos(theta)*sin(phi)*vel; 
        electron_velocity[array_index + 1] = sin(theta)*sin(phi)*vel;
        electron_velocity[array_index + 2] = cos(phi)*vel; 
    }

    for(int e = 0; e < conduction_electrons; e++) {
      int array_index = 3*e;
      electron_velocity_output << e << ", " << electron_potential[e] << ", " << electron_velocity[array_index] << ", " << electron_velocity[array_index+1] << ", " << electron_velocity[array_index+2] << std::endl; // ", " << 1e10*constants::K*electron_potential[e] << ", " << 1e10*vel*vel*constants::m_e*0.5 << ", " << 1e10*(electron_potential[e]*constants::K + vel*vel*constants::m_e*0.5) << std::endl;
    //  atom_phonon_output << e << ", " << atom_potential[e] << std::endl;
    }

    electron_velocity_output.close();
    atom_phonon_output.close();
            if (err::check) std::cout << "Electron velocity ready..." << std::endl;
}

double M_B_distrib(const double& epsilon, const double& beta) {

  return (exp(epsilon / beta) / (beta*(exp(epsilon / beta) + 1.0)*(exp(epsilon / beta) + 1.0)));
  
}

void create_phonon_distribution(std::vector<double>& distribution, const double& beta) {

  const double step_size = 1.0 / double(conduction_electrons);
  const double offset  = beta*3.0;
  int count = 0;

 // std::cout << step_size << ", " << offset << std::endl;
  while(count < conduction_electrons) {
    if(beta == 0) {
      distribution[count] = E_f_A;
      count++;
      continue;
    }
    double electron = double(omp_int_random[omp_get_thread_num()]() % conduction_electrons);
    double epsilon = (step_size *electron)  - offset;
    if(omp_uniform_random[omp_get_thread_num()]() < ((0.01 / 16.0)*(epsilon+(3.0*beta))*(epsilon+(3.0*beta))*exp(-0.5*(epsilon+(3.0*beta)) / beta) / (beta*beta*beta))) {
   //  if(E_f_A*(epsilon+1.0) < E_f_A - 3.0*constants::kB_r*Tp) std::cout << E_f_A*(epsilon+1.0) << ", " << E_f_A - (3.0*constants::kB_r*Tp) << ", " << epsilon << ", " << beta  << std::endl;
      distribution[count] = E_f_A*(epsilon + 1.0);
      count++;
    }
  }
}

void create_phonon_distribution(const std::string& name, std::vector<double>& distribution, const double& beta) {

  char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog." << std::endl;
      }
  std::ofstream distrib;
  distrib.open(string(directory) +"/"+ name);
  distrib.precision(20);

  const double step_size = 1.0 / double(conduction_electrons);
  const double offset  = beta*3.0;
  int count = 0;

 // std::cout << step_size << ", " << offset << std::endl;
  while(count < conduction_electrons) {
    if(beta == 0) {
      distribution[count] = E_f_A;
      count++;
      continue;
    }
    double electron = double(omp_int_random[omp_get_thread_num()]() % conduction_electrons);
    double epsilon = (step_size *electron)  - offset;
    if(omp_uniform_random[omp_get_thread_num()]() < ((0.01 / 16.0)*(epsilon+(3.0*beta))*(epsilon+(3.0*beta))*exp(-0.5*(epsilon+(3.0*beta)) / beta) / (beta*beta*beta))) {
     // if(E_f_A*(epsilon+1.0) < E_f_A - 3.0*constants::kB_r*Te) std::cout << E_f_A*(epsilon+1.0) << ", " << E_f_A - (3.0*constants::kB_r*Te) << ", " << epsilon << ", " << beta  << std::endl;
      distribution[count] = E_f_A*(epsilon + 1.0);
      distrib << count << ", " << distribution[count] << "\n";
      count++;
    }
  }
  distrib.close();
}

void create_fermi_distribution(const std::string& name, std::vector<double>& distribution, const double& beta) {

  char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog." << std::endl;
      }
  std::ofstream distrib;
  distrib.open(string(directory) +"/"+ name);
  distrib.precision(20);

  const double step_size = 1.0 / double(conduction_electrons);
  const double offset  = E_f_A + 8.0*constants::kB_r*Te;
  int count = 0;

 // std::cout << step_size << ", " << offset << std::endl;
  while(count < conduction_electrons) {
    // if(beta == 0) {
    //   distribution[count] = E_f_A;
    //   count++;
    //   continue;
    // }
   // double electron = double(count);//double(omp_int_random[omp_get_thread_num()]() % conduction_electrons);
    //double epsilon = step_size *electron*offset/E_f_A;
  //  if(omp_uniform_random[omp_get_thread_num()]() < (1.0/(exp(epsilon/beta) + 1.0))) {
      distribution[count] = count*(E_f_A + 4.0*constants::kB_r*Te)/double(conduction_electrons);
      distrib << count << ", " << distribution[count] << "\n";
      count++;
 //   }
  }
  distrib.close();
}

double return_phonon_distribution(const double& epsilon, const double& beta ) {
    return (1.0/(exp((epsilon-E_f_A)/beta) + 1.0));
}


/*
double e_a_scattering(int e, int a, const double& l_x, const double& l_y, const double& l_z) {

    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> range_distrib(-1,1);
    std::uniform_real_distribution<double> scattering_prob_distrib(0,1);
    std::normal_distribution<double> velocity_gaussian_distrib(0.25,0.1);

    int array_index = 3*e;
    bool collision = false;
   // double atom_energy = atomic_phonon_energy[2*a];

    
    double excitation_constant = velocity_gaussian_distrib(gen);
    
    double Px = electron_velocity[array_index];
    double Py = electron_velocity[array_index+1];
    double Pz = electron_velocity[array_index+2];
    double P = sqrt(Px*Px+Py*Py+Pz*Pz);
    double P_p = sqrt(Px*Px+Py*Py);
    double P_i = sqrt(Pz*Pz);
    double scattering_velocity = P;
    double electron_TE = 1e10*(electron_potential[e]*constants::K + P*P*constants::m_e*0.5);
   
   

    double d_p = sqrt(l_x*l_x+l_y*l_y);
    double d_i = sqrt(l_z*l_z);
    double d_r = sqrt((l_x*l_x)+(l_y*l_y)+(l_z*l_z));
    double polar_value = M_PI * range_distrib(gen);
    double incline_value = M_PI * range_distrib(gen);
    double normal_polar_angle = acos((l_x*Px+l_y*Py)/(P_p*d_p));
    double normal_incline_angle = acos((l_z*Pz) / (P_i*d_i));

    double excitation_energy = abs(excitation_constant)*(electron_TE - atomic_phonon_energy[2*a]);

    
    if (excitation_energy > 0) {
      //  std::cout << electron_TE << ", " << atomic_phonon_energy[2*a] << ", " << excitation_energy << ", " << electron_TE - excitation_energy << std::endl;
        excitation_energy = electron_TE - excitation_energy;
        
        scattering_velocity = 1e-5*sqrt(2*excitation_energy/constants::m_e);
       // std::cout << P - scattering_velocity << std::endl;
      // if(scattering_velocity < 0) std::cout << P << ", " << scattering_velocity << ", " << P - scattering_velocity <<  std::endl;
    

        if(normal_polar_angle > M_PI) normal_polar_angle = -1*(normal_polar_angle-M_PI);
        else if (normal_polar_angle < -1*M_PI) normal_polar_angle = -1*(normal_polar_angle+M_PI);

        if(normal_incline_angle > M_PI) normal_incline_angle = -1*(normal_incline_angle-M_PI);
        else if (normal_incline_angle < -1*M_PI) normal_incline_angle = -1*(normal_incline_angle+M_PI);

    
        double polar_scattering_angle = atanl(Py/Px);
        if(polar_scattering_angle < 0) polar_scattering_angle += M_PI;
        double incline_scattering_angle = acos(Pz/P);

        double polar_prob = velocity_gaussian_distrib(gen)*polar_value*polar_value* ( ( ((M_PI/2)+normal_polar_angle) * exp( (polar_value-1)*(polar_value-1)/(-8)) )+( ((M_PI/2)-normal_polar_angle) * exp( (polar_value+1)*(polar_value+1)/(-8)) ) )/(d_r*2*M_PI*sqrt(2*M_PI));
        double incline_prob = velocity_gaussian_distrib(gen)*incline_value*incline_value* ( ( ((M_PI/2)+normal_incline_angle)  * exp( (incline_value-1)*(incline_value-1)/(-8)) ) + ( ((M_PI/2)-normal_incline_angle) * exp( (incline_value-1)*(incline_value-1)/(-8)) ) )/(d_r*2*M_PI*sqrt(2*M_PI));
    
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
        
            if(scattering_velocity < P) {
              //  if(electron_TE - atomic_phonon_energy[2*a] < 0) std::cout << "TE: " << electron_TE << ", E_f: " << atomic_phonon_energy[2*a] << std::endl;
                P = scattering_velocity;
            }
            
        
            #pragma omp critical
            {
        
            atomic_phonon_energy[2*a] += excitation_energy;
            TLE += excitation_energy;
            }
        } 
        //#pragma omp critical
     //   std::cout << "Scattering Velocity: " << scattering_velocity << ", incoming_velocity" << P << std::endl; //", polar_probability " << polar_prob << ", incline_prob " << incline_prob << ", normal_polar_angle " << normal_polar_angle << ", " << ", polar_scattering_angle " << polar_scattering_angle << ", incline_scattering_angle" << incline_scattering_angle << std::endl;
        //#pragma omp critical
      //  TLE += P*P - scattering_velocity*scattering_velocity;
       // #pragma omp critical
        //if(TLE != 0) std::cout << TLE << std::endl;
     //   std::cout << scattering_velocity*1e5 << std::endl;
        electron_velocity[array_index]   = P * cos(polar_scattering_angle)*sin(incline_scattering_angle);
        electron_velocity[array_index+1] = P * sin(polar_scattering_angle)*sin(incline_scattering_angle);
        electron_velocity[array_index+2] = P * cos(incline_scattering_angle); 
    }

    //TLE += atomic_phonon_energy[2*a];
    return 0;//TLE;
}

double e_p_scattering(int e, int a, const double& x_distance, const double& y_distance, const double& z_distance) {
   
    int b;
    
    double d_x,d_y,d_z;
   
    std::srand(std::time(nullptr));
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> phonon_jump_distrib(0,1);
    std::uniform_int_distribution<> phonon_random_walk(0,26);

    double atomic_energy = atomic_phonon_energy[2*a];
    int p_p_coupling;
    int p_e_reverse_coupling; // = atomic_phonon_energy[2*a + 1];
    int array_index = 3*e;
    double scattering_velocity = sqrt((electron_velocity[array_index]*electron_velocity[array_index]) + (electron_velocity[array_index+1]*electron_velocity[array_index+1]) + (electron_velocity[array_index+2]*electron_velocity[array_index+2]));
    double length = sqrt((x_distance*x_distance)+(y_distance*y_distance)+(z_distance*z_distance));
    double electron_energy = 0.5*constants::m_e*1e10*scattering_velocity*scattering_velocity;
    double excitation_energy = exp(-1*length)*(electron_energy - atomic_energy);
    if(excitation_energy > 0) {
        if(phonon_jump_distrib(gen) > 0.2*exp(-1*excitation_energy)) {
            #pragma omp critical 
            {
            atomic_phonon_energy[2*a] += excitation_energy;
            atomic_phonon_energy[2*a +1]++;
            TLE += excitation_energy;
            }
            scattering_velocity -= 1e-5*sqrt(2*excitation_energy/constants::m_e);
            
        }
    } else {
        if(phonon_jump_distrib(gen) < 0.1*p_e_reverse_coupling*exp(-1*length)) {
            #pragma omp critical 
            {
            atomic_phonon_energy[2*a] += excitation_energy;
            atomic_phonon_energy[2*a+1] = 0;
            TLE += excitation_energy;
            }
            scattering_velocity += 1e-5*sqrt(-2*excitation_energy/constants::m_e);
        
        } else {
            #pragma omp critical
            atomic_phonon_energy[2*a+1]++;
        }
    }
    return 0;
} 
*/
void output_data() {
  
    //=========
    // Output equilibration step data
    //=========

    if((current_time_step % (CASTLE_output_rate * CASTLE_MD_rate)) == 0) {
      time_stamp = std::to_string(current_time_step);
    
      char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog." << std::endl;
      }

      std::ofstream atomic_phonon_output;
      atomic_phonon_output.open(string(directory) + "/Atom_Energy/" + time_stamp);
      atomic_phonon_output.precision(10);
      atomic_phonon_output << std::scientific;

      // std::ofstream atomic_phonon_output_hist;
      // atomic_phonon_output_hist.open(string(directory) + "/Atom_Energy/" + time_stamp + ".hist");
      // atomic_phonon_output_hist.precision(10);
      // atomic_phonon_output_hist << std::scientific;
      // electron_position_output_down.open(string(directory) + "/Electron_Position/" + time_stamp + ".xyz");
      // electron_position_output_down << conduction_electrons << "\n";
      // electron_position_output_down << time_stamp << "\n";
      // electron_position_output_down.precision(10);
      // electron_position_output_down << std::scientific;

      electron_velocity_output.open(string(directory) + "/Electron_Velocity/" + time_stamp + ".csv");
      electron_velocity_output << "Electron number,    x-component,     y-component,    z-component,     length, energy" << "\n";
      electron_velocity_output.precision(10);
      electron_velocity_output << std::scientific;
    
      std::ofstream temp_map_list [8];
     
     
      double tempMap[8];
      for(int i = 0; i < 8; i++) {
         tempMap[i] = 0.0;
         temp_map_list[i].open(string(directory) + "/Temp_Map" + std::to_string(i) + "/" + time_stamp);
      //electron_temperature_map.precision(10);
     // electron_temperature_map << std::scientific;
      }
    double  velocity_length; 
    double x_vel, y_vel, z_vel, x_pos, y_pos, z_pos;
    int array_index, array_index_y, array_index_z;

    for(int e = 0; e < conduction_electrons; e++) {
      array_index   = 3*e;
      array_index_y = array_index + 1;
      array_index_z = array_index + 2;
      x_pos = electron_position[array_index];
      y_pos = electron_position[array_index_y]; 
      z_pos = electron_position[array_index_z];
       x_vel = 1e5*electron_velocity[array_index];
       y_vel = 1e5*electron_velocity[array_index+1];
       z_vel = 1e5*electron_velocity[array_index+2];
      velocity_length = electron_potential[e];
      if(x_pos < (lattice_width * 0.5) && y_pos < (lattice_depth*0.5) && z_pos < (lattice_height * 0.5)) {
        tempMap[0] += velocity_length;
        temp_map_list[0] << e << ", " << x_pos << ", " << y_pos << ", " << z_pos << ", " << velocity_length << "\n";
      }
      if(x_pos > lattice_width * 0.5 && y_pos < lattice_depth*0.5 && z_pos < lattice_height * 0.5) {
        tempMap[1] += velocity_length;
        temp_map_list[1] << e << ", " << x_pos << ", " << y_pos << ", " << z_pos << ", " << velocity_length << "\n";
      }
      if(x_pos < lattice_width * 0.5 && y_pos > lattice_depth*0.5 && z_pos < lattice_height * 0.5) {
        tempMap[2] += velocity_length;
        temp_map_list[2] << e << ", " << x_pos << ", " << y_pos << ", " << z_pos << ", " << velocity_length << "\n";
      }
      if(x_pos > lattice_width * 0.5 && y_pos > lattice_depth*0.5 && z_pos < lattice_height * 0.5) {
        tempMap[3] += velocity_length;
        temp_map_list[3] << e << ", " << x_pos << ", " << y_pos << ", " << z_pos << ", " << velocity_length << "\n";
      }
      if(x_pos < lattice_width * 0.5 && y_pos < lattice_depth*0.5 && z_pos > lattice_height * 0.5) {
        tempMap[4] += velocity_length;
        temp_map_list[4] << e << ", " << x_pos << ", " << y_pos << ", " << z_pos << ", " << velocity_length << "\n";
      }
      if(x_pos > lattice_width * 0.5 && y_pos < lattice_depth*0.5 && z_pos > lattice_height * 0.5) {
        tempMap[5] += velocity_length;
        temp_map_list[5] << e << ", " << x_pos << ", " << y_pos << ", " << z_pos << ", " << velocity_length << "\n";
      }
      if(x_pos < lattice_width * 0.5 && y_pos > lattice_depth*0.5 && z_pos > lattice_height * 0.5) {
        tempMap[6] += velocity_length;
        temp_map_list[6] << e << ", " << x_pos << ", " << y_pos << ", " << z_pos << ", " << velocity_length << "\n";
      }
      if(x_pos > lattice_width * 0.5 && y_pos > lattice_depth*0.5 && z_pos > lattice_height * 0.5) {
        tempMap[7] += velocity_length;
        temp_map_list[7] << e << ", " << x_pos << ", " << y_pos << ", " << z_pos << ", " << velocity_length << "\n";
      }

      electron_velocity_output << e << ", " << velocity_length << ", " << x_vel << ", " << y_vel << ", " << z_vel << "\n";
      atomic_phonon_output << e << ", " << atom_potential[e] << "\n";
    }
    // electron_position_output_down.close();
    electron_velocity_output.close();
    atomic_phonon_output.close();
    
    for(int i = 0; i < 8; i++) {
      temp_map_list[i].close();
      tempMap[i] = sqrt(1.5*e_heat_capacity_i*1e-30*(tempMap[i] - (E_f_A*8.0/double(conduction_electrons)))* n_f*300.0*8.0/double(conduction_electrons));
    }
    temp_data << current_time_step << ", " << tempMap[0] << ", " << tempMap[1] << ", " << tempMap[2] << ", " << tempMap[3] << ", "\
              << tempMap[4] << ", " << tempMap[5] << ", " << tempMap[6] << ", " << tempMap[7] << std::endl;
  
      std::cout << "  " << current_time_step / total_time_steps * 100 << "%. " << std::endl; 
    // for(int e = 0; e < conduction_electrons; e++) {
    //   atomic_phonon_output_hist << w*double(e)<< ", " << double(count[e]) / double(conduction_electrons) << "\n";
    // }
    // atomic_phonon_output_hist.close();


  //  const std::string name = "E_distrib";
    //create_phonon_distribution(name, electron_potential, Te*e_specific_heat / E_f_A);
  // if(hot_electrons > 0) {
  //   std::ofstream electron_hot_output;
  //   electron_hot_output.open(string(directory) + "/Hot_Electrons/" + time_stamp + ".xyz");
  //   electron_hot_output << hot_electrons << "\n";
  //   electron_hot_output << time_stamp << "\n";
  //   electron_hot_output.precision(10);
  //   electron_hot_output << std::fixed;
  //   std::ofstream atom_hot_output;
  //   atom_hot_output.open(string(directory) + "/Hot_Atoms/" + time_stamp + ".xyz");
  //   atom_hot_output << hot_atoms << "\n";
  //   atom_hot_output << time_stamp << "\n";
  //   atom_hot_output.precision(10);
  //   atom_hot_output << std::fixed;
  //   for(int e = 0; e < hot_electrons; e++) {
  //     array_index = 3*hot_e_list.front();
  //     hot_e_list.pop_front();
  //     electron_hot_output << "H" << ", " << new_electron_position[array_index] << ", " << new_electron_position[array_index+1] << ", " << new_electron_position[array_index+2] << "\n";
  //   }  
  //   for(int e = 0; e < hot_atoms; e++) {
  //     array_index = 3*hot_a_list.front();
  //     hot_a_list.pop_front();
  //     atom_hot_output << "H" << ", " << atom_position[array_index] << ", " << atom_position[array_index+1] << ", " << atom_position[array_index+2] << "\n";
  //   }  
  //   electron_hot_output.close();
  //   atom_hot_output.close();
  // }
  } 
    double mean_ea_rad = 0.0;
    double mean_ee_rad = 0.0;

    // #pragma omp parallel for reduction(+:mean_ee_rad, mean_ea_rad)
    // for (int e = 0; e < conduction_electrons; e++) {
    //  // mean_ea_rad += sqrt(mean_radius[2*e]);
    //   mean_ee_rad += sqrt(mean_radius[2*e + 1]);
    //  // mean_radius[2*e] = 1.0;
    //   mean_radius[2*e + 1] = 1.0;
    // }
  //  mean_ea_rad /= conduction_electrons;
  //  mean_ee_rad /= conduction_electrons;
  
    mean_data.precision(10);
    mean_data << std::scientific;

    double j = x_flux * constants::e * 1e20 / (1600 * CASTLE_output_rate * dt); //current density
    double nu = j / (n_f * constants::e); //drift velocity
    double I = n_f * 1600 * 1e-20 * nu * constants::e; //current

    if(!current_time_step) {
    mean_data << CASTLE_real_time << ", " << current_time_step << ", " 
      << Te*Te*e_heat_capacity * 1e-20/300.0 << ", " << Tp*a_heat_capacity*1e-20 << ", " 
      << Te << ", " << Tp << ", "  \
      << TTMe << ", " << TTMp << ", "
      << mean_ea_rad << ", " << mean_ee_rad << ", " << a_a_scattering_count << ", " << e_a_scattering_count << ", " << e_e_scattering_count  << ", " << x_flux << ", " << y_flux << ", " << z_flux  << ", " \
       << std::endl;
    }

    else {
    mean_data << CASTLE_real_time << ", " << current_time_step << ", " 
      << Te*Te*e_heat_capacity * 1e-20/300.0 << ", " << Tp*a_heat_capacity*1e-20 << ", "  
      << Te << ", " << Tp << ", " //<< TEKE << ", " << TLE << ", " 
      << TTMe << ", " << TTMp << ", "
      << mean_ea_rad << ", " << mean_ee_rad << ", " << std::fixed; mean_data.precision(1); mean_data << double(a_a_scattering_count) / CASTLE_output_rate << ", " << double(e_a_scattering_count) / CASTLE_output_rate << ", " << double(e_e_scattering_count) / double(CASTLE_output_rate) << ", " << double(x_flux) / double(CASTLE_output_rate) << ", " << double(y_flux) / CASTLE_output_rate << ", " << double(z_flux) / double(CASTLE_output_rate)  << ", " \
      << std::endl;
    }
   
    x_flux = 0;
    y_flux = 0;
    z_flux = 0;
    e_a_scattering_count = 0;
    e_e_scattering_count = 0;
    a_a_scattering_count = 0;
}




} //end of CASTLE namespace



