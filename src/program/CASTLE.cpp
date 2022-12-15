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

// double return_fermi_distribution(const double& epsilon, const double& beta );
void create() {

            if (err::check) std::cout << "Creating CASTLE..." << std::endl; 
            if (err::check) std::cout << "  ";
  
  /*  double x,y,z,r,r_min, d_x,d_y,d_z, p_x,p_y,p_z, d_p_x,d_p_y,d_p_z, p, force, f_x,f_y,f_z, d_f_x,d_f_y,d_f_z, theta, phi, potential;
    double a_x,a_y, a_d_x,a_d_y,a_d_z, a_p_x,a_p_y, a_f_x,a_f_y, a_d_f_x,a_d_f_y;
    std::ofstream electron_output;
    std::ofstream ballistic_data;

    dt = 1e-1; //fs

    a_x = 0;
    a_y = 0;

    x = -7;
    y = 0;

    p_x = 15.0; //A/fs 
    p_y = 0;

ballistic_data.open("Ballistic_Data");
for(int a = 0; a < 1200; a++) {
    if(a == 600) continue;
    r = 7;
    r_min = r;
    int t = 0;
    std::string height  = std::to_string(a);
    electron_output.open("Ballistic_Electron/" + height + ".xyz");
    electron_output << "Ballistic Electron" << "\n" << "\n";

    theta = 1.25*M_PI*(1.0 - (0.5*(a/1500.0)));
    p_x = 15.0; //A/fs 
    p_y = 0;
    x = r*cos(theta);
    y = r*sin(theta);

   force = exp(-0.5*r)*(0.5*r + 1.0) / (r*r); //N -> Kg A/fs**2 
    //if(r < 1.6) force = (sqrt(3.0)*exp(-0.5*r)*(0.5*r+1)/(r*r)) - (sqrt(3)*(exp(-0.5*r)*(2.97182 - 2.47652*r))/((r - 3.2)*(r-3.2)));
    // else force = 0.0;
    ballistic_data << theta << ", " << theta / M_PI << ", ";
    f_x = force * cos(theta);
    f_y = force * sin(theta);
    electron_output << t << " " << x << " " << y << "\n";
    while(r <= 7.01) {
        
        d_x = x + (p_x*dt) ;//+ (f_x*dt*dt / (2*constants::m_e_r));
        d_y = y + (p_y*dt) ;//+ (f_y*dt*dt / (2*constants::m_e_r));
       
        //update forces
        r = sqrt(((d_x -a_x)*(d_x - a_x)) + ((d_y - a_y)*(d_y - a_y)));
        if(r < r_min) r_min = r;

        force = exp(-0.5*r)*(0.5*r + 1.0) / (r*r); //N -> Kg A/fs**2 
        // if(r < 1.6) force = (sqrt(3.0)*exp(-0.5*r)*(0.5*r+1)/(r*r)) - (sqrt(3.0)*(exp(-0.5*r)*(2.97182 - 2.47652*r))/((r - 3.2)*(r-3.2)));
        // else force = 0.0;

        electron_output << t+1 << ", " << d_x << ", " << d_y << " " <<  (force)*1e-1*constants::m_e_r_i << " "\
          << sqrt((p_x*p_x)+(p_y*p_y))*1e-1 + force*1e-1*1e-1/(2*constants::m_e_r)<< "\n";

        theta = atan2(d_y,  d_x);
        if(theta != theta) theta = 0.0;
  
        d_f_x = force * cos(theta);
        d_f_y = force * sin(theta);

        //update velocity
        d_p_x =p_x+ ((d_f_x + f_x)*dt/ (2*constants::m_e_r));
        d_p_y =p_y+ ((d_f_y + f_y)*dt/ (2*constants::m_e_r));

        phi = atan2(d_p_y, d_p_x);
          if(phi!=phi) phi = 0.0;

        p_x = 15.0*cos(phi);
        p_y = 15.0*sin(phi);

        x = d_x;
        y = d_y;

        f_x = d_f_x;
        f_y = d_f_y;

        t++;

        phi = atan2(p_y, p_x);
        if(phi != phi) phi = 0.0;
       // electron_output << t << " " << x << " " << y << "\n";
    }
    electron_output.close();
    ballistic_data << theta / M_PI << ", " << r_min << ", " << phi / M_PI << ", " << sqrt(p_x*p_x + p_y*p_y) <<  "\n";

} ballistic_data.close();   */

        std::cout << "Building CASTLE..." <<std::endl;

    stopwatch_t castle_watch;
    castle_watch.start();
     omp_threads = sim::CASTLE_omp_threads;
    equilibrium_step = true;
    //========
    // Initialize the lattice, electrons, and starting forces
    //========
            if (err::check) std::cout << "Prepare to initialize..." << std::endl;
       omp_set_dynamic(0);
       omp_set_num_threads(omp_threads);
      #pragma omp parallel 
            if(omp_get_thread_num() == omp_threads-1) std::cout << "OpenMP capability detected. Parallelizing integration for " << omp_get_thread_num() +1<<  " threads" << std::endl;
        
    initialize();
    
    
        std::cout << "CASTLE build time.at(s): " << castle_watch.elapsed_seconds() << std::endl;
  
        std::cout << "Storming CASTLE..." << std::endl;
   
 
    //========
    // Integrate total time steps
    //========
    castle_watch.start();
    sim::integrate(total_time_steps);
    
        std::cout << "Average time step.at(s):  " << (castle_watch.elapsed_seconds()) / total_time_steps << std::endl;

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

    std::cout << "Averaging complete. " << castle_watch.elapsed_seconds() << " .at(s) elapsed. Step-time: " << castle_watch.elapsed_seconds() / total_time_steps << std::endl;
    
    mean_data.close();
    //temp_data.close();
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
    conduction_electrons = atoms::num_atoms;//int(round(0.6*double(atoms::num_atoms)));
    lattice_atoms = atoms::num_atoms; //Better lattice creation will come from VAMPIRE in future
    CASTLE_output_rate = sim::partial_time;
    dt = mp::dt_SI * 1e15;//-4; //S -> femptoSeconds
    TTMe = TTMp = d_TTMe = d_TTMp = Tp = Te = d_Tp = d_Te = sim::temperature;
    // d_TTMp = TTMp = Tp = d_Tp = 300.0;
    
    // Te += 1.0;
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
    lattice_width = cs::system_dimensions[0]+x_unit_size; //A
    lattice_depth  = cs::system_dimensions[1]+y_unit_size; // A
    lattice_height  = cs::system_dimensions[2]+z_unit_size; // A

    // if( lattice_height < 26 || lattice_depth < 26 || lattice_width < 26) {
    //   std::cerr << "Lattice dimension Error: castle must be 4.5x4.5x4.5 nm^3 minimum" << std::endl;
    // }

    x_flux = 0;
    y_flux = 0;
    z_flux = 0;
    current_time_step = 0;
    CASTLE_real_time = 0;
    // uniform_random.seed(2137082040);
    // int_random.seed(2137082040);

    mu_r = (atomic_mass + constants::m_e_r) / (atomic_mass * constants::m_e_r );
    combined_mass = 1 / (atomic_mass + constants::m_e_r);
    applied_voltage = sqrt(1.60218e-19*2.0*sim::applied_voltage/constants::eps_0/(1e-30*lattice_width * lattice_height * lattice_depth)); //eV -> V/m
    n_f = 1.0e3 * conduction_electrons / (lattice_width * lattice_height * lattice_depth); // e- / A**3 -> e-/m**3
    E_f = constants::h * constants::h * pow(3.0 * M_PI * M_PI * n_f*1e27, 0.66666666666666666666666667) / (8.0 * M_PI * M_PI * constants::m_e); //Fermi-energy
    E_f_A = E_f*1e20; //AJ
    mu_f = 5.0 * E_f_A / (3.0 * conduction_electrons);//Fermi-level
    v_f = sqrt(2.0 * E_f_A * constants::m_e_r_i); // A/fs

    a_specific_heat = 25.0 /6.02e3; // //sim::TTCl; //J/K/mol -> .at(e20/Na) AJ/K/e-   
    a_specific_heat_i = 1.0 / a_specific_heat;
    a_heat_capacity = 1.14*a_specific_heat * n_f; //AJ/K/particle .at() AJ/K/nm**3
    a_heat_capacity_i = 1.0 / a_heat_capacity;

    e_specific_heat = 0.5*M_PI*M_PI*constants::kB*constants::kB/E_f; // gamma; //J/K**2/e- 
    e_specific_heat_i = 1.0 / e_specific_heat;
    e_heat_capacity = (10.4)*(1.0+0.084)*1e20*e_specific_heat * n_f; //J/K/e- -> AJ/K**2/nm**3
    e_heat_capacity_i = 1.0 / e_heat_capacity;

    ea_coupling_strength = 1e-6*sim::ea_coupling_strength*constants::eV_to_AJ*constants::eV_to_AJ/constants::hbar_r; // meV**2 -> AJ/fs
    phonon_energy = 1e-3*sqrt(sim::ea_coupling_strength/0.084)*constants::eV_to_AJ; // meV [e-3] AJ
    std::cout << "phonon occupation test: " << return_BE_integrand(phonon_energy, Te) << ", " << return_BE_integrand(phonon_energy, Tp) << std::endl;
    // std::cout << "phonon energy " << phonon_energy << std::endl;

    atomic_mass = 58.69 * 1.6726219e3; // kg * 1e30 for reduced units
    power_density = 1e1*sim::fluence; // mJ/cm**2 -> .at(e17/e14/e2(fs)) AJ/fs/nm**2
    
    const static double tau = E_f_A /(M_PI*ea_coupling_strength); // fs/AJ
    G = 300.0*e_heat_capacity*E_f_A*2.0/tau; //AJ/fs/K/nm**3 [e-20*e27*e15 = e22]  
    //G = sim::TTG*1e-23;
    //G=Ce/Te-p = pihbar constant (meV**2)Ef*n_f*(1/eps)**2
    ea_rate = -30.0*dt*E_f_A/tau;  //AJ(ready for E_i)  AJfs/fs
    ee_rate = -1.0*dt*sim::ee_coupling_strength/(constants::eV_to_AJ*constants::eV_to_AJ); //eV^-2 fs^-1 -> fs**-1 AJ**-2

    omp_set_num_threads(omp_threads);
  // omp_uniform_random(omp_threads);
  //  int_random.seed(omp_threads);
  //  omp_int_random.resize(omp_threads);
  //  // std::cout << omp_get_num_threads() << std::endl;
  //   std::srand(std::time(nullptr));
  //   std::random_device rd;  //Will be used to obtain a seed for the random number engine
  //   std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  //   std::uniform_int_distribution<> test_int(0, conduction_electrons -1);
  //   std::uniform_real_distribution<double> test_uniform;
    #pragma omp parallel
    {
    for(uint32_t i = 0; i < omp_threads; i++) {
    // //  // std::cout << omp_get_num_threads() << std::endl;
     omp_int_random.at(i).seed(i);
     omp_uniform_random.at(i).seed(i);
    //   omp_int_random.at(i).seed(i);
    }
    }
    /*std::cout << "testing grnd() bounds..." << std::endl;
    #pragma omp parallel 
    {
    const int thread = omp_get_thread_num();
    for( int e = 0; e < conduction_electrons*conduction_electrons; e++) {
      if(int(omp_uniform_random.at(thread)()*2147483647)%conduction_electrons < 0 || int(omp_uniform_random.at(thread)()*2147483647)%conduction_electrons >= conduction_electrons) std::cout << "grnd() error " << std::endl;
    }
    for( int e = 0; e < conduction_electrons*conduction_electrons; e++) double test = omp_uniform_random.at(thread)();
    }
    std::cout << "grnd passed." << std::endl; */
   // int test = omp_int_random[0]() % conduction_electrons;
    //test = omp_int_random.at(24)MTRand_int32::rand_int32() %conduction_electrons ;
    // Initialize lattice
    //========
    initialize_positions();
            if (err::check) std::cout << "Lattice built " << std::endl;
    
    //========
    // initialize electrons: lattice and conduction bands, velocity, spin, etc.
    //=======
   // initialize_forces();
            if (err::check) std::cout << "Electrons ready..." << std::endl;
    //=======
    // Calls forces set up
    //=======
            if (err::check) std::cout << "Forces ready..." << std::endl;

   initialize_velocities();

             if (err::check) std::cout << "Particles a movin" << std::endl;
    dos_occ  = round(phonon_energy*(conduction_electrons/(E_f_A-core_cutoff)));
    local_dos_occ = round(0.9*phonon_energy*(ee_density/3.0/(E_f_A-core_cutoff)));
    std::cout << "global dos occupation: "<< dos_occ  << ", local dos occupation: " << local_dos_occ << std::endl;

    std::cout << "E_f(AJ): " << E_f_A << std::scientific << ", gamma(J/m**3/K**2): " << e_heat_capacity*1e7 << ", C_l(J/K/m**3): " << a_heat_capacity*1e7 << ", G@300K(J/K/s/m**3): " <<  G*1e22  << \
    ", ea_rate@300K(J/s/K/m**3): " << -1e22*ea_rate*n_f/300.0 <<  ", tau(fs/AJ): " << tau/E_f_A << ", photon max rate: " << 1e-2*power_density*lattice_width*lattice_depth/(sim::applied_voltage*constants::eV_to_AJ) << std::fixed << std::endl;
   //G = -1.0*ea_rate*n_f/300.0; //J/s/K/m**3 []
     initialize_cell_omp(); 
  // else std::cout << "CASTLE lattice integration is most efficient at greater than 4 15 Angstrom unit cells wide. Generating standard OpenMP lattice." << std::endl;
    // std::cout << "electron DoS summation: " << int(round(4*ee_density/3.0/(3*constants::kB_r*300.0+E_f_A - core_cutoff)/4.0))  << std::endl;
    
   if (err::check)  std::cout << "Total lattice cells: " << total_cells << ", maximum cell size: " << cell_integration_lists[0].size() << ", maximum integration list size: " << electron_integration_list[0].size() << std::endl;
    char directory [256];
    if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog. data" << std::endl;
    }
    
   // temp_data.open(string(directory) + "/temp_data.csv");
    mean_data.open(string(directory) + "/mean_data.csv");
    mean_data << "time, step, mean-EKE, mean-LE, mean-Te, mean-Tp,  mean-radius, mean-e-a-collisions, mean-e-e-collisions, mean-x_flux, mean-y_flux, mean-z_flux" << "\n";

}


void initialize_cell_omp() {

  // double 
  // x_omp_cells = int(floor(lattice_width / 20.0));
  // y_omp_cells = int(floor(lattice_depth / 20.0));
  // z_omp_cells = int(floor(lattice_height/ 20.0));
  x_omp_cells = 8;
  y_omp_cells = 8;
  z_omp_cells = 8;

  total_cells = x_omp_cells*y_omp_cells*z_omp_cells;

  x_step_size = lattice_width / double(x_omp_cells);
  y_step_size = lattice_depth / double(y_omp_cells);
  z_step_size = lattice_height/ double(z_omp_cells);
  boundary_conditions_cutoff = fmax(x_step_size, fmax(y_step_size, z_step_size));
  cell_integration_lists.resize(total_cells);
  old_cell_integration_lists.resize(total_cells);

  for(int i=0; i < total_cells; i++) {
    cell_integration_lists[i].resize(int(4*double(conduction_electrons) / double(total_cells)), 0);
    old_cell_integration_lists[i].resize(int(4*double(conduction_electrons) / double(total_cells)), 0);
    cell_integration_lists[i][0] = 1;
    old_cell_integration_lists[i][0] = 1;
  }

  escaping_electrons.resize( int(double(conduction_electrons) / double(total_cells)));
  escaping_electrons[0] = 1;
  lattice_cell_coordinate.resize(x_omp_cells);
  cell_lattice_coordinate.resize(total_cells);
  for(int c = 0; c < total_cells; c++) {
    cell_lattice_coordinate.at(c).resize(3);
  }
  int cell_count = 0;
    if(err::check)   std::cout << "cell integration arrays generated." << std::endl;

  //omp lattice coordinates map
  for(int i = 0; i < x_omp_cells; i++) {
    lattice_cell_coordinate.at(i).resize(y_omp_cells);
    for(int k = 0; k < y_omp_cells; k++) {
      lattice_cell_coordinate.at(i).at(k).resize(z_omp_cells);
      for(int j = 0; j < z_omp_cells; j++) {
        lattice_cell_coordinate.at(i).at(k).at(j) = cell_count;
        cell_lattice_coordinate.at(cell_count)[0] = i;
        cell_lattice_coordinate.at(cell_count)[1] = k;
        cell_lattice_coordinate.at(cell_count)[2] = j;
        cell_count++;
      }
    }
  }
  if(err::check)   std::cout << "lattice coordinate arrays initiated." << std::endl;
 
  //lattice cell division
  //int current_lattice_current_end;
//  #pragma omp parallel
 
  // if(omp_get_num_threads() != omp_threads) std::cout << "omp pragma error " << omp_get_num_threads() << " != " << omp_threads << std::endl;
  // #pragma omp for 
  for(int e = 0; e < conduction_electrons; e++) {
    const int array_index = 3*e;
     int x_cell = int(floor(electron_position.at(array_index) / x_step_size));
      //if (x_cell < 0 || x_cell > x_omp_cells) std::cout << electron_position.at(array_index) << ", " << x_step_size << ", " << floor(electron_position.at(array_index) / x_step_size) << std::endl;

     int y_cell = int(floor(electron_position.at(array_index+1) / y_step_size));
      //if (y_cell < 0 || y_cell > y_omp_cells) std::cout << electron_position.at(array_index+1) << ", " << y_step_size << ", " << floor(electron_position.at(array_index+1) / y_step_size) << std::endl;
    
     int z_cell = int(floor(electron_position.at(array_index+2) / z_step_size));
      //if (z_cell < 0 || z_cell > z_omp_cells) std::cout << electron_position.at(array_index+2) << ", " << z_step_size << ", " << floor(electron_position.at(array_index+2) / z_step_size) << std::endl;
    
    if(x_cell != x_cell || y_cell != y_cell || z_cell != z_cell ) std::cout << x_cell << ", " << y_cell << ", " << z_cell << ", " <<\
    int(floor(electron_position.at(array_index) / x_step_size)) << ", " << int(floor(electron_position.at(array_index+1) / y_step_size)) << ", " << \
    int(floor(electron_position.at(array_index+2) / z_step_size)) << std::endl;

       if(x_cell >= x_omp_cells || y_cell >= y_omp_cells || z_cell >= z_omp_cells) {std::cout << "cell sorting for ee integration exceeds bounds " << \
      x_cell << ", " << y_cell << ", " << z_cell << std::endl;
      x_cell --;
      y_cell --;
      z_cell --; }
          if(x_cell < 0 || y_cell < 0 || z_cell < 0) {std::cout << "cell sorting for ee integration less than zero " << \
      x_cell << ", " << y_cell << ", " << z_cell << std::endl;
      x_cell = 0;
      y_cell = 0;
      z_cell = 0; }
     const int omp_cell = lattice_cell_coordinate.at(x_cell).at(y_cell).at(z_cell); 
     int current_lattice_current_end = cell_integration_lists.at(omp_cell)[0];

    cell_integration_lists.at(omp_cell).at(current_lattice_current_end) = e;
    old_cell_integration_lists.at(omp_cell).at(current_lattice_current_end) = e;
  
    cell_integration_lists.at(omp_cell)[0]++;
    old_cell_integration_lists.at(omp_cell)[0]++;
  
          if(current_lattice_current_end >= cell_integration_lists.at(omp_cell).size()) std::cout << \
            omp_cell << ", " << current_lattice_current_end << ", " << cell_integration_lists.at(omp_cell).size() << ", " << e << std::endl;
      
  }
 
  cell_nearest_neighbor_list.resize(total_cells);
  std::vector<int> spiral_cell_counter;
  spiral_cell_counter.resize(total_cells, 0);
      if(err::check) std::cout << "electrons in cells" << std::endl;
  //spiral lattice integration
  for(int c = 0; c < total_cells; c++) {
    cell_nearest_neighbor_list.at(c).resize(27); 
    const int x_cell = cell_lattice_coordinate.at(c)[0];
    const int y_cell = cell_lattice_coordinate.at(c)[1];
    const int z_cell = cell_lattice_coordinate.at(c)[2];
   // cell_nearest_neighbor_list.at(c)[0] = c;
      //if(c == 0) std::cout << "starting spiral for cell 0" << std::endl;
   
    for(int s = 0; s < 27; s++) {
    //  if(s == 14) continue;
         // if(s == 0) std::cout << "nn 1 for spiral 0. Expecting <4,4,4>" << std::endl;

      int cell_spiral_xmod = x_cell + (s % 3) - 1;
      if(cell_spiral_xmod < 0) cell_spiral_xmod = x_omp_cells - 1;
      else if (cell_spiral_xmod > x_omp_cells-1) cell_spiral_xmod = 0;
         // if(s == 0) std::cout << "<" << cell_spiral_xmod;

      int cell_spiral_ymod = y_cell + (int(floor(s / 3)) % 3) -1;
      if(cell_spiral_ymod < 0) cell_spiral_ymod = y_omp_cells - 1;
      else if (cell_spiral_ymod > y_omp_cells-1) cell_spiral_ymod = 0;
         // if(s == 0) std::cout << "," << cell_spiral_ymod;

      int cell_spiral_zmod = z_cell + int(floor(s / 9)) - 1;
      if(cell_spiral_zmod < 0) cell_spiral_zmod = z_omp_cells - 1;
      else if (cell_spiral_zmod > z_omp_cells-1) cell_spiral_zmod = 0;
         // if(s == 0) std::cout << "," << cell_spiral_zmod << ">" << std::endl;

      int cell = lattice_cell_coordinate.at(cell_spiral_xmod).at(cell_spiral_ymod).at(cell_spiral_zmod);
      cell_nearest_neighbor_list.at(c).at(s) = cell;
      spiral_cell_counter[cell]++;
      //if(s == 13) std::cout << cell_nearest_neighbor_list.at(c).at(s) << ", " << c << ", " << x_cell << ", " << y_cell << ", " << z_cell << std::endl;
   
    }
  }

  for(int c = 0; c < total_cells; c++) {
    if(spiral_cell_counter[c] != 27) std::cout << "spiral cell error on cell " << c << " with " << spiral_cell_counter[c] << std::endl;
  }
  
  
    if(err::check) std::cout << "spiral integration coordiantes initialized." << std::endl;

    const int max_x_threads = 4;
    const int max_y_threads = 4;
    const int max_z_threads = 4;  

    int max_total_threads = (x_omp_cells/max_x_threads) *(y_omp_cells/ max_y_threads) * (z_omp_cells/ max_z_threads);
   if(max_total_threads != omp_threads) std::cout << "maximum omp threads based on given lattice parameters: " << max_total_threads << "\n Given threads: " << omp_threads << "\n Reducing to max threads" << std::endl;
   omp_threads = int(std::min(double(max_total_threads), double(omp_threads)));

    cells_per_thread = total_cells / omp_threads;
    if (err::check) std::cout << "thread count " << omp_threads << ", " << cells_per_thread << ", " << total_cells << ", " << x_omp_cells << ", " << y_omp_cells << ", " << z_omp_cells << std::endl;
    lattice_cells_per_omp.resize(omp_threads);

    for(int t=0; t< omp_threads;t++){
      lattice_cells_per_omp.at(t).resize(cells_per_thread);
    }
    // int omp_checkerboard_scheme = 0;
    // if(omp_threads == max_x_threads) omp_checkerboard_scheme = 1;
    // if(omp_threads == max_y_threads) omp_checkerboard_scheme = 2;
   // std::cout << x_omp_cells << ", " << y_omp_cells << ", " << z_omp_cells << ", " <<  cells_per_thread << ", " << max_x_threads << ", " << max_y_threads << ", " << max_z_threads << std::endl;

    #pragma omp parallel
    {
    const int thread = omp_get_thread_num();
    for(int l = 0; l < cells_per_thread; l++) {
      // #pragma omp critical
      // std::cout << omp_get_thread_num() << ", " << (omp_get_thread_num()*2)%x_omp_cells + floor(l/z_omp_cells) << ", " << (int(floor(omp_get_thread_num()*2/x_omp_cells)*2)%(y_omp_cells)) + int(floor(l/(z_omp_cells/2)))%2 << ", " << 4*floor(omp_get_thread_num()/(2*y_omp_cells)) + (l % (z_omp_cells/2)) << std::endl;
      
      //double decker cells
    
     lattice_cells_per_omp[thread].at(l) = lattice_cell_coordinate.at((thread*max_x_threads)%x_omp_cells + int(floor(l/(max_y_threads*max_z_threads)))).at((max_y_threads)*(int(floor(thread*max_x_threads/(x_omp_cells)))%(y_omp_cells/ max_y_threads)) + int(floor(l/(max_z_threads)))%(max_y_threads)).at((max_z_threads)*floor(thread/(x_omp_cells*y_omp_cells/(max_x_threads*max_y_threads))) + (l % (max_z_threads)));
      
      //single decker cells
     // lattice_cells_per_omp.at(omp_get_thread_num()).at(l) = lattice_cell_coordinate.at((omp_get_thread_num()*2)%x_omp_cells + floor(l/(2*z_omp_cells))).at((int(floor(omp_get_thread_num()*2/x_omp_cells)*2)%y_omp_cells) + int(floor(l/(z_omp_cells)))%2).at( l % (z_omp_cells));
      
      ////serial cells
     // lattice_cells_per_omp.at(omp_get_thread_num()).at(l) = l;  
    }
    }

    std::vector<int> cell_check;
    cell_check.resize(total_cells, 0);

    for(int l = 0; l < omp_threads; l++) {
      for (int e = 0; e < cells_per_thread; e++) {
          cell_check[lattice_cells_per_omp.at(l).at(e)]++;
      }
    }
    for(int c = 0; c < total_cells-1; c++) {
      if(cell_check[c] == cell_check[c+1] == 1 ) continue;
      std::cout << "unintegrated cell: " << cell_check[c] << ", " << cell_check[c+1] << std::endl;
    }
    for(int c = 0; c < total_cells; c++) {
      if(cell_check[c] == 1 ) continue;
      std::cout << "over-integrated cell: " << cell_check[c]  << std::endl;
    }
        omp_set_dynamic(0);
       omp_set_num_threads(omp_threads);
   
  std::vector<int> x_comp;
  std::vector<int> y_comp;
  std::vector<int> z_comp;
  x_comp.resize(omp_threads);
  y_comp.resize(omp_threads);
  z_comp.resize(omp_threads);

  for (int c = 0; c < cells_per_thread; c++) {
    for(int l = 0; l < omp_threads; l++) {
      int current_cell = lattice_cells_per_omp.at(l).at(c);
      int current_x = cell_lattice_coordinate.at(current_cell)[0];
      int current_y = cell_lattice_coordinate.at(current_cell)[1];
      int current_z = cell_lattice_coordinate.at(current_cell)[2];

      x_comp[l] = current_x;
      y_comp[l] = current_y;
      z_comp[l] = current_z;
    }
    for(int l = 0; l < omp_threads-1; l++) {
      if(abs(x_comp[l] - x_comp[l+1]) < 1 && abs(y_comp[l] - y_comp[l+1]) < 1) std::cout << x_comp[l] << ", " << x_comp[l+1] << ", " << y_comp[l] << ", " << y_comp[l+1] << ", " << l << ", " << l+1 << std::endl;
      if(abs(z_comp[l] - z_comp[l+1]) < 1 && y_comp[l] == y_comp[l+1] && x_comp[l] == x_comp[l+1]) std::cout <<  z_comp[l] << ", " <<  z_comp[l+1] << ", " << l << ", " << l+1 << std::endl;
    }
  }
    if(err::check)  std::cout << " checkerboard omp parallelisation scheme tested" << std::endl;
    // switch (omp_checkerboard_scheme)
    // {
    // case omp_checkerboard_scheme==1:    
    //   int initial = lattice_cell_coordinate.at(omp_get_thread_num() - 1)[0][0];
    //   for(int t = 0; t < cells_per_thread; t++) {
    //     lattice_cells_per_omp.at(t) = lattice_cell_coordinate.at(omp_get_thread_num() - 1).at(t % y_omp_cells).at(floor(t/z_omp_celss))
    //   }
    //   break;
    // default:
    //   break;
    // }
  /*    omp_set_dynamic(0);
       omp_set_num_threads(omp_threads);
    #pragma omp paralell
    {
      int cell = 
      for( int e = 1; e < size; e++) {
    for(int electron = 0; electron < conduction_electrons; electron++) {
      const int array_index = 3*electron;
      const int x_cell = int(floor(electron_position.at(array_index) / x_step_size));
      //if (x_cell < 0 || x_cell > x_omp_cells) std::cout << electron_position.at(array_index) << ", " << x_step_size << ", " << floor(electron_position.at(array_index) / x_step_size) << std::endl;

      const int y_cell = int(floor(electron_position.at(array_index+1) / y_step_size));
      //if (y_cell < 0 || y_cell > y_omp_cells) std::cout << electron_position.at(array_index+1) << ", " << y_step_size << ", " << floor(electron_position.at(array_index+1) / y_step_size) << std::endl;
    
      const int z_cell = int(floor(electron_position.at(array_index+2) / z_step_size));
      //if (z_cell < 0 || z_cell > z_omp_cells) std::cout << electron_position.at(array_index+2) << ", " << z_step_size << ", " << floor(electron_position.at(array_index+2) / z_step_size) << std::endl;
    
      const int c = lattice_cell_coordinate.at(x_cell).at(y_cell).at(z_cell);

      int ee_scattering_count = 2;
      int ee_dos_count = 1;
      int ee_integration_count = 1;
      for(int s = 0; s < 27; s++) {
        const int cell = cell_nearest_neighbor_list.at(c).at(s);

        const int size = cell_integration_lists.at(cell)[0];
       // if(electron == 0) std::cout << size << ", " << cell << ", " << omp_get_thread_num() << std::endl;

        for(int i = 1; i < size; i++) {
            const int array_index_i = 3*cell_integration_lists.at(cell).at(i);
            if (array_index_i == array_index) continue; //no self repulsion

            double x_distance = electron_position.at(array_index)   - electron_position.at(array_index_i);
            double y_distance = electron_position.at(array_index+1) - electron_position.at(array_index_i + 1);
            double z_distance = electron_position.at(array_index+2) - electron_position.at(array_index_i + 2); 

            if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance = x_distance + lattice_width;
            else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance = x_distance - lattice_width;

            if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance = y_distance + lattice_depth;
            else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance = y_distance - lattice_depth;

            if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance = z_distance + lattice_height;
            else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance = z_distance - lattice_height;

            const double length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
                  if (length != length) std::cout << x_distance << ", " << y_distance << ", " << z_distance <<  ", " << \
        electron_position.at(array_index) << ", " << electron_position.at(array_index_i) << ", " << electron_position.at(array_index+1) << ", " << electron_position.at(array_index_i + 1) << ", " <<\
        electron_position.at(array_index+2) << ", " <<  electron_position.at(array_index_i + 2) << std::endl;
            if(length > e_e_integration_cutoff) continue;
            electron_integration_list.at(electron).at(ee_integration_count) = array_index_i;
            ee_integration_count++;
              if(ee_integration_count >= electron_integration_list.at(electron).size() - 3) {std::cout << electron << ", " << ee_integration_count << " > " << electron_integration_list.at(electron).size() << ", " << length << std::endl;
              break; }

            if (length > e_e_neighbor_cutoff) continue;
            electron_nearest_electron_list.at(electron).at(ee_dos_count) = array_index_i/3;
            ee_dos_count++;

              if(ee_dos_count >= electron_nearest_electron_list.at(electron).size() - 3) {std::cout << electron << ", " << ee_dos_count << " > " << electron_nearest_electron_list.at(electron).size() << ", " << length << std::endl;
              break; }

            if(length > e_e_coulomb_cutoff) continue;
            electron_ee_scattering_list.at(electron).at(ee_scattering_count) = array_index_i/3;
              //  if(electron==0)  std::cout << ee_scattering_count << ", " << length << ", " << cell << ", " << array_index_i/3 << std::endl;
            
               //  if(ee_scattering_count >= electron_ee_scattering_list.at(e).size() - 3) {std::cout << e << ", " << ee_scattering_count << " > " << electron_ee_scattering_list.at(e).size() << ", " << length << std::endl;
                // break; }
            ee_scattering_count++;
        }
      }
      electron_integration_list.at(electron)[0] = ee_integration_count;
      electron_nearest_electron_list.at(electron)[0] = ee_dos_count;
      electron_ee_scattering_list.at(electron)[1] = ee_scattering_count;
  //   std::cout << ee_integration_count << std::endl;
  }
    
    std::cout << "nearest neighbor integration list complete." << std::endl;
*/
} 


//====================================
// Creates and outputs atomic lattice
//      Currently static lattice
//====================================
void initialize_positions() {

   // initialize_lattice();

     if (err::check)  std::cout << "Lattice" << std::endl;

    initialize_electrons();

     if (err::check)  std::cout << "Electrons" << std::endl;
}


void initialize_lattice() {
    // atom_anchor_position.resize(lattice_atoms*3,0);
    //   omp_set_dynamic(0);
    //    omp_set_num_threads(omp_threads);
    // #pragma omp parallel for schedule(static) 
    // for (int a = 0; a < lattice_atoms; ++a) {  
    //     int array_index = 3*a;
    //     atom_anchor_position.at(array_index)     = atoms::x_coord_array.at(a) + 0.5*x_unit_size;
    //     atom_anchor_position.at(array_index + 1) = atoms::y_coord_array.at(a) + 0.5*y_unit_size;
    //     atom_anchor_position.at(array_index + 2) = atoms::z_coord_array.at(a) + 0.5*z_unit_size;
        
    // }
}

//====================================
// Function call to initialize basic electron structures
//      Will set up lattice and conduction electrons and velocites
//      Any additional features, e.g. spin, will need to be added here through function calls
//====================================
void initialize_electrons() {

    // char directory .at(omp_threads6);
    // if(getcwd(directory, sizeof(directory)) == NULL){
    //         std::cerr << "Fatal getcwd error in datalog." << std::endl;
    // }
    // electron_position_output_down.open(string(directory) + "/Electron_Position_Init.xyz");
    // electron_position_output_down << conduction_electrons << "\n";  
    // electron_position_output_down << "Initial positions for electrons" << "\n";  

            if (err::check) std::cout << "Lattice output file and electron position file opened..." << std::endl;
    //=========
    // Initialize arrays for holding electron variables
    //      Arrays in super array format to take advantage of caching
    //========
    electron_position.resize(conduction_electrons * 3, 0); // ""'Memory is cheap. Time is expensive' -Steve Jobs; probably" -Michael Scott." -Headcannon.
    electron_velocity.resize(conduction_electrons * 3, 0); //Angstroms
    electron_potential.resize(conduction_electrons, 0);
    ee_dos_hist.resize(conduction_electrons);
    relaxation_time_hist_ee.resize(3*conduction_electrons);
    // relaxation_time_hist_ea.resize(3*conduction_electrons);

    dos_size = int(floor((100.0/phonon_energy)/1.0))+1;
    
    global_e_dos.resize(dos_size);
    // temp_Map.resize(8);
    flux_index.resize(100,0); 
    //const static double step_size = 8.0*((8.0*constants::kB_r*Te) + ((1.0 - 0.9817)*E_f_A)) / double(conduction_electrons);
    // for(int i = 0; i < 8; i++) {
    //     temp_Map.at(i).resize(dos_size+10,0);
    //     // std::cout << temp_Map.at(i)[0] << std::endl;
    //     // std::cout << temp_Map.at(i).at(round(conduction_electrons*0.3)-1) << std::endl;
    // }
    //std::cout << temp_Map.at(7).at(3335) << std::endl;
    e_a_scattering_count = 0;
    e_e_scattering_count = 0;
    ee_transport_scattering_count = 0;
    ee_core_scattering_count = 0;
    ea_transport_scattering_count = 0;
    ea_core_scattering_count = 0;
    ee_scattering_angle = sim::ee_scattering_angle;
    e_e_neighbor_cutoff = pow((lattice_width/8.0)-1.0,2.0);
    
    half_int_var =  2;//(e_e_integration_cutoff - e_e_neighbor_cutoff) / (dt*v_f);
    // full_int_var = 4;//2*half_int_var;
 //   boundary_conditions_cutoff = 18.0; //_e_integration_cutoff - 2;
   // e_e_neighbor_cutoff *= e_e_neighbor_cutoff;
    e_e_integration_cutoff = pow(lattice_width/8.0,2.0);
    e_e_coulomb_cutoff = pow(1.4*1.4*1.4, 2.0);
    
   // std::cout << half_int_var << ", " << full_int_var << ", " << boundary_conditions_cutoff << ", " << e_e_integration_cutoff << std::endl;
    // electron_transport_list.resize(conduction_electrons, false);
    electron_integration_list.resize(conduction_electrons);
    electron_nearest_electron_list.resize(conduction_electrons);
  
    electron_ee_scattering_list.resize(conduction_electrons);
    electron_ea_scattering_list.resize(conduction_electrons);

        if (err::check) std::cout << "Prepare to set position: " << std::endl;
    const int e_density =   int(round(3*int(round(pow(e_e_integration_cutoff,1.5)*1.25*M_PI * 1.0*n_f * 1e-3))));
     ee_density =  3*int(round(pow(e_e_neighbor_cutoff, 1.5)*1.25*M_PI * 1.0*n_f * 1e-3));
    const int ee_scattering = int(6*round(pow(e_e_coulomb_cutoff,   3.0)*1.25*M_PI * 1.0*n_f * 1e-3));
        if (err::check)  std::cout << e_density << ", " << ee_density << ", " << ee_scattering << std::endl;
      
    for(int h = 0; h < dos_size; h++) {
      global_e_dos[h].resize(2,0);
    }

    omp_set_dynamic(0);
    omp_set_num_threads(omp_threads);
    #pragma omp parallel for schedule(static) 
    for (int e = 0; e < conduction_electrons; e++) {

        electron_integration_list.at(e).resize(e_density, 0);
        electron_nearest_electron_list.at(e).resize(ee_density, 0);
        electron_ee_scattering_list.at(e).resize(ee_scattering*2, 0);
        electron_ea_scattering_list.at(e).resize(2,0);
        ee_dos_hist.at(e).resize(dos_size,0);
        
        relaxation_time_hist_ee[3*e].resize(4*100,0);
        relaxation_time_hist_ee[3*e+1].resize(4*100,0);
        relaxation_time_hist_ee[3*e+2].resize(4*100,0);

        // relaxation_time_hist_ea[3*e].resize(4*70,0);
        // relaxation_time_hist_ea[3*e+1].resize(4*70,0);
        // relaxation_time_hist_ea[3*e+2].resize(4*70,0);
        
        const int array_index = 3*e;
        electron_position.at(array_index)     = atoms::x_coord_array.at((e)%lattice_atoms) + 0.5*x_unit_size;
        electron_position.at(array_index + 1) = atoms::y_coord_array.at((e)%lattice_atoms) + 0.5*y_unit_size;
        electron_position.at(array_index + 2) = atoms::z_coord_array.at((e)%lattice_atoms) + 0.5*z_unit_size;
        //initialize and output electron posititons
      //  = atom_anchor_position.at(3*(e%lattice_atoms));//   + cos(theta)*sin(phi)*screening_depth;//*radius_mod(gen)); //Angstroms
       // electron_position.at(array_index + 2) = atom_anchor_position.at(3*(e%lattice_atoms)+2);// + cos(phi)*screening_depth;//*radius_mod(gen);
    }
    // int array_index;
    // for(int e = 0; e < conduction_electrons; e++) {
    //   array_index = 3*e;
    //   electron_position_output_down << "H" << "    " << electron_position.at(array_index)  << "    " << electron_position.at(array_index + 1) << "    " <<  electron_position.at(array_index + 2)  << "\n";    
    // }
    
    // electron_position_output_down.close(); 

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
        electron_num = captured_electron_list.at(2*e);
     //   std::cout << captured_electron_list.at(2*e+1) << std::endl;
        radius = captured_electron_list.at(2*e + 1)*radius_distrib(gen);
        if(radius > 1) radius = 2 - radius;
     //   std::cout << radius << std::endl;
        theta = M_PI * Theta_pos_distrib(gen);
        phi = M_PI * Phi_pos_distrib(gen);

        d_x = cos(theta)*sin(phi)*radius;
        d_y = sin(theta)*sin(phi)*radius;
        d_z = cos(phi)*radius;

        new_electron_position.at(3*electron_num)     = d_x + ((atomic_size * round((new_electron_position.at(electron_num*3)-1) / atomic_size)) + 1); //closest x atom index
        new_electron_position.at(3*electron_num + 1) = d_y + ((atomic_size * round((new_electron_position.at(electron_num*3+1)-1) / atomic_size)) + 1); //closest x atom index;
        new_electron_position.at(3*electron_num + 2) = d_z + ((atomic_size * round((new_electron_position.at(electron_num*3+2)-1) / atomic_size)) + 1); //closest x atom index;
    }

    for(int e = 0; e < captured_electron_list.size()/2; e++) {
        new_potential = 0;
        electron_num = captured_electron_list.at(2*e);
    
        Px = electron_velocity.at(electron_num*3);
        Py = electron_velocity.at(electron_num*3 + 1);
        Pz = electron_velocity.at(electron_num*3 + 2);
        P = (Px*Px)+(Py*Py)+(Pz*Pz);
        
        d_x = new_electron_position.at(electron_num*3) - ((atomic_size * round((new_electron_position.at(electron_num*3)-1) / atomic_size)) + 1); //closest x atom index
        d_y = new_electron_position.at(electron_num*3+1) - ((atomic_size * round((new_electron_position.at(electron_num*3+1)-1) / atomic_size)) + 1); //closest y atom index
        d_z = new_electron_position.at(electron_num*3+2) - ((atomic_size * round((new_electron_position.at(electron_num*3+2)-1) / atomic_size)) + 1); //closest z atom index
        
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

            
            electron_force.at(electron_num)     += force * cos(theta)*sin(phi);
            electron_force.at(electron_num + 1) += force * sin(theta)*sin(phi);
            electron_force.at(electron_num + 2) += force * cos(phi);

        }
        // std::cout << new_potential << std::endl;
        for (int i = 1; i < electron_nearest_electron_list.at(electron_num)[0]+1; i++) {
            if (i == electron_num) continue; //no self repulsion
            //if (symmetry_list.at(e).at(i)) continue;  //make use of symmetry

            //   electron_spin_two = conduction_electron_spin.at(i);
            array_index_i = 3*electron_nearest_electron_list.at(electron_num).at(i);
            x_distance = new_electron_position.at(electron_num*3) - new_electron_position.at(array_index_i);
            y_distance = new_electron_position.at(electron_num*3+1) - new_electron_position.at(array_index_i + 1);
            z_distance = new_electron_position.at(electron_num*3+2) -new_electron_position.at(array_index_i + 2);

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

            
            electron_force.at(electron_num*3)     += force * cos(theta)*sin(phi);
            electron_force.at(electron_num*3 + 1) += force * sin(theta)*sin(phi);
            electron_force.at(electron_num*3 + 2) += force * cos(phi);
        }
       //KE1+PE1 = KE2+PE2
       //v2 = sqrt(2KE2/m_e) = sqrt(2(KE1+PE1-PE2)/m_e)
        P = 1e-6*sqrt(abs(2*((electron_potential.at(electron_num)*constants::K*1e10) + (P*constants::m_e*1e10/2)-(new_potential*constants::K*1e10))/constants::m_e));
      //  std::cout  << electron_potential.at(electron_num)*constants::K*1e10 + ((Px*Px)+(Py*Py)+(Pz*Pz))*constants::m_e*1e10/2 - new_potential*constants::K*1e10 - P*P*constants::m_e*1e10/2 <<  std::endl;
        theta = M_PI * Theta_Vel_distrib(gen);
        phi = M_PI * Phi_Vel_distrib(gen);
        
              //vel = sqrt(2*KE/m_e) =               TE     -     PE
       // if (e ==0 ) std::cout << electron_potential.at(e) << std::endl;
      
      //  if(e == 0) std::cout << "KE: " << 0.5*constants::m_e*v_f*v_f << ", PE: " << electron_potential.at(e)*1e10*constants::K << std::endl;
        electron_velocity.at(electron_num)     = cos(theta)*sin(phi)*P; 
        electron_velocity.at(electron_num + 1) = sin(theta)*sin(phi)*P;
        electron_velocity.at(electron_num + 2) = cos(phi)*P;
       // std::cout << constants::K*1e10*(electron_potential.at(electron_num)) << ", " << ((Px*Px)+(Py*Py)+(Pz*Pz))*constants::m_e*1e10/2 << ", " << new_potential*constants::K*1e10 << ", " << P*P*constants::m_e*1e10/2 << std::endl;
        TLE += constants::K*1e10*electron_potential.at(electron_num)+ ((Px*Px)+(Py*Py)+(Pz*Pz))*constants::m_e*1e10/2 - new_potential*constants::K*1e10 - P*P*constants::m_e*1e10/2;
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

   // initialize_atomic_interactions();
    if (err::check)  std::cout << "Atomic interactions" << std::endl;

   // initialize_electron_atom_interactions();
    if (err::check)  std::cout << "Atomic electron interactions" << std::endl;
}

void initialize_electron_interactions() {
  omp_set_dynamic(0);
       omp_set_num_threads(omp_threads);
    #pragma omp parallel for schedule(static)
    for (int e = 0; e < conduction_electrons; e++) {
      int array_index_i;
     
      double x_distance,y_distance,z_distance, length;
      int array_index = 3*e;
      int integration_count = 1;
      int neighbor_count = 1;
      int scattering_count = 2;
                if (err::check) if(e ==0) std::cout << "Calculating conduction electron repulsion" << std::endl;

        for (int i = 0; i < conduction_electrons; i++) {
            if (i == e) continue; //no self repulsion

            array_index_i = 3*i;
            x_distance = electron_position.at(array_index)     - electron_position.at(array_index_i);
            y_distance = electron_position.at(array_index + 1) - electron_position.at(array_index_i + 1);
            z_distance = electron_position.at(array_index + 2) - electron_position.at(array_index_i + 2);
            
            if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance = x_distance + lattice_width;
            else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance = x_distance - lattice_width;

            if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance = y_distance + lattice_depth;
            else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance = y_distance - lattice_depth;

            if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance = z_distance + lattice_height;
            else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance = z_distance - lattice_height;

            length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);

            if(length > e_e_integration_cutoff) continue;
            electron_integration_list.at(e).at(integration_count) = array_index_i;
            integration_count++;

            if (length > e_e_neighbor_cutoff) continue;
            electron_nearest_electron_list.at(e).at(neighbor_count) = array_index_i/3;
            neighbor_count++;

            if(length > e_e_coulomb_cutoff) {
              
              continue; }
              //std::cout << e << ", " << i << ", " << length << std::endl;
            electron_ee_scattering_list.at(e).at(scattering_count) = array_index_i/3;
            scattering_count++;
        }
        electron_integration_list.at(e)[0] = integration_count;
        electron_nearest_electron_list.at(e)[0] = neighbor_count;
        electron_ee_scattering_list.at(e)[1] = scattering_count;
     }
  }
/*
// void initialize_atomic_interactions() {

//     #pragma omp parallel for schedule(static)
//     for (int e = 0; e < lattice_atoms; e++) {

//       int array_index_i;
//       double x_distance,y_distance,z_distance, length;
//       int array_index = 3*e;

//       int neighbor_count = 2;

//       for (int i = 0; i < lattice_atoms; i++) {
//           if (i == e) continue; //no self repulsion

//         array_index_i = i*3;

//         x_distance = atom_position.at(array_index)     - atom_position.at(array_index_i);
//         y_distance = atom_position.at(array_index + 1) - atom_position.at(array_index_i + 1);
//         z_distance = atom_position.at(array_index + 2) - atom_position.at(array_index_i + 2);

//         if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
//         else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

//         if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
//         else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

//         if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
//         else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

//         length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
           
//         if(length > a_a_neighbor_cutoff) continue;
       
//         atomic_nearest_atom_list.at(e).at(neighbor_count) = i;
//         neighbor_count++;      
//       }
//       atomic_nearest_atom_list.at(e)[1] = neighbor_count; 
    
//             if(err::check) if(e ==0) std::cout << "Calculating conduction-lattice repulsion" << std::endl;
//     }

// }

// void initialize_electron_atom_interactions() { //we'll need a more developed algorithmn for #electrons != #atoms
    

//     #pragma omp parallel for schedule(static)
//     for (int e = 0; e < conduction_electrons; e++) {

//       int array_index, array_index_a;
//       double x_distance,y_distance,z_distance, length;

//       array_index = 3*e;
//       int nearest_electron_count = 1;
      
//         for (int a = 0; a < lattice_atoms; a++) {
            
//             array_index_a = 3*a;
//             x_distance = electron_position.at(array_index)     - atom_position.at(array_index_a);
//             y_distance = electron_position.at(array_index + 1) - atom_position.at(array_index_a+1);
//             z_distance = electron_position.at(array_index + 2) - atom_position.at(array_index_a+2); 

//             if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
//             else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

//             if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
//             else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

//             if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
//             else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

//             length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
      
//             if(length > e_a_neighbor_cutoff) continue;
       
//             atomic_nearest_electron_list.at(e).at(nearest_electron_count) = array_index_a;
//             nearest_electron_count++;
//         }
//         atomic_nearest_electron_list.at(e)[0] = nearest_electron_count;
//     }
// }
*/
void initialize_velocities() {
    
    // char directory .at(omp_threads6);
    // if(getcwd(directory, sizeof(directory)) == NULL){
    //         std::cerr << "Fatal getcwd error in datalog." << std::endl;
    // }
    // electron_velocity_output.open(string(directory) + "/Electron_Velocity/init.csv");
    // electron_velocity_output << "electron number, x-component, y-component, z-component, length" << std::endl;  
    // std::ofstream atom_phonon_output;
    // atom_phonon_output.open(string(directory) + "/Atom_Energy/init.csv");
    // atom_phonon_output << "atom number, energy" << std::endl;

    const std::string n = "Init_E_distrib";
    std::cout << "conduction electrons " << conduction_electrons << std::endl;
   create_fermi_distribution(n, electron_potential,constants::kB_r*Te/E_f_A);
        if (err::check) std::cout << "distribution generated" << std::endl;
  //  const std::string na = "P_distrib";
    //create_phonon_distribution(na, atom_potential,constants::kB_r*Te/E_f_A);
   
      omp_set_dynamic(0);
       omp_set_num_threads(omp_threads);
       p_x = 0.0;
       p_y = 0.0;
       p_z = 0.0;

    // std::srand(std::time(nullptr));
    // std::random_device rd;  //Will be used to obtain a seed for the random number engine
    // std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    // std::uniform_int_distribution<> test_int(0, conduction_electrons -1);
    // std::uniform_real_distribution<double> test_uniform;

    double minimum = E_f_A;
   // transport_cutoff = 106.978;
    
    #pragma omp parallel 
    {
      
      #pragma omp for schedule(dynamic, 8) reduction(+:p_x,p_y,p_z) 
      for(int e = 0; e < conduction_electrons; e++) {
        double phi;
        double theta;
        double sign;
        const int array_index = 3*e;
   
        const double energy = electron_potential[int(omp_uniform_random.at(omp_get_thread_num())()*2147483647)%conduction_electrons];
        //  if(energy > transport_cutoff) electron_transport_list.at(e) = true;
        
        const double vel = sqrt(2.0*energy*constants::m_e_r_i);

        if(sim::CASTLE_x_vector < 0.0 && sim::CASTLE_y_vector < 0.0 && sim::CASTLE_z_vector < 0.0) {
          theta = 2.0*M_PI*omp_uniform_random[omp_get_thread_num()]();
          phi = M_PI*omp_uniform_random[omp_get_thread_num()]();
         
        } else {
          const double unit = sqrt((sim::CASTLE_x_vector*sim::CASTLE_x_vector)+(sim::CASTLE_y_vector*sim::CASTLE_y_vector)+(sim::CASTLE_z_vector*sim::CASTLE_z_vector));
          theta = atan2(sim::CASTLE_y_vector , sim::CASTLE_x_vector);
            if(sim::CASTLE_y_vector == sim::CASTLE_x_vector == 0.0) theta = 0.0;
          phi = acos(sim::CASTLE_z_vector / unit);
         // if(sim::CASTLE_z_vector < 0.0) theta += M_PI;
        }

          if(theta != theta || phi != phi || vel != vel) std::cout << theta << ", " << phi << ", " << vel << ", " << energy << std::endl;
      
            if (err::check) if(e==0) std::cout << "Electron velocity ready..." << std::endl;
        electron_velocity.at(array_index)     = cos(theta)*sin(phi)*vel; 
        electron_velocity.at(array_index + 1) = sin(theta)*sin(phi)*vel;
        electron_velocity.at(array_index + 2) = cos(phi)*vel; 
        

        p_x += electron_velocity.at(array_index);
        p_y += electron_velocity.at(array_index+1);
        p_z += electron_velocity.at(array_index+2);
      }
    }

   // std::cout << count << std::endl;
        if (err::check) std::cout << "distribution assigned" << std::endl;
      std::ofstream Init_E_vel;
      Init_E_vel.open("Init_E_vel");
    #pragma omp parallel for schedule(static) reduction(min:minimum)
    for(int e = 0; e < conduction_electrons; e++) {
      const int array_index = 3*e;
      const double energy = 0.5*constants::m_e_r*(electron_velocity.at(array_index)*electron_velocity.at(array_index) + electron_velocity.at(array_index+1)*electron_velocity.at(array_index+1) + electron_velocity.at(array_index+2)*electron_velocity.at(array_index+2));
      electron_potential.at(e) = energy;
      if (energy < minimum) minimum = energy;
    }

    core_cutoff = minimum;
    // transport_cutoff = core_cutoff + 44.01;

  #pragma omp parallel 
  {
    std::vector<int> local_e_dos;
    local_e_dos.resize(dos_size,0);
    #pragma omp for nowait
    for(int e = 0; e < conduction_electrons; e++) {
      local_e_dos[int(std::min(dos_size-1.0, std::max(0.0, floor((electron_potential[e] - core_cutoff)/phonon_energy))))]++;
    }
    
    #pragma omp critical 
    {
      for(int h = 0; h < dos_size; h++) {
        global_e_dos[h][1] += local_e_dos[h];
      }
    }
  }
    char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog. time stamps" << std::endl;
      }

    std::ofstream global_e_dos_out;
    global_e_dos_out.open(string(directory) + "/Temp_Map/init");
    lattice_output.open(string(directory) + "/ee_scattering");
    int total = 0;
    for(int h = 0 ; h < dos_size; h++) {
      global_e_dos_out << h + int(core_cutoff) << ", " <<  global_e_dos[h][1] << std::endl;
      total += global_e_dos[h][1];
    }

    if(total != conduction_electrons) std::cout <<"hist problem " << total << ", " << conduction_electrons << std::endl;

    std::cout << "core cutoff: " << core_cutoff << ", transport cutoff: " << transport_cutoff << std::endl;
    if (err::check) std::cout << p_x/double(conduction_electrons) << ", " << p_y/double(conduction_electrons) << ", " << p_z/double(conduction_electrons) << std::endl;
    for(int e = 0; e < conduction_electrons; e++) {
      // electron_velocity.at(3*e)   -= (p_x/double(conduction_electrons));
      // electron_velocity.at(3*e+1) -= (p_y/double(conduction_electrons));  
      // electron_velocity.at(3*e+2) -= (p_z/double(conduction_electrons));

      Init_E_vel << e << ", " << electron_potential.at(e)<< ", " <<  electron_velocity.at(3*e) << ", " << electron_velocity.at(3*e+1) << ", " << electron_velocity.at(3*e+2) << "\n";
    }
   if (err::check)  std::cout << "distribution output" << std::endl;
    
            if (err::check) std::cout << "Electron velocity ready..." << std::endl;
}

double M_B_distrib(const double& epsilon, const double& beta) {

  return (exp(epsilon / beta) / (beta*(exp(epsilon / beta) + 1.0)*(exp(epsilon / beta) + 1.0)));
  
}

void create_phonon_distribution(std::vector<double>& distribution, const double& beta) {

//   const double step_size = 1.0 / double(conduction_electrons);
//   const double offset  = beta*3.0;
//   int count = 0;

//  // std::cout << step_size << ", " << offset << std::endl;
//   while(count < conduction_electrons) {
//     if(beta == 0) {
//       distribution.at(count) = E_f_A;
//       count++;
//       continue;
//     }
//     double electron = double(omp_int_random.at(omp_get_thread_num())() % conduction_electrons);
//     double epsilon = (step_size *electron)  - offset;
//     if(omp_uniform_random.at(omp_get_thread_num())() < ((0.01 / 16.0)*(epsilon+(3.0*beta))*(epsilon+(3.0*beta))*exp(-0.5*(epsilon+(3.0*beta)) / beta) / (beta*beta*beta))) {
//    //  if(E_f_A*(epsilon+1.0) < E_f_A - 3.0*constants::kB_r*Tp) std::cout << E_f_A*(epsilon+1.0) << ", " << E_f_A - (3.0*constants::kB_r*Tp) << ", " << epsilon << ", " << beta  << std::endl;
//       distribution.at(count) = E_f_A*(epsilon + 1.0);
//       count++;
//     }
//   }
}

void create_phonon_distribution(const std::string& name, std::vector<double>& distribution, const double& beta) {

  // char directory [256];
  //     if(getcwd(directory, sizeof(directory)) == NULL){
  //           std::cerr << "Fatal getcwd error in datalog." << std::endl;
  //     }
  // std::ofstream distrib;
  // distrib.open(string(directory) +"/"+ name);
  // distrib.precision(20);

  // const double step_size = 1.0 / double(conduction_electrons);
//   const double offset  = beta*3.0;
//   int count = 0;

//  // std::cout << step_size << ", " << offset << std::endl;
//   while(count < conduction_electrons) {
//     if(beta == 0) {
//       distribution.at(count) = E_f_A;
//       count++;
//       continue;
//     }
//     double electron = double(omp_int_random.at(omp_get_thread_num())() % conduction_electrons);
//     double epsilon = (step_size *electron)  - offset;
//     if(omp_uniform_random.at(omp_get_thread_num())() < ((0.01 / 16.0)*(epsilon+(3.0*beta))*(epsilon+(3.0*beta))*exp(-0.5*(epsilon+(3.0*beta)) / beta) / (beta*beta*beta))) {
//      // if(E_f_A*(epsilon+1.0) < E_f_A - 3.0*constants::kB_r*Te) std::cout << E_f_A*(epsilon+1.0) << ", " << E_f_A - (3.0*constants::kB_r*Te) << ", " << epsilon << ", " << beta  << std::endl;
//       distribution.at(count) = E_f_A*(epsilon + 1.0);
//       distrib << count << ", " << distribution.at(count) << "\n";
//       count++;
//     }
//   }
//   distrib.close();
}

void create_fermi_distribution(const std::string& name, std::vector<double>& distribution, const double beta) {

  char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog. fermi dist" << std::endl;
      }
  std::ofstream distrib;
  distrib.open(string(directory) +"/"+ name);
  distrib.precision(20);

  double step_size = 0.002;
  const double max  = E_f_A + 3.0*constants::kB_r*Te;
  const double min = E_f_A - 3.0*constants::eV_to_AJ;
  if (err::check)  std::cout << "min: " << min << ", max: " << max << std::endl;
  step_size *= (max-min);
 // else step_size = 1.0;

  int count = 0;
  int subCount = 0;
  int energy_step = 0;
  if(beta <= 0.000001 ) {
    while(count < conduction_electrons) {
      distribution.at(count) = E_f_A;
      distrib << count << ", " << distribution.at(count) << "\n";
      count++;
    }
  } else {
    while(count < conduction_electrons) { 

      double epsilon = max - step_size*energy_step;
      int occupation = int(round((conduction_electrons/(max-min))*0.1*(return_fermi_distribution((epsilon-E_f_A)/E_f_A, beta))));//+return_fermi_distribution((epsilon - 0.5*step_size - E_f_A)/E_f_A, beta))));
     
     // int steps = round(1.0 / double(occupation));
      for(int o = 0; o < occupation; o++) {
        distribution.at(count) = epsilon - 0.5*step_size;
        distrib << count << ", " << distribution.at(count) << "\n";
        count++;
      
        if(count == conduction_electrons) break;
      }
      if(1 - return_fermi_distribution((epsilon-E_f_A)/E_f_A, beta) > 1e-3) transport_cutoff = epsilon;
      energy_step++;
    }
  }

 // std::cout << "Total atoms to fill " << count << " electrons: " << subCount << std::endl;
  distrib.close();
}
/*
void create_gaussian_distribution(const std::string& name, std::vector<double>& distribution, const double& beta) {

  char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog. fermi dist" << std::endl;
      }
  std::ofstream distrib;
  distrib.open(string(directory) +"/"+ name);
  distrib.precision(20);

  double step_size = 0.002;
  const double max  = E_f_A + 3.0*constants::kB_r*Te;
  const double min = E_f_A - 3.0*constants::eV_to_AJ;
  std::cout << "min: " << min << ", max: " << max << std::endl;
  step_size *= (max-min);
 // else step_size = 1.0;

  int count = 0;
  int subCount = 0;
  int energy_step = 0;
  if(beta <= 0.000001 ) {
    while(count < conduction_electrons) {
      distribution.at(count) = E_f_A;
      distrib << count << ", " << distribution.at(count) << "\n";
      count++;
    }
  } else {
    while(count < conduction_electrons) { 

      double epsilon = max - step_size*energy_step;
      int occupation = int(round((conduction_electrons/(max-min))*0.05*(return_fermi_distribution((epsilon-E_f_A)/E_f_A, beta)+return_fermi_distribution((epsilon - 0.5*step_size - E_f_A)/E_f_A, beta))));
     
     // int steps = round(1.0 / double(occupation));
      for(int o = 0; o < occupation; o++) {
        distribution.at(count) = epsilon - 0.25*step_size;
        distrib << count << ", " << distribution.at(count) << "\n";
        count++;
        //if (count == conduction_electrons/2) transport_cutoff = epsilon - step_size*0.5;
 
        if(count == conduction_electrons) break;
      }
      if(1 - return_fermi_distribution((epsilon-E_f_A)/E_f_A, beta) > 1e-4) transport_cutoff = epsilon;
      energy_step++;
     // if (epsilon <= min) break;

   // subCount++;
    }
  }
  //transport_cutoff = 173.78;
  //transport_cutoff = core_cutoff;
  
 // std::cout << "Total atoms to fill " << count << " electrons: " << subCount << std::endl;
  distrib.close();
} */
double return_fermi_distribution(const double epsilon, const double beta) 
{
  if(beta <= 0.00001) return 1.0;
  else return (1.0/(exp(epsilon/beta) + 1.0));
}

// double return_gaussian_distribution(const double epsilon, const double beta) 
// {
//   if(beta <= 0.00001) return 1.0;
//   else return exp(-0.5*((epsilon*epsilon) - E_f_A)/(beta*beta));
// }

double return_BE_integrand(const double phonon_e, const double temperature) {
  return 1.0/(exp(phonon_e/(temperature*constants::kB_r)) - 1.0);
}

void output_data() {
  
    //=========
    // Output equilibration step data
    //=========
    double scat_size = 0.0;
    double scat_stddev = 0.0;
    double e_stddev = 0.0;
    double e_size = 0.0;

    omp_set_dynamic(0);
    omp_set_num_threads(omp_threads);
    double global_d_U = 0.0;
    double local_d_U = 0.0;
    int electron_counter = 0;
    #pragma omp parallel for reduction(+:e_size, scat_size, global_d_U, local_d_U, electron_counter)
    for(int e = 0; e < conduction_electrons; e++) {
      double e_energy = electron_potential[e];
      int index = std::min(dos_size-1.0, std::max(0.0, floor((e_energy-core_cutoff)/phonon_energy)));
      if(e_energy > 0.8*E_f_A) {
        e_size += electron_nearest_electron_list[e][0];
        scat_size += electron_ee_scattering_list[e][1];
        electron_counter++;
      }
      if(e_energy < E_f_A) {
        if(e_energy > transport_cutoff) {
          local_d_U += (E_f_A - e_energy)*std::max(0.0, 1.0 - double(ee_dos_hist[e][index])/local_dos_occ);
          global_d_U+= (E_f_A - e_energy)*std::max(0.0, 1.0 - (double(global_e_dos[index][0])/dos_occ)); 
        } else {
          local_d_U += (E_f_A - e_energy)*std::max(0.0, 1.0 - (double(global_e_dos[index][0])/double(global_e_dos[index][1])));
          global_d_U+= (E_f_A - e_energy)*std::max(0.0, 1.0 - (double(global_e_dos[index][0])/double(global_e_dos[index][1])));  
        }
      } else {
        global_d_U += (e_energy - E_f_A)*double(global_e_dos[index][0])/dos_occ;
        local_d_U += (e_energy - E_f_A)*double(ee_dos_hist[e][index])/local_dos_occ;
      }
    }
    e_size /= electron_counter;
    scat_size /= electron_counter;

    #pragma omp parallel for reduction(+:e_stddev,scat_stddev)
    for(int e = 0; e < conduction_electrons; e++) {
      if(electron_potential[e] > 0.8*E_f_A){
        e_stddev += (electron_nearest_electron_list[e][0] - e_size)   *(electron_nearest_electron_list[e][0] - e_size);
        scat_stddev += (electron_ee_scattering_list[e][1] - scat_size)*(electron_ee_scattering_list[e][1] - scat_size);
      }
    }
   
    e_stddev = sqrt(e_stddev/electron_counter);
    scat_stddev = sqrt(scat_stddev/electron_counter);

    if((current_time_step % (CASTLE_output_rate * CASTLE_MD_rate)) == 0) {
      time_stamp = std::to_string(current_time_step);
    
      char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog. time stamps" << std::endl;
      }

      std::ofstream temp_map;
      std::ofstream flux_hist;
      std::ofstream relaxation_time;
    
      temp_map.open(string(directory) + "/Temp_Map/" + time_stamp);
      flux_hist.open(string(directory) + "/flux_hist/" + time_stamp );   
      relaxation_time.open(string(directory)  + "/relaxation_time/" + time_stamp);

      for(int e = 0; e < flux_index.size(); e++) {
        flux_hist << e << ", " << flux_index[e] << "\n";
        flux_index[e] = 0;
      }

    double ee_avg;
    double ea_avg;
    int ee_total;
    int ea_total;
   
    // std::cout << relaxation_time_hist_ee[0].size() << std::endl;
    for(int h = 0; h < relaxation_time_hist_ee[0].size(); h++) {
      ee_avg = 0.0;
      // ea_avg = 0.0;
      ee_total = 0;
      // ea_total = 0;
      for(int e = 0; e < conduction_electrons; e++) { 
        ee_avg += double(relaxation_time_hist_ee[3*e + 2][h])/std::max(1.0, double(relaxation_time_hist_ee[3*e][h]) );
        // ea_avg += double(relaxation_time_hist_ea[3*e + 2][h])/std::max(1.0, double(relaxation_time_hist_ea[3*e][h]) );
        ee_total += relaxation_time_hist_ee[3*e][h];
        // if(relaxation_time_hist_ea[3*e + 2][h] > 0) ea_total++;
      }

      relaxation_time << double(h)/4.0 + int(round(core_cutoff)) << ", " << ee_avg/std::max(1.0,double(ee_total)) << ", " << ee_total  << "\n";
      // << ", " << ea_avg/std::max(1.0,double(ea_total)) << ", " << ea_total
    }
    // }
  //  const int output_count_lr = int(round(transport_cutoff-core_cutoff));
    const int output_count_hr = ee_dos_hist[0].size();
    // #pragma omp parallel
    // {
    // const int thread = omp_get_thread_num();
    // for(int e = 1; e < electron_nearest_electron_list[thread*10000][0]; e++) {
    //   const   int electron = electron_nearest_electron_list[thread*10000][e];
    //   const double energy = electron_potential.at(electron);
    //   if(energy < transport_cutoff) {
    //     const   int hist = int(std::min(double(output_count_lr-1), std::max(0.0, floor((energy- core_cutoff)/step_size_lr))));
    //     temp_Map[thread].at(hist)++;
    //   } else {
    //     const   int hist = int(std::min(double(output_count_hr-1), std::max(0.0, floor((energy- transport_cutoff)/step_size_hr))));
    //     temp_Map[thread+4].at(hist)++;
    //   }
    //  // if(x_pos < (lattice_width * 0.5) && y_pos < (lattice_depth*0.5) && z_pos < (lattice_height * 0.5)) {        
    //        // double x_pos = electron_position.at(array_index);
    //   // double y_pos = electron_position.at(array_index+1); 
    //   // double z_pos = electron_position.at(array_index+2);
    //  // }
    //   // if(x_pos > (lattice_width * 0.5) && y_pos < lattice_depth*0.5 && z_pos < lattice_height * 0.5) {
    //   //    #pragma omp atomic update
    //   //   temp_Map[1].at(hist)++;
    //   // }
    //   // if(x_pos < (lattice_width * 0.5) && y_pos > lattice_depth*0.5 && z_pos < lattice_height * 0.5) {
    //   //    #pragma omp atomic update
    //   //   temp_Map[2].at(hist)++;
    //   // }
    //   // if(x_pos > lattice_width * 0.5 && y_pos > lattice_depth*0.5 && z_pos < lattice_height * 0.5) {
    //   //    #pragma omp atomic update
    //   //   temp_Map[3].at(hist)++;
    //   // }
    //   // if(x_pos < lattice_width * 0.5 && y_pos < lattice_depth*0.5 && z_pos > lattice_height * 0.5) {
    //   //    #pragma omp atomic update
    //   //   temp_Map.at(4).at(hist)++;
    //   // }
    //   // if(x_pos > (lattice_width * 0.5) && y_pos < lattice_depth*0.5 && z_pos > lattice_height * 0.5) {
    //   //   #pragma omp atomic update
    //   //   temp_Map.at(5).at(hist)++;
    //   // }
    //   // if(x_pos < (lattice_width * 0.5) && y_pos > lattice_depth*0.5 && z_pos > lattice_height * 0.5) {
    //   //    #pragma omp atomic update
    //   //   temp_Map.at(6).at(hist)++;
    //   // }
    //   // if(x_pos > (lattice_width * 0.5) && y_pos > lattice_depth*0.5 && z_pos > lattice_height * 0.5) {
    //   //    #pragma omp atomic update
    //   //   temp_Map.at(7).at(hist)++; 
    //   // }    
    // }
    // }  
    
    int count = 0;
    int electrons[4];
    while(count < 3) {
      int selection = int(omp_uniform_random[0]()*2147483647) % (conduction_electrons);
       if(electron_potential[selection] > E_f_A*0.8) {electrons[count] = selection; count++;}
    }
    
    for(int i = 0; i < output_count_hr; i++) {
     // if(i == 11) temp_map_0 << i << ", " << temp_Map[0].at(i) << "\n";
     // if(i < output_count_lr) { 
      temp_map << i*phonon_energy + int(round(core_cutoff)) << ", " << global_e_dos.at(i)[0] << ", " << global_e_dos[i][1] 
                      << ", " << ee_dos_hist[electrons[0]].at(i) \
                      << ", " << ee_dos_hist[electrons[1]].at(i) \
                      << ", " << ee_dos_hist[electrons[2]].at(i) \
                      << ", " << ee_dos_hist[electrons[3]].at(i) << "\n";
      // temp_map_0 << 4*i+1 + int(round(core_cutoff)) << ", " << ee_dos_hist[electrons[0]].at(i) \
      //                 << ", " << ee_dos_hist[electrons[1]].at(i) \
      //                 << ", " << ee_dos_hist[electrons[2]].at(i) \
      //                 << ", " << ee_dos_hist[electrons[3]].at(i) << "\n";
      // temp_map_0 << 4*i +2 + int(round(core_cutoff))<< ", " << ee_dos_hist[electrons[0]].at(i) \
      //                 << ", " << ee_dos_hist[electrons[1]].at(i) \
      //                 << ", " << ee_dos_hist[electrons[2]].at(i) \
      //                 << ", " << ee_dos_hist[electrons[3]].at(i) << "\n";
      // temp_map_0 << 4*i +3 + int(round(core_cutoff))<< ", " << ee_dos_hist[electrons[0]].at(i) \
      //                 << ", " << ee_dos_hist[electrons[1]].at(i) \
      //                 << ", " << ee_dos_hist[electrons[2]].at(i) \
      //                 << ", " << ee_dos_hist[electrons[3]].at(i) << "\n";
      // } else {
      // temp_map_0 << i + int(round(core_cutoff))
      //   << ", " << ee_dos_hist[electrons[0]].at(i) \
      //   << ", " << ee_dos_hist[electrons[1]].at(i) \
      //   << ", " << ee_dos_hist[electrons[2]].at(i) \
      //   << ", " << ee_dos_hist[electrons[3]].at(i) << "\n"; 
      // }
    }
    temp_map.close();

   /* for(int i = 0; i < output_count_lr; i++) {
      //if(i == 11) temp_map_1 << i << ", " << temp_Map[1].at(i) << "\n";
       temp_map_1 << i << ", " << ee_dos_hist[electrons[1]].at(i) << "\n";
     
    }
    temp_map_1.close();

    for(int i = 0; i < output_count_lr; i++) {
     // if(i == 11) temp_map_2 << i << ", " << temp_Map[2].at(i) << "\n";
       temp_map_2 << i << ", " << ee_dos_hist[electrons[2]].at(i) << "\n";
     
    }
    temp_map_2.close();
    
    for(int i = 0; i < output_count_lr; i++) {
     // if(i == 11) temp_map_3 << i << ", " << temp_Map[3].at(i) << "\n";
       temp_map_3 << i << ", " << ee_dos_hist[electrons[3]].at(i) << "\n";
      
    }
    temp_map_3.close();
    */
   
   /* for(int i = output_count_lr; i < output_count_hr; i++) {
    //  if(i == 11) temp_map_4 << i << ", " << temp_Map.at(4).at(i) << "\n";
       temp_map_4 << i-output_count_lr << ", " << ee_dos_hist.at(electrons[0]).at(i) << "\n";
     
    } 
    temp_map_4.close();

    for(int i = output_count_lr; i < output_count_hr; i++) {
      temp_map_5<< i-output_count_lr << ", " << ee_dos_hist.at(electrons[1]).at(i) << "\n";
  
    }
    temp_map_5.close();

    for(int i = output_count_lr; i < output_count_hr; i++) {
      temp_map_6 << i-output_count_lr << ", " << ee_dos_hist.at(electrons[2]).at(i) << "\n";
     
    }
    temp_map_6.close();

    for(int i = output_count_lr; i < output_count_hr; i++) {
      temp_map_7 << i-output_count_lr << ", " << ee_dos_hist.at(electrons[3]).at(i) << "\n";
 
    }
    temp_map_7.close(); */
      // std::ofstream E_vel;
      // E_vel.open("velocity/"+time_stamp);
      // for(int e = 0; e<conduction_electrons; e++) {
      //   E_vel <<  e << ", " << electron_potential.at(e) << ", " << electron_velocity.at(3*e) << ", " << electron_velocity.at(3*e+1) << ", " << electron_velocity.at(3*e+2) << "\n";
      // }
      // E_vel.close();
      // std::ofstream E_pos;
      // E_pos.open("position/"+time_stamp);
      // for(int e = 0; e<conduction_electrons; e++) {
      //   E_pos <<  e << ", " << electron_potential.at(e) << ", " << electron_position.at(3*e) << ", " << electron_position.at(3*e+1) << ", " << electron_position.at(3*e+2) << "\n";
      // }
      // E_pos.close();

      std::cout << std::fixed; std::cout.precision(3); std::cout << "  " << current_time_step / total_time_steps * 100 << "%. " << std::endl; 
    }
    mean_data.precision(10);
    mean_data << std::scientific;

    const double I = double(x_flux) * constants::e * 1e35 / double(CASTLE_output_rate) / dt / (lattice_height*lattice_depth); //current density/m**2

    if(!current_time_step) {
    mean_data << CASTLE_real_time << ", " << current_time_step << ", " 
      << d_Te*d_Te*e_heat_capacity << ", "  << global_d_U << ", " << local_d_U << ", "// (d_Te*d_Te*e_heat_capacity +Tp*a_heat_capacity) << ", " 
      << Te << ", " << Tp << ", " << TEKE << ", " << TLE << ", " 
      << d_TTMe << ", " << d_TTMp << ", " <<  I*double(CASTLE_output_rate) << ", " << e_size << ", " << e_stddev << ", " << scat_size << ", " << scat_stddev << ", " << p_x/double(conduction_electrons) << ", " << p_y/double(conduction_electrons) << ", " << p_z/double(conduction_electrons) << ", " 
      << std::fixed; mean_data.precision(1); mean_data << double(e_a_scattering_count) / 1 << ", " << double(e_e_scattering_count) / double(1) << ", " << \
      double(ee_core_scattering_count) / double(1) << ", " << double(ee_transport_scattering_count) / double(1) << ", " <<\
      double(ea_core_scattering_count) / double(1) << ", " << double(ea_transport_scattering_count) / double(1) << ", " <<\
      double(x_flux) / double(1) << ", " << double(y_flux) / 1 << ", " << double(z_flux) / double(1)  << ", " \
      << std::endl;
    } else {
    mean_data << CASTLE_real_time << ", " << current_time_step << ", " 
      << d_Te*d_Te*e_heat_capacity << ", "  << global_d_U << ", " << local_d_U << ", " //<<  (d_Te*d_Te*e_heat_capacity +Tp*a_heat_capacity) << ", " 
      << d_Te << ", " << d_Tp << ", " << TEKE << ", " << TLE << ", " 
      << d_TTMe << ", " << d_TTMp << ", " <<  I << ", " << e_size << ", " << e_stddev << ", " << scat_size << ", " << scat_stddev << ", " << p_x/double(conduction_electrons) << ", " << p_y/double(conduction_electrons) << ", " << p_z/double(conduction_electrons) << ", " 
      << std::fixed; mean_data.precision(1); mean_data << double(e_a_scattering_count) / CASTLE_output_rate << ", " << double(e_e_scattering_count) / double(CASTLE_output_rate) << ", " << \
      double(ee_core_scattering_count) / double(CASTLE_output_rate) << ", " << double(ee_transport_scattering_count) / double(CASTLE_output_rate) << ", " <<\
      double(ea_core_scattering_count) / double(CASTLE_output_rate) << ", " << double(ea_transport_scattering_count) / double(CASTLE_output_rate) << ", " <<\
      double(x_flux) / double(CASTLE_output_rate) << ", " << double(y_flux) / CASTLE_output_rate << ", " << double(z_flux) / double(CASTLE_output_rate)  << ", " \
      << std::endl;
    }
  //  << std::accumulate(electron_potential.begin(), electron_potential.end(), 0.0) << ", "
    x_flux = 0;
    y_flux = 0;
    z_flux = 0;
    e_a_scattering_count = 0;
    ee_core_scattering_count = 0;
    ee_transport_scattering_count = 0;
    ea_core_scattering_count = 0;
    ea_transport_scattering_count = 0;
    e_e_scattering_count = 0;
    if(transport_cutoff > core_cutoff+45.61 - 0.5*floor((d_Te - 300.0)/100.0)) std::cout << "transport cutoff shift from " << transport_cutoff << " to " << core_cutoff+45.61 - 0.5*floor((d_Te - 300.0)/100.0) << std::endl;
    transport_cutoff = core_cutoff+45.61 - 0.5*floor((d_Te - 300.0)/100.0);
    
    // if(Te > 900) transport_cutoff = 102.951;
  // if(current_time_step > (sim::equilibration_time/2)) ee_rate = -1.0*dt*sim::ee_coupling_strength/(constants::eV_to_AJ*constants::eV_to_AJ);
   // a_a_scattering_count = 0;
}




} //end of CASTLE namespace



