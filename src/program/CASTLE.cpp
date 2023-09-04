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
    lattice_atoms = atoms::num_atoms; //Better lattice creation will come from VAMPIRE in future
    phonon_energy = 1e-3*sqrt(sim::ea_coupling_strength/0.12)*constants::eV_to_AJ;
    dos_en_step = phonon_energy/4.0;
    i_dos_en_step = 1.0/dos_en_step;
     // meV [e-3] AJ
    i_phonon_energy = 1.0/phonon_energy;
       
    // std::cout << "phonon energy " << phonon_energy << std::endl;
    min_as = sim::CASTLE_min_as*constants::eV_to_AJ;
    max_as = sim::CASTLE_max_as*constants::eV_to_AJ;
    DoS_cutoff = dos_en_step* int(floor(-1.0*min_as*i_dos_en_step)); // E_f_A - min
    dos_size = int(floor(((max_as-min_as)*i_dos_en_step)/1.0))+1;

     E_f = 8.0*constants::eV_to_AJ*1e-20;// constants::h * constants::h * pow(3.0 * M_PI * M_PI * n_f*1e27, 0.66666666666666666666666667) / (8.0 * M_PI * M_PI * constants::m_e); //Fermi-energy
    E_f_A =  E_f*1e20; //AJ
      mu_f = 0.00; //mu adjust at 300K
    lattice_width = cs::system_dimensions[0]+x_unit_size; //A
    lattice_depth  = cs::system_dimensions[1]+y_unit_size; // A
    lattice_height  = cs::system_dimensions[2]+z_unit_size; // A
    lattice_volume = lattice_depth*lattice_height*lattice_width;
    cell_volume = x_unit_size*y_unit_size*z_unit_size;

    double l_f = lattice_atoms/lattice_volume; //atoms/A^3
           std::cout << "Initialising lattice: " << lattice_depth*lattice_height*lattice_width*1e-3 << " nm^3; active space range (AJ): " <<E_f_A+ min_as << " <> " << E_f_A+ max_as << ", (" << max_as-min_as << " AJ). Resultant resolution: " << dos_en_step << " for " << dos_size << " bins" << std::endl;
    
    dt = mp::dt_SI * 1e15;//-4; //S -> femptoSeconds
    TTMe = TTMp = d_TTMe = d_TTMp = Tp = Te = d_Tp = d_Te = sim::temperature;

         if(err::check) std::cout << "phonon occupation test and DoS offset discretisation: " << return_BE_distribution(phonon_energy, Te) << ", " << int(floor(1.5*constants::eV_to_AJ*i_dos_en_step)) << std::endl;

    
      char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog. fermi dist" << std::endl;
      }
        // std::cout << min_as << "2< " << max_as << ", " << DoS_cutoff << ", " << dos_size << std::endl;
    std::string dos_name = "Ag-DoS-sorted.csv";
    std::vector<std::vector< double> > dos_scale;
    int dos_intput_size = 281;
    dos_scale.resize(dos_intput_size);
    std::fstream dos_file (dos_name, std::ios::in);
    std::string line;

      double value;
      if(dos_file.is_open()) {
        for (int l = 0; l < dos_intput_size; l++) {

          dos_scale[l].resize(2, 0.0);
          std::getline(dos_file, line);
          std::stringstream s(line);
          std::string s_value;
          int count = 0;
          while(std::getline(s, s_value, ',')) {
             if(count == 0) dos_scale[l][count] = (std::stod(s_value))*constants::eV_to_AJ;
             else dos_scale[l][count] = std::stod(s_value)/constants::eV_to_AJ;// e48 (e-/J/m^3) -> [1/e20/e30] -> e-/AJ/A^3
             count++;
          }
        }
      } else std::cout << dos_name << " did not open" << std::endl;
     
      double min = -1.0*DoS_cutoff;
      
      dos_standard.resize(dos_size, dos_scale[dos_intput_size-1][1]); //states/eV/ -> eV/AJ m^3 -> states/AJ/atom
     // int l = int(ceil(min - dos_scale[0][0]));
      // if(l < 0) std::cout << min << ", " << dos_scale[0][0] << std::endl;
     // std::cout << l << ", " << min << ", " << dos_scale[l][0] << ", " << min - dos_scale[0][0] << std::endl;
      int d_count = 0;
      
      for(int l=0; l < dos_intput_size && d_count < (dos_size-2); l++) {
        if(dos_scale[l][0] < min) continue;
        int count = 0;
        double avg = 0.0;
      
        while(dos_scale[l][0] > min && d_count < (dos_size-2) ) {
          avg += dos_scale[l][1]; //e-/eV/m^3
          min += dos_en_step;
          count++;
          d_count++;
          // avg /= count;
          dos_standard.at(d_count-1) = avg / count; //avg e-/AJ/atom
          // std::cout << d_count-1 << ", " << avg / count << std::endl;
        }
        if(count == 0)  {
          dos_standard.at(d_count) = dos_scale[l][1]; //e-/AJ/atom
          min += dos_en_step;
          d_count++;
        }
      }
      dos_file.close();
      std::ofstream dos_standard_output;
      dos_standard_output.open(string(directory) + "/dos_standard.csv");
      if(!dos_standard_output.is_open()) std::cout << "dos_standard_output.csv is not open" << std::endl;
        // std::cout << min_as << "< " << max_as << ", " << DoS_cutoff << ", " << dos_size << std::endl;
      // min = DoS_cutoff;
      
      for(int d = 1; d < dos_size-2; d++) {
        // dos_standard[d] = return_fermi_distribution(d*dos_en_step+87.5561-E_f_A, Te)*
        
        // if(d*dos_en_step -DoS_cutoff >  dos_scale[0][0]) 
        total_e_scaling += dos_en_step*dos_standard[d]*return_fermi_distribution(d*dos_en_step-DoS_cutoff, 0);//+4.0*return_fermi_distribution(d*dos_en_step-DoS_cutoff + 0.5*dos_en_step , 0) + return_fermi_distribution((d+1)*dos_en_step-DoS_cutoff, 0))/6.0; //e-/A^3 for Fermi window
        dos_standard[d] = lattice_atoms*dos_standard[d]; //e-/AJ for Fermi window for whole lattice 
        // std::cout << d*dos_en_step-DoS_cutoff << ", " << dos_standard[d]*dos_en_step << ", " << 0.5*(dos_standard[d+1]-dos_standard[d-1])/lattice_atoms/dos_en_step << '\n';
        dos_standard_output << d*dos_en_step-DoS_cutoff << ", " << dos_standard[d]*dos_en_step << ", " << (dos_standard[d+1]-dos_standard[d-1]/lattice_atoms)/dos_en_step << '\n';
      }
      dos_standard_output.close();
      conduction_electrons = int(round(total_e_scaling*lattice_atoms));

      // total_e /= constants::eV_to_AJ;
    // DoS_cutoff *= 1.01;
    
     std::cout << "total e- dos: " << total_e_scaling << " (e-/Fermi window/atom); full lattice: " << conduction_electrons << " @ 0K from " << -DoS_cutoff << " to 0 " << std::endl;

    CASTLE_output_rate = sim::partial_time;
    
    // d_TTMp = TTMp = Tp = d_Tp = 300.0;
    
    // Te += 1.0;
    total_time_steps = sim::equilibration_time; //100
    CASTLE_MD_rate = sim::CASTLE_MD_rate;
  
    applied_voltage_sim = sim::applied_voltage_sim;
    heat_pulse_sim = sim::heat_pulse_sim;

    photon_energy = sim::photon_energy * constants::eV_to_AJ;

    ee_coupling = sim::ee_coupling;
    ea_coupling = sim::ea_coupling;

    //========
    // Initialzie variables used by all functions
    //========
    
    x_flux = 0;
    y_flux = 0;
    z_flux = 0;
    current_time_step = 0;
    CASTLE_real_time = 0;
    // uniform_random.seed(2137082040);
    // int_random.seed(2137082040);

    mu_r = (atomic_mass + constants::m_e_r) / (atomic_mass * constants::m_e_r );
    combined_mass = 1 / (atomic_mass + constants::m_e_r);
    applied_voltage = sim::applied_voltage;
    n_f = 1.0e3 * conduction_electrons / (lattice_width * lattice_height * lattice_depth)/total_e_scaling; // e- / A**3 -> atoms/nm**3
   
  
    v_f = sqrt(2.0 * E_f_A * constants::m_e_r_i); // A/fs
 
    a_specific_heat = 25.0 /6.02e3; // //sim::TTCl; //J/K/mol -> .at(e20/Na) AJ/K/e-   
    a_specific_heat_i = 1.0 / a_specific_heat;
    a_heat_capacity = 1.14*a_specific_heat * n_f; //AJ/K/particle .at() AJ/K/nm**3
    a_heat_capacity_i = 1.0 / a_heat_capacity;

    e_specific_heat = 0.5*M_PI*M_PI*constants::kB*constants::kB/E_f; // gamma; //J/K**2/e- 
    e_specific_heat_i = 1.0 / e_specific_heat;
    e_heat_capacity = 63.3*1e-7;// (13.5)*(1.0+0.084)*1e20*e_specific_heat * n_f; //J/K/e- -> AJ/K**2/nm**3 -> [e27/e20] -> J/K**2/m**3
    e_heat_capacity_i = 1.0 / e_heat_capacity;

    ea_coupling_strength = 1e-6*sim::ea_coupling_strength*constants::eV_to_AJ*constants::eV_to_AJ/constants::hbar_r; // meV**2 -> AJ/fs

    atomic_mass = 58.69 * 1.6726219e3; // kg * 1e30 for reduced units
    fluence = 1e3*sim::fluence; // mJ/cm**2 -> 1e4/1e3 J/m^2 -> 1e20/1e18 AJ/nm**2 
    if(photon_energy == 0.0) std::cout << "sim error: photon energy gives infinite intensity " << std::endl;
    const double tau = 3.0*E_f_A /(M_PI*ea_coupling_strength); // fs/AJ
    G = 300.0*e_heat_capacity*E_f_A*3.0/tau; //AJ/fs/K/nm**3 [e-20*e27*e15 = e22]  
    //G = sim::TTG*1e-23;
    //G=Ce/Te-p = pihbar constant (meV**2)Ef*n_f*(1/eps)**2
    ea_rate = -30.0*dt*E_f_A/tau;  //AJ(ready for E_i)  AJfs/fs
    double p = phonon_energy/(constants::hbar_r*3650.0*1e5); // A^-1;  E_D = p * c_s * hbar //  1e15 / A1e10 = 1/A
    
    ea_rate = -sim::ee_coupling_strength*dt*2.0*M_PI/ constants::hbar_r; //(4.16*4.16*4.16*(p*p + 2.5*2.5))/constants::hbar_r;
    // modifier: 0.02
    // ee_rate = -1.0*dt*sim::ee_coupling_strength/(constants::eV_to_AJ*constants::eV_to_AJ); //eV^-2 fs^-1 -> fs**-1 AJ**-2
        // 2pi/hbar e^4 / eps_o^2 (1.25pi 2.75^3)^2 
    ee_rate = -1.0*dt*2.0*M_PI*pow(constants::e*constants::e/(constants::eps_0_A*4.16*4.16*4.16), 2.0)/constants::hbar_r; // AJ/A^4
    E_f_A -= ((E_f_A - 1.5*constants::eV_to_AJ)*i_dos_en_step - floor((E_f_A - 1.5*constants::eV_to_AJ)*i_dos_en_step))*dos_en_step;
    // core_cutoff = E_f_A - DoS_cutoff;
     std::cout << "E_f(AJ): " << E_f*1e20 << ", discretised (AJ): " << E_f_A << std::scientific << ", gamma(J/m**3/K**2): " << e_heat_capacity*1e7 << ", C_l(J/K/m**3): " << a_heat_capacity*1e7 << ", G@300K(J/K/s/m**3): " <<  G*1e22  << \
    ", ea_rate@300K(J/s/K/m**3): " << -1e22*ea_rate*n_f/300.0 <<  ", tau_ep(fs/AJ): " << 1.0/(-1.0*ea_rate/dt) << ", ee_rate(J/s/K/m^3): " << -1e22*ee_rate*n_f/dt/300.0 << ", tau_ee(fs/AJ): " << 1.0/((-1.0/(11.9*11.9))*ee_rate/dt) << \
    ", phonon energy: " << phonon_energy << std::fixed << std::endl;
        
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
       
   
    // co/nduction_electrons = round(0.5*conduction_electrons);
     double K_avg = 0.0;
     dWdE_standard.resize(int(ceil(12.0*constants::eV_to_AJ*i_dos_en_step)), 0.0);
     double int_dos = 0.0;
    for(int e = 1; e < dos_size-1; e++){
        int_dos += dos_en_step*(1.0*dos_standard[e+1]+4.0*dos_standard[e]+1.0*dos_standard[e-1])*n_f/lattice_atoms/6.0;
        // e-/AJ/atom * atom/nm^3 * 1e-3 nm^3/A^3
        //
        K_avg = pow(3.0e-3*M_PI*M_PI*int_dos, 0.333333333333); //1/nm -> 1/A
        // (3 pi^2 int_0^E  dE D(E) )^1/3
        dWdE_standard[e] = K_avg;
    }


      std::ofstream dWdE_standard_output;
      dWdE_standard_output.open(string(directory) + "/dWdE_standard.csv");
      // min = E_f_A -1.5*constants::eV_to_AJ;
     // double total_e = 0.0;
      for(int d = 2; d < dWdE_standard.size()-2; d++) {
        // dos_standard[d] = return_fermi_distribution(d*dos_en_step+87.5561-E_f_A, Te)*
        if(dWdE_standard[d] ==0) continue;
        double d_dW = (-1.0*dWdE_standard[d+2]+8.0*dWdE_standard[d+1]-8.0*dWdE_standard[d-1]+dWdE_standard[d-2])*i_dos_en_step/12.0;
      
        dWdE_standard_output << (d*dos_en_step + E_f_A - DoS_cutoff) << ", " << dWdE_standard[d] << ", " << \
                                return_dWdE(d*dos_en_step + E_f_A -DoS_cutoff) << ", " << \
                                return_vel(d*dos_en_step + E_f_A - DoS_cutoff) << ", " <<\
                                1.0/(constants::hbar_r*d_dW) << ", " << 1.0/(constants::hbar_r*(dWdE_standard[d+1]-dWdE_standard[d-1])*i_dos_en_step/2.0) <<  '\n';
      }
      dWdE_standard_output.close();


  double q = 0.0;
  double q_1 = 0.0;
  double q_2 = 0.0;
  double f_0;
  // double f_2;
  double k;

  for(int h = 2; h < dos_size-2; h++) {
  
    // if( ((h)*dos_en_step+DoS_cutoff) > transport_cutoff) f_0 = double(global_e_dos[h][0])/(dos_standard[h]*dos_en_step);
    // else f_0 = double(global_e_dos[h][0])/double(global_e_dos[h][1]);
    f_0 = return_fermi_distribution(h*dos_en_step - DoS_cutoff, 0.0);
    k = dWdE_standard[h];
    if(k <= 0.0) continue;
    q += dos_en_step*dos_standard[h]*f_0*1.0/(k*k*lattice_atoms*x_unit_size*y_unit_size*z_unit_size);
    k = return_dWdE(h*dos_en_step-DoS_cutoff);
    q_1 += dos_en_step*dos_standard[h]*f_0*1.0/(k*k*lattice_atoms*x_unit_size*y_unit_size*z_unit_size);
            // s^2 /  fs^2
    q_2 += dos_en_step*dos_standard[h]*(return_fermi_distribution((h+1)*dos_en_step - DoS_cutoff, 0.0)-return_fermi_distribution((h-1)*dos_en_step - DoS_cutoff, 0.0))/(2.0*dos_en_step);
  }
	DoS_cutoff += dos_en_step;
  q = sqrt(constants::K*q);
  q_1 = sqrt(constants::K*q_1);
  q_2 = sqrt(constants::K*constants::K*abs(q_2));
  q_offset = 11.9 - q;

  //  q_sq = k_sq();
   if(err::check) std::cout << "q: " << q << ", q_1 " <<  q_1 << ", q_2: " << q_2 << std::endl;
  q_sq = 2.5*2.5;//11.9*11.9;
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
    dos_occ  = dos_en_step * dos_standard[int(floor((E_f_A-DoS_cutoff)*i_dos_en_step))];
    // local_dos_occ = round(dos_en_step*(ee_density/6.0/(E_f_A-core_cutoff)));
     if(err::check)  std::cout << "DoS(E_f) : " << dos_occ  << std::endl;

   
   //G = -1.0*ea_rate*n_f/300.0; //J/s/K/m**3 []
     initialize_cell_omp(); 
  // else std::cout << "CASTLE lattice integration is most efficient at greater than 4 15 Angstrom unit cells wide. Generating standard OpenMP lattice." << std::endl;
    // std::cout << "electron DoS summation: " << int(round(4*ee_density/3.0/(3*constants::kB_r*300.0+E_f_A - core_cutoff)/4.0))  << std::endl;
    output_verbose_initialisation_data(string(directory));
   if (err::check)  std::cout << "Total lattice cells: " << total_cells << ", maximum cell size: " << cell_integration_lists[0].size() << ", maximum integration list size: " << electron_integration_list[0].size() << std::endl;
   
    /*
    // std::string dWdE_name = "Ni-dWdE.csv";
    // std::vector<std::vector< double> > dWdE_scale;
    // dWdE_scale.resize(25);
    // std::fstream dWdE_file (dWdE_name, std::ios::in);
    //  min = E_f_A - DoS_cutoff;
    //  double offset = 0.0;
    //   if(dWdE_file.is_open()) {
    //     for (int l = 0; l < 25; l++) {

    //       dWdE_scale[l].resize(2, 0.0);
    //       std::getline(dWdE_file, line);
    //       std::stringstream s(line);
    //       std::string s_value;
    //       int count = 0;
    //       while(std::getline(s, s_value, ',')) {
    //         if(count == 0) {
    //           if(l == 0) {
    //             offset = std::stod(s_value)*constants::eV_to_AJ;
    //             dWdE_scale[l][count] = min;
    //           }
    //           else dWdE_scale[l][count] = std::stod(s_value)*constants::eV_to_AJ - offset + min;
    //           }
    //         else dWdE_scale[l][count] = std::stod(s_value);
    //           count++;
    //         }
    //     }
    //   } else std::cout << dWdE_name << " did not open" << std::endl;
    //   dWdE_file.close();
    //   dWdE_standard.resize(dos_size, dWdE_scale[24][1]);
    //   dWdE_standard_i.resize(dos_size, 0.0);
    //    d_count = 0;
    //   for(int l = 0; l < 25; l++) {

    //     int count = 0;
    //     double avg = 0.0;

    //     while(dWdE_scale[l][0] > min) {
    //       avg += dWdE_scale[l][1];
    //       min += dos_en_step;
    //       count++;
    //       d_count++;
    //       // avg /= count;
    //       dWdE_standard[d_count-1] = avg / count;
    //       if(d_count == dos_size) {l = 25; break;}
    //       // std::cout << d_count-1 << ", " << avg / count << std::endl;
    //     }
    //     if(count == 0)  {
    //       dWdE_standard[d_count] = dWdE_scale[l][1];
    //       min += dos_en_step;
    //       d_count++;
    //     }
    //   }
 */
   
    
    // double k_ffset = 3.4
    // kg = kg fs / A  1e30 
    // m  = hbar k/v
    //  m_eff_ratio = return_dWdE(E_f_A+5.0)*return_vel(E_f_A)/return_vel(E_f_A+5.0)/return_dWdE(E_f_A);
    //  std::cout << "m_ee ratio: " << m_eff_ratio << std::endl;
   // temp_data.open(string(directory) + "/temp_data.csv");
    mean_data.open(string(directory) + "/mean_data.csv");
    mean_data << "time, step, mean-EKE, mean-LE, mean-Te, mean-Tp,  mean-radius, mean-e-a-collisions, mean-e-e-collisions, mean-x_flux, mean-y_flux, mean-z_flux" << "\n";

}

void initialize_cell_omp() {

  //number of cells in each lattice direction
  x_omp_cells = 6;
  y_omp_cells = 6;
  z_omp_cells = 6;

  total_cells = x_omp_cells*y_omp_cells*z_omp_cells;

  x_step_size = lattice_width / double(x_omp_cells);
  y_step_size = lattice_depth / double(y_omp_cells);
  z_step_size = lattice_height/ double(z_omp_cells);
      if(err::check) std::cout << "step sizes(A): " << x_step_size << ", " << y_step_size << ", " << z_step_size << std::endl;
  boundary_conditions_cutoff = 1.0*fmax(x_step_size, fmax(y_step_size, z_step_size));
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
  cell_lr_neighbor_list.resize(total_cells);
  std::vector<int> spiral_cell_counter;
  spiral_cell_counter.resize(total_cells, 0);
  std::vector<int> lr_cell_counter;
  lr_cell_counter.resize(total_cells, 0);
      if(err::check) std::cout << "electrons in cells; checking occupancy" << std::endl;

  int total_c = 0;
  for(int c = 0; c < total_cells; c++) {
    total_c += cell_integration_lists[c][0]-1;
  }
    if(total_c != conduction_electrons) std::cout << "not enough electrons: " << total_c << " in " << total_cells << " instead of " << conduction_electrons <<std::endl;
  
  //spiral lattice integration
  for(int c = 0; c < total_cells; c++) {
    cell_nearest_neighbor_list.at(c).resize(27); 
    cell_lr_neighbor_list.at(c).resize(125);
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

    for(int s = 0; s < 125; s++) {

      int cell_spiral_xmod = x_cell + (s % 5) - 2;
      if(cell_spiral_xmod < -1) cell_spiral_xmod = x_omp_cells - 2;
      else if(cell_spiral_xmod < 0) cell_spiral_xmod = x_omp_cells - 1;
      else if (cell_spiral_xmod > x_omp_cells) cell_spiral_xmod = 1;
      else if (cell_spiral_xmod > x_omp_cells-1) cell_spiral_xmod = 0;
        // if(s == 0) std::cout << "<" << cell_spiral_xmod;

      int cell_spiral_ymod = y_cell + (int(floor(s / 5)) % 5) - 2;
      if(cell_spiral_ymod < -1) cell_spiral_ymod = y_omp_cells - 2;
      else if(cell_spiral_ymod < 0) cell_spiral_ymod = y_omp_cells - 1;
      else if (cell_spiral_ymod > y_omp_cells) cell_spiral_ymod = 1;
      else if (cell_spiral_ymod > y_omp_cells-1) cell_spiral_ymod = 0;
         // if(s == 0) std::cout << "," << cell_spiral_ymod;

      int cell_spiral_zmod = z_cell + int(floor(s / 25)) - 2;
      if(cell_spiral_zmod < -1) cell_spiral_zmod = z_omp_cells - 2;
      else if(cell_spiral_zmod < 0) cell_spiral_zmod = z_omp_cells - 1;
      else if (cell_spiral_zmod > z_omp_cells) cell_spiral_zmod = 1;
      else if (cell_spiral_zmod > z_omp_cells-1) cell_spiral_zmod = 0;
         // if(s == 0) std::cout << "," << cell_spiral_zmod << ">" << std::endl;

      int cell = lattice_cell_coordinate.at(cell_spiral_xmod).at(cell_spiral_ymod).at(cell_spiral_zmod);
      cell_lr_neighbor_list.at(c).at(s) = cell;
      lr_cell_counter[cell]++;
      //if(s == 13) std::cout << cell_nearest_neighbor_list.at(c).at(s) << ", " << c << ", " << x_cell << ", " << y_cell << ", " << z_cell << std::endl;
   
    }
  }

  for(int c = 0; c < total_cells; c++) {
    if(spiral_cell_counter[c] != 27) std::cout << "spiral cell error on cell " << c << " with " << spiral_cell_counter[c] << std::endl;
    if(lr_cell_counter[c] != 125) std::cout << "lr cell error on cell " << c << " with " << lr_cell_counter[c] << std::endl;
  }
  
  
    if(err::check) std::cout << "spiral integration coordiantes initialized." << std::endl;

    //number of cells each thread takes in each lattice direction
    // (12/6)^3 = n_threads; 24 -> 2*2*6
    const int max_x_threads = 2; //3 
    const int max_y_threads = 2; //3 
    const int max_z_threads = 3; //2 -> 18

    int max_total_threads = (x_omp_cells/max_x_threads) *(y_omp_cells/ max_y_threads) * (z_omp_cells/ max_z_threads);
   if(max_total_threads != omp_threads) std::cout << "maximum omp threads based on given lattice parameters: " << max_total_threads << "\n Given threads: " << omp_threads << "\n Reducing to max threads" << std::endl;
   omp_threads = int(std::min(double(max_total_threads), double(omp_threads)));

    cells_per_thread = total_cells / omp_threads;
    if (err::check) std::cout << "thread count " << omp_threads << ", " << cells_per_thread << ", " << total_cells << ", " << x_omp_cells << ", " << y_omp_cells << ", " << z_omp_cells << std::endl;
    lattice_cells_per_omp.resize(omp_threads);

    for(int t=0; t< omp_threads;t++){
      lattice_cells_per_omp.at(t).resize(cells_per_thread);
    }

    #pragma omp parallel
    {
    const int thread = omp_get_thread_num();
    for(int l = 0; l < cells_per_thread; l++) {
      if(err::check && l==0) {
        #pragma omp critical 
        std::cout << "thread: " << thread << ", x_offset: " << (thread*max_x_threads)%x_omp_cells << ", y offset: " << (max_y_threads)*(int(floor(thread*max_x_threads/(x_omp_cells)))%(y_omp_cells/ max_y_threads)) << ", z offset: " << (max_z_threads)*floor(thread/(x_omp_cells*y_omp_cells/(max_x_threads*max_y_threads))) << std::endl;
      }
      // #pragma omp critical
      // std::cout << omp_get_thread_num() << ", " << (omp_get_thread_num()*2)%x_omp_cells + floor(l/z_omp_cells) << ", " << (int(floor(omp_get_thread_num()*2/x_omp_cells)*2)%(y_omp_cells)) + int(floor(l/(z_omp_cells/2)))%2 << ", " << 4*floor(omp_get_thread_num()/(2*y_omp_cells)) + (l % (z_omp_cells/2)) << std::endl;
      
      //double decker cells
    
     lattice_cells_per_omp[thread].at(l) = lattice_cell_coordinate.at((thread*max_x_threads)%x_omp_cells + int(floor(l/(max_y_threads*max_z_threads))))\
                                                                  .at((max_y_threads)*(int(floor(thread*max_x_threads/(x_omp_cells)))%(y_omp_cells/ max_y_threads)) + int(floor(l/(max_z_threads)))%(max_y_threads))\
                                                                  .at((max_z_threads)*floor(thread/(x_omp_cells*y_omp_cells/(max_x_threads*max_y_threads))) + (l % (max_z_threads)));
      
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
   
} 


//====================================
// Creates and outputs atomic lattice
//      Currently static lattice
//====================================
void initialize_positions() {

   // initialize_lattice();

    // if (err::check)  std::cout << "Lattice" << std::endl;

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

            if (err::check) std::cout << "Lattice output file and electron position file opened..." << std::endl;
    //=========
    // Initialize arrays for holding electron variables
    //      Arrays in super array format to take advantage of caching
    //========
      
    e_a_scattering_count = 0;
    e_e_scattering_count = 0;
    ee_transport_scattering_count = 0;
    ee_core_scattering_count = 0;
    ea_transport_scattering_count = 0;
    ea_core_scattering_count = 0;
    // ee_scattering_angle = sim::ee_scattering_angle;
    // e_e_neighbor_cutoff = pow((lattice_width/4.0)-1.0,2.0);
    half_int_var.resize(2,0);
    half_int_var[0] =  3;
    half_int_var[1] =  3;
    
    e_e_integration_cutoff = lattice_width/6.0; //
    e_e_coulomb_cutoff = 4.16*2.0; //
    double deltaX = e_e_integration_cutoff-e_e_coulomb_cutoff;
    std::cout << "band 1 velocity(A/fs): " << return_vel(E_f_A) << ", minimun separation criteria(dt): " << floor(deltaX/(2.0*return_vel(E_f_A)*dt)) << "...";
    if( (4*return_vel(E_f_A)*dt*half_int_var[0]) < deltaX) {
      terminaltextcolor(GREEN);
       std::cout << "criteria exceeded. Consider increasing stride to " << floor(deltaX/(2*return_vel(E_f_A)*dt)) << std::endl;
      terminaltextcolor(WHITE);
    }
    else if ( (2*return_vel(E_f_A)*dt*half_int_var[0]) < deltaX)  std::cout << "criteria met: " << deltaX/(2*return_vel(E_f_A)*dt*half_int_var[0]) << std::endl;
    else {
      terminaltextcolor(RED);
      std::cerr << "criteria not met: " << deltaX << " < " << 2*return_vel(E_f_A)*dt*half_int_var[0] << std::endl;
    } 
   //  deltaX = 2*e_e_integration_cutoff - e_e_coulomb_cutoff;
   //  std::cout << "band 2 velocity(A/fs): " << return_vel(E_f_A+5.0) << ", minimun separation criteria(dt): " << floor(deltaX/(2*return_vel(E_f_A+5.0)*dt)) << "...";
   //  if( (4*return_vel(E_f_A+5.0)*dt*half_int_var[1]) < deltaX) {
   //    terminaltextcolor(GREEN);
   //     std::cout << "criteria exceeded. Consider increasing stride to " << floor(deltaX/(2*return_vel(E_f_A+5.0)*dt)) << std::endl;
   //    terminaltextcolor(WHITE);
   //  }
   //  else if ( (2*return_vel(E_f_A+5.0)*dt*half_int_var[1]) < deltaX)  std::cout << "criteria met" << std::endl;
   //  else {
   //    terminaltextcolor(RED);
   //    std::cerr << "criteria not met: " << deltaX << " < " << 2*return_vel(E_f_A+5.0)*dt*half_int_var[1] << std::endl;
   //  }

    e_e_integration_cutoff = pow(e_e_integration_cutoff,2.0);
    e_e_coulomb_cutoff = pow(e_e_coulomb_cutoff, 2.0);
    
      if (err::check) std::cout << "Prepare to set arrays: " << std::endl;
          
    electron_position.resize(conduction_electrons * 3, 0); // ""'Memory is cheap. Time is expensive' -Steve Jobs; probably" -Michael Scott." -Headcannon.
    electron_velocity.resize(conduction_electrons * 3, 0); //Angstroms
    electron_potential.resize(conduction_electrons, 0);
            
    relaxation_time_hist_ee.resize(3*conduction_electrons);

    electron_integration_list.resize(conduction_electrons);
   
    electron_ee_scattering_list.resize(conduction_electrons);
    electron_ea_scattering_list.resize(conduction_electrons);
    flux_index.resize(dos_size,0.0); 
    global_tau_ep.resize(dos_size*2, 0.0);
    global_tau_ee.resize(dos_size, 0.0);

        if (err::check) std::cout << "Prepare to set position: " << std::endl;
     ee_density =   3*int(round(pow(e_e_integration_cutoff,1.5)*1.25*M_PI * total_e_scaling*n_f * 1e-3));
    //  ee_density =  3*int(round(pow(e_e_neighbor_cutoff, 1.5)*1.25*M_PI * 3.8*n_f * 1e-3));
    const int ee_scattering = 12*round(pow(e_e_coulomb_cutoff,   1.5)*1.25*M_PI * total_e_scaling*n_f * 1e-3);
        if (err::check)  std::cout << ee_density << ", " << ee_scattering << std::endl;
    
    omp_set_dynamic(0);
    omp_set_num_threads(omp_threads);
    #pragma omp parallel for schedule(dynamic,4) 
    for (int e = 0; e < conduction_electrons; e++) {

        electron_integration_list[e].resize(ee_density, 0);
        // electron_nearest_electron_list.at(e).resize(ee_density, 0);
        electron_ee_scattering_list[e].resize(ee_scattering, 0);
        electron_ea_scattering_list[e].resize(2,0);
        
         const int array_index = 3*e;
        relaxation_time_hist_ee[array_index].resize(4*60,0);
        relaxation_time_hist_ee[array_index+1].resize(4*60,0);
        relaxation_time_hist_ee[array_index+2].resize(4*60,0);

        // ee_dos_hist.at(e).resize(dos_size,0);
        // relaxation_time_hist_ea[3*e].resize(4*70,0);
        // relaxation_time_hist_ea[3*e+1].resize(4*70,0);
        // relaxation_time_hist_ea[3*e+2].resize(4*70,0);
        
       
        electron_position[array_index]     = atoms::x_coord_array.at((2*e)%int(lattice_atoms)) + 0.5*x_unit_size;
        electron_position[array_index + 1] = atoms::y_coord_array.at((2*e)%int(lattice_atoms)) + 0.5*y_unit_size;
        electron_position[array_index + 2] = atoms::z_coord_array.at((2*e)%int(lattice_atoms)) + 0.5*z_unit_size;
        //initialize and output electron posititons
      //  = atom_anchor_position.at(3*(e%lattice_atoms));//   + cos(theta)*sin(phi)*screening_depth;//*radius_mod(gen)); //Angstroms
       // electron_position.at(array_index + 2) = atom_anchor_position.at(3*(e%lattice_atoms)+2);// + cos(phi)*screening_depth;//*radius_mod(gen);
    }
}

void initialize_forces() {

    initialize_electron_interactions();
    if (err::check)  std::cout << "Electron interactions" << std::endl;

   // initialize_atomic_interactions();
    if (err::check)  std::cout << "Atomic interactions" << std::endl;

   // initialize_electron_atom_interactions();
    if (err::check)  std::cout << "Atomic electron interactions" << std::endl;
}

void initialize_electron_interactions() {
  // omp_set_dynamic(0);
  //      omp_set_num_threads(omp_threads);
  //   #pragma omp parallel for schedule(static)
  //   for (int e = 0; e < conduction_electrons; e++) {
  //     int array_index_i;
     
  //     double x_distance,y_distance,z_distance, length;
  //     int array_index = 3*e;
  //     int integration_count = 1;
  //     int neighbor_count = 1;
  //     int scattering_count = 2;
  //               if (err::check) if(e ==0) std::cout << "Calculating conduction electron repulsion" << std::endl;

  //       for (int i = 0; i < conduction_electrons; i++) {
  //           if (i == e) continue; //no self repulsion

  //           array_index_i = 3*i;
  //           x_distance = electron_position.at(array_index)     - electron_position.at(array_index_i);
  //           y_distance = electron_position.at(array_index + 1) - electron_position.at(array_index_i + 1);
  //           z_distance = electron_position.at(array_index + 2) - electron_position.at(array_index_i + 2);
            
  //           if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance = x_distance + lattice_width;
  //           else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance = x_distance - lattice_width;

  //           if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance = y_distance + lattice_depth;
  //           else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance = y_distance - lattice_depth;

  //           if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance = z_distance + lattice_height;
  //           else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance = z_distance - lattice_height;

  //           length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);

  //           if(length > e_e_integration_cutoff) continue;
  //           electron_integration_list.at(e).at(integration_count) = array_index_i;
  //           integration_count++;

  //           // if (length > e_e_neighbor_cutoff) continue;
  //           // electron_nearest_electron_list.at(e).at(neighbor_count) = array_index_i/3;
  //           // neighbor_count++;

  //           if(length > e_e_coulomb_cutoff) {
              
  //             continue; }
  //             //std::cout << e << ", " << i << ", " << length << std::endl;
  //           electron_ee_scattering_list.at(e).at(scattering_count) = array_index_i/3;
  //           scattering_count++;
  //       }
  //       electron_integration_list.at(e)[0] = integration_count;
  //       // electron_nearest_electron_list.at(e)[0] = neighbor_count;
  //       electron_ee_scattering_list.at(e)[1] = scattering_count;
  //    }
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
    

    const std::string n = "Init_E_distrib";
    // std::cout << "conduction electrons " << conduction_electrons << std::endl;

   create_defined_fermi_distribution(n, electron_potential, Te);
        if (err::check) std::cout << "distribution generated" << std::endl;
  
      omp_set_dynamic(0);
       omp_set_num_threads(omp_threads);
       p_x = 0.0;
       p_y = 0.0;
       p_z = 0.0;

    // double minimum = E_f_A;
   // transport_cutoff = 106.978;
    std::random_device rd;
    std::mt19937 g(rd());
    for(int r = 0; r < 7; r++) {
      std::shuffle(electron_potential.begin(), electron_potential.end(), g);
    }
    // double min = 40.0;
    #pragma omp parallel 
    {
      int thread = omp_get_thread_num();
      #pragma omp for schedule(dynamic, 4) reduction(+:p_x,p_y,p_z) 
      for(int e = 0; e < conduction_electrons; e++) {
        double phi;
        double theta;
        const int array_index = 3*e;
   
        const double energy = electron_potential[e];
        
        double mom = return_dWdE(energy);
        // if(vel < min) min = vel;

        if(sim::CASTLE_x_vector < 0.0 && sim::CASTLE_y_vector < 0.0 && sim::CASTLE_z_vector < 0.0) {
          theta = 2.0*M_PI*omp_uniform_random[thread]();
          phi = M_PI*omp_uniform_random[thread]();
         
        } else {
          const double unit = sqrt((sim::CASTLE_x_vector*sim::CASTLE_x_vector)+(sim::CASTLE_y_vector*sim::CASTLE_y_vector)+(sim::CASTLE_z_vector*sim::CASTLE_z_vector));
          theta = atan2(sim::CASTLE_y_vector, sim::CASTLE_x_vector);
            if(theta != theta) theta = 0.0;
          phi = acos(sim::CASTLE_z_vector / unit);
         // if(sim::CASTLE_z_vector < 0.0) theta += M_PI;
        }

          if(theta != theta || phi != phi || mom != mom) std::cout << theta << ", " << phi << ", " << mom << ", " << energy << std::endl;
      
            if (err::check) if(e==0) std::cout << "Electron velocity ready..." << std::endl;
        electron_velocity[array_index]     = cos(theta)*sin(phi); 
        electron_velocity[array_index + 1] = sin(theta)*sin(phi);
        electron_velocity[array_index + 2] = cos(phi); 
        
        p_x += electron_velocity[array_index]*mom;
        p_y += electron_velocity[array_index+1]*mom;
        p_z += electron_velocity[array_index+2]*mom;
      }
    }

    p_x /= double(conduction_electrons);
    p_y /= double(conduction_electrons);
    p_z /= double(conduction_electrons);

  //  std::cout << min << std::endl;
        if (err::check) std::cout << "distribution assigned" << std::endl;
      std::ofstream Init_E_vel;
      Init_E_vel.open("Init_E_vel");

  double n_p_x = 0.0;
  double n_p_y = 0.0;
  double n_p_z = 0.0;
  #pragma omp parallel 
  {
    std::vector<int> local_e_dos;
    local_e_dos.resize(dos_size,0);
    #pragma omp for schedule(dynamic, 4) nowait reduction(+:n_p_x,n_p_y,n_p_z)
    for(int e = 0; e < conduction_electrons; e++) {
      int array_index = 3*e;
      double energy = electron_potential[e];
      double mom = return_dWdE(energy);
      local_e_dos[int(std::min(dos_size-1.0, std::max(0.0, floor((energy - DoS_cutoff)*i_dos_en_step))))]++;
     
      double v_x = electron_velocity[array_index] - p_x/mom;
      double v_y = electron_velocity[array_index+1] - p_y/mom;
      double v_z = electron_velocity[array_index+2] - p_z/mom;

      double v_mod = 1/sqrt(v_x*v_x + v_y*v_y + v_z*v_z);

      electron_velocity[array_index]   = v_x*v_mod;
      electron_velocity[array_index+1] = v_y*v_mod;
      electron_velocity[array_index+2] = v_z*v_mod;

      n_p_x += v_x*v_mod;
      n_p_y += v_y*v_mod;
      n_p_z += v_z*v_mod;
    }
    
    #pragma omp critical 
    {
      for(int h = 0; h < dos_size; h++) {
        global_e_dos[h][0] += local_e_dos[h];
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
      // if(h == 0) std::cout << "DoS vs core offset " << h*dos_en_step + core_cutoff - DoS_cutoff << ", Discrete: " << int(floor(h*dos_en_step + core_cutoff - DoS_cutoff)*i_phonn_energy) << std::endl;
      global_e_dos_out << h*dos_en_step + DoS_cutoff << ", " <<  global_e_dos[h][1] << ", " << dos_standard[int(floor((h*dos_en_step) *i_dos_en_step))]*dos_en_step << std::endl;
      total += global_e_dos[h][1];
    }

    if(total != conduction_electrons) std::cout <<"hist problem " << total << ", " << conduction_electrons << std::endl;

    std::cout << "core cutoff: " << core_cutoff << ", transport cutoff: " << transport_cutoff << ", ballistic cutoff: " << E_f_A+4.86166 << std::endl;
    std::cout << "center of mass adjustment: <" << p_x << "," << p_y << "," << p_z << "> -> <" << n_p_x/double(conduction_electrons) << "," << n_p_y/double(conduction_electrons) << "," << n_p_z/double(conduction_electrons) << ">" <<  std::endl;
   

       if (err::check)  std::cout << "distribution output" << std::endl;
  
            if (err::check) std::cout << "Electron velocity ready..." << std::endl;
}

void output_verbose_initialisation_data(string directory) {

    //mu_F
    std::ofstream d_U_output_Ni;
    d_U_output_Ni.open(string(directory) + "/dU_Ni_calibration.txt");
    d_U_output_Ni.precision(10);
    std::ofstream d_U_output_Ag;
    d_U_output_Ag.open(string(directory) + "/dU_Ag_calibration.txt");
    d_U_output_Ag.precision(10);
    std::ofstream d_U_output_feg;
    d_U_output_feg.open(string(directory) + "/dU_feg_calibration.txt");
    d_U_output_feg.precision(10);

    double U_int_min_Ni = 0.0;
    double U_min_Ni = 0.0;
    double U_int_zero_Ni = 0.0;

    double U_int_min_Ag = 0.0;
    double U_min_Ag = 0.0;
    double U_int_zero_Ag = 0.0;

    double U_int_min_feg = 0.0;
    double U_min_feg = 0.0;
    double U_int_zero_feg = 0.0;

    double u_Ni  = 0.0;
    double u_Ag = 0.0;
    double u_feg = 0.0;
    double dos_mu;//
    double dos_e_f = dos_standard[int(floor((E_f_A-DoS_cutoff)*i_dos_en_step))];
    double delta_dos_Ni = 0.0;
    double delta_dos_feg = 0.0;// 0.75*conduction_electrons/E_f_A;
    double Ag_dos = 0.33333 * lattice_atoms;
    double feg_dos = 1.5*conduction_electrons/E_f_A; // *sqrt(e/E_f)

    for( double t = 0.0; t < 2000.0; t++) {
      u_Ni = (t < 400.0) ? 2.42*t/1250.0 : 2.42*(1.0+ t-400.0)/225.0;
      u_feg =  -1.0*pow(M_PI*constants::kB_r*t/E_f_A, 2.0)/12.0;
      double dU_Ag = 0.0;
      double U_int_Ag = 0.0;
      double U_int_temp_Ag = 0.0;

      double dU_feg = 0.0;
      double U_int_feg = 0.0;
      double U_int_temp_feg = 0.0;

      double dU_Ni = 0.0;
      double U_int_Ni = 0.0;
      double U_int_temp_Ni = 0.0;

      double transient_entropy = 0.0;
      for( double e = core_cutoff; e < E_f_A+10.0*constants::kB_r*t; e += dos_en_step) {
        int index  = int(floor((e-DoS_cutoff)*i_dos_en_step));
        //Ni
        if(e < E_f_A) dU_Ni += dos_standard[index]*dos_en_step*(E_f_A - e)*abs(1.0-(return_fermi_distribution(e-E_f_A - u_Ni, t)+4.0*return_fermi_distribution(e + 0.5*dos_en_step-E_f_A-u_Ni,t)+return_fermi_distribution(e + dos_en_step-E_f_A-u_Ni,t))/6.0);
        else dU_Ni += dos_standard[index]*dos_en_step*(e - E_f_A)*(return_fermi_distribution(e -E_f_A - u_Ni,t)+4.0*return_fermi_distribution(e + 0.5*dos_en_step-E_f_A- u_Ni,t)+return_fermi_distribution(e + dos_en_step-E_f_A-u_Ni,t))/6.0;
        U_int_Ni += dos_standard[index]*e*dos_en_step*(return_fermi_distribution(e -E_f_A-u_Ni,t)+4.0*return_fermi_distribution(e +0.5*dos_en_step-E_f_A-u_Ni,t)+return_fermi_distribution(e +dos_en_step-E_f_A-u_Ni,t))/6.0;

        //Ag
        if(e < E_f_A) dU_Ag += Ag_dos*dos_en_step*(E_f_A-e)*abs(1.0-(return_fermi_distribution(e-E_f_A + u_Ag, t)+4.0*return_fermi_distribution(e + 0.5*dos_en_step-E_f_A+u_Ag,t)+return_fermi_distribution(e + dos_en_step-E_f_A+u_Ag,t))/6.0);
        else dU_Ag += Ag_dos*dos_en_step*(e-E_f_A)*(return_fermi_distribution(e-E_f_A + u_Ag, t)+4.0*return_fermi_distribution(e + 0.5*dos_en_step-E_f_A+u_Ag,t)+return_fermi_distribution(e + dos_en_step-E_f_A+u_Ag,t))/6.0;
        U_int_Ag += Ag_dos*e*dos_en_step*(return_fermi_distribution(e -E_f_A+u_Ag,t)+4.0*return_fermi_distribution(e +0.5*dos_en_step-E_f_A+u_Ag,t)+return_fermi_distribution(e +dos_en_step-E_f_A+u_Ag,t))/6.0;
      
        //feg
        if(e < E_f_A) dU_feg += feg_dos*sqrt(e/E_f_A)*dos_en_step*(E_f_A-e)*abs(1.0-(return_fermi_distribution(e-E_f_A -u_feg, t)+4.0*return_fermi_distribution(e + 0.5*dos_en_step-E_f_A-u_feg,t)+return_fermi_distribution(e + dos_en_step-E_f_A-u_feg,t))/6.0);
        else dU_feg += feg_dos*sqrt(e/E_f_A)*dos_en_step*(e-E_f_A)*(return_fermi_distribution(e-E_f_A -u_feg, t)+4.0*return_fermi_distribution(e + 0.5*dos_en_step-E_f_A-u_feg,t)+return_fermi_distribution(e + dos_en_step-E_f_A-u_feg,t))/6.0;
        U_int_feg += lattice_atoms*feg_dos*sqrt(e/E_f_A)*e*dos_en_step*(return_fermi_distribution(e-E_f_A -u_feg, t)+4.0*return_fermi_distribution(e + 0.5*dos_en_step-E_f_A-u_feg,t)+return_fermi_distribution(e + dos_en_step-E_f_A-u_feg,t))/6.0;

        double dos = dos_standard[index]/lattice_atoms;
        double fermi_dist = return_fermi_distribution(e -E_f_A - u_Ni,t);
        double onemin_fermi_dist = 1 - fermi_dist;
     
        transient_entropy += dos_en_step*dos*(fermi_dist*log(fermi_dist) + ((onemin_fermi_dist == 0) ? 0.0 : onemin_fermi_dist*log(onemin_fermi_dist)));
          
      }
      //Ni
      int mu_index  = int(floor((E_f_A + u_Ni -DoS_cutoff)*i_dos_en_step));
      delta_dos_Ni = 0.0;//0.5*(dos_standard[mu_index+1]-dos_standard[mu_index-1])/dos_en_step; 
      dos_mu = dos_standard[int(floor((E_f_A + u_Ni-DoS_cutoff)*i_dos_en_step))];
      dU_Ni /= lattice_atoms;
      dU_Ni /= abs(E_f_A*delta_dos_Ni + dos_e_f)/lattice_atoms;
      U_int_Ni /=  conduction_electrons;//*(dos_e_f/lattice_atoms);
      if(t == 0.0) U_int_zero_Ni = U_int_Ni;
      U_int_min_Ni = U_int_zero_Ni - u_Ni*E_f_A*dos_e_f/lattice_atoms;

      U_int_temp_Ni = sqrt(6.0*std::max(0.0, U_int_Ni-U_int_min_Ni)/(M_PI*M_PI*constants::kB_r*constants::kB_r*abs(E_f_A*delta_dos_Ni + dos_e_f)/lattice_atoms));
      d_U_output_Ni << t << ", " << dU_Ni <<  ", " << sqrt(std::max(0.0,6.0*dU_Ni))/M_PI/constants::kB_r << ", "  << U_int_Ni << ", " << U_int_temp_Ni <<  ", " << transient_entropy << std::endl;
      
      //Ag
      dU_Ag /= lattice_atoms;
      dU_Ag /= Ag_dos/lattice_atoms;
      U_int_Ag /= conduction_electrons;
      if(t == 0.0) U_int_zero_Ag = U_int_Ag;
      U_int_min_Ag = U_int_zero_Ag - u_Ag*E_f_A*Ag_dos/lattice_atoms;
      
      U_int_temp_Ag = sqrt(6.0*std::max(0.0, U_int_Ag-U_int_min_Ag)/(M_PI*M_PI*constants::kB_r*constants::kB_r*abs(Ag_dos)/lattice_atoms));
      d_U_output_Ag << t << ", " << dU_Ag <<  ", " << sqrt(std::max(0.0,6.0*dU_Ag))/M_PI/constants::kB_r << ", "  << U_int_Ag << ", " << U_int_temp_Ag << std::endl;
      
      //feg
      dU_feg /= abs(E_f_A*delta_dos_feg/sqrt(E_f_A*(E_f_A-u_feg)) + feg_dos);
      U_int_feg /= conduction_electrons;
      if(t == 0.0) U_int_zero_feg = U_int_feg;
      U_int_min_feg = U_int_zero_feg - u_feg*E_f_A*feg_dos;
      U_int_temp_feg = sqrt(6.0*std::max(0.0, U_int_feg-U_int_min_feg)/(M_PI*M_PI*constants::kB_r*constants::kB_r*abs(E_f_A*delta_dos_feg/sqrt(E_f_A*(E_f_A-u_feg)) + feg_dos)));
      d_U_output_feg << t << ", " << dU_feg <<  ", " << sqrt(std::max(0.0,6.0*dU_feg))/M_PI/constants::kB_r << ", "  << U_int_feg << ", " << U_int_temp_feg << ", " << U_int_min_feg << std::endl;
     
    }
    d_U_output_Ni.close();
    d_U_output_Ag.close();
    d_U_output_feg.close();
            if(err::check) std::cout << "normalised Uint for 0K: " << U_int_zero_Ni << std::endl;


    //Au e-p BTE
    const double e_energy = E_f_A;
    const double e_mom = return_dWdE(e_energy); // 1/A
    const double c_s = 3560.0; //m/s
    const double E_D = 32.3e-3*constants::eV_to_AJ; //eV -> AJ
    const double T_D = 225.0; //K
    const double V = 4.16*4.16*4.16; //A^3
    const double k = 2.5; //1/A
    const double dos_p = 9.0/(pow(constants::kB_r*T_D, 3.0)*V); // 9/V E_D

    const double M_ep_coeff = 1e30*constants::e*constants::e/(2.0*constants::eps_0*V); // 1e30 kg/ s^2 / A^2
    const double coeff = V*M_PI*M_PI*M_PI/(constants::hbar_r*e_mom);
    
    std::ofstream BTE_out;
    double integrand_p = 0.0;
    double integrand_m = 0.0;
    double temp = 300.0;
    BTE_out.open("BTE-output.txt");
    for (double e = dos_en_step; e < E_D; e += dos_en_step) {
        double d_e_p = e_energy + e;
        double d_e_m = e_energy - e;
        double d_e_mom_p = return_dWdE(d_e_p);
        double d_e_mom_m = return_dWdE(d_e_m);
        double q = e*1e5/(constants::hbar_r*c_s); //  1e5/ A
        double M_ep = M_ep_coeff/(q*q + k*k);
        double f_p = return_fermi_distribution(d_e_p + DoS_cutoff - E_f_A, temp);
        double f_m = return_fermi_distribution(d_e_m + DoS_cutoff - E_f_A, temp);
        double f = return_fermi_distribution(e + DoS_cutoff - E_f_A, temp);
        double g = return_BE_distribution(e, temp);
        double F_p = f_p*(1.0-f)*(g+1.0) - f*(1.0-f_p)*g;
        double F_m = f_m*(1.0-f)*g - f*(1.0-f_m)*(g+1.0);
        double chi_p = (q*q + e_mom*e_mom - d_e_mom_p*d_e_mom_p)/(2.0*q*e_mom);
        double chi_m = (q*q + e_mom*e_mom - d_e_mom_m*d_e_mom_m)/(2.0*q*e_mom);
        double scat_p = coeff * dos_en_step * (dos_standard[int(floor((d_e_p-DoS_cutoff)*i_dos_en_step))]/(e_mom+q))*(dos_p*e*e/q)*F_p;
        double scat_m = coeff * dos_en_step * (dos_standard[int(floor((d_e_m-DoS_cutoff)*i_dos_en_step))]/(e_mom-q))*(dos_p*e*e/q)*F_m;
        integrand_p += scat_p;
        integrand_m += scat_m;
        std::cout << q << " + " << e_mom << " = " << q+e_mom << "; " << d_e_mom_p << std::endl;// << ", " << chi_p << ", " << scat_m << ", " << chi_m << std::endl;
    } 
    std::cout << " 0, 0, 0, 0, 0, " << integrand_p << ", " << integrand_m << std::endl;
}
void create_fermi_distribution(const std::string& name, std::vector<double>& distribution, const double temp) {

  char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog. fermi dist" << std::endl;
      }
  std::ofstream distrib;
  distrib.open(string(directory) +"/"+ name);
  distrib.precision(20);

  double step_size = dos_en_step/100.0;
  const double max  = E_f_A;// + 3.0*constants::kB_r*Te;
  const double min = E_f_A - 1.5*constants::eV_to_AJ;
  std::cout << "DoS E_f_A discrete offset " << int(DoS_cutoff) << ". Active space size " << dos_size << std::endl;
 
  global_e_dos.resize(dos_size);
  for(int h = 0; h < dos_size; h++) {
      global_e_dos[h].resize(2,0);
  }
  // for(int e = 0; e < conduction_electrons; e++) {
  //   
  // }

  int count = 0;
  int subCount = 0;
  int energy_step = 0;
  double epsilon = E_f_A;
  if(temp <= 1) {
    transport_cutoff = E_f_A;
    while(count < conduction_electrons && epsilon > min) {
      epsilon = max - step_size*energy_step;
      if(epsilon < min) epsilon = min;
      DoS_cutoff = int(floor((epsilon-min)*i_dos_en_step));

      // if(DoS_cutoff > 123) std::cout << DoS_cutoff << std::endl;
      int occupation = int(floor(dos_standard[DoS_cutoff]*step_size));

      for(int o = 0; o < occupation; o++) {
        distribution.at(count) = epsilon;
        distrib << count << ", " << distribution.at(count) << "\n";
        count++;

        if(count == conduction_electrons) { 
          std::cout << "potential break early at: " << epsilon << " instead of: " << min << std::endl;
           break;}
      }
      energy_step++;
    }
    if(count != conduction_electrons) std::cout << "potential early break: " << count << " rather than " << conduction_electrons << std::endl;
    core_cutoff = epsilon;
    // E_f_A = epsilon;
  } else {
    while(count < conduction_electrons && epsilon > min) { 

      epsilon = max + 30.0*temp*constants::kB_r - step_size*energy_step;
      DoS_cutoff = int(floor((epsilon-min)*i_dos_en_step));
      int occupation = int(round(dos_standard[DoS_cutoff]*step_size*(return_fermi_distribution(epsilon-E_f_A, temp)+return_fermi_distribution(epsilon-step_size-E_f_A, temp)+4.0*return_fermi_distribution(epsilon-0.5*step_size-E_f_A, temp))/6.0));
      
     // int steps = round(1.0 / double(occupation));
      for(int o = 0; o < occupation; o++) {
        distribution.at(count) = epsilon;// + 0.5*step_size;
        distrib << count << ", " << distribution.at(count) << "\n";
        count++;
      
        if(count == conduction_electrons) { 
          std::cout << "potential break early at: " << epsilon << " instead of: " << min << std::endl;
           break;}
      }
      if(1.0 - return_fermi_distribution(epsilon-E_f_A, temp) > 1e-3) transport_cutoff = epsilon;
      energy_step++;
    }
    if(count != conduction_electrons) std::cout << "potential early break: " << count << " rather than " << conduction_electrons << std::endl;
    core_cutoff = epsilon;
  }
  DoS_cutoff =  min;
  // transport_cutoff = E_f_A;
  std::cout << "Discrete DoS_cutoff at " << DoS_cutoff << ", offset " << floor((core_cutoff-DoS_cutoff)*i_dos_en_step) << std::endl;
  distrib.close();
}

void create_defined_fermi_distribution(const std::string& name, std::vector<double>& distribution, const double temp) {

  char directory [256];
      if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog. fermi dist" << std::endl;
      }
  std::ofstream distrib;
  distrib.open(string(directory) +"/"+ name);
  distrib.precision(20);
  int steps = 100;
  double step_size = dos_en_step/double(steps);
  const double max  = E_f_A;// + 3.0*constants::kB_r*Te;
  const double min = E_f_A - DoS_cutoff;
    // std::cout << "DoS E_f_A discrete cutoff " << DoS_cutoff << ". offset " << core_cutoff << std::endl;
  
  global_e_dos.resize(dos_size);
  for(int h = 0; h < dos_size; h++) {
      global_e_dos[h].resize(2,0);
  }

  int count = 0;
  int energy_step = 0;
  double epsilon = E_f_A - step_size*0.001;
  int hist_sub_E_f = floor(( DoS_cutoff)*i_dos_en_step);
  // std::cout << hist_sub_E_f << std::endl;
  int hist_super_E_f = dos_size - hist_sub_E_f;
  transport_cutoff = E_f_A;
  bool full = false;
  double zero_temp = 84.319184;
  int finite_electron_number = 0;
  double mu = 0.0;
  double mu_Ni = 0.0;
  // epsilon = max;
  if(temp <= 1) {
    for(int h = 0; h < hist_sub_E_f && !full; h++) {
      // epsilon = max - step_size*energy_step;
      const double DoS_occupation = dos_en_step* dos_standard[int(floor((epsilon-min)*i_dos_en_step))];
      const int step_occupation = int(floor(DoS_occupation/double(steps)));
        // if(step_occupation == 0) std::cout << DoS_occupation << ", " << step_occupation << std::endl;
      int small_adjust = int(round(DoS_occupation - step_occupation*steps)) % steps;
        // std::cout << small_adjust << ", " << DoS_occupation << ", " << step_occupation << std::endl;
      int offset = 0;
      for(int s = 0; s < steps && !full; s++) {
        if(small_adjust > 0) offset = 1;
        else offset = 0;
        small_adjust--;
        for(int o = 0; o < (step_occupation + offset); o++) {
          distribution.at(count) = epsilon;
          distrib << count << ", " << epsilon << "\n";
          count++;
          // if(h == 0 && small_adjust > 0) std::cout << small_adjust << ", " << offset << ", " << step_occupation << std::endl;
            if(count == conduction_electrons) { 
            std::cout << "potential break early at: " << epsilon << " instead of: " << min << std::endl;
            full = true;
            break;}
        }
        // energy_step++;
        epsilon -= step_size;
      }
    }
    core_cutoff = epsilon;
    if(count != conduction_electrons) {
      terminaltextcolor(RED);
      std::cout << "potential early break: " << count << " rather than " << conduction_electrons << std::endl;
    }
    // E_f_A = epsilon;
  } else {
    // double u_Ni = 0.0;
    mu_Ni = 0.0;//(temp < 400.0) ? 2.42*temp/1250.0 : 2.42*(1.0+ temp-400.0)/225.0;
    epsilon = E_f_A + 0.0000001;// + step_size*energy_step;
    for(int h = 0; h < hist_super_E_f; h++) {
     
      const double DoS_occupation = dos_en_step* dos_standard[int(floor((epsilon-min)*i_dos_en_step))]*(return_fermi_distribution(epsilon-E_f_A-mu_Ni,temp)+return_fermi_distribution(epsilon+step_size-E_f_A-mu_Ni,temp)+4.0*return_fermi_distribution(epsilon+0.5*step_size-E_f_A-mu_Ni,temp))/6.0;
      if(DoS_occupation < 1) break;
      const int step_occupation = int(floor(DoS_occupation/double(steps)));
      int small_adjust = int(round(DoS_occupation - step_occupation*steps)) % steps;
      int offset = 0;
      for(int s = 0; s < steps; s++) {
        if(small_adjust > 0) offset = 1;
        else offset = 0;
        small_adjust--;
        for(int o = 0; o < (step_occupation + offset); o++) {
          distribution.at(count) = epsilon;
          distrib << count << ", " << epsilon << "\n";
          count++;
         
          if(count == conduction_electrons) { 
            std::cout << "potential break early at: " << epsilon << " instead of: " << min << std::endl;
         
            break;}
        }
        // energy_step++;
        epsilon += step_size;
      }
    }
    epsilon = E_f_A - 0.000001;
    for(int h = 0; h < hist_sub_E_f && !full; h++) {
     
      const double DoS_occupation = dos_en_step* dos_standard[int(floor((epsilon-min)*i_dos_en_step))]*(return_fermi_distribution(epsilon-E_f_A-mu_Ni,temp)+return_fermi_distribution(epsilon-step_size-E_f_A-mu_Ni,temp)+4.0*return_fermi_distribution(epsilon-0.5*step_size-E_f_A-mu_Ni,temp))/6.0;
      const int step_occupation = int(floor(DoS_occupation/double(steps)));
      int small_adjust = int(round(DoS_occupation - step_occupation*steps)) % steps;
      int offset = 0;
      for(int s = 0; s < steps && !full; s++) {
        if(small_adjust > 0) offset = 1;
        else offset = 0;
        small_adjust--;
        for(int o = 0; o < (step_occupation + offset); o++) {
          distribution.at(count) = epsilon;
          distrib << count << ", " << epsilon << "\n";
          count++;
           if(count == conduction_electrons) { 
              std::cout << "potential break early at: " << epsilon << " instead of: " << min << "; conduction electrons: " << count << std::endl;
           full = true;
           break;}
        }
        // energy_step++;
        if(return_fermi_distribution(epsilon - E_f_A-mu_Ni, temp) < 1.0-1e-3) transport_cutoff = epsilon;
        epsilon -= step_size;
        if(epsilon >= zero_temp) finite_electron_number = count;
      }
    }
  }
  transport_cutoff -= 2.0;
  DoS_cutoff =  min;
  core_cutoff = epsilon;

  // transport_cutoff = E_f_A;
  std::cout << "Discrete DoS_cutoff at " << DoS_cutoff << ", offset " << floor((core_cutoff-DoS_cutoff)*i_dos_en_step) << std::endl;
  std::cout << "mu cutoff error gives " << conduction_electrons - count << " electrons. % = " << 100.0*double(conduction_electrons-count)/double(conduction_electrons) << \
  ". Estimated mu_feg (AJ): " <<  -1.0* pow(M_PI*constants::kB_r*temp/E_f_A, 2.0)/12.0 << ", used: " << mu_Ni << std::endl;
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
double return_fermi_distribution(const double energy, const double temp) 
{
  if(temp <= 0.00001) return (1.0/(exp(energy/(constants::kB_r*1.0)) + 1.0));
  else return (1.0/(exp(energy/(constants::kB_r*temp)) + 1.0));
}

double return_BE_distribution(const double phonon_e, const double temperature) {
  return 1.0/(exp(phonon_e/(temperature*constants::kB_r)) - 1.0);
}

double return_dWdE(const double e_energy) {
	 return 0.0013937*(e_energy-E_f_A) + 2.3893;
//   if(e_energy > (E_f_A+4.86166)) return 0.00275775*(e_energy-E_f_A) + 3.81605;
//   else return 0.0244961*(e_energy-E_f_A) + 3.71037;
}

double return_dWdE_i(const double momentum) {
	 return E_f_A + ((momentum-2.3893)/0.0013937);
//   if(momentum > 3.82946) return E_f_A + ((momentum-3.81605)/0.00275775);
//   else return E_f_A + ((momentum- 3.71037)/0.0244961);
}

double return_vel(const double energy) {
	 return 1.0/(constants::hbar_r*0.0013937);
//   if(energy > (E_f_A+4.86166)) return 1.0/(constants::hbar_r*0.00275775);
//   else return 1.0/(constants::hbar_r*0.0244961);
}

double return_m_e_r(const double energy) {
  if(energy > (E_f_A+4.86166)) return constants::hbar_r*constants::hbar_r*0.00275775*(0.00275775*(energy-E_f_A) + 3.81605);
  else return  (0.0244961*(energy-E_f_A) + 3.71037)*constants::hbar_r*0.0244961*constants::hbar_r;
}

double return_DoS_phonon(const double energy) {
  return 9.0*energy*energy/(phonon_energy*phonon_energy*phonon_energy);
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
    // double local_d_U = 0.0;
    double d_U_avg = 0.0;
    double transient_entropy = 0.0;
    int ballistic_count = 0;
    int count = 0;
    #pragma omp parallel for schedule(dynamic,4)  reduction(+:e_size, scat_size, d_U_avg, count)
    for(int e = 0; e < conduction_electrons; e++) {
      double e_energy = electron_potential[e];
      
      if(e_energy >= 84.181114) {d_U_avg += e_energy*dos_en_step; count++;}
      // int index = std::min(dos_size-1.0, std::max(0.0, floor((e_energy-DoS_cutoff)*i_dos_en_step)));
      // if(e_energy < E_f_A)  local_d_U += (E_f_A - e_energy)*dos_en_step*abs(1.0 - double(ee_dos_hist[e][index])/(dos_standard[index]*dos_en_step));
      // else local_d_U += (e_energy - E_f_A)*dos_en_step*double(ee_dos_hist[e][index])/(dos_standard[index]*dos_en_step);

      e_size += electron_integration_list[e][0];
      scat_size += electron_ee_scattering_list[e][1];

      if(e_energy  > (E_f_A+4.865)) {
        #pragma omp atomic 
        ballistic_count++;
      }
    }

    if(ballistic_count > 0) std::cout << "Ballistic count: " << ballistic_count << std::endl;
    // local_d_U /= lattice_atoms* (dos_standard[int(floor((E_f_A-DoS_cutoff)*i_dos_en_step))]/lattice_atoms); // normalise by volume (e- number) then D(E_f) and DoS width for sum
    // global_d_U /= conduction_electrons/(E_f_A-core_cutoff);
    e_size /= conduction_electrons;
    scat_size /= conduction_electrons;
    d_U_avg /= count; // summation scaled by DoS resolution; divided by D(E_f) after zero point reduction
    // d_U_avg *= 3.0*constants::eV_to_AJ; //(lattice_depth*lattice_width*lattice_height);

      if(err::check)  std::cout << "local dU calcualted " << std::endl;

    #pragma omp parallel for schedule(dynamic, 4)  reduction(+:e_stddev,scat_stddev)
    for(int e = 0; e < conduction_electrons; e++) {
      // if(electron_potential[e] > 0.8*E_f_A){
        e_stddev += (electron_integration_list[e][0] - e_size)   *(electron_integration_list[e][0] - e_size);
        scat_stddev += (electron_ee_scattering_list[e][1] - scat_size)*(electron_ee_scattering_list[e][1] - scat_size);
      // }
    }
     if(err::check)  std::cout << "stdDev occ calcualted " << std::endl;
    // int|E_f window : (1-f_i)(E_f - E_i)D(e_i)de -> Reiman sum with step size: dos_en_step

    double tau_ep = 0.0;
    double tau_ee = 0.0;
    for(double u_e = core_cutoff; u_e < (E_f_A+max_as); u_e += dos_en_step) {
      
      const int index = std::min(dos_size-1.0, floor((u_e-DoS_cutoff)*i_dos_en_step));
      tau_ep += global_tau_ep[index];
      tau_ee += global_tau_ee[index];
      // global_tau_ep[index] = 0.0;
      // global_tau_ee[index] = 0.0;
      if(u_e < E_f_A) {
        if(u_e > transport_cutoff) {
          double dos = dos_standard[index]/lattice_atoms;
          double fermi_dist = std::min(1.0, double(global_e_dos[index][0])/(dos_standard[index]*dos_en_step));
          double onemin_fermi_dist = 1 - fermi_dist;
          global_d_U += (E_f_A - u_e)*dos_en_step*abs(1.0 - fermi_dist)*dos; 
          transient_entropy += dos_en_step*dos*(fermi_dist*log(fermi_dist) + ((onemin_fermi_dist == 0) ? 0.0 : onemin_fermi_dist*log(onemin_fermi_dist)));
            // std::cout << dos_en_step*dos*(fermi_dist*log(fermi_dist) + ((onemin_fermi_dist == 0) ? 0.0 : onemin_fermi_dist*log(onemin_fermi_dist))) << ", " << (fermi_dist*log(fermi_dist) + ((onemin_fermi_dist == 0) ? 0.0 : onemin_fermi_dist*log(onemin_fermi_dist))) << ", " << fermi_dist*log(fermi_dist) << ", " << ((onemin_fermi_dist == 0) ? 0.0 : onemin_fermi_dist*log(onemin_fermi_dist)) << std::endl; 
        } else {
          double fermi_dist = std::min(1.0, (double(global_e_dos[index][0])/std::max(1.0,double(global_e_dos[index][1]))));
          double dos = double(global_e_dos[index][1])/lattice_atoms/dos_en_step;
          double onemin_fermi_dist = 1 - fermi_dist;
          global_d_U += (E_f_A - u_e)*dos_en_step*abs(1.0 - fermi_dist)*dos; 
          transient_entropy += dos_en_step*dos*(fermi_dist*log(fermi_dist) + ((onemin_fermi_dist == 0) ? 0.0 : onemin_fermi_dist*log(onemin_fermi_dist)));
          //  std::cout << dos_en_step*dos*(fermi_dist*log(fermi_dist) + ((onemin_fermi_dist == 0) ? 0.0 : onemin_fermi_dist*log(onemin_fermi_dist))) << ", " << (fermi_dist*log(fermi_dist) + ((onemin_fermi_dist == 0) ? 0.0 : onemin_fermi_dist*log(onemin_fermi_dist))) << ", " << fermi_dist*log(fermi_dist) << ", " << ((onemin_fermi_dist == 0) ? 0.0 : onemin_fermi_dist*log(onemin_fermi_dist)) << std::endl;
        }
      } else {
        double fermi_dist = std::min(1.0, double(global_e_dos[index][0])/(dos_standard[index]*dos_en_step));
        double dos = dos_standard[index]/lattice_atoms;
        double onemin_fermi_dist = 1 - fermi_dist;
        global_d_U += (u_e - E_f_A)*dos_en_step*fermi_dist*dos;
        transient_entropy += dos_en_step*dos*(((fermi_dist == 0) ? 0.0 : fermi_dist*log(fermi_dist)) + ((onemin_fermi_dist == 0) ? 0.0 : onemin_fermi_dist*log(onemin_fermi_dist)));
        //  std::cout << fermi_dist << ", " << dos_en_step*dos*(fermi_dist*log(fermi_dist) + ((onemin_fermi_dist == 0) ? 0.0 : onemin_fermi_dist*log(onemin_fermi_dist))) << ", " << (fermi_dist*log(fermi_dist) + ((onemin_fermi_dist == 0) ? 0.0 : onemin_fermi_dist*log(onemin_fermi_dist))) << ", " << fermi_dist*log(fermi_dist) << ", " << ((onemin_fermi_dist == 0) ? 0.0 : onemin_fermi_dist*log(onemin_fermi_dist)) << std::endl;
      } 
    } 


    tau_ep /= double(conduction_electrons);
    tau_ee /= double(conduction_electrons);
    e_stddev = sqrt(e_stddev/double(conduction_electrons));
    scat_stddev = sqrt(scat_stddev/double(conduction_electrons));

          if(err::check)  std::cout << "global dU occ done " << std::endl;
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

    // std::cout << "why?" << std::endl;
     double normalise = 1.0/(double(CASTLE_output_rate * CASTLE_MD_rate)*lattice_depth*lattice_height*lattice_width);
      for(int e = 0; e < flux_index.size(); e++) {
        flux_hist << e*dos_en_step + DoS_cutoff << ", " << flux_index[e]*normalise << "\n";
        flux_index[e] = 0;
      }
      flux_hist.close();

    // stopwatch_t omp_time;
    // omp_time.start();
    double ee_avg = 0.0;
    int ee_total = 0;
    #pragma omp parallel 
    {
      std::vector<double> ee_avg_array;
      std::vector<int> ee_total_array;
      ee_avg_array.resize(relaxation_time_hist_ee[0].size(), 0.0);
      ee_total_array.resize(relaxation_time_hist_ee[0].size(), 0);
      #pragma omp for 
      for(int e = 0; e < conduction_electrons; e++) { 
        for(int h = 0; h < relaxation_time_hist_ee[0].size(); h++) {
        
          ee_avg_array[h] += double(relaxation_time_hist_ee[3*e + 2][h])/std::max(1.0, double(relaxation_time_hist_ee[3*e][h]) );
          ee_total_array[h] += relaxation_time_hist_ee[3*e][h];
        }
      }

      for(int h = 0; h < dos_size ; h++) {

      if(h < relaxation_time_hist_ee[0].size()) {
        #pragma omp critical
        {
        ee_avg += ee_avg_array[h];
        ee_total += ee_total_array[h];
        ee_avg_array[h] = 0;
        ee_total_array[h] = 0;
        }
      }
        #pragma omp barrier

        #pragma omp single 
        {
        if(h < dos_size) {
          relaxation_time << double(h)/4.0 + int(round(DoS_cutoff)) << ", " << ee_avg/std::max(1.0,double(ee_total)) << ", " << ee_total << ", ";
          ee_avg = 0.0;
          ee_total = 0;
        }
        relaxation_time << h*dos_en_step + DoS_cutoff << ", " << global_tau_ep[2*h] << ", " << global_tau_ep[2*h+1] << ", " << global_tau_ee[h] << std::endl;
        global_tau_ep[h] = 0.0;
        global_tau_ee[h] = 0.0;
        }
      }
    }
     if(err::check)  std::cout << "relaxation time output " << std::endl;
    // std::cout << "relaxatime output timing: " << omp_time.elapsed_seconds() << std::endl;
    relaxation_time.close();
      // std::cout << "why not" << std::endl;
    // }
  //  const int output_count_lr = int(round(transport_cutoff-core_cutoff));
    const int output_count_hr = global_e_dos[0].size();
   
    for(int i = 0; i < dos_size; i++) {
     // if(i == 11) temp_map_0 << i << ", " << temp_Map[0].at(i) << "\n";
     // if(i < output_count_lr) { 
      temp_map << i*dos_en_step + DoS_cutoff << ", " << global_e_dos[i][0] << ", " << global_e_dos[i][1] << ", " << dos_standard[i]*dos_en_step << "\n";
      if(global_e_dos[i][0]/(dos_standard[i]*dos_en_step) > 1.02) {
        terminaltextcolor(RED);
        std::cout << "DoS violation at " << i*dos_en_step+DoS_cutoff << ", " << global_e_dos[i][0] << " over " << dos_standard[i]*dos_en_step << std::endl;
        terminaltextcolor(WHITE);
      }
    }
    // std::cout << "why not" << std::endl;
    temp_map.close();
     if(err::check)  std::cout << "occupation output " << std::endl;

      std::cout << std::fixed; std::cout.precision(3); std::cout << "  " << current_time_step / total_time_steps * 100 << "%. Lifetime energy leak: " << full_int_var << std::endl; 
    }
    mean_data.precision(10);
    mean_data << std::scientific;

    const double I = double(x_flux) * constants::e * 1e35 / double(CASTLE_output_rate) / dt / (lattice_height*lattice_depth); //current density/m**2
    // double u_Ni = (Te < 400.0) ? 2.42*Te/1250.0 : 2.42*(1.0+ Te-400.0)/225.0;
    global_d_U /= (dos_standard[int(floor((E_f_A -DoS_cutoff)*i_dos_en_step))]/lattice_atoms);
    if(!current_time_step) {
    mean_data << CASTLE_real_time << ", " << current_time_step << ", " 
      << sqrt(6.0*global_d_U)/M_PI/constants::kB_r << ", "  << sqrt(6.0*global_d_U)/M_PI/constants::kB_r << ", " <<  global_d_U << ", "// (d_Te*d_Te*e_heat_capacity +Tp*a_heat_capacity) << ", " 
      << Te << ", " << Tp << ", " << -1.0*transient_entropy << ", " << total_TEKE << ", " 
      << d_TTMe << ", " << d_TTMp << ", " <<  I*double(CASTLE_output_rate) << ", " << e_size << ", " << e_stddev << ", " << scat_size << ", " << scat_stddev << ", " << p_x/double(conduction_electrons) << ", " << p_y/double(conduction_electrons) << ", " << p_z/double(conduction_electrons) << ", " 
      << std::fixed; mean_data.precision(1); mean_data << double(e_a_scattering_count) / 1 << ", " << double(e_e_scattering_count) / double(1) << ", " << \
      double(ee_core_scattering_count) / double(1) << ", " << double(ee_transport_scattering_count) / double(1) << ", " <<\
      double(ea_core_scattering_count) / double(1) << ", " << double(ea_transport_scattering_count) / double(1) << ", " <<\
      double(x_flux) / double(1) << ", " << double(y_flux) / 1 << ", " << double(z_flux) / double(1)  << ", " \
      << std::endl;
    } else {
    mean_data << CASTLE_real_time << ", " << current_time_step << ", " 
      << sqrt(6.0*global_d_U)/M_PI/constants::kB_r << ", " << sqrt(6.0*global_d_U)/M_PI/constants::kB_r << ", " << global_d_U << ", " //<<  (d_Te*d_Te*e_heat_capacity +Tp*a_heat_capacity) << ", " 
      << d_Te << ", " << d_Tp << ", " << -1.0*transient_entropy << ", " << total_TEKE << ", " 
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
    // if(!equilibrium_step) half_int_var = 1;
    
    // if(equilibrium_step) transport_cutoff = E_f_A - (E_f_A - core_cutoff)*current_time_step/sim::equilibration_time;
    if(transport_cutoff > (E_f_A - 0.5*floor((d_TTMe - 300.0)/100.0))) {
      std::cout << "transport cutoff shift from " << transport_cutoff << ", to: ";
      transport_cutoff = (E_f_A - 0.5*floor((d_TTMe - 300.0)/100.0));
      std::cout << transport_cutoff << std::endl;
    }
    // if(Te > 900) transport_cutoff = 102.951;
  
}




} //end of CASTLE namespace



