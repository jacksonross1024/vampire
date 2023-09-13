// ==============================================
//   Coupled Atomistic Spin Thermalized Lattice Environment
//
//  =========    ========     ========  ============  ||           =========
// ||          ||        ||  ||              ||       ||          ||
// ||          ||        ||  ||              ||       ||          ||
// ||          ||        ||  ||              ||       ||          ||
// ||          || ====== ||   ========       ||       ||           =========
// ||          ||        ||          ||      ||       ||          ||
// ||          ||        ||          ||      ||       ||          ||
// ||          ||        ||          ||      ||       ||          ||
//  =========                 ========                  =========   =========
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
 //ouble return_phonon_distribution(const double& epsilon, const double& beta );


int velocity_verlet_step(double time_step) {
  
    TEKE = 0.0;
    TLE  = 0.0;
    p_x = 0.0;
    p_y = 0.0;
    p_z = 0.0;


    int total_1 = 0;
    int total_2 = 0;
    for(int h = 0 ; h < dos_size; h++) {
      //  std::cout << h << ", " << global_e_dos[h][0] << ", " << global_e_dos[h][1] << std::endl;
      total_1 += global_e_dos[h][0];
      total_2 += global_e_dos[h][1];
      global_e_dos[h][0] = 0;
    }
          if(err::check)  std::cout << "Updating new electron position." << std::endl;
   update_position();

         if(err::check) std::cout << " positions updated. updating dos " << std::endl;
    update_dynamics();
         if(err::check)  std::cout << "Output mean data" << std::endl;
    
    if (current_time_step % CASTLE_output_rate == 0) output_data(); //std::cout << "x_flux: " << x_flux / CASTLE_output_rate << "\n"; x_flux = 0;
          if(err::check)  std::cout << "resetting global dos " << std::endl;

   // const int size = ee_dos_hist[0].size();
    #pragma omp parallel for schedule(dynamic, 4)
    for(int e = 0; e < conduction_electrons; e++) {
      // if(electron_potential[e] < E_f_A-24.25) continue;
      electron_ee_scattering_list[e][0] = 0;
      // for(int h = 0; h < size; h++) {
      //   ee_dos_hist[e][h] = 0;
      // }
    }

    
    if(total_1 != conduction_electrons) std::cout << "hist problem " << total_1 << " != " << conduction_electrons << " != " << total_2 << std::endl;
    
    //reset integration
    CASTLE_real_time += dt;
    current_time_step++;
   
    return EXIT_SUCCESS;
}


void update_position() {

    omp_set_dynamic(0);
    omp_set_num_threads(omp_threads);
    
    if(current_time_step % half_int_var[1] == 0) {
      cell_integration_lists.swap(old_cell_integration_lists);
     
      for(int c = 0; c < total_cells; c++) {
        cell_integration_lists[c][0] = 1;
      }
        if(err::check) std::cout << "cell lists swapped.";
    }
     if(err::check) std::cout << " updating position..." << std::endl;
    int total = 0;
    int total_1 = 0;
  
  #pragma omp parallel reduction(+:x_flux,y_flux,z_flux, p_x,p_y,p_z)
  {
    std::vector<int> local_e_dos;
    local_e_dos.resize(dos_size,0);
    std::vector<double> local_x_dos;
    local_x_dos.resize(dos_size,0.0);
    const int thread = omp_get_thread_num();
    //if(omp_get_num_threads() != omp_threads) std::cout << "omp pragma error " << omp_get_num_threads() << " != " << omp_threads << std::endl;
    int local = 0;
    for (int l = 0 ; l < cells_per_thread; l++) {
      const int cell = lattice_cells_per_omp[thread][l];
      // #pragma omp critical
      // std::cout << "cell " << cell << " done by thread " << thread << std::endl;
      int size;
      if(current_time_step % half_int_var[1] == 0) size = old_cell_integration_lists[cell][0];
      else  size = cell_integration_lists[cell][0];

      std::vector<int> cell_integration_list;
      cell_integration_list = (current_time_step % half_int_var[1] == 0) ? old_cell_integration_lists[cell] : cell_integration_lists[cell];
      local += size - 1;
  //  if((current_time_step % half_int_var) == 0 && ( (l % ((z_omp_cells/2)-1) )== 0 || (l % (z_omp_cells/2)) == 0) ) {
      #pragma omp barrier
  
      for (int e = 1; e < size; e++) { 
         int electron = cell_integration_list[e];
         const double energy = electron_potential[electron];
         int index = int(std::min(dos_size-1.0, std::max(0.0, floor((energy - DoS_cutoff)*i_dos_en_step))));
         
   
         const double vel = return_vel(energy);
         const double mom = return_dWdE(energy);
         local_e_dos[index]++;
         const int array_index = 3*electron;

         double pos[3] = {electron_position[array_index], electron_position[array_index+1], electron_position[array_index+2]};

         double v[3] = {electron_velocity[array_index], electron_velocity[array_index+1], electron_velocity[array_index+2]};

        // if(electron_potential[e] > 0.98*E_f_A) {
          p_x += v[0]*mom;
          p_y += v[1]*mom;
          p_z += v[2]*mom;
          pos[0] += v[0]*vel * dt;// + (electron_force[array_index]   * dt * dt * constants::K_A / 2); // x superarray component
          pos[1] += v[1]*vel * dt;// + (electron_force[array_index_y) * dt * dt * constants::K_A / 2); // y superarray component
          pos[2] += v[2]*vel * dt;// + (electron_force[array_index_z) * dt * dt * constants::K_A / 2); // z superarray component
        // }
         //  if(v_x != v_x || v_y != v_y || v_z != v_z) std::cout << energy << ", " << vel << ", " << v_x << ", " << v_y << ", " << v_z << std::endl;
        if (pos[0] < 0.0) {
          pos[0] += lattice_width; 
          x_flux--;
        } else if (pos[0] > lattice_width) {
          pos[0] -= lattice_width; 
          x_flux++; 
        }
      
        local_x_dos[index] += electron_velocity[array_index]*return_dWdE(energy);
          
        if (pos[1] < 0.0) {pos[1] += lattice_depth; y_flux--;}
        else if (pos[1] > lattice_depth) {pos[1] -= lattice_depth; y_flux++;}

        if (pos[2] < 0.0) {pos[2] += lattice_height; z_flux--;}
        else if (pos[2] > lattice_height) {pos[2] -= lattice_height; z_flux++;}

        electron_position[array_index]   = pos[0];
        electron_position[array_index+1] = pos[1];
        electron_position[array_index+2] = pos[2];
    
      if(current_time_step % half_int_var[1] == 0) {
          
           int x_cell = int(floor(pos[0] / x_step_size));

           int y_cell = int(floor(pos[1] / y_step_size));
    
           int z_cell = int(floor(pos[2] / z_step_size));

          const int omp_cell = lattice_cell_coordinate[x_cell][y_cell][z_cell];
          
          if ((abs(x_cell - cell_lattice_coordinate[cell][0]) > 1 &&\
               abs(x_cell - cell_lattice_coordinate[cell][0]) <= (x_omp_cells-2)) ||
              (abs(y_cell - cell_lattice_coordinate[cell][1]) > 1  &&\
               abs(y_cell - cell_lattice_coordinate[cell][1]) <= (y_omp_cells-2)) ||\
              (abs(z_cell - cell_lattice_coordinate[cell][2]) > 1  &&\
               abs(z_cell - cell_lattice_coordinate[cell][2]) <= (z_omp_cells-2) )) {
            #pragma omp critical(halo)
            {
            escaping_electrons[escaping_electrons[0]] = electron;
            escaping_electrons[0]++;
            }
          } else if (omp_cell != cell) {
              #pragma omp critical(shell)
              {
              cell_integration_lists[omp_cell][cell_integration_lists[omp_cell][0]] = electron;
              cell_integration_lists[omp_cell][0]++;
              }
          } else {
            cell_integration_lists[cell][cell_integration_lists[cell][0]] = electron;
            cell_integration_lists[cell][0]++;
           // if(cell_integration_lists[cell][0] >= cell_integration_lists[cell).size() - 1) std::cerr << "too big? " << cell_integration_lists[cell][0] << " > " << cell_integration_lists[cell).size() << std::endl;
           }
        }
      }
    }
        
        if(err::check && omp_get_thread_num()==0) std::cout << "positions updated. updating global dos..." << std::endl;
   #pragma omp critical 
   {
      total += local;

      for (int h = 0; h < dos_size; h++) {
         flux_index[h] += local_x_dos[h];
         global_e_dos[h][0] += local_e_dos[h];
         total_1 += local_e_dos[h];
      }
   } 
   }
      if(err::check) std::cout << "global dos updated. catching escaped electrons: " << escaping_electrons[0] -1<< std::endl;
  if(total != conduction_electrons || total_1 != conduction_electrons) std::cout << total << ", " << conduction_electrons << std::endl;
   if(current_time_step % half_int_var[1] == 0) {    
    const int size = escaping_electrons[0];
  //  std::cout << "escaping electrons " << size-1 << ", " << escaping_electrons.size() << std::endl;
    
    for (int e = 1; e < size; e++) {
      const   int electron = escaping_electrons[e];
      const   int array_index = 3*electron;
      int x_cell = int(floor(electron_position[array_index] / x_step_size));

      int y_cell = int(floor(electron_position[array_index+1] / y_step_size));
    
      int z_cell = int(floor(electron_position[array_index+2] / z_step_size));

      //      if(x_cell >= x_omp_cells || y_cell >= y_omp_cells || z_cell >= z_omp_cells) {std::cout << "cell sorting for ee integration exceeds bounds" << \
      // x_cell << ", " << y_cell << ", " << z_cell << std::endl;
      // x_cell --;
      // y_cell --;
      // z_cell --; }
      //     if(x_cell < 0 || y_cell < 0 || z_cell < 0) {std::cout << "cell sorting for ee integration less than zero" << \
      // x_cell << ", " << y_cell << ", " << z_cell << std::endl;
      // x_cell = 0;
      // y_cell = 0;
      // z_cell = 0; }
      
        const int omp_cell = lattice_cell_coordinate[x_cell][y_cell][z_cell];
        cell_integration_lists[omp_cell][cell_integration_lists[omp_cell][0]] = electron;
        cell_integration_lists[omp_cell][0]++;

      //  escaping_electrons[e] = 0;
    }
    escaping_electrons[0] = 1;
  }
}

void update_dynamics() {
  

    //I(t) = 0.94 F/(t_p * L_z) AJ/nm^2 / fs / A
    const double intesnity = 0.94e-14*fluence/sim::pump_time/lattice_height; //AJ/nm^2/A/s -> e1/e15 AJ/nm^3/fs
    int photons_at_dt = 0; //hv*dt
    double pump = 0.0; // AJ / fs
    double external_potential = 0.0; //AJ/e-
    const static double sigma = 1e-15 / sim::pump_time; // 1/s -> 1e-15 -> 1/fs
	  const double four_ln_2 = 2.77258872224;

    int count = 0;
    // std::cout << sigma << std::endl;
    if(!equilibrium_step) {
      if(heat_pulse_sim) {
        //hv(dt)/fs
        pump = intesnity*exp(-four_ln_2*sigma*sigma*((double(current_time_step) - ((40.0/dt)+sim::equilibration_time))*(double(current_time_step) - ((40.0 / dt)+sim::equilibration_time)))); // AJ/fs/nm**3
        photons_at_dt = 1e-3*pump*dt*lattice_depth*lattice_height*lattice_width/photon_energy; //AJ/fs/nm^3 * A^3(1e-3 nm^3/A^3) * dt -> AJ / AJ/gamma
        //external_potential = photon_energy; //AJ/hv     ;//1e27*pump*dt/ n_f; // AJ / particle
     //   std::cout << photons_at_dt << std::endl;
        TTMe = d_TTMe;
        TTMp = d_TTMp;
          //AJ/fs/K -> g(T-T)=C/t
        d_TTMe = ((G*(TTMp - TTMe)+pump)*dt*e_heat_capacity_i/TTMe) + TTMe;
      //  std::cout <<" TTMe: " << pump*dt*e_heat_capacity_i*300.0/TTMe << ", " << pump*dt*e_heat_capacity_i*300.0/Te;
        d_TTMp = ( G*(TTMe - TTMp)      *dt*a_heat_capacity_i)      + TTMp;
      }
    }
   double potential = 0.0;
      if(!equilibrium_step && applied_voltage_sim) {
         //
         potential = 1e10*dt*constants::e*applied_voltage *std::min(((current_time_step-sim::equilibration_time) /10000.0), 1.0)/constants::hbar_r;
         // std::cout << potential << ", " << (current_time_step-sim::equilibration_time) /100.0 << std::endl;
      } 

    omp_set_dynamic(0);
    omp_set_num_threads(omp_threads);
    std::vector<int> chosen;
    chosen.resize(photons_at_dt);
    while(count < photons_at_dt) {
      int select = int(omp_uniform_random[omp_get_thread_num()]()*2147483647) % (conduction_electrons);
      double e_energy = electron_potential[select];
      if(e_energy > (DoS_cutoff+2.0*dos_en_step) && return_fermi_distribution(e_energy+photon_energy-E_f_A, Te) < (1.0-1e-2) && (e_energy+photon_energy < E_f_A+ max_as)) {
        chosen[count] = select;
        count++;
      }
    }
    count = 0;
      //  pump = 0.0;
    Tp = d_Tp;
    Te = d_Te;
      if(err::check) std::cout << "conditions set; updating scattering lists" << std::endl;

  //  pump = 0.0;

    #pragma omp parallel for schedule(dynamic, 4) reduction(+:external_potential)
    for (int e = 0; e < conduction_electrons; e++) {
      // if(electron_potential[e] < (E_f_A-24.25)) continue;
      const int array_index = 3*e;
      
      if(current_time_step % half_int_var[(electron_potential[e]<(E_f_A+4.8) ? 0 : 1)] == 0) e_e_coulomb(e, array_index);
      // else if (current_time_step % half_int_ == 0 && electron_transport_list[e]) e_e_coulomb(e, array_index);
      else  neighbor_e_e_coulomb(e, array_index);

      if(photons_at_dt > 0 && std::end(chosen) != std::find(chosen.begin(), chosen.end(), e)) {
        // #pragma omp atomic
        // count++;
        
        // pump += external_potential;
        electron_thermal_field(e, array_index, photon_energy, omp_get_thread_num());
      }
      
      if(!equilibrium_step && applied_voltage_sim) external_potential += electron_applied_voltage(e, array_index, potential);
      //if(equilibrium_step && heat_pulse_sim) ea_scattering(e, array_index, omp_get_thread_num());
      elastic_ea_scattering(e, array_index, omp_get_thread_num());
    }

    // if(count != photons_at_dt) std::cout << photons_at_dt << ", " << count<< std::endl;
       TEKE += external_potential;
        if(err::check) std::cout << "lists updated. ee scattering step..." << std::endl;
     // q_sq = k_sq();
     
    ee_scattering();
   // pump /= 1e-3*lattice_depth*lattice_height*lattice_width;  
   //  if(!equilibrium_step) 
    d_Tp = a_heat_capacity_i*TLE *n_f*total_e_scaling/conduction_electrons + Tp;
   //  else d_Tp = sim::temperature;

    d_Te = (e_heat_capacity_i*TEKE*n_f*total_e_scaling/conduction_electrons/abs(Te)) + (dt*e_heat_capacity_i*pump/abs(Te)) + Te;
 
   total_TEKE += TEKE;
        if (err::check) std::cout << "reset scattering." << std::endl;
       
}

void electron_thermal_field(const int e, const int array_index, const double EKE, const int thread) {       
    
   const double theta = omp_uniform_random[thread]() * 2.0 * M_PI;
   const double phi   = omp_uniform_random[thread]() * M_PI; 

   electron_potential[e] += EKE;
  
   electron_velocity[array_index]   = cos(theta)*sin(phi);
   electron_velocity[array_index+1] = sin(theta)*sin(phi);
   electron_velocity[array_index+2] = cos(phi);

   electron_ea_scattering_list[e][0] = 1;
   return;
}

void e_e_coulomb(const int e, const int array_index) {

   const int cell = lattice_cell_coordinate[int(floor(electron_position[array_index] / x_step_size))][int(floor(electron_position[array_index+1] / y_step_size))][int(floor(electron_position[array_index+2] / z_step_size))];
   int ee_dos_count = 1;
   int ee_integration_count = 1;
   int ee_scattering_list = 1;

   double integration_range = e_e_integration_cutoff; //(electron_potential[e] > E_f_A+4.8) ? e_e_integration_cutoff*2.0: e_e_integration_cutoff;
   int i_size = 3*round(pow(integration_range,1.5)*1.25*M_PI * total_e_scaling*n_f * 1e-3);
   // if(electron_ee_scattering_list[e].size() < e_size) electron_ee_scattering_list[e].resize(e_size, 0);
   // if(electron_integration_list[e].size() < i_size) {
   //    // std::cout << "resize. energy better be > " << electron_potential[e] << " > " << E_f_A+4.8 << std::endl;
   //    electron_integration_list[e].resize(i_size, 0);
   // }
   int cells = (i_size > ee_density) ? 125 : 27;
  //  if(cells > 27) std::cout << i_size << ", " << ee_density << ", " << electron_potential[e] << std::endl;
   std::vector<int> integration_list;
   integration_list = (cells > 27) ? cell_lr_neighbor_list[cell] : cell_nearest_neighbor_list[cell];

    for(int s = 0; s < cells; s++) {
      // if(s > 27) std::cout << "energy better be > " << E_f_A+4.8 << ": " << electron_potential[e] << std::endl;
      const int int_cell = integration_list[s];
      const int size = cell_integration_lists[int_cell][0];
      if(size != size || size < 1  ) std::cout << size << std::endl;
      for(int i = 1; i < size; i++) {
        const int electron = cell_integration_lists[int_cell][i];
        if(e == electron) continue;
        // if(electron_potential[electron] < (E_f_A-24.25)) continue;
        const int array_index_i = 3*electron;

        double distance[3] = {electron_position[array_index]  - electron_position[array_index_i], \
                         electron_position[array_index+1]- electron_position[array_index_i+1],\
                         electron_position[array_index+2]- electron_position[array_index_i+2]};
    
      	if (distance[0] < -1.0*boundary_conditions_cutoff)  distance[0]  += lattice_width;
   		else if (distance[0]  > boundary_conditions_cutoff) distance[0]  -= lattice_width;
   
   		if (distance[1]  < -1.0*boundary_conditions_cutoff) distance[1] += lattice_depth;
   		else if (distance[1] > boundary_conditions_cutoff)  distance[1] -= lattice_depth;    
   
   		if (distance[2] <  -1.0*boundary_conditions_cutoff) distance[2] += lattice_height;
   		else if (distance[2] > boundary_conditions_cutoff)  distance[2] -= lattice_height;
        
        double length  = (distance[0] *distance[0]  + distance[1]*distance[1] + distance[2]*distance[2]);
   		if(length == 0.0) return ;
        
      
		  if(length != length) std::cout << distance[0] << ", " << distance[1] << ", " << distance[2] <<  ", " << \
        electron_position[array_index] << ", " << electron_position[array_index_i] << ", " << electron_position[array_index+1] << ", " << electron_position[array_index_i + 1] << ", " <<\
        electron_position[array_index+2] << ", " <<  electron_position[array_index_i + 2] << std::endl;
        
        if(length > integration_range) continue;
        electron_integration_list[e][ee_integration_count] = array_index_i;
          if(ee_integration_count >= electron_integration_list[e].size() - 2) {std::cout << e << " ee integration " << ee_integration_count << " > " << electron_integration_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
          break; }
        ee_integration_count++;

      //   if(length > scattering_range) continue;
      //   electron_nearest_electron_list[e][ee_dos_count] = electron;
      //   ee_dos_count++;
      //   ee_dos_hist[e].at(int(std::max(0.0, std::min(dos_size-1.0, floor((electron_potential[electron]-DoS_cutoff)*i_dos_en_step)))))++;
      //   if(ee_dos_count >= electron_nearest_electron_list[e].size() - 2) {std::cout << e << " ee dos " << ee_dos_count << " > " << electron_nearest_electron_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
      //     break; }

        if (length > e_e_coulomb_cutoff) continue; 
        electron_ee_scattering_list[e][ee_scattering_list*2] = electron;
        electron_ee_scattering_list[e][ee_scattering_list*2 +1] = length;
        if(ee_scattering_list >= electron_ee_scattering_list[e].size() - 1) {std::cout << e << " ee scattering " << ee_scattering_list << " > " << electron_ee_scattering_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
          break; }

        ee_scattering_list++;
      }
    }
    electron_integration_list[e][0] = ee_integration_count;
   //  electron_nearest_electron_list[e][0] = ee_dos_count;
    electron_ee_scattering_list[e][1] = ee_scattering_list;
  //  if(e == 0 ) std::cout << ee_integration_count << ", " << ee_dos_count << ", " << ee_scattering_list << std::endl;
}

void neighbor_e_e_coulomb(const int e, const int array_index) {
    // if(electron_potential[e] < 0.8*E_f_A) return;
   const int size = electron_integration_list[e][0];
   int neighbor_count = 1;
   int scattering_count = 1;
   // int s_size = 3*round(pow(e_e_integration_cutoff,1.5)*1.25*M_PI * 3.8*n_f * 1e-3);
   // int l_size = 3*round(pow(400.0,1.5)*1.25*M_PI * 3.8*n_f * 1e-3);

   // double scattering_range = (electron_integration_list[e].size() > s_size) ? 400.0 : 49.0;
    for (int i = 1; i < size; i++) {
        const int array_index_i = electron_integration_list[e][i];
        const int electron = array_index_i/3;
        // if(electron_potential[electron] < (E_f_A-24.25)) continue;

        double distance[3] = {electron_position[array_index]  - electron_position[array_index_i], \
                         electron_position[array_index+1]- electron_position[array_index_i+1],\
                         electron_position[array_index+2]- electron_position[array_index_i+2]};
    
      	if (distance[0] < -1.0*boundary_conditions_cutoff)  distance[0]  += lattice_width;
   		else if (distance[0]  > boundary_conditions_cutoff) distance[0]  -= lattice_width;
   
   		if (distance[1]  < -1.0*boundary_conditions_cutoff) distance[1] += lattice_depth;
   		else if (distance[1] > boundary_conditions_cutoff)  distance[1] -= lattice_depth;    
   
   		if (distance[2] <  -1.0*boundary_conditions_cutoff) distance[2] += lattice_height;
   		else if (distance[2] > boundary_conditions_cutoff)  distance[2] -= lattice_height;
        
        double length  = (distance[0] *distance[0]  + distance[1]*distance[1] + distance[2]*distance[2]);
   		if(length == 0.0) return ;
        
       if (length != length) std::cout << distance[0] << ", " << distance[1] << ", " << distance[2] <<  ", " << \
        electron_position[array_index] << ", " << electron_position[array_index_i] << ", " << electron_position[array_index+1] << ", " << electron_position[array_index_i + 1] << ", " <<\
        electron_position[array_index+2] << ", " <<  electron_position[array_index_i + 2] << std::endl;
        
      //   if(length > scattering_range) continue;
      //   electron_nearest_electron_list[e][neighbor_count] = array_index_i/3;
      //   neighbor_count++;
      //   ee_dos_hist[e][int(std::max(0.0, std::min(dos_size-1.0, floor((electron_potential[electron]-DoS_cutoff)*i_dos_en_step))))]++;
      //   if(neighbor_count >= electron_nearest_electron_list[e].size() - 2) {std::cout << e << ", " << neighbor_count << " > " << electron_nearest_electron_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
      //   break; }

        if (length > e_e_coulomb_cutoff) continue; 
        electron_ee_scattering_list[e][scattering_count*2] = array_index_i/3;
        electron_ee_scattering_list[e][scattering_count*2 +1] = length;
        scattering_count++;
    }

   //  electron_nearest_electron_list[e][0] = neighbor_count;
    electron_ee_scattering_list[e][1] = scattering_count;
    //  if(e == 0 ) std::cout << size << ", " << neighbor_count << ", " << scattering_count << std::endl;
}

double electron_applied_voltage(const int e, const int array_index, const double potential) {

   const double energy = electron_potential[e];
   double k = return_dWdE(energy);
   double k_1 = electron_velocity[array_index]*k;
   double k_2 = electron_velocity[array_index+1]*k;
   double k_3 = electron_velocity[array_index+2]*k;

   double d_k_1 = k_1-potential; 
   // e-4  / A
   double d_k = sqrt(d_k_1*d_k_1 + k_2*k_2 + k_3*k_3);

   double d_energy = return_dWdE_i(d_k);
   double deltaE = d_energy - energy;
   
   if(d_energy > E_f_A+ max_as) return 0.0;
   else if(d_energy <= transport_cutoff) {
    
         electron_velocity[array_index]   *= -1.0;//cos(theta)*sin(phi);
         electron_velocity[array_index+1] *= -1.0;// sin(theta)*sin(phi);
         electron_velocity[array_index+2] *= -1.0;//cos(phi);
        //  electron_potential[e] = d_energy;
    return 0.0;
   } 
  // else if (d_energy > transport_cutoff) {
      const int d_index   = int(std::min( dos_size-2.0, std::max(1.0, floor((d_energy - DoS_cutoff)*i_dos_en_step))));
      double e_occupation   = std::min(1.0, double(global_e_dos[d_index][0]) / (dos_standard[d_index]*dos_en_step));
     // if(e_occupation < 0.99) {
         double theta = atan2(k_2, d_k_1);
         if(theta != theta) theta = 0.0;
         double phi = acos(k_3/d_k);
         electron_velocity[array_index]   = cos(theta)*sin(phi);
         electron_velocity[array_index+1] = sin(theta)*sin(phi);
         electron_velocity[array_index+2] = cos(phi);
         electron_potential[e] = d_energy;

         return deltaE;   
     // }
  // } 
  //  return 0.0;
}

void elastic_ea_scattering(const int e, const int array_index, const int thread) {

    const double e_energy = electron_potential[e];
    const double k = return_dWdE(e_energy);
        if(e_energy < core_cutoff+1.5) return;
    double e_occupation;
    double f_e_occupation;
    double r_e_occupation;
    const double phonon_factor = phonon_energy + (1.0*phonon_energy*mtrandom::gaussianc(omp_uniform_random[thread])/8.0); 
    double p = phonon_factor*1e5/(constants::hbar_r*3650.0); // 1e5/ A

    // #pragma omp critical 
    // std::cout << "phonon energy: " << phonon_energy << ", factor: " << phonon_factor << std::endl;
    //Gamma_ep =( V*pi^3 / hbar k) * \pm DoS_standard[e']/k' DoS_p (E_q)/q F M_ep^2
   // const double coeff = cell_volume*M_PI*M_PI*M_PI/(constants::hbar_r*k); //A^3 A/ AJ fs 
    const double p_energy = e_energy + phonon_factor;
    const double m_energy = e_energy - phonon_factor;
    const double p_k = return_dWdE(p_energy);
    const double m_k = return_dWdE(m_energy);
    // const double p_dos = dos_standard[int((p_energy-DoS_cutoff)*i_dos_en_step)]/p_k; //states A/AJ
    // const double m_dos = dos_standard[int((m_energy-DoS_cutoff)*i_dos_en_step)]/m_k;
    const double pho_s = dos_en_step*return_DoS_phonon(phonon_factor);///sqrt(2.0*M_PI); //
    const double M_ep = 0.5*constants::K* phonon_factor/(q_sq*q_sq + p*p);
    const double h_p = 0.5 + (p_k*p_k - p*p)/(2.0*k*k);
    const double h_m = 0.5 + (m_k*m_k - p*p)/(2.0*k*k);
    const double r_p = sqrt(p_k*p_k - h_p*h_p*k*k);
    const double r_m = sqrt(m_k*m_k - h_m*h_m*k*k);

    const double p_mom_ratio = 1.5*r_p/(p*p*p);
    const double m_mom_ratio = 1.5*r_m/(p*p*p);

    const int e_index   = int(std::min( dos_size-2.0, std::max(1.0, floor((e_energy - DoS_cutoff)*i_dos_en_step))));
    const int f_e_index = int(std::min( dos_size-2.0, std::max(1.0, floor((p_energy - DoS_cutoff)*i_dos_en_step))));
    const int r_e_index = int(std::min( dos_size-2.0, std::max(1.0, floor((m_energy - DoS_cutoff)*i_dos_en_step))));
    //  std::cout << e_index << ", " << f_e_index << ", " << r_e_index << std::endl;
    // const double local_d_dos = std::min(local_dos_occ, dos_en_step*double(electron_nearest_electron_list[e][0])/(E_f_A-DoS_cutoff));
    if( e_energy > transport_cutoff) {
      e_occupation   = std::min(1.0, double(global_e_dos[e_index  ][0]) / (dos_standard[e_index]*dos_en_step)); 
      f_e_occupation = std::min(1.0, double(global_e_dos[f_e_index][0]) / (dos_standard[f_e_index]*dos_en_step)); 
      r_e_occupation = std::min(1.0, double(global_e_dos[r_e_index][0]) / (dos_standard[r_e_index]*dos_en_step)); 
      // e_occupation   =  std::min(1.0,double(ee_dos_hist[e][e_index  ]) / local_d_dos); 
      // f_e_occupation =  std::min(1.0,double(ee_dos_hist[e][e_index+1]) / local_d_dos); 
      // r_e_occupation =  std::min(1.0,double(ee_dos_hist[e][e_index-1]) / local_d_dos); 
    } else {
     e_occupation   = std::min(1.0, double(global_e_dos[e_index  ][0]) / std::max(1.0, double(global_e_dos[e_index  ][1])));    
     f_e_occupation = std::min(1.0, double(global_e_dos[f_e_index][0]) / std::max(1.0, double(global_e_dos[f_e_index][1]))); 
     r_e_occupation = std::min(1.0, double(global_e_dos[r_e_index][0]) / std::max(1.0, double(global_e_dos[r_e_index][1]))); 
    }
    
    const double thermal_factor = return_BE_distribution(phonon_factor, Tp);
    const double f_factor = thermal_factor*(1.0 - f_e_occupation);//*e_occupation;// - f_e_occupation*(1.0-e_occupation);
    const double r_factor = (thermal_factor + 1.0) *(1.0-r_e_occupation);//*e_occupation;
    const double rta_p = ea_rate*M_ep*pho_s*f_factor;
    const double rta_m = ea_rate*M_ep*pho_s*r_factor;
    global_tau_ep[2*e_index] += rta_p;
    global_tau_ep[2*e_index+1] -= rta_m;

    // #pragma omp critical 
    // std::cout << exp(rta_p) << ", " << exp(rta_m) << std::endl;

    if( (f_factor > 0.0) &&  ( omp_uniform_random[thread]() > std::exp(rta_p)) ) {
//    #pragma omp critical 
//    std::cout << std::exp(rta_p) << ", " << ea_rate << ", " << phonon_factor/(q_sq*q_sq + p*p) << ", " << pho_s << ", " << f_factor << ", " << r_factor << std::endl;
      if( p_energy > (E_f_A+ max_as) ) return;

        if(p_k + p <= k ) std::cout << "weird stuff " << p_k << " + " << p << " < " << k << std::endl;
    if(m_k + p <= k) std::cout << "normal stuff " << m_k << " + " << p << " < " << k << std::endl;
    if((k + std::min(p_k, p)) <= std::max(p_k, p)) std::cout << "normal stuff " << (k + std::min(p_k, p)) << " < " <<  std::max(p_k, p) << "; " << k << ", " << p_k << ", " << p << std::endl;
    if((k + std::min(m_k, p)) <= std::max(m_k, p)) std::cout << "weird? stuff " << (k + std::min(m_k, p)) << " < " <<  std::max(m_k, p) << "; " << k << ", " << m_k << ", " << p << std::endl;

    const double n[3] = {electron_velocity[array_index], electron_velocity[array_index+1], electron_velocity[array_index+2]};
    double O[3] = {1.0/sqrt(3.0),1.0/sqrt(3.0),1.0/sqrt(3.0)};
    if(n[0]*O[0] + n[1]*O[1] + n[2]*O[2] > 0.99) {
     //   std::cout << "change basis: " << n[0]*O[0] << ", " << n[1]*O[1] << ", " << n[2]*O[2] << std::endl;
        O[0] = 1.0;
        O[1] = 0.0;
        O[2] = 0.0;
    }

    double t[3] = {
                        O[1]*n[2] - O[2]*n[1],n[0]*O[2]-O[0]*n[2],n[1]*O[0]-n[0]*O[1]
    };
    //  #pragma omp critical 
    // {
    //     std::cout << "t_i: <" << O[1]*n[2] - O[2]*n[1] << ", " << n[0]*O[2]-O[0]*n[2] << ", " << n[1]*O[0]-n[0]*O[1] << std::endl; 
    //   // std::cout << "n_i: <" << n[0] << ", " << n[1] << ", " << n[2] << ">, t_i:  <" << t[0] << ", "<< t[1] << ", " << t[2] << std::endl; 
    // } 
    const double t_n = 1.0/sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
    t[0] *= t_n;
    t[1] *= t_n;
    t[2] *= t_n;
    // #pragma omp critical 
    // {
    //    std::cout << "n_i: <" << n[0] << ", " << n[1] << ", " << n[2] << ">, t_i:  <" << t[0] << ", "<< t[1] << ", " << t[2] << std::endl; 
    // } 
   
    double b[3] = {
                        t[1]*n[2] - t[2]*n[1],n[0]*t[2]-t[0]*n[2],n[1]*t[0]-n[0]*t[1]
    };
    
    double d_p[3] = {
                        k*n[0]*h_p, k*n[1]*h_p, k*n[2]*h_p
    };
    double theta_p = omp_uniform_random[thread]()*M_PI*2.0;

    d_p[0] += r_p*(t[0]*cos(theta_p)+b[0]*sin(theta_p));
    d_p[1] += r_p*(t[1]*cos(theta_p)+b[1]*sin(theta_p));
    d_p[2] += r_p*(t[2]*cos(theta_p)+b[2]*sin(theta_p));
      
      int e_index = int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25)))) ;
      relaxation_time_hist_ee[array_index][e_index]++;
      relaxation_time_hist_ee[array_index + 2][e_index] += current_time_step - relaxation_time_hist_ee[array_index + 1][e_index ];
      
      electron_potential[e] = p_energy;
      
      relaxation_time_hist_ee[array_index + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((electron_potential[e]-DoS_cutoff)/0.25)))) ] = current_time_step;
      
      electron_velocity[array_index]   = d_p[0]/p_k;
      electron_velocity[array_index+1] = d_p[1]/p_k;
      electron_velocity[array_index+2] = d_p[2]/p_k;
    
      electron_ee_scattering_list[e][0] = 1;

      #pragma omp critical(eascattering)
      {
     //  std::cout << p_k << " vs " << sqrt(d_p[0]*d_p[0] + d_p[1]*d_p[1] + d_p[2]*d_p[2]) << ", " << p_energy << ", " <<  return_dWdE_i(p_k) << std::endl;
      ea_transport_scattering_count++;
      TEKE += p_energy-e_energy;
      TLE -= p_energy-e_energy;
      e_a_scattering_count++;
      }
      return;
    }  

    if(r_factor > 0.0 && omp_uniform_random[thread]() > exp(rta_m)) {
      
//       #pragma omp critical 
//    std::cout << ea_rate << ", " << M_ep << ", " << pho_s << ", " << p_mom_ratio/(1.5) << ", " << m_mom_ratio/(1.5) << ", " << f_factor << ", " << r_factor << std::endl;
   
      if( m_energy < core_cutoff) return;

    if(p_k + p <= k ) std::cout << "weird stuff " << p_k << " + " << p << " < " << k << std::endl;
    if(m_k + p <= k) std::cout << "normal stuff " << m_k << " + " << p << " < " << k << std::endl;
    if((k + std::min(p_k, p)) <= std::max(p_k, p)) std::cout << "normal stuff " << (k + std::min(p_k, p)) << " < " <<  std::max(p_k, p) << "; " << k << ", " << p_k << ", " << p << std::endl;
    if((k + std::min(m_k, p)) <= std::max(m_k, p)) std::cout << "weird? stuff " << (k + std::min(m_k, p)) << " < " <<  std::max(m_k, p) << "; " << k << ", " << m_k << ", " << p << std::endl;
    
    const double n[3] = {electron_velocity[array_index], electron_velocity[array_index+1], electron_velocity[array_index+2]};
    double O[3] = {1.0/sqrt(3.0),1.0/sqrt(3.0),1.0/sqrt(3.0)};
    if(n[0]*O[0] + n[1]*O[1] + n[2]*O[2] > 0.99) {
     //   std::cout << "change basis: " << n[0]*O[0] << ", " << n[1]*O[1] << ", " << n[2]*O[2] << std::endl;
        O[0] = 1.0;
        O[1] = 0.0;
        O[2] = 0.0;
    }

    double t[3] = {
                        O[1]*n[2] - O[2]*n[1],n[0]*O[2]-O[0]*n[2],n[1]*O[0]-n[0]*O[1]
    };
   
    const double t_n = 1.0/sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
    t[0] *= t_n;
    t[1] *= t_n;
    t[2] *= t_n;
   
    double b[3] = {
                        t[1]*n[2] - t[2]*n[1],n[0]*t[2]-t[0]*n[2],n[1]*t[0]-n[0]*t[1]
    };
    
    double d_p[3] = {
                        k*n[0]*h_m, k*n[1]*h_m, k*n[2]*h_m
    };
    double theta_p = omp_uniform_random[thread]()*M_PI*2.0;

    d_p[0] += r_m*(t[0]*cos(theta_p)+b[0]*sin(theta_p));
    d_p[1] += r_m*(t[1]*cos(theta_p)+b[1]*sin(theta_p));
    d_p[2] += r_m*(t[2]*cos(theta_p)+b[2]*sin(theta_p));
      
      int e_index = int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25)))) ;
      relaxation_time_hist_ee[array_index][e_index]++;
      relaxation_time_hist_ee[array_index + 2][e_index] += current_time_step - relaxation_time_hist_ee[array_index + 1][e_index ];
      
      electron_potential[e] = m_energy;
      
      relaxation_time_hist_ee[array_index + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((electron_potential[e]-DoS_cutoff)/0.25)))) ] = current_time_step;
      
      electron_velocity[array_index]   = d_p[0]/m_k;
      electron_velocity[array_index+1] = d_p[1]/m_k;
      electron_velocity[array_index+2] = d_p[2]/m_k;
    
      electron_ee_scattering_list[e][0] = 1;

      #pragma omp critical(eascattering)
      {
     //  std::cout << p_k << " vs " << sqrt(d_p[0]*d_p[0] + d_p[1]*d_p[1] + d_p[2]*d_p[2]) << ", " << p_energy << ", " <<  return_dWdE_i(p_k) << std::endl;
      ea_core_scattering_count++;
      TEKE += m_energy-e_energy;
      TLE -= m_energy-e_energy;
      e_a_scattering_count++;
      }
      return;
    }
}


void ea_scattering(const int e, const int array_index, const int thread) {

    const double e_energy = electron_potential[e];
   // if(e_energy < 0.0) std::cout << "bigger problem " << std::endl;
    // double deltaE = 1.0; //double(e_index+1)+DoS_cutoff - e_energy;
   if(e_energy < core_cutoff+1.5) return;
    double e_occupation;
    double f_e_occupation;
    double r_e_occupation;
    const double phonon_factor = phonon_energy*omp_uniform_random[thread](); + (1.0*phonon_energy*mtrandom::gaussianc(omp_uniform_random[thread])/8.0);// ;//+ 0.16*abs(mtrandom::gaussianc(omp_uniform_random[thread]));// ;//
    double p = phonon_factor*1e5/(constants::hbar_r*3650.0); // 1e5/ A
   //if(phonon_factor < 0.0) std::cout << "problem " << phonon_factor << std::endl;
    const double theta = 2.0*M_PI*omp_uniform_random[thread]();
    const double phi   = M_PI*omp_uniform_random[thread]();
    double p_x = p*cos(theta)*sin(phi);
    double p_y = p*sin(theta)*sin(phi);
    double p_z = p*cos(phi);

    const double k = return_dWdE(e_energy);
    // std::cout << k << ", " << p << std::endl;

    double d_k_x_p = k*electron_velocity[array_index]+p_x;
    double d_k_y_p = k*electron_velocity[array_index + 1] + p_y;
    double d_k_z_p = k*electron_velocity[array_index+2] + p_z;
    double k_p = sqrt(d_k_x_p*d_k_x_p + d_k_y_p*d_k_y_p + d_k_z_p*d_k_z_p);
    double p_energy = return_dWdE_i(k_p);

    
    double d_k_x_m = k*electron_velocity[array_index]-p_x;
    double d_k_y_m = k*electron_velocity[array_index + 1] - p_y;
    double d_k_z_m = k*electron_velocity[array_index+2] - p_z;
    double k_m = sqrt(d_k_x_m*d_k_x_m + d_k_y_m*d_k_y_m + d_k_z_m*d_k_z_m);
    double m_energy = return_dWdE_i(k_m);
    
    // if(phonon_factor == (p_energy-e_energy) || phonon_factor == (e_energy-m_energy)){
    // #pragma omp critical
    // {
    // std::cout << e_energy << ", " << p_energy << ", " << e_energy-phonon_factor << ", " << m_energy << std::endl;
    
    // }
    // }
    const int e_index   = int(std::min( dos_size-2.0, std::max(1.0, floor((e_energy - DoS_cutoff)*i_dos_en_step))));
    const int f_e_index = int(std::min( dos_size-2.0, std::max(1.0, floor((p_energy - DoS_cutoff)*i_dos_en_step))));
    const int r_e_index = int(std::min( dos_size-2.0, std::max(1.0, floor((m_energy - DoS_cutoff)*i_dos_en_step))));
    //  std::cout << e_index << ", " << f_e_index << ", " << r_e_index << std::endl;
    // const double local_d_dos = std::min(local_dos_occ, dos_en_step*double(electron_nearest_electron_list[e][0])/(E_f_A-DoS_cutoff));
    if( e_energy > transport_cutoff) {
      e_occupation   = std::min(1.0, double(global_e_dos[e_index  ][0]) / (dos_standard[e_index]*dos_en_step)); 
      f_e_occupation = std::min(1.0, double(global_e_dos[f_e_index][0]) / (dos_standard[f_e_index]*dos_en_step)); 
      r_e_occupation = std::min(1.0, double(global_e_dos[r_e_index][0]) / (dos_standard[r_e_index]*dos_en_step)); 
      // e_occupation   =  std::min(1.0,double(ee_dos_hist[e][e_index  ]) / local_d_dos); 
      // f_e_occupation =  std::min(1.0,double(ee_dos_hist[e][e_index+1]) / local_d_dos); 
      // r_e_occupation =  std::min(1.0,double(ee_dos_hist[e][e_index-1]) / local_d_dos); 
    } else {
     e_occupation   = std::min(1.0, double(global_e_dos[e_index  ][0]) / std::max(1.0, double(global_e_dos[e_index  ][1])));    
     f_e_occupation = std::min(1.0, double(global_e_dos[f_e_index][0]) / std::max(1.0, double(global_e_dos[f_e_index][1]))); 
     r_e_occupation = std::min(1.0, double(global_e_dos[r_e_index][0]) / std::max(1.0, double(global_e_dos[r_e_index][1]))); 
    }
    
    const double thermal_factor = return_BE_distribution(phonon_factor, Tp);
    const double f_factor = thermal_factor*(1.0 - f_e_occupation);// - f_e_occupation*(1.0-e_occupation);
    const double r_factor = (thermal_factor + 1.0) *(1.0-r_e_occupation);
    global_tau_ep[2*e_index] += ea_rate*(f_factor)/(p*p+q_sq*q_sq);
    global_tau_ep[2*e_index+1] -=  ea_rate*(r_factor)/(p*p+q_sq*q_sq);
    if(f_factor > 0.0 && omp_uniform_random[thread]() > exp(ea_rate*f_factor/(p*p+q_sq*q_sq))) {
      
      
      if( p_energy > (E_f_A+ max_as) ) return;
      // if(omp_uniform_random[thread]()*factor < r_factor) deltaE *= -1.0;
      //  std::cout << e_occupation << ", " << d_e_occupation << ", " << double(e_index+1)+DoS_cutoff - e_energy << ", " << thermal_factor << ", "<< return_BE_integrand(abs(deltaE),Tp)*(e_occupation-d_e_occupation) << ", " << d_e_occupation*(1.0-e_occupation) << ", " << exp(ea_rate*abs(thermal_factor)) << std::endl;  
      relaxation_time_hist_ee[array_index].at(int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25)))) )++;
      relaxation_time_hist_ee[array_index + 2].at(int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25)))) ) += current_time_step - relaxation_time_hist_ee[array_index + 1].at(int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25)))) );
      
      electron_potential[e] = p_energy;
      
      relaxation_time_hist_ee[array_index + 1].at(int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((electron_potential[e]-DoS_cutoff)/0.25)))) ) = current_time_step;
      
      // electron_transport_list[e] = true;
      // if(electron_potential[e] < transport_cutoff) electron_transport_list[e] = false;
      // const double theta = 2.0*M_PI*omp_uniform_random[thread]();
      // const double phi   = M_PI*omp_uniform_random[thread]();
     // if(electron_velocity[array_index+2] < 0.0) theta += M_PI;

      // const double scattering_velocity = return_vel(electron_potential[e]);
      electron_velocity[array_index]   = d_k_x_p /k_p;//cos(theta)*sin(phi);
      electron_velocity[array_index+1] = d_k_y_p/k_p;//sin(theta)*sin(phi);
      electron_velocity[array_index+2] = d_k_z_p/k_p;//cos(phi); 
    
      electron_ee_scattering_list[e][0] = 1;

      #pragma omp critical(eascattering)
      {
      //lattice_output << phonon_factor*i_dos_en_step/2.0 << ", " << f_factor << ", " << r_factor << ", " << thermal_factor << ", " << e_occupation << ", " <<\
      f_e_occupation << ", " << r_e_occupation << ", " << \
      thermal_factor*(return_fermi_distribution((e_index*dos_en_step+core_cutoff-E_f_A), Te)-return_fermi_distribution((e_index*dos_en_step+core_cutoff+phonon_factor-E_f_A), Te))\
    - return_fermi_distribution((e_index*dos_en_step+core_cutoff-E_f_A), Te)*(1.0 - return_fermi_distribution((e_index*dos_en_step+core_cutoff-E_f_A), Te))  << ", " << \
      thermal_factor*(return_fermi_distribution((e_index*dos_en_step+core_cutoff-E_f_A), Te)-return_fermi_distribution((e_index*dos_en_step+core_cutoff-phonon_factor-E_f_A), Te))\
    + return_fermi_distribution((e_index*dos_en_step+core_cutoff-E_f_A), Te)*(1.0 - return_fermi_distribution((e_energy-phonon_factor-E_f_A), Te)) << ", " << \
      return_fermi_distribution((e_index*dos_en_step+core_cutoff-E_f_A), Te) << ", " << \
      return_fermi_distribution((e_index*dos_en_step+core_cutoff+phonon_factor-E_f_A), Te) << ", " <<\
      return_fermi_distribution((e_index*dos_en_step+core_cutoff-phonon_factor-E_f_A), Te) << std::endl;
      //f(deltaE < 0.0) ea_core_scattering_count++;
      ea_transport_scattering_count++;
      TEKE += p_energy-e_energy;
      TLE -= p_energy-e_energy;
      e_a_scattering_count++;
      }
      return;
    }  
    if(r_factor > 0.0 && omp_uniform_random[thread]() > exp(ea_rate*r_factor/(p*p+q_sq*q_sq))) {
     // double deltaE = phonon_factor;
      if(m_energy   < core_cutoff ) return;
     // if(omp_uniform_random[thread]()*factor < r_factor) deltaE *= -1.0;
      //  std::cout << e_occupation << ", " << d_e_occupation << ", " << double(e_index+1)+core_cutoff - e_energy << ", " << thermal_factor << ", "<< return_BE_integrand(abs(deltaE),Tp)*(e_occupation-d_e_occupation) << ", " << d_e_occupation*(1.0-e_occupation) << ", " << exp(ea_rate*abs(thermal_factor)) << std::endl;  
      relaxation_time_hist_ee[array_index][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))))]++;
      relaxation_time_hist_ee[array_index + 2][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))))] += current_time_step - relaxation_time_hist_ee[array_index + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))))];
      
      electron_potential[e] = m_energy;
      
      relaxation_time_hist_ee[array_index + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((electron_potential[e]-DoS_cutoff)/0.25))))] = current_time_step;
      
      // electron_transport_list[e] = true;
      // if(electron_potential[e] < transport_cutoff) electron_transport_list[e] = false;
    //   const double theta = 2.0*M_PI*omp_uniform_random[thread]();
    //   const double phi   = M_PI*omp_uniform_random[thread]();
     // if(electron_velocity[array_index+2] < 0.0) theta += M_PI;

      // const double scattering_velocity = return_vel(electron_potential[e]);
      electron_velocity[array_index]   = d_k_x_m /k_m;//cos(theta)*sin(phi);
      electron_velocity[array_index+1] = d_k_y_m /k_m;//sin(theta)*sin(phi);
      electron_velocity[array_index+2] = d_k_z_m /k_m;//cos(phi); 
    
      electron_ee_scattering_list[e][0] = 1;

      #pragma omp critical(eascattering)
      {
     //lattice_output << phonon_factor*i_dos_en_step/2.0 << ", " << f_factor << ", " << r_factor << ", " << thermal_factor << ", " << e_occupation << ", " <<\
      f_e_occupation << ", " << r_e_occupation << ", " << \
      thermal_factor*(return_fermi_distribution((e_index*dos_en_step+core_cutoff-E_f_A), Te)-return_fermi_distribution((e_index*dos_en_step+core_cutoff+phonon_factor-E_f_A), Te))\
    - return_fermi_distribution((e_index*dos_en_step+core_cutoff-E_f_A), Te)*(1.0 - return_fermi_distribution((e_index*dos_en_step+core_cutoff-E_f_A), Te))  << ", " << \
      thermal_factor*(return_fermi_distribution((e_index*dos_en_step+core_cutoff-E_f_A), Te)-return_fermi_distribution((e_index*dos_en_step+core_cutoff-phonon_factor-E_f_A), Te))\
    + return_fermi_distribution((e_index*dos_en_step+core_cutoff-E_f_A), Te)*(1.0 - return_fermi_distribution((e_energy-phonon_factor-E_f_A), Te)) << ", " << \
      return_fermi_distribution((e_index*dos_en_step+core_cutoff-E_f_A), Te) << ", " << \
      return_fermi_distribution((e_index*dos_en_step+core_cutoff+phonon_factor-E_f_A), Te) << ", " <<\
      return_fermi_distribution((e_index*dos_en_step+core_cutoff-phonon_factor-E_f_A), Te) << std::endl;
      ea_core_scattering_count++;
     // else ea_transport_scattering_count++;
      TEKE += m_energy-e_energy;
      TLE -= m_energy-e_energy;
      e_a_scattering_count++;
    //   std::cout << std::scientific;
     // std::cout << TLE << ", " << m_energy-e_energy << ", " << k_m << ", " << p  << ", " << phonon_factor << std::endl;
     // std::cout << std::fixed;
      }
      return;
    }  
}

void ee_scattering() {
         if (err::check) std::cout << "ee_scattering." << std::endl;

  omp_set_dynamic(0);
  omp_set_num_threads(omp_threads);
 // const static double q_sq = 1.4*1.4;//*constants::hbar_over_me_sqrt*constants::hbar_over_me_sqrt;
  #pragma omp parallel 
  {
    
  const int thread = omp_get_thread_num();

  for(int l = 0; l < cells_per_thread; l++) {
    const   int cell = lattice_cells_per_omp[thread][l];
    const   int size = cell_integration_lists[cell][0];
    
    #pragma omp barrier 

    for(int e = 1; e < size; e++) {
      
      const int electron = cell_integration_lists[cell][e];
      if (electron_ee_scattering_list[electron][0] == 1) continue;
      const int scattering_size = electron_ee_scattering_list[electron][1];
  
      for(int a = 1; a < scattering_size; a++) {

         int electron_collision = electron_ee_scattering_list[electron][a*2];
         if (electron_ee_scattering_list[electron_collision][0] == 1) continue;
       
         const double e_energy = electron_potential[electron];
         const int array_index = 3*electron;
         const double d_e_energy = electron_potential[electron_collision];
        
         const int array_index_i = 3*electron_collision;
         // bool e_ballistic = (e_energy > (E_f_A+4.8)) ? true : false;
         // bool d_e_ballistic = (d_e_energy > (E_f_A+4.8)) ? true : false;

         
         // if(!
         elastic_scattering(thread, electron, array_index, electron_collision, array_index_i, e_energy, d_e_energy);// && (e_ballistic != d_e_ballistic) && (e_ballistic || d_e_ballistic) ) \
         //   inelastic_scattering(thread, electron, array_index, electron_collision, array_index_i, e_energy, d_e_energy);

      }
    } 
  }   
  } if (err::check) std::cout << "ee_scattering done." << std::endl;
  
}

bool elastic_scattering(int thread, int e, int array_index, int d_e, int array_index_i, double e_energy, double d_e_energy ) {
   
   double distance[3] = {electron_position[array_index]  - electron_position[array_index_i], \
                         electron_position[array_index+1]- electron_position[array_index_i+1],\
                         electron_position[array_index+2]- electron_position[array_index_i+2]};

   if (distance[0] < -1.0*boundary_conditions_cutoff)  distance[0]  += lattice_width;
   else if (distance[0]  > boundary_conditions_cutoff) distance[0]  -= lattice_width;
   
   if (distance[1]  < -1.0*boundary_conditions_cutoff) distance[1] += lattice_depth;
   else if (distance[1] > boundary_conditions_cutoff)  distance[1] -= lattice_depth;    
   
   if (distance[2] <  -1.0*boundary_conditions_cutoff) distance[2] += lattice_height;
   else if (distance[2] > boundary_conditions_cutoff)  distance[2] -= lattice_height;

   double length  = sqrt(distance[0] *distance[0]  + distance[1]*distance[1] + distance[2]*distance[2]);
   if(length == 0.0) return false;
        
  double l_theta = atan2(distance[1],distance[0] );
    if(l_theta!=l_theta) l_theta = 0.0;
  double l_phi = acos(distance[2]/length);

  distance[0]  = 4.16*cos(l_theta)*sin(l_phi);
  distance[1] = 4.16*sin(l_theta)*sin(l_phi);
  distance[2] = 4.16*cos(l_phi);
  length = 4.16*4.16;
  
   const double k1 = return_dWdE(e_energy);
   const double k2 = return_dWdE(d_e_energy);

    double k_1[3] = {electron_velocity[array_index  ]*k1,\
                          electron_velocity[array_index+1]*k1,\
                          electron_velocity[array_index+2]*k1};

    double k_2[3] = {electron_velocity[array_index_i  ]*k2,\
                          electron_velocity[array_index_i+1]*k2,\
                          electron_velocity[array_index_i+2]*k2};

   const double v_x_dot_product =(( k_1[0] - k_2[0])*distance[0] \
                                 + (k_1[1] - k_2[1])*distance[1]\
                                 + (k_1[2] - k_2[2])*distance[2]) /length;

          if(v_x_dot_product != v_x_dot_product) {std::cout << "ee error: " << v_x_dot_product << std::endl; return false;}
        
   double d_k_1[3] = {k_1[0] - distance[0] *v_x_dot_product,\
                      k_1[1] - distance[1]*v_x_dot_product,\
                       k_1[2] - distance[2]*v_x_dot_product};

   double d_k_2[3] = {k_2[0] + distance[0] *v_x_dot_product,\
                      k_2[1] + distance[1]*v_x_dot_product,\
                      k_2[2] + distance[2]*v_x_dot_product};

   double d_k1 = sqrt(d_k_1[0]*d_k_1[0] + d_k_1[1]*d_k_1[1] + d_k_1[2]*d_k_1[2]);
   double d_k2 = sqrt(d_k_2[0]*d_k_2[0] + d_k_2[1]*d_k_2[1] + d_k_2[2]*d_k_2[2]);
      
   double d_e_1 = return_dWdE_i(d_k1);
   double d_e_2 = return_dWdE_i(d_k2);

   double deltaE = e_energy+d_e_energy - d_e_1 - d_e_2;
   const double deltaK = (k1 - d_k1)*(k1 - d_k1);
   // if(d_e_1 > (E_f_A+4.8) || d_e_2 > (E_f_A+4.8)) std::cout << d_e_1 << ", " << d_k_1 << ", " << d_e_2 << ", "  << d_k_2 << ", " << deltaE << ", " << deltaK << std::endl;
   // double a_factor = abs((k_1*k_1 + k_2*k_2 - d_k_1*d_k_1 - d_k_2*d_k_2)/(2.0*d_k_1*d_k_2));
   // double b_factor = k_1*k_2/(d_k_1*d_k_2);

      if(d_e_1 < core_cutoff || d_e_2 < core_cutoff) return false;
      if( d_e_1 > (E_f_A+ max_as) || d_e_2 > (E_f_A+ max_as)) return false;
      if(abs(deltaE) > 0.01) return false; //std::cout << "ee collision deltaE: " << deltaE << ", " << e_energy << ", " << d_e_energy << ", " << d_e_1 << ", " << d_e_2 << std::endl;
     // if(a_factor > (b_factor+1)) {
      //#//pragma omp critical 
      //    std::cout << " a/b factor " << a_factor << ", " <<  b_factor << ", " << k_1  +k_2 << ", " << d_k_1 + d_k_2 << std::endl;
        //  return false;}
       
      double d_e_occupation;  
      double d_d_occupation;
      double e_occupation;  
      double d_occupation;

      const int e_index = int(std::min( dos_size-1.0, std::max(0.0, floor((e_energy - DoS_cutoff)*i_dos_en_step))));
      if ( e_energy > transport_cutoff) e_occupation = std::min(1.0, double(global_e_dos[e_index][0])/(dos_standard[e_index]*dos_en_step));  
      else  e_occupation = std::min(1.0, double(global_e_dos[e_index][0]) / double(global_e_dos[e_index][1]));  

      const int d_index = int(std::min( dos_size-1.0, std::max(0.0, floor((d_e_energy - DoS_cutoff)*i_dos_en_step))));
      if ( d_e_energy > transport_cutoff)  d_occupation = std::min(1.0, double(global_e_dos[d_index][0])/(dos_standard[d_index]*dos_en_step));  
      else d_occupation = std::min(1.0, double(global_e_dos[d_index][0]) / double(global_e_dos[d_index][1])); 
     
      const int d_e_index = int(std::min( dos_size-1.0, std::max(0.0, floor((d_e_1 - DoS_cutoff)*i_dos_en_step))));
      if ( d_e_1 > transport_cutoff) d_e_occupation = std::max(0.0, 1.0 - double(global_e_dos[d_e_index][0])/(dos_standard[d_e_index]*dos_en_step));  
      else  d_e_occupation = std::max(0.0, 1.0 - double(global_e_dos[d_e_index][0]) / double(global_e_dos[d_e_index][1]));  

      const int d_d_index = int(std::min( dos_size-1.0, std::max(0.0, floor((d_e_2 - DoS_cutoff)*i_dos_en_step))));
      if ( d_e_2 > transport_cutoff)  d_d_occupation = std::max(0.0, 1.0 - double(global_e_dos[d_d_index][0])/(dos_standard[d_d_index]*dos_en_step));  
      else d_d_occupation = std::max(0.0, 1.0  - double(global_e_dos[d_d_index][0]) / double(global_e_dos[d_d_index][1])); 

      double occupation_factor = ee_rate*d_d_occupation*d_e_occupation/((q_sq+deltaK)*(q_sq+deltaK));

      // if(e_energy > E_f_A+4.8 || d_e_energy > E_f_A+4.8) occupation_factor /= q_sq*q_sq;
      // else
      // occupation_factor /= (q_sq+deltaK)*(q_sq+deltaK);

      global_tau_ee[d_e_index] += occupation_factor;
      global_tau_ee[d_d_index] += occupation_factor;

   if(omp_uniform_random[thread]() > exp(occupation_factor)) {
          
          int e_relax_index = int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))));
          int d_relax_index = int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_energy-DoS_cutoff)/0.25))));
          // if(a_factor > (b_factor+1.0) || abs(deltaE) < 0.1) std::cout << a_factor << " > " << b_factor+1 << ", dE: "<< deltaE << ", " << k_1 << ", " << k_2 << ", " << d_k_1 << ", " << d_k_2 << ", " << e_energy << ", " << d_e_energy << ", " << d_e_1 << ", " << d_e_2 << std::endl;
          relaxation_time_hist_ee[array_index][e_relax_index]++;
          relaxation_time_hist_ee[array_index_i][d_relax_index]++;
        
          relaxation_time_hist_ee[array_index + 2][e_relax_index] += current_time_step - relaxation_time_hist_ee[array_index + 1][e_relax_index];
          relaxation_time_hist_ee[array_index_i + 2][d_relax_index] += current_time_step - relaxation_time_hist_ee[array_index_i + 1][d_relax_index];
        
          electron_potential[e] = d_e_1;
          electron_potential[d_e] = d_e_2;

          relaxation_time_hist_ee[array_index + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_1-DoS_cutoff)/0.25))))] = current_time_step;
          relaxation_time_hist_ee[array_index_i + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_2-DoS_cutoff)/0.25))))] = current_time_step;

          electron_ee_scattering_list[e][0] = 1;
          electron_ee_scattering_list[d_e][0] = 1;

          double theta = atan2(d_k_1[1], d_k_1[0]);
          if(theta != theta) theta == 0.0;
          double phi = acos(d_k_1[2]/d_k1);
         
          electron_velocity[array_index]    = sin(phi)*cos(theta);
          electron_velocity[array_index+1]  = sin(phi)*sin(theta);
          electron_velocity[array_index+2]  = cos(phi); 

           theta = atan2(d_k_2[1], d_k_2[0]);
              if(theta != theta) theta == 0.0;
           phi = acos(d_k_2[2]/d_k2);

          electron_velocity[array_index_i]  = sin(phi)*cos(theta);
          electron_velocity[array_index_i+1]= sin(phi)*sin(theta);
          electron_velocity[array_index_i+2]= cos(phi); 

          #pragma omp critical(eescattering)
          {
          if(deltaE > 0.001) full_int_var++;
          else if (deltaE < -0.001) full_int_var--;
          
          if (electron_potential[e] < E_f_A + 4.8) ee_core_scattering_count++;
          else ee_transport_scattering_count++;
          if (electron_potential[d_e] < E_f_A + 4.8) ee_core_scattering_count++;
          else ee_transport_scattering_count++;
          // lattice_output << e_e_scattering_count << ", " << electron << ", " << electron_collision << ", " \
                         << x_distance*v_x_dot_product << ", " << y_distance*v_x_dot_product << ", " << z_distance*v_x_dot_product << ", " \
                         << deltaE << ", " << deltaK << ", "  << sqrt(length) << ", " << v_x_dot_product << ", " \
                         << ee_rate*e_occupation*d_e_occupation/((0.25+(deltaK))*(0.25+(deltaK))) << ", " << ee_rate*e_occupation*d_e_occupation/((0.25+(deltaE*deltaE))*(0.25+(deltaE*deltaE))) << ", " << rspace_Hsr << std::endl;
          e_e_scattering_count += 2;
          // ee_core_scattering_count += 2;
          //if(electron_potential[electron] > transport_cutoff || electron_potential[electron_collision] > transport_cutoff )lattice_output <<\
            electron << ", " << e_occupation << ", " << d_e_occupation << ", " << 1.0/((q_sq+(deltaK*deltaK))*(q_sq+(deltaK*deltaK))) << ", " << deltaK << ", " << deltaE << ", " << e_energy<< ", " << d_e_energy << ", " <<\
            1-return_fermi_distribution(e_index*dos_en_step+deltaE+core_cutoff-E_f_A,constants::kB_r*Te) << ", " << 1-return_fermi_distribution(d_index*dos_en_step+core_cutoff-E_f_A-deltaE,constants::kB_r*Te) << ", " <<\
            sqrt(2.0*constants::m_e_r_i*e_energy/constants::hbar_r/constants::hbar_r) << ", " << sqrt(2.0*constants::m_e_r_i*(e_energy+deltaE)/constants::hbar_r/constants::hbar_r) << ", " << constants::m_over_hbar_sqrt*sqrt(v_x*v_x + v_y*v_y + v_z*v_z) << ", " << constants::m_over_hbar_sqrt*sqrt(d_v_x_2*d_v_x_2 + d_v_y_2*d_v_y_2 + d_v_z_2*d_v_z_2) << std::endl;
          // if(e_back_pressure > 1.0) e_e_scattering_count += 2;
          }
      return true;
   } else return false;
}

void inelastic_scattering(int thread, int e, int array_index, int d_e, int array_index_i, double e_energy, double d_e_energy) {
   
   double deltaE = 0.1*(e_energy-d_e_energy);
   double d_e_occupation;  
   double d_d_occupation;
   double d_e_1 = e_energy-deltaE;
   double d_e_2 = d_e_energy+deltaE;
   
   if(d_e_1 < core_cutoff || d_e_2 < core_cutoff) return;
   if( d_e_1 > (E_f_A+ max_as) || d_e_2 > (E_f_A+ max_as)) return;

   const int e_index = int(std::min( dos_size-1.0, std::max(0.0, floor((d_e_1 - DoS_cutoff)*i_dos_en_step))));
   if ( d_e_1 > transport_cutoff) d_e_occupation = std::max(0.0, 1.0 - double(global_e_dos[e_index][0])/(dos_standard[e_index]*dos_en_step));  
   else  d_e_occupation = std::max(0.0, 1.0 - double(global_e_dos[e_index][0]) / double(global_e_dos[e_index][1]));  

   const int d_index = int(std::min( dos_size-1.0, std::max(0.0, floor((d_e_2 - DoS_cutoff)*i_dos_en_step))));
   if ( d_e_2 > transport_cutoff)  d_d_occupation = std::max(0.0, 1.0 - double(global_e_dos[d_index][0])/(dos_standard[d_index]*dos_en_step));  
   else d_d_occupation = std::max(0.0, 1.0  - double(global_e_dos[d_index][0]) / double(global_e_dos[d_index][1])); 
   double occupation_factor = sim::ee_coupling_strength* ee_rate*d_e_occupation*d_d_occupation/(q_sq*q_sq);///((q_sq+(deltaK))*(q_sq+(deltaK)));//*exp(0.15*(d_occupation*e_occupation-1.0));

      global_tau_ee[e_index] += occupation_factor;
      global_tau_ee[d_index] += occupation_factor;

   if(omp_uniform_random[thread]() > exp(occupation_factor)) {
          
      int e_relax_index = int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))));
          int d_relax_index = int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_energy-DoS_cutoff)/0.25))));
          // if(a_factor > (b_factor+1.0) || abs(deltaE) < 0.1) std::cout << a_factor << " > " << b_factor+1 << ", dE: "<< deltaE << ", " << k_1 << ", " << k_2 << ", " << d_k_1 << ", " << d_k_2 << ", " << e_energy << ", " << d_e_energy << ", " << d_e_1 << ", " << d_e_2 << std::endl;
          relaxation_time_hist_ee[array_index][e_relax_index]++;
          relaxation_time_hist_ee[array_index_i][d_relax_index]++;
        
          relaxation_time_hist_ee[array_index + 2][e_relax_index] += current_time_step - relaxation_time_hist_ee[array_index + 1][e_relax_index];
          relaxation_time_hist_ee[array_index_i + 2][d_relax_index] += current_time_step - relaxation_time_hist_ee[array_index_i + 1][d_relax_index];
        
          electron_potential[e] = d_e_1;
          electron_potential[d_e] = d_e_2;

          relaxation_time_hist_ee[array_index + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_1-DoS_cutoff)/0.25))))] = current_time_step;
          relaxation_time_hist_ee[array_index_i + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_2-DoS_cutoff)/0.25))))] = current_time_step;

      electron_ee_scattering_list[e][0] = 1;
      electron_ee_scattering_list[d_e][0] = 1;

      double theta = 2.0*M_PI*omp_uniform_random[thread]();
      double phi = M_PI*omp_uniform_random[thread]();

          electron_velocity[array_index]    = sin(phi)*cos(theta);
          electron_velocity[array_index+1]  = sin(phi)*sin(theta);
          electron_velocity[array_index+2]  = cos(phi); 

          electron_velocity[array_index_i]  = sin(-phi)*cos(-theta);
          electron_velocity[array_index_i+1]= sin(-phi)*sin(-theta);
          electron_velocity[array_index_i+2]= cos(-phi); 

          #pragma omp critical(eescattering)
          {
          if(deltaE > 0.001) full_int_var++;
          else if (deltaE < -0.001) full_int_var--;
         //  if (electron_potential[e] < transport_cutoff) ee_core_scattering_count++;
         //  else ee_transport_scattering_count++;
         //  if (electron_potential[d_e] < transport_cutoff) ee_core_scattering_count++;
         //  else ee_transport_scattering_count++;
          // lattice_output << e_e_scattering_count << ", " << electron << ", " << electron_collision << ", " \
                         << x_distance*v_x_dot_product << ", " << y_distance*v_x_dot_product << ", " << z_distance*v_x_dot_product << ", " \
                         << deltaE << ", " << deltaK << ", "  << sqrt(length) << ", " << v_x_dot_product << ", " \
                         << ee_rate*e_occupation*d_e_occupation/((0.25+(deltaK))*(0.25+(deltaK))) << ", " << ee_rate*e_occupation*d_e_occupation/((0.25+(deltaE*deltaE))*(0.25+(deltaE*deltaE))) << ", " << rspace_Hsr << std::endl;
          e_e_scattering_count += 2;
          ee_transport_scattering_count += 2;
          //if(electron_potential[electron] > transport_cutoff || electron_potential[electron_collision] > transport_cutoff )lattice_output <<\
            electron << ", " << e_occupation << ", " << d_e_occupation << ", " << 1.0/((q_sq+(deltaK*deltaK))*(q_sq+(deltaK*deltaK))) << ", " << deltaK << ", " << deltaE << ", " << e_energy<< ", " << d_e_energy << ", " <<\
            1-return_fermi_distribution(e_index*dos_en_step+deltaE+core_cutoff-E_f_A,constants::kB_r*Te) << ", " << 1-return_fermi_distribution(d_index*dos_en_step+core_cutoff-E_f_A-deltaE,constants::kB_r*Te) << ", " <<\
            sqrt(2.0*constants::m_e_r_i*e_energy/constants::hbar_r/constants::hbar_r) << ", " << sqrt(2.0*constants::m_e_r_i*(e_energy+deltaE)/constants::hbar_r/constants::hbar_r) << ", " << constants::m_over_hbar_sqrt*sqrt(v_x*v_x + v_y*v_y + v_z*v_z) << ", " << constants::m_over_hbar_sqrt*sqrt(d_v_x_2*d_v_x_2 + d_v_y_2*d_v_y_2 + d_v_z_2*d_v_z_2) << std::endl;
          // if(e_back_pressure > 1.0) e_e_scattering_count += 2;
          }
      return;
   }
}

double k_sq() {
  double q = 0.0;
  double f_0;
  // double f_2;
  double k;
 
  for(int h = int(floor((core_cutoff-DoS_cutoff)*i_dos_en_step)); h < dos_size-2; h++) {
  
    if( ((h)*dos_en_step+DoS_cutoff) > transport_cutoff) f_0 = double(global_e_dos[h][0])/(dos_standard[h]*dos_en_step);
    else f_0 = double(global_e_dos[h][0])/double(global_e_dos[h][1]);

    k = return_dWdE(h*dos_en_step + DoS_cutoff);
    q += dos_en_step*dos_standard[h]*f_0*1.0/(k*k*lattice_atoms*x_unit_size*y_unit_size*z_unit_size);
            // s^2 /  fs^2
  }

  return pow(sqrt(constants::K* q) + q_offset, 2.0);
}

} //end CASTLE namespace

