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

    int total_1 = 0;
    int total_2 = 0;
    for(int h = 0 ; h < dos_size; h++) {
       //std::cout << h << ", " << global_e_dos[h][0] << ", " << global_e_dos[h][1] << std::endl;
      total_1 += global_e_dos[h][0];
      total_2 += global_e_dos[h][1];
      global_e_dos[h][0] = 0;
    }
    // if(total_1 != conduction_electrons) std::cout << "hist problem " << total_1 << " != " << conduction_electrons << " != " << total_2 << std::endl;
    
    //reset integration
    CASTLE_real_time += dt;
    current_time_step++;
   
    return EXIT_SUCCESS;
}


void update_position() {

    omp_set_dynamic(0);
    omp_set_num_threads(omp_threads);
    
    if(current_time_step % half_int_var == 0) {
      cell_integration_lists.swap(old_cell_integration_lists);
     
      for(int c = 0; c < total_cells; c++) {
        cell_integration_lists[c][0] = 1;
      }
        if(err::check) std::cout << "cell lists swapped.";
    }
     if(err::check) std::cout << " updating position..." << std::endl;
    int total = 0;
  
  #pragma omp parallel reduction(+:x_flux,y_flux,z_flux, p_x,p_y,p_z)
  {
    std::vector<int> local_e_dos;
    local_e_dos.resize(dos_size,0);
    const int thread = omp_get_thread_num();
    //if(omp_get_num_threads() != omp_threads) std::cout << "omp pragma error " << omp_get_num_threads() << " != " << omp_threads << std::endl;
    int local = 0;
    for (int l = 0 ; l < cells_per_thread; l++) {
      const int cell = lattice_cells_per_omp[thread][l];
      int size;
      if(current_time_step % half_int_var == 0) size = old_cell_integration_lists[cell][0];
      else  size = cell_integration_lists[cell][0];
      std::vector<int> cell_integration_list;
      cell_integration_list = (current_time_step % half_int_var == 0) ? old_cell_integration_lists[cell] : cell_integration_lists[cell];
      local += size - 1;
  //  if((current_time_step % half_int_var) == 0 && ( (l % ((z_omp_cells/2)-1) )== 0 || (l % (z_omp_cells/2)) == 0) ) {
      #pragma omp barrier
  
      for (int e = 1; e < size; e++) { 
        int electron = cell_integration_list[e];
      //   if(current_time_step % half_int_var == 0) electron = old_cell_integration_lists[cell][e];
      //   else electron = cell_integration_lists[cell][e];
        const double energy = electron_potential[electron];
   
        const double vel = return_vel(energy);
        local_e_dos[int(std::min(dos_size-1.0, std::max(0.0, floor((energy - DoS_cutoff)*i_phonon_energy))))]++;
        const int array_index = 3*electron;

        double x_pos = electron_position[array_index];
        double y_pos = electron_position[array_index+1];
        double z_pos = electron_position[array_index+2];

        const double v_x = vel*electron_velocity[array_index];
        const double v_y = vel*electron_velocity[array_index+1];
        const double v_z = vel*electron_velocity[array_index+2];

        // if(electron_potential[e] > 0.98*E_f_A) {
          p_x += v_x;
          p_y += v_y;
          p_z += v_z;
          x_pos += v_x * dt;// + (electron_force[array_index]   * dt * dt * constants::K_A / 2); // x superarray component
          y_pos += v_y * dt;// + (electron_force[array_index_y) * dt * dt * constants::K_A / 2); // y superarray component
          z_pos += v_z * dt;// + (electron_force[array_index_z) * dt * dt * constants::K_A / 2); // z superarray component
        // }
          if(v_x != v_x || v_y != v_y || v_z != v_z) std::cout << energy << ", " << vel << ", " << v_x << ", " << v_y << ", " << v_z << std::endl;
        if (x_pos < 0.0) {
          x_pos += lattice_width; 
          x_flux--;
          if (current_time_step % CASTLE_output_rate == 0) {
            #pragma omp atomic
            flux_index[int(std::max(0.0, std::min(99.0, floor((energy - DoS_cutoff)/1.0))))]--;
          }
        } else if (x_pos > lattice_width) {
          x_pos -= lattice_width; 
          x_flux++;
          if (current_time_step % CASTLE_output_rate == 0) {
            #pragma omp atomic
            flux_index[int(std::max(0.0, std::min(99.0, floor((energy - DoS_cutoff)/1.0))))]++; 
          } 
        }

        if (y_pos < 0.0) {y_pos += lattice_depth; y_flux--;}
        else if (y_pos > lattice_depth) {y_pos -= lattice_depth; y_flux++;}

        if (z_pos < 0.0) {z_pos += lattice_height; z_flux--;}
        else if (z_pos > lattice_height) {z_pos -= lattice_height; z_flux++;}

        electron_position[array_index]   = x_pos;
        electron_position[array_index+1] = y_pos;
        electron_position[array_index+2] = z_pos;
    
      if(current_time_step % half_int_var == 0) {
          
           int x_cell = int(floor(x_pos / x_step_size));

           int y_cell = int(floor(y_pos / y_step_size));
    
           int z_cell = int(floor(z_pos / z_step_size));

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
         global_e_dos[h][0] += local_e_dos[h];
      }
   } 
   }
      if(err::check) std::cout << "global dos updated. catching escaped electrons: " << escaping_electrons[0] << std::endl;
  // if(total != conduction_electrons) std::cout << total << ", " << conduction_electrons << std::endl;
   if(current_time_step % half_int_var == 0) {    
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
  
  //  int array_index;
    double photon_energy = sim::applied_voltage*constants::eV_to_AJ; //AJ/hv
    //AJ/fs/nm**2 -> [1e-3volume(A**3)/(AJ/hv)) hv/fs 
    const static double photon_rate = 1e-2*power_density*lattice_width*lattice_depth/photon_energy; //hv/fs
    int photons_at_dt = 0; //hv*dt
    double pump = 0.0; // AJ / fs
    double external_potential = 0.0; //AJ/e-
    const static double sigma = dt * 0.1;
    int count = 0;

    if(!equilibrium_step) {
      if(heat_pulse_sim) {
        //hv(dt)/fs
        photons_at_dt = int(round(photon_rate*dt*exp(-0.5*sigma*sigma*((double(current_time_step) - ((40.0/dt)+sim::equilibration_time))*(double(current_time_step) - ((40.0 / dt)+sim::equilibration_time)))))); // AJ/fs/nm**3
        pump = 1e3*photons_at_dt*photon_energy/(dt*lattice_depth*lattice_height*lattice_width); //AJ/fs/nm**3
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
      
    omp_set_dynamic(0);
    omp_set_num_threads(omp_threads);
    std::vector<int> chosen;
    chosen.resize(photons_at_dt);
    while(count < photons_at_dt) {
      int select = int(omp_uniform_random[omp_get_thread_num()]()*2147483647) % (conduction_electrons);
      double e_energy = electron_potential[select];
      if(e_energy > (DoS_cutoff+2.0*phonon_energy) && return_fermi_distribution(e_energy+photon_energy-E_f_A, Te) < (1.0-1e-2) && (e_energy+photon_energy < core_cutoff+59.0)) {
        chosen[count] = select;
        count++;
      }
    }
    count = 0;
      //  pump = 0.0;
    Tp = d_Tp;
    Te = d_Te;
      if(err::check) std::cout << "conditions set; updating dos " << std::endl;


    #pragma omp parallel for schedule(dynamic, 4)
    for (int e = 0; e < conduction_electrons; e++) {
      // if(electron_potential[e] < (E_f_A-24.25)) continue;
      const int array_index = 3*e;
      
      if(current_time_step % half_int_var == 0) e_e_coulomb(e, array_index);
      // else if (current_time_step % half_int_ == 0 && electron_transport_list[e]) e_e_coulomb(e, array_index);
      else  neighbor_e_e_coulomb(e, array_index);

      if(photons_at_dt > 0 && std::end(chosen) != std::find(chosen.begin(), chosen.end(), e)) {
        #pragma omp atomic
        count++;
        // pump += external_potential;
        electron_thermal_field(e, array_index, photon_energy, omp_get_thread_num());
      }
      
      // if(!equilibrium_step) external_potential += electron_applied_voltage(e, array_index, photon_energy);
      // if(equilibrium_step && electron_potential[e] > 0.8*E_f_A) ea_scattering(e, array_index, omp_get_thread_num());
      ea_scattering(e, array_index, omp_get_thread_num());
    }

    if(count != photons_at_dt) std::cout << photons_at_dt << ", " << count<< std::endl;
      //  TEKE += external_potential;
        if(err::check) std::cout << "dos updated. ee scattering step..." << std::endl;
      q_sq = k_sq();
     
    ee_scattering();
   // pump /= 1e-3*lattice_depth*lattice_height*lattice_width;  
    if(!equilibrium_step) d_Tp =  a_heat_capacity_i*TLE *n_f*4.087312/conduction_electrons + Tp;
    else d_Tp = sim::temperature;

    d_Te = (e_heat_capacity_i*TEKE*n_f*4.087312/conduction_electrons/abs(Te)) + (dt*e_heat_capacity_i*pump/abs(Te)) + Te;
 
        if (err::check) std::cout << "reset scattering." << std::endl;
       
}

void electron_thermal_field(const int e, const int array_index, const double EKE, const int thread) {       
    
    const double theta = omp_uniform_random[thread]() * 2.0 * M_PI;
    const double phi   = omp_uniform_random[thread]() * M_PI; 

    electron_potential[e] += EKE;
  
    electron_velocity[array_index]   = cos(theta)*sin(phi);
    electron_velocity[array_index+1] = sin(theta)*sin(phi);
    electron_velocity[array_index+2] = cos(phi);

  //  electron_ee_scattering_list[e][0] = 1;
    electron_ea_scattering_list[e][0] = 1;
 
}

void e_e_coulomb(const int e, const int array_index) {

   const int cell = lattice_cell_coordinate[int(floor(electron_position[array_index] / x_step_size))][int(floor(electron_position[array_index+1] / y_step_size))][int(floor(electron_position[array_index+2] / z_step_size))];
   int ee_dos_count = 1;
   int ee_integration_count = 1;
   int ee_scattering_list = 1;
   double scattering_range = (electron_potential[e] > E_f_A+4.8) ? 400.0 : 49.0;
   double integration_range = (electron_potential[e] > E_f_A+4.8) ? 400.0 : 100.0;
   int e_size = 6*round(pow(scattering_range, 1.5)*1.25*M_PI * 3.8*n_f * 1e-3);
   int i_size = 3*round(pow(e_e_integration_cutoff,1.5)*1.25*M_PI * 3.8*n_f * 1e-3);
   if(electron_ee_scattering_list[e].size() < e_size) electron_ee_scattering_list[e].resize(e_size, 0);
   if(electron_integration_list[e].size() < i_size) electron_integration_list[e].resize(i_size, 0);
   int cells = (i_size > e_e_integration_cutoff) ? 125 : 27;
   std::vector<int> integration_list;
   integration_list = (cells > 27) ? cell_lr_neighbor_list[cell] : cell_nearest_neighbor_list[cell];

    for(int s = 0; s < cells; s++) {
      const int int_cell = integration_list[s];
      const int size = cell_integration_lists[int_cell][0];
      if(size != size || size < 1  ) std::cout << size << std::endl;
      for(int i = 1; i < size; i++) {
        const int electron = cell_integration_lists[int_cell][i];
        if(e == electron) continue;
        // if(electron_potential[electron] < (E_f_A-24.25)) continue;
        const int array_index_i = 3*electron;

        double x_distance = electron_position[array_index]   - electron_position[array_index_i];
        double y_distance = electron_position[array_index+1] - electron_position[array_index_i + 1];
        double z_distance = electron_position[array_index+2] - electron_position[array_index_i + 2]; 
    
        if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance += lattice_width;
        else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance -= lattice_width;

        if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance += lattice_depth;
        else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance -= lattice_depth;

        if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance += lattice_height;
        else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance -= lattice_height;
        
        const double length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
            
        if(length == 0.0) continue;
       
            if (length != length) std::cout << x_distance << ", " << y_distance << ", " << z_distance <<  ", " << \
        electron_position[array_index] << ", " << electron_position[array_index_i] << ", " << electron_position[array_index+1] << ", " << electron_position[array_index_i + 1] << ", " <<\
        electron_position[array_index+2] << ", " <<  electron_position[array_index_i + 2] << std::endl;
        
        if(length > e_e_integration_cutoff) continue;
        electron_integration_list[e][ee_integration_count] = array_index_i;
          if(ee_integration_count >= electron_integration_list[e].size() - 2) {std::cout << e << " ee integration " << ee_integration_count << " > " << electron_integration_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
          break; }
        ee_integration_count++;

      //   if(length > scattering_range) continue;
      //   electron_nearest_electron_list[e][ee_dos_count] = electron;
      //   ee_dos_count++;
      //   ee_dos_hist[e].at(int(std::max(0.0, std::min(dos_size-1.0, floor((electron_potential[electron]-DoS_cutoff)*i_phonon_energy)))))++;
      //   if(ee_dos_count >= electron_nearest_electron_list[e].size() - 2) {std::cout << e << " ee dos " << ee_dos_count << " > " << electron_nearest_electron_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
      //     break; }

        if (length > scattering_range) continue; 
        electron_ee_scattering_list[e][ee_scattering_list*2] = electron;
        electron_ee_scattering_list[e][ee_scattering_list*2 +1] = length;
        if(ee_scattering_list >= electron_ee_scattering_list[e].size() - 3) {std::cout << e << " ee scattering " << ee_scattering_list << " > " << electron_ee_scattering_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
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
   int s_size = 3*round(pow(e_e_integration_cutoff,1.5)*1.25*M_PI * 3.8*n_f * 1e-3);
   // int l_size = 3*round(pow(400.0,1.5)*1.25*M_PI * 3.8*n_f * 1e-3);

   double scattering_range = (electron_integration_list[e].size() > s_size) ? 400.0 : 49.0;
    for (int i = 1; i < size; i++) {
        const int array_index_i = electron_integration_list[e][i];
        const int electron = array_index_i/3;
        // if(electron_potential[electron] < (E_f_A-24.25)) continue;

        double x_distance = electron_position[array_index]   - electron_position[array_index_i];
        double y_distance = electron_position[array_index+1] - electron_position[array_index_i + 1];
        double z_distance = electron_position[array_index+2] - electron_position[array_index_i + 2]; 
    
        if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance = x_distance + lattice_width;
        else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance = x_distance - lattice_width;

        if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance = y_distance + lattice_depth;
        else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance = y_distance - lattice_depth;

        if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance = z_distance + lattice_height;
        else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance = z_distance - lattice_height;
        
        const double length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
        if(length == 0.0) continue;
        
       if (length != length) std::cout << x_distance << ", " << y_distance << ", " << z_distance <<  ", " << \
        electron_position[array_index] << ", " << electron_position[array_index_i] << ", " << electron_position[array_index+1] << ", " << electron_position[array_index_i + 1] << ", " <<\
        electron_position[array_index+2] << ", " <<  electron_position[array_index_i + 2] << std::endl;
        
      //   if(length > scattering_range) continue;
      //   electron_nearest_electron_list[e][neighbor_count] = array_index_i/3;
      //   neighbor_count++;
      //   ee_dos_hist[e][int(std::max(0.0, std::min(dos_size-1.0, floor((electron_potential[electron]-DoS_cutoff)*i_phonon_energy))))]++;
      //   if(neighbor_count >= electron_nearest_electron_list[e].size() - 2) {std::cout << e << ", " << neighbor_count << " > " << electron_nearest_electron_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
      //   break; }

        if (length > scattering_range) continue; 
        electron_ee_scattering_list[e][scattering_count*2] = array_index_i/3;
        electron_ee_scattering_list[e][scattering_count*2 +1] = length;
        scattering_count++;
    }

   //  electron_nearest_electron_list[e][0] = neighbor_count;
    electron_ee_scattering_list[e][1] = scattering_count;
    //  if(e == 0 ) std::cout << size << ", " << neighbor_count << ", " << scattering_count << std::endl;
}

double electron_applied_voltage(const int e, const int array_index, double& external_potential) {

  electron_velocity[array_index] += 5e-11*dt*applied_voltage*constants::e_A*constants::m_e_r_i; //V/m * q = F_x/kg = 0.5*a_x*dt = d_v_x
  const double deltaE = 0.5*constants::m_e_r*((electron_velocity[array_index]*electron_velocity[array_index])\
  +(electron_velocity[array_index+1]*electron_velocity[array_index+1])+(electron_velocity[array_index+2]*electron_velocity[array_index+2])) - electron_potential[e];

  electron_potential[e] += deltaE;
 
  return deltaE;
}

void ea_scattering(const int e, const int array_index, const int thread) {

    const double e_energy = electron_potential[e];
    // double deltaE = 1.0; //double(e_index+1)+DoS_cutoff - e_energy;
   if(e_energy < core_cutoff+1.5) return;
    double e_occupation;
    double f_e_occupation;
    double r_e_occupation;
    const double phonon_factor = 4.0*phonon_energy + (phonon_energy*mtrandom::gaussianc(omp_uniform_random[thread])/1.0);// ;//+ 0.16*abs(mtrandom::gaussianc(omp_uniform_random[thread]));// ;//
    const int e_index   = int(std::min( dos_size-2.0, std::max(1.0, floor((e_energy - DoS_cutoff)*i_phonon_energy))));
    const int f_e_index = int(std::min( dos_size-2.0, std::max(1.0, floor((e_energy + phonon_factor - DoS_cutoff)*i_phonon_energy))));
    const int r_e_index = int(std::min( dos_size-2.0, std::max(1.0, floor((e_energy - phonon_factor - DoS_cutoff)*i_phonon_energy))));
    // const double local_d_dos = std::min(local_dos_occ, phonon_energy*double(electron_nearest_electron_list[e][0])/(E_f_A-DoS_cutoff));
    if( e_energy > transport_cutoff) {
      e_occupation   = std::min(1.0, double(global_e_dos[e_index  ][0]) / (dos_standard[e_index]*phonon_energy)); 
      f_e_occupation = std::min(1.0, double(global_e_dos[f_e_index][0]) / (dos_standard[f_e_index]*phonon_energy)); 
      r_e_occupation = std::min(1.0, double(global_e_dos[r_e_index][0]) / (dos_standard[r_e_index]*phonon_energy)); 
      // e_occupation   =  std::min(1.0,double(ee_dos_hist[e][e_index  ]) / local_d_dos); 
      // f_e_occupation =  std::min(1.0,double(ee_dos_hist[e][e_index+1]) / local_d_dos); 
      // r_e_occupation =  std::min(1.0,double(ee_dos_hist[e][e_index-1]) / local_d_dos); 
    } else {
     e_occupation   = std::min(1.0, double(global_e_dos[e_index  ][0]) / std::max(1.0, double(global_e_dos[e_index  ][1])));    
     f_e_occupation = std::min(1.0, double(global_e_dos[f_e_index][0]) / std::max(1.0, double(global_e_dos[f_e_index][1]))); 
     r_e_occupation = std::min(1.0, double(global_e_dos[r_e_index][0]) / std::max(1.0, double(global_e_dos[r_e_index][1]))); 
    }
    
   const double thermal_factor = return_BE_distribution(phonon_factor, TTMp);
   const double f_factor = thermal_factor*(1.0 - f_e_occupation);// - f_e_occupation*(1.0-e_occupation);
   const double r_factor = (thermal_factor + 1.0) *(1.0-r_e_occupation);
   global_tau_ep[2*e_index] += ea_rate*(f_factor);
   global_tau_ep[2*e_index+1] -=  ea_rate*(r_factor );
   if(f_factor > 0.0 && omp_uniform_random[thread]() > exp(ea_rate*f_factor)) {
      double deltaE = phonon_factor;
      
      if( e_energy + deltaE > (core_cutoff+60.0) ) return;
      // if(omp_uniform_random[thread]()*factor < r_factor) deltaE *= -1.0;
      //  std::cout << e_occupation << ", " << d_e_occupation << ", " << double(e_index+1)+DoS_cutoff - e_energy << ", " << thermal_factor << ", "<< return_BE_integrand(abs(deltaE),Tp)*(e_occupation-d_e_occupation) << ", " << d_e_occupation*(1.0-e_occupation) << ", " << exp(ea_rate*abs(thermal_factor)) << std::endl;  
      relaxation_time_hist_ee[array_index].at(int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25)))) )++;
      relaxation_time_hist_ee[array_index + 2].at(int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25)))) ) += current_time_step - relaxation_time_hist_ee[array_index + 1].at(int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25)))) );
      
      electron_potential[e] += deltaE;
      
      relaxation_time_hist_ee[array_index + 1].at(int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((electron_potential[e]-DoS_cutoff)/0.25)))) ) = current_time_step;
      
      // electron_transport_list[e] = true;
      // if(electron_potential[e] < transport_cutoff) electron_transport_list[e] = false;

      const double theta = 2.0*M_PI*omp_uniform_random[thread]();
      const double phi   = M_PI*omp_uniform_random[thread]();
     // if(electron_velocity[array_index+2] < 0.0) theta += M_PI;

      // const double scattering_velocity = return_vel(electron_potential[e]);
      electron_velocity[array_index]   = cos(theta)*sin(phi);
      electron_velocity[array_index+1] = sin(theta)*sin(phi);
      electron_velocity[array_index+2] = cos(phi); 
    
      electron_ee_scattering_list[e][0] = 1;

      #pragma omp critical(eascattering)
      {
      //lattice_output << phonon_factor*i_phonon_energy/2.0 << ", " << f_factor << ", " << r_factor << ", " << thermal_factor << ", " << e_occupation << ", " <<\
      f_e_occupation << ", " << r_e_occupation << ", " << \
      thermal_factor*(return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), Te)-return_fermi_distribution((e_index*phonon_energy+core_cutoff+phonon_factor-E_f_A), Te))\
    - return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), Te)*(1.0 - return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), Te))  << ", " << \
      thermal_factor*(return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), Te)-return_fermi_distribution((e_index*phonon_energy+core_cutoff-phonon_factor-E_f_A), Te))\
    + return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), Te)*(1.0 - return_fermi_distribution((e_energy-phonon_factor-E_f_A), Te)) << ", " << \
      return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), Te) << ", " << \
      return_fermi_distribution((e_index*phonon_energy+core_cutoff+phonon_factor-E_f_A), Te) << ", " <<\
      return_fermi_distribution((e_index*phonon_energy+core_cutoff-phonon_factor-E_f_A), Te) << std::endl;
      //f(deltaE < 0.0) ea_core_scattering_count++;
      ea_transport_scattering_count++;
      TEKE += deltaE;
      TLE -= deltaE;
      e_a_scattering_count++;
      }
      return;
    }  
    if(r_factor > 0.0 && omp_uniform_random[thread]() > exp(ea_rate*r_factor)) {
      double deltaE = phonon_factor;
      if(e_energy - deltaE  < core_cutoff ) return;
     // if(omp_uniform_random[thread]()*factor < r_factor) deltaE *= -1.0;
      //  std::cout << e_occupation << ", " << d_e_occupation << ", " << double(e_index+1)+core_cutoff - e_energy << ", " << thermal_factor << ", "<< return_BE_integrand(abs(deltaE),Tp)*(e_occupation-d_e_occupation) << ", " << d_e_occupation*(1.0-e_occupation) << ", " << exp(ea_rate*abs(thermal_factor)) << std::endl;  
      relaxation_time_hist_ee[array_index][int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))))]++;
      relaxation_time_hist_ee[array_index + 2][int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))))] += current_time_step - relaxation_time_hist_ee[array_index + 1][int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))))];
      
      electron_potential[e] -= deltaE;
      
      relaxation_time_hist_ee[array_index + 1][int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((electron_potential[e]-DoS_cutoff)/0.25))))] = current_time_step;
      
      // electron_transport_list[e] = true;
      // if(electron_potential[e] < transport_cutoff) electron_transport_list[e] = false;

      const double theta = 2.0*M_PI*omp_uniform_random[thread]();
      const double phi   = M_PI*omp_uniform_random[thread]();
     // if(electron_velocity[array_index+2] < 0.0) theta += M_PI;

      // const double scattering_velocity = return_vel(electron_potential[e]);
      electron_velocity[array_index]   = cos(theta)*sin(phi);
      electron_velocity[array_index+1] = sin(theta)*sin(phi);
      electron_velocity[array_index+2] = cos(phi); 
    
      electron_ee_scattering_list[e][0] = 1;

      #pragma omp critical(eascattering)
      {
     //lattice_output << phonon_factor*i_phonon_energy/2.0 << ", " << f_factor << ", " << r_factor << ", " << thermal_factor << ", " << e_occupation << ", " <<\
      f_e_occupation << ", " << r_e_occupation << ", " << \
      thermal_factor*(return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), Te)-return_fermi_distribution((e_index*phonon_energy+core_cutoff+phonon_factor-E_f_A), Te))\
    - return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), Te)*(1.0 - return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), Te))  << ", " << \
      thermal_factor*(return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), Te)-return_fermi_distribution((e_index*phonon_energy+core_cutoff-phonon_factor-E_f_A), Te))\
    + return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), Te)*(1.0 - return_fermi_distribution((e_energy-phonon_factor-E_f_A), Te)) << ", " << \
      return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), Te) << ", " << \
      return_fermi_distribution((e_index*phonon_energy+core_cutoff+phonon_factor-E_f_A), Te) << ", " <<\
      return_fermi_distribution((e_index*phonon_energy+core_cutoff-phonon_factor-E_f_A), Te) << std::endl;
      ea_core_scattering_count++;
     // else ea_transport_scattering_count++;
      TEKE -= deltaE;
      TLE += deltaE;
      e_a_scattering_count++;
      }
      return;
    }  
}

void ee_scattering() {
         if (err::check) std::cout << "ee_scattering." << std::endl;

  omp_set_dynamic(0);
  omp_set_num_threads(omp_threads);
  int e_excess_p = 0;
  int e_excess_n = 0;
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
         bool e_ballistic = (e_energy > (E_f_A+4.8)) ? true : false;
         bool d_e_ballistic = (d_e_energy > (E_f_A+4.8)) ? true : false;

         if((e_ballistic != d_e_ballistic) && (e_ballistic || d_e_ballistic) ) inelastic_scattering(thread, electron, array_index, electron_collision, array_index_i, e_energy, d_e_energy);
         else elastic_scattering(thread, electron, array_index, electron_collision, array_index_i, e_energy, d_e_energy);

      }
    } 
  }   
  } if (err::check) std::cout << "ee_scattering done." << std::endl;
  
}

void elastic_scattering(int thread, int e, int array_index, int d_e, int array_index_i, double e_energy, double d_e_energy ) {
   
   double x_distance = electron_position[array_index]  - electron_position[array_index_i];
   double y_distance = electron_position[array_index+1]- electron_position[array_index_i+1];
   double z_distance = electron_position[array_index+2]- electron_position[array_index_i+2]; 

   if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance += lattice_width;
   else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance -= lattice_width;
   
   if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance += lattice_depth;
   else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance -= lattice_depth;    
   
   if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance += lattice_height;
   else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance -= lattice_height;

   double length  = (x_distance*x_distance + y_distance*y_distance + z_distance*z_distance);
   if(length == 0.0) return;
        
   const double k_1 = return_dWdE(e_energy);
   const double k_2 = return_dWdE(d_e_energy);

   const double k_x_1 = electron_velocity[array_index]*k_1;
   const double k_y_1 = electron_velocity[array_index+1]*k_1;
   const double k_z_1 = electron_velocity[array_index+2]*k_1;

   const double k_x_2 = electron_velocity[array_index_i]*k_2;
   const double k_y_2 = electron_velocity[array_index_i+1]*k_2;
   const double k_z_2 = electron_velocity[array_index_i+2]*k_2;

   const double v_x_dot_product =(( k_x_1 - k_x_2)*x_distance\
                                 + (k_y_1 - k_y_2)*y_distance\
                                 + (k_z_1 - k_z_2)*z_distance) /length;

          if(v_x_dot_product != v_x_dot_product) {std::cout << "ee error: " << v_x_dot_product << std::endl; return;}
        
   const double d_k_x_1 = k_x_1 - x_distance*v_x_dot_product;
   const double d_k_y_1 = k_y_1 - y_distance*v_x_dot_product;
   const double d_k_z_1 = k_z_1 - z_distance*v_x_dot_product;
   const double d_k_x_2 = k_x_2 + x_distance*v_x_dot_product;
   const double d_k_y_2 = k_y_2 + y_distance*v_x_dot_product;
   const double d_k_z_2 = k_z_2 + z_distance*v_x_dot_product;

   double d_k_1 = sqrt(d_k_x_1*d_k_x_1 + d_k_y_1*d_k_y_1 + d_k_z_1*d_k_z_1);
   double d_k_2 = sqrt(d_k_x_2*d_k_x_2 + d_k_y_2*d_k_y_2 + d_k_z_2*d_k_z_2);
        
   double d_e_1 = return_dWdE_i(d_k_1);
   double d_e_2 = return_dWdE_i(d_k_2);

   double deltaE = e_energy+d_e_energy - d_e_1 - d_e_2;
   const double deltaK = (k_1 - d_k_1)*(k_1 - d_k_1);
   double a_factor = abs((k_1*k_1 + k_2*k_2 - d_k_1*d_k_1 - d_k_2*d_k_2)/(2.0*d_k_1*d_k_2));
   double b_factor = k_1*k_2/(d_k_1*d_k_2);

      if(d_e_1 < core_cutoff || d_e_2 < core_cutoff) return;
      if( d_e_1 > (core_cutoff+60.0) || d_e_2 > (core_cutoff+60.0)) return;
      if(abs(deltaE) > 0.01) return; //std::cout << "ee collision deltaE: " << deltaE << ", " << e_energy << ", " << d_e_energy << ", " << d_e_1 << ", " << d_e_2 << std::endl;
      if(a_factor > (b_factor+1)) {
          std::cout << a_factor << ", " <<  b_factor << std::endl;
          return;}
       
      double d_e_occupation;  
      double d_d_occupation;
      const int e_index = int(std::min( dos_size-1.0, std::max(0.0, floor((d_e_1 - DoS_cutoff)*i_phonon_energy))));
      if ( d_e_1 > transport_cutoff) d_e_occupation = std::max(0.0, 1.0 - double(global_e_dos[e_index][0])/(dos_standard[e_index]*phonon_energy));  
      else  d_e_occupation = std::max(0.0, 1.0 - double(global_e_dos[e_index][0]) / double(global_e_dos[e_index][1]));  

      const int d_index = int(std::min( dos_size-1.0, std::max(0.0, floor((d_e_2 - DoS_cutoff)*i_phonon_energy))));
      if ( d_e_2 > transport_cutoff)  d_d_occupation = std::max(0.0, 1.0 - double(global_e_dos[d_index][0])/(dos_standard[d_index]*phonon_energy));  
      else d_d_occupation = std::max(0.0, 1.0  - double(global_e_dos[d_index][0]) / double(global_e_dos[d_index][1])); 
      double occupation_factor = ee_rate*d_e_occupation*d_d_occupation;///((q_sq+(deltaK))*(q_sq+(deltaK)));//*exp(0.15*(d_occupation*e_occupation-1.0));

      if(e_energy > E_f_A+4.8 || d_e_energy > E_f_A+4.8) occupation_factor /= q_sq*q_sq;
      else occupation_factor /= (q_sq+deltaK)*(q_sq+deltaK);

      global_tau_ee[e_index] += occupation_factor;
      global_tau_ee[d_index] += occupation_factor;

   if(omp_uniform_random[thread]() > exp(occupation_factor)) {
          
          // if(a_factor > (b_factor+1.0) || abs(deltaE) < 0.1) std::cout << a_factor << " > " << b_factor+1 << ", dE: "<< deltaE << ", " << k_1 << ", " << k_2 << ", " << d_k_1 << ", " << d_k_2 << ", " << e_energy << ", " << d_e_energy << ", " << d_e_1 << ", " << d_e_2 << std::endl;
          relaxation_time_hist_ee[array_index][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))))]++;
          relaxation_time_hist_ee[array_index_i][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_energy-DoS_cutoff)/0.25))))]++;
        
          relaxation_time_hist_ee[array_index + 2][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))))] += current_time_step - relaxation_time_hist_ee[array_index + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))))];
          relaxation_time_hist_ee[array_index_i + 2][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_energy-DoS_cutoff)/0.25))))] += current_time_step - relaxation_time_hist_ee[array_index_i + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_energy-DoS_cutoff)/0.25))))];
        
          electron_potential[e] = d_e_1;
          electron_potential[d_e] = d_e_2;

          relaxation_time_hist_ee[array_index + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_1-DoS_cutoff)/0.25))))] = current_time_step;
          relaxation_time_hist_ee[array_index_i + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_2-DoS_cutoff)/0.25))))] = current_time_step;

          electron_ee_scattering_list[e][0] = 1;
          electron_ee_scattering_list[d_e][0] = 1;

          double theta = atan2(d_k_y_1, d_k_x_1);
          if(theta != theta) theta == 0.0;
          double phi = acos(d_k_z_1/d_k_1);
         
          electron_velocity[array_index]    = sin(phi)*cos(theta);
          electron_velocity[array_index+1]  = sin(phi)*sin(theta);
          electron_velocity[array_index+2]  = cos(phi); 

           theta = atan2(d_k_y_2, d_k_x_2);
              if(theta != theta) theta == 0.0;
           phi = acos(d_k_z_2/d_k_2);

          electron_velocity[array_index_i]  = sin(phi)*cos(theta);
          electron_velocity[array_index_i+1]= sin(phi)*sin(theta);
          electron_velocity[array_index_i+2]= cos(phi); 

          #pragma omp critical(eescattering)
          {
          if(deltaE > 0.001) full_int_var++;
          else if (deltaE < -0.001) full_int_var--;
          
          if (electron_potential[e] < transport_cutoff) ee_core_scattering_count++;
          else ee_transport_scattering_count++;
          if (electron_potential[d_e] < transport_cutoff) ee_core_scattering_count++;
          else ee_transport_scattering_count++;
          // lattice_output << e_e_scattering_count << ", " << electron << ", " << electron_collision << ", " \
                         << x_distance*v_x_dot_product << ", " << y_distance*v_x_dot_product << ", " << z_distance*v_x_dot_product << ", " \
                         << deltaE << ", " << deltaK << ", "  << sqrt(length) << ", " << v_x_dot_product << ", " \
                         << ee_rate*e_occupation*d_e_occupation/((0.25+(deltaK))*(0.25+(deltaK))) << ", " << ee_rate*e_occupation*d_e_occupation/((0.25+(deltaE*deltaE))*(0.25+(deltaE*deltaE))) << ", " << rspace_Hsr << std::endl;
          e_e_scattering_count += 2;
          //if(electron_potential[electron] > transport_cutoff || electron_potential[electron_collision] > transport_cutoff )lattice_output <<\
            electron << ", " << e_occupation << ", " << d_e_occupation << ", " << 1.0/((q_sq+(deltaK*deltaK))*(q_sq+(deltaK*deltaK))) << ", " << deltaK << ", " << deltaE << ", " << e_energy<< ", " << d_e_energy << ", " <<\
            1-return_fermi_distribution(e_index*phonon_energy+deltaE+core_cutoff-E_f_A,constants::kB_r*Te) << ", " << 1-return_fermi_distribution(d_index*phonon_energy+core_cutoff-E_f_A-deltaE,constants::kB_r*Te) << ", " <<\
            sqrt(2.0*constants::m_e_r_i*e_energy/constants::hbar_r/constants::hbar_r) << ", " << sqrt(2.0*constants::m_e_r_i*(e_energy+deltaE)/constants::hbar_r/constants::hbar_r) << ", " << constants::m_over_hbar_sqrt*sqrt(v_x*v_x + v_y*v_y + v_z*v_z) << ", " << constants::m_over_hbar_sqrt*sqrt(d_v_x_2*d_v_x_2 + d_v_y_2*d_v_y_2 + d_v_z_2*d_v_z_2) << std::endl;
          // if(e_back_pressure > 1.0) e_e_scattering_count += 2;
          }
      return;
   }
}

void inelastic_scattering(int thread, int e, int array_index, int d_e, int array_index_i, double e_energy, double d_e_energy) {
   
   double deltaE = 0.1*(e_energy-d_e_energy);
   double d_e_occupation;  
   double d_d_occupation;
   double d_e_1 = e_energy-deltaE;
   double d_e_2 = d_e_energy+deltaE;
   const int e_index = int(std::min( dos_size-1.0, std::max(0.0, floor((d_e_1 - DoS_cutoff)*i_phonon_energy))));
   if ( d_e_1 > transport_cutoff) d_e_occupation = std::max(0.0, 1.0 - double(global_e_dos[e_index][0])/(dos_standard[e_index]*phonon_energy));  
   else  d_e_occupation = std::max(0.0, 1.0 - double(global_e_dos[e_index][0]) / double(global_e_dos[e_index][1]));  

   const int d_index = int(std::min( dos_size-1.0, std::max(0.0, floor((d_e_2 - DoS_cutoff)*i_phonon_energy))));
   if ( d_e_2 > transport_cutoff)  d_d_occupation = std::max(0.0, 1.0 - double(global_e_dos[d_index][0])/(dos_standard[d_index]*phonon_energy));  
   else d_d_occupation = std::max(0.0, 1.0  - double(global_e_dos[d_index][0]) / double(global_e_dos[d_index][1])); 
   double occupation_factor = ee_rate*d_e_occupation*d_d_occupation/(q_sq*q_sq);///((q_sq+(deltaK))*(q_sq+(deltaK)));//*exp(0.15*(d_occupation*e_occupation-1.0));

      global_tau_ee[e_index] += occupation_factor;
      global_tau_ee[d_index] += occupation_factor;

   if(omp_uniform_random[thread]() > exp(occupation_factor)) {
          
      relaxation_time_hist_ee[array_index][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))))]++;
      relaxation_time_hist_ee[array_index_i][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_energy-DoS_cutoff)/0.25))))]++;
        
      relaxation_time_hist_ee[array_index + 2][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))))] += current_time_step - relaxation_time_hist_ee[array_index + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((e_energy-DoS_cutoff)/0.25))))];
      relaxation_time_hist_ee[array_index_i + 2][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_energy-DoS_cutoff)/0.25))))] += current_time_step - relaxation_time_hist_ee[array_index_i + 1][int(std::max(0.0, std::min( 4.0*60.0 - 1.0, floor((d_e_energy-DoS_cutoff)/0.25))))];
        
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
          if (electron_potential[e] < transport_cutoff) ee_core_scattering_count++;
          else ee_transport_scattering_count++;
          if (electron_potential[d_e] < transport_cutoff) ee_core_scattering_count++;
          else ee_transport_scattering_count++;
          // lattice_output << e_e_scattering_count << ", " << electron << ", " << electron_collision << ", " \
                         << x_distance*v_x_dot_product << ", " << y_distance*v_x_dot_product << ", " << z_distance*v_x_dot_product << ", " \
                         << deltaE << ", " << deltaK << ", "  << sqrt(length) << ", " << v_x_dot_product << ", " \
                         << ee_rate*e_occupation*d_e_occupation/((0.25+(deltaK))*(0.25+(deltaK))) << ", " << ee_rate*e_occupation*d_e_occupation/((0.25+(deltaE*deltaE))*(0.25+(deltaE*deltaE))) << ", " << rspace_Hsr << std::endl;
          e_e_scattering_count += 2;
          //if(electron_potential[electron] > transport_cutoff || electron_potential[electron_collision] > transport_cutoff )lattice_output <<\
            electron << ", " << e_occupation << ", " << d_e_occupation << ", " << 1.0/((q_sq+(deltaK*deltaK))*(q_sq+(deltaK*deltaK))) << ", " << deltaK << ", " << deltaE << ", " << e_energy<< ", " << d_e_energy << ", " <<\
            1-return_fermi_distribution(e_index*phonon_energy+deltaE+core_cutoff-E_f_A,constants::kB_r*Te) << ", " << 1-return_fermi_distribution(d_index*phonon_energy+core_cutoff-E_f_A-deltaE,constants::kB_r*Te) << ", " <<\
            sqrt(2.0*constants::m_e_r_i*e_energy/constants::hbar_r/constants::hbar_r) << ", " << sqrt(2.0*constants::m_e_r_i*(e_energy+deltaE)/constants::hbar_r/constants::hbar_r) << ", " << constants::m_over_hbar_sqrt*sqrt(v_x*v_x + v_y*v_y + v_z*v_z) << ", " << constants::m_over_hbar_sqrt*sqrt(d_v_x_2*d_v_x_2 + d_v_y_2*d_v_y_2 + d_v_z_2*d_v_z_2) << std::endl;
          // if(e_back_pressure > 1.0) e_e_scattering_count += 2;
          }
      return;
   }
}

double k_sq() {
  double q_sq = 0.0;
  double f_0;
  // double f_2;
  double k;
  double df_dE;
  for(int h = 15; h < dos_size-2; h++) {
  
    if( ((h)*phonon_energy+DoS_cutoff) > transport_cutoff) f_0 = double(global_e_dos[h][0])/(dos_standard[h]*phonon_energy);
    else f_0 = double(global_e_dos[h][0])/double(global_e_dos[h][1]);

    k = return_dWdE(h*phonon_energy + DoS_cutoff);
    q_sq += phonon_energy*dos_standard[h]*f_0*1.0/(k*k*lattice_atoms*x_unit_size*y_unit_size*z_unit_size);
            // s^2 /  fs^2
  }

  return pow(sqrt(constants::K* q_sq) + (11.9-4.569471), 2.0);
}

} //end CASTLE namespace

