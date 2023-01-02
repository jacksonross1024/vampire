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
  
    TEKE_1 = 0.0;
    TLE_1  = 0.0;

    TEKE_2 = 0.0;
    TLE_2  = 0.0;

    layer_1 = 0;
    layer_2 = 0;
    
    p_x = 0.0;
    p_y = 0.0;
    p_z = 0.0;

          if(err::check)  std::cout << "Updating new electron position." << std::endl;
   update_position();

    n_f_1 = 2.0e3 * layer_1 / ((lattice_width) * (lattice_height-120.0) * (lattice_depth)); // e- / A**3 -> e-/m**3
    E_f_1 = constants::h * constants::h * pow(3.0 * M_PI * M_PI * n_f_1*1e27, 0.66666666666666666666666667) / (8.0 * M_PI * M_PI * constants::m_e); //Fermi-energy
    E_f_A_1 = E_f_1*1e20; //AJ
    n_f_2 = 2.0e3 * layer_2 / ((lattice_width) * (lattice_height-60.0) * (lattice_depth)); // e- / A**3 -> e-/m**3
    E_f_2 = constants::h * constants::h * pow(3.0 * M_PI * M_PI * n_f_2*1e27, 0.66666666666666666666666667) / (8.0 * M_PI * M_PI * constants::m_e); //Fermi-energy
    E_f_A_2 = E_f_2*1e20; //AJ

    e_specific_heat_1 = 0.5*M_PI*M_PI*constants::kB*constants::kB/E_f_1; // gamma; //J/K**2/e- 
    e_specific_heat_i_1 = 1.0 / e_specific_heat_1;
    e_heat_capacity_1 = (10.4)*(1.0+0.084)*1e20*e_specific_heat_1 * n_f_1; //J/K/e- -> AJ/K**2/nm**3
    e_heat_capacity_i_1 = 1.0 / e_heat_capacity_1;

    e_specific_heat_2 = 0.5*M_PI*M_PI*constants::kB*constants::kB/E_f_2; // gamma; //J/K**2/e- 
    e_specific_heat_i_2 = 1.0 / e_specific_heat_2;
    e_heat_capacity_2 = (10.4)*(1.0+0.084)*1e20*e_specific_heat_2 * n_f_2; //J/K/e- -> AJ/K**2/nm**3
    e_heat_capacity_i_2 = 1.0 / e_heat_capacity_2;

         if(err::check) std::cout << " positions updated. updating dos " << std::endl;
    update_dynamics();
         if(err::check)  std::cout << "Output mean data" << std::endl;
    
    if (current_time_step % CASTLE_output_rate == 0) output_data(); //std::cout << "x_flux: " << x_flux / CASTLE_output_rate << "\n"; x_flux = 0;
          if(err::check)  std::cout << "resetting global dos " << std::endl;

    const int size = ee_dos_hist[0].size();
    #pragma omp parallel for schedule(dynamic, 8)
    for(int e = 0; e < conduction_electrons; e++) {
      electron_ee_scattering_list[e][0] = 0;
      for(int h = 0; h < size; h++) {
        ee_dos_hist[e][h] = 0;
      }
    }

    int total_1 = 0;
    int total_2 = 0;
    for(int h = 0 ; h < dos_size; h++) {
      total_1 += global_e_dos_1[h][0];
      total_2 += global_e_dos_2[h][0];
      global_e_dos_1[h][0] = 0;
      global_e_dos_2[h][0] = 0;
    }

    if( total_1+total_2 != conduction_electrons) std::cout << "problem " << total_1+total_2 << std::endl;

    CASTLE_real_time += dt;
    current_time_step++;
   
    return EXIT_SUCCESS;
}

void setup_output() {

    CASTLE_output_data = true;
}


void update_position(){

    omp_set_dynamic(0);
    omp_set_num_threads(omp_threads);
    cell_integration_lists.swap(old_cell_integration_lists);
     
    for(int c = 0; c < total_cells; c++) {
      cell_integration_lists[c][0] = 1; 
    }
    
    int total_1 = 0;
    int total_2 = 0;
      if(err::check) std::cout << "cell lists swapped. updating position..." << std::endl;
  #pragma omp parallel reduction(+:x_flux,y_flux,z_flux, p_x,p_y,p_z, layer_1,layer_2)
  {
    std::vector<int> local_e_dos_1;
    local_e_dos_1.resize(dos_size,0);
    std::vector<int> local_e_dos_2;
    local_e_dos_2.resize(dos_size,0);
    const int thread = omp_get_thread_num();
   
    int local_1 = 0;
    int local_2 = 0;
    for (int l = 0 ; l < cells_per_thread; l++) {
      const int cell = lattice_cells_per_omp[thread][l];
      const int size = old_cell_integration_lists[cell][0];
   
      local_1 += size - 1;
      local_2 += size - 1;
  
      #pragma omp barrier
  
      for (int e = 1; e < size; e++) { 
        const int electron = old_cell_integration_lists[cell][e];
        
        const int array_index = 3*electron;

        double x_pos = electron_position[array_index];
        double y_pos = electron_position[array_index+1];
        double z_pos = electron_position[array_index+2];

        const double v_x = electron_velocity[array_index];
        const double v_y = electron_velocity[array_index+1];
        const double v_z = electron_velocity[array_index+2];

        if(z_pos > (lattice_height-60.0)) { local_e_dos_1[int(std::min(dos_size-1.0, std::max(0.0, floor((electron_potential[electron] - core_cutoff)/phonon_energy))))]++; layer_1++; }
        else { local_e_dos_2[int(std::min(dos_size-1.0, std::max(0.0, floor((electron_potential[electron] - core_cutoff)/phonon_energy))))]++; layer_2++; }
       
          p_x += v_x;
          p_y += v_y;
          p_z += v_z;
          x_pos += v_x * dt;// + (electron_force[array_index]   * dt * dt * constants::K_A / 2); // x superarray component
          y_pos += v_y * dt;// + (electron_force[array_index_y) * dt * dt * constants::K_A / 2); // y superarray component
          z_pos += v_z * dt;// + (electron_force[array_index_z) * dt * dt * constants::K_A / 2); // z superarray component

        if (x_pos < 0.0) {
          x_pos += lattice_width; 
          x_flux--;
          if (current_time_step % CASTLE_output_rate == 0) {
            #pragma omp atomic
            flux_index[int(std::max(0.0, std::min(99.0, floor((electron_potential[electron] - core_cutoff)/1.0))))]--;
          }
        } else if (x_pos > lattice_width) {
          x_pos -= lattice_width; 
          x_flux++;
          if (current_time_step % CASTLE_output_rate == 0) {
            #pragma omp atomic
            flux_index[int(std::max(0.0, std::min(99.0, floor((electron_potential[electron] - core_cutoff)/1.0))))]++; 
          } 
        }

        if (y_pos < 0.0) {y_pos += lattice_depth; y_flux--;}
        else if (y_pos > lattice_depth) {y_pos -= lattice_depth; y_flux++;}

        // if (z_pos < 0.0) {z_pos += lattice_height; z_flux--;}
        // else if (z_pos > lattice_height) {z_pos -= lattice_height; z_flux++;}

        if(z_pos < 0.0) { z_pos = 0.001; electron_velocity[array_index+2] *= -1.0;}
        else if (z_pos > lattice_height) { z_pos = lattice_height-0.001; electron_velocity[array_index+2] *= -1.0;}
        
        else if (electron_position[array_index+2] > (lattice_height-60.0) && z_pos < (lattice_height-60.0)) {
          if( omp_uniform_random[thread]() > exp((E_f_A_2-E_f_A_1)*108.099/electron_potential[e]) ) {
            // const double theta = omp_uniform_random[thread]() * 2.0 * M_PI;
            // const double phi   = omp_uniform_random[thread]() * M_PI; 
            // const double vel = sqrt(2.0*electron_potential[e]*constants::m_e_r_i);
            // electron_velocity[array_index]   = vel*cos(theta)*sin(phi);
            // electron_velocity[array_index+1] = vel*sin(theta)*sin(phi);
            // electron_velocity[array_index+2] = vel*cos(phi);
            // if( cos(phi) < 0.0)
              z_flux--; 
              layer_1--;
              layer_2++;
          } else { 
            z_pos = lattice_height-59.99; 
            electron_velocity[array_index+2] *= -1.0;
            }
        }
        else if (electron_position[array_index+2] < (lattice_height-60.0) && z_pos > (lattice_height-60.0)) {
          if( omp_uniform_random[thread]() > exp((E_f_A_1-E_f_A_2)*108.099/electron_potential[e]) ) {
            // const double theta = omp_uniform_random[thread]() * 2.0 * M_PI;
            // const double phi   = omp_uniform_random[thread]() * M_PI; 
            // const double vel = sqrt(2.0*electron_potential[e]*constants::m_e_r_i);
            // electron_velocity[array_index]   = vel*cos(theta)*sin(phi);
            // electron_velocity[array_index+1] = vel*sin(theta)*sin(phi);
            // electron_velocity[array_index+2] = vel*cos(phi);
            // if( cos(phi) > 0.0) 
             z_flux++; 
             layer_2--;
              layer_1++;
          }
          else { 
            z_pos = lattice_height-60.01; 
            electron_velocity[array_index+2] *= -1.0;
            }
        }

        electron_position[array_index]   = x_pos;
        electron_position[array_index+1] = y_pos;
        electron_position[array_index+2] = z_pos;
    
      // if(current_time_step % half_int_var == 0) {
          
           int x_cell = int(floor(x_pos / x_step_size));

           int y_cell = int(floor(y_pos / y_step_size));
    
           int z_cell = int(floor(z_pos / z_step_size));

          const int omp_cell = lattice_cell_coordinate[x_cell][y_cell][z_cell];
          
          if ((abs(x_cell - cell_lattice_coordinate.at(cell)[0]) > 1 &&\
               abs(x_cell - cell_lattice_coordinate.at(cell)[0]) <= (x_omp_cells-2)) ||
              (abs(y_cell - cell_lattice_coordinate.at(cell)[1]) > 1  &&\
               abs(y_cell - cell_lattice_coordinate.at(cell)[1]) <= (y_omp_cells-2)) ||\
              (abs(z_cell - cell_lattice_coordinate.at(cell)[2]) > 1  &&\
               abs(z_cell - cell_lattice_coordinate.at(cell)[2]) <= (z_omp_cells-2) )) {
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
        
        if(err::check && omp_get_thread_num()==0) std::cout << "positions updated. updating global dos..." << std::endl;
    #pragma omp critical 
    {
    total_1 += local_1;
    total_2 += local_2;
    for (int h = 0; h < dos_size; h++) {
      global_e_dos_1[h][0] += local_e_dos_1[h];
      global_e_dos_2[h][0] += local_e_dos_2[h];
    }
    } 
  }
      if(err::check) std::cout << "global dos updated. catching escaped electrons: " << escaping_electrons[0] << std::endl;
  // if(total != conduction_electrons) std::cout << total << ", " << conduction_electrons << std::endl;
  //if(current_time_step % half_int_var == 0) {    
    const int size = escaping_electrons[0];
  //  std::cout << "escaping electrons " << size-1 << ", " << escaping_electrons.size() << std::endl;
    
    for (int e = 1; e < size; e++) {
      const   int electron = escaping_electrons[e];
      const   int array_index = 3*electron;
      int x_cell = int(floor(electron_position[array_index] / x_step_size));

      int y_cell = int(floor(electron_position[array_index+1] / y_step_size));
    
      int z_cell = int(floor(electron_position[array_index+2] / z_step_size));
      
        const int omp_cell = lattice_cell_coordinate[x_cell][y_cell][z_cell];
        cell_integration_lists[omp_cell][cell_integration_lists[omp_cell][0]] = electron;
        cell_integration_lists[omp_cell][0]++;

      //  escaping_electrons[e] = 0;
    }
    escaping_electrons[0] = 1;
  
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
        pump = 1e3*photons_at_dt*photon_energy/(dt*lattice_depth*(lattice_height-120.0)*lattice_width); //AJ/fs/nm**3
        //external_potential = photon_energy; //AJ/hv     ;//1e27*pump*dt/ n_f; // AJ / particle
     //   std::cout << photons_at_dt << std::endl;
        TTMe_1 = d_TTMe_1;
        TTMp_1 = d_TTMp_1;
        TTMe_2 = d_TTMe_2;
        TTMp_2 = d_TTMp_2;
          //AJ/fs/K -> g(T-T)=C/t
        d_TTMe_1 = ((0.01*G*(TTMe_2 - TTMe_1)+G*(TTMp_1 - TTMe_1)+pump)*dt*e_heat_capacity_i_1/TTMe_1)+ TTMe_1;
        d_TTMp_1 = (  0.1*G*(TTMp_2 - TTMp_1)+G*(TTMe_1 - TTMp_1)      *dt*a_heat_capacity_i)       + TTMp_1;

          //AJ/fs/K -> g(T-T)=C/t
        d_TTMe_2 = ((0.01*G*(TTMe_1 - TTMe_2)+G*(TTMp_2 - TTMe_2)+pump)*dt*e_heat_capacity_i_2/TTMe_2)+ TTMe_2;
        d_TTMp_2 = (  0.1*G*(TTMp_1 - TTMp_2)+G*(TTMe_2 - TTMp_2)      *dt*a_heat_capacity_i)       + TTMp_2;
      }
    }
      

    dos_occ_1  = round(phonon_energy*(layer_1/(E_f_A_1-core_cutoff)));
    dos_occ_2  = round(phonon_energy*(layer_2/(E_f_A_2-core_cutoff)));

    omp_set_dynamic(0);
    omp_set_num_threads(omp_threads);
    std::vector<int> chosen;
    chosen.resize(photons_at_dt);
    while(count < photons_at_dt) {
      int select = int(omp_uniform_random[omp_get_thread_num()]()*2147483647) % (conduction_electrons);
      double e_energy = electron_potential[select];
      if( return_fermi_distribution((e_energy+photon_energy-E_f_A_1)/E_f_A_1, constants::kB_r*Te_1/E_f_A_1) < (1.0-1e-2) && electron_position[select*3 +2] > (lattice_height-60.0)) {
        chosen[count] = select;
        count++;
      }
    }
    count = 0;
      //  pump = 0.0;
    Tp_1 = d_Tp_1;
    Te_1 = d_Te_1;
    Tp_2 = d_Tp_2;
    Te_2 = d_Te_2;
      if(err::check) std::cout << "conditions set; updating dos " << std::endl;
    #pragma omp parallel for schedule(dynamic, 8)
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
      // if(!equilibrium_step) 
      // if(equilibrium_step && electron_potential[e] > 0.8*E_f_A) ea_scattering(e, array_index, omp_get_thread_num());
    //  else if (!equilibrium_step)  
      ea_scattering(e, array_index, omp_get_thread_num());
    }
    if(count != photons_at_dt) std::cout << photons_at_dt << ", " << count<< std::endl;
      //  TEKE += external_potential;
        if(err::check) std::cout << "dos updated. ee scattering step..." << std::endl;
    ee_scattering();
   // pump /= 1e-3*lattice_depth*lattice_height*lattice_width;
    d_Tp_1 = (a_heat_capacity_i*TLE_1 *n_f_1/layer_1) + Tp_1;
    d_Te_1 = (e_heat_capacity_i_1*TEKE_1*n_f_1/layer_1/abs(Te_1)) + (dt*e_heat_capacity_i_1*pump/Te_1) + Te_1;

    d_Tp_2 = (a_heat_capacity_i*TLE_2 *n_f_2/layer_2) + Tp_2;
    d_Te_2 = (e_heat_capacity_i_2*TEKE_2*n_f_2/layer_2/abs(Te_2)) + (dt*e_heat_capacity_i_2*pump/Te_2) + Te_2;
 
        if (err::check) std::cout << "reset scattering." << std::endl;
}

void electron_thermal_field(const int e, const int array_index, const double EKE, const int thread) {       
    
    const double theta = omp_uniform_random[thread]() * 2.0 * M_PI;
    const double phi   = omp_uniform_random[thread]() * M_PI; 

    electron_potential[e] += EKE;
    
    const double vel = sqrt(2.0*electron_potential[e]*constants::m_e_r_i);
    electron_velocity[array_index]   = vel*cos(theta)*sin(phi);
    electron_velocity[array_index+1] = vel*sin(theta)*sin(phi);
    electron_velocity[array_index+2] = vel*cos(phi);

    electron_ee_scattering_list[e][0] = 1;
}

/*
// void e_a_coulomb(const int& e, const int& array_index]{

//     double x_distance;
//     double y_distance;
//     double z_distance;
   
//     double length;
 
//     int array_index_a, nearest_electron_count = 1;
//     int count = 2;

//   for (int a = 0; a < lattice_atoms; a++) {

//         array_index_a = 3*a;
//         x_distance = electron_position[array_index]     - atom_position[array_index_a);
//         y_distance = electron_position[array_index + 1) - atom_position[array_index_a + 1);
//         z_distance = electron_position[array_index + 2) - atom_position[array_index_a + 2); 
       
//         if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
//         else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

//         if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
//         else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

//         if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
//         else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

//         length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
        
//         if(length > e_a_neighbor_cutoff) continue;

//         atomic_nearest_electron_list[e][nearest_electron_count) = array_index_a;
//         nearest_electron_count++;

//         if (length > e_a_coulomb_cutoff) continue;
     
//         if(mean_radius[2*e) > length) mean_radius[2*e) = length;

//         if(ea_coupling) {
//             electron_ea_scattering_list[e][count) = a;
//             count++;
//         }
//   }

//  atomic_nearest_electron_list[e][0] = nearest_electron_count;
//  electron_ea_scattering_list[e][1] = count;
 
// }

// void neighbor_e_a_coulomb(const int& e, const int& array_index]{
                     
//     double x_distance;
//     double y_distance;
//     double z_distance;
    
//     double length;
//     int size = atomic_nearest_electron_list[e][0];
//     int count = 2;
  
   
//     for (int a = 1; a < size; a++) {
      
//         int array_index_a = atomic_nearest_electron_list[e][a);
        
//         x_distance = electron_position[array_index]     - atom_position[array_index_a);
//         y_distance = electron_position[array_index + 1) - atom_position[array_index_a + 1);
//         z_distance = electron_position[array_index + 2) - atom_position[array_index_a + 2); 
      
//         if (x_distance < (10.0 - lattice_width))       x_distance = x_distance + lattice_width;
//         else if (x_distance > (lattice_width - 10.0))  x_distance = x_distance - lattice_width;

//         if (y_distance < (10.0 - lattice_depth))       y_distance = y_distance + lattice_depth;
//         else if (y_distance > (lattice_depth - 10.0))  y_distance = y_distance - lattice_depth;

//         if (z_distance <  (10.0 - lattice_height))     z_distance = z_distance + lattice_height;
//         else if (z_distance > (lattice_height - 10.0)) z_distance = z_distance - lattice_height;

//       length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance); //Angstroms
   
//       if (length > e_a_coulomb_cutoff) continue;
     
//       if(mean_radius[2*e) > length) mean_radius[2*e) = length;
      
//       if(ea_coupling) {
//         electron_ea_scattering_list[e][count) = array_index_a / 3;
//         count++;
//       }
//     }
//    electron_ea_scattering_list[e][1] = count;
// }
*/
void e_e_coulomb(const int e, const int array_index) {
    // if(electron_potential[e] < 0.8*E_f_A) return;

    int x_cell = int(floor(electron_position[array_index] / x_step_size));
      //if (x_cell < 0 || x_cell > x_omp_cells) std::cout << electron_position[array_index] << ", " << x_step_size << ", " << floor(electron_position[array_index] / x_step_size) << std::endl;

    int y_cell = int(floor(electron_position[array_index+1] / y_step_size));
      //if (y_cell < 0 || y_cell > y_omp_cells) std::cout << electron_position[array_index+1] << ", " << y_step_size << ", " << floor(electron_position[array_index+1] / y_step_size) << std::endl;
    
    int z_cell = int(floor(electron_position[array_index+2] / z_step_size));
      //if (z_cell < 0 || z_cell > z_omp_cells) std::cout << electron_position[array_index+2] << ", " << z_step_size << ", " << floor(electron_position[array_index+2] / z_step_size) << std::endl;

    const int cell = lattice_cell_coordinate[x_cell][y_cell][z_cell];
    int ee_dos_count = 1;
    int ee_integration_count = 1;
    int ee_scattering_list = 1;
   
    for(int s = 0; s < 27; s++) {
      const int int_cell = cell_nearest_neighbor_list[cell][s];
      const int size = cell_integration_lists[int_cell][0];
    
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

        if (z_distance <  (boundary_conditions_cutoff - lattice_height))     continue;
        else if (z_distance > (lattice_height - boundary_conditions_cutoff)) continue;
        
        const double length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
            
        if(length == 0.0) continue;
       
            if (length != length) std::cout << x_distance << ", " << y_distance << ", " << z_distance <<  ", " << \
        electron_position[array_index] << ", " << electron_position[array_index_i] << ", " << electron_position[array_index+1] << ", " << electron_position[array_index_i + 1] << ", " <<\
        electron_position[array_index+2] << ", " <<  electron_position[array_index_i + 2] << std::endl;
        
        if(length > e_e_integration_cutoff) continue;
        electron_integration_list[e][ee_integration_count] = array_index_i;
          if(ee_integration_count >= electron_integration_list[e].size() - 2) {std::cout << e << ", " << ee_integration_count << " > " << electron_integration_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
          break; }
        ee_integration_count++;

        if(length > e_e_neighbor_cutoff) continue;
        electron_nearest_electron_list[e][ee_dos_count] = electron;
        ee_dos_count++;
   
        ee_dos_hist[e].at(int(std::max(0.0, std::min(dos_size-1.0, floor((electron_potential[electron]-core_cutoff)/phonon_energy)))))++;
        
        if(ee_dos_count >= electron_nearest_electron_list[e].size() - 2) {std::cout << e << ", " << ee_dos_count << " > " << electron_nearest_electron_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
          break; }

        if (length > e_e_coulomb_cutoff) continue; 
        electron_ee_scattering_list[e][ee_scattering_list*2] = electron;
        electron_ee_scattering_list[e][ee_scattering_list*2 +1] = length;
        if(ee_scattering_list >= electron_ee_scattering_list[e].size() - 3) {std::cout << e << ", " << ee_scattering_list << " > " << electron_ee_scattering_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
          break; }

        ee_scattering_list++;
      }
    }
    electron_integration_list[e][0] = ee_integration_count;
    electron_nearest_electron_list[e][0] = ee_dos_count;
    electron_ee_scattering_list[e][1] = ee_scattering_list;
  //  if(e == 0 ) std::cout << ee_integration_count << ", " << ee_dos_count << ", " << ee_scattering_list << std::endl;
  }
void neighbor_e_e_coulomb(const int e, const int array_index) {
    // if(electron_potential[e] < 0.8*E_f_A) return;
    const int size = electron_integration_list[e][0];
      int neighbor_count = 1;
      int scattering_count = 1;
 
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

        if (z_distance <  (boundary_conditions_cutoff - lattice_height))     continue;
        else if (z_distance > (lattice_height - boundary_conditions_cutoff)) continue;
        
        const double length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
        if(length == 0.0) continue;
        
       if (length != length) std::cout << x_distance << ", " << y_distance << ", " << z_distance <<  ", " << \
        electron_position[array_index] << ", " << electron_position[array_index_i] << ", " << electron_position[array_index+1] << ", " << electron_position[array_index_i + 1] << ", " <<\
        electron_position[array_index+2] << ", " <<  electron_position[array_index_i + 2] << std::endl;
        
        if(length > e_e_neighbor_cutoff) continue;
        electron_nearest_electron_list[e][neighbor_count] = array_index_i/3;
        neighbor_count++;
      
        ee_dos_hist[e].at(int(std::max(0.0, std::min(dos_size-1.0, floor((electron_potential[electron]-core_cutoff)/phonon_energy)))))++;
      
        if(neighbor_count >= electron_nearest_electron_list[e].size() - 2) {std::cout << e << ", " << neighbor_count << " > " << electron_nearest_electron_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
        break; }

        if (length > e_e_coulomb_cutoff) continue; 
        electron_ee_scattering_list[e][scattering_count*2] = array_index_i/3;
        electron_ee_scattering_list[e][scattering_count*2 +1] = length;
        scattering_count++;
        if(scattering_count >= electron_ee_scattering_list[e].size() - 3) {std::cout << e << ", " << scattering_count << " > " << electron_ee_scattering_list[e].size() << ", " << length << ", " << electron_potential[e]  << std::endl;
           std::cout << x_distance << ", " << y_distance << ", " << z_distance << std::endl;
           break; }
      //   std::cout << e << ", " << i << ", " << length << std::endl;
    }

    electron_nearest_electron_list[e][0] = neighbor_count;
    electron_ee_scattering_list[e][1] = scattering_count;
    //  if(e == 0 ) std::cout << size << ", " << neighbor_count << ", " << scattering_count << std::endl;
}
/*
// void a_a_coulomb(const int a, const int array_index, \
//                         double& a_x_force, double& a_y_force, double& a_z_force, double& LPE) {
   
//     double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
//  //   int neighbor_count = 1;

//         x_distance = atom_position[array_index];
//         y_distance = atom_position[array_index+1];
//         z_distance = atom_position[array_index+2];
      
//                // if (err::check)  std::cout << "Calculating atomic interactions" << std::endl;

//        // neighbor_count = 1;

//         for (int i = 0; i < lattice_atoms; i++) {
//             if (i == a) continue; //no self repulsion

//             x_distance -= atom_position[array_index];
//             y_distance -= atom_position[array_index + 1);
//             z_distance -= atom_position[array_index + 2);

//             length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
           
//             if(length > 6) continue;

//             //length = sqrt(length);
//          //   force = -2000*(length - 2);
//            // PE = 1000*(length - 2)*(length - 2);
//         //   if(e == 0) std::cout << force << ", " << length <<  std::endl;
 
//             //phi   = acos(z_distance / length);
//             //theta = atan(y_distance / x_distance);
//             if (x_distance < 0) theta += M_PI;

//             atom_position[array_index]     += force * cos(theta)*sin(phi) / atomic_mass;
//             atom_position[array_index + 1) += force * sin(theta)*sin(phi) / atomic_mass;
//             atom_position[array_index + 2) += force * cos(phi) / atomic_mass;
//         }
//    // atomic_nearest_atom_list[a][0] = neighbor_count;
//     LPE += PE/2;
//    // if(a == 100) std::cout << neighbor_count << std::endl;
// }

// void neighbor_a_a_coulomb(const int a, const int array_index, \
//                         double& a_x_force, double& a_y_force, double& a_z_force, double& LPE) {
    
//     int array_index_i;
//     double x_distance,y_distance,z_distance,length, force, theta,phi, PE = 0;
//     double size = atomic_nearest_atom_list[a][0];

//         x_distance = atom_position[array_index];
//         y_distance = atom_position[array_index+1];
//         z_distance = atom_position[array_index+2];

//         for (int i = 1; i < size; i++) {
//            // if (i == a) continue; //no self repulsion
           
//             array_index_i = atomic_nearest_atom_list[a][i);

//             x_distance -= atom_position[array_index_i];
//             y_distance -= atom_position[array_index_i + 1);
//             z_distance -= atom_position[array_index_i + 2); 

//             if (x_distance < -30)     x_distance = x_distance + 40;
//             else if (x_distance > 30) x_distance = x_distance - 40;
//             if (y_distance < -30)     y_distance = y_distance + 40;
//             else if (y_distance > 30) y_distance = y_distance - 40;
//             if (z_distance <  -30)    z_distance = z_distance + 40;
//             else if (z_distance > 30) z_distance = z_distance - 40;

//             length = (x_distance*x_distance) + (y_distance*y_distance) + (z_distance*z_distance);
             
//             // if (length > a_a_coulomb_cutoff) continue; 

//             length = sqrt(length);
//             force = 4*(-2*exp(4 - 2*length) + 2*exp(length + 2));
//             PE = 4*(exp(4 - 2*length) - (2*exp(2 -length)) + 1);
            
//             phi   = acos(z_distance / length);
//             theta = atan(y_distance / x_distance);
//             if (x_distance < 0) theta += M_PI;

//             a_x_force += force * cos(theta)*sin(phi) / atomic_mass;
//             a_y_force += force * sin(theta)*sin(phi) / atomic_mass;
//             a_z_force += force * cos(phi) / atomic_mass;

//         }
//     LPE += PE/2;
//    // if(a == 100) std::cout << size << std::endl;
// }
*/
double electron_applied_voltage(const int e, const int array_index, double& external_potential) {

  electron_velocity[array_index] += 5e-11*dt*applied_voltage*constants::e_A*constants::m_e_r_i; //V/m * q = F_x/kg = 0.5*a_x*dt = d_v_x
  const double deltaE = 0.5*constants::m_e_r*((electron_velocity[array_index]*electron_velocity[array_index])\
  +(electron_velocity[array_index+1]*electron_velocity[array_index+1])+(electron_velocity[array_index+2]*electron_velocity[array_index+2])) - electron_potential[e];

  electron_potential[e] += deltaE;
 
  return deltaE;
}
/*
// void aa_scattering() {
//         if (err::check) std::cout << "aa_scattering." << std::endl;
//   double scattering_prob;
  
//   const static double aa_rate = -10.0*dt/sim::ee_coupling_strength;
//     for(int e = 0; e < conduction_electrons; e++) {
  
//       if(atomic_nearest_atom_list[e][0]) continue;
      
//       double deltaE = atom_potential[e] - E_f_A; - 3.0*constants::kB_r*Tp;

//      // if(deltaE < 8.0) scattering_prob = 8.0;
//        scattering_prob = deltaE*deltaE;

//       if(omp_uniform_random[omp_get_thread_num())() > exp(aa_rate*scattering_prob)) {
        
//         deltaE *= 0.5-0.05*exp(-2.0*deltaE/(E_f_A-3.0*constants::kB_r*Tp));
//         int size = atomic_nearest_atom_list[e][1] - 2;
//         size = 2+ (omp_int_random[omp_get_thread_num())() % size);
//       //  if(size >= atomic_nearest_atom_list[e].size() -1) std::cout << "bus error: " << size << " < " << atomic_nearest_atom_list[e].size() << std::endl; 
//         int a = atomic_nearest_atom_list[e][size);
      
//        // if(a >= conduction_electrons -1) std::cout << "bus error: " << a << " < " << conduction_electrons-1 << std::endl; 
//         if(atomic_nearest_atom_list[a][0]) continue;

//         deltaE *= 0.5;
//         atom_potential[e] -= deltaE;
//         atom_potential[a) += deltaE;  

//         atomic_nearest_atom_list[e][0] = 1;
//         atomic_nearest_atom_list[a][0] = 1;
        
//         a_a_scattering_count++;
//       }
//     }
//     if (err::check) std::cout << "aa_scattering done." << std::endl;
// }
*/
void ea_scattering(const int e, const int array_index, const int thread) {

    const double e_energy = electron_potential[e];
    if( e_energy + phonon_energy > (core_cutoff+100.0) ) return;
    double e_occupation;
    double f_e_occupation;
    double r_e_occupation;
    const int e_index = int(std::min( dos_size-2.0, std::max(1.0, floor((e_energy - core_cutoff)/phonon_energy))));
    const double phonon_factor = phonon_energy*0.666667*(1.5-0.25*mtrandom::gaussianc(omp_uniform_random[thread]));
    const double local_d_dos = std::min(local_dos_occ, phonon_energy*double(electron_nearest_electron_list[e][0])/(E_f_A_1-core_cutoff));
    double thermal_factor;
    if(electron_position[array_index+2] > (lattice_height-60.0) ) {
      thermal_factor = return_BE_integrand(phonon_factor, Tp_1);
      if( e_energy > transport_cutoff_1) {
        e_occupation   = std::min(1.0,double(global_e_dos_1[e_index  ][0]) / dos_occ_1); 
        f_e_occupation = std::min(1.0,double(global_e_dos_1[e_index+1][0]) / dos_occ_1); 
        r_e_occupation = std::min(1.0,double(global_e_dos_1[e_index-1][0]) / dos_occ_1); 
      // e_occupation   =  std::min(1.0,double(ee_dos_hist[e][e_index  ]) / local_d_dos); 
      // f_e_occupation = std::min(1.0,double(ee_dos_hist[e][e_index+1]) / local_d_dos); 
      // r_e_occupation =  std::min(1.0,double(ee_dos_hist[e][e_index-1]) / local_d_dos); 
      } else {
      e_occupation   = std::min(1.0, double(global_e_dos_1[e_index  ][0]) / std::max(1.0, double(global_e_dos_1[e_index  ][1])));    
      f_e_occupation = std::min(1.0, double(global_e_dos_1[e_index+1][0]) / std::max(1.0, double(global_e_dos_1[e_index+1][1]))); 
      r_e_occupation = std::min(1.0, double(global_e_dos_1[e_index-1][0]) / std::max(1.0, double(global_e_dos_1[e_index-1][1]))); 
      // e_occupation   = std::min(1.0,double(global_e_dos_1[e_index  ][0]) / dos_occ_1); 
      // f_e_occupation = std::min(1.0,double(global_e_dos_1[e_index+1][0]) / dos_occ_1); 
      // r_e_occupation = std::min(1.0,double(global_e_dos_1[e_index-1][0]) / dos_occ_1); 
      }
    } else {
      thermal_factor = return_BE_integrand(phonon_factor, Tp_2);
      if( e_energy > transport_cutoff_2) {
        e_occupation   = std::min(1.0,double(global_e_dos_2[e_index  ][0]) / dos_occ_2); 
        f_e_occupation = std::min(1.0,double(global_e_dos_2[e_index+1][0]) / dos_occ_2); 
        r_e_occupation = std::min(1.0,double(global_e_dos_2[e_index-1][0]) / dos_occ_2); 
      // e_occupation   =  std::min(1.0, double(ee_dos_hist[e][e_index  ]) / local_d_dos); 
      // f_e_occupation = std::min(1.0, double(ee_dos_hist[e][e_index+1]) / local_d_dos); 
      // r_e_occupation =  std::min(1.0, double(ee_dos_hist[e][e_index-1]) / local_d_dos); 
      } else {
        e_occupation   = std::min(1.0, double(global_e_dos_2[e_index  ][0]) / std::max(1.0, double(global_e_dos_2[e_index  ][1])));    
        f_e_occupation = std::min(1.0, double(global_e_dos_2[e_index+1][0]) / std::max(1.0, double(global_e_dos_2[e_index+1][1]))); 
        r_e_occupation = std::min(1.0, double(global_e_dos_2[e_index-1][0]) / std::max(1.0, double(global_e_dos_2[e_index-1][1]))); 
        // e_occupation   = std::min(1.0,double(global_e_dos_2[e_index  ][0]) / dos_occ_2); 
        // f_e_occupation = std::min(1.0,double(global_e_dos_2[e_index+1][0]) / dos_occ_2); 
        // r_e_occupation = std::min(1.0,double(global_e_dos_2[e_index-1][0]) / dos_occ_2); 
      }
    }
    //  else {
    //   thermal_factor = return_BE_integrand(phonon_factor, Tp_2);
    //   e_occupation   = return_fermi_distribution(e_energy-E_f_A, constants::kB_r*Te_2);
    //   f_e_occupation = return_fermi_distribution(e_energy + phonon_factor - E_f_A, constants::kB_r*Te_2);
    //   r_e_occupation = return_fermi_distribution(e_energy - phonon_factor - E_f_A, constants::kB_r*Te_2);
    // } 
    const double f_factor = thermal_factor*(e_occupation - f_e_occupation) - f_e_occupation*(1.0-e_occupation);
    const double r_factor = thermal_factor*(e_occupation - r_e_occupation) + e_occupation*(1.0-r_e_occupation);
    // const double factor = f_factor + r_factor;
    
    if(f_factor > 0.0 && omp_uniform_random[thread]() > exp(ea_rate*f_factor)) {
      double deltaE = phonon_factor;
      // if(omp_uniform_random[thread]()*factor < r_factor) deltaE *= -1.0;
      //  std::cout << e_occupation << ", " << d_e_occupation << ", " << double(e_index+1)+core_cutoff - e_energy << ", " << thermal_factor << ", "<< return_BE_integrand(abs(deltaE),Tp)*(e_occupation-d_e_occupation) << ", " << d_e_occupation*(1.0-e_occupation) << ", " << exp(ea_rate*abs(thermal_factor)) << std::endl;  
      const int index = int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((electron_potential[e]-core_cutoff)/0.25))));
      if(electron_position[array_index+2] > (lattice_height-60.0)) {
        relaxation_time_hist_ee_1[3*e][index]++;
        relaxation_time_hist_ee_1[3*e + 2][index] += current_time_step - relaxation_time_hist_ee_1[3*e + 1][index];
        electron_potential[e] += deltaE; 
        relaxation_time_hist_ee_1[3*e + 1][int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((electron_potential[e]-core_cutoff)/0.25))))] = current_time_step;
      }  else {
        relaxation_time_hist_ee_2[3*e][index]++;
        relaxation_time_hist_ee_2[3*e + 2][index] += current_time_step - relaxation_time_hist_ee_2[3*e + 1][index];
        electron_potential[e] += deltaE; 
        relaxation_time_hist_ee_2[3*e + 1][int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((electron_potential[e]-core_cutoff)/0.25))))] = current_time_step;
      } 

      const double theta = 2.0*M_PI*omp_uniform_random[thread]();
      const double phi   = M_PI*omp_uniform_random[thread]();

      const double scattering_velocity = sqrt(2.0*electron_potential[e]*constants::m_e_r_i);
      electron_velocity[array_index]   = scattering_velocity * cos(theta)*sin(phi);
      electron_velocity[array_index+1] = scattering_velocity * sin(theta)*sin(phi);
      electron_velocity[array_index+2] = scattering_velocity * cos(phi); 
    
      electron_ee_scattering_list[e][0] = 1;

      #pragma omp critical(eascattering)
      {
     // lattice_output << e << ", " << f_factor << ", " << r_factor << ", " << thermal_factor << ", " << e_occupation << ", " <<\
      f_e_occupation << ", " << r_e_occupation << ", " << \
      thermal_factor*(return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), constants::kB_r*Te)-return_fermi_distribution((e_index*phonon_energy+core_cutoff+phonon_factor-E_f_A), constants::kB_r*Te))\
    - return_fermi_distribution((e_index*phonon_energy+core_cutoff+phonon_factor-E_f_A), constants::kB_r*Te)*(1.0 - return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), constants::kB_r*Te))  << ", " << \
    - 1.0*(thermal_factor*(return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), constants::kB_r*Te)-return_fermi_distribution((e_index*phonon_energy+core_cutoff-phonon_factor-E_f_A), constants::kB_r*Te))\
    + return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), constants::kB_r*Te)*(1.0 - return_fermi_distribution((e_index*phonon_energy+core_cutoff-phonon_factor-E_f_A), constants::kB_r*Te))) << ", " << \
      return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), constants::kB_r*Te) << ", " << \
      return_fermi_distribution((e_index*phonon_energy+core_cutoff+phonon_factor-E_f_A), constants::kB_r*Te) << ", " <<\
      return_fermi_distribution((e_index*phonon_energy+core_cutoff-phonon_factor-E_f_A), constants::kB_r*Te) << std::endl;
      //f(deltaE < 0.0) ea_core_scattering_count++;
      if(electron_position[array_index+2] > (lattice_height-60.0) ) { TEKE_1 += deltaE; TLE_1 -= deltaE; ea_transport_scattering_count_1++; }
      else {TEKE_2 += deltaE; TLE_2 -= deltaE; ea_transport_scattering_count_2++;}
      }
      return;
    }

    if(r_factor > 0.0 && omp_uniform_random[thread]() > exp(ea_rate*r_factor)) {
      double deltaE = phonon_factor;
      const int index = int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((electron_potential[e]-core_cutoff)/0.25))));
      if(electron_position[array_index+2] > (lattice_height-60.0)) {
        relaxation_time_hist_ee_1[3*e][index]++;
        relaxation_time_hist_ee_1[3*e + 2][index] += current_time_step - relaxation_time_hist_ee_1[3*e + 1][index];
        electron_potential[e] -= deltaE; 
        relaxation_time_hist_ee_1[3*e + 1][int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((electron_potential[e]-core_cutoff)/0.25))))] = current_time_step;
      } else {
        relaxation_time_hist_ee_2[3*e][index]++;
        relaxation_time_hist_ee_2[3*e + 2][index] += current_time_step - relaxation_time_hist_ee_2[3*e + 1][index];
        electron_potential[e] -= deltaE; 
        relaxation_time_hist_ee_2[3*e + 1][int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((electron_potential[e]-core_cutoff)/0.25))))] = current_time_step;
      } 
      const double theta = 2.0*M_PI*omp_uniform_random[thread]();
      const double phi   = M_PI*omp_uniform_random[thread]();
     // if(electron_velocity[array_index+2] < 0.0) theta += M_PI;

      const double scattering_velocity = sqrt(2.0*electron_potential[e]*constants::m_e_r_i);
      electron_velocity[array_index]   = scattering_velocity * cos(theta)*sin(phi);
      electron_velocity[array_index+1] = scattering_velocity * sin(theta)*sin(phi);
      electron_velocity[array_index+2] = scattering_velocity * cos(phi); 
    
      electron_ee_scattering_list[e][0] = 1;

      #pragma omp critical(eascattering)
      {
    // lattice_output << e << ", " << f_factor << ", " << r_factor << ", " << thermal_factor << ", " << e_occupation << ", " <<\
      f_e_occupation << ", " << r_e_occupation << ", " << \
      thermal_factor*(return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), constants::kB_r*Te)-return_fermi_distribution((e_index*phonon_energy+core_cutoff+phonon_factor-E_f_A), constants::kB_r*Te))\
    - return_fermi_distribution((e_index*phonon_energy+core_cutoff+phonon_factor-E_f_A), constants::kB_r*Te)*(1.0 - return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), constants::kB_r*Te))  << ", " << \
    -1.0*(thermal_factor*(return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), constants::kB_r*Te)-return_fermi_distribution((e_index*phonon_energy+core_cutoff-phonon_factor-E_f_A), constants::kB_r*Te))\
    + return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), constants::kB_r*Te)*(1.0 - return_fermi_distribution((e_index*phonon_energy+core_cutoff-phonon_factor-E_f_A), constants::kB_r*Te))) << ", " << \
      return_fermi_distribution((e_index*phonon_energy+core_cutoff-E_f_A), constants::kB_r*Te) << ", " << \
      return_fermi_distribution((e_index*phonon_energy+core_cutoff+phonon_factor-E_f_A), constants::kB_r*Te) << ", " <<\
      return_fermi_distribution((e_index*phonon_energy+core_cutoff-phonon_factor-E_f_A), constants::kB_r*Te) << std::endl;
      if(electron_position[array_index+2] > (lattice_height-60.0) ) {TEKE_1 -= deltaE; TLE_1 += deltaE; ea_core_scattering_count_1++;}
      else {TEKE_2 -= deltaE; TLE_2 += deltaE; ea_core_scattering_count_2++;}
      }
      return;
    }  
}

void ee_scattering() {
         if (err::check) std::cout << "ee_scattering." << std::endl;

  omp_set_dynamic(0);
  omp_set_num_threads(omp_threads);
  const static double q_sq = 1.4*1.4;
  #pragma omp parallel reduction(+: ee_core_scattering_count_1, ee_transport_scattering_count_1, ee_core_scattering_count_2, ee_transport_scattering_count_2) 
  {
    
  const int thread = omp_get_thread_num();
  for(int l = 0; l < cells_per_thread; l++) {
    const int cell = lattice_cells_per_omp[thread][l];
    const int size = cell_integration_lists[cell][0];
    
    #pragma omp barrier 

    for(int e = 1; e < size; e++) {
      const int electron = cell_integration_lists[cell][e];
      if (electron_ee_scattering_list[electron][0] == 1) continue;
      
      const int scattering_size = electron_ee_scattering_list[electron][1];
      const double e_energy = electron_potential[electron];
      
      const int array_index = 3*electron;

      for(int a = 1; a < scattering_size; a++) {
        int electron_collision = electron_ee_scattering_list[electron][a*2];
        if (electron_ee_scattering_list[electron_collision][0] == 1) continue;
        const double d_e_energy = electron_potential[electron_collision];
        
        const int array_index_i = 3*electron_collision;

        const double length = electron_ee_scattering_list[electron][2*a + 1];
        double x_distance = electron_position[array_index]  - electron_position[array_index_i];
        double y_distance = electron_position[array_index+1]- electron_position[array_index_i+1];
        double z_distance = electron_position[array_index+2]- electron_position[array_index_i+2]; 

        if (x_distance < (boundary_conditions_cutoff - lattice_width))       x_distance += lattice_width;
        else if (x_distance > (lattice_width - boundary_conditions_cutoff))  x_distance -= lattice_width;
        if (y_distance < (boundary_conditions_cutoff - lattice_depth))       y_distance += lattice_depth;
        else if (y_distance > (lattice_depth - boundary_conditions_cutoff))  y_distance -= lattice_depth;    
        if (z_distance <  (boundary_conditions_cutoff - lattice_height))     z_distance += lattice_height;
        else if (z_distance > (lattice_height - boundary_conditions_cutoff)) z_distance -= lattice_height;

        if(length == 0.0) continue;
        const double v_x = electron_velocity[array_index];
        const double v_y = electron_velocity[array_index+1];
        const double v_z = electron_velocity[array_index+2];
        const double v_x_2 = electron_velocity[array_index_i];
        const double v_y_2 = electron_velocity[array_index_i+1];
        const double v_z_2 = electron_velocity[array_index_i+2];

        const double v_x_dot_product = (( v_x - v_x_2)*x_distance\
                                       + (v_y - v_y_2)*y_distance\
                                       + (v_z - v_z_2)*z_distance) /length;

        const double d_v_x_1 = v_x_2 + x_distance*v_x_dot_product;
        const double d_v_y_1 = v_y_2 + y_distance*v_x_dot_product;
        const double d_v_z_1 = v_z_2 + z_distance*v_x_dot_product;
        const double d_v_x_2 = v_x   - x_distance*v_x_dot_product;
        const double d_v_y_2 = v_y   - y_distance*v_x_dot_product;
        const double d_v_z_2 = v_z   - z_distance*v_x_dot_product;
  
        double deltaE = (0.5*constants::m_e_r*(d_v_x_2*d_v_x_2 + d_v_y_2*d_v_y_2 + d_v_z_2*d_v_z_2) -\
                           0.5*constants::m_e_r*(v_x*v_x + v_y*v_y + v_z*v_z));
        if(electron_position[array_index+2] > (lattice_height-120.0) && (e_energy + deltaE < transport_cutoff_1)) continue;
        else if( e_energy + deltaE < transport_cutoff_2 ) continue;
        
        if(electron_position[array_index_i+2] > (lattice_height-120.0) && (d_e_energy - deltaE < transport_cutoff_1)) continue;
        else if(d_e_energy - deltaE < transport_cutoff_2 ) continue;

        if( (e_energy + deltaE) > (core_cutoff+100.0) || (d_e_energy - deltaE > (core_cutoff+100.0))) continue;
        
        double e_occupation;  
        double d_e_occupation; 
        const int e_index = int(std::min( dos_size-1.0, std::max(0.0, floor((e_energy + deltaE - core_cutoff)/phonon_energy))));
        const int d_index = int(std::min( dos_size-1.0, std::max(0.0, floor((d_e_energy - deltaE - core_cutoff)/phonon_energy))));
        
         //**global density of states**
        if(electron_position[array_index+2] > (lattice_height-60.0)) { 
          if ( (e_energy + deltaE) > transport_cutoff_1) e_occupation = std::max(0.0, 1.0 - double(global_e_dos_1[e_index][0])/dos_occ_1);  
          else e_occupation = std::max(0.0, 1.0 - double(global_e_dos_1[e_index][0]) / double(global_e_dos_1[e_index][1]));    
        } else {
          if ( (e_energy + deltaE) > transport_cutoff_2) e_occupation = std::max(0.0, 1.0 - double(global_e_dos_2[e_index][0])/dos_occ_2);  
          else e_occupation = std::max(0.0, 1.0 - double(global_e_dos_2[e_index][0]) / double(global_e_dos_2[e_index][1]));    
        } 
        // else e_occupation = 1.0-return_fermi_distribution(e_energy+deltaE-E_f_A, constants::kB_r*Te_2)
        if(electron_position[array_index_i+2] > (lattice_height-60.0)) {
          if ( (d_e_energy - deltaE) > transport_cutoff_1) d_e_occupation = std::max(0.0, 1.0 - double(global_e_dos_1[d_index][0])/dos_occ_1);  
          else d_e_occupation = std::max(0.0, 1.0  - double(global_e_dos_1[d_index][0]) / double(global_e_dos_1[d_index][1])); 
        } else {
          if ( (d_e_energy - deltaE) > transport_cutoff_2) d_e_occupation = std::max(0.0, 1.0 - double(global_e_dos_2[d_index][0])/dos_occ_2);  
          else d_e_occupation = std::max(0.0, 1.0  - double(global_e_dos_2[d_index][0]) / double(global_e_dos_2[d_index][1])); 
        } 
        // else d_e_occupation = 1.0-return_fermi_distribution(d_e_energy-deltaE-E_f_A, constants::kB_r*Te_2);
            
        //// **local density of states**
        // const double local_d_dos = std::min(local_dos_occ, phonon_energy*double(electron_nearest_electron_list[electron_collision][0])/(E_f_A_1-core_cutoff));
        // const double local_e_dos = std::min(local_dos_occ, phonon_energy*double(electron_nearest_electron_list[electron][0])/(E_f_A_1-core_cutoff));
        
        // e_occupation   = std::max(0.0, 1.0 - (double(ee_dos_hist[electron][e_index])/local_e_dos));     
        // d_e_occupation = std::max(0.0, 1.0 - (double(ee_dos_hist[electron_collision][d_index])/local_d_dos));
      
        const double deltaK = constants::m_over_hbar_sqrt*(sqrt(d_v_x_2*d_v_x_2 + d_v_y_2*d_v_y_2 + d_v_z_2*d_v_z_2)-sqrt(v_x*v_x+v_y*v_y+v_z*v_z));
          
        const double occupation_factor = e_occupation*d_e_occupation;
        if(occupation_factor < 0.0) continue;
        if(omp_uniform_random[thread]() > exp(ee_rate*occupation_factor/((q_sq+(deltaK*deltaK))*(q_sq+(deltaK*deltaK))))) {
          
          int relax_index = int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((electron_potential[electron]-core_cutoff)/0.25))));
          if(electron_position[array_index+2] > (lattice_height-60.0)){
            relaxation_time_hist_ee_1[3*electron][relax_index]++;
            relaxation_time_hist_ee_1[3*electron + 2][relax_index] += current_time_step - relaxation_time_hist_ee_1[3*electron + 1][relax_index]; 
            electron_potential[electron] += deltaE;     
            relaxation_time_hist_ee_1[3*electron + 1][int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((electron_potential[electron]-core_cutoff)/0.25))))] = current_time_step;
          } else {
            relaxation_time_hist_ee_2[3*electron][relax_index]++;
            relaxation_time_hist_ee_2[3*electron + 2][relax_index] += current_time_step - relaxation_time_hist_ee_2[3*electron + 1][relax_index]; 
            electron_potential[electron] += deltaE;     
            relaxation_time_hist_ee_2[3*electron + 1][int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((electron_potential[electron]-core_cutoff)/0.25))))] = current_time_step;
          }

          relax_index = int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((electron_potential[electron_collision]-core_cutoff)/0.25))));
          if(electron_position[array_index_i+2] > (lattice_height-60.0)) {
            relaxation_time_hist_ee_1[3*electron_collision][relax_index]++;
            relaxation_time_hist_ee_1[3*electron_collision + 2][relax_index] += current_time_step - relaxation_time_hist_ee_1[3*electron_collision + 1][relax_index];
            electron_potential[electron_collision]   -= deltaE;
            relaxation_time_hist_ee_1[3*electron_collision + 1][int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((electron_potential[electron_collision]-core_cutoff)/0.25))))] = current_time_step;
          } else {
            relaxation_time_hist_ee_2[3*electron_collision][relax_index]++;
            relaxation_time_hist_ee_2[3*electron_collision + 2][relax_index] += current_time_step - relaxation_time_hist_ee_2[3*electron_collision + 1][relax_index];
            electron_potential[electron_collision]   -= deltaE;
            relaxation_time_hist_ee_2[3*electron_collision + 1][int(std::max(0.0, std::min( 4.0*100.0 - 1.0, floor((electron_potential[electron_collision]-core_cutoff)/0.25))))] = current_time_step;
          }

          electron_ee_scattering_list[electron][0] = 1;
          electron_ee_scattering_list[electron_collision][0] = 1;
          electron_velocity[array_index_i]  = d_v_x_1;
          electron_velocity[array_index_i+1]= d_v_y_1; 
          electron_velocity[array_index_i+2]= d_v_z_1;
          electron_velocity[array_index]    = d_v_x_2;
          electron_velocity[array_index+1]  = d_v_y_2;
          electron_velocity[array_index+2]  = d_v_z_2; 

          #pragma omp critical(eescattering)
          {
          if(electron_position[array_index+2] > (lattice_height-60.0)) {
            if (electron_potential[electron] < transport_cutoff_1) ee_core_scattering_count_1++;
            else ee_transport_scattering_count_1++;
            if (electron_potential[electron_collision] < transport_cutoff_1) ee_core_scattering_count_1++;
            else ee_transport_scattering_count_1++;
          // lattice_output << e_e_scattering_count << ", " << electron << ", " << electron_collision << ", " \
                         << x_distance*v_x_dot_product << ", " << y_distance*v_x_dot_product << ", " << z_distance*v_x_dot_product << ", " \
                         << deltaE << ", " << deltaK << ", "  << sqrt(length) << ", " << v_x_dot_product << ", " \
                         << ee_rate*e_occupation*d_e_occupation/((0.25+(deltaK))*(0.25+(deltaK))) << ", " << ee_rate*e_occupation*d_e_occupation/((0.25+(deltaE*deltaE))*(0.25+(deltaE*deltaE))) << ", " << rspace_Hsr << std::endl;
          }  else {
            if (electron_potential[electron] < transport_cutoff_2) ee_core_scattering_count_2++;
            else ee_transport_scattering_count_2++;
            if (electron_potential[electron_collision] < transport_cutoff_2) ee_core_scattering_count_2++;
            else ee_transport_scattering_count_2++;
          // lattice_output << e_e_scattering_count << ", " << electron << ", " << electron_collision << ", " \
                         << x_distance*v_x_dot_product << ", " << y_distance*v_x_dot_product << ", " << z_distance*v_x_dot_product << ", " \
                         << deltaE << ", " << deltaK << ", "  << sqrt(length) << ", " << v_x_dot_product << ", " \
                         << ee_rate*e_occupation*d_e_occupation/((0.25+(deltaK))*(0.25+(deltaK))) << ", " << ee_rate*e_occupation*d_e_occupation/((0.25+(deltaE*deltaE))*(0.25+(deltaE*deltaE))) << ", " << rspace_Hsr << std::endl;
          }
          }
          break;
        }
      }
    } 
  }   if (err::check) std::cout << "ee_scattering done." << std::endl;
  }
}


} //end CASTLE namespace

