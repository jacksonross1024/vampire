//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <cmath>
#include <algorithm>
#include <iostream>

// Vampire headers
#include "internal.hpp"
#include "atoms.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "ltmp.hpp"
#include "create.hpp"

// Forward declaration - output functions are in output.cpp
namespace st { namespace internal {
   void output_thermal_microcell_data();
   void open_thermal_temperature_profile_file();
}}

namespace st {
namespace internal {

//-----------------------------------------------------------------------------
// Initialize Debye lookup table
//-----------------------------------------------------------------------------
// Debye kernel: x^4 exp(x) / (exp(x)-1)^2
static double debye_kernel(const double x) {
   if(x <= 0.0) return 0.0;
   const double exp_x = std::exp(x);
   const double exp_m1 = exp_x - 1.0;
   return x*x*x*x*exp_x/(exp_m1*exp_m1);
}

void initialise_debye_table() {
   sc1d_debye_table.resize(12001, 0.0); // T_D/T from 0 to 12.0 with 0.001 resolution

   // Accumulate the integral using Simpson's 3/8 rule (matching ltmp implementation)
   volatile double integrand = 0.0;
   for(double T_D_over_T = 0.0; T_D_over_T < 12.0; T_D_over_T += 0.001) {
      int int_resolution = static_cast<int>(std::round(T_D_over_T * 1000.0));
      if(int_resolution >= 0 && int_resolution < 12000) {
         double left = T_D_over_T + 0.00050;
         double right = T_D_over_T - 0.00050;
         // Simpson's 3/8 rule: for 4 points (right, mid2, mid1, left) spanning 0.001
         const double h = 0.001 / 3.0;
         const double simpson_coeff = 3.0 * h / 8.0;
         integrand += simpson_coeff*(debye_kernel(left)
            + 3.0*debye_kernel(2.0*left*0.333+right*0.333)
            + 3.0*debye_kernel(left*0.333+2.0*right*0.333)
            + debye_kernel(right));
         
         // Store Debye heat capacity coefficient: 3.0 * integral / (T_D/T)^3
         // At high T (T >> T_D), this approaches 1.0, giving Cp(T) → phonon_heat_capacity (constant)
         if(T_D_over_T > 1e-10) {
            sc1d_debye_table[int_resolution] = 3.0*integrand/(T_D_over_T*T_D_over_T*T_D_over_T);
         } else {
            sc1d_debye_table[int_resolution] = 1.0; // High temperature limit
         }
      }
   }
}

//-----------------------------------------------------------------------------
// Initialize thermal gradients
//-----------------------------------------------------------------------------
void initialise_thermal_gradients() {
   if(sc1d_thermal_gradients_initialised) {
      zlog << zTs() << "Warning: Thermal gradients already initialized. Continuing." << std::endl;
      return;
   }
   
   if(!sc1d_thermal_gradients_enable) {
      return;
   }
   
   zlog << zTs() << "Initializing thermal gradients for spin-currents module." << std::endl;
   
   // Initialize Debye lookup table
   initialise_debye_table();
   
   const int ncz = num_microcells_per_stack;
   const double dz = micro_cell_thickness * 1e-10; // Angstroms to meters
   
   // Allocate temperature arrays (per stack, per cell)
   // Initialize to simulation temperature (not reference temperature)
   // Reference temperature is only for parameter scaling
   const int total_stacks = num_stacks_y;
   const size_t total_cells = total_stacks * ncz;
   const double initial_T = sim::temperature; // Use actual simulation temperature
   
   sc1d_Te.resize(total_cells, initial_T);
   sc1d_Tp.resize(total_cells, initial_T);
   sc1d_sqrt_Te.resize(total_cells, std::sqrt(initial_T));
   sc1d_sqrt_Tp.resize(total_cells, std::sqrt(initial_T));
   
   // Allocate material property arrays (per cell, same for all stacks)
   sc1d_Ce.resize(ncz, 0.0);
   sc1d_Cp.resize(ncz, 0.0);
   sc1d_kappa_e.resize(ncz, 0.0);
   sc1d_kappa_p.resize(ncz, 0.0);
   sc1d_G.resize(ncz, 0.0);
   sc1d_T_Debye.resize(ncz, 0.0);
   
   // Allocate atom mapping arrays
   #ifdef MPICF
      const int num_local_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
   #else
      const int num_local_atoms = atoms::num_atoms;
   #endif
   
   sc1d_atom_cell_idx.resize(num_local_atoms, -1);
   sc1d_atom_use_phonon.resize(num_local_atoms, false);
   
   // Count atoms per cell
   std::vector<int> num_atoms_in_cell(ncz, 0);
   
   // Map atoms to coarse cells and accumulate material properties
   for(int atom = 0; atom < num_local_atoms; ++atom) {
      const double z_angstrom = atoms::z_coord_array[atom];
      const int cell = static_cast<int>(std::floor(z_angstrom / micro_cell_thickness));
      
      // Clamp cell index to valid range
      const int cell_idx = std::max(0, std::min(cell, ncz - 1));
      sc1d_atom_cell_idx[atom] = cell_idx;
      
      // Determine if atom couples to phonon temperature
      const int mat = atoms::type_array[atom];
      sc1d_atom_use_phonon[atom] = mp::material[mat].couple_to_phonon_temperature;
      
      // Accumulate material properties from spin-torque material structure
      if(mat < (int)st::internal::mp.size()) {
         sc1d_Ce[cell_idx] += st::internal::mp[mat].electron_heat_capacity;
         sc1d_Cp[cell_idx] += st::internal::mp[mat].phonon_heat_capacity;
         sc1d_kappa_e[cell_idx] += st::internal::mp[mat].electron_thermal_conductivity;
         sc1d_kappa_p[cell_idx] += st::internal::mp[mat].phonon_thermal_conductivity;
         sc1d_G[cell_idx] += st::internal::mp[mat].electron_phonon_coupling;
         sc1d_T_Debye[cell_idx] += st::internal::mp[mat].einstein_temperature;
         num_atoms_in_cell[cell_idx]++;
      } else {
         // Material index out of range, just count atom
         num_atoms_in_cell[cell_idx]++;
      }
   }

   // Removed non-magnetic atoms seed thermal cell averages (Ce, Cp, kappa, G, T).
   // They are not mapped into sc1d_atom_cell_idx. Keep-nonmagnetic atoms are
   // already in the local atom list. Each rank adds the removed atoms it holds;
   // the Allreduce sums occupation and constants.
   for(size_t atom=0; atom<cs::non_magnetic_atoms_array.size(); ++atom){
      const cs::nm_atom_t& nm = cs::non_magnetic_atoms_array[atom];
      const int cell = static_cast<int>(std::floor(nm.z / micro_cell_thickness));
      const int cell_idx = std::max(0, std::min(cell, ncz - 1));
      const int mat = nm.mat;
      if(mat >= 0 && mat < (int)st::internal::mp.size()) {
         sc1d_Ce[cell_idx] += st::internal::mp[mat].electron_heat_capacity;
         sc1d_Cp[cell_idx] += st::internal::mp[mat].phonon_heat_capacity;
         sc1d_kappa_e[cell_idx] += st::internal::mp[mat].electron_thermal_conductivity;
         sc1d_kappa_p[cell_idx] += st::internal::mp[mat].phonon_thermal_conductivity;
         sc1d_G[cell_idx] += st::internal::mp[mat].electron_phonon_coupling;
         sc1d_T_Debye[cell_idx] += st::internal::mp[mat].einstein_temperature;
         num_atoms_in_cell[cell_idx]++;
      } else {
         num_atoms_in_cell[cell_idx]++;
      }
   }
   
   // Reduce across MPI ranks if needed
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &sc1d_Ce[0], ncz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &sc1d_Cp[0], ncz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &sc1d_kappa_e[0], ncz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &sc1d_kappa_p[0], ncz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &sc1d_G[0], ncz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &sc1d_T_Debye[0], ncz, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &num_atoms_in_cell[0], ncz, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif
   
   // Average material properties per cell
   for(int cell = 0; cell < ncz; ++cell) {
      if(num_atoms_in_cell[cell] > 0) {
         sc1d_Ce[cell] /= num_atoms_in_cell[cell];
         sc1d_Cp[cell] /= num_atoms_in_cell[cell];
         sc1d_kappa_e[cell] /= num_atoms_in_cell[cell];
         sc1d_kappa_p[cell] /= num_atoms_in_cell[cell];
         sc1d_G[cell] /= num_atoms_in_cell[cell];
         sc1d_T_Debye[cell] /= num_atoms_in_cell[cell];
      } else {
         // Empty cell: copy from previous cell
         if(cell > 0) {
            sc1d_Ce[cell] = sc1d_Ce[cell-1];
            sc1d_Cp[cell] = sc1d_Cp[cell-1];
            sc1d_kappa_e[cell] = sc1d_kappa_e[cell-1];
            sc1d_kappa_p[cell] = sc1d_kappa_p[cell-1];
            sc1d_G[cell] = sc1d_G[cell-1];
            sc1d_T_Debye[cell] = sc1d_T_Debye[cell-1];
         } else {
            terminaltextcolor(RED);
            std::cerr << "Error: Cell 0 has no atoms! Cannot initialize thermal gradients." << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
   }
   
   sc1d_thermal_gradients_initialised = true;
   zlog << zTs() << "Thermal gradients initialized successfully." << std::endl;
   
   // Output microcell configuration
   output_thermal_microcell_data();
   
   // Open temperature profile output file
   open_thermal_temperature_profile_file();
}

} // end of namespace internal
} // end of namespace st
