//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>

// Vampire headers
#include "spintorque.hpp"
#include "atoms.hpp"
#include "create.hpp"
#include "material.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace {

static constexpr double kMuB = 9.27400968e-24;
static constexpr double kE = 1.60217662e-19;
static constexpr double kHbar = 1.05457162e-34;

// Same atomic volume as the laser-electrical Eq. 6 conversion: cell in Å³.
static double eq6_atom_volume(){
   const double natom = std::max(1.0, static_cast<double>(cs::unit_cell.atom.size()));
   return cs::unit_cell.dimensions[0] * cs::unit_cell.dimensions[1] * cs::unit_cell.dimensions[2]
          / natom * 1.0e-30;
}

// Serban Eq. (6) as a field the LLG can cross, same algebra as calculate_full_spin_fields.
// T_S = -(De/λ_J²) m×S − (De/λ_φ²) m×(m×S), added as (V/μ_S) T_S.
// De/λ_J² is J_sd/(2ℏ), the precession rate already used for S.
// The solved S is a charge density; (μ_B/e) S is the moment density in Eq. (6).
// H = (r − α p)(m × ŝ) + (p + α r) ŝ, with p and r the field-like and damping-like Tesla amplitudes.
static void apply_sc1d_atom_spin_torque_fields(const std::vector<double>& x_spin_array,
                                               const std::vector<double>& y_spin_array,
                                               const std::vector<double>& z_spin_array,
                                               const std::vector<int>& atom_type_array,
                                               const std::vector<double>& mu_s_array){
   (void)mu_s_array;
   const int natoms = st::internal::num_local_atoms;
   const int nf = st::internal::sc1d_nf;
   const bool in_equilibration = (st::internal::sc1d_step_counter <= st::internal::sc1d_relax_steps);
   const double V_atom = eq6_atom_volume();
   const double muB_over_e = kMuB / kE;

   for(int atom=0; atom<natoms; ++atom){
      st::internal::x_field_array[atom] = 0.0;
      st::internal::y_field_array[atom] = 0.0;
      st::internal::z_field_array[atom] = 0.0;

      if(nf <= 0 || st::internal::sc1d_Sfine.empty()) continue;
      if(atom >= (int)st::internal::sc1d_atom_ls.size()) continue;

      const int ls = st::internal::sc1d_atom_ls[atom];
      const int i_lo = st::internal::sc1d_atom_fine_lo[atom];
      const int i_hi = st::internal::sc1d_atom_fine_hi[atom];
      if(ls < 0 || i_lo < 0 || i_hi < i_lo || i_hi >= nf) continue;

      const int mat = atom_type_array[atom];
      if(mat < 0 || mat >= (int)st::internal::mp.size()) continue;
      const double Jsd = st::internal::mp[mat].sd_exchange;
      if(Jsd <= 0.0) continue;
      if(in_equilibration) continue;

      const double* Sbase = st::internal::sc1d_Sfine.data() + (std::size_t)ls * (std::size_t)nf * 3u;

      double sx = 0.0;
      double sy = 0.0;
      double sz = 0.0;
      int ncell = 0;
      for(int i=i_lo; i<=i_hi; ++i){
         sx += Sbase[3*i+0];
         sy += Sbase[3*i+1];
         sz += Sbase[3*i+2];
         ncell++;
      }
      if(ncell <= 0) continue;
      const double inv = 1.0 / (double)ncell;
      sx *= inv * muB_over_e;
      sy *= inv * muB_over_e;
      sz *= inv * muB_over_e;

      const double mu = mp::material[mat].mu_s_SI;
      const double gamma = mp::gamma_SI * mp::material[mat].gamma_rel;
      if(!(mu > 0.0) || !(gamma > 0.0) || !(V_atom > 0.0)) continue;

      const double De = st::internal::mp[mat].diffusion;
      const double lphi = st::internal::mp[mat].lambda_phi;
      const double c_field = Jsd / (2.0 * kHbar);
      const double c_damp = (lphi > 0.0 && De > 0.0) ? (De / (lphi * lphi)) : 0.0;
      const double pref = V_atom / (gamma * mu);
      const double alpha = mp::material[mat].alpha;
      const double along = pref * (c_field + alpha * c_damp);
      const double cross = pref * (c_damp - alpha * c_field);

      const double mx = x_spin_array[atom];
      const double my = y_spin_array[atom];
      const double mz = z_spin_array[atom];
      const double cx = my * sz - mz * sy;
      const double cy = mz * sx - mx * sz;
      const double cz = mx * sy - my * sx;

      st::internal::x_field_array[atom] = along * sx + cross * cx;
      st::internal::y_field_array[atom] = along * sy + cross * cy;
      st::internal::z_field_array[atom] = along * sz + cross * cz;
   }
}

}

namespace st{


   //-----------------------------------------------------------------------------
   // Function for updating spin torque fields
   //-----------------------------------------------------------------------------
   void update_spin_torque_fields(const std::vector<double>& x_spin_array,
                                  const std::vector<double>& y_spin_array,
                                  const std::vector<double>& z_spin_array,
                                  const std::vector<int>& atom_type_array,
                                  const std::vector<double>& mu_s_array){
       
      if(st::internal::enabled==false) return;


      // update magnetisations
      st::internal::update_cell_magnetisation(x_spin_array, y_spin_array, z_spin_array, atom_type_array, mu_s_array);

      // calculate spin_accumulation
      if(st::internal::sot_sa) st::internal::calculate_sot_accumulation();
      else {
         if(st::internal::sc1d_enable) st::internal::calculate_spin_accumulation_1d();
         else st::internal::calculate_spin_accumulation();
      }

      if(st::internal::sc1d_enable){
         apply_sc1d_atom_spin_torque_fields(x_spin_array, y_spin_array, z_spin_array,
                                            atom_type_array, mu_s_array);
      } else {
         for(int atom=0; atom<st::internal::num_local_atoms; ++atom) {
            const int cell3 = 3*st::internal::atom_st_index[atom];
            const double i_mu_s = 1.0/(mu_s_array[atom_type_array[atom]]);
            st::internal::x_field_array[atom] = st::internal::spin_torque[cell3+0]*i_mu_s;
            st::internal::y_field_array[atom] = st::internal::spin_torque[cell3+1]*i_mu_s;
            st::internal::z_field_array[atom] = st::internal::spin_torque[cell3+2]*i_mu_s;
         }
      }
   }

   //-----------------------------------------------------------------------------
   // Function for adding atomic spin torque fields to external field array
   //-----------------------------------------------------------------------------
   void get_spin_torque_fields(std::vector<double>& x_total_external_field_array,
                               std::vector<double>& y_total_external_field_array,
                               std::vector<double>& z_total_external_field_array,
                               const int start_index,
                               const int end_index){
     
      if(st::internal::enabled==false) return;


      // Add spin torque fields
      for(int i=start_index; i<end_index; ++i) x_total_external_field_array[i] += st::internal::x_field_array[i];
      for(int i=start_index; i<end_index; ++i) y_total_external_field_array[i] += st::internal::y_field_array[i];
      for(int i=start_index; i<end_index; ++i) z_total_external_field_array[i] += st::internal::z_field_array[i];

      return;
   }

   //-----------------------------------------------------------------------------
   // Spin-current equilibration, before ASD equilibration starts the LLG.
   // Magnetisation is not integrated. Each call steps the 1D solver once;
   // while the counter is inside spin-currents-1d-relax-steps the
   // charge current stays pinned at zero and the torque field stays zero.
   //-----------------------------------------------------------------------------
   void equilibrate_spin_currents_1d(){
      if(!internal::sc1d_enable || internal::sc1d_relax_steps == 0UL) return;
      static bool done = false;
      if(done) return;
      done = true;

      for(unsigned long n = 0; n < internal::sc1d_relax_steps; ++n){
         update_spin_torque_fields(atoms::x_spin_array,
                                   atoms::y_spin_array,
                                   atoms::z_spin_array,
                                   atoms::type_array,
                                   mp::mu_s_array);
      }
   }

} // end of st namespace
