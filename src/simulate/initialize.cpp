//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2020. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <iostream>

// Vampire headers
#include "constants.hpp"
#include "create.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "vio.hpp"
#include "internal.hpp"

namespace sim{

   //-------------------------------------------------------------------------------
   // initialise sim namespace variables
   //-------------------------------------------------------------------------------
   void initialize(int num_materials){

      // unroll slonczewski spin transfer torque arrays
      sim::internal::stt_asm.resize(num_materials,0.0);
      sim::internal::stt_rj.resize(num_materials,0.0);
      sim::internal::stt_pj.resize(num_materials,0.0);

      // unroll spin orbit torque arrays
      sim::internal::sot_asm.resize(num_materials,0.0);
      sim::internal::sot_asm_2nd_order.resize(num_materials,0.0);
      sim::internal::sot_rj.resize(num_materials,0.0);
      sim::internal::sot_pj.resize(num_materials,0.0);

      sim::internal::sot_asm2.resize(num_materials,0.0);
      sim::internal::sot_rj2.resize(num_materials,0.0);
      sim::internal::sot_pj2.resize(num_materials,0.0);

      sim::internal::let_stt_fl.resize(num_materials,0.0);
      sim::internal::let_stt_dl.resize(num_materials,0.0);
      sim::internal::let_fl_torkance.resize(num_materials,0.0);
      sim::internal::let_dl_torkance.resize(num_materials,0.0);
      sim::internal::let_eq6_fl.resize(num_materials,0.0);
      sim::internal::let_eq6_dl.resize(num_materials,0.0);
      sim::internal::let_fl_mode.resize(num_materials,0);
      sim::internal::let_dl_mode.resize(num_materials,0);
      sim::internal::let_P.resize(num_materials,1.0);
      sim::internal::let_px.resize(num_materials,0.0);
      sim::internal::let_py.resize(num_materials,0.0);
      sim::internal::let_pz.resize(num_materials,0.0);

      sim::internal::lot_lt_x.resize(num_materials, 0.0);
      sim::internal::lot_lt_y.resize(num_materials, 0.0);
      sim::internal::lot_lt_z.resize(num_materials, 0.0);

      sim::internal::vcmak.resize(num_materials);

      sim::STDspin_parallel_initialized = false;
      sim::c_octants.resize(8);
      sim::b_octants.resize(8);
      // loop over materials set by user
      for(unsigned int m=0; m < sim::internal::mp.size(); ++m){
         // copy values set by user to arrays
         if(sim::internal::mp[m].stt_asm.is_set()) sim::internal::stt_asm[m] = sim::internal::mp[m].stt_asm.get();
         if(sim::internal::mp[m].stt_rj.is_set())  sim::internal::stt_rj[m]  = sim::internal::mp[m].stt_rj.get();
         if(sim::internal::mp[m].stt_pj.is_set())  sim::internal::stt_pj[m]  = sim::internal::mp[m].stt_pj.get();

         if(sim::internal::mp[m].sot_asm.is_set()) sim::internal::sot_asm[m] = sim::internal::mp[m].sot_asm.get();
         if(sim::internal::mp[m].sot_asm_2nd_order.is_set()) sim::internal::sot_asm_2nd_order[m] = sim::internal::mp[m].sot_asm_2nd_order.get();
         if(sim::internal::mp[m].sot_rj.is_set())  sim::internal::sot_rj[m]  = sim::internal::mp[m].sot_rj.get();
         if(sim::internal::mp[m].sot_pj.is_set())  sim::internal::sot_pj[m]  = sim::internal::mp[m].sot_pj.get();

         if(sim::internal::mp[m].sot_asm2.is_set()) sim::internal::sot_asm2[m] = sim::internal::mp[m].sot_asm2.get();
         if(sim::internal::mp[m].sot_rj2.is_set())  sim::internal::sot_rj2[m]  = sim::internal::mp[m].sot_rj2.get();
         if(sim::internal::mp[m].sot_pj2.is_set())  sim::internal::sot_pj2[m]  = sim::internal::mp[m].sot_pj2.get();

         if(sim::internal::mp[m].let_stt_fl.is_set()) sim::internal::let_stt_fl[m] = sim::internal::mp[m].let_stt_fl.get();
         if(sim::internal::mp[m].let_stt_dl.is_set()) sim::internal::let_stt_dl[m] = sim::internal::mp[m].let_stt_dl.get();
         if(sim::internal::mp[m].let_fl_torkance.is_set()) sim::internal::let_fl_torkance[m] = sim::internal::mp[m].let_fl_torkance.get();
         if(sim::internal::mp[m].let_dl_torkance.is_set()) sim::internal::let_dl_torkance[m] = sim::internal::mp[m].let_dl_torkance.get();
         if(sim::internal::mp[m].let_P.is_set()) sim::internal::let_P[m] = sim::internal::mp[m].let_P.get();

         if(m < sim::internal::let_px.size() && sim::internal::mp[m].let_pvec_set){
            sim::internal::let_px[m] = sim::internal::mp[m].let_pvec[0];
            sim::internal::let_py[m] = sim::internal::mp[m].let_pvec[1];
            sim::internal::let_pz[m] = sim::internal::mp[m].let_pvec[2];
         }

         // Serban Eq. 6: T_S = -(De/λJ²) m×S − (De/λφ²) m×(m×S). LLG wants Tesla.
         // Js = P (μB/e) Jc. Eq. 4 De∇S ~ Js over the decay length k λ (default k=1):
         // S_FL = (k λJ / De) Js, S_DL = (k λφ / De) Js. Then H = [V_atom/(γ μ_S)] (De/λ²) S.
         const bool have_de = sim::internal::mp[m].let_diffusion.is_set();
         const bool have_lj = sim::internal::mp[m].let_lambda_J.is_set();
         const bool have_lp = sim::internal::mp[m].let_lambda_phi.is_set();
         if(have_de && (have_lj || have_lp)){
            const double natom = std::max<double>(1.0, static_cast<double>(cs::unit_cell.atom.size()));
            const double V_atom = cs::unit_cell.dimensions[0]*cs::unit_cell.dimensions[1]*cs::unit_cell.dimensions[2] / natom * 1.0e-30;
            const double mu = mp::material[m].mu_s_SI;
            const double g = mp::gamma_SI * mp::material[m].gamma_rel;
            const double pref = V_atom / (g * mu);
            const double De = sim::internal::mp[m].let_diffusion.get();
            const double kabs = sim::internal::mp[m].let_abs_factor.is_set() ? sim::internal::mp[m].let_abs_factor.get() : 1.0;
            if(have_lj){
               const double lJ = sim::internal::mp[m].let_lambda_J.get();
               const double S_per_Js = kabs * lJ / De;
               sim::internal::let_eq6_fl[m] = pref * (De / (lJ * lJ)) * S_per_Js;
            }
            if(have_lp){
               const double lphi = sim::internal::mp[m].let_lambda_phi.get();
               const double S_per_Js = kabs * lphi / De;
               sim::internal::let_eq6_dl[m] = pref * (De / (lphi * lphi)) * S_per_Js;
            }
            zlog << zTs() << "LET STT Eq. 6 for material " << m+1
                 << ": P " << sim::internal::let_P[m]
                 << " k " << kabs
                 << " p " << sim::internal::let_px[m] << " "
                 << sim::internal::let_py[m] << " " << sim::internal::let_pz[m]
                 << " FL " << sim::internal::let_eq6_fl[m]
                 << " DL " << sim::internal::let_eq6_dl[m] << " T s / A" << std::endl;
            const double p2 = sim::internal::let_px[m]*sim::internal::let_px[m]
                            + sim::internal::let_py[m]*sim::internal::let_py[m]
                            + sim::internal::let_pz[m]*sim::internal::let_pz[m];
            if(p2 < 1.0e-18){
               zlog << zTs() << "Warning: LET polarisation for material " << m+1
                    << " is zero; STT field will be zero." << std::endl;
            }
         }

         // Priority: torkance > Eq. 6 De/λ > Tesla at Jc_ref
         if(sim::internal::mp[m].let_fl_torkance.is_set()) sim::internal::let_fl_mode[m] = 1;
         else if(have_de && have_lj) sim::internal::let_fl_mode[m] = 2;
         else sim::internal::let_fl_mode[m] = 0;

         if(sim::internal::mp[m].let_dl_torkance.is_set()) sim::internal::let_dl_mode[m] = 1;
         else if(have_de && have_lp) sim::internal::let_dl_mode[m] = 2;
         else sim::internal::let_dl_mode[m] = 0;

         if(sim::internal::let_fl_mode[m] != 0 || sim::internal::let_dl_mode[m] != 0)
            sim::internal::enable_spin_torque_fields = true;

         if(sim::internal::mp[m].lt_x.is_set())  sim::internal::lot_lt_x[m] = sim::internal::mp[m].lt_x.get();
         if(sim::internal::mp[m].lt_y.is_set())  sim::internal::lot_lt_y[m] = sim::internal::mp[m].lt_y.get();
         if(sim::internal::mp[m].lt_z.is_set())  sim::internal::lot_lt_z[m] = sim::internal::mp[m].lt_z.get();

         // set vcma coefficients (requires sim::internal::enable_vcma_fields == true) but this should be default
         if(sim::internal::mp[m].vcmak.is_set()){
            const double imu_s = 1.0 / mp::material[m].mu_s_SI; // calculate inverse moment
            sim::internal::vcmak[m] = imu_s * sim::internal::mp[m].vcmak.get();
         }
      }

      return;
   }

} // end of namespace gpu
