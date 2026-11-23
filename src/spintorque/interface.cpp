//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

// Vampire headers
#include "errors.hpp"
#include "spintorque.hpp"
#include "units.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "sim.hpp"


// Spin Torque headers
#include "internal.hpp"

namespace st{

// New material slots are value-initialised to 0. Conductivity uses negative
// as the unset sentinel (Einstein σ from De and Te). An explicit 0 is an insulator.
static void extend_spintorque_materials(const std::size_t n){
   if(n > 100) return;
   const std::size_t old = st::internal::mp.size();
   if(n <= old) return;
   st::internal::mp.resize(n);
   for(std::size_t i = old; i < n; ++i){
      st::internal::mp[i].conductivity = -1.0;
   }
}

// Directed pair, mirrored so cell k -> k+1 and the reverse share one value.
// Each keyword writes only its own fields.
static void store_sheet_pair(const int a, const int b,
                             const bool have_g, const double gup, const double gdn,
                             const bool have_mix, const double re, const double im,
                             const bool have_delta, const double delta){
   const std::size_t nmat = st::internal::mp.size();
   st::internal::ensure_sheet_matrix_size(nmat);
   const int ij[2][2] = {{a, b}, {b, a}};
   for(int n=0; n<2; ++n){
      const std::size_t idx = static_cast<std::size_t>(ij[n][0]) * nmat + static_cast<std::size_t>(ij[n][1]);
      if(have_g){
         st::internal::sheet_g_up_pair[idx] = gup;
         st::internal::sheet_g_down_pair[idx] = gdn;
      }
      if(have_mix){
         st::internal::sheet_re_pair[idx] = re;
         st::internal::sheet_im_pair[idx] = im;
      }
      if(have_delta) st::internal::sheet_delta_pair[idx] = delta;
      st::internal::sheet_pair_set[idx] = 1;
   }
   st::internal::sheet_enabled = true;
}

static bool parse_sheet_partner(const std::string& word, const std::string& value, const int line,
                                const int super_index, const int nparam, int& other_index,
                                std::vector<double>& params){
   params = vin::doubles_from_string(value);
   if((int)params.size() < nparam){
      terminaltextcolor(RED);
      std::cerr << "Error - expected " << nparam << " comma-separated values for " << word
                << " at line " << line << std::endl;
      terminaltextcolor(WHITE);
      err::vexit();
   }
   other_index = static_cast<int>(params[0]) - 1;
   if(other_index < 0 || other_index > 99){
      terminaltextcolor(RED);
      std::cerr << "Error - invalid partner material index for " << word << " at line " << line << std::endl;
      terminaltextcolor(WHITE);
      err::vexit();
   }
   const std::size_t want = static_cast<std::size_t>(std::max(super_index, other_index) + 1);
   if(want > st::internal::mp.size() && want < 101) extend_spintorque_materials(want);
   st::internal::ensure_sheet_matrix_size(st::internal::mp.size());
   return true;
}

bool match_material(string const word, string const value, string const unit, int const line, int const super_index){

   // add prefix string
   std::string prefix="material:";

   // Check for material id > current array size and if so dynamically expand mp array
   if((unsigned int) super_index + 1 > st::internal::mp.size() && super_index + 1 < 101){
      extend_spintorque_materials(static_cast<std::size_t>(super_index) + 1);
   }

   //STT constants------------------------------------------------------------
   std::string test="spin-diffusion-length"; 
   /*
      float spin-diffusion-length 
         Details
      */
   if(word==test){
      double lsdl=atof(value.c_str());
      vin::check_for_valid_value(lsdl, word, line, prefix, unit, "length", 0.01, 1.0e10,"material"," 0.01 - 1e10 Angstroms");
      st::internal::mp[super_index].lambda_sdl=lsdl*1.e-10; // defined in metres
      return true;
   }
   //--------------------------------------------------------------------
   test="spin-polarisation-conductivity"; //beta (eq 2)
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double betac=atof(value.c_str());
      vin::check_for_valid_value(betac, word, line, prefix, unit, "none", 0.0, 1.0,"material"," -1 - 1");
      st::internal::mp[super_index].beta_cond=betac;
      return true;
   }
   //--------------------------------------------------------------------
   test="spin-polarisation-diffusion"; //beta'
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double betad=atof(value.c_str());
      vin::check_for_valid_value(betad, word, line, prefix, unit, "none", 0.0, 1.0,"material"," -1 - 1");
      st::internal::mp[super_index].beta_diff=betad;
      return true;
   }
   //--------------------------------------------------------------------
   test="spin-accumulation"; //m(infinity)
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double sa=atof(value.c_str());
      vin::check_for_valid_value(sa, word, line, prefix, unit, "none", 0.0, 1.0e10,"material"," 1.0e-6 - 1.0e10");
      st::internal::mp[super_index].sa_infinity=sa;
      return true;
   }
   //--------------------------------------------------------------------
   test="diffusion-constant"; //D_0 (eq 2)
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double dc=atof(value.c_str());
      // Zero is an insulator: the face harmonic mean closes and no flux crosses.
      vin::check_for_valid_value(dc, word, line, prefix, unit, "none", 0.0, 100.0,"material"," 0 - 100"); //m^2/s
      st::internal::mp[super_index].diffusion=dc;
      return true;
   }
   //--------------------------------------------------------------------
   test="sd-exchange-constant";// 
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double sd=atof(value.c_str());
      vin::check_for_valid_value(sd, word, line, prefix, unit, "energy", 0.0, 1.e-17,"material"," 1.0e-30 - 1.e-17 J");
      st::internal::mp[super_index].sd_exchange=sd;
      return true;
   }

   //--------------------------------------------------------------------
   test="spin-dephasing-length"; // lambda_phi
   if(word==test){
      double lphi=atof(value.c_str());
      vin::check_for_valid_value(lphi, word, line, prefix, unit, "length", 0.0, 1.0e10,"material"," 0.0 - 1e10 Angstroms");
      st::internal::mp[super_index].lambda_phi=lphi*1.e-10;
      return true;
   }

   //--------------------------------------------------------------------
   test="demag-spin-coupling"; // chi (C/m^3): source -chi*d(|m|/msat)/dt on S
   if(word==test){
      double chi=atof(value.c_str());
      vin::check_for_valid_value(chi, word, line, prefix, unit, "none", 0.0, 1.0e30,"material"," 0.0 - 1e30");
      st::internal::mp[super_index].chi_demag=chi;
      return true;
   }

   //--------------------------------------------------------------------
   test="seebeck-coefficient"; // Seebeck coefficient S (V/K) for charge Seebeck effect
   if(word==test){
      double S=atof(value.c_str());
      vin::check_for_valid_value(S, word, line, prefix, unit, "none", -1.0e-3, 1.0e-3,"material"," -1e-3 - 1e-3 V/K");
      st::internal::mp[super_index].seebeck_coefficient=S;
      return true;
   }

   // P_S = (S_up - S_down) / S. beta_c_seebeck is computed at init unless beta-seebeck is set.
   test="seebeck-polarization";
   if(word==test){
      double ps=atof(value.c_str());
      vin::check_for_valid_value(ps, word, line, prefix, unit, "none", -5.0, 5.0,"material"," -5 - 5");
      st::internal::mp[super_index].seebeck_polarization=ps;
      return true;
   }

   test="beta-seebeck";
   if(word==test){
      double bs=atof(value.c_str());
      vin::check_for_valid_value(bs, word, line, prefix, unit, "none", -5.0, 5.0,"material"," -5 - 5");
      st::internal::mp[super_index].beta_seebeck=bs;
      st::internal::mp[super_index].beta_seebeck_set=true;
      return true;
   }

   // Polarization of the transported hot current, (j_up - j_down)/(j_up + j_down).
   // Same sign as beta_c: positive is majority. Default 0.
   test="superdiffusion-polarization";
   if(word==test){
      double ph=atof(value.c_str());
      vin::check_for_valid_value(ph, word, line, prefix, unit, "none", -5.0, 5.0,"material"," -5 - 5");
      st::internal::mp[super_index].superdiffusion_polarization=ph;
      return true;
   }

   test="electrical-conductivity";
   if(word==test){
      double sig=atof(value.c_str());
      vin::check_for_valid_value(sig, word, line, prefix, unit, "none", 0.0, 1.0e10,"material"," 0.0 - 1e10 S/m (0 = insulator; omit for Einstein from D,T)");
      st::internal::mp[super_index].conductivity=sig;
      return true;
   }

   //--------------------------------------------------------------------
   // Pair-wise interfacial spin conductance (Robin coupling)
   // Syntax: material[i]:spin-interface-conductance = j, <value>
   //    where j = partner material index (1-based), value = conductance in m/s
   // Example: material[1]:spin-interface-conductance = 2, 1.0e14
   // Units: m/s. Internally stored as resistance R_int = 1/G_int (s/m)
   //--------------------------------------------------------------------
   test="spin-interface-conductance";
   if(word == test){
      // Parse partner index and conductance value from comma-separated value string
      // Format: "partner_index, conductance_value"
      std::vector<double> params = vin::doubles_from_string(value);
      if(params.size() < 2){
         terminaltextcolor(RED);
         std::cerr << "Error - expected syntax 'material[i]:spin-interface-conductance = j, value' at line " << line << std::endl;
         std::cerr << "        where j is partner material index and value is conductance in m/s" << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      
      const int other_index = static_cast<int>(params[0]) - 1; // convert 1-based to 0-based
      const double gint = params[1];
      
      if(other_index < 0 || other_index > 99){
         terminaltextcolor(RED);
         std::cerr << "Error - invalid partner material index " << (other_index+1) << " for interface conductance at line " << line << std::endl;
         std::cerr << "        Partner index must be between 1 and 100" << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }

      // ensure material array sizes can accommodate both indices
      const std::size_t want = static_cast<std::size_t>(std::max(super_index, other_index) + 1);
      if(want > st::internal::mp.size() && want < 101) extend_spintorque_materials(want);
      st::internal::ensure_interface_matrix_size(st::internal::mp.size());

      if(gint < 0.0 || gint > 1.0e30){
         terminaltextcolor(RED);
         std::cerr << "Error - interface conductance value " << gint << " out of range (0 - 1e30 m/s) at line " << line << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      
      // store resistance; gint==0 treated as infinite resistance (disabled)
      if(gint > 0.0){
         const std::size_t nmat = st::internal::mp.size();
         const double rint = 1.0/gint;
         st::internal::r_int_pair[static_cast<std::size_t>(super_index)*nmat + static_cast<std::size_t>(other_index)] = rint;
         // default symmetric if reverse direction not yet set
         const std::size_t rev = static_cast<std::size_t>(other_index)*nmat + static_cast<std::size_t>(super_index);
         if(st::internal::r_int_pair[rev] == 0.0) st::internal::r_int_pair[rev] = rint;
         st::internal::interface_coupling_enabled = true;
      }
      return true;
   }
   //--------------------------------------------------------------------
   // Channel conductances in S/m^2. Unset pairs stay transparent.
   // Syntax: material[i]:spin-interface-conductance-up-down = j, G_up, G_down
   test="spin-interface-conductance-up-down";
   if(word == test){
      int other_index = 0;
      std::vector<double> params;
      parse_sheet_partner(word, value, line, super_index, 3, other_index, params);
      const double gup = params[1];
      const double gdn = params[2];
      if(gup < 0.0 || gdn < 0.0 || gup > 1.0e20 || gdn > 1.0e20){
         terminaltextcolor(RED);
         std::cerr << "Error - channel conductance out of range (0 - 1e20 S/m^2) at line " << line << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      store_sheet_pair(super_index, other_index, true, gup, gdn, false, 0.0, 0.0, false, 0.0);
      return true;
   }
   //--------------------------------------------------------------------
   // Tunnel conductance of an insulating segment or a zero-thickness face.
   // Not written into the bulk conductivity. Syntax: material[i]:tunnel-conductance = j, G_up, G_down
   test="tunnel-conductance";
   if(word == test){
      int other_index = 0;
      std::vector<double> params;
      parse_sheet_partner(word, value, line, super_index, 3, other_index, params);
      const double gup = params[1];
      const double gdn = params[2];
      if(gup < 0.0 || gdn < 0.0 || gup > 1.0e20 || gdn > 1.0e20){
         terminaltextcolor(RED);
         std::cerr << "Error - tunnel conductance out of range (0 - 1e20 S/m^2) at line " << line << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      const std::size_t nmat = st::internal::mp.size();
      st::internal::ensure_tunnel_matrix_size(nmat);
      const int ij[2][2] = {{super_index, other_index}, {other_index, super_index}};
      for(int n = 0; n < 2; ++n){
         const std::size_t idx = static_cast<std::size_t>(ij[n][0]) * nmat + static_cast<std::size_t>(ij[n][1]);
         st::internal::tunnel_g_up_pair[idx] = gup;
         st::internal::tunnel_g_down_pair[idx] = gdn;
         st::internal::tunnel_pair_set[idx] = 1;
      }
      st::internal::tunnel_enabled = true;
      return true;
   }
   //--------------------------------------------------------------------
   // Decay length of each tunnel channel, in the length unit (Angstroms).
   // Syntax: material[i]:tunnel-decay = j, lambda_up, lambda_down
   // The conductance used by the oxide is G_sigma * exp(-t / lambda_sigma),
   // with t the thickness of the insulating segment. A missing keyword leaves
   // that channel at its zero-thickness conductance.
   test="tunnel-decay";
   if(word == test){
      int other_index = 0;
      std::vector<double> params;
      parse_sheet_partner(word, value, line, super_index, 3, other_index, params);
      const double lup_A = params[1];
      const double ldn_A = params[2];
      if(lup_A < 0.01 || ldn_A < 0.01 || lup_A > 1.0e6 || ldn_A > 1.0e6){
         terminaltextcolor(RED);
         std::cerr << "Error - tunnel decay length out of range (0.01 - 1e6 Angstroms) at line " << line << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      const std::size_t nmat = st::internal::mp.size();
      st::internal::ensure_tunnel_matrix_size(nmat);
      const int ij[2][2] = {{super_index, other_index}, {other_index, super_index}};
      for(int n = 0; n < 2; ++n){
         const std::size_t idx = static_cast<std::size_t>(ij[n][0]) * nmat + static_cast<std::size_t>(ij[n][1]);
         st::internal::tunnel_l_up_pair[idx] = lup_A * 1.0e-10;
         st::internal::tunnel_l_down_pair[idx] = ldn_A * 1.0e-10;
      }
      return true;
   }
   //--------------------------------------------------------------------
   // Unpolarised inelastic conductance in parallel with the elastic channels.
   // Syntax: material[i]:tunnel-inelastic = j, G_e, G_p, lambda, p_e, p_p
   // G_e and G_p are the 300 K prefactors (S/m^2) of the Te and Tp pieces.
   // Each is multiplied by (T/300)^p and by exp(-t/lambda). Neither carries spin.
   test="tunnel-inelastic";
   if(word == test){
      int other_index = 0;
      std::vector<double> params;
      parse_sheet_partner(word, value, line, super_index, 6, other_index, params);
      const double ge = params[1];
      const double gp = params[2];
      const double lam_A = params[3];
      const double pe = params[4];
      const double pp = params[5];
      if(ge < 0.0 || gp < 0.0 || ge > 1.0e20 || gp > 1.0e20){
         terminaltextcolor(RED);
         std::cerr << "Error - inelastic tunnel conductance out of range at line " << line << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      if(lam_A < 0.01 || lam_A > 1.0e6){
         terminaltextcolor(RED);
         std::cerr << "Error - inelastic tunnel decay length out of range (0.01 - 1e6 Angstroms) at line " << line << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      if(pe < 0.0 || pp < 0.0 || pe > 10.0 || pp > 10.0){
         terminaltextcolor(RED);
         std::cerr << "Error - inelastic tunnel temperature power out of range (0 - 10) at line " << line << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      const std::size_t nmat = st::internal::mp.size();
      st::internal::ensure_tunnel_matrix_size(nmat);
      const int ij[2][2] = {{super_index, other_index}, {other_index, super_index}};
      for(int n = 0; n < 2; ++n){
         const std::size_t idx = static_cast<std::size_t>(ij[n][0]) * nmat + static_cast<std::size_t>(ij[n][1]);
         st::internal::tunnel_gi_e_pair[idx] = ge;
         st::internal::tunnel_gi_p_pair[idx] = gp;
         st::internal::tunnel_li_pair[idx] = lam_A * 1.0e-10;
         st::internal::tunnel_pe_pair[idx] = pe;
         st::internal::tunnel_pp_pair[idx] = pp;
      }
      return true;
   }
   //--------------------------------------------------------------------
   // Real and imaginary mixing conductance in S/m^2.
   // Syntax: material[i]:spin-mixing-conductance = j, Re, Im
   test="spin-mixing-conductance";
   if(word == test){
      int other_index = 0;
      std::vector<double> params;
      parse_sheet_partner(word, value, line, super_index, 3, other_index, params);
      const double re = params[1];
      const double im = params[2];
      if(re < -1.0e20 || re > 1.0e20 || im < -1.0e20 || im > 1.0e20){
         terminaltextcolor(RED);
         std::cerr << "Error - mixing conductance out of range at line " << line << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      store_sheet_pair(super_index, other_index, false, 0.0, 0.0, true, re, im, false, 0.0);
      return true;
   }
   //--------------------------------------------------------------------
   // Spin memory loss, dimensionless.
   // Syntax: material[i]:spin-memory-loss = j, delta
   test="spin-memory-loss";
   if(word == test){
      int other_index = 0;
      std::vector<double> params;
      parse_sheet_partner(word, value, line, super_index, 2, other_index, params);
      const double delta = params[1];
      if(delta < 0.0 || delta > 100.0){
         terminaltextcolor(RED);
         std::cerr << "Error - spin memory loss out of range (0 - 100) at line " << line << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      store_sheet_pair(super_index, other_index, false, 0.0, 0.0, false, 0.0, 0.0, true, delta);
      return true;
   }
   //--------------------------------------------------------------------
   test="spin-torque-free-layer"; //
   /*
    * Spin torque free layer flag
    */
   if(word==test){
      st::internal::free_layer = super_index;
      return true;
   }
   //--------------------------------------------------------------------
   test="spin-torque-reference-layer";
   /*
    * Spin torque free layer flag
    */
   if(word==test){
      st::internal::reference_layer = super_index;
      return true;
   }

   //SOT-SA constants------------------------------------------------------------
   test="sot-spin-diffusion-length"; 
   /*
      float spin-diffusion-length 
         Details
      */
   if(word==test){
      double lsdl=atof(value.c_str());
      vin::check_for_valid_value(lsdl, word, line, prefix, unit, "length", 0.01, 1.0e10,"material"," 0.01 - 1e10 Angstroms");
      st::internal::mp[super_index].sot_lambda_sdl=lsdl*1.e-10; // defined in metres
      return true;
   }
   //--------------------------------------------------------------------
   test="sot-spin-polarisation-conductivity"; //beta (eq 2)
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double betac=atof(value.c_str());
      vin::check_for_valid_value(betac, word, line, prefix, unit, "none", 1.0e-3, 1.0e3,"material"," 0.001 - 1000");
      st::internal::mp[super_index].sot_beta_cond=betac;
      return true;
   }
   //--------------------------------------------------------------------
   test="sot-spin-polarisation-diffusion"; //beta'
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double betad=atof(value.c_str());
      vin::check_for_valid_value(betad, word, line, prefix, unit, "none", 1.0e-9, 1.0e3,"material"," 0.001 - 1000");
      st::internal::mp[super_index].sot_beta_diff=betad;
      return true;
   }
   //--------------------------------------------------------------------
   test="sot-spin-accumulation"; //m(infinity)
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double sa=atof(value.c_str());
      vin::check_for_valid_value(sa, word, line, prefix, unit, "none", 0.0, 1.0e10,"material"," 1.0e-6 - 1.0e10");
      st::internal::mp[super_index].sot_sa_infinity=sa;
      return true;
   }
   //--------------------------------------------------------------------
   test="sot-diffusion-constant"; //D_0 (eq 2)
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double dc=atof(value.c_str());
      vin::check_for_valid_value(dc, word, line, prefix, unit, "none", 1.0e-9, 100.0,"material"," 1.0e-9 - 100"); //m^2/s
      st::internal::mp[super_index].sot_diffusion=dc;
      return true;
   }
   //--------------------------------------------------------------------
   test="sot-sd-exchange-constant";// 
   /*
      float spin-diffusion-length
         Details
      */
   if(word==test){
      double sd=atof(value.c_str());
      vin::check_for_valid_value(sd, word, line, prefix, unit, "energy", 0.0, 1.e-17,"material"," 1.0e-30 - 1.e-17 J");
      st::internal::mp[super_index].sot_sd_exchange=sd;
      return true;
   }
   //--------------------------------------------------------------------
   // Thermal gradient material parameters
   //--------------------------------------------------------------------
   test="electron-heat-capacity";
   if(word==test){
      double num=atof(value.c_str());
      vin::check_for_valid_value(num, word, line, prefix, unit, "J/(K^2*m^3)", 1.0e-3, 1.0e10,"material"," 0.001 - 1e10");
      st::internal::mp[super_index].electron_heat_capacity = num;
      return true;
   }

   test="phonon-heat-capacity";
   if(word==test){
      double num=atof(value.c_str());
      vin::check_for_valid_value(num, word, line, prefix, unit, "J/(K*m^3)", 1.0e-3, 1.0e10,"material"," 0.001 - 1e10");
      st::internal::mp[super_index].phonon_heat_capacity = num;
      return true;
   }

   test="electron-phonon-coupling";
   if(word==test){
      double num=atof(value.c_str());
      vin::check_for_valid_value(num, word, line, prefix, unit, "W/(m^3*K)", 0.0, 1.0e20,"material"," 0.0 - 1e20");
      st::internal::mp[super_index].electron_phonon_coupling = num;
      return true;
   }

   test="electron-thermal-conductivity";
   if(word==test){
      double num=atof(value.c_str());
      vin::check_for_valid_value(num, word, line, prefix, unit, "W/(m*K)", 0.0, 1.0e3,"material"," 0.0 - 1000");
      st::internal::mp[super_index].electron_thermal_conductivity = num;
      return true;
   }

   test="phonon-thermal-conductivity";
   if(word==test){
      double num=atof(value.c_str());
      vin::check_for_valid_value(num, word, line, prefix, unit, "W/(m*K)", 1.0e-3, 1.0e3,"material"," 0.001 - 1000");
      st::internal::mp[super_index].phonon_thermal_conductivity = num;
      return true;
   }

   test="einstein-temp";
   if(word==test){
      double num=atof(value.c_str());
      vin::check_for_valid_value(num, word, line, prefix, unit, "K", 0.0, 1.0e4,"material"," 0 - 10,000");
      // Store as Einstein temperature (ltmp multiplies by 1.25 internally, but we'll store as-is)
      st::internal::mp[super_index].einstein_temperature = num;
      return true;
   }

   //--------------------------------------------------------------------
   // keyword not found
   //--------------------------------------------------------------------
   return false;

}



   //-----------------------------------------------------------------------------
   // Function to process input file parameters for ST module
   //-----------------------------------------------------------------------------
   bool match_input_parameter(string const key, string const word, string const value, string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="spin-torque";


      if(key!=prefix) return false;

      //----------------------------------
      // Now test for all valid options
      //----------------------------------

      std::string test="enable-ST-fields";

       if(word==test){
         st::internal::enabled = true;
         return true;
       }

       //-------------------------------------------------

       test="TMRenable"; //set false

       if(word==test){
         st::internal::TMRenable = true;
         return true;
       }


      //-------------------------------------------------
      test="current-density"; //
       if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 1.0e15,"input","0.0 - 1.0e13 A/m2");
         st::internal::je =T;
         return true;
        }

      //-------------------------------------------------

       test="current-direction";
        if(word==test){
         std::string je_dir="x";
             if(value==je_dir){
                st::internal::current_direction=0;
                return true;
             }
             else
               je_dir ="y";
             if(value==je_dir){
              st::internal::current_direction=1;
              return true;
            }
            else{
              st::internal::current_direction=2;
              return true;
            }
        }

      //-------------------------------------------------

       test="micro-cell-size";
       if(word==test){
         std::vector<double> u(3);
         u=vin::doubles_from_string(value);
         //vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
         st::internal::micro_cell_size[0]=u.at(0);
         st::internal::micro_cell_size[1]=u.at(1);
         st::internal::micro_cell_size[2]=u.at(2);
         return true;
      }
      test = "micro-cell-decomp";
      if(word == test){
         std::string decomp_value = "mpi";
         if(decomp_value == value) {
            st::internal::microcell_decomp_type = decomp_value;
            return true;
         }
         decomp_value = "A";
         if(decomp_value == value) {
            st::internal::microcell_decomp_type = decomp_value;
            return true;
         }
         terminaltextcolor(RED);
            std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
            std::cerr << "\t\"mpi\"" << std::endl;
            std::cerr << "\t\"A\"" << std::endl;
         terminaltextcolor(WHITE);
            return false;

      }
      //-------------------------------------------------

      test="micro-cell-thickness";
       if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 20,"input","0.0 - 20 A");
         st::internal::micro_cell_thickness =T;
         return true;
        }

     //-------------------------------------------------

       test="initial-spin-polarisation"; //
       if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 20,"input","0.0 - 20 C/m^3");
         st::internal::initial_beta =T;
         return true;
        }

      test="spin-hall-angle"; //
       if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 10,"input","Spin Hall Angle: 0.0 - 1.0");
         st::internal::initial_theta =T;
         return true;
        }

      test="remove-material-type"; //
       if(word==test){
         double T=atof(value.c_str());
         vin::check_for_valid_value(T, word, line, prefix, unit, "none", 0.0, 20,"input","material type");
         st::internal::remove_nm =T;
         return true;
        }
    //--------------------------------------------------------------------
      test="initial-mag-direction";//might not matter
      if(word==test){
      std::vector<double> u(3);
      u=vin::doubles_from_string(value);
      vin::check_for_valid_unit_vector(u, word, line, prefix, "input");
      st::internal::initial_m[0]=u.at(0);
      st::internal::initial_m[1]=u.at(1);
      st::internal::initial_m[2]=u.at(2);
      return true;
   }
   test="finite-boundary-condition";
   if(word==test){
      if(value == "true"){
         st::internal::fbc = true;
      }  else if(value == "false"){
         st::internal::fbc = false;
      } else {
         terminaltextcolor(RED);
                std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
                std::cerr << "\t\"true\"" << std::endl;
                std::cerr << "\t\"false\"" << std::endl;
            terminaltextcolor(WHITE);
               return false;
      }

      return true;
   }

   test = "output-torque-data";
   if(word == test) {
      std::string output_type = "init";
      if(value == output_type){
         st::internal::output_torque_data = output_type;
         return true;
      }
      output_type = "int";
      if(value == output_type) {
         st::internal::output_torque_data = output_type;
         return true;
      }
      output_type = "final";
      if(value == output_type) {
         st::internal::output_torque_data = output_type;
         return true;
      }
      
      terminaltextcolor(RED);
      std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
      std::cerr << "\t\"init\"" << std::endl;
      std::cerr << "\t\"int\"" << std::endl;
      std::cerr << "\t\"final\"" << std::endl;
      terminaltextcolor(WHITE);
      return false;
      
   }

   test = "sot-check-only";
   if(word == test) {
      st::internal::sot_check = true;

      return true;
   }
      //-------------------------------------------------

   test="ST-output-rate";
   if(word==test){
      int T=atoi(value.c_str());
      vin::check_for_valid_int(T, word, line, prefix, 0, 2000000000,"input","0 - 2,000,000,000");
      st::internal::ST_output_rate =T;
      return true;
   }

   //-------------------------------------------------

   test="spin-currents-1d-enable";
   if(word==test){
      st::internal::sc1d_enable = true;
      st::internal::enabled = true;  // Automatically enable spin torque fields
      return true;
   }

   // Accumulation, currents, and output stay on. The atom spin-torque field is not applied.
   test="spin-currents-1d-calibrate";
   if(word==test){
      if(value == "true"){
         st::internal::sc1d_calibrate = true;
      } else if(value == "false"){
         st::internal::sc1d_calibrate = false;
      } else {
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"true\"" << std::endl;
         std::cerr << "\t\"false\"" << std::endl;
         terminaltextcolor(WHITE);
         return false;
      }
      return true;
   }

   // coarse: unit-cell bins, columns summed. atomic: one row per field window, columns summed.
   // verbose: every fine cell of every column.
   // The fine-grid constants file is written at initialisation for every mode.
   test="spin-currents-1d-output";
   if(word==test){
      if(value == "coarse"){
         st::internal::sc1d_output_mode = st::internal::sc1d_out_coarse;
      } else if(value == "atomic"){
         st::internal::sc1d_output_mode = st::internal::sc1d_out_atomic;
      } else if(value == "verbose"){
         st::internal::sc1d_output_mode = st::internal::sc1d_out_verbose;
      } else {
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"coarse\"" << std::endl;
         std::cerr << "\t\"atomic\"" << std::endl;
         std::cerr << "\t\"verbose\"" << std::endl;
         terminaltextcolor(WHITE);
         return false;
      }
      return true;
   }

   test="spin-currents-1d-fine-dz";
   if(word==test){
      double dz=atof(value.c_str());
      vin::check_for_valid_value(dz, word, line, prefix, unit, "length", 0.01, 1.0e10,"input"," 0.01 - 1e10 Angstroms");
      st::internal::sc1d_fine_dz = dz;
      return true;
   }

   test="spin-currents-1d-diffusion";
   if(word==test){
      if(value == "backward-euler" || value == "be"){
         st::internal::sc1d_diffusion = st::internal::sc1d_diff_be;
      } else if(value == "bdf2"){
         st::internal::sc1d_diffusion = st::internal::sc1d_diff_bdf2;
      } else {
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"backward-euler\"" << std::endl;
         std::cerr << "\t\"bdf2\"" << std::endl;
         terminaltextcolor(WHITE);
         return false;
      }
      return true;
   }

   test="spin-currents-1d-prolongation";
   if(word==test){
      if(value == "step"){
         st::internal::sc1d_prolongation = st::internal::sc1d_prolong_step;
      } else if(value == "linear"){
         st::internal::sc1d_prolongation = st::internal::sc1d_prolong_linear;
      } else if(value == "smoothstep"){
         st::internal::sc1d_prolongation = st::internal::sc1d_prolong_smooth;
      } else {
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"linear\"" << std::endl;
         std::cerr << "\t\"step\"" << std::endl;
         std::cerr << "\t\"smoothstep\"" << std::endl;
         terminaltextcolor(WHITE);
         return false;
      }
      return true;
   }

   // Width of the linear and smoothstep blend, in Angstroms. Step ignores it.
   // Omitted or non-positive values use the crystallographic unit-cell height.
   test="spin-currents-1d-prolongation-decay";
   if(word==test){
      double decay=atof(value.c_str());
      vin::check_for_valid_value(decay, word, line, prefix, unit, "length", 0.01, 1.0e10,"input"," 0.01 - 1e10 Angstroms");
      st::internal::sc1d_prolongation_decay = decay;
      return true;
   }

   // analytic: broaden the material-file interface with the prolongation kernel.
   // diffuse: average that same interface by the atomic populations, then prolong.
   test="spin-currents-1d-interface";
   if(word==test){
      if(value == "analytic"){
         st::internal::sc1d_interface = st::internal::sc1d_interface_analytic;
      } else if(value == "diffuse"){
         st::internal::sc1d_interface = st::internal::sc1d_interface_diffuse;
      } else {
         terminaltextcolor(RED);
         std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
         std::cerr << "\t\"analytic\"" << std::endl;
         std::cerr << "\t\"diffuse\"" << std::endl;
         terminaltextcolor(WHITE);
         return false;
      }
      return true;
   }

   test="spin-currents-1d-spin-stride";
   if(word==test){
      int stride = atoi(value.c_str());
      if(stride < 1) stride = 1;
      st::internal::sc1d_spin_stride = stride;
      return true;
   }

   test="spin-currents-1d-charge-stride";
   if(word==test){
      int stride = atoi(value.c_str());
      if(stride < 1) stride = 1;
      st::internal::sc1d_charge_stride = stride;
      return true;
   }

   test="spin-currents-1d-temperature";
   if(word==test){
      double Tk = atof(value.c_str());
      vin::check_for_valid_value(Tk, word, line, prefix, unit, "none", 0.0, 1.0e9,"input"," 0.0 - 1e9 K");
      st::internal::sc1d_temperature = Tk;
      return true;
   }

   test="spin-currents-1d-relax-steps";
   if(word==test){
      long steps_long = atol(value.c_str());
      double steps_double = static_cast<double>(steps_long);
      vin::check_for_valid_value(steps_double, word, line, prefix, unit, "none", 0.0, 1.0e9,"input"," 0 - 1e9 steps");
      if(steps_long < 0) steps_long = 0;
      st::internal::sc1d_relax_steps = static_cast<unsigned long>(steps_long);
      return true;
   }

test="spin-currents-1d-laser-enable";
if(word==test){
   st::internal::sc1d_laser_enable = true;
   return true;
}

test="spin-currents-1d-laser-Q0";
if(word==test){
   double Q0 = atof(value.c_str());
   vin::check_for_valid_value(Q0, word, line, prefix, unit, "none", 0.0, 1.0e30,"input"," 0.0 - 1e30 W/m^3");
   st::internal::sc1d_laser_Q0 = Q0;
   return true;
}

test="spin-currents-1d-optical-absorption-length";
if(word==test){
   double d = atof(value.c_str());
   vin::check_for_valid_value(d, word, line, prefix, unit, "none", 1.0e-12, 1.0e-3,"input"," 1e-12 - 1e-3 m");
   st::internal::sc1d_optical_absorption_length = d;
   return true;
}

test="spin-currents-1d-laser-wavelength";
if(word==test){
   double lam = atof(value.c_str());
   vin::check_for_valid_value(lam, word, line, prefix, unit, "none", 1.0e-9, 1.0e-3,"input"," 1e-9 - 1e-3 m");
   st::internal::sc1d_laser_wavelength = lam;
   return true;
}

test="spin-currents-1d-laser-eta";
if(word==test){
   double eta = atof(value.c_str());
   vin::check_for_valid_value(eta, word, line, prefix, unit, "none", 0.0, 100.0,"input"," 0.0 - 100.0");
   st::internal::sc1d_laser_eta = eta;
   return true;
}

test="spin-currents-1d-laser-t0";
if(word==test){
   double t0 = atof(value.c_str());
   vin::check_for_valid_value(t0, word, line, prefix, unit, "none", -1.0e3, 1.0e3,"input"," -1e3 - 1e3 s");
   st::internal::sc1d_laser_t0 = t0;
   return true;
}

test="spin-currents-1d-laser-fwhm";
if(word==test){
   double fwhm = atof(value.c_str());
   vin::check_for_valid_value(fwhm, word, line, prefix, unit, "none", 0.0, 1.0e3,"input"," 0.0 - 1e3 s");
   st::internal::sc1d_laser_fwhm = fwhm;
   return true;
}

test="spin-currents-1d-tau-s";
if(word==test){
   double tau = atof(value.c_str());
   vin::check_for_valid_value(tau, word, line, prefix, unit, "none", 0.0, 1.0e3,"input"," 0.0 - 1e3 s");
   st::internal::sc1d_tau_s = tau;
   return true;
}

test="spin-currents-1d-tau-demag";
if(word==test){
   double tau_d = atof(value.c_str());
   vin::check_for_valid_value(tau_d, word, line, prefix, unit, "none", 0.0, 1.0e3,"input"," 0.0 - 1e3 s");
   st::internal::sc1d_tau_demag = tau_d;
   return true;
}

//--------------------------------------------------------------------
// Superdiffusive transport parameters
//--------------------------------------------------------------------
test="spin-currents-1d-superdiffusive-enable";
if(word==test){
   if(value == "true"){
      st::internal::sc1d_superdiffusive_enable = true;
   } else if(value == "false"){
      st::internal::sc1d_superdiffusive_enable = false;
   } else {
      terminaltextcolor(RED);
      std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
      std::cerr << "\t\"true\"" << std::endl;
      std::cerr << "\t\"false\"" << std::endl;
      terminaltextcolor(WHITE);
      return false;
   }
   return true;
}

// gradient, a velocity (constant speed), or a time (decay lifetime).
// Velocity is stored in m/s. Vampire's velocity unit is angstrom/second,
// and m/s converts by 1e10, so divide that factor back out here.
// Time is already seconds after units::convert.
test="spin-currents-1d-superdiffusive-type";
if(word==test){
   if(value == "gradient"){
      st::internal::sc1d_superdiffusive_type = 0;
      return true;
   }
   if(unit.size() == 0){
      terminaltextcolor(RED);
      std::cerr << "Error - value for \'spin-torque:" << word << "\' must be \"gradient\", a velocity, or a time." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error - spin-torque:" << word << " needs \"gradient\", a velocity unit, or a time unit." << std::endl;
      err::vexit();
   }
   double num = std::atof(value.c_str());
   std::string utype;
   if(units::convert(unit, num, utype) == EXIT_FAILURE){
      terminaltextcolor(RED);
      std::cerr << "Error: Unit \'" << unit << "\' specified on line " << line << " of input file is not a valid unit." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error: Unit \'" << unit << "\' specified on line " << line << " of input file is not a valid unit." << std::endl;
      err::vexit();
   }
   if(utype == "velocity"){
      const double v_mps = num / 1.0e10;
      if(!std::isfinite(v_mps) || std::fabs(v_mps) > 1.0e8){
         terminaltextcolor(RED);
         std::cerr << "Error: spin-torque:" << word << " velocity on line " << line << " must be within 1e8 m/s." << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Error: spin-torque:" << word << " velocity on line " << line << " must be within 1e8 m/s." << std::endl;
         err::vexit();
      }
      // The number is a speed. signed_hot_electron_velocity points it toward the substrate.
      st::internal::sc1d_superdiffusive_type = 1;
      st::internal::sc1d_superdiffusive_speed = v_mps;
      return true;
   }
   if(utype == "time"){
      if(!std::isfinite(num) || num < 1.0e-18 || num > 1.0){
         terminaltextcolor(RED);
         std::cerr << "Error: spin-torque:" << word << " lifetime on line " << line << " must be in the range 1e-18 - 1 s." << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Error: spin-torque:" << word << " lifetime on line " << line << " must be in the range 1e-18 - 1 s." << std::endl;
         err::vexit();
      }
      st::internal::sc1d_superdiffusive_type = 2;
      st::internal::sc1d_superdiffusive_lifetime = num;
      return true;
   }
   terminaltextcolor(RED);
   std::cerr << "Error: Unit \'" << unit << "\' of type \'" << utype << "\' on line " << line
             << " is invalid for spin-torque:" << word << ". Use a velocity or a time." << std::endl;
   terminaltextcolor(WHITE);
   zlog << zTs() << "Error: Unit \'" << unit << "\' of type \'" << utype << "\' on line " << line
        << " is invalid for spin-torque:" << word << "." << std::endl;
   err::vexit();
}

test="spin-currents-1d-tau-e";
if(word==test){
   double tau_e = atof(value.c_str());
   vin::check_for_valid_value(tau_e, word, line, prefix, unit, "none", 1.0e-18, 1.0e-9,"input"," 1e-18 - 1e-9 s (typically 50-200 fs)");
   st::internal::sc1d_tau_e = tau_e;
   return true;
}

test="spin-currents-1d-v-e0";
if(word==test){
   double v_e0 = atof(value.c_str());
   vin::check_for_valid_value(v_e0, word, line, prefix, unit, "none", 0.0, 1.0e7,"input"," 0.0 - 1e7 m/s (0 = auto from De/dopt)");
   st::internal::sc1d_v_e0 = v_e0;
   return true;
}

//--------------------------------------------------------------------
// ltmp integration flags for temperature-dependent scaling
//--------------------------------------------------------------------
// Electron and phonon temperatures are read from ltmp. The spin-current
// module does not step its own temperature profile.

test="spin-currents-1d-use-ltmp-temperatures";
if(word==test){
   if(value == "true"){
      st::internal::sc1d_use_ltmp_temperatures = true;
   } else if(value == "false"){
      st::internal::sc1d_use_ltmp_temperatures = false;
   } else {
      terminaltextcolor(RED);
      std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
      std::cerr << "\t\"true\"" << std::endl;
      std::cerr << "\t\"false\"" << std::endl;
      terminaltextcolor(WHITE);
      return false;
   }
   return true;
}

test="spin-currents-1d-reference-temperature";
if(word==test){
   double Tref = atof(value.c_str());
   vin::check_for_valid_value(Tref, word, line, prefix, unit, "none", 0.0, 1.0e6,"input"," 0.0 - 1e6 K");
   st::internal::sc1d_reference_temperature = Tref;
   return true;
}

// Transport constants (D, lambda, tau, S, sigma) use Te and Tp only when this
// is true. The temperature laws themselves are not filled in yet.
test="spin-currents-1d-thermal-effects";
if(word==test){
   if(value == "true"){
      st::internal::sc1d_thermal_effects = true;
   } else if(value == "false"){
      st::internal::sc1d_thermal_effects = false;
   } else {
      terminaltextcolor(RED);
      std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
      std::cerr << "\t\"true\"" << std::endl;
      std::cerr << "\t\"false\"" << std::endl;
      terminaltextcolor(WHITE);
      return false;
   }
   return true;
}

//--------------------------------------------------------------------
// Seebeck effect parameters
//--------------------------------------------------------------------
test="spin-currents-1d-bias";
if(word==test){
   if(value == "true"){
      st::internal::sc1d_bias_enable = true;
   } else if(value == "false"){
      st::internal::sc1d_bias_enable = false;
   } else {
      terminaltextcolor(RED);
      std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
      std::cerr << "\t\"true\"" << std::endl;
      std::cerr << "\t\"false\"" << std::endl;
      terminaltextcolor(WHITE);
      return false;
   }
   return true;
}

test="spin-currents-1d-voltage-bottom";
if(word==test){
   double volts = atof(value.c_str());
   vin::check_for_valid_value(volts, word, line, prefix, unit, "none", -1.0e3, 1.0e3,"input"," -1000 - 1000 V");
   st::internal::sc1d_voltage_bottom = volts;
   return true;
}

test="spin-currents-1d-voltage-top";
if(word==test){
   double volts = atof(value.c_str());
   vin::check_for_valid_value(volts, word, line, prefix, unit, "none", -1.0e3, 1.0e3,"input"," -1000 - 1000 V");
   st::internal::sc1d_voltage_top = volts;
   return true;
}

test="spin-currents-1d-seebeck-enable";
if(word==test){
   if(value == "true"){
      st::internal::sc1d_seebeck_enable = true;
   } else if(value == "false"){
      st::internal::sc1d_seebeck_enable = false;
   } else {
      terminaltextcolor(RED);
      std::cerr << "Error - value for \'spin-torque:" << word << "\' must be one of:" << std::endl;
      std::cerr << "\t\"true\"" << std::endl;
      std::cerr << "\t\"false\"" << std::endl;
      terminaltextcolor(WHITE);
      return false;
   }
   return true;
}

   test="SOT-spin-accumulation";
   if(word==test){
      
      st::internal::sot_sa =true;
      return true;
   }

      //--------------------------------------------------------------------
      // input parameter not found here
      return false;
   }

} // end of namespace st
