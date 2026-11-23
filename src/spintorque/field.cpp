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
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Vampire headers
#include "spintorque.hpp"
#include "atoms.hpp"
#include "create.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

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

// One atom that receives Eq. (6). Window, scale, and the two Tesla coefficients
// are fixed after the fine grid is built. S and m are the step quantities.
struct TorqueTerm {
   int atom;
   int ncell;
   std::size_t s_off;
   double s_scale;
   double along;
   double cross;
};

std::vector<TorqueTerm> torque_terms;

// Serban Eq. (6) as a field the LLG can cross, same algebra as calculate_full_spin_fields.
// T_S = -(De/λ_J²) m×S − (De/λ_φ²) m×(m×S), added as (V/μ_S) T_S.
// De/λ_J² is J_sd/(2ℏ), the precession rate already used for S.
// The solved S is a charge density; (μ_B/e) S is the moment density in Eq. (6).
// H = (r − α p)(m × ŝ) + (p + α r) ŝ, with p and r the field-like and damping-like Tesla amplitudes.
// Relax holds m fixed and does not store H. A few samples are still evaluated
// for sc1d-heff-tesla.log; those rows are marked stored=0.
static int pick_named_torque_atom(const std::string& name){
   double sx = 0.0;
   double sy = 0.0;
   int n = 0;
   for(std::size_t t = 0; t < torque_terms.size(); ++t){
      const int atom = torque_terms[t].atom;
      const int mat = atoms::type_array[atom];
      if(mat < 0 || mat >= static_cast<int>(mp::material.size())) continue;
      if(mp::material[static_cast<std::size_t>(mat)].name != name) continue;
      sx += atoms::x_coord_array[atom];
      sy += atoms::y_coord_array[atom];
      ++n;
   }
   if(n == 0) return -1;
   const double cx = sx / static_cast<double>(n);
   const double cy = sy / static_cast<double>(n);
   int best = -1;
   double best_d = 0.0;
   for(std::size_t t = 0; t < torque_terms.size(); ++t){
      const int atom = torque_terms[t].atom;
      const int mat = atoms::type_array[atom];
      if(mat < 0 || mat >= static_cast<int>(mp::material.size())) continue;
      if(mp::material[static_cast<std::size_t>(mat)].name != name) continue;
      const double dx = atoms::x_coord_array[atom] - cx;
      const double dy = atoms::y_coord_array[atom] - cy;
      const double d = dx * dx + dy * dy;
      if(best < 0 || d < best_d){
         best = atom;
         best_d = d;
      }
   }
   return best;
}

static bool heff_sample_due(const unsigned long counter){
   const unsigned long nrelax = st::internal::sc1d_relax_steps;
   if(counter == 1UL || counter == nrelax || counter == nrelax + 1UL) return true;
   if(counter <= nrelax) return false;
   if(!(mp::dt_SI > 0.0)) return false;
   const double t = mp::dt_SI * static_cast<double>(sim::time - sim::equilibration_time);
   const double half = 0.51 * mp::dt_SI;
   const double targets[4] = {0.34e-12, 0.68e-12, 2.34e-12, 2.68e-12};
   for(int i = 0; i < 4; ++i){
      if(std::fabs(t - targets[i]) <= half) return true;
   }
   return false;
}

static void apply_sc1d_atom_spin_torque_fields(const std::vector<double>& x_spin_array,
                                               const std::vector<double>& y_spin_array,
                                               const std::vector<double>& z_spin_array,
                                               const std::vector<int>& atom_type_array,
                                               const std::vector<double>& mu_s_array){
   (void)atom_type_array;
   (void)mu_s_array;
   // Relax holds m fixed and the torque at zero. The field arrays are already
   // zero from initialisation, and this function is the only writer.
   const bool store_field = st::internal::sc1d_step_counter > st::internal::sc1d_relax_steps;
   const bool log_now = heff_sample_due(st::internal::sc1d_step_counter);
   if(!store_field && !log_now) return;
   static int sample_co = -2;
   static int sample_cfb = -2;
   if(sample_co == -2){
      sample_co = pick_named_torque_atom("Co-5");
      sample_cfb = pick_named_torque_atom("CoFeB");
   }

   const double* Sfine = st::internal::sc1d_Sfine.data();
   const double* mx = x_spin_array.data();
   const double* my = y_spin_array.data();
   const double* mz = z_spin_array.data();
   double* Hx = st::internal::x_field_array.data();
   double* Hy = st::internal::y_field_array.data();
   double* Hz = st::internal::z_field_array.data();
   const TorqueTerm* terms = torque_terms.data();
   const std::size_t nterm = torque_terms.size();

   struct Held {
      bool hit;
      double fx, fy, fz;
      double sx, sy, sz;
      double px, py, pz;
   };
   Held held_co = {false, 0, 0, 0, 0, 0, 0, 0, 0, 0};
   Held held_cfb = {false, 0, 0, 0, 0, 0, 0, 0, 0, 0};

   for(std::size_t t = 0; t < nterm; ++t){
      const TorqueTerm& term = terms[t];
      const double* S = Sfine + term.s_off;
      double sx = 0.0;
      double sy = 0.0;
      double sz = 0.0;
      const int ncell = term.ncell;
      for(int i = 0; i < ncell; ++i){
         sx += S[3 * i];
         sy += S[3 * i + 1];
         sz += S[3 * i + 2];
      }
      const double scale = term.s_scale;
      sx *= scale;
      sy *= scale;
      sz *= scale;

      const int atom = term.atom;
      const double px = mx[atom];
      const double py = my[atom];
      const double pz = mz[atom];
      const double cx = py * sz - pz * sy;
      const double cy = pz * sx - px * sz;
      const double cz = px * sy - py * sx;
      const double fx = term.along * sx + term.cross * cx;
      const double fy = term.along * sy + term.cross * cy;
      const double fz = term.along * sz + term.cross * cz;
      if(store_field){
         Hx[atom] = fx;
         Hy[atom] = fy;
         Hz[atom] = fz;
      }
      if(log_now && atom == sample_co){
         held_co.hit = true;
         held_co.fx = fx; held_co.fy = fy; held_co.fz = fz;
         held_co.sx = sx; held_co.sy = sy; held_co.sz = sz;
         held_co.px = px; held_co.py = py; held_co.pz = pz;
      } else if(log_now && atom == sample_cfb){
         held_cfb.hit = true;
         held_cfb.fx = fx; held_cfb.fy = fy; held_cfb.fz = fz;
         held_cfb.sx = sx; held_cfb.sy = sy; held_cfb.sz = sz;
         held_cfb.px = px; held_cfb.py = py; held_cfb.pz = pz;
      }
   }

   if(!log_now) return;
#ifdef MPICF
   if(vmpi::my_rank != 0) return;
#endif
   static std::ofstream out("sc1d-heff-tesla.log");
   static bool header = false;
   if(!header){
      header = true;
      out << "# H in Tesla. formula_* is Eq. (H-sc1d). llg_* is the vector stored for the next LLG step.\n";
      out << "# stored=0 during spin-current relax: the formula is evaluated and not written, so llg_* stays 0.\n";
      out << "# t_drive_ps is the drive time at this solver call (0 during relax). The laser envelope peaks at pump_time, not at half width.\n";
      out << "# S_* is the window-averaged spin accumulation in C/m^3. s = (mu_B/e) S is what enters H.\n";
      out << "phase\tstep\tsim_time\tt_drive_ps\tstored\tatom\tmaterial\tz_A"
          << "\tmx\tmy\tmz"
          << "\tllg_Hx_T\tllg_Hy_T\tllg_Hz_T"
          << "\tformula_Hx_T\tformula_Hy_T\tformula_Hz_T"
          << "\tformula_Hpar_T\tformula_Hperp_T"
          << "\tS_x\tS_y\tS_z\n";
      out.flush();
   }
   const bool relax = !store_field;
   const double t_drive_ps = relax ? 0.0
      : (mp::dt_SI * static_cast<double>(sim::time - sim::equilibration_time) * 1.0e12);
   const char* phase = relax ? "relax" : "drive";
   const double e_over_muB = kE / kMuB;
   const std::ios::fmtflags cout_flags = std::cout.flags();
   const std::streamsize cout_prec = std::cout.precision();
   auto emit = [&](const char* tag, const int atom, const Held& held){
      const int mat = (atom >= 0) ? atoms::type_array[atom] : -1;
      const std::string name = (mat >= 0 && mat < static_cast<int>(mp::material.size()))
         ? mp::material[static_cast<std::size_t>(mat)].name : tag;
      const double z = (atom >= 0) ? atoms::z_coord_array[atom] : 0.0;
      const double lx = (atom >= 0) ? Hx[atom] : 0.0;
      const double ly = (atom >= 0) ? Hy[atom] : 0.0;
      const double lz = (atom >= 0) ? Hz[atom] : 0.0;
      const double mdotr = held.px * held.px + held.py * held.py + held.pz * held.pz;
      const double hdotm = held.fx * held.px + held.fy * held.py + held.fz * held.pz;
      const double hpar = (mdotr > 0.0) ? (hdotm / std::sqrt(mdotr)) : 0.0;
      const double pxn = (mdotr > 0.0) ? held.px / std::sqrt(mdotr) : 0.0;
      const double pyn = (mdotr > 0.0) ? held.py / std::sqrt(mdotr) : 0.0;
      const double pzn = (mdotr > 0.0) ? held.pz / std::sqrt(mdotr) : 0.0;
      const double tx = held.fx - hpar * pxn;
      const double ty = held.fy - hpar * pyn;
      const double tz = held.fz - hpar * pzn;
      const double hperp = std::sqrt(tx * tx + ty * ty + tz * tz);
      std::ostringstream row;
      row << std::scientific << std::setprecision(8);
      row << phase << '\t' << st::internal::sc1d_step_counter << '\t' << sim::time
          << '\t' << t_drive_ps
          << '\t' << (store_field ? 1 : 0)
          << '\t' << atom << '\t' << name << '\t' << z
          << '\t' << held.px << '\t' << held.py << '\t' << held.pz
          << '\t' << lx << '\t' << ly << '\t' << lz
          << '\t' << held.fx << '\t' << held.fy << '\t' << held.fz
          << '\t' << hpar << '\t' << hperp
          << '\t' << (held.sx * e_over_muB) << '\t' << (held.sy * e_over_muB) << '\t' << (held.sz * e_over_muB)
          << '\n';
      out << row.str();
      std::ostringstream msg;
      msg << std::scientific << std::setprecision(6);
      msg << "sc1d Heff Tesla " << phase
          << " step " << st::internal::sc1d_step_counter
          << " t_drive_ps " << t_drive_ps
          << " " << name
          << " stored " << (store_field ? 1 : 0)
          << " llg " << lx << ' ' << ly << ' ' << lz
          << " formula " << held.fx << ' ' << held.fy << ' ' << held.fz
          << " Hperp " << hperp;
      std::cout << msg.str() << std::endl;
   };
   emit("Co-5", sample_co, held_co);
   emit("CoFeB", sample_cfb, held_cfb);
   out.flush();
   std::cout.flags(cout_flags);
   std::cout.precision(cout_prec);
}

}

namespace st{
namespace internal{

void cache_sc1d_torque_atoms(){
   torque_terms.clear();
   const int natoms = num_local_atoms;
   const int nf = sc1d_nf;
   const int nmat = static_cast<int>(mp.size());
   const double V_atom = eq6_atom_volume();
   if(natoms <= 0 || nf <= 0 || !(V_atom > 0.0) || nmat <= 0) return;

   const double muB_over_e = kMuB / kE;
   const bool use_fine_rates = !sc1d_Jsd_fine.empty()
      && !sc1d_D_fine.empty()
      && !sc1d_lphi_fine.empty();
   const int* stack_of = sc1d_atom_ls.data();
   const int* lo_of = sc1d_atom_fine_lo.data();
   const int* hi_of = sc1d_atom_fine_hi.data();
   const int* mat_of = atoms::type_array.data();
   const double* Jf_all = sc1d_Jsd_fine.data();
   const double* Df_all = sc1d_D_fine.data();
   const double* Lf_all = sc1d_lphi_fine.data();
   const std::size_t nf_sz = static_cast<std::size_t>(nf);

   torque_terms.reserve(static_cast<std::size_t>(natoms));
   for(int atom = 0; atom < natoms; ++atom){
      const int ls = stack_of[atom];
      const int i_lo = lo_of[atom];
      const int i_hi = hi_of[atom];
      if(ls < 0 || i_lo < 0 || i_hi < i_lo || i_hi >= nf) continue;

      const int mat = mat_of[atom];
      if(mat < 0 || mat >= nmat) continue;
      const double Jsd = mp[static_cast<std::size_t>(mat)].sd_exchange;
      // Zero marks no s-d coupling. A negative material exchange is kept.
      if(Jsd == 0.0) continue;

      const double mu = mp::material[static_cast<std::size_t>(mat)].mu_s_SI;
      const double gamma = mp::gamma_SI * mp::material[static_cast<std::size_t>(mat)].gamma_rel;
      if(!(mu > 0.0) || !(gamma > 0.0)) continue;

      const int ncell = i_hi - i_lo + 1;
      double c_field = Jsd / (2.0 * kHbar);
      double c_damp = 0.0;
      const double De = mp[static_cast<std::size_t>(mat)].diffusion;
      const double lphi = mp[static_cast<std::size_t>(mat)].lambda_phi;
      c_damp = (lphi > 0.0 && De > 0.0) ? (De / (lphi * lphi)) : 0.0;
      // The fine grid already holds the prolonged rates, including an analytic
      // sheet or a diffuse in-cell mix.
      if(use_fine_rates){
         const std::size_t base = static_cast<std::size_t>(ls) * nf_sz;
         const double* Jf = Jf_all + base;
         const double* Df = Df_all + base;
         const double* Lf = Lf_all + base;
         double sum_j = 0.0;
         double sum_d = 0.0;
         for(int i = i_lo; i <= i_hi; ++i){
            sum_j += Jf[i];
            if(Lf[i] > 0.0 && Df[i] > 0.0) sum_d += Df[i] / (Lf[i] * Lf[i]);
         }
         const double invn = 1.0 / static_cast<double>(ncell);
         c_field = (sum_j * invn) / (2.0 * kHbar);
         c_damp = sum_d * invn;
      }

      const double pref = V_atom / (gamma * mu);
      const double alpha = mp::material[static_cast<std::size_t>(mat)].alpha;
      TorqueTerm term;
      term.atom = atom;
      term.ncell = ncell;
      term.s_off = (static_cast<std::size_t>(ls) * nf_sz + static_cast<std::size_t>(i_lo)) * 3u;
      term.s_scale = muB_over_e / static_cast<double>(ncell);
      term.along = pref * (c_field + alpha * c_damp);
      term.cross = pref * (c_damp - alpha * c_field);
      torque_terms.push_back(term);
   }
}

} // namespace internal

   //-----------------------------------------------------------------------------
   // Function for updating spin torque fields
   //-----------------------------------------------------------------------------
   void update_spin_torque_fields(const std::vector<double>& x_spin_array,
                                  const std::vector<double>& y_spin_array,
                                  const std::vector<double>& z_spin_array,
                                  const std::vector<int>& atom_type_array,
                                  const std::vector<double>& mu_s_array){
       
      if(st::internal::enabled==false) return;

      if(st::internal::sc1d_enable){
         st::internal::calculate_spin_accumulation_1d();
         // Calibrate keeps the magnetisation source and every output, and
         // leaves the atom field at the zero it was given at initialisation.
         if(!st::internal::sc1d_calibrate){
            apply_sc1d_atom_spin_torque_fields(x_spin_array, y_spin_array, z_spin_array,
                                               atom_type_array, mu_s_array);
         }
         return;
      }

      // update magnetisations
      st::internal::update_cell_magnetisation(x_spin_array, y_spin_array, z_spin_array, atom_type_array, mu_s_array);

      // calculate spin_accumulation
      if(st::internal::sot_sa) st::internal::calculate_sot_accumulation();
      else st::internal::calculate_spin_accumulation();
      for(int atom=0; atom<st::internal::num_local_atoms; ++atom) {
         const int cell3 = 3*st::internal::atom_st_index[atom];
         const double i_mu_s = 1.0/(mu_s_array[atom_type_array[atom]]);
         st::internal::x_field_array[atom] = st::internal::spin_torque[cell3+0]*i_mu_s;
         st::internal::y_field_array[atom] = st::internal::spin_torque[cell3+1]*i_mu_s;
         st::internal::z_field_array[atom] = st::internal::spin_torque[cell3+2]*i_mu_s;
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
      // Calibrate: the solver still ran, but the atom does not feel the field.
      if(st::internal::sc1d_enable && st::internal::sc1d_calibrate) return;


      // Add spin torque fields
      for(int i=start_index; i<end_index; ++i) x_total_external_field_array[i] += st::internal::x_field_array[i];
      for(int i=start_index; i<end_index; ++i) y_total_external_field_array[i] += st::internal::y_field_array[i];
      for(int i=start_index; i<end_index; ++i) z_total_external_field_array[i] += st::internal::z_field_array[i];

      return;
   }

   //-----------------------------------------------------------------------------
   // Spin-current relax, after Monte-Carlo preconditioning and before ASD
   // equilibration starts the LLG. Order: precondition → relax → eq → drive.
   // The first transport step reseeds m-hat and S from the preconditioned
   // atoms (S = sa_inf * unit m-hat). Magnetisation is then held fixed for
   // spin-currents-1d-relax-steps. Torque stays zero and the applied bias
   // stays off until the drive phase.
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
