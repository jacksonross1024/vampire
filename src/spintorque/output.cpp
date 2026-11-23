//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>

// Vampire headers
#include "spintorque.hpp"
#include "vmpi.hpp"
#include "vio.hpp"
#include "sim.hpp"
#include "material.hpp"
#include "ltmp.hpp"
#include "create.hpp"

// Spin Torque headers
#include "internal.hpp"

namespace st{
   namespace internal{

      static void ensure_sc1d_output_directory(){
         if(vmpi::my_rank == 0){
            const int rc = mkdir("spin-acc", 0775);
            (void)rc;
         }
#ifdef MPICF
         MPI_Barrier(MPI_COMM_WORLD);
#endif
      }

      static void stack_xy(const int stack, int& ix, int& iy){
         const int ny = std::max(1, num_y_stacks);
         ix = stack / ny;
         iy = stack - ix * ny;
      }

      static bool sc1d_rank_owns_column(const int stack){
#ifdef MPICF
         const int nrank = std::max(1, vmpi::num_processors);
         return (stack % nrank) == vmpi::my_rank;
#else
         (void)stack;
         return true;
#endif
      }

      // Owners hold their columns. Sum onto rank 0 without leaving another
      // rank's values in the root buffer for the next step.
      static void reduce_owned_sum(std::vector<double>& field, const int stride){
#ifdef MPICF
         if(vmpi::num_processors <= 1 || field.empty() || stride <= 0) return;
         std::vector<double> send = field;
         for(int s = 0; s < num_stacks_y; ++s){
            if(sc1d_rank_owns_column(s)) continue;
            const std::size_t off = (std::size_t)s * (std::size_t)stride;
            const std::size_t n = std::min(send.size(), off + (std::size_t)stride);
            for(std::size_t i = off; i < n; ++i) send[i] = 0.0;
         }
         std::vector<double> recv(vmpi::my_rank == 0 ? field.size() : 0);
         MPI_Reduce(send.data(), vmpi::my_rank == 0 ? recv.data() : nullptr,
                    static_cast<int>(send.size()), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
         if(vmpi::my_rank == 0) field.swap(recv);
#else
         (void)field;
         (void)stride;
#endif
      }

      void output_sc1d_fine_constants(){
         if(!sc1d_enable) return;
         ensure_sc1d_output_directory();
         if(vmpi::my_rank != 0) return;
         const int nf = sc1d_nf;
         const int nsub = std::max(1, sc1d_nsub);
         if(nf <= 0 || num_stacks_y <= 0) return;
         const double dzA = micro_cell_thickness / (double)nsub;
         const double cell_xy = micro_cell_size[0];

         zlog << zTs() << "Fine-grid constants go to spin-acc/fine-constants." << std::endl;
         std::ofstream ofile("spin-acc/fine-constants");
         if(!ofile){
            zlog << zTs() << "Could not open spin-acc/fine-constants. Create the spin-acc directory before the run." << std::endl;
            return;
         }
         ofile << "# Fine-grid constants after initialisation, before the laser changes conductivity.\n";
         ofile << "# z is Angstroms from the substrate. D and De include the mixing sheet through sigma_up+sigma_down.\n";
         ofile << "# beta_c is the up/down imbalance. sigma_bulk is the conductivity before the sheet.\n";
         ofile << "# The hot-electron speed keeps the pre-sheet diffusivity and is not a column.\n";
         ofile << "# ix, iy are the column. iz is the fine cell along that column. Each column is written.\n";
         ofile << "# ix\tiy\tiz\tx(A)\ty(A)\tz(A)\tbeta_c\tbeta_d\tD(m^2/s)\tDe(m^2/s)\tlambda_sdl(m)\tlambda_phi(m)\tJsd(J)\tchi\tsa_inf(C/m^3)\tP_N\talpha\tsigma(S/m)\tsigma_bulk(S/m)\tS(V/K)\tbeta_seebeck\tP_hot\toccupied\tconnected\tspacing\tnatom\n";
         ofile << std::scientific << std::setprecision(6);
         for(int stack = 0; stack < num_stacks_y; ++stack){
            const int start = stack_index_y[stack];
            int ix = 0;
            int iy = 0;
            stack_xy(stack, ix, iy);
            const double xA = (pos[3*start+0] + 0.5) * cell_xy;
            const double yA = (pos[3*start+1] + 0.5) * cell_xy;
            const std::size_t base = (std::size_t)stack * (std::size_t)nf;
            for(int i = 0; i < nf; ++i){
               const std::size_t c = base + (std::size_t)i;
               const double zA = (i + 0.5) * dzA;
               ofile << ix << "\t" << iy << "\t" << i << "\t";
               ofile << std::fixed << std::setprecision(2);
               ofile << xA << "\t" << yA << "\t" << zA << "\t";
               ofile << std::scientific << std::setprecision(6);
               const unsigned char mask = (c < sc1d_cell_mask.size()) ? sc1d_cell_mask[c] : 0;
               const double natom = (c < sc1d_natom_fine.size()) ? sc1d_natom_fine[c] : 0.0;
               const double sigma_bulk = (c < sc1d_sigma_bulk_fine.size()) ? sc1d_sigma_bulk_fine[c] : 0.0;
               ofile << sc1d_Bc_fine[c] << "\t" << sc1d_Bd_fine[c] << "\t"
                     << sc1d_D_fine[c] << "\t" << sc1d_D_fine[c] << "\t"
                     << sc1d_lsf_fine[c] << "\t" << sc1d_lphi_fine[c] << "\t"
                     << sc1d_Jsd_fine[c] << "\t" << sc1d_chi_fine[c] << "\t"
                     << sc1d_sa_inf_fine[c] << "\t"
                     << sc1d_PN_fine[c] << "\t" << sc1d_alpha_fine[c] << "\t"
                     << sc1d_sigma_fine[c] << "\t"
                     << sigma_bulk << "\t"
                     << sc1d_seebeck_fine[c] << "\t" << sc1d_beta_seebeck_fine[c] << "\t"
                     << sc1d_phot_fine[c] << "\t"
                     << ((mask & sc1d_mask_occupied) ? 1 : 0) << "\t"
                     << ((mask & sc1d_mask_connected) ? 1 : 0) << "\t"
                     << ((mask & sc1d_mask_spacing) ? 1 : 0) << "\t"
                     << natom << "\n";
            }
         }
      }

      static void output_sc1d_verbose(const std::string& filename,
                                      const std::string& stage,
                                      const unsigned long index){
         const int nf = sc1d_nf;
         const int nsub = std::max(1, sc1d_nsub);
         const int nface = nf + 1;
         if(nf <= 0) return;
         reduce_owned_sum(sc1d_ns_fine, nf);
         reduce_owned_sum(sc1d_ne_fine, nf);
         reduce_owned_sum(sc1d_V_fine, nf);
         reduce_owned_sum(sc1d_Jc_edge_fine, nface);
         reduce_owned_sum(sc1d_J_restore_edge, nface);
         reduce_owned_sum(sc1d_J_cd_edge, nface);
         reduce_owned_sum(sc1d_J_sd_edge, nface);
         reduce_owned_sum(sc1d_J_seebeck_edge, nface);
         reduce_owned_sum(sc1d_J_hot_edge, nface);
         reduce_owned_sum(sc1d_Js_fine, 3 * nf);
         reduce_owned_sum(sc1d_mhat_fine, 3 * nf);
         if(vmpi::my_rank != 0) return;

         const double dzA = micro_cell_thickness / (double)nsub;
         const double cell_xy = micro_cell_size[0];
         std::ofstream ofile(filename.c_str());
         if(!ofile) return;
         ofile << "# stage " << stage << " index " << index
               << " step " << sc1d_step_counter << " sim_time " << sim::time << "\n";
         ofile << "# One row per fine cell of each column. z is the cell centre.\n";
         ofile << "# ix, iy are the column. iz is the fine cell. Verbose keeps the columns separate.\n";
         ofile << "# Charge components, Jc and Js are the average of the two faces of that cell, in A/m^2.\n";
         ofile << "# Jc is the stored series current, including the tunnel current on an insulating face.\n";
         ofile << "# Js is the face current the step carries. Diffusion is dz*A*(S_left-S_right),\n";
         ofile << "# which is zero on a closed face (either neighbour has D = 0). Drift is the edge\n";
         ofile << "# used by the divergence, including the tunnel value. Walls add beta*Jc*m and no diffusion.\n";
         ofile << "# V is the fine-cell potential in volts.\n";
         ofile << "# beta_c through sigma_bulk are the live fine-grid constants used by this step.\n";
         ofile << "# ix\tiy\tiz\tx(A)\ty(A)\tz(A)\tns(m^-3)\tne(m^-3)\tV(V)\tJ_field\tJ_cd\tJ_sd\tJ_seebeck\tJ_hot\tJc\tJs_x(A/m^2)\tJs_y(A/m^2)\tJs_z(A/m^2)\tSa_x(C/m^3)\tSa_y(C/m^3)\tSa_z(C/m^3)\tm_x\tm_y\tm_z\tTe(K)\tTp(K)\tbeta_c\tbeta_d\tbeta_seebeck\tP_hot\tD(m^2/s)\tDe(m^2/s)\tlambda_sdl(m)\tlambda_phi(m)\tJsd(J)\tsa_inf(C/m^3)\tP_N\talpha\tsigma(S/m)\tsigma_bulk(S/m)\n";
         for(int stack = 0; stack < num_stacks_y; ++stack){
            const int start = stack_index_y[stack];
            int ix = 0;
            int iy = 0;
            stack_xy(stack, ix, iy);
            const double xA = (pos[3*start+0] + 0.5) * cell_xy;
            const double yA = (pos[3*start+1] + 0.5) * cell_xy;
            const std::size_t base = (std::size_t)stack * (std::size_t)nf;
            const std::size_t ebase = (std::size_t)stack * (std::size_t)nface;
            const double* ns = sc1d_ns_fine.data() + base;
            const double* ne = sc1d_ne_fine.data() + base;
            const double* Vp = sc1d_V_fine.data() + base;
            const double* Jc = sc1d_Jc_edge_fine.data() + ebase;
            const double* S = sc1d_Sfine.data() + base * 3u;
            const double* Js = sc1d_Js_fine.data() + base * 3u;
            const double* mh = sc1d_mhat_fine.data() + base * 3u;
            const double* mlv = sc1d_ml_fine.data() + base;
            const double* Jf = sc1d_J_restore_edge.data() + ebase;
            const double* Jcd = sc1d_J_cd_edge.data() + ebase;
            const double* Jsd = sc1d_J_sd_edge.data() + ebase;
            const double* Jsee = sc1d_J_seebeck_edge.data() + ebase;
            const double* Jhot = sc1d_J_hot_edge.data() + ebase;
            const double* Bc = sc1d_Bc_fine.data() + base;
            const double* Bd = sc1d_Bd_fine.data() + base;
            const double* Bs = sc1d_beta_seebeck_fine.data() + base;
            const double* Ph = sc1d_phot_fine.data() + base;
            const double* D = sc1d_D_fine.data() + base;
            const double* De = sc1d_Dmat_fine.data() + base;
            const double* lsf = sc1d_lsf_fine.data() + base;
            const double* lphi = sc1d_lphi_fine.data() + base;
            const double* Jsd_fine = sc1d_Jsd_fine.data() + base;
            const double* sa_inf = sc1d_sa_inf_fine.data() + base;
            const double* PN = sc1d_PN_fine.data() + base;
            const double* alpha = sc1d_alpha_fine.data() + base;
            const double* sigma = sc1d_sigma_fine.data() + base;
            const double* sigma_bulk = sc1d_sigma_bulk_fine.data() + base;
            for(int i = 0; i < nf; ++i){
               const double zA = (static_cast<double>(i) + 0.5) * dzA;
               ofile << ix << "\t" << iy << "\t" << i << "\t";
               const double jf = 0.5 * (Jf[i] + Jf[i + 1]);
               const double jcd = 0.5 * (Jcd[i] + Jcd[i + 1]);
               const double jsd = 0.5 * (Jsd[i] + Jsd[i + 1]);
               const double jsee = 0.5 * (Jsee[i] + Jsee[i + 1]);
               const double jhot = 0.5 * (Jhot[i] + Jhot[i + 1]);
               const double jc = 0.5 * (Jc[i] + Jc[i + 1]);
               ofile << std::fixed << std::setprecision(2);
               ofile << xA << "\t" << yA << "\t" << zA << "\t";
               ofile << std::scientific << std::setprecision(6);
               ofile << ns[i] << "\t" << ne[i] << "\t" << Vp[i] << "\t"
                     << jf << "\t" << jcd << "\t" << jsd << "\t"
                     << jsee << "\t" << jhot << "\t"
                     << jc << "\t"
                     << Js[3*i+0] << "\t" << Js[3*i+1] << "\t" << Js[3*i+2] << "\t"
                     << S[3*i+0] << "\t" << S[3*i+1] << "\t" << S[3*i+2] << "\t"
                     << mlv[i] * mh[3*i+0] << "\t" << mlv[i] * mh[3*i+1] << "\t" << mlv[i] * mh[3*i+2] << "\t"
                     << sc1d_Te_fine[static_cast<std::size_t>(i)] << "\t"
                     << sc1d_Tp_fine[static_cast<std::size_t>(i)] << "\t"
                     << Bc[i] << "\t" << Bd[i] << "\t" << Bs[i] << "\t" << Ph[i] << "\t"
                     << D[i] << "\t" << De[i] << "\t"
                     << lsf[i] << "\t" << lphi[i] << "\t" << Jsd_fine[i] << "\t"
                     << sa_inf[i] << "\t" << PN[i] << "\t"
                     << alpha[i] << "\t" << sigma[i] << "\t" << sigma_bulk[i] << "\n";
            }
         }
      }

      static void reduce_sc1d_dynamics(){
         const int nf = sc1d_nf;
         const int nface = nf + 1;
         if(nf <= 0) return;
         reduce_owned_sum(sc1d_ns_fine, nf);
         reduce_owned_sum(sc1d_ne_fine, nf);
         reduce_owned_sum(sc1d_J_restore_edge, nface);
         reduce_owned_sum(sc1d_J_cd_edge, nface);
         reduce_owned_sum(sc1d_J_sd_edge, nface);
         reduce_owned_sum(sc1d_J_seebeck_edge, nface);
         reduce_owned_sum(sc1d_J_hot_edge, nface);
         reduce_owned_sum(sc1d_Js_fine, 3 * nf);
         reduce_owned_sum(sc1d_mhat_fine, 3 * nf);
      }

      struct sc1d_sample_t {
         double ns, ne, jf, jcd, jsd, jsee, jhot;
         double jsx, jsy, jsz, sax, say, saz, mx, my, mz;
         int n;
      };

      static void add_fine_cell(const int stack, const int i, sc1d_sample_t& s){
         const int nf = sc1d_nf;
         if(i < 0 || i >= nf) return;
         const std::size_t base = (std::size_t)stack * (std::size_t)nf + (std::size_t)i;
         const std::size_t ebase = (std::size_t)stack * (std::size_t)(nf + 1);
         const double* S = sc1d_Sfine.data() + (std::size_t)stack * (std::size_t)nf * 3u;
         const double* Js = sc1d_Js_fine.data() + (std::size_t)stack * (std::size_t)nf * 3u;
         const double* mh = sc1d_mhat_fine.data() + (std::size_t)stack * (std::size_t)nf * 3u;
         s.ns += sc1d_ns_fine[base];
         s.ne += sc1d_ne_fine[base];
         s.jf += 0.5 * (sc1d_J_restore_edge[ebase + (std::size_t)i] + sc1d_J_restore_edge[ebase + (std::size_t)i + 1u]);
         s.jcd += 0.5 * (sc1d_J_cd_edge[ebase + (std::size_t)i] + sc1d_J_cd_edge[ebase + (std::size_t)i + 1u]);
         s.jsd += 0.5 * (sc1d_J_sd_edge[ebase + (std::size_t)i] + sc1d_J_sd_edge[ebase + (std::size_t)i + 1u]);
         s.jsee += 0.5 * (sc1d_J_seebeck_edge[ebase + (std::size_t)i] + sc1d_J_seebeck_edge[ebase + (std::size_t)i + 1u]);
         s.jhot += 0.5 * (sc1d_J_hot_edge[ebase + (std::size_t)i] + sc1d_J_hot_edge[ebase + (std::size_t)i + 1u]);
         s.jsx += Js[3 * i + 0];
         s.jsy += Js[3 * i + 1];
         s.jsz += Js[3 * i + 2];
         s.sax += S[3 * i + 0];
         s.say += S[3 * i + 1];
         s.saz += S[3 * i + 2];
         s.mx += mh[3 * i + 0];
         s.my += mh[3 * i + 1];
         s.mz += mh[3 * i + 2];
         s.n += 1;
      }

      static void add_sample(sc1d_sample_t& dst, const sc1d_sample_t& src){
         dst.ns += src.ns;
         dst.ne += src.ne;
         dst.jf += src.jf;
         dst.jcd += src.jcd;
         dst.jsd += src.jsd;
         dst.jsee += src.jsee;
         dst.jhot += src.jhot;
         dst.jsx += src.jsx;
         dst.jsy += src.jsy;
         dst.jsz += src.jsz;
         dst.sax += src.sax;
         dst.say += src.say;
         dst.saz += src.saz;
         dst.mx += src.mx;
         dst.my += src.my;
         dst.mz += src.mz;
         dst.n += src.n;
      }

      static void write_sc1d_sample(std::ofstream& ofile, const int ix, const int iy, const int iz,
                                    const double xA, const double yA, const double zA, const sc1d_sample_t& s){
         if(s.n <= 0) return;
         const double inv = 1.0 / static_cast<double>(s.n);
         const double sax = s.sax * inv;
         const double say = s.say * inv;
         const double saz = s.saz * inv;
         double mx = s.mx * inv;
         double my = s.my * inv;
         double mz = s.mz * inv;
         const double mm = std::sqrt(mx * mx + my * my + mz * mz);
         if(mm > 1.0e-12){
            mx /= mm;
            my /= mm;
            mz /= mm;
         }
         const double sa_mag = std::sqrt(sax * sax + say * say + saz * saz);
         ofile << ix << "\t" << iy << "\t" << iz << "\t";
         ofile << std::fixed << std::setprecision(2);
         ofile << xA << "\t" << yA << "\t" << zA << "\t";
         ofile << std::scientific << std::setprecision(6);
         ofile << (s.ns * inv) << "\t" << (s.ne * inv) << "\t"
               << (s.jf * inv) << "\t" << (s.jcd * inv) << "\t" << (s.jsd * inv) << "\t"
               << (s.jsee * inv) << "\t" << (s.jhot * inv) << "\t"
               << ((s.jf + s.jcd + s.jsd + s.jsee + s.jhot) * inv) << "\t"
               << (s.jsx * inv) << "\t" << (s.jsy * inv) << "\t" << (s.jsz * inv) << "\t"
               << sax << "\t" << say << "\t" << saz << "\t" << sa_mag << "\t"
               << mx << "\t" << my << "\t" << mz << "\n";
      }

      static void output_sc1d_profile(const std::string& filename, const std::string& stage, const unsigned long index){
         const int nf = sc1d_nf;
         const int nsub = std::max(1, sc1d_nsub);
         if(nf <= 0) return;
         reduce_sc1d_dynamics();
         if(vmpi::my_rank != 0) return;
         const double dzA = micro_cell_thickness / static_cast<double>(nsub);
         std::ofstream ofile(filename.c_str());
         if(!ofile) return;
         ofile << "# stage " << stage << " index " << index
               << " step " << sc1d_step_counter << " sim_time " << sim::time << "\n";
         ofile << "# One row per output cell. The columns are summed. ix and iy are -1.\n";
         ofile << "# Coarse bins average a unit cell. Atomic bins average a field window. The sum is over columns.\n";
         ofile << "# Jc and Js are the face currents the step carries, averaged with the same cells.\n";
         ofile << "# Js diffusion is zero on a closed face. It is not a centred difference of S.\n";
         ofile << "# ix\tiy\tiz\tx(A)\ty(A)\tz(A)\tns(m^-3)\tne(m^-3)\tJ_field\tJ_cd\tJ_sd\tJ_seebeck\tJ_hot\tJc\tJs_x(A/m^2)\tJs_y(A/m^2)\tJs_z(A/m^2)\tSa_x(C/m^3)\tSa_y(C/m^3)\tSa_z(C/m^3)\tSa_mag(C/m^3)\tm_x\tm_y\tm_z\n";
         const double xA = 0.5 * cs::system_dimensions[0];
         const double yA = 0.5 * cs::system_dimensions[1];
         if(sc1d_output_mode == sc1d_out_atomic){
            for(std::size_t p = 0; p < sc1d_plane_z.size(); ++p){
               const int lo = sc1d_plane_lo[p];
               const int hi = sc1d_plane_hi[p];
               if(lo < 0 || hi < lo) continue;
               sc1d_sample_t sum = {};
               int per = 0;
               for(int stack = 0; stack < num_stacks_y; ++stack){
                  sc1d_sample_t s = {};
                  for(int i = lo; i <= hi; ++i) add_fine_cell(stack, i, s);
                  if(per == 0) per = s.n;
                  add_sample(sum, s);
               }
               sum.n = per;
               write_sc1d_sample(ofile, -1, -1, static_cast<int>(p), xA, yA, sc1d_plane_z[p], sum);
            }
         } else {
            double uc = cs::unit_cell.dimensions[2];
            if(!(uc > 0.0)) uc = micro_cell_thickness;
            const double z_max = static_cast<double>(nf) * dzA;
            const int nbin = std::max(1, static_cast<int>(std::ceil(z_max / uc)));
            for(int b = 0; b < nbin; ++b){
               const double z0 = static_cast<double>(b) * uc;
               const double z1 = z0 + uc;
               sc1d_sample_t sum = {};
               int per = 0;
               for(int stack = 0; stack < num_stacks_y; ++stack){
                  sc1d_sample_t s = {};
                  for(int i = 0; i < nf; ++i){
                     const double z = (static_cast<double>(i) + 0.5) * dzA;
                     if(z < z0 || z >= z1) continue;
                     add_fine_cell(stack, i, s);
                  }
                  if(per == 0) per = s.n;
                  add_sample(sum, s);
               }
               sum.n = per;
               write_sc1d_sample(ofile, -1, -1, b, xA, yA, z0 + 0.5 * uc, sum);
            }
         }
      }

      // Stage of a completed transport step.
      // relax: sc1d_step_counter is 1..relax-steps and sim::time is still 0.
      // eq: the update that runs as sim::time advances to equilibration-time-steps
      // is the last equilibration step, so the test is <=, not <.
      // drive: sim::time has moved past equilibration. The index is the count
      // inside that stage. sim_time in the header is always sim::time, which
      // stays 0 through relax.
      static void sc1d_completed_stage(std::string& stage, unsigned long& index){
         if(sc1d_step_counter <= sc1d_relax_steps){
            stage = "relax";
            index = sc1d_step_counter;
            return;
         }
         if(sim::time <= sim::equilibration_time){
            stage = "eq";
            index = static_cast<unsigned long>(sim::time);
            return;
         }
         stage = "drive";
         index = static_cast<unsigned long>(sim::time - sim::equilibration_time);
      }

      // Index 0 is the stage-entry frame, written once by the origin path.
      // A later step writes only when its stage index is a positive multiple
      // of ST-output-rate, so a zero-length stage cannot emit index 0 twice.
      static bool sc1d_periodic_due(const unsigned long index){
         const unsigned long rate = static_cast<unsigned long>(std::max(ST_output_rate, 1));
         return (index > 0UL) && (index % rate == 0UL);
      }

      static void write_sc1d_frame_body(const std::string& stage, const unsigned long index){
         std::stringstream filename;
         filename << "spin-acc/" << stage << "-"
                  << std::setw(7) << std::setfill('0') << index;
         if(sc1d_output_mode == sc1d_out_verbose) output_sc1d_verbose(filename.str(), stage, index);
         else output_sc1d_profile(filename.str(), stage, index);
         config_file_counter++;
      }

      static void write_sc1d_frame(const std::string& stage, const unsigned long index){
         ensure_sc1d_output_directory();
         #ifdef MPICF
            MPI_Barrier(MPI_COMM_WORLD);
         #endif
         write_sc1d_frame_body(stage, index);
      }

      //-----------------------------------------------------------------------------
      // Function to output 1D spin currents data
      //-----------------------------------------------------------------------------
      void output_sc1d_data(){
         if(!sc1d_enable) return;
         ensure_sc1d_output_directory();

         #ifdef MPICF
            MPI_Barrier(MPI_COMM_WORLD);
         #endif

         std::string stage;
         unsigned long index = 0UL;
         sc1d_completed_stage(stage, index);
         if(!sc1d_periodic_due(index)) return;
         write_sc1d_frame_body(stage, index);
      }

      // One index-0 file when an LLG stage is entered. A stage whose step
      // count is zero is not entered by the integrator, so it is not created.
      void output_sc1d_stage_origin(const std::string& stage){
         if(!sc1d_enable) return;
         static bool eq_done = false;
         static bool drive_done = false;
         if(stage == "eq"){
            if(eq_done) return;
            eq_done = true;
            if(sim::equilibration_time == 0) return;
            if(sim::time != 0) return;
            write_sc1d_frame("eq", 0UL);
            return;
         }
         if(stage == "drive"){
            if(drive_done) return;
            drive_done = true;
            if(sim::time != sim::equilibration_time) return;
            write_sc1d_frame("drive", 0UL);
         }
      }

      // Starting relax frame: S has been reseeded from the preconditioned
      // moments and the first diffusion step has not run. The transport
      // counter is already 1, because that seed lives in the first step.
      // Rank 0's gather is undone so the diffusion step still sees owner data.
      void output_sc1d_relax_origin(){
         if(!sc1d_enable || sc1d_relax_steps == 0UL) return;
         static bool wrote = false;
         if(wrote) return;
         wrote = true;

         struct kept_field {
            std::vector<double>* live;
            std::vector<double> saved;
         };
         std::vector<kept_field> kept;
         if(vmpi::my_rank == 0){
            std::vector<double>* fields[] = {
               &sc1d_Sfine, &sc1d_ns_fine, &sc1d_ne_fine, &sc1d_V_fine,
               &sc1d_Jc_edge_fine, &sc1d_J_restore_edge, &sc1d_J_cd_edge,
               &sc1d_J_sd_edge, &sc1d_J_seebeck_edge, &sc1d_J_hot_edge,
               &sc1d_Js_fine, &sc1d_mhat_fine
            };
            const int nfields = static_cast<int>(sizeof(fields) / sizeof(fields[0]));
            kept.reserve(static_cast<std::size_t>(nfields));
            for(int f = 0; f < nfields; ++f){
               kept_field item;
               item.live = fields[f];
               item.saved = *fields[f];
               kept.push_back(item);
            }
         }
         if(sc1d_nf > 0) reduce_owned_sum(sc1d_Sfine, 3 * sc1d_nf);
         write_sc1d_frame("relax", 0UL);
         for(std::size_t i = 0; i < kept.size(); ++i) kept[i].live->swap(kept[i].saved);
      }

      //-----------------------------------------------------------------------------
      // One-line cumulative timing for the 1D solver. Same output cadence as
      // output_sc1d_data. Rank 0 only; other ranks keep their own accumulators.
      //-----------------------------------------------------------------------------
#if ST_SC1D_TIMINGS
      static void append_sc1d_bucket(std::ostringstream& line,
                                     const char* name,
                                     const double seconds,
                                     const double inv_pct){
         std::ostringstream sec;
         sec << std::scientific << std::setprecision(3) << seconds;
         std::ostringstream pcts;
         pcts << std::fixed << std::setprecision(1) << (seconds * inv_pct);
         line << " " << name << "=" << sec.str() << " (" << pcts.str() << "%)";
      }
#endif

      void report_sc1d_timing(){
#if !ST_SC1D_TIMINGS
         return;
#else
         if(!sc1d_enable) return;

         std::string stage_name;
         unsigned long index = 0UL;
         sc1d_completed_stage(stage_name, index);
         (void)stage_name;
         if(!sc1d_periodic_due(index)) return;

         double total = sc1d_time_total;
         double io = sc1d_time_io;
         double charge = sc1d_time_charge;
         double spin_current = sc1d_time_spin_current;
         double spin_acc = sc1d_time_spin_acc;
         double interp = sc1d_time_interp;
         double thermal = sc1d_time_thermal;
         double bcast = sc1d_time_bcast;
         double other = sc1d_time_other;
#ifdef MPICF
         // Print the rank that did the most column work. Idle ranks still
         // sit in the broadcast, so a plain rank-0 sample can miss the solve.
         double local_times[9] = {total, io, charge, spin_current, spin_acc, interp, thermal, bcast, other};
         const int nrank = std::max(1, vmpi::num_processors);
         std::vector<double> all_times((std::size_t)nrank * 9u, 0.0);
         MPI_Allgather(local_times, 9, MPI_DOUBLE, all_times.data(), 9, MPI_DOUBLE, MPI_COMM_WORLD);
         if(vmpi::my_rank != 0) return;
         int busiest = 0;
         double busiest_work = -1.0;
         for(int r=0; r<nrank; ++r){
            const double* row = all_times.data() + (std::size_t)r * 9u;
            const double work = row[2] + row[3] + row[4] + row[5] + row[6];
            if(work > busiest_work){
               busiest_work = work;
               busiest = r;
            }
         }
         const double* best = all_times.data() + (std::size_t)busiest * 9u;
         total = best[0];
         io = best[1];
         charge = best[2];
         spin_current = best[3];
         spin_acc = best[4];
         interp = best[5];
         thermal = best[6];
         bcast = best[7];
         other = best[8];
#endif
         const double inv_pct = (total > 0.0) ? (100.0 / total) : 0.0;

         std::ostringstream line;
         std::ostringstream tot;
         tot << std::scientific << std::setprecision(3) << total;
         line << "sc1d timing [s] total=" << tot.str();
         append_sc1d_bucket(line, "io", io, inv_pct);
         append_sc1d_bucket(line, "charge", charge, inv_pct);
         append_sc1d_bucket(line, "spin_current", spin_current, inv_pct);
         append_sc1d_bucket(line, "spin_acc", spin_acc, inv_pct);
         append_sc1d_bucket(line, "interp_coarse", interp, inv_pct);
         append_sc1d_bucket(line, "thermal_refresh", thermal, inv_pct);
         append_sc1d_bucket(line, "bcast", bcast, inv_pct);
         append_sc1d_bucket(line, "other", other, inv_pct);

         std::cout << line.str() << std::endl;
         zlog << zTs() << line.str() << std::endl;
#endif
      }

      void report_sc1d_profile(const double run_s){
         if(!sc1d_enable) return;

         const char* scheme = (sc1d_diffusion == sc1d_diff_bdf2) ? "bdf2" : "backward-euler";
         double diffusion = sc1d_time_diffusion;
         double charge = sc1d_time_charge;
         double eom = sc1d_time_eom;
         double prolong = sc1d_time_interp;
         double thermal = sc1d_time_thermal;
         double io = sc1d_time_io;
         double drift = sc1d_time_spin_current;
         double bcast = sc1d_time_bcast;
         double coarsen = sc1d_time_coarsen;
         double other = sc1d_time_other;
         double total = sc1d_time_total;
#ifdef MPICF
         double local_times[11] = {
            diffusion, charge, eom, prolong, thermal, io, drift, bcast, coarsen, other, total
         };
         const int nrank = std::max(1, vmpi::num_processors);
         std::vector<double> all_times((std::size_t)nrank * 11u, 0.0);
         MPI_Allgather(local_times, 11, MPI_DOUBLE, all_times.data(), 11, MPI_DOUBLE, MPI_COMM_WORLD);
         if(vmpi::my_rank != 0) return;
         int busiest = 0;
         double busiest_work = -1.0;
         for(int r = 0; r < nrank; ++r){
            const double* row = all_times.data() + (std::size_t)r * 11u;
            const double work = row[0] + row[1] + row[2] + row[3] + row[4];
            if(work > busiest_work){
               busiest_work = work;
               busiest = r;
            }
         }
         const double* best = all_times.data() + (std::size_t)busiest * 11u;
         diffusion = best[0];
         charge = best[1];
         eom = best[2];
         prolong = best[3];
         thermal = best[4];
         io = best[5];
         drift = best[6];
         bcast = best[7];
         coarsen = best[8];
         other = best[9];
         total = best[10];
#endif
         const double inv = (total > 0.0) ? (100.0 / total) : 0.0;
         const double of_run = (run_s > 0.0) ? (100.0 * total / run_s) : 0.0;

         std::ostringstream out;
         out << std::fixed << std::setprecision(3);
         out << "Spin-current profile, initialisation excluded. Scheme " << scheme << "\n";
         out << "  diffusion_solve        " << diffusion << " s  " << (diffusion * inv) << "% of spin steps\n";
         out << "  charge_poisson         " << charge << " s  " << (charge * inv) << "% of spin steps\n";
         out << "  spin_eom               " << eom << " s  " << (eom * inv) << "% of spin steps\n";
         out << "  m_prolongation         " << prolong << " s  " << (prolong * inv) << "% of spin steps\n";
         out << "  thermal_refresh        " << thermal << " s  " << (thermal * inv) << "% of spin steps\n";
         out << "  io                     " << io << " s  " << (io * inv) << "% of spin steps\n";
         out << "  drift_assembly         " << drift << " s  " << (drift * inv) << "% of spin steps\n";
         out << "  coarsen_torque         " << coarsen << " s  " << (coarsen * inv) << "% of spin steps\n";
         out << "  broadcast              " << bcast << " s  " << (bcast * inv) << "% of spin steps\n";
         out << "  other                  " << other << " s  " << (other * inv) << "% of spin steps\n";
         out << "  spin_steps_total       " << total << " s  " << of_run << "% of simulation run time\n";

         std::cout << out.str() << std::flush;
         zlog << zTs() << out.str() << std::flush;
      }


      //-----------------------------------------------------------------------------
      // Function to output base microcell properties
      //-----------------------------------------------------------------------------
      void output_base_microcell_data(){
         ensure_sc1d_output_directory();

         using st::internal::beta_cond;
         using st::internal::beta_diff;
         using st::internal::sa_infinity;
         using st::internal::lambda_sdl;
         using st::internal::pos;

         const int num_cells = beta_cond.size();

         // only output on root process
         if(vmpi::my_rank==0){
            if(sim::time%(ST_output_rate) ==0){
               zlog << zTs() << "Outputting ST base microcell data" << std::endl;
               std::ofstream ofile;
               ofile.open("spin-acc/st-microcells-base.cfg");
               ofile << num_cells << std::endl;
               ofile << num_stacks_y << std::endl;
               for(int cell=0; cell < num_cells; ++cell){
                  // if(cell_natom[cell] == 0) continue;
                  if(sot_sa) {
                        ofile << cell_stack_index[cell] << "\t" << "\t" << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] \
                        << "\t" << beta_cond[cell] << "\t" << beta_diff[cell] << "\t" << sa_infinity[cell] << "\t" << lambda_sdl[cell] << "\t" << \
                        st::internal::a[cell] << "\t" << st::internal::b[cell] << "\t" \
                        << "\t" << st::internal::sot_beta_cond[cell] << "\t" << st::internal::sot_beta_diff[cell] << "\t" << st::internal::sot_sa_infinity[cell] << "\t" << st::internal::sot_lambda_sdl[cell] << "\t" \
                        << st::internal::sot_a[cell] << "\t" << st::internal::sot_b[cell] << "\t" << st::internal::spin_acc_sign[cell] << "\t" << st::internal::sot_sa_source[cell] << "\t" << st::internal::cell_natom[cell] << std::endl;
                  } else {
                        ofile << cell_stack_index[cell] << "\t" << "\t" << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] \
                        << "\t" << beta_cond[cell] << "\t" << beta_diff[cell] << "\t" << sa_infinity[cell] << "\t" << lambda_sdl[cell] << "\t" << 
                        st::internal::a[cell] << "\t" << st::internal::b[cell] << "\t" << st::internal::cell_natom[cell] << std::endl;
                     }
               }    
		         ofile.close();
            }
         }

      return;
      }

      //-----------------------------------------------------------------------------
      // Function to output base microcell properties
      //-----------------------------------------------------------------------------
      void output_microcell_sa_data(){
         
         // If 1D spin currents solver is enabled, skip non-1D output
         if(sc1d_enable) return;

         // only output on root process
         #ifdef MPICF
           MPI_Barrier(MPI_COMM_WORLD);

            if(sim::time%(ST_output_rate) ==0){

               const int size = sa_final.size();
               const int num_cells = size/3;
               // determine file name
               std::stringstream filename;
               filename << "spin-acc/" << config_file_counter;
          
               MPI_Reduce(&sa_final[0], &sa_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&j_final_up_x[0], &j_final_up_x_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&j_final_up_y[0], &j_final_up_y_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&j_final_down_y[0], &j_final_down_y_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&coeff_ast[0], &coeff_ast_sum[0], num_cells, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&coeff_nast[0], &coeff_nast_sum[0], num_cells, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&ast[0], &ast_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&nast[0], &nast_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&total_ST[0], &total_ST_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                          
            if(vmpi::my_rank == 0) {
               std::ofstream ofile;
               ofile.open(std::string(filename.str()).c_str());
	          
               ofile << "pos_x \t pos_y \t pos_z \t m_x \t m_y \t m_z \t spin_acc_x \t spin_acc_y \t spin_acc_z \t " << \
               "jup_x \t jup_y \t jup_z \t jdown_x \t jdown_y \t jdown_z \t ast_x \t ast_y \t ast_z \t nast_x \t nast_y \t nast_z \t torque_x \t torque_y \t torque_z \t num_atom" << std::endl;
               for(int cell=0; cell<num_cells; ++cell){
                  //   if( (st::internal::cell_stack_index[cell]-1)%3 == 0) continue;
                  if(cell_natom[cell] == 0) continue;
                  double mag = sqrt(m[3*cell+0]*m[3*cell+0] + m[3*cell+1]*m[3*cell+1] + m[3*cell+2]*m[3*cell+2]);
                  mag = (mag == 0.0) ? 0.0: 1/mag;
                  ofile << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] << "\t";
                  ofile << m[3*cell+0]*mag << "\t" << m[3*cell+1]*mag << "\t" << m[3*cell+2]*mag << "\t";
                 // if(st::internal::sot_check) ofile << (sa_sum[3*cell+0]-sa_infinity[cell]*m[3*cell]*mag)/sa_infinity[cell] << "\t" << (sa_sum[3*cell+1]-sa_infinity[cell]*m[3*cell+1]*mag)/sa_infinity[cell] << "\t" << (sa_sum[3*cell+2]-m[3*cell+2]*mag)/sa_infinity[cell] << "\t";
                  //else 
                  ofile << std::fixed << std::setprecision(10);
                  ofile << sa_sum[3*cell+0] << "\t" << sa_sum[3*cell+1] << "\t" << sa_sum[3*cell+2] << "\t";
                  ofile << std::setprecision(6);
                  ofile << j_final_up_x_sum[3*cell+0] << "\t" << j_final_up_x_sum[3*cell+1] << "\t" << j_final_up_x_sum[3*cell+2] << "\t";
                  ofile << j_final_up_y_sum[3*cell+0] << "\t" << j_final_up_y_sum[3*cell+1] << "\t" << j_final_up_y_sum[3*cell+2] << "\t";
                  ofile << j_final_down_y_sum[3*cell+0] << "\t" << j_final_down_y_sum[3*cell+1] << "\t" << j_final_down_y_sum[3*cell+2] << "\t";
                  ofile << coeff_ast_sum[cell] << "\t";
                  ofile << coeff_nast_sum[cell] << "\t";
                  ofile << ast_sum[3*cell+0] << "\t" << ast_sum[3*cell+1] << "\t" << ast_sum[3*cell+2] << "\t";
                  ofile << nast_sum[3*cell+0] << "\t" << nast_sum[3*cell+1] << "\t" << nast_sum[3*cell+2] << "\t";
                  ofile << total_ST_sum[3*cell+0] << "\t" << total_ST_sum[3*cell+1] << "\t" << total_ST_sum[3*cell+2] << "\t" << sqrt(total_ST_sum[3*cell+0]*total_ST_sum[3*cell+0] + total_ST_sum[3*cell+1]*total_ST_sum[3*cell+1] + total_ST_sum[3*cell+2]*total_ST_sum[3*cell+2]);
                  ofile << "\t" << cell_natom[cell] << "\n";

                  sa_sum[3*cell+0] = 0.0;
                  sa_sum[3*cell+1] = 0.0;
                  sa_sum[3*cell+2] = 0.0;
                  j_final_up_x_sum[3*cell+0] = j_final_up_x_sum[3*cell+1] = j_final_up_x_sum[3*cell+2] = 0.0;
                  j_final_up_y_sum[3*cell+0] = j_final_up_y_sum[3*cell+1] = j_final_up_y_sum[3*cell+2] = 0.0;
                  j_final_down_y_sum[3*cell+0] = j_final_down_y_sum[3*cell+1] = j_final_down_y_sum[3*cell+2] = 0.0;
                  coeff_ast_sum[cell] = 0.0;
                  coeff_nast_sum[cell]  = 0.0;
                  ast_sum[3*cell+0] = ast_sum[3*cell+1] = ast_sum[3*cell+2] = 0.0;
                  nast_sum[3*cell+0] = nast_sum[3*cell+1] = nast_sum[3*cell+2] = 0.0;
                  total_ST_sum[3*cell+0] = total_ST_sum[3*cell+1] = total_ST_sum[3*cell+2] = 0.0;
           
               }

            ofile.close();
       
            // update config_file_counter
           }
         config_file_counter++;
         }
         #else
      if(sim::time%(ST_output_rate) ==0){ 

      // using namespace st::internal;
         const int size = m.size();
         const int num_cells = size/3;

            // determine file name
            std::stringstream filename;
            filename << "spin-acc/" << config_file_counter;
         
         
                  std::ofstream ofile;
            ofile.open(std::string(filename.str()).c_str());
            //  ofile<<"Time:"<< "\t" << sim::time*mp::dt_SI<< std::endl;
            ofile << "pos_x \t pos_y \t pos_z \t m_x \t m_y \t m_z \t spin_acc_x \t spin_acc_y \t spin_acc_z \t " << \
            "j_x \t j_y \t j_z \t ast_x \t ast_y \t ast_z \t nast_x \t nast_y \t nast_z \t torque_x \t torque_y \t torque_z \t num_atom" << std::endl;
            for(int cell=0; cell<num_cells; ++cell){
               //   if( (st::internal::cell_stack_index[cell]-1)%3 == 0) continue;
               if(cell_natom[cell] == 0) continue;
               double mag = sqrt(m[3*cell+0]*m[3*cell+0] + m[3*cell+1]*m[3*cell+1] + m[3*cell+2]*m[3*cell+2]);
               mag = (mag == 0.0) ? 0.0: 1/mag;
               ofile << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] << "\t";
               ofile << m[3*cell+0] << "\t" << m[3*cell+1] << "\t" << m[3*cell+2] << "\t";
               if(st::internal::sot_check) ofile << (sa_final[3*cell+0]-sa_infinity[cell]*m[3*cell]*mag)/sa_infinity[cell] << "\t" << (sa_final[3*cell+1]-sa_infinity[cell]*m[3*cell+1]*mag)/sa_infinity[cell] << "\t" << (sa_final[3*cell+2]-m[3*cell+2]*mag)/sa_infinity[cell] << "\t";
               else ofile << sa_final[3*cell+0] << "\t" << sa_final[3*cell+1] << "\t" << sa_final[3*cell+2] << "\t";
               ofile << j_final_up_x[3*cell+0] << "\t" << j_final_up_x[3*cell+1] << "\t" << j_final_up_x[3*cell+2] << "\t";
               ofile << j_final_up_y[3*cell+0] << "\t" << j_final_up_y[3*cell+1] << "\t" << j_final_up_y[3*cell+2] << "\t";
               ofile << j_final_down_y[3*cell+0] << "\t" << j_final_down_y[3*cell+1] << "\t" << j_final_down_y[3*cell+2] << "\t";
               ofile << coeff_ast[cell] << "\t";
               ofile << coeff_nast[cell] << "\t";
               ofile << ast[3*cell+0] << "\t" << ast[3*cell+1] << "\t" << ast[3*cell+2] << "\t";
               ofile << nast[3*cell+0] << "\t" << nast[3*cell+1] << "\t" << nast[3*cell+2] << "\t";
               ofile << total_ST[3*cell+0] << "\t" << total_ST[3*cell+1] << "\t" << total_ST[3*cell+2];
               ofile << "\t" << cell_natom[cell] << "\n";
               
            }

            ofile.close();
       
            // update config_file_counter
         config_file_counter++;
      }
           
   #endif 
   
   return;
   }

      void output_microcell_data(){

         // If 1D spin currents solver is enabled, only output 1D data and return
         if(sc1d_enable){
            output_sc1d_data();
            return;
         }

         // only output on root process
         #ifdef MPICF 
           MPI_Barrier(MPI_COMM_WORLD);

            if(sim::time%(ST_output_rate) ==0){

               const int size = sa_final.size();
               const int num_cells = size/3;
               // determine file name
               std::stringstream filename;
               filename << "spin-acc/" << config_file_counter;
          
               MPI_Reduce(&sa_final[0], &sa_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&j_final_up_x[0], &j_final_up_x_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               // MPI_Reduce(&j_final_up_y[0], &j_final_up_y_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               // MPI_Reduce(&j_final_down_y[0], &j_final_down_y_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&coeff_ast[0], &coeff_ast_sum[0], num_cells, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&coeff_nast[0], &coeff_nast_sum[0], num_cells, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&ast[0], &ast_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&nast[0], &nast_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&total_ST[0], &total_ST_sum[0], size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                          
            if(vmpi::my_rank == 0) {
               std::ofstream ofile;
               ofile.open(std::string(filename.str()).c_str());
	          
               ofile << "pos_x \t pos_y \t pos_z \t m_x \t m_y \t m_z \t spin_acc_x \t spin_acc_y \t spin_acc_z \t " << \
               "jup_x \t jup_y \t jup_z \t jdown_x \t jdown_y \t jdown_z \t ast_x \t ast_y \t ast_z \t nast_x \t nast_y \t nast_z \t torque_x \t torque_y \t torque_z \t num_atom" << std::endl;
               for(int cell=0; cell<num_cells; ++cell){
                  //   if( (st::internal::cell_stack_index[cell]-1)%3 == 0) continue;
                  if(cell_natom[cell] == 0) continue;
                  double mag = sqrt(m[3*cell+0]*m[3*cell+0] + m[3*cell+1]*m[3*cell+1] + m[3*cell+2]*m[3*cell+2]);
                  mag = 1;//(mag == 0.0) ? 0.0: 1/mag;
                  ofile << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] << "\t";
                  ofile << m[3*cell+0]*mag << "\t" << m[3*cell+1]*mag << "\t" << m[3*cell+2]*mag << "\t";
                 // if(st::internal::sot_check) ofile << (sa_sum[3*cell+0]-sa_infinity[cell]*m[3*cell]*mag)/sa_infinity[cell] << "\t" << (sa_sum[3*cell+1]-sa_infinity[cell]*m[3*cell+1]*mag)/sa_infinity[cell] << "\t" << (sa_sum[3*cell+2]-m[3*cell+2]*mag)/sa_infinity[cell] << "\t";
                  //else 
                  ofile << std::fixed << std::setprecision(10);
                  ofile << sa_sum[3*cell+0] << "\t" << sa_sum[3*cell+1] << "\t" << sa_sum[3*cell+2] << "\t";
                  ofile << std::setprecision(6);
                  ofile << j_final_up_x_sum[3*cell+0] << "\t" << j_final_up_x_sum[3*cell+1] << "\t" << j_final_up_x_sum[3*cell+2] << "\t";
                  // ofile << j_final_up_y_sum[3*cell+0] << "\t" << j_final_up_y_sum[3*cell+1] << "\t" << j_final_up_y_sum[3*cell+2] << "\t";
                  // ofile << j_final_down_y_sum[3*cell+0] << "\t" << j_final_down_y_sum[3*cell+1] << "\t" << j_final_down_y_sum[3*cell+2] << "\t";
                  ofile << coeff_ast_sum[cell] << "\t";
                  ofile << coeff_nast_sum[cell] << "\t";
                  ofile << ast_sum[3*cell+0] << "\t" << ast_sum[3*cell+1] << "\t" << ast_sum[3*cell+2] << "\t";
                  ofile << nast_sum[3*cell+0] << "\t" << nast_sum[3*cell+1] << "\t" << nast_sum[3*cell+2] << "\t";
                  ofile << total_ST_sum[3*cell+0] << "\t" << total_ST_sum[3*cell+1] << "\t" << total_ST_sum[3*cell+2] << "\t" << sqrt(total_ST_sum[3*cell+0]*total_ST_sum[3*cell+0] + total_ST_sum[3*cell+1]*total_ST_sum[3*cell+1] + total_ST_sum[3*cell+2]*total_ST_sum[3*cell+2]);
                  ofile << "\t" << cell_natom[cell] << "\n";

                  sa_sum[3*cell+0] = 0.0;
                  sa_sum[3*cell+1] = 0.0;
                  sa_sum[3*cell+2] = 0.0;
                  j_final_up_x_sum[3*cell+0] = j_final_up_x_sum[3*cell+1] = j_final_up_x_sum[3*cell+2] = 0.0;
                  // j_final_up_y_sum[3*cell+0] = j_final_up_y_sum[3*cell+1] = j_final_up_y_sum[3*cell+2] = 0.0;
                  // j_final_down_y_sum[3*cell+0] = j_final_down_y_sum[3*cell+1] = j_final_down_y_sum[3*cell+2] = 0.0;
                  coeff_ast_sum[cell] = 0.0;
                  coeff_nast_sum[cell]  = 0.0;
                  ast_sum[3*cell+0] = ast_sum[3*cell+1] = ast_sum[3*cell+2] = 0.0;
                  nast_sum[3*cell+0] = nast_sum[3*cell+1] = nast_sum[3*cell+2] = 0.0;
                  total_ST_sum[3*cell+0] = total_ST_sum[3*cell+1] = total_ST_sum[3*cell+2] = 0.0;
           
               }

            ofile.close();
       
            // update config_file_counter
           }
         config_file_counter++;
         }
         #else
      if(sim::time%(ST_output_rate) ==0){ 

      // using namespace st::internal;
         const int size = m.size();
         const int num_cells = size/3;

            // determine file name
            std::stringstream filename;
            filename << "spin-acc/" << config_file_counter;
         
         
                  std::ofstream ofile;
            ofile.open(std::string(filename.str()).c_str());
            //  ofile<<"Time:"<< "\t" << sim::time*mp::dt_SI<< std::endl;
            ofile << "pos_x \t pos_y \t pos_z \t m_x \t m_y \t m_z \t spin_acc_x \t spin_acc_y \t spin_acc_z \t " << \
            "j_x \t j_y \t j_z \t ast_x \t ast_y \t ast_z \t nast_x \t nast_y \t nast_z \t torque_x \t torque_y \t torque_z \t num_atom" << std::endl;
            for(int cell=0; cell<num_cells; ++cell){
               //   if( (st::internal::cell_stack_index[cell]-1)%3 == 0) continue;
               if(cell_natom[cell] == 0) continue;
               double mag = sqrt(m[3*cell+0]*m[3*cell+0] + m[3*cell+1]*m[3*cell+1] + m[3*cell+2]*m[3*cell+2]);
               mag = (mag == 0.0) ? 0.0: 1/mag;
               ofile << pos[3*cell+0] << "\t" << pos[3*cell+1] << "\t" << pos[3*cell+2] << "\t";
               ofile << m[3*cell+0] << "\t" << m[3*cell+1] << "\t" << m[3*cell+2] << "\t";
               // if(st::internal::sot_check) ofile << (sa_final[3*cell+0]-sa_infinity[cell]*m[3*cell]*mag)/sa_infinity[cell] << "\t" << (sa_final[3*cell+1]-sa_infinity[cell]*m[3*cell+1]*mag)/sa_infinity[cell] << "\t" << (sa_final[3*cell+2]-m[3*cell+2]*mag)/sa_infinity[cell] << "\t";
               // else 
               ofile << sa_final[3*cell+0] << "\t" << sa_final[3*cell+1] << "\t" << sa_final[3*cell+2] << "\t";
               ofile << j_final_up_x[3*cell+0] << "\t" << j_final_up_x[3*cell+1] << "\t" << j_final_up_x[3*cell+2] << "\t";
               // ofile << j_final_up_y[3*cell+0] << "\t" << j_final_up_y[3*cell+1] << "\t" << j_final_up_y[3*cell+2] << "\t";
               // ofile << j_final_down_y[3*cell+0] << "\t" << j_final_down_y[3*cell+1] << "\t" << j_final_down_y[3*cell+2] << "\t";
               ofile << coeff_ast[cell] << "\t";
               ofile << coeff_nast[cell] << "\t";
               ofile << ast[3*cell+0] << "\t" << ast[3*cell+1] << "\t" << ast[3*cell+2] << "\t";
               ofile << nast[3*cell+0] << "\t" << nast[3*cell+1] << "\t" << nast[3*cell+2] << "\t";
               ofile << total_ST[3*cell+0] << "\t" << total_ST[3*cell+1] << "\t" << total_ST[3*cell+2];
               ofile << "\t" << cell_natom[cell] << "\n";
               
            }

            ofile.close();
       
            // update config_file_counter
         config_file_counter++;
      }
           
   #endif 
   
   return;
   }

   //-----------------------------------------------------------------------------
   // Thermal gradient output functions
   //-----------------------------------------------------------------------------
   
   // Local variables for thermal gradient output
   std::ofstream thermal_microcell_file;
   std::ofstream thermal_temperature_file;
   int thermal_output_counter = 0;
   bool thermal_output_initialised = false;

   //-----------------------------------------------------------------------------
   // Function to output thermal gradient microcell configuration
   //-----------------------------------------------------------------------------
   void output_thermal_microcell_data(){
      if(!sc1d_thermal_gradients_enable || !ltmp::is_enabled()) return;
      
      // only output on root process
      #ifdef MPICF
      if(vmpi::my_rank==0){
      #else
      {
      #endif
         std::ofstream ofile;
         ofile.open("spin-acc/thermal_microcell_config.cfg");
         
         // Header: cell position (Angstroms), attenuation, material constants
         ofile << "# z(A)\tattenuation\tC_e(J/m^3/K)\tC_p(J/m^3/K)\tG_ep(J/s/m^3/K)\tk_e(J/s/m/K)\tk_p(J/s/m/K)\tT_Debye(K)\n";
         
         const int ncz = num_microcells_per_stack;
         const double dz = micro_cell_thickness; // Angstroms
         const double d = std::max(sc1d_optical_absorption_length * 1e10, 1e-8); // m to Angstroms
         const double z_max = (double)ncz * dz; // total stack height in Angstroms
         
         for(int cell = 0; cell < ncz; ++cell){
            const double z = (cell + 0.5) * dz; // cell center position
            const double distance_from_top = z_max - z; // distance from top surface
            const double attenuation = std::exp(-distance_from_top / d);
            
            ofile << std::fixed << std::setprecision(2);
            ofile << z << "\t";
            ofile << std::scientific << std::setprecision(6);
            ofile << attenuation << "\t";
            if(cell < (int)sc1d_Ce.size()){
               ofile << sc1d_Ce[cell] << "\t";
               ofile << sc1d_Cp[cell] << "\t";
               ofile << sc1d_G[cell] << "\t";
               ofile << sc1d_kappa_e[cell] << "\t";
               ofile << sc1d_kappa_p[cell] << "\t";
               ofile << sc1d_T_Debye[cell] << "\n";
            } else {
               ofile << "0\t0\t0\t0\t0\t0\n";
            }
         }
         
         ofile.close();
      }
      
      return;
   }

   //-----------------------------------------------------------------------------
   // Function to open temperature profile output file
   //-----------------------------------------------------------------------------
   void open_thermal_temperature_profile_file(){
      if(!sc1d_thermal_gradients_enable || !ltmp::is_enabled()) return;
      
      thermal_output_counter = 0;
      thermal_output_initialised = true;
      
      #ifdef MPICF
      if(vmpi::my_rank==0){
      #else
      {
      #endif
         thermal_temperature_file.open("spin-acc/thermal_temperature_profile.dat");
         // Header: time step, then Te and Tp for each cell (stack 0 only)
         thermal_temperature_file << "# step\t";
         const int ncz = num_microcells_per_stack;
         for(int cell = 0; cell < ncz; ++cell){
            thermal_temperature_file << "Te_" << cell << "\t" << "Tp_" << cell << "\t";
         }
         thermal_temperature_file << "\n";
      }
      
      return;
   }

   //-----------------------------------------------------------------------------
   // Function to write temperature profile data
   //-----------------------------------------------------------------------------
   void write_thermal_temperature_data(){
      if(!sc1d_thermal_gradients_enable || !ltmp::is_enabled()) return;
      if(!thermal_output_initialised) return;
      
      if(sim::time % ST_output_rate != 0) return;
      
      #ifdef MPICF
      if(vmpi::my_rank==0){
      #else
      {
      #endif
         const int ncz = num_microcells_per_stack;
         const double dz = micro_cell_thickness;
         
         thermal_temperature_file << thermal_output_counter << "\t";
         
         for(int cell = 0; cell < ncz; ++cell){
            const double zA = (cell + 0.5) * dz;
            const double Te = ltmp::get_electron_temperature_at_z(zA);
            const double Tp = ltmp::get_phonon_temperature_at_z(zA);
            thermal_temperature_file << std::fixed << std::setprecision(1);
            thermal_temperature_file << Te << "\t" << Tp << "\t";
         }
         
         thermal_temperature_file << std::endl;
      }
      
      thermal_output_counter++;
      
      return;
   }

   } // end of internal namespace

   void report_sc1d_profile(const double run_s){
      internal::report_sc1d_profile(run_s);
   }

   // Called before each LLG iteration. Writes the index-0 frame once, at the
   // first step of equilibration and the first step of the drive.
   void begin_sc1d_llg_stage_output(){
      if(!internal::sc1d_enable) return;
      if(sim::time < sim::equilibration_time) internal::output_sc1d_stage_origin("eq");
      else internal::output_sc1d_stage_origin("drive");
   }

} // end of st namespace
