//------------------------------------------------------------------------------
// fine_grid.cpp
//
// Build the 1D fine-grid transport constants from the atoms, not from the
// coarse microcell averages. Occupied cells keep the alloy mean. Empty cells
// are filled from those occupied anchors. Analytic broadening then spreads
// each material-file interface across the prolongation window. Diffuse uses
// the same interface constants, averaged by the atomic populations, and does
// not add that algebraic tail. Magnetisation is prolonged once at the start of
// relax from the Monte-Carlo-preconditioned atoms, held fixed through relax,
// then prolonged again on equilibration and drive steps.
//------------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "atoms.hpp"
#include "create.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

#include "co_bulk_4k_fit.hpp"
#include "internal.hpp"

namespace st {
namespace internal {

namespace {

static constexpr double kHbarFine = 1.05457162e-34;
static constexpr double kEFine = 1.60217662e-19;
static constexpr double kBoltzmannFine = 1.380649e-23;
static constexpr double kMetalNFine = 2.856e27; // Co DoS window n = g kT; match spincurrents.cpp
static constexpr double kEpsFine = 1e-30;

static double beta_c_seebeck_of(const mp_t& mat){
   if(mat.beta_seebeck_set) return mat.beta_seebeck;
   const double b = mat.beta_cond;
   return b + 0.5 * mat.seebeck_polarization * (1.0 - b * b);
}

static double einstein_sigma(const double D, const double Te){
   const double T = std::max(Te, 1.0);
   return (kMetalNFine * kEFine * kEFine * std::max(D, 0.0)) / (kBoltzmannFine * T);
}

static double blend_t(const double z, const double zL, const double zR, const double decay){
   const double span = zR - zL;
   if(!(span > 0.0)) return 0.0;
   const double zmid = 0.5 * (zL + zR);
   double s = 0.0;
   if(sc1d_prolongation == sc1d_prolong_step){
      return (z >= zmid) ? 1.0 : 0.0;
   }
   if(!(decay > 0.0) || 0.5 * span <= decay){
      s = (z - zL) / span;
   } else if(z <= zmid - decay){
      return 0.0;
   } else if(z >= zmid + decay){
      return 1.0;
   } else {
      s = (z - (zmid - decay)) / (2.0 * decay);
   }
   if(s < 0.0) s = 0.0;
   if(s > 1.0) s = 1.0;
   if(sc1d_prolongation == sc1d_prolong_smooth) s = s * s * (3.0 - 2.0 * s);
   return s;
}

static double kernel_shape(const double z, const double zif, const double half){
   if(!(half > 0.0)) return 0.0;
   const double x = std::fabs(z - zif) / half;
   if(!(x < 1.0)) return 0.0;
   if(sc1d_prolongation == sc1d_prolong_linear) return 1.0 - x;
   const double t = 1.0 - x;
   return t * t * (3.0 - 2.0 * t);
}

static double hermite(const double s){
   if(!(s > 0.0)) return 0.0;
   if(!(s < 1.0)) return 1.0;
   return s * s * (3.0 - 2.0 * s);
}

// Smoothstep coefficients, built once. A cell outside every ±decay window keeps
// the blend already written: that blend is the nearest anchor, which is what the
// old overlap wrote there. block[i] != 0 is not a prolongation target.
struct StackSmooth {
   std::vector<int> cell;
   std::vector<int> term_begin;
   std::vector<int> left;
   std::vector<int> right;
   std::vector<double> cL;
   std::vector<double> cR;
   std::vector<double> w;
};

// Overlapping ±decay windows are averaged. Each window is a convex combination
// of its two anchors, so the average cannot exceed the larger anchor: 7/10/7
// stays at or below 10, and 10/7/10 can lift the middle without passing 10.
static StackSmooth make_smooth_cache(const std::vector<int>& anchors,
                                     const int nf,
                                     const double dzA,
                                     const double decay,
                                     const char* block){
   StackSmooth cache;
   if(anchors.size() < 2u || nf <= 0 || !(decay > 0.0) || !(dzA > 0.0)) return cache;
   struct Term { int left; int right; double cL; double cR; double w; };
   std::vector<std::vector<Term>> terms(static_cast<std::size_t>(nf));
   for(std::size_t a = 0; a + 1u < anchors.size(); ++a){
      const double zL = (static_cast<double>(anchors[a]) + 0.5) * dzA;
      const double zR = (static_cast<double>(anchors[a + 1u]) + 0.5) * dzA;
      const double zif = 0.5 * (zL + zR);
      for(int i = 0; i < nf; ++i){
         const double z = (static_cast<double>(i) + 0.5) * dzA;
         if(z >= zif + decay) break;
         if(!(std::fabs(z - zif) < decay)) continue;
         if(block != nullptr && block[i] != 0) continue;
         const double s = (z - (zif - decay)) / (2.0 * decay);
         const double h = hermite(s);
         const double w = kernel_shape(z, zif, decay);
         if(!(w > 0.0)) continue;
         Term t;
         t.left = anchors[a];
         t.right = anchors[a + 1u];
         t.cL = w * (1.0 - h);
         t.cR = w * h;
         t.w = w;
         terms[static_cast<std::size_t>(i)].push_back(t);
      }
   }
   cache.term_begin.push_back(0);
   for(int i = 0; i < nf; ++i){
      const std::vector<Term>& row = terms[static_cast<std::size_t>(i)];
      if(row.empty()) continue;
      cache.cell.push_back(i);
      for(std::size_t k = 0; k < row.size(); ++k){
         cache.left.push_back(row[k].left);
         cache.right.push_back(row[k].right);
         cache.cL.push_back(row[k].cL);
         cache.cR.push_back(row[k].cR);
         cache.w.push_back(row[k].w);
      }
      cache.term_begin.push_back(static_cast<int>(cache.left.size()));
   }
   return cache;
}

static std::vector<double> smooth_src;

static void apply_smooth_cache(double* field, const int nf, const StackSmooth& cache){
   if(field == nullptr || nf <= 0 || cache.cell.empty()) return;
   smooth_src.resize(static_cast<std::size_t>(nf));
   std::copy(field, field + nf, smooth_src.begin());
   for(std::size_t c = 0; c < cache.cell.size(); ++c){
      const int i = cache.cell[c];
      if(i < 0 || i >= nf) continue;
      double value = 0.0;
      double weight = 0.0;
      const int b0 = cache.term_begin[c];
      const int b1 = cache.term_begin[c + 1u];
      for(int t = b0; t < b1; ++t){
         const std::size_t u = static_cast<std::size_t>(t);
         const int iL = cache.left[u];
         const int iR = cache.right[u];
         const double sL = (iL >= 0 && iL < nf) ? smooth_src[static_cast<std::size_t>(iL)] : 0.0;
         const double sR = (iR >= 0 && iR < nf) ? smooth_src[static_cast<std::size_t>(iR)] : 0.0;
         value += cache.cL[u] * sL + cache.cR[u] * sR;
         weight += cache.w[u];
      }
      if(weight > 0.0) field[i] = value / weight;
   }
}

// Occupied cells of a sigma = 0 material, and empty cells inside that
// material's height. Both are the insulator, so neither is a prolongation target.
static bool sigma_blocks_prolongation(const std::size_t slot, const int fine, const double dzA){
   const bool occupied = (sc1d_cell_mask[slot] & sc1d_mask_occupied) != 0;
   if(occupied){
      const int mat = sc1d_mat_fine[slot];
      return mat >= 0 && mat < static_cast<int>(mp.size())
         && mp[static_cast<std::size_t>(mat)].conductivity == 0.0;
   }
   const double z = (static_cast<double>(fine) + 0.5) * dzA;
   const double film_z = cs::system_dimensions[2];
   for(int mat = 0; mat < static_cast<int>(mp.size()); ++mat){
      const mp_t& p = mp[static_cast<std::size_t>(mat)];
      if(p.conductivity != 0.0) continue;
      const double z_lo = create::get_material_height_min(mat) * film_z;
      const double z_hi = create::get_material_height_max(mat) * film_z;
      if(z_hi > z_lo && z >= z_lo && z < z_hi) return true;
   }
   return false;
}

// Smoothstep and the gap fill both leave sigma = 0 untouched. Empty cells in
// the slab still need the material's own constants, which the blend never wrote.
static void stamp_zero_sigma(const std::size_t slot, const int mat){
   if(mat < 0 || mat >= static_cast<int>(mp.size())) return;
   const mp_t& p = mp[static_cast<std::size_t>(mat)];
   if(p.conductivity != 0.0) return;
   sc1d_Bc_fine[slot] = p.beta_cond;
   sc1d_Bd_fine[slot] = p.beta_diff;
   sc1d_sa_inf_fine[slot] = p.sa_infinity;
   sc1d_lsf_fine[slot] = p.lambda_sdl;
   sc1d_D_fine[slot] = p.diffusion;
   sc1d_Dmat_fine[slot] = p.diffusion;
   sc1d_Jsd_fine[slot] = p.sd_exchange;
   sc1d_lphi_fine[slot] = p.lambda_phi;
   sc1d_chi_fine[slot] = p.chi_demag;
   sc1d_seebeck_fine[slot] = p.seebeck_coefficient;
   sc1d_beta_seebeck_fine[slot] = beta_c_seebeck_of(p);
   sc1d_phot_fine[slot] = p.superdiffusion_polarization;
   sc1d_PN_fine[slot] = p.spin_dos_polarization;
   sc1d_alpha_fine[slot] = p.spin_dos_alpha;
   sc1d_sigma_fine[slot] = 0.0;
   sc1d_sigma_bulk_fine[slot] = 0.0;
}

static void require_finite(const char* name, const std::vector<double>& field){
   for(std::size_t i = 0; i < field.size(); ++i){
      if(std::isfinite(field[i])) continue;
      terminaltextcolor(RED);
      std::cerr << "Error - non-finite " << name << " on the fine grid at slot " << i
                << " = " << field[i] << ". This is not capped." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error - non-finite " << name << " at fine slot " << i
           << " = " << field[i] << "." << std::endl;
      err::vexit();
   }
}

// Series addition of one interface density onto the prolonged bulk. acc_up and
// acc_dn are resistivities (ohm m). acc_re and acc_im are S/m^3. acc_d2 is 1/m^2.
// The integral of each density over the window is the lumped interface value.
static void apply_sheet_density(const int i,
                                const double acc_up,
                                const double acc_dn,
                                const double acc_re,
                                const double acc_im,
                                const double acc_d2,
                                const double beta_b,
                                const double beta_d,
                                const double D_b,
                                const double lsf_b,
                                const double lphi_b,
                                const double Jsd_b,
                                const double sigma_b,
                                double* Bc,
                                double* D,
                                double* lsf,
                                double* lphi,
                                double* Jsd,
                                double* sigma){
   const double g_e = 2.0;
   double sig = sigma_b;
   double beta = beta_b;
   if(sigma_b > 0.0){
      const double sup_b = sigma_b * (1.0 + beta_b) / 2.0;
      const double sdn_b = sigma_b * (1.0 - beta_b) / 2.0;
      double sup = sup_b;
      double sdn = sdn_b;
      if(sup_b > 0.0 && acc_up > 0.0) sup = 1.0 / (1.0 / sup_b + acc_up);
      if(sdn_b > 0.0 && acc_dn > 0.0) sdn = 1.0 / (1.0 / sdn_b + acc_dn);
      sig = sup + sdn;
      if(sig > 0.0) beta = (sup - sdn) / sig;
      if(beta < 0.0) beta = 0.0;
   }

   double Diff = D_b;
   if(sigma_b > 0.0) Diff = D_b * sig / sigma_b;
   else Diff = 0.0;

   const double prod = spin_gamma(beta, beta_d);
   const double denom = 1.0 - prod;
   const double BBp = 1.0 / std::sqrt(denom);
   double inv_sf = acc_d2;
   if(lsf_b > 0.0){
      const double lam = lsf_b * BBp;
      inv_sf += 1.0 / (lam * lam);
   }
   double length_sf = lsf_b;
   if(inv_sf > 0.0 && BBp > 0.0) length_sf = 1.0 / (std::sqrt(inv_sf) * BBp);

   double inv_phi = 0.0;
   if(lphi_b > 0.0) inv_phi = 1.0 / (lphi_b * lphi_b);
   if(sigma_b > 0.0) inv_phi += g_e * acc_re / sigma_b;
   if(!(inv_phi >= 0.0) || !std::isfinite(inv_phi)){
      terminaltextcolor(RED);
      std::cerr << "Error - real mixing conductance made 1/lambda_phi^2 = " << inv_phi
                << " at fine cell " << i << ". A negative dephasing rate is not a length." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error - 1/lambda_phi^2 = " << inv_phi << " at fine cell " << i << "." << std::endl;
      err::vexit();
   }
   const double length_phi = (inv_phi > 0.0) ? (1.0 / std::sqrt(inv_phi)) : 0.0;

   double exchange = Jsd_b;
   if(sigma_b > 0.0) exchange += 2.0 * kHbarFine * Diff * (g_e * acc_im / sigma_b);

   Bc[i] = beta;
   D[i] = Diff;
   lsf[i] = length_sf;
   if(lphi != nullptr) lphi[i] = length_phi;
   if(Jsd != nullptr) Jsd[i] = exchange;
   sigma[i] = sig;
}

struct SheetMix {
   double gup;
   double gdn;
   double re;
   double im;
   double delta;
   double rint;
   bool sheet;
};

struct TransportSlot {
   double bc;
   double bd;
   double sa;
   double lsf;
   double diff;
   double jsd;
   double pn;
   double alpha;
   double sigma;
};

static TransportSlot material_transport(const mp_t& p,
                                        const double Te,
                                        const double Tp,
                                        const double ml){
   TransportSlot q;
   q.bc = p.beta_cond;
   q.bd = p.beta_diff;
   q.sa = p.sa_infinity;
   q.lsf = p.lambda_sdl;
   q.diff = p.diffusion;
   q.jsd = p.sd_exchange;
   q.pn = p.spin_dos_polarization;
   q.alpha = p.spin_dos_alpha;
   const bool thermal = sc1d_thermal_effects && p.spin_currents_thermal;
   if(thermal){
      const co_bulk_4k_fit::Result fit = co_bulk_4k_fit::evaluate(Te, Tp, ml);
      const double beta_c0 = (co_bulk_4k_fit::kSigmaUp0 - co_bulk_4k_fit::kSigmaDown0)
         / co_bulk_4k_fit::kSigma0;
      const double beta_d0 = (co_bulk_4k_fit::kDUp0 - co_bulk_4k_fit::kDDown0)
         / (co_bulk_4k_fit::kDUp0 + co_bulk_4k_fit::kDDown0);
      if(std::fabs(beta_c0) > 1.0e-30) q.bc = p.beta_cond * fit.beta_c / beta_c0;
      if(q.bc < 0.0) q.bc = 0.0;
      if(std::fabs(beta_d0) > 1.0e-30) q.bd = p.beta_diff * fit.beta_d / beta_d0;
      q.diff = p.diffusion * fit.diffusion_ratio;
      q.lsf = p.lambda_sdl * fit.lambda_sdl_ratio;
      q.pn = p.spin_dos_polarization * fit.pn_ratio;
      q.alpha = p.spin_dos_alpha * fit.alpha_ratio;
      // File constants are the m_l = 1 values. The fine cell supplies m_l.
      // Stiffness is m_l^1.7. J_sd is not a polynomial and is not replaced by
      // a Curie-Bloch m(Ts): the spin equation uses J_sd/(2 hbar), so m_l = 0
      // is a zero rate.
      const double rel_m = std::max(ml / co_bulk_4k_fit::kM0, 0.0);
      const double stiffness = std::pow(rel_m, 1.7);
      q.sa = p.sa_infinity * stiffness * (Te / co_bulk_4k_fit::kAnchorT);
      q.jsd = p.sd_exchange * stiffness;
      if(p.conductivity > 0.0) q.sigma = p.conductivity * fit.sigma_ratio;
      else if(p.conductivity < 0.0) q.sigma = einstein_sigma(q.diff, Te);
      else q.sigma = 0.0;
   } else {
      const double sigma_T = co_bulk_4k_fit::kAnchorT;
      if(p.conductivity > 0.0) q.sigma = p.conductivity;
      else if(p.conductivity < 0.0) q.sigma = einstein_sigma(q.diff, sigma_T);
      else q.sigma = 0.0;
   }
   return q;
}

// Directed pair from the lower composition to the upper one. A set sheet
// replaces that pair's lumped Robin resistance, as the coarse edge did.
static SheetMix mix_pair(const std::size_t slotL,
                         const std::size_t slotR,
                         const std::vector<double>& count,
                         const std::vector<double>& mat_count,
                         const std::size_t nslot,
                         const std::size_t nmat){
   SheetMix mix;
   mix.gup = 0.0;
   mix.gdn = 0.0;
   mix.re = 0.0;
   mix.im = 0.0;
   mix.delta = 0.0;
   mix.rint = 0.0;
   mix.sheet = false;
   const double n1 = count[slotL];
   const double n2 = count[slotR];
   if(n1 < 0.5 || n2 < 0.5 || nmat == 0) return mix;
   for(std::size_t a = 0; a < nmat; ++a){
      const double fa = mat_count[a * nslot + slotL] / n1;
      if(fa <= 0.0) continue;
      for(std::size_t b = 0; b < nmat; ++b){
         const double fb = mat_count[b * nslot + slotR] / n2;
         if(fb <= 0.0) continue;
         const std::size_t ip = a * nmat + b;
         const double w = fa * fb;
         const bool sheet = sheet_enabled && ip < sheet_pair_set.size() && sheet_pair_set[ip];
         if(sheet){
            mix.sheet = true;
            mix.gup += w * sheet_g_up_pair[ip];
            mix.gdn += w * sheet_g_down_pair[ip];
            mix.re += w * sheet_re_pair[ip];
            mix.im += w * sheet_im_pair[ip];
            mix.delta += w * sheet_delta_pair[ip];
         } else if(interface_coupling_enabled && ip < r_int_pair.size() && r_int_pair[ip] != 0.0){
            mix.rint += w * r_int_pair[ip];
         }
      }
   }
   return mix;
}

static void add_density(const SheetMix& mix,
                        const double sd,
                        double& acc_up,
                        double& acc_dn,
                        double& acc_re,
                        double& acc_im,
                        double& acc_d2){
   if(!mix.sheet || !(sd > 0.0)) return;
   if(mix.gup > 0.0) acc_up += sd / mix.gup;
   if(mix.gdn > 0.0) acc_dn += sd / mix.gdn;
   acc_re += sd * mix.re;
   acc_im += sd * mix.im;
   acc_d2 += (mix.delta * sd) * (mix.delta * sd);
}

// Flat window from the left occupied cell through the right one, including the
// empty cells between them and nothing outside. sd is chosen so the integral
// of the window is the full material-file interface, even where windows meet.
static void deposit_step_window(int iL,
                                int iR,
                                const SheetMix& mix,
                                const double dz_m,
                                const int nf,
                                double* rint,
                                std::vector<double>& acc_up,
                                std::vector<double>& acc_dn,
                                std::vector<double>& acc_re,
                                std::vector<double>& acc_im,
                                std::vector<double>& acc_d2){
   if(iR < iL) std::swap(iL, iR);
   if(iL < 0) iL = 0;
   if(iR >= nf) iR = nf - 1;
   const int ncell = iR - iL + 1;
   if(ncell <= 0 || !(dz_m > 0.0)) return;
   if(mix.sheet){
      const double sd = 1.0 / (static_cast<double>(ncell) * dz_m);
      for(int i = iL; i <= iR; ++i){
         add_density(mix, sd,
                     acc_up[static_cast<std::size_t>(i)], acc_dn[static_cast<std::size_t>(i)],
                     acc_re[static_cast<std::size_t>(i)], acc_im[static_cast<std::size_t>(i)],
                     acc_d2[static_cast<std::size_t>(i)]);
      }
   }
   if(mix.rint != 0.0 && rint != nullptr && iR > iL){
      const double share = mix.rint / static_cast<double>(iR - iL);
      for(int i = iL; i < iR; ++i) rint[i] += share;
   }
}

// lambda_sf = lambda_sdl / sqrt(1 - beta_c beta_d). The diffusion operator
// decays with this length. Memory loss shortens it directly.
static double lambda_sf_of(const double lsf, const double bc, const double bd){
   if(!(lsf > 0.0)) return 0.0;
   const double gap = 1.0 - spin_gamma(bc, bd);
   if(!(gap > 0.0)) return 0.0;
   return lsf / std::sqrt(gap);
}

static double lambda_gap_of(const double bc, const double bd){
   const double gap = 1.0 - spin_gamma(bc, bd);
   return (gap > 0.0) ? gap : 0.0;
}

// Hermite bump of one memory-loss interface. k is in [0, 1]. The half-width
// is at least four fine cells. The equivalent interface used for transmission
// is 3.5 A thick, independent of that window.
struct MemorySample {
   int cell;
   double k;
   double lam;
};

struct MemoryFace {
   double zif;
   double z_lo;
   double z_hi;
   int left;
   int right;
   double delta_file;
};

static std::vector<StackSmooth> smooth_cache;
static std::vector<std::vector<MemoryFace>> memory_faces;
// Hermite weight, fixed at build. Temperature changes the interior lengths
// this weight joins. The weight does not mix those interiors with the
// prolonged blend: that blend still carries the long-metal tail.
static std::vector<double> lambda_sf_hermite;
// Occupied-cell lambda_sf before prolongation. Interface loss is sized from
// these interiors and replaces the prolonged value inside the window.
static std::vector<double> lambda_sf_interior;
// Owning memory-loss face within the stack, or -1 outside every window.
static std::vector<int> lambda_sf_face;

static constexpr double kInterfaceThicknessA = 3.5;

static double memory_half_width(const double decay, const double dzA){
   const double four_cells = 4.0 * dzA;
   return (decay > four_cells) ? decay : four_cells;
}

static void fill_memory_window(const MemoryFace& face,
                               const double* lsf,
                               const double* bc,
                               const double* bd,
                               const int nf,
                               const double dzA,
                               const double half,
                               std::vector<MemorySample>& samples){
   samples.clear();
   for(int i = 0; i < nf; ++i){
      const double z = (static_cast<double>(i) + 0.5) * dzA;
      if(z >= face.zif + half) break;
      // Same ±half window as the prolongation. The stored interiors are kept
      // separately, so cells inside the window, including the anchors, can be
      // replaced.
      if(!(std::fabs(z - face.zif) < half)) continue;
      const double k = kernel_shape(z, face.zif, half);
      if(!(k > 0.0)) continue;
      const double lam = lambda_sf_of(lsf[i], bc[i], bd[i]);
      if(!(lam > 0.0)) continue;
      MemorySample sample;
      sample.cell = i;
      sample.k = k;
      sample.lam = lam;
      samples.push_back(sample);
   }
}

// Occupied cells still hold the material lambda_sf here. Prolongation and the
// sheet both run later and must not be what the interface is sized from.
static void snapshot_interior_lambda(){
   const std::size_t nslot = sc1d_lsf_fine.size();
   if(lambda_sf_interior.size() != nslot) lambda_sf_interior.assign(nslot, 0.0);
   for(std::size_t slot = 0; slot < nslot; ++slot){
      lambda_sf_interior[slot] = 0.0;
      if(slot >= sc1d_cell_mask.size() || slot >= sc1d_Bc_fine.size() || slot >= sc1d_Bd_fine.size()) continue;
      if((sc1d_cell_mask[slot] & sc1d_mask_occupied) == 0) continue;
      lambda_sf_interior[slot] = lambda_sf_of(sc1d_lsf_fine[slot], sc1d_Bc_fine[slot], sc1d_Bd_fine[slot]);
   }
}

// Center length from the two interior lambda_sf values. The material-file delta
// is the total spin-memory optical depth through an L-thick interface (order
// unity transmission exp(-delta)). Prolongation already attributes chi_pro =
// L/lambda_bar of that budget to diffusion at the mean interior length
// lambda_bar; the join centre supplies the remainder chi_c = delta - chi_pro
// via lambda_c = L/chi_c. The centre stays at or below the shorter interior
// so the cobalt side of a wide-spacer join cannot exceed the cobalt interior.
static double interface_center_length(const double lamL,
                                      const double lamR,
                                      const double thickness_m,
                                      const double delta_file){
   if(!(lamL > 0.0) || !(lamR > 0.0) || !(thickness_m > 0.0)) return 0.0;
   const double shorter = (lamL < lamR) ? lamL : lamR;
   const double lambda_bar = 0.5 * (lamL + lamR);
   if(!(lambda_bar > 0.0)) return 0.0;
   const double chi_pro = thickness_m / lambda_bar;
   const double chi_center = std::max(delta_file, 0.0) - chi_pro;
   if(!(chi_center > 0.0)) return shorter;
   const double eff = thickness_m / chi_center;
   return (eff < shorter) ? eff : shorter;
}

// Freeze which face owns each cell and the hermite weight. The lengths are
// applied later from the current interiors.
static void accumulate_lambda_delta(const int stack,
                                    const int nf,
                                    const double dzA,
                                    const double dz_m,
                                    const double decay){
   (void)dz_m;
   if(stack < 0 || static_cast<std::size_t>(stack) >= memory_faces.size()) return;
   const std::vector<MemoryFace>& faces = memory_faces[static_cast<std::size_t>(stack)];
   if(faces.empty()) return;
   const std::size_t base = static_cast<std::size_t>(stack) * static_cast<std::size_t>(nf);
   if(lambda_sf_hermite.size() < base + static_cast<std::size_t>(nf)) return;
   if(lambda_sf_face.size() < base + static_cast<std::size_t>(nf)) return;
   const double* lsf = sc1d_lsf_fine.data() + base;
   const double* bc = sc1d_Bc_fine.data() + base;
   const double* bd = sc1d_Bd_fine.data() + base;
   const double half = memory_half_width(decay, dzA);
   std::vector<MemorySample> samples;
   std::vector<int> owner(static_cast<std::size_t>(nf), -1);
   std::vector<double> nearest(static_cast<std::size_t>(nf), 1.0e99);
   for(std::size_t f = 0; f < faces.size(); ++f){
      fill_memory_window(faces[f], lsf, bc, bd, nf, dzA, half, samples);
      for(std::size_t s = 0; s < samples.size(); ++s){
         const int cell = samples[s].cell;
         const double z = (static_cast<double>(cell) + 0.5) * dzA;
         const double dist = std::fabs(z - faces[f].zif);
         if(!(dist < nearest[static_cast<std::size_t>(cell)])) continue;
         nearest[static_cast<std::size_t>(cell)] = dist;
         owner[static_cast<std::size_t>(cell)] = static_cast<int>(f);
      }
   }
   for(std::size_t f = 0; f < faces.size(); ++f){
      fill_memory_window(faces[f], lsf, bc, bd, nf, dzA, half, samples);
      for(std::size_t s = 0; s < samples.size(); ++s){
         const int cell = samples[s].cell;
         if(owner[static_cast<std::size_t>(cell)] != static_cast<int>(f)) continue;
         const std::size_t slot = base + static_cast<std::size_t>(cell);
         lambda_sf_hermite[slot] = samples[s].k;
         lambda_sf_face[slot] = static_cast<int>(f);
      }
   }
}

// Replace the prolonged lambda_sf inside the window. On each side the profile
// is (1 - k) * that side's interior + k * center. k is 1 at the interface and
// 0 at ±half, so the window meets the interiors and the long-metal prolongation
// does not leak through the shoulder.
static void apply_fixed_lambda_delta(){
   const int nf = sc1d_nf;
   const int nstack = std::max(0, num_stacks_y);
   if(nf <= 0 || nstack <= 0 || lambda_sf_hermite.empty() || lambda_sf_face.empty()) return;
   const std::size_t nslot = static_cast<std::size_t>(nstack) * static_cast<std::size_t>(nf);
   if(lambda_sf_hermite.size() < nslot || lambda_sf_face.size() < nslot) return;
   if(lambda_sf_interior.size() < nslot) return;
   if(sc1d_lsf_fine.size() < nslot || sc1d_Bc_fine.size() < nslot || sc1d_Bd_fine.size() < nslot) return;
   if(static_cast<int>(memory_faces.size()) < nstack) return;
   const int nsub = std::max(1, sc1d_nsub);
   const double dzA = micro_cell_thickness / static_cast<double>(nsub);
   if(!(dzA > 0.0)) return;
   const double thickness_m = kInterfaceThicknessA * 1.0e-10;

   for(int stack = 0; stack < nstack; ++stack){
      const std::vector<MemoryFace>& faces = memory_faces[static_cast<std::size_t>(stack)];
      const std::size_t base = static_cast<std::size_t>(stack) * static_cast<std::size_t>(nf);
      std::vector<double> center(faces.size(), 0.0);
      std::vector<double> side_left(faces.size(), 0.0);
      std::vector<double> side_right(faces.size(), 0.0);
      for(std::size_t f = 0; f < faces.size(); ++f){
         const int iL = faces[f].left;
         const int iR = faces[f].right;
         if(iL < 0 || iR < 0 || iL >= nf || iR >= nf) continue;
         const double lamL = lambda_sf_interior[base + static_cast<std::size_t>(iL)];
         const double lamR = lambda_sf_interior[base + static_cast<std::size_t>(iR)];
         side_left[f] = lamL;
         side_right[f] = lamR;
         center[f] = interface_center_length(lamL, lamR, thickness_m, faces[f].delta_file);
      }
      for(int i = 0; i < nf; ++i){
         const std::size_t slot = base + static_cast<std::size_t>(i);
         const int face = lambda_sf_face[slot];
         if(face < 0 || static_cast<std::size_t>(face) >= faces.size()) continue;
         double k = lambda_sf_hermite[slot];
         if(!(k > 0.0)) continue;
         if(k > 1.0) k = 1.0;
         const std::size_t uf = static_cast<std::size_t>(face);
         const double lamL = side_left[uf];
         const double lamR = side_right[uf];
         const double lamC = center[uf];
         if(!(lamL > 0.0) || !(lamR > 0.0) || !(lamC > 0.0)) continue;
         const double z = (static_cast<double>(i) + 0.5) * dzA;
         const double side = (z < faces[uf].zif) ? lamL : lamR;
         const double lam = (1.0 - k) * side + k * lamC;
         const double gap = lambda_gap_of(sc1d_Bc_fine[slot], sc1d_Bd_fine[slot]);
         if(!(lam > 0.0) || !(gap > 0.0)) continue;
         sc1d_lsf_fine[slot] = lam * std::sqrt(gap);
      }
   }
}

// Temperature replaces the occupied anchors. Every constant that temperature
// rewrote is smoothstepped again, including conductivity, diffusion and sa_inf.
static void reapply_smooth_bulk(){
   const int nf = sc1d_nf;
   const int nstack = std::max(0, num_stacks_y);
   if(nf <= 0 || nstack <= 0) return;
   if(sc1d_interface == sc1d_interface_diffuse) return;
   if(sc1d_prolongation != sc1d_prolong_smooth) return;
   if(static_cast<int>(smooth_cache.size()) != nstack) return;
   std::vector<double>* fields[8] = {
      &sc1d_Bc_fine, &sc1d_Bd_fine, &sc1d_D_fine, &sc1d_lsf_fine,
      &sc1d_sa_inf_fine, &sc1d_PN_fine, &sc1d_alpha_fine, &sc1d_sigma_fine
   };
   for(int stack = 0; stack < nstack; ++stack){
      const std::size_t base = static_cast<std::size_t>(stack) * static_cast<std::size_t>(nf);
      const StackSmooth& cache = smooth_cache[static_cast<std::size_t>(stack)];
      for(int f = 0; f < 8; ++f){
         if(fields[f]->size() < base + static_cast<std::size_t>(nf)) continue;
         apply_smooth_cache(fields[f]->data() + base, nf, cache);
      }
   }
}

static bool density_set(const double up, const double dn, const double re, const double im, const double d2){
   return up > 0.0 || dn > 0.0 || re != 0.0 || im != 0.0 || d2 > 0.0;
}

static int coarse_id_of(const double x, const double y, const double z){
   const int ncx = num_x_stacks;
   const int ncy = num_y_stacks;
   const int ncz = num_microcells_per_stack;
   const int d[3] = {ncx, ncy, ncz};
   const double cs[3] = {micro_cell_size[stx], micro_cell_size[sty], micro_cell_thickness};
   double c[3];
   c[stx] = x;
   c[sty] = y;
   c[stz] = z;
   int scc[3] = {0, 0, 0};
   for(int dim = 0; dim < 3; ++dim){
      scc[dim] = static_cast<int>(std::floor(c[dim] / cs[dim]));
      if(scc[dim] < 0 || scc[dim] >= d[dim]) return -1;
   }
   return ((scc[0] * ncy) + scc[1]) * ncz + scc[2];
}

static int fine_index_of(const double z, const double dzA, const int nf){
   int fine = static_cast<int>(std::floor(z / dzA));
   if(fine < 0) fine = 0;
   if(fine >= nf) fine = nf - 1;
   return fine;
}

} // namespace

void build_fine_material_grid(){
   const int nf = sc1d_nf;
   const int nstack = std::max(0, num_stacks_y);
   if(nf <= 0 || nstack <= 0) return;
   const int nsub = std::max(1, sc1d_nsub);
   const double dzA = micro_cell_thickness / static_cast<double>(nsub);
   const double dz_m = dzA * 1.0e-10;
   if(!(dzA > 0.0) || !(dz_m > 0.0)){
      terminaltextcolor(RED);
      std::cerr << "Error - fine-grid spacing is not positive." << std::endl;
      terminaltextcolor(WHITE);
      err::vexit();
   }

   double decay = sc1d_prolongation_decay;
   if(!(decay > 0.0)){
      decay = cs::unit_cell.dimensions[2];
      if(!(decay > 0.0)) decay = micro_cell_thickness;
      sc1d_prolongation_decay = decay;
   }

   const std::size_t nslot = static_cast<std::size_t>(nstack) * static_cast<std::size_t>(nf);
   const std::size_t nmat = mp.size();

   std::vector<double> count(nslot, 0.0);
   std::vector<double> bc(nslot, 0.0);
   std::vector<double> bd(nslot, 0.0);
   std::vector<double> sa(nslot, 0.0);
   std::vector<double> lsdl(nslot, 0.0);
   std::vector<double> diff(nslot, 0.0);
   std::vector<double> jsd(nslot, 0.0);
   std::vector<double> lphi(nslot, 0.0);
   std::vector<double> chi(nslot, 0.0);
   std::vector<double> seebeck(nslot, 0.0);
   std::vector<double> beta_s(nslot, 0.0);
   std::vector<double> phot(nslot, 0.0);
   std::vector<double> pn(nslot, 0.0);
   std::vector<double> alpha(nslot, 0.0);
   std::vector<double> sigma(nslot, 0.0);
   std::vector<double> msat(nslot, 0.0);
   std::vector<double> mat_count(nmat * nslot, 0.0);

   const int nconst = 16;
   auto add_material = [&](const int stack, const int fine, const int mat){
      if(mat < 0 || mat >= static_cast<int>(nmat)){
         terminaltextcolor(RED);
         std::cerr << "Error - spin-torque material parameters missing for material index " << mat << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      const std::size_t slot = static_cast<std::size_t>(stack) * static_cast<std::size_t>(nf) + static_cast<std::size_t>(fine);
      const mp_t& p = mp[static_cast<std::size_t>(mat)];
      bc[slot] += p.beta_cond;
      bd[slot] += p.beta_diff;
      sa[slot] += p.sa_infinity;
      lsdl[slot] += p.lambda_sdl;
      diff[slot] += p.diffusion;
      jsd[slot] += p.sd_exchange;
      lphi[slot] += p.lambda_phi;
      chi[slot] += p.chi_demag;
      seebeck[slot] += p.seebeck_coefficient;
      beta_s[slot] += beta_c_seebeck_of(p);
      phot[slot] += p.superdiffusion_polarization;
      pn[slot] += p.spin_dos_polarization;
      alpha[slot] += p.spin_dos_alpha;
      sigma[slot] += p.conductivity;
      if(mat < static_cast<int>(mp::mu_s_array.size())) msat[slot] += mp::mu_s_array[static_cast<std::size_t>(mat)];
      count[slot] += 1.0;
      mat_count[static_cast<std::size_t>(mat) * nslot + slot] += 1.0;
   };

   for(int atom = 0; atom < num_local_atoms; ++atom){
      if(atom >= static_cast<int>(sc1d_atom_ls.size()) || atom >= static_cast<int>(sc1d_atom_fine.size())) continue;
      const int stack = sc1d_atom_ls[static_cast<std::size_t>(atom)];
      const int fine = sc1d_atom_fine[static_cast<std::size_t>(atom)];
      if(stack < 0 || fine < 0 || fine >= nf) continue;
      add_material(stack, fine, atoms::type_array[static_cast<std::size_t>(atom)]);
   }

   for(std::size_t atom = 0; atom < cs::non_magnetic_atoms_array.size(); ++atom){
      const cs::nm_atom_t& nm = cs::non_magnetic_atoms_array[atom];
      const int id = coarse_id_of(nm.x, nm.y, nm.z + 0.0001);
      if(id < 0 || id >= static_cast<int>(cell_stack_index.size())){
         terminaltextcolor(RED);
         std::cerr << "Error - removed non-magnetic atom out of range in the fine-grid occupation." << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }
      const int stack = cell_stack_index[static_cast<std::size_t>(id)] - 1;
      if(stack < 0 || stack >= nstack) continue;
      add_material(stack, fine_index_of(nm.z, dzA, nf), nm.mat);
   }

#ifdef MPICF
   if(vmpi::num_processors > 1 && nslot > 0){
      std::vector<double> packed(static_cast<std::size_t>(nconst) * nslot + mat_count.size(), 0.0);
      std::vector<double>* parts[nconst] = {
         &count, &bc, &bd, &sa, &lsdl, &diff, &jsd, &lphi, &chi, &seebeck, &beta_s, &phot, &pn, &alpha, &sigma, &msat
      };
      std::size_t off = 0;
      for(int p = 0; p < nconst; ++p){
         std::copy(parts[p]->begin(), parts[p]->end(), packed.begin() + static_cast<std::ptrdiff_t>(off));
         off += nslot;
      }
      std::copy(mat_count.begin(), mat_count.end(), packed.begin() + static_cast<std::ptrdiff_t>(off));
      MPI_Allreduce(MPI_IN_PLACE, packed.data(), static_cast<int>(packed.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      off = 0;
      for(int p = 0; p < nconst; ++p){
         std::copy(packed.begin() + static_cast<std::ptrdiff_t>(off),
                   packed.begin() + static_cast<std::ptrdiff_t>(off + nslot),
                   parts[p]->begin());
         off += nslot;
      }
      std::copy(packed.begin() + static_cast<std::ptrdiff_t>(off), packed.end(), mat_count.begin());
   }
#endif

   sc1d_cell_mask.assign(nslot, 0);
   sc1d_natom_fine.assign(nslot, 0.0);
   sc1d_msat_fine.assign(nslot, 0.0);
   sc1d_anchor_left.assign(nslot, 0);
   sc1d_anchor_right.assign(nslot, 0);
   sc1d_anchor_t.assign(nslot, 0.0);
   std::fill(sc1d_Rint_face.begin(), sc1d_Rint_face.end(), 0.0);
   std::fill(sc1d_sheet_acc_up_fine.begin(), sc1d_sheet_acc_up_fine.end(), 0.0);
   std::fill(sc1d_sheet_acc_dn_fine.begin(), sc1d_sheet_acc_dn_fine.end(), 0.0);
   std::fill(sc1d_sheet_acc_re_fine.begin(), sc1d_sheet_acc_re_fine.end(), 0.0);
   std::fill(sc1d_sheet_acc_im_fine.begin(), sc1d_sheet_acc_im_fine.end(), 0.0);
   std::fill(sc1d_sheet_acc_d2_fine.begin(), sc1d_sheet_acc_d2_fine.end(), 0.0);
   std::fill(sc1d_Rint_base_face.begin(), sc1d_Rint_base_face.end(), 0.0);
   sc1d_m_prev_fine.assign(nslot, 0.0);
   sc1d_m0_fine.assign(nslot, 0.0);
   sc1d_ml_fine.assign(nslot, 0.0);
   // Non-magnetic cells stay at zero length. A unit +z filler becomes a
   // positive sa_inf*m peak wherever the scalar sa_inf tail enters Ta or MgO.
   sc1d_mhat_ref_fine.assign(nslot * 3u, 0.0);

   int n_occupied = 0;
   sc1d_mix_fine.assign(nslot, std::vector<fine_mix_t>());
   for(std::size_t slot = 0; slot < nslot; ++slot){
      const int fine = static_cast<int>(slot % static_cast<std::size_t>(nf));
      if(!(count[slot] > 0.0)) continue;
      const double inv = 1.0 / count[slot];
      bc[slot] *= inv;
      bd[slot] *= inv;
      sa[slot] *= inv;
      lsdl[slot] *= inv;
      diff[slot] *= inv;
      jsd[slot] *= inv;
      lphi[slot] *= inv;
      chi[slot] *= inv;
      seebeck[slot] *= inv;
      beta_s[slot] *= inv;
      phot[slot] *= inv;
      pn[slot] *= inv;
      alpha[slot] *= inv;
      sigma[slot] *= inv;
         if(sigma[slot] < 0.0) sigma[slot] = einstein_sigma(diff[slot], sc1d_Te_fine[static_cast<std::size_t>(fine)]);
      sc1d_Bc_fine[slot] = bc[slot];
      sc1d_Bd_fine[slot] = bd[slot];
      sc1d_sa_inf_fine[slot] = sa[slot];
      sc1d_lsf_fine[slot] = lsdl[slot];
      sc1d_D_fine[slot] = diff[slot];
      sc1d_Jsd_fine[slot] = jsd[slot];
      sc1d_lphi_fine[slot] = lphi[slot];
      sc1d_chi_fine[slot] = chi[slot];
      sc1d_seebeck_fine[slot] = seebeck[slot];
      sc1d_beta_seebeck_fine[slot] = beta_s[slot];
      sc1d_phot_fine[slot] = phot[slot];
      sc1d_PN_fine[slot] = pn[slot];
      sc1d_alpha_fine[slot] = alpha[slot];
      sc1d_sigma_fine[slot] = sigma[slot];
      int best_mat = -1;
      double best_n = 0.0;
      for(std::size_t mat_i = 0; mat_i < nmat; ++mat_i){
         const double nc = mat_count[mat_i * nslot + slot];
         if(nc > best_n){
            best_n = nc;
            best_mat = static_cast<int>(mat_i);
         }
      }
      sc1d_mat_fine[slot] = best_mat;
      std::vector<fine_mix_t> mix;
      for(std::size_t mat_i = 0; mat_i < nmat; ++mat_i){
         const double nc = mat_count[mat_i * nslot + slot];
         if(nc > 0.0) mix.push_back(fine_mix_t{static_cast<int>(mat_i), nc});
      }
      if(mix.size() > 1) sc1d_mix_fine[slot] = std::move(mix);
      sc1d_natom_fine[slot] = count[slot];
      sc1d_msat_fine[slot] = msat[slot];
      sc1d_cell_mask[slot] = static_cast<unsigned char>(sc1d_cell_mask[slot] | sc1d_mask_occupied);
      ++n_occupied;
   }
   if(n_occupied == 0){
      terminaltextcolor(RED);
      std::cerr << "Error - the fine grid has no occupied cells, so transport constants cannot be prolonged." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error - fine grid occupation is empty." << std::endl;
      err::vexit();
   }

   // The field window is per atom, hence per column. Each rank paints only the
   // atoms it holds. The column mark is the union of those windows.
   std::vector<int> connected(nslot, 0);
   for(int atom = 0; atom < num_local_atoms; ++atom){
      if(atom >= static_cast<int>(sc1d_atom_ls.size())) continue;
      const int stack = sc1d_atom_ls[static_cast<std::size_t>(atom)];
      const int lo = sc1d_atom_fine_lo[static_cast<std::size_t>(atom)];
      const int hi = sc1d_atom_fine_hi[static_cast<std::size_t>(atom)];
      if(stack < 0 || lo < 0 || hi < lo) continue;
      for(int i = lo; i <= hi && i < nf; ++i){
         const std::size_t slot = static_cast<std::size_t>(stack) * static_cast<std::size_t>(nf) + static_cast<std::size_t>(i);
         connected[slot] = 1;
      }
   }
#ifdef MPICF
   if(vmpi::num_processors > 1 && !connected.empty()){
      MPI_Allreduce(MPI_IN_PLACE, connected.data(), static_cast<int>(connected.size()),
                    MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   }
#endif
   int n_connected = 0;
   int n_spacing = 0;
   for(std::size_t slot = 0; slot < nslot; ++slot){
      if(connected[slot] != 0){
         sc1d_cell_mask[slot] = static_cast<unsigned char>(sc1d_cell_mask[slot] | sc1d_mask_connected);
         ++n_connected;
      } else {
         sc1d_cell_mask[slot] = static_cast<unsigned char>(sc1d_cell_mask[slot] | sc1d_mask_spacing);
         ++n_spacing;
      }
   }

   std::vector<double>* fields[14] = {
      &sc1d_Bc_fine, &sc1d_Bd_fine, &sc1d_D_fine, &sc1d_lsf_fine, &sc1d_lphi_fine, &sc1d_Jsd_fine,
      &sc1d_chi_fine, &sc1d_sa_inf_fine, &sc1d_PN_fine, &sc1d_alpha_fine, &sc1d_sigma_fine,
      &sc1d_seebeck_fine, &sc1d_beta_seebeck_fine, &sc1d_phot_fine
   };
   smooth_cache.assign(static_cast<std::size_t>(nstack), StackSmooth());
   memory_faces.assign(static_cast<std::size_t>(nstack), std::vector<MemoryFace>());
   lambda_sf_hermite.assign(nslot, 0.0);
   lambda_sf_interior.assign(nslot, 0.0);
   lambda_sf_face.assign(nslot, -1);
   snapshot_interior_lambda();

   for(int stack = 0; stack < nstack; ++stack){
      const std::size_t base = static_cast<std::size_t>(stack) * static_cast<std::size_t>(nf);
      std::vector<int> anchors;
      anchors.reserve(static_cast<std::size_t>(nf));
      for(int i = 0; i < nf; ++i){
         if(sc1d_cell_mask[base + static_cast<std::size_t>(i)] & sc1d_mask_occupied) anchors.push_back(i);
      }
      if(anchors.empty()){
         for(int i = 0; i < nf; ++i){
            sc1d_anchor_left[base + static_cast<std::size_t>(i)] = 0;
            sc1d_anchor_right[base + static_cast<std::size_t>(i)] = 0;
            sc1d_anchor_t[base + static_cast<std::size_t>(i)] = 0.0;
         }
         continue;
      }
      for(int i = 0; i < nf; ++i){
         const std::size_t slot = base + static_cast<std::size_t>(i);
         if(sc1d_cell_mask[slot] & sc1d_mask_occupied){
            sc1d_anchor_left[slot] = i;
            sc1d_anchor_right[slot] = i;
            sc1d_anchor_t[slot] = 0.0;
            continue;
         }
         int L = anchors.front();
         int R = anchors.back();
         if(i < anchors.front()){
            R = anchors.front();
            L = R;
         } else if(i > anchors.back()){
            L = anchors.back();
            R = L;
         } else {
            for(std::size_t a = 0; a + 1u < anchors.size(); ++a){
               if(anchors[a] < i && i < anchors[a + 1u]){
                  L = anchors[a];
                  R = anchors[a + 1u];
                  break;
               }
            }
         }
         const double z = (static_cast<double>(i) + 0.5) * dzA;
         const double zL = (static_cast<double>(L) + 0.5) * dzA;
         const double zR = (static_cast<double>(R) + 0.5) * dzA;
         sc1d_anchor_left[slot] = L;
         sc1d_anchor_right[slot] = R;
         // Diffuse follows the atomic populations. Those are sharp on this
         // stack, so the blend is a step. Analytic uses the selected kernel.
         if(sc1d_interface == sc1d_interface_diffuse){
            sc1d_anchor_t[slot] = (z >= 0.5 * (zL + zR)) ? 1.0 : 0.0;
         } else {
            sc1d_anchor_t[slot] = blend_t(z, zL, zR, decay);
         }
      }
      std::vector<char> block(static_cast<std::size_t>(nf), 0);
      for(int i = 0; i < nf; ++i){
         if(sigma_blocks_prolongation(base + static_cast<std::size_t>(i), i, dzA)) block[static_cast<std::size_t>(i)] = 1;
      }
      for(int f = 0; f < 14; ++f){
         if(fields[f]->size() < base + static_cast<std::size_t>(nf)) continue;
         double* field = fields[f]->data() + base;
         for(int i = 0; i < nf; ++i){
            if(block[static_cast<std::size_t>(i)] != 0) continue;
            const std::size_t slot = base + static_cast<std::size_t>(i);
            const int L = sc1d_anchor_left[slot];
            const int R = sc1d_anchor_right[slot];
            const double t = sc1d_anchor_t[slot];
            if(L == i && R == i) continue;
            field[i] = (1.0 - t) * field[L] + t * field[R];
         }
      }
      if(sc1d_interface != sc1d_interface_diffuse && sc1d_prolongation == sc1d_prolong_smooth){
         smooth_cache[static_cast<std::size_t>(stack)] = make_smooth_cache(anchors, nf, dzA, decay, block.data());
         const StackSmooth& cache = smooth_cache[static_cast<std::size_t>(stack)];
         for(int f = 0; f < 14; ++f){
            if(fields[f]->size() < base + static_cast<std::size_t>(nf)) continue;
            apply_smooth_cache(fields[f]->data() + base, nf, cache);
         }
      }
   }

   // Kept for the hot-electron speed. Charge diffusion uses sc1d_D_fine, which the
   // sheet then scales with σ↑+σ↓. The up/down imbalance is beta_c, not this copy.
   std::copy(sc1d_D_fine.begin(), sc1d_D_fine.begin() + static_cast<std::ptrdiff_t>(nslot), sc1d_Dmat_fine.begin());
   sc1d_sigma_bulk_fine = sc1d_sigma_fine;

   for(int stack = 0; stack < nstack; ++stack){
      const std::size_t base = static_cast<std::size_t>(stack) * static_cast<std::size_t>(nf);
      std::vector<int> anchors;
      for(int i = 0; i < nf; ++i){
         if(sc1d_cell_mask[base + static_cast<std::size_t>(i)] & sc1d_mask_occupied) anchors.push_back(i);
      }
      std::vector<double> acc_up(static_cast<std::size_t>(nf), 0.0);
      std::vector<double> acc_dn(static_cast<std::size_t>(nf), 0.0);
      std::vector<double> acc_re(static_cast<std::size_t>(nf), 0.0);
      std::vector<double> acc_im(static_cast<std::size_t>(nf), 0.0);
      std::vector<double> acc_d2(static_cast<std::size_t>(nf), 0.0);

      if(sc1d_interface == sc1d_interface_diffuse){
         for(std::size_t a = 0; a < anchors.size(); ++a){
            const int i = anchors[a];
            const std::size_t slot = base + static_cast<std::size_t>(i);
            std::vector<int> mats;
            for(std::size_t m = 0; m < nmat; ++m){
               if(mat_count[m * nslot + slot] > 0.0) mats.push_back(static_cast<int>(m));
            }
            if(mats.size() < 2u) continue;
            SheetMix mix;
            mix.gup = mix.gdn = mix.re = mix.im = mix.delta = mix.rint = 0.0;
            mix.sheet = false;
            const double ncell = count[slot];
            for(std::size_t ia = 0; ia < mats.size(); ++ia){
               for(std::size_t ib = ia + 1u; ib < mats.size(); ++ib){
                  const std::size_t ma = static_cast<std::size_t>(mats[ia]);
                  const std::size_t mb = static_cast<std::size_t>(mats[ib]);
                  const double fa = mat_count[ma * nslot + slot] / ncell;
                  const double fb = mat_count[mb * nslot + slot] / ncell;
                  std::size_t ip = ma * nmat + mb;
                  bool set = sheet_enabled && ip < sheet_pair_set.size() && sheet_pair_set[ip];
                  if(!set){
                     ip = mb * nmat + ma;
                     set = sheet_enabled && ip < sheet_pair_set.size() && sheet_pair_set[ip];
                  }
                  if(!set) continue;
                  const double w = fa * fb;
                  mix.sheet = true;
                  mix.gup += w * sheet_g_up_pair[ip];
                  mix.gdn += w * sheet_g_down_pair[ip];
                  mix.re += w * sheet_re_pair[ip];
                  mix.im += w * sheet_im_pair[ip];
                  mix.delta += w * sheet_delta_pair[ip];
               }
            }
            add_density(mix, 1.0 / dz_m, acc_up[static_cast<std::size_t>(i)], acc_dn[static_cast<std::size_t>(i)],
                        acc_re[static_cast<std::size_t>(i)], acc_im[static_cast<std::size_t>(i)], acc_d2[static_cast<std::size_t>(i)]);
         }
      }
      if(anchors.size() >= 2u){
         const int nface = static_cast<int>(anchors.size()) - 1;
         const bool step_window = (sc1d_interface == sc1d_interface_diffuse)
            || (sc1d_prolongation == sc1d_prolong_step);
         for(int face = 0; face < nface; ++face){
            const int iL = anchors[static_cast<std::size_t>(face)];
            const int iR = anchors[static_cast<std::size_t>(face + 1)];
            const SheetMix mix = mix_pair(base + static_cast<std::size_t>(iL), base + static_cast<std::size_t>(iR),
                                          count, mat_count, nslot, nmat);
            const double zL = (static_cast<double>(iL) + 0.5) * dzA;
            const double zR = (static_cast<double>(iR) + 0.5) * dzA;
            const double zif = 0.5 * (zL + zR);
            // Step, and every diffuse case, puts the sheet on the cells from the
            // left occupied plane through the right one. Nothing outside that
            // pair is averaged. Diffuse uses this rather than the algebraic
            // kernel: the curve has to come from the atomic populations.
            if(step_window){
               deposit_step_window(iL, iR, mix, dz_m, nf,
                                   sc1d_Rint_face.data() + base,
                                   acc_up, acc_dn, acc_re, acc_im, acc_d2);
               continue;
            }
            // Smoothstep keeps the full ±decay window, including occupied planes
            // and the overlap with the next interface. The kernel is normalised
            // on that window, so the integrated conductance is still the full 1/G.
            // Memory loss is a later hermite mix of lambda_sf, not an acc_d2 sheet.
            // Linear stays clipped to the neighbouring interface midpoints.
            double z_lo = (face == 0) ? zL : 0.5 * ((static_cast<double>(anchors[static_cast<std::size_t>(face - 1)]) + 0.5) * dzA + zL);
            double z_hi = (face + 1 == nface) ? zR : 0.5 * (zR + (static_cast<double>(anchors[static_cast<std::size_t>(face + 2)]) + 0.5) * dzA);
            const bool wide = (sc1d_prolongation == sc1d_prolong_smooth && decay > 0.0);
            if(wide){
               z_lo = zif - decay;
               z_hi = zif + decay;
            } else if(sc1d_prolongation != sc1d_prolong_step && decay > 0.0){
               z_lo = std::max(z_lo, zif - decay);
               z_hi = std::min(z_hi, zif + decay);
            }
            if(!(z_hi > z_lo)) continue;

            std::vector<int> cells;
            std::vector<double> shape;
            double shape_sum = 0.0;
            const bool last_face = (face + 1 == nface);
            for(int i = 0; i < nf; ++i){
               const double z = (static_cast<double>(i) + 0.5) * dzA;
               const bool inside = (wide || last_face) ? (z >= z_lo && z <= z_hi) : (z >= z_lo && z < z_hi);
               if(!inside) continue;
               const double side = (z >= zif) ? (z_hi - zif) : (zif - z_lo);
               const double s = kernel_shape(z, zif, side);
               if(!(s > 0.0)) continue;
               cells.push_back(i);
               shape.push_back(s);
               shape_sum += s;
            }
            if(cells.empty()){
               int nearest = iL;
               double best_d = 1.0e99;
               for(int i = 0; i < nf; ++i){
                  const double d = std::fabs((static_cast<double>(i) + 0.5) * dzA - zif);
                  if(d < best_d){
                     best_d = d;
                     nearest = i;
                  }
               }
               cells.push_back(nearest);
               shape.push_back(1.0);
               shape_sum = 1.0;
            }
            if(mix.sheet && shape_sum > 0.0){
               const bool match_memory = wide && mix.delta > 0.0;
               SheetMix cond = mix;
               if(match_memory) cond.delta = 0.0;
               if(match_memory){
                  MemoryFace recorded;
                  recorded.zif = zif;
                  recorded.z_lo = (face == 0) ? (zL - decay) : ((static_cast<double>(anchors[static_cast<std::size_t>(face - 1)]) + 0.5) * dzA);
                  recorded.z_hi = (face + 1 == nface) ? (zR + decay) : ((static_cast<double>(anchors[static_cast<std::size_t>(face + 2)]) + 0.5) * dzA);
                  recorded.left = iL;
                  recorded.right = iR;
                  recorded.delta_file = mix.delta;
                  memory_faces[static_cast<std::size_t>(stack)].push_back(recorded);
               }
               for(std::size_t c = 0; c < cells.size(); ++c){
                  const int i = cells[c];
                  const double sd = shape[c] / (shape_sum * dz_m);
                  add_density(cond, sd, acc_up[static_cast<std::size_t>(i)], acc_dn[static_cast<std::size_t>(i)],
                              acc_re[static_cast<std::size_t>(i)], acc_im[static_cast<std::size_t>(i)], acc_d2[static_cast<std::size_t>(i)]);
               }
            }
            if(mix.rint != 0.0){
               std::vector<int> faces;
               std::vector<double> fshape;
               double fsum = 0.0;
               for(int i = 0; i < nf - 1; ++i){
                  const double zface = static_cast<double>(i + 1) * dzA;
                  const bool inside = (wide || last_face) ? (zface >= z_lo && zface <= z_hi) : (zface >= z_lo && zface < z_hi);
                  if(!inside) continue;
                  const double side = (zface >= zif) ? (z_hi - zif) : (zif - z_lo);
                  const double s = kernel_shape(zface, zif, side);
                  if(!(s > 0.0)) continue;
                  faces.push_back(i);
                  fshape.push_back(s);
                  fsum += s;
               }
               if(faces.empty()){
                  int best = std::min(iL, nf - 2);
                  if(best < 0) best = 0;
                  faces.push_back(best);
                  fshape.push_back(1.0);
                  fsum = 1.0;
               }
               for(std::size_t f = 0; f < faces.size(); ++f){
                  sc1d_Rint_face[base + static_cast<std::size_t>(faces[f])] += mix.rint * fshape[f] / fsum;
               }
            }
         }
      }

      for(int i = 0; i < nf; ++i){
         const std::size_t slot = base + static_cast<std::size_t>(i);
         sc1d_sheet_acc_up_fine[slot] = acc_up[static_cast<std::size_t>(i)];
         sc1d_sheet_acc_dn_fine[slot] = acc_dn[static_cast<std::size_t>(i)];
         sc1d_sheet_acc_re_fine[slot] = acc_re[static_cast<std::size_t>(i)];
         sc1d_sheet_acc_im_fine[slot] = acc_im[static_cast<std::size_t>(i)];
         sc1d_sheet_acc_d2_fine[slot] = acc_d2[static_cast<std::size_t>(i)];
      }

      double* Bc = sc1d_Bc_fine.data() + base;
      double* Bd = sc1d_Bd_fine.data() + base;
      double* D = sc1d_D_fine.data() + base;
      double* Dmat = sc1d_Dmat_fine.data() + base;
      double* lsf = sc1d_lsf_fine.data() + base;
      double* lphi = sc1d_lphi_fine.data() + base;
      double* Jsd = sc1d_Jsd_fine.data() + base;
      double* sigma_b = sc1d_sigma_bulk_fine.data() + base;
      double* sigma_f = sc1d_sigma_fine.data() + base;
      std::vector<double> Bc0(Bc, Bc + nf);
      std::vector<double> lsf0(lsf, lsf + nf);
      std::vector<double> lphi0(lphi, lphi + nf);
      std::vector<double> Jsd0(Jsd, Jsd + nf);
      for(int i = 0; i < nf; ++i){
         if(!density_set(acc_up[static_cast<std::size_t>(i)], acc_dn[static_cast<std::size_t>(i)],
                         acc_re[static_cast<std::size_t>(i)], acc_im[static_cast<std::size_t>(i)],
                         acc_d2[static_cast<std::size_t>(i)])) continue;
         // An explicit insulator stays at zero conductivity. The sheet is a metallic series
         // resistance and is not smeared into that cell.
         if(!(sigma_b[i] > 0.0)) continue;
         apply_sheet_density(i, acc_up[static_cast<std::size_t>(i)], acc_dn[static_cast<std::size_t>(i)],
                             acc_re[static_cast<std::size_t>(i)], acc_im[static_cast<std::size_t>(i)],
                             acc_d2[static_cast<std::size_t>(i)],
                             Bc0[static_cast<std::size_t>(i)], Bd[i], Dmat[i], lsf0[static_cast<std::size_t>(i)],
                             lphi0[static_cast<std::size_t>(i)], Jsd0[static_cast<std::size_t>(i)], sigma_b[i],
                             Bc, D, lsf, lphi, Jsd, sigma_f);
      }
      accumulate_lambda_delta(stack, nf, dzA, dz_m, decay);
   }
   apply_fixed_lambda_delta();
   // sigma = 0 is not prolonged. Occupied cells already kept their material
   // values; this stamps the same constants onto any cell the blend did touch.
   // A conductor with diffusion 0 still has only D cleared.
   for(std::size_t slot = 0; slot < nslot; ++slot){
      const int mat = sc1d_mat_fine[slot];
      if(mat < 0 || mat >= static_cast<int>(mp.size())) continue;
      if(mp[static_cast<std::size_t>(mat)].conductivity == 0.0){
         stamp_zero_sigma(slot, mat);
      } else if(mp[static_cast<std::size_t>(mat)].diffusion == 0.0){
         sc1d_D_fine[slot] = 0.0;
         sc1d_Dmat_fine[slot] = 0.0;
      }
   }
   const double film_z = cs::system_dimensions[2];
   for(int mat = 0; mat < static_cast<int>(mp.size()); ++mat){
      const mp_t& p = mp[static_cast<std::size_t>(mat)];
      const bool kill_sigma = p.conductivity == 0.0;
      const bool kill_diff = p.diffusion == 0.0;
      if(!kill_sigma && !kill_diff) continue;
      const double z_lo = create::get_material_height_min(mat) * film_z;
      const double z_hi = create::get_material_height_max(mat) * film_z;
      if(!(z_hi > z_lo)) continue;
      for(std::size_t slot = 0; slot < nslot; ++slot){
         const int fine = static_cast<int>(slot % static_cast<std::size_t>(nf));
         const double z = (static_cast<double>(fine) + 0.5) * dzA;
         if(z < z_lo || !(z < z_hi)) continue;
         const bool occupied = (sc1d_cell_mask[slot] & sc1d_mask_occupied) != 0;
         if(occupied && sc1d_mat_fine[slot] != mat) continue;
         if(!occupied) sc1d_mat_fine[slot] = mat;
         if(kill_sigma){
            stamp_zero_sigma(slot, mat);
         } else if(kill_diff){
            sc1d_D_fine[slot] = 0.0;
            sc1d_Dmat_fine[slot] = 0.0;
         }
      }
   }
   (void)nconst;
   sc1d_Rint_base_face = sc1d_Rint_face;

   require_finite("beta_c", sc1d_Bc_fine);
   require_finite("beta_d", sc1d_Bd_fine);
   require_finite("D", sc1d_D_fine);
   require_finite("De", sc1d_Dmat_fine);
   require_finite("lambda_sdl", sc1d_lsf_fine);
   require_finite("lambda_phi", sc1d_lphi_fine);
   require_finite("Jsd", sc1d_Jsd_fine);
   require_finite("chi", sc1d_chi_fine);
   require_finite("sa_inf", sc1d_sa_inf_fine);
   require_finite("P_N", sc1d_PN_fine);
   require_finite("alpha", sc1d_alpha_fine);
   require_finite("sigma", sc1d_sigma_fine);
   require_finite("sigma_bulk", sc1d_sigma_bulk_fine);
   require_finite("seebeck", sc1d_seebeck_fine);
   require_finite("beta_seebeck", sc1d_beta_seebeck_fine);
   require_finite("P_hot", sc1d_phot_fine);
   require_finite("Rint", sc1d_Rint_face);

#ifdef MPICF
   if(vmpi::my_rank == 0)
#endif
   {
      const char* shape = "smoothstep";
      if(sc1d_prolongation == sc1d_prolong_step) shape = "step";
      else if(sc1d_prolongation == sc1d_prolong_linear) shape = "linear";
      const char* face = (sc1d_interface == sc1d_interface_diffuse) ? "diffuse" : "analytic";
      zlog << zTs() << "Fine grid " << nf << " cells, " << nstack << " columns. Occupied " << n_occupied
           << ", connected " << n_connected << ", spacing " << n_spacing
           << ". Prolongation " << shape << ", decay " << sc1d_prolongation_decay
           << " A, interface " << face << ".";
      if(sc1d_interface == sc1d_interface_diffuse){
         zlog << " Diffuse averages the material-file interface by population and prolongs it as a step.";
      }
      zlog << std::endl;
   }
}

void refresh_fine_transport_constants(){
   const int nf = sc1d_nf;
   const int nstack = std::max(0, num_stacks_y);
   if(nf <= 0 || nstack <= 0) return;
   const std::size_t nslot = static_cast<std::size_t>(nstack) * static_cast<std::size_t>(nf);
   if(sc1d_Bc_fine.size() < nslot || sc1d_mat_fine.size() < nslot) return;

   for(std::size_t slot = 0; slot < nslot; ++slot){
      if(!(sc1d_cell_mask[slot] & sc1d_mask_occupied)) continue;
      const int fine = static_cast<int>(slot % static_cast<std::size_t>(nf));
      const double Te = sc1d_Te_fine[static_cast<std::size_t>(fine)];
      const double Tp = sc1d_Tp_fine[static_cast<std::size_t>(fine)];
      const double ml = (slot < sc1d_ml_fine.size()) ? sc1d_ml_fine[slot] : co_bulk_4k_fit::kM0;

      TransportSlot sum;
      sum.bc = sum.bd = sum.sa = sum.lsf = sum.diff = sum.jsd = sum.pn = sum.alpha = sum.sigma = 0.0;
      double count = 0.0;
      if(slot < sc1d_mix_fine.size() && !sc1d_mix_fine[slot].empty()){
         const std::vector<fine_mix_t>& mix = sc1d_mix_fine[slot];
         for(std::size_t k = 0; k < mix.size(); ++k){
            const int mat = mix[k].mat;
            if(mat < 0 || mat >= static_cast<int>(mp.size())) continue;
            const double w = mix[k].count;
            const TransportSlot q = material_transport(mp[static_cast<std::size_t>(mat)], Te, Tp, ml);
            sum.bc += w * q.bc;
            sum.bd += w * q.bd;
            sum.sa += w * q.sa;
            sum.lsf += w * q.lsf;
            sum.diff += w * q.diff;
            sum.jsd += w * q.jsd;
            sum.pn += w * q.pn;
            sum.alpha += w * q.alpha;
            sum.sigma += w * q.sigma;
            count += w;
         }
      } else {
         const int mat = sc1d_mat_fine[slot];
         if(mat >= 0 && mat < static_cast<int>(mp.size())){
            const TransportSlot q = material_transport(mp[static_cast<std::size_t>(mat)], Te, Tp, ml);
            sum = q;
            count = 1.0;
         }
      }
      if(!(count > 0.0)) continue;
      const double inv = 1.0 / count;
      sc1d_Bc_fine[slot] = sum.bc * inv;
      sc1d_Bd_fine[slot] = sum.bd * inv;
      sc1d_sa_inf_fine[slot] = sum.sa * inv;
      sc1d_lsf_fine[slot] = sum.lsf * inv;
      sc1d_D_fine[slot] = sum.diff * inv;
      sc1d_Jsd_fine[slot] = sum.jsd * inv;
      sc1d_PN_fine[slot] = sum.pn * inv;
      sc1d_alpha_fine[slot] = sum.alpha * inv;
      sc1d_sigma_fine[slot] = sum.sigma * inv;
   }
   snapshot_interior_lambda();

   std::vector<double>* fields[9] = {
      &sc1d_Bc_fine, &sc1d_Bd_fine, &sc1d_D_fine, &sc1d_lsf_fine,
      &sc1d_sa_inf_fine, &sc1d_PN_fine, &sc1d_alpha_fine, &sc1d_sigma_fine,
      &sc1d_Jsd_fine
   };
   for(int stack = 0; stack < nstack; ++stack){
      const std::size_t base = static_cast<std::size_t>(stack) * static_cast<std::size_t>(nf);
      for(int f = 0; f < 9; ++f){
         if(fields[f]->size() < base + static_cast<std::size_t>(nf)) continue;
         double* field = fields[f]->data() + base;
         for(int i = 0; i < nf; ++i){
            const std::size_t slot = base + static_cast<std::size_t>(i);
            if(sc1d_cell_mask[slot] & sc1d_mask_occupied) continue;
            const int L = sc1d_anchor_left[slot];
            const int R = sc1d_anchor_right[slot];
            const double t = sc1d_anchor_t[slot];
            field[i] = (1.0 - t) * field[L] + t * field[R];
         }
      }
   }

   reapply_smooth_bulk();

   std::copy(sc1d_D_fine.begin(), sc1d_D_fine.begin() + static_cast<std::ptrdiff_t>(nslot), sc1d_Dmat_fine.begin());
   sc1d_sigma_bulk_fine = sc1d_sigma_fine;
   if(sc1d_Rint_base_face.size() == sc1d_Rint_face.size()) sc1d_Rint_face = sc1d_Rint_base_face;

   for(std::size_t slot = 0; slot < nslot; ++slot){
      const int mat = sc1d_mat_fine[slot];
      if(mat < 0 || mat >= static_cast<int>(mp.size())) continue;
      if(mp[static_cast<std::size_t>(mat)].conductivity == 0.0){
         stamp_zero_sigma(slot, mat);
      } else if(mp[static_cast<std::size_t>(mat)].diffusion == 0.0){
         sc1d_D_fine[slot] = 0.0;
         sc1d_Dmat_fine[slot] = 0.0;
      }
   }

   for(int stack = 0; stack < nstack; ++stack){
      const std::size_t base = static_cast<std::size_t>(stack) * static_cast<std::size_t>(nf);
      double* Bc = sc1d_Bc_fine.data() + base;
      double* Bd = sc1d_Bd_fine.data() + base;
      double* D = sc1d_D_fine.data() + base;
      double* Dmat = sc1d_Dmat_fine.data() + base;
      double* lsf = sc1d_lsf_fine.data() + base;
      double* lphi = sc1d_lphi_fine.data() + base;
      double* Jsd = sc1d_Jsd_fine.data() + base;
      double* sigma_b = sc1d_sigma_bulk_fine.data() + base;
      double* sigma_f = sc1d_sigma_fine.data() + base;
      std::vector<double> Bc0(Bc, Bc + nf);
      std::vector<double> lsf0(lsf, lsf + nf);
      for(int i = 0; i < nf; ++i){
         const std::size_t slot = base + static_cast<std::size_t>(i);
         const double acc_up = sc1d_sheet_acc_up_fine[slot];
         const double acc_dn = sc1d_sheet_acc_dn_fine[slot];
         const double acc_re = sc1d_sheet_acc_re_fine[slot];
         const double acc_im = sc1d_sheet_acc_im_fine[slot];
         const double acc_d2 = sc1d_sheet_acc_d2_fine[slot];
         if(!density_set(acc_up, acc_dn, acc_re, acc_im, acc_d2)) continue;
         if(!(sigma_b[i] > 0.0)) continue;
         apply_sheet_density(i, acc_up, acc_dn, acc_re, acc_im, acc_d2,
                             Bc0[static_cast<std::size_t>(i)], Bd[i], Dmat[i],
                             lsf0[static_cast<std::size_t>(i)], lphi[i], Jsd[i],
                             sigma_b[i], Bc, D, lsf, nullptr, Jsd, sigma_f);
      }
   }

   apply_fixed_lambda_delta();

   require_finite("beta_c", sc1d_Bc_fine);
   require_finite("beta_d", sc1d_Bd_fine);
   require_finite("D", sc1d_D_fine);
   require_finite("De", sc1d_Dmat_fine);
   require_finite("lambda_sdl", sc1d_lsf_fine);
   require_finite("Jsd", sc1d_Jsd_fine);
   require_finite("sa_inf", sc1d_sa_inf_fine);
   require_finite("P_N", sc1d_PN_fine);
   require_finite("alpha", sc1d_alpha_fine);
   require_finite("sigma", sc1d_sigma_fine);
   require_finite("sigma_bulk", sc1d_sigma_bulk_fine);
   require_finite("Rint", sc1d_Rint_face);
}

namespace {

constexpr std::size_t kNoMbase = static_cast<std::size_t>(-1);
std::size_t cached_nslot = 0;
std::size_t cached_natoms = 0;
std::vector<std::size_t> atom_mbase;
std::vector<double> atom_mu;
std::vector<std::size_t> anchor_slot_left;
std::vector<std::size_t> anchor_slot_right;
std::vector<double> macc_scratch;
std::vector<double> dmdt_anchor;
std::vector<double> vx_scratch;
std::vector<double> vy_scratch;
std::vector<double> vz_scratch;
std::vector<double> thresh_scratch;
std::vector<char> have_floor;

} // namespace

static double norm3_fine(const double x, const double y, const double z){
   return std::sqrt(x * x + y * y + z * z);
}

// Ta and MgO have no moment. Occupied cells report that through msat.
// Empty cells are classified by the material height, because only sigma = 0
// slabs receive a material id.
static bool fine_cell_nonmagnetic(const std::size_t slot){
   const int nf = sc1d_nf;
   if(nf <= 0 || slot >= sc1d_cell_mask.size()) return false;
   const bool occupied = (sc1d_cell_mask[slot] & sc1d_mask_occupied) != 0;
   if(occupied && !(sc1d_msat_fine[slot] > 0.0)) return true;
   const int nmat = static_cast<int>(mp::mu_s_array.size());
   const int mat = (slot < sc1d_mat_fine.size()) ? sc1d_mat_fine[slot] : -1;
   if(mat >= 0 && mat < nmat) return !(mp::mu_s_array[static_cast<std::size_t>(mat)] > 0.0);
   const int nsub = std::max(1, sc1d_nsub);
   const double dzA = micro_cell_thickness / static_cast<double>(nsub);
   const int fine = static_cast<int>(slot % static_cast<std::size_t>(nf));
   const double z = (static_cast<double>(fine) + 0.5) * dzA;
   const double film_z = cs::system_dimensions[2];
   const int nmp = static_cast<int>(mp.size());
   for(int m = 0; m < nmat && m < nmp; ++m){
      if(mp::mu_s_array[static_cast<std::size_t>(m)] > 0.0) continue;
      const double z_lo = create::get_material_height_min(m) * film_z;
      const double z_hi = create::get_material_height_max(m) * film_z;
      if(z_hi > z_lo && z >= z_lo && z < z_hi) return true;
   }
   return false;
}

void cache_prolong_maps(){
   const int nf = sc1d_nf;
   const int nstack = std::max(0, num_stacks_y);
   cached_natoms = static_cast<std::size_t>(std::max(0, num_local_atoms));
   if(nf <= 0 || nstack <= 0){
      cached_nslot = 0;
      atom_mbase.clear();
      atom_mu.clear();
      anchor_slot_left.clear();
      anchor_slot_right.clear();
      macc_scratch.clear();
      dmdt_anchor.clear();
      vx_scratch.clear();
      vy_scratch.clear();
      vz_scratch.clear();
      thresh_scratch.clear();
      have_floor.clear();
      return;
   }

   const std::size_t nf_sz = static_cast<std::size_t>(nf);
   const std::size_t nslot = static_cast<std::size_t>(nstack) * nf_sz;
   cached_nslot = nslot;
   atom_mbase.assign(cached_natoms, kNoMbase);
   atom_mu.assign(cached_natoms, 0.0);

   const int* stack_of = sc1d_atom_ls.data();
   const int* fine_of = sc1d_atom_fine.data();
   const int* mat_of = atoms::type_array.data();
   const double* mus = mp::mu_s_array.data();
   const int nmat = static_cast<int>(mp::mu_s_array.size());
   for(std::size_t atom = 0; atom < cached_natoms; ++atom){
      const int stack = stack_of[atom];
      const int fine = fine_of[atom];
      if(stack < 0 || stack >= nstack || fine < 0 || fine >= nf) continue;
      const int mat = mat_of[atom];
      if(mat >= 0 && mat < nmat) atom_mu[atom] = mus[mat];
      atom_mbase[atom] = (static_cast<std::size_t>(stack) * nf_sz + static_cast<std::size_t>(fine)) * 3u;
   }

   anchor_slot_left.resize(nslot);
   anchor_slot_right.resize(nslot);
   for(std::size_t slot = 0; slot < nslot; ++slot){
      const std::size_t base = (slot / nf_sz) * nf_sz;
      anchor_slot_left[slot] = base + static_cast<std::size_t>(sc1d_anchor_left[slot]);
      anchor_slot_right[slot] = base + static_cast<std::size_t>(sc1d_anchor_right[slot]);
   }

   macc_scratch.assign(nslot * 3u, 0.0);
   dmdt_anchor.assign(nslot, 0.0);
   vx_scratch.assign(nslot, 0.0);
   vy_scratch.assign(nslot, 0.0);
   vz_scratch.assign(nslot, 0.0);
   thresh_scratch.assign(nslot, 0.0);
   have_floor.assign(nslot, 0);
}

void prolong_fine_magnetisation(const double dt){
   const std::size_t nslot = cached_nslot;
   if(nslot == 0) return;

   double* macc = macc_scratch.data();
   std::fill(macc_scratch.begin(), macc_scratch.end(), 0.0);

   const double* sx = atoms::x_spin_array.data();
   const double* sy = atoms::y_spin_array.data();
   const double* sz = atoms::z_spin_array.data();
   const std::size_t* base_of = atom_mbase.data();
   const double* mu_of = atom_mu.data();
   const std::size_t natoms = cached_natoms;
   for(std::size_t atom = 0; atom < natoms; ++atom){
      const std::size_t b = base_of[atom];
      if(b == kNoMbase) continue;
      const double mu = mu_of[atom];
      macc[b] += sx[atom] * mu;
      macc[b + 1u] += sy[atom] * mu;
      macc[b + 2u] += sz[atom] * mu;
   }
#ifdef MPICF
   if(vmpi::num_processors > 1 && !macc_scratch.empty()){
      MPI_Allreduce(MPI_IN_PLACE, macc, static_cast<int>(macc_scratch.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
#endif

   for(std::size_t slot = 0; slot < nslot; ++slot){
      if(!(sc1d_cell_mask[slot] & sc1d_mask_occupied)) continue;
      if(!(sc1d_msat_fine[slot] > 0.0)){
         sc1d_mhat_ref_fine[slot * 3u + 0u] = 0.0;
         sc1d_mhat_ref_fine[slot * 3u + 1u] = 0.0;
         sc1d_mhat_ref_fine[slot * 3u + 2u] = 0.0;
         vx_scratch[slot] = 0.0;
         vy_scratch[slot] = 0.0;
         vz_scratch[slot] = 0.0;
         dmdt_anchor[slot] = 0.0;
         thresh_scratch[slot] = 0.0;
         have_floor[slot] = 0;
         continue;
      }
      const double mx = macc[slot * 3u];
      const double my = macc[slot * 3u + 1u];
      const double mz = macc[slot * 3u + 2u];
      const double mm = norm3_fine(mx, my, mz);
      // m0 is the cell saturation. The gate uses the reduced length |m|/m0,
      // which runs from 0 to 1. The random-walk scale 1/sqrt(2N) is capped
      // at 1 so a fully aligned few-atom cell still supplies a direction.
      const double m0 = sc1d_msat_fine[slot];
      const double reduced = (m0 > 0.0) ? mm / m0 : 0.0;
      sc1d_ml_fine[slot] = reduced;
      const bool first = !(sc1d_m0_fine[slot] > 0.0);
      if(first && mm > 0.0){
         sc1d_m0_fine[slot] = (m0 > 0.0) ? m0 : mm;
         sc1d_mhat_ref_fine[slot * 3u + 0u] = mx / mm;
         sc1d_mhat_ref_fine[slot * 3u + 1u] = my / mm;
         sc1d_mhat_ref_fine[slot * 3u + 2u] = mz / mm;
      }
      const double N = sc1d_natom_fine[slot];
      bool ok = N >= 1.0 && m0 > 0.0;
      double floor_red = 0.0;
      if(ok){
         floor_red = 1.0 / std::sqrt(2.0 * N);
         if(floor_red > 1.0) floor_red = 1.0;
      }
      const double floor_m = floor_red * m0;
      const bool above = ok && reduced >= floor_red && mm > kEpsFine;
      if(!first && above){
         sc1d_mhat_ref_fine[slot * 3u + 0u] = mx / mm;
         sc1d_mhat_ref_fine[slot * 3u + 1u] = my / mm;
         sc1d_mhat_ref_fine[slot * 3u + 2u] = mz / mm;
      }
      if(sc1d_step_counter <= 1UL){
         dmdt_anchor[slot] = 0.0;
         sc1d_m_prev_fine[slot] = mm;
      } else {
         const double sat = sc1d_msat_fine[slot];
         dmdt_anchor[slot] = (sat > 0.0) ? ((mm - sc1d_m_prev_fine[slot]) / (dt * sat)) : 0.0;
         sc1d_m_prev_fine[slot] = mm;
      }
      if((above || first) && mm > kEpsFine){
         vx_scratch[slot] = mx;
         vy_scratch[slot] = my;
         vz_scratch[slot] = mz;
      } else {
         const double scale = mm;
         vx_scratch[slot] = sc1d_mhat_ref_fine[slot * 3u + 0u] * scale;
         vy_scratch[slot] = sc1d_mhat_ref_fine[slot * 3u + 1u] * scale;
         vz_scratch[slot] = sc1d_mhat_ref_fine[slot * 3u + 2u] * scale;
      }
      thresh_scratch[slot] = floor_m;
      have_floor[slot] = ok ? 1 : 0;
   }

   const double* vx = vx_scratch.data();
   const double* vy = vy_scratch.data();
   const double* vz = vz_scratch.data();
   const double* thresh = thresh_scratch.data();
   const char* floor_ok = have_floor.data();
   const std::size_t* sL_of = anchor_slot_left.data();
   const std::size_t* sR_of = anchor_slot_right.data();
   for(std::size_t slot = 0; slot < nslot; ++slot){
      const std::size_t sL = sL_of[slot];
      const std::size_t sR = sR_of[slot];
      const double t = sc1d_anchor_t[slot];
      if(!(sc1d_cell_mask[slot] & sc1d_mask_occupied)){
         sc1d_ml_fine[slot] = (1.0 - t) * sc1d_ml_fine[sL] + t * sc1d_ml_fine[sR];
      }
      if(fine_cell_nonmagnetic(slot)){
         sc1d_mhat_fine[slot * 3u + 0u] = 0.0;
         sc1d_mhat_fine[slot * 3u + 1u] = 0.0;
         sc1d_mhat_fine[slot * 3u + 2u] = 0.0;
         sc1d_dmdt_fine[slot] = 0.0;
         continue;
      }
      const double px = (1.0 - t) * vx[sL] + t * vx[sR];
      const double py = (1.0 - t) * vy[sL] + t * vy[sR];
      const double pz = (1.0 - t) * vz[sL] + t * vz[sR];
      const double hm = norm3_fine(px, py, pz);
      const double hx = (1.0 - t) * sc1d_mhat_ref_fine[sL * 3u + 0u] + t * sc1d_mhat_ref_fine[sR * 3u + 0u];
      const double hy = (1.0 - t) * sc1d_mhat_ref_fine[sL * 3u + 1u] + t * sc1d_mhat_ref_fine[sR * 3u + 1u];
      const double hz = (1.0 - t) * sc1d_mhat_ref_fine[sL * 3u + 2u] + t * sc1d_mhat_ref_fine[sR * 3u + 2u];
      const double f0 = thresh[sL];
      const double f1 = thresh[sR];
      const bool ok0 = floor_ok[sL] != 0;
      const bool ok1 = floor_ok[sR] != 0;
      bool above = false;
      if(ok0 && ok1) above = hm >= ((1.0 - t) * f0 + t * f1);
      else if(ok0) above = (t < 1.0) && (hm >= f0);
      else if(ok1) above = (t > 0.0) && (hm >= f1);

      double ux = sc1d_mhat_ref_fine[sL * 3u + 0u];
      double uy = sc1d_mhat_ref_fine[sL * 3u + 1u];
      double uz = sc1d_mhat_ref_fine[sL * 3u + 2u];
      if(above && hm > kEpsFine){
         ux = px / hm;
         uy = py / hm;
         uz = pz / hm;
      } else {
         const double hr = norm3_fine(hx, hy, hz);
         if(hr > kEpsFine){
            ux = hx / hr;
            uy = hy / hr;
            uz = hz / hr;
         }
      }
      sc1d_mhat_fine[slot * 3u + 0u] = ux;
      sc1d_mhat_fine[slot * 3u + 1u] = uy;
      sc1d_mhat_fine[slot * 3u + 2u] = uz;
      sc1d_dmdt_fine[slot] = (1.0 - t) * dmdt_anchor[sL] + t * dmdt_anchor[sR];
   }
}

} // namespace internal
} // namespace st
