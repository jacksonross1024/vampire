#include "moire_coarse.hpp"

#include <algorithm>
#include <cmath>

namespace {

void levi_dm(const double J[3][3], double d[3]) {
   d[0] = 0.5 * (J[1][2] - J[2][1]);
   d[1] = 0.5 * (J[2][0] - J[0][2]);
   d[2] = 0.5 * (J[0][1] - J[1][0]);
}

double isotropic_exchange(const double J[3][3]) {
   return (J[0][0] + J[1][1] + J[2][2]) / 3.0;
}

void add_vec3(std::vector<double>& arr, int base, const double v[3]) {
   for (int a = 0; a < 3; ++a)
      arr[static_cast<size_t>(base + a)] += v[a];
}

void add_outer(std::vector<double>& arr, int base, const double r[3],
               const double d[3]) {
   for (int i = 0; i < 3; ++i) {
      for (int a = 0; a < 3; ++a)
         arr[static_cast<size_t>(base + 3 * i + a)] += r[i] * d[a];
   }
}

}  // namespace

MoireAverages coarseGrainMoireBonds(const std::vector<spin_interaction>& bonds,
                                    int nx_moire, int ny_moire, double moire_Lx,
                                    double moire_Ly, double /*cell_volume*/,
                                    double nn_tolerance,
                                    double interlayer_tolerance) {
   MoireAverages out;
   const int ncell = nx_moire * ny_moire;
   out.J_cell.assign(static_cast<size_t>(ncell), 0.0);
   out.J_nn.assign(static_cast<size_t>(ncell), 0.0);
   out.J_homo_mean.assign(static_cast<size_t>(ncell), 0.0);
   out.J_homo_Jr2.assign(static_cast<size_t>(ncell), 0.0);
   out.D_cell_vec.assign(static_cast<size_t>(ncell) * 3u, 0.0);
   out.M_homo.assign(static_cast<size_t>(ncell) * 9u, 0.0);
   out.M_inhomo.assign(static_cast<size_t>(ncell) * 9u, 0.0);
   out.counts_cell.assign(static_cast<size_t>(ncell), 0);
   out.counts_nn.assign(static_cast<size_t>(ncell), 0);
   out.counts_intra.assign(static_cast<size_t>(ncell), 0);

   auto wrap_delta = [](double d, double L) {
      while (d >= 0.5 * L)
         d -= L;
      while (d < -0.5 * L)
         d += L;
      return d;
   };

   auto cell_index_from_position = [&](double x, double y) -> int {
      double fx = x / moire_Lx;
      double fy = y / moire_Ly;
      fx -= std::floor(fx);
      fy -= std::floor(fy);
      int ix = std::min(nx_moire - 1, static_cast<int>(fx * nx_moire));
      int iy = std::min(ny_moire - 1, static_cast<int>(fy * ny_moire));
      return iy * nx_moire + ix;
   };

   for (const auto& b : bonds) {
      // Microcell averages: only interfacial Cr sheets (h_id == 1 on both sites).
      if (b.h_i != 1 || b.h_j != 1)
         continue;

      int c = cell_index_from_position(b.xi, b.yi);

      double r[3];
      r[0] = wrap_delta(b.xj - b.xi, moire_Lx);
      r[1] = wrap_delta(b.yj - b.yi, moire_Ly);
      r[2] = b.zj - b.zi;
      double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];

      double Jiso = isotropic_exchange(b.J);
      double d[3];
      levi_dm(b.J, d);

      bool same_layer = (b.layer_i == b.layer_j);
      bool same_sublattice = (b.si == b.sj);
      bool opposite_sublattice = !same_sublattice;

      bool intracell_lateral =
          (std::abs(r[0]) < nn_tolerance && std::abs(r[1]) < nn_tolerance);
      bool interlayer_bond = (std::abs(r[2]) > interlayer_tolerance);

      //intracell aex
      if (interlayer_bond && opposite_sublattice && intracell_lateral) {
         out.J_cell[static_cast<size_t>(c)] += Jiso;
         add_vec3(out.D_cell_vec, 3 * c, d);
         out.counts_cell[static_cast<size_t>(c)] += 1;
         continue;
      }

      //inhomogeneous exchange intercell
      if (interlayer_bond && opposite_sublattice && !intracell_lateral) {
         out.J_nn[static_cast<size_t>(c)] += Jiso;
         add_outer(out.M_inhomo, 9 * c, r, d);
         out.counts_nn[static_cast<size_t>(c)] += 1;
         continue;
      }

      if (same_layer && same_sublattice) {
         out.J_homo_mean[static_cast<size_t>(c)] += Jiso;
         out.J_homo_Jr2[static_cast<size_t>(c)] += Jiso * r2;
         add_outer(out.M_homo, 9 * c, r, d);
         out.counts_intra[static_cast<size_t>(c)] += 1;
      }
   }

   for (int c = 0; c < ncell; ++c) {
      size_t ci = static_cast<size_t>(c);
      if (out.counts_cell[ci] > 0) {
         out.J_cell[ci] /= static_cast<double>(out.counts_cell[ci]);
         for (int a = 0; a < 3; ++a)
            out.D_cell_vec[ci * 3u + static_cast<size_t>(a)] /=
                static_cast<double>(out.counts_cell[ci]);
      }
      if (out.counts_nn[ci] > 0) {
         out.J_nn[ci] /= static_cast<double>(out.counts_nn[ci]);
         for (int q = 0; q < 9; ++q)
            out.M_inhomo[ci * 9u + static_cast<size_t>(q)] /=
                static_cast<double>(out.counts_nn[ci]);
      }
      if (out.counts_intra[ci] > 0) {
         const double inv = 1.0 / static_cast<double>(out.counts_intra[ci]);
         out.J_homo_mean[ci] *= inv;
         out.J_homo_Jr2[ci] *= inv;
         for (int q = 0; q < 9; ++q)
            out.M_homo[ci * 9u + static_cast<size_t>(q)] *= inv;
      }
   }

   return out;
}

TwistedDoubleBilayerMoireAverages coarseGrainTwistedDoubleBilayerMoireBonds(
    const std::vector<spin_interaction>& bonds, int nx_moire, int ny_moire,
    double moire_Lx, double moire_Ly, double /*cell_volume*/,
    double nn_tolerance, double interlayer_tolerance) {
   TwistedDoubleBilayerMoireAverages out;
   const int ncell = nx_moire * ny_moire;
   out.J_23_cell.assign(static_cast<size_t>(ncell), 0.0);
   out.J_12_nn.assign(static_cast<size_t>(ncell), 0.0);
   out.J_34_nn.assign(static_cast<size_t>(ncell), 0.0);
   out.J_layer1_homo_mean.assign(static_cast<size_t>(ncell), 0.0);
   out.J_layer1_homo_Jr2.assign(static_cast<size_t>(ncell), 0.0);
   out.J_layer2_homo_mean.assign(static_cast<size_t>(ncell), 0.0);
   out.J_layer2_homo_Jr2.assign(static_cast<size_t>(ncell), 0.0);
   out.J_layer3_homo_mean.assign(static_cast<size_t>(ncell), 0.0);
   out.J_layer3_homo_Jr2.assign(static_cast<size_t>(ncell), 0.0);
   out.J_layer4_homo_mean.assign(static_cast<size_t>(ncell), 0.0);
   out.J_layer4_homo_Jr2.assign(static_cast<size_t>(ncell), 0.0);
   out.counts_23_cell.assign(static_cast<size_t>(ncell), 0);
   out.counts_12_nn.assign(static_cast<size_t>(ncell), 0);
   out.counts_34_nn.assign(static_cast<size_t>(ncell), 0);
   out.counts_layer1_homo.assign(static_cast<size_t>(ncell), 0);
   out.counts_layer2_homo.assign(static_cast<size_t>(ncell), 0);
   out.counts_layer3_homo.assign(static_cast<size_t>(ncell), 0);
   out.counts_layer4_homo.assign(static_cast<size_t>(ncell), 0);
   out.D_23_cell_vec.assign(static_cast<size_t>(ncell) * 3u, 0.0);
   out.M_inhomo_12.assign(static_cast<size_t>(ncell) * 9u, 0.0);
   out.M_inhomo_34.assign(static_cast<size_t>(ncell) * 9u, 0.0);
   out.M_homo_layer1.assign(static_cast<size_t>(ncell) * 9u, 0.0);
   out.M_homo_layer2.assign(static_cast<size_t>(ncell) * 9u, 0.0);
   out.M_homo_layer3.assign(static_cast<size_t>(ncell) * 9u, 0.0);
   out.M_homo_layer4.assign(static_cast<size_t>(ncell) * 9u, 0.0);

   auto wrap_delta = [](double d, double L) {
      while (d >= 0.5 * L)
         d -= L;
      while (d < -0.5 * L)
         d += L;
      return d;
   };

   auto cell_index_from_position = [&](double x, double y) -> int {
      double fx = x / moire_Lx;
      double fy = y / moire_Ly;
      fx -= std::floor(fx);
      fy -= std::floor(fy);
      int ix = std::min(nx_moire - 1, static_cast<int>(fx * nx_moire));
      int iy = std::min(ny_moire - 1, static_cast<int>(fy * ny_moire));
      return iy * nx_moire + ix;
   };

   auto accumulate_layer_homo = [&](int layer, int cell, double Jiso,
                                    double r2) {
      const size_t ci = static_cast<size_t>(cell);
      switch (layer) {
         case 1:
            out.J_layer1_homo_mean[ci] += Jiso;
            out.J_layer1_homo_Jr2[ci] += Jiso * r2;
            out.counts_layer1_homo[ci] += 1;
            return;
         case 2:
            out.J_layer2_homo_mean[ci] += Jiso;
            out.J_layer2_homo_Jr2[ci] += Jiso * r2;
            out.counts_layer2_homo[ci] += 1;
            return;
         case 3:
            out.J_layer3_homo_mean[ci] += Jiso;
            out.J_layer3_homo_Jr2[ci] += Jiso * r2;
            out.counts_layer3_homo[ci] += 1;
            return;
         case 4:
            out.J_layer4_homo_mean[ci] += Jiso;
            out.J_layer4_homo_Jr2[ci] += Jiso * r2;
            out.counts_layer4_homo[ci] += 1;
            return;
         default:
            return;
      }
   };

   auto normalize_scalar_family = [&](std::vector<double>& values,
                                      const std::vector<int>& counts) {
      for (int c = 0; c < ncell; ++c) {
         const size_t ci = static_cast<size_t>(c);
         if (counts[ci] > 0)
            values[ci] /= static_cast<double>(counts[ci]);
      }
   };

   auto normalize_homo_family = [&](std::vector<double>& means,
                                    std::vector<double>& jr2,
                                    const std::vector<int>& counts) {
      for (int c = 0; c < ncell; ++c) {
         const size_t ci = static_cast<size_t>(c);
         if (counts[ci] <= 0)
            continue;
         const double inv = 1.0 / static_cast<double>(counts[ci]);
         means[ci] *= inv;
         jr2[ci] *= inv;
      }
   };

   auto normalize_homo_dmi_layer = [&](std::vector<double>& tens,
                                       const std::vector<int>& counts) {
      for (int c = 0; c < ncell; ++c) {
         const size_t ci = static_cast<size_t>(c);
         if (counts[ci] <= 0)
            continue;
         const double inv = 1.0 / static_cast<double>(counts[ci]);
         for (int q = 0; q < 9; ++q)
            tens[ci * 9u + static_cast<size_t>(q)] *= inv;
      }
   };

   auto add_homo_dmi_for_layer = [&](int layer, int cell, const double r[3],
                                     const double d[3]) {
      const size_t base = 9u * static_cast<size_t>(cell);
      std::vector<double>* tgt = nullptr;
      switch (layer) {
         case 1:
            tgt = &out.M_homo_layer1;
            break;
         case 2:
            tgt = &out.M_homo_layer2;
            break;
         case 3:
            tgt = &out.M_homo_layer3;
            break;
         case 4:
            tgt = &out.M_homo_layer4;
            break;
         default:
            return;
      }
      add_outer(*tgt, static_cast<int>(base), r, d);
   };

   for (const auto& b : bonds) {
      // h_id: 0 = outer Cr sheets (layers 1, 4 → mumax cells 0, 2); 1 = interface
      // sheets (layers 2, 3 → middle cell). Include both so L1/L4 intralayer and
      // 1–2 / 3–4 couplings are not dropped.
      if (b.layer_i < 1 || b.layer_i > 4 || b.layer_j < 1 || b.layer_j > 4)
         continue;

      const int c = cell_index_from_position(b.xi, b.yi);
      double r[3];
      r[0] = wrap_delta(b.xj - b.xi, moire_Lx);
      r[1] = wrap_delta(b.yj - b.yi, moire_Ly);
      r[2] = b.zj - b.zi;
      const double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
      const double Jiso = isotropic_exchange(b.J);
      double d_dm[3];
      levi_dm(b.J, d_dm);

      if (b.layer_i == b.layer_j) {
         accumulate_layer_homo(b.layer_i, c, Jiso, r2);
         add_homo_dmi_for_layer(b.layer_i, c, r, d_dm);
         continue;
      }

      const bool interlayer_bond = (std::abs(r[2]) > interlayer_tolerance);
      if (!interlayer_bond)
         continue;

      const bool intracell_lateral =
          (std::abs(r[0]) < nn_tolerance && std::abs(r[1]) < nn_tolerance);
      const int lo = std::min(b.layer_i, b.layer_j);
      const int hi = std::max(b.layer_i, b.layer_j);
      const size_t ci = static_cast<size_t>(c);

      if (lo == 2 && hi == 3 && intracell_lateral) {
         out.J_23_cell[ci] += Jiso;
         add_vec3(out.D_23_cell_vec, 3 * c, d_dm);
         out.counts_23_cell[ci] += 1;
         continue;
      }

      // The outer bilayer interfaces map onto neighbouring z cells in the
      // 3-cell micromagnetic stack, so keep them as separate intercell families.
      if (lo == 1 && hi == 2) {
         out.J_12_nn[ci] += Jiso;
         add_outer(out.M_inhomo_12, 9 * c, r, d_dm);
         out.counts_12_nn[ci] += 1;
         continue;
      }
      if (lo == 3 && hi == 4) {
         out.J_34_nn[ci] += Jiso;
         add_outer(out.M_inhomo_34, 9 * c, r, d_dm);
         out.counts_34_nn[ci] += 1;
      }
   }

   normalize_scalar_family(out.J_23_cell, out.counts_23_cell);
   normalize_scalar_family(out.J_12_nn, out.counts_12_nn);
   normalize_scalar_family(out.J_34_nn, out.counts_34_nn);
   normalize_homo_family(out.J_layer1_homo_mean, out.J_layer1_homo_Jr2,
                         out.counts_layer1_homo);
   normalize_homo_family(out.J_layer2_homo_mean, out.J_layer2_homo_Jr2,
                         out.counts_layer2_homo);
   normalize_homo_family(out.J_layer3_homo_mean, out.J_layer3_homo_Jr2,
                         out.counts_layer3_homo);
   normalize_homo_family(out.J_layer4_homo_mean, out.J_layer4_homo_Jr2,
                         out.counts_layer4_homo);

   for (int c = 0; c < ncell; ++c) {
      const size_t ci = static_cast<size_t>(c);
      if (out.counts_23_cell[ci] > 0) {
         const double inv = 1.0 / static_cast<double>(out.counts_23_cell[ci]);
         for (int a = 0; a < 3; ++a)
            out.D_23_cell_vec[ci * 3u + static_cast<size_t>(a)] *= inv;
      }
      if (out.counts_12_nn[ci] > 0) {
         for (int q = 0; q < 9; ++q)
            out.M_inhomo_12[ci * 9u + static_cast<size_t>(q)] /=
                static_cast<double>(out.counts_12_nn[ci]);
      }
      if (out.counts_34_nn[ci] > 0) {
         for (int q = 0; q < 9; ++q)
            out.M_inhomo_34[ci * 9u + static_cast<size_t>(q)] /=
                static_cast<double>(out.counts_34_nn[ci]);
      }
   }
   normalize_homo_dmi_layer(out.M_homo_layer1, out.counts_layer1_homo);
   normalize_homo_dmi_layer(out.M_homo_layer2, out.counts_layer2_homo);
   normalize_homo_dmi_layer(out.M_homo_layer3, out.counts_layer3_homo);
   normalize_homo_dmi_layer(out.M_homo_layer4, out.counts_layer4_homo);

   return out;
}
