#ifndef MOIRE_COARSE_HPP
#define MOIRE_COARSE_HPP

#include <vector>
#include "positions.hpp"

/// Coarse-grained Moiré-cell quantities in atomistic units (meV, meV·Å, etc.),
/// produced from a list of bond tensors. See planning doc for .npz mapping.
struct MoireAverages {
   std::vector<double> J_cell;
   std::vector<double> J_nn;
   /// Mean isotropic J (meV), same-layer same-sublattice (no r² weight).
   std::vector<double> J_homo_mean;
   /// Mean J·|r|² (meV·Å²), same bond class — for stiffness / cell_volume conversion.
   std::vector<double> J_homo_Jr2;
   std::vector<double> D_cell_vec;
   std::vector<double> M_homo;
   std::vector<double> M_inhomo;
   std::vector<int> counts_cell;
   std::vector<int> counts_nn;
   std::vector<int> counts_intra;
};

/// Twisted-double-bilayer coarse-grained exchange families in atomistic units.
/// The four physical layers are retained explicitly while the twisted 2-3
/// interface is reduced onto the middle AFM cell. Uses all bonds with
/// layer_i/j in 1–4 (from atom S); h_id 0 and 1 are both kept (outer layers vs
/// interface), unlike coarseGrainMoireBonds which restricts to h_id == 1.
struct TwistedDoubleBilayerMoireAverages {
   std::vector<double> J_23_cell;
   std::vector<double> J_12_nn;
   std::vector<double> J_34_nn;
   std::vector<double> J_layer1_homo_mean;
   std::vector<double> J_layer1_homo_Jr2;
   std::vector<double> J_layer2_homo_mean;
   std::vector<double> J_layer2_homo_Jr2;
   std::vector<double> J_layer3_homo_mean;
   std::vector<double> J_layer3_homo_Jr2;
   std::vector<double> J_layer4_homo_mean;
   std::vector<double> J_layer4_homo_Jr2;
   std::vector<int> counts_23_cell;
   std::vector<int> counts_12_nn;
   std::vector<int> counts_34_nn;
   std::vector<int> counts_layer1_homo;
   std::vector<int> counts_layer2_homo;
   std::vector<int> counts_layer3_homo;
   std::vector<int> counts_layer4_homo;

   /// DMI (antisymmetric part of bond tensor): lab-frame bond vector ``r`` with
   /// axial vector ``d = levi(J)`` (see ``coarseGrainMoireBonds``). Same family
   /// routing as exchange: L2–L3 intracell lateral → ``D_23_cell_vec``; L1–L2 /
   /// L3–L4 intercell → ``M_inhomo_*``; same-layer → per-layer ``M_homo_layer*``.
   std::vector<double> D_23_cell_vec;
   std::vector<double> M_inhomo_12;
   std::vector<double> M_inhomo_34;
   std::vector<double> M_homo_layer1;
   std::vector<double> M_homo_layer2;
   std::vector<double> M_homo_layer3;
   std::vector<double> M_homo_layer4;
};

/// Bin bonds onto an nx_moire × ny_moire grid with periodic xy wrapping.
/// \param cell_volume  Reserved for future normalization; unused in current implementation.
MoireAverages coarseGrainMoireBonds(
   const std::vector<spin_interaction>& bonds,
   int nx_moire,
   int ny_moire,
   double moire_Lx,
   double moire_Ly,
   double cell_volume,
   double nn_tolerance,
   double interlayer_tolerance);

/// Twisted-double-bilayer reduction that keeps four physical-layer stiffness
/// families plus the three interlayer exchange families used by the 3-cell AFM
/// mumaxplus stack.
TwistedDoubleBilayerMoireAverages coarseGrainTwistedDoubleBilayerMoireBonds(
   const std::vector<spin_interaction>& bonds,
   int nx_moire,
   int ny_moire,
   double moire_Lx,
   double moire_Ly,
   double cell_volume,
   double nn_tolerance,
   double interlayer_tolerance);

#endif
