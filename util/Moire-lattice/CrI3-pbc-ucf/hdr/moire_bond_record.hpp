#ifndef MOIRE_BOND_RECORD_HPP
#define MOIRE_BOND_RECORD_HPP

#include <array>
#include "positions.hpp"

/// Build a bond tensor matching the UCF normalised-tensorial layout (exchange.cpp
/// otext lines): xx=ex[0], xy=ex[3], xz=-ex[2], yx=-ex[3], yy=ex[0], yz=ex[1],
/// zx=ex[2], zy=-ex[1], zz=ex[0]. Isotropic branch uses ex[0] on the diagonal.
inline spin_interaction make_spin_interaction_from_exchange(
    const spin& atom_i, const spin& atom_j,
    const std::array<float, 4>& ex,
    bool tensorial_dmi) {
   spin_interaction b;
   b.xi = atom_i.x;
   b.yi = atom_i.y;
   b.zi = atom_i.z;
   b.xj = atom_j.x;
   b.yj = atom_j.y;
   b.zj = atom_j.z;
   b.si = static_cast<int>(atom_i.S);
   b.sj = static_cast<int>(atom_j.S);
   b.layer_i = b.si;
   b.layer_j = b.sj;
   b.h_i = atom_i.h_id;
   b.h_j = atom_j.h_id;
   if (tensorial_dmi) {
      b.J[0][0] = ex[0];
      b.J[0][1] = ex[3];
      b.J[0][2] = -ex[2];
      b.J[1][0] = -ex[3];
      b.J[1][1] = ex[0];
      b.J[1][2] = ex[1];
      b.J[2][0] = ex[2];
      b.J[2][1] = -ex[1];
      b.J[2][2] = ex[0];
   } else {
      for (int r = 0; r < 3; ++r) {
         for (int c = 0; c < 3; ++c) {
            b.J[r][c] = 0.0;
         }
      }
      b.J[0][0] = b.J[1][1] = b.J[2][2] = ex[0];
   }
   return b;
}

#endif
