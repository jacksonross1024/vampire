#ifndef MOIRE_WIRE_HPP
#define MOIRE_WIRE_HPP

#include <vector>
#include "positions.hpp"

/// Global list filled during interaction generation (see exchange.cpp).
extern std::vector<spin_interaction> g_moire_spin_interactions;

void moire_spin_interactions_clear();
void moire_spin_interactions_merge(std::vector<spin_interaction>& thread_chunk);

/// Run coarse-grain helpers on ``g_moire_spin_interactions`` and write binaries
/// in *cwd* according to ``moire_coarse_write_twisted_*`` (default: MOTD only;
/// optional ``--moire-export`` in main: ``both``, ``bilayer_only``, ``double_only``).
/// *moire_phys* is ``{ x0, y0, Lx, Ly }`` in the same length units as atom
/// coordinates (Å), after any in-place update from build_spin_moire_cell.
/// Grid: ``nx_m = floor(Lx/(2*a0x))``, ``ny_m = floor(Ly/(2*a1y))`` (independent of
/// microcell_Nx/Ny used elsewhere). nn_tol from ``moire_coarse_nn_tol_scale``.
void moire_spin_interactions_finalize_and_write(const char* cwd,
                                                 const double moire_phys[4]);

#endif
