#ifndef INITALISE_HPP
#define INITALISE_HPP

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>

#include <unistd.h>
extern double angle;
extern double twist_loction;

extern double system_size_x;
extern double system_size_y;
extern double system_size_z;
extern int number_of_unit_cells_z;

// Microcell resolution: each microcell spans microcell_scale_x*a0x by microcell_scale_y*a1y (use 1, 2, or 3).
// unit_x_lr = floor(x / (microcell_scale_x*a0x)), unit_y_lr = floor(y / (microcell_scale_y*a1y)).
extern int microcell_scale_x;
extern int microcell_scale_y;
extern int microcell_Nx;  // number_of_unit_cells_x / microcell_scale_x
extern int microcell_Ny;  // number_of_unit_cells_y / microcell_scale_y

/// Dimensionless factor for moiré micromagnetic coarse grain: nn_tol (Å) =
/// moire_coarse_nn_tol_scale * min(cell_wx, cell_wy). If this exceeds ~9.9 Å
/// relative to in-plane bond extent, interlayer terms fall in the intracell bin.
/// Optional CLI argv[11]; default matches former hard-coded 0.51.
extern double moire_coarse_nn_tol_scale;

/// Which moiré coarse-grain binaries ``moire_spin_interactions_finalize_and_write`` emits.
/// Defaults: twisted-double-bilayer only (MOTD). Pass optional argv[12] ``both`` to also
/// write legacy twisted-bilayer ``moire_coarse_v2.bin``; ``bilayer_only`` for legacy only.
extern bool moire_coarse_write_twisted_bilayer;
extern bool moire_coarse_write_twisted_double_bilayer;

extern double dmi12;
extern double dmi23;
extern double dmi34;
extern double dmi_decay;
extern bool DMI;
extern double exchange12;
extern double exchange23;
extern double exchange34;

extern double separation;

void initialise_variables();
void resize_arrays(std::vector < std::vector < double > > &A, int sizex, int sizey);

#endif
