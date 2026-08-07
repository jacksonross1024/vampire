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
// Binning uses a half-cell grid offset: unit_x_lr = floor((x + ax/2)/ax), unit_y_lr = floor((y + ay/2)/ay).
extern int microcell_scale_x;
extern int microcell_scale_y;
extern int microcell_Nx;  // number_of_unit_cells_x / microcell_scale_x
extern int microcell_Ny;  // number_of_unit_cells_y / microcell_scale_y

extern double dmi12;
extern double dmi23;
extern double dmi34;
extern double dmi_decay;
extern bool DMI;
extern bool BQ;
extern bool bake_only; // --bake-only: detect moiré period, write bake file, skip exchange
extern bool write_config_output; // --config-atoms: write config_energy_atomic/cells.txt
extern int moire_force_candidate; // --candidate N: force ranked AA quartet index (-1 = try all)
extern double exchange12;
extern double exchange23;
extern double exchange34;

extern double separation;

void initialise_variables();
void resize_arrays(std::vector < std::vector < double > > &A, int sizex, int sizey);

#endif
