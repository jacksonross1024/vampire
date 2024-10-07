#ifndef INITALISE_HPP
#define INITALISE_HPP

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>

extern double angle;
extern double twist_loction;

extern double system_size_x;
extern double system_size_y;
extern double system_size_z;
extern int number_of_unit_cells_z;

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
