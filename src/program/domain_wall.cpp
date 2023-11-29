//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and Andrea Meo 2014-2018. All rights reserved.
//
//-----------------------------------------------------------------------------
//

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"
#include "../exchange/internal.hpp"
#include "anisotropy.hpp"
#include "../simulate/internal.hpp"
#include "../program/internal.hpp"
namespace program{

	//------------------------------------------------------------------------------
	// Program to calculate a simple time series
	//------------------------------------------------------------------------------
	void domain_wall(){

		// check calling of routine if error checking is activated
		if(err::check==true) std::cout << "program::domain walls has been called" << std::endl;

		char directory [256];
      	if(getcwd(directory, sizeof(directory)) == NULL){
            std::cerr << "Fatal getcwd error in datalog." << std::endl;
      	}

		if(sim::domain_wall_angle == -1) std::cerr << "no domain angle specified" << std::endl;

		double temp=sim::temperature;
		#ifdef MPICF
		int num_local_atoms =  vmpi::num_core_atoms+vmpi::num_bdry_atoms;
		#else
		int num_local_atoms = atoms::num_atoms;
		#endif

		int num_averages = 50;
		int num_mag_cat = (3 + 3);
		int num_mag_types = 2;
		num_averages = num_averages/2.0;
		sim::domain_wall_centre = 0.1*sim::domain_wall_position; //nm
		sim::domain_wall_velocity = 0.0;
		//if(sim::domain_wall_discretisation == 1) 
		//sim::domain_wall_discretisation[] = 1.0*sim::unit_cell_x;// + 0.0001; 
		
		sim::domain_wall_discretisation[0] *= sim::unit_cell_x;
		sim::domain_wall_discretisation[1] *= sim::unit_cell_y;
		sim::domain_wall_discretisation[2] *= sim::unit_cell_z/6.0;
		
		int num_dw_cells_x = (cs::system_dimensions[0]/sim::domain_wall_discretisation[0]) + 1;
		int num_dw_cells_y = (cs::system_dimensions[1]/sim::domain_wall_discretisation[1]) + 1;
		int num_dw_cells_z = (cs::system_dimensions[2]/sim::domain_wall_discretisation[2]) + 1;
			if(num_dw_cells_z < 0) num_dw_cells_z = 6;
			std::cout << "system dimension x: " << cs::system_dimensions[0] << ", discretisation cells x (A): " << sim::domain_wall_discretisation[0] << ", num cells(per type): " << num_dw_cells_x << ", types: " << num_mag_types << std::endl;
			std::cout << "system dimension y: " << cs::system_dimensions[1] << ", discretisation cells y (A): " << sim::domain_wall_discretisation[1] << ", num cells(per type): " << num_dw_cells_y << ", types: " << num_mag_types << std::endl;
			std::cout << "system dimension z: " << cs::system_dimensions[2] << ", discretisation cells z (A): " << sim::domain_wall_discretisation[2] << ", num cells(per type): " << num_dw_cells_z << ", types: " << num_mag_types << std::endl;

		int num_dw_cells = num_mag_types*num_dw_cells_z*num_dw_cells_y*num_dw_cells_x;
		// int num_categories = 6;
		std::vector <double > atom_to_cell_array(num_local_atoms,0.0);
		std::vector <int > cell_to_lattice_array(3*num_dw_cells, 0);
		std::vector  < double > mag(num_dw_cells*num_mag_cat, 0.0);

		std::vector < int > num_atoms_in_cell(num_dw_cells,0);
		// std::vector<double> topological_charge(num_categories*num_dw_cells,0.0);
		// for (int mat = 0; mat < mp::num_materials; mat ++){
		// 	// if(sim::domain_wall_angle == 0) std::cout <<mat << " 90 degree wall " <<  sim::domain_wall_second_vector_x[mat] << "\t" <<sim::domain_wall_second_vector_y[mat] << "\t" <<sim::domain_wall_second_vector_z[mat] << "\t" << std::endl;
		// 	// if(sim::domain_wall_angle == 1) std::cout <<mat << " 180 degree wall " <<  sim::domain_wall_second_vector_x[mat] << "\t" <<sim::domain_wall_second_vector_y[mat] << "\t" <<sim::domain_wall_second_vector_z[mat] << "\t" << std::endl;
	
		// }
	//	int avg_atom_num = 0;
		//reverses the magentisation of atoms further away than the domain wall distance.
		if (!sim::load_checkpoint_flag){
			if (sim::domain_wall_axis == 0) {
				//90 degree 
				if(sim::domain_wall_angle == 0) {
					for(int atom=0;atom<num_local_atoms;atom++) {
					//		std::cout <<atom << '\t' <<  atoms::x_coord_array[atom] << "\t" << cs::system_dimensions[0]*sim::domain_wall_position -sim::domain_wall_width/2.0 << std::endl;
						// if (atoms::x_coord_array[atom] > cs::system_dimensions[0]*sim::domain_wall_position -sim::domain_wall_width*30.0) {
							int mat = atoms::type_array[atom]-1;
							double mod = 1.0;///sqrt(pos_x*pos_x + pos_y*pos_y + 1.0);

							double theta = std::atan(exp(-1.0*(atoms::x_coord_array[atom] - sim::domain_wall_position)/sim::domain_wall_width)) -M_PI*0.25;
							double pos_x = std::cos(theta)*mod;
							double pos_y = std::sin(theta)*mod; //std::tanh(pos-M_PI/4.0)/(std::cosh(pos-M_PI/4.0)*std::cosh(pos-M_PI/4.0));
					
							atoms::x_spin_array[atom] =  (mat==1?-1.0:1.0)* pos_x;
							atoms::y_spin_array[atom] =  (mat==1?-1.0:1.0)* pos_y;
							atoms::z_spin_array[atom] = 0.0;// (sim::domain_wall_second_vector_z[mat] - atoms::z_spin_array[atom])*theta;
						// int cat = atoms::sublayer_array[atom];
						// if(cat > 5) std::cout << cat << std::endl;
						int x_cell = (atoms::x_coord_array[atom]+0.0001)/sim::domain_wall_discretisation[0];
						int y_cell = (atoms::y_coord_array[atom]+0.0001)/sim::domain_wall_discretisation[1];
						int z_cell = (atoms::z_coord_array[atom]+0.0001)/sim::domain_wall_discretisation[2];
				//				std::cout << atom << '\t' << mat << '\t' << cell << "\t" << atom_to_cell_array.size() << "\t" << num_atoms_in_cell.size() << '\t' << num_dw_cells*mat + cell << std::endl;
						int cell = z_cell*num_dw_cells_x*num_dw_cells_y*num_mag_types + y_cell*num_dw_cells_x*num_mag_types + x_cell*num_mag_types + mat;
						//          1* + 
					//	std::cout << cell << "\t";
						atom_to_cell_array[atom] = cell;
						cell_to_lattice_array[3*cell + 0] = x_cell;
						cell_to_lattice_array[3*cell + 1] = y_cell;
						cell_to_lattice_array[3*cell + 2] = z_cell;

						num_atoms_in_cell[cell]++;
					//	avg_atom_num ++;
						// if( z_cell > 5) {
						// 	std::cout << cell << ", " << z_cell << ", " << atoms::z_coord_array[atom]+0.001 << ", " << sim::domain_wall_discretisation[2] << std::endl;
						// }
						// if (num_atoms_in_cell[num_dw_cells*cat + cell] < 0) std::cout << cell << ", " << atoms::x_coord_array[atom] << ", " << atoms::y_coord_array[atom] << ", " << atoms::z_coord_array[atom] << std::endl;
					}
				}
				
				//180 degree:
				else if(sim::domain_wall_angle == 1) {
					for(int atom=0;atom<num_local_atoms;atom++) {
					//		std::cout <<atom << '\t' <<  atoms::x_coord_array[atom] << "\t" << cs::system_dimensions[0]*sim::domain_wall_position -sim::domain_wall_width/2.0 << std::endl;
						// if (atoms::x_coord_array[atom] > cs::system_dimensions[0]*sim::domain_wall_position -sim::domain_wall_width*3.0) {
						int mat = atoms::type_array[atom];
						double mod = 1.0;///sqrt(pos_x*pos_x + pos_y*pos_y + 1.0);

					
						double theta = -std::atan(std::sinh((atoms::x_coord_array[atom] - sim::domain_wall_position)/sim::domain_wall_width));
						
						double pos_x = std::cos(theta)*mod; 
						double pos_y = std::sin(theta)*mod;
					
						atoms::x_spin_array[atom] =  (mat==2?-1.0:1.0)* pos_x;
						atoms::y_spin_array[atom] =  (mat==2?-1.0:1.0)* pos_y;
						atoms::z_spin_array[atom] +=  (sim::domain_wall_second_vector_z[mat] - atoms::z_spin_array[atom])*theta;
						
						int x_cell = (atoms::x_coord_array[atom]+0.001)/sim::domain_wall_discretisation[0];
						int y_cell = (atoms::y_coord_array[atom]+0.001)/sim::domain_wall_discretisation[1];
						int z_cell = (atoms::z_coord_array[atom]+0.001)/sim::domain_wall_discretisation[2];
				
						int cell =  z_cell*num_dw_cells_x*num_dw_cells_y + y_cell*num_dw_cells_x + x_cell;
						atom_to_cell_array[atom] = cell;
						num_atoms_in_cell[cell] ++;
					}
				}
			}
			if (sim::domain_wall_axis == 1){
				for(int atom=0;atom<num_local_atoms;atom++) {
					//		std::cout <<atom << '\t' <<  atoms::x_coord_array[atom] << "\t" << cs::system_dimensions[0]*sim::domain_wall_position -sim::domain_wall_width/2.0 << std::endl;
					if (atoms::y_coord_array[atom] > sim::domain_wall_position -sim::domain_wall_width*3.0) {
						int mat = atoms::type_array[atom];
						double pos = std::tanh((atoms::y_coord_array[atom] - sim::domain_wall_position)*M_PI/sim::domain_wall_width);
						//volatile double pos_x = std::cos(pos);
						//volatile double pos_y = std::sin(pos); //std::tanh(pos-M_PI/4.0)/(std::cosh(pos-M_PI/4.0)*std::cosh(pos-M_PI/4.0));
					
						double mod = 1.0;///sqrt(pos_x*pos_x + pos_y*pos_y + 1.0);
					
						atoms::x_spin_array[atom] =  (mat==2?-1.0:1.0)* sqrt(1.0-pos*pos)*mod;
						atoms::y_spin_array[atom] =  (mat==2?-1.0:1.0)* pos*mod;
						atoms::z_spin_array[atom] +=  (sim::domain_wall_second_vector_z[mat] - atoms::z_spin_array[atom])*pos;
			
					}
				}
			}
			// if (sim::domain_wall_axis == 2){
			// 	for(int atom=0;atom<num_local_atoms;atom++){
			// 		//			std::cout <<" here" << "\t" << atom << '\t' <<  atoms::z_coord_array[atom] << "\t" << cs::system_dimensions[2] << "\t" << sim::domain_wall_position  << '\t'<< sim::domain_wall_width/2.0 << "\t" << cs::system_dimensions[2]*sim::domain_wall_position -sim::domain_wall_width/2.0 << std::endl;
			// 		if (atoms::z_coord_array[atom] > cs::system_dimensions[2]*sim::domain_wall_position -sim::domain_wall_width/2.0){
			// 			//	 std::cout << "entered" << std::endl;
			// 			int mat = atoms::type_array[atom];
			// 			double pos = (atoms::z_coord_array[atom] - cs::system_dimensions[2]*sim::domain_wall_position + sim::domain_wall_width/2.0)/sim::domain_wall_width;
			// 			if (pos > 1) pos = 1;
			// 			if (pos < 0) pos = 0;
			// 			//				std::cout << atoms::z_coord_array[atom] << '\t' << pos << std::endl;
			// 			double dx = -atoms::x_spin_array[atom] + sim::domain_wall_second_vector_x[mat];
			// 			double dy = -atoms::y_spin_array[atom] + sim::domain_wall_second_vector_y[mat];
			// 			double dz = -atoms::z_spin_array[atom] + sim::domain_wall_second_vector_z[mat];
			// 			//	std::cout << dx << '\t' << dy << "\t" << dz << std::endl;
			// 			double mx = atoms::x_spin_array[atom] + dx*pos;
			// 			double my = atoms::y_spin_array[atom] + dy*pos;
			// 			double mz = atoms::z_spin_array[atom] + dz*pos;
			// 			//	std::cout << dx << '\t' << dy << "\t" << dz << std::endl;
			// 			atoms::x_spin_array[atom] = mx;
			// 			atoms::y_spin_array[atom] = my;
			// 			atoms::z_spin_array[atom] = mz;
			// 			//		if (mat ==0) std::cout << atom << '\t' << mat << "\t" << atoms::z_coord_array[atom] << '\t' << atoms::x_spin_array[atom] << '\t' << atoms::y_spin_array[atom] << '\t' << atoms::z_spin_array[atom] << '\t' << sim::domain_wall_second_vector_x[mat] << '\t' << sim::domain_wall_second_vector_y[mat] << '\t' << sim::domain_wall_second_vector_z[mat] << '\t' <<std::endl;
			// 		}
			// 	}
			// }
			// if (sim::domain_wall_axis == 1){
			//    for(int atom=0;atom<num_local_atoms;atom++){
			//       if (atoms::y_coord_array[atom] > cs::system_dimensions[1]*sim::domain_wall_position){
			// 			int mat = atoms::type_array[atom];
			// 			atoms::x_spin_array[atom] = sim::domain_wall_second_vector_x[mat];
			// 			atoms::y_spin_array[atom] = sim::domain_wall_second_vector_y[mat];
			// 		 atoms::z_spin_array[atom] = sim::domain_wall_second_vector_z[mat];
			//        }
			//    }
			// }
			//
			// if (sim::domain_wall_axis == 2){
			//    for(int atom=0;atom<num_local_atoms;atom++){
			//       if (atoms::z_coord_array[atom] > cs::system_dimensions[2]*sim::domain_wall_position){
			// 			int mat = atoms::type_array[atom];
			// 			atoms::x_spin_array[atom] = sim::domain_wall_second_vector_x[mat];
			// 			atoms::y_spin_array[atom] = sim::domain_wall_second_vector_y[mat];
			// 		  atoms::z_spin_array[atom] = sim::domain_wall_second_vector_z[mat];
			//        }
			//    }
			// }
		}
		else {
			if (sim::domain_wall_axis == 0) {
				//90 degree 
				if(sim::domain_wall_angle == 0) {
					for(int atom=0;atom<num_local_atoms;atom++) {
					//		std::cout <<atom << '\t' <<  atoms::x_coord_array[atom] << "\t" << cs::system_dimensions[0]*sim::domain_wall_position -sim::domain_wall_width/2.0 << std::endl;
						// if (atoms::x_coord_array[atom] > cs::system_dimensions[0]*sim::domain_wall_position -sim::domain_wall_width*30.0) {
							int mat = atoms::type_array[atom]-1;
						//	double mod = 1.0;///sqrt(pos_x*pos_x + pos_y*pos_y + 1.0);

						//	double theta = std::atan(exp(-1.0*(atoms::x_coord_array[atom] - sim::domain_wall_position)/sim::domain_wall_width)) -M_PI*0.25;
						//	double pos_x = std::cos(theta)*mod;
						//	double pos_y = std::sin(theta)*mod; //std::tanh(pos-M_PI/4.0)/(std::cosh(pos-M_PI/4.0)*std::cosh(pos-M_PI/4.0));
					
						//	atoms::x_spin_array[atom] =  (mat==1?-1.0:1.0)* pos_x;
						//	atoms::y_spin_array[atom] =  (mat==1?-1.0:1.0)* pos_y;
						//	atoms::z_spin_array[atom] = 0.0;// (sim::domain_wall_second_vector_z[mat] - atoms::z_spin_array[atom])*theta;
						// int cat = atoms::sublayer_array[atom];
						// if(cat > 5) std::cout << cat << std::endl;
						int x_cell = (atoms::x_coord_array[atom]+0.0001)/sim::domain_wall_discretisation[0];
						int y_cell = (atoms::y_coord_array[atom]+0.0001)/sim::domain_wall_discretisation[1];
						int z_cell = (atoms::z_coord_array[atom]+0.0001)/sim::domain_wall_discretisation[2];
				//				std::cout << atom << '\t' << mat << '\t' << cell << "\t" << atom_to_cell_array.size() << "\t" << num_atoms_in_cell.size() << '\t' << num_dw_cells*mat + cell << std::endl;
						int cell = z_cell*num_dw_cells_x*num_dw_cells_y*num_mag_types + y_cell*num_dw_cells_x*num_mag_types + x_cell*num_mag_types + mat;
						//          1* + 
					//	std::cout << cell << "\t";
						atom_to_cell_array[atom] = cell;
						cell_to_lattice_array[3*cell + 0] = x_cell;
						cell_to_lattice_array[3*cell + 1] = y_cell;
						cell_to_lattice_array[3*cell + 2] = z_cell;

						num_atoms_in_cell[cell]++;
					//	avg_atom_num ++;
						// if( z_cell > 5) {
						// 	std::cout << cell << ", " << z_cell << ", " << atoms::z_coord_array[atom]+0.001 << ", " << sim::domain_wall_discretisation[2] << std::endl;
						// }
						// if (num_atoms_in_cell[num_dw_cells*cat + cell] < 0) std::cout << cell << ", " << atoms::x_coord_array[atom] << ", " << atoms::y_coord_array[atom] << ", " << atoms::z_coord_array[atom] << std::endl;
					}
				}
			}
		}
		/*
		int count = 0;
		// std::cout << sim::anti_PBC[0] << '\t' <<  sim::anti_PBC[1] << '\t' <<  sim::anti_PBC[2] <<std::endl;
		if (sim::anti_PBC[0] || sim::anti_PBC[1] || sim::anti_PBC[2]){
		  for (int atom = 0; atom <num_local_atoms; atom++){
		    const int start = atoms::neighbour_list_start_index[atom];
		    const int end   = atoms::neighbour_list_end_index[atom] + 1;
		    for(int nn=start;nn<end;nn++){
		 	 const int natom = atoms::neighbour_list_array[nn];
		 	 bool edge = false;
		 	 if (sim::anti_PBC[0] == true) {
		 	   #ifdef MPICF
		 	   if (atoms::x_coord_array[natom] < -0.01 || atoms::x_coord_array[natom] > (cs::system_dimensions[0]-0.01) || atoms::x_coord_array[atom] < -0.01 || atoms::x_coord_array[atom] > (cs::system_dimensions[0]-0.01)){
					//  std::cout << atom << '\t' << natom  << "\t" << atoms::x_coord_array[atom] << '\t' <<atoms::x_coord_array[natom] << '\t' <<  std::endl;
		
					 edge = true;
					 count++;
		 	   }
		 	   #else
		 	   const double dx = (atoms::x_coord_array[atom] - atoms::x_coord_array[natom])*(atoms::x_coord_array[atom] - atoms::x_coord_array[natom]);
		 	   if (dx > (cs::system_dimensions[0]-10)*(cs::system_dimensions[0]-10)){
		 	     edge = true;
		 	   }
		 	   #endif
		 	 }
		 	 if (sim::anti_PBC[1] == true){
		 	   //std::cout << "y" << std::endl;
		 	   #ifdef MPICF
		 	   if (atoms::y_coord_array[natom] < 0 || atoms::y_coord_array[natom] > (cs::system_dimensions[1]-1.665) || atoms::y_coord_array[atom] < 0 || atoms::y_coord_array[atom] > (cs::system_dimensions[1]-1.665)){
		 	     edge = true;
		 	   }
		 	   #else
		 	   const double dy = (atoms::y_coord_array[atom] - atoms::y_coord_array[natom])*(atoms::y_coord_array[atom] - atoms::y_coord_array[natom]);
		 	   if (dy > (cs::system_dimensions[1]-5)*(cs::system_dimensions[1]-5)){
		 	     edge = true;
		 	   }
		 	   #endif
		
		 	 }
		 	 if (sim::anti_PBC[2] == true){
		 	   //std::cout << "z" << std::endl;
		 	   #ifdef MPICF
		 	   if (atoms::z_coord_array[natom] < 0 || atoms::z_coord_array[natom] > cs::system_dimensions[2] || atoms::z_coord_array[atom] < 0 || atoms::z_coord_array[atom] > cs::system_dimensions[2]){
		 	     edge = true;
		 	   }
		 	   #else
		 	   const double dz = (atoms::z_coord_array[atom] - atoms::z_coord_array[natom])*(atoms::z_coord_array[atom] - atoms::z_coord_array[natom]);
		 	   if (dz > (cs::system_dimensions[2]-5)*(cs::system_dimensions[2]-5)){
		 	     edge = true;
		 	   }
		 	   #endif
		 	 }
			// 	std::cout << edge << std::endl;
		
		 	 	 if (edge == true){
		 		   	//   std::cout <<  exchange::internal::exchange_type << '\t' << edge << '\t' << nn << "\t" << exchange::tensorial << std::endl;
		 		     switch(exchange::internal::exchange_type){
		 	   case exchange::isotropic:
		 	     //std::cout << "iso" <<std::endl;
		 	     	atoms::i_exchange_list[nn].Jij = -1.0*atoms::i_exchange_list[nn].Jij;
						break;
		 	   case exchange::vectorial:
				 	// std::cout << atom << '\t' << natom << '\t' << atoms::v_exchange_list[nn].Jij[0] << "\t" << atoms::x_coord_array[atom] << '\t' <<atoms::x_coord_array[natom] << '\t' <<  std::endl;
		 	     atoms::v_exchange_list[nn].Jij[0] = -1.0*atoms::v_exchange_list[nn].Jij[0];
		 	     atoms::v_exchange_list[nn].Jij[1] = -1.0*atoms::v_exchange_list[nn].Jij[1];
		 	     atoms::v_exchange_list[nn].Jij[2] = -1.0*atoms::v_exchange_list[nn].Jij[2];
					 break;
		 	   case exchange::tensorial:
		 	    //std::cout << "tensor"<<std::endl;
			if(atoms::t_exchange_list[nn].Jij[0][0] < -300 ) {	
				// std::cout << atoms::t_exchange_list[nn].Jij[0][0] << std::endl;
		 	    //  if(atoms::x_coord_array[atom] > cs::system_dimensions[0]-1.665) atoms::t_exchange_list[nn].Jij[1][0] = -1.0*atoms::t_exchange_list[nn].Jij[0][0];
		 	    //  else atoms::t_exchange_list[nn].Jij[1][0] = 1.0*atoms::t_exchange_list[nn].Jij[0][0];
				//  if(atoms::x_coord_array[atom] > cs::system_dimensions[0]-1.665) atoms::t_exchange_list[nn].Jij[0][1] = -1.0*atoms::t_exchange_list[nn].Jij[0][0];
		 	    //  else atoms::t_exchange_list[nn].Jij[0][1] = 1.0*atoms::t_exchange_list[nn].Jij[0][0];
			 	atoms::t_exchange_list[nn].Jij[0][1] = 1.5*atoms::t_exchange_list[nn].Jij[0][0];
				atoms::t_exchange_list[nn].Jij[1][0] = 1.5*atoms::t_exchange_list[nn].Jij[0][0];
		 	     atoms::t_exchange_list[nn].Jij[0][2] = 0.0;//-1.0*atoms::t_exchange_list[nn].Jij[0][2];
		 	     atoms::t_exchange_list[nn].Jij[1][2] = 0.0;//-1.0*atoms::t_exchange_list[nn].Jij[1][2];
		 	     atoms::t_exchange_list[nn].Jij[2][0] = 0.0;//-1.0*atoms::t_exchange_list[nn].Jij[2][0];
		 	     atoms::t_exchange_list[nn].Jij[2][1] = 0.0;//-1.0*atoms::t_exchange_list[nn].Jij[2][1];
				 atoms::t_exchange_list[nn].Jij[0][0] = 1.0*atoms::t_exchange_list[nn].Jij[0][0];
				 atoms::t_exchange_list[nn].Jij[1][1] = 1.0*atoms::t_exchange_list[nn].Jij[1][1];
		 	     atoms::t_exchange_list[nn].Jij[2][2] = 1.0*atoms::t_exchange_list[nn].Jij[2][2];
				 if(atoms::x_coord_array[atom] > -0.01) {
				atoms::t_exchange_list[nn].Jij[1][0] *= -0.50;//*atoms::t_exchange_list[nn].Jij[1][1];
				atoms::t_exchange_list[nn].Jij[0][0] *= 1.0;//*atoms::t_exchange_list[nn].Jij[0][0];
				// atoms::t_exchange_list[nn].Jij[2][2] *= 1.0;//*atoms::t_exchange_list[nn].Jij[2][2];
				atoms::t_exchange_list[nn].Jij[1][1] *= 1.0;//*atoms::t_exchange_list[nn].Jij[2][2];
				atoms::t_exchange_list[nn].Jij[0][1] *= -0.50;//*atoms::t_exchange_list[nn].Jij[0][0];
				}
				if(atoms::x_coord_array[natom] < -0.01) {
				atoms::t_exchange_list[nn].Jij[1][0] *= -0.5;//*atoms::t_exchange_list[nn].Jij[1][1];
				atoms::t_exchange_list[nn].Jij[0][0] *= -1.0;//*atoms::t_exchange_list[nn].Jij[0][0];
				// atoms::t_exchange_list[nn].Jij[2][2] *= 1.0;//*atoms::t_exchange_list[nn].Jij[2][2];
				atoms::t_exchange_list[nn].Jij[1][1] *= -1.0;//*atoms::t_exchange_list[nn].Jij[2][2];
				atoms::t_exchange_list[nn].Jij[0][1] *= -0.500;//*atoms::t_exchange_list[nn].Jij[0][0];
				}
					 break;
			} else if(atoms::t_exchange_list[nn].Jij[0][0] > 90 ) {
				// std::cout << atoms::t_exchange_list[nn].Jij[0][0] << std::endl;
				// if(atoms::x_coord_array[atom] > cs::system_dimensions[0]-1.665) atoms::t_exchange_list[nn].Jij[0][1] = -1.0*atoms::t_exchange_list[nn].Jij[0][0];
				// else atoms::t_exchange_list[nn].Jij[0][1] = 1.0*atoms::t_exchange_list[nn].Jij[0][0];
				atoms::t_exchange_list[nn].Jij[0][1] = 1.0*atoms::t_exchange_list[nn].Jij[0][0];
		 	     atoms::t_exchange_list[nn].Jij[0][2] = 0.0;//-1.0*atoms::t_exchange_list[nn].Jij[0][2];
		 	    // atoms::t_exchange_list[nn].Jij[1][0] = -1.0*atoms::t_exchange_list[nn].Jij[1][1];
		 	     atoms::t_exchange_list[nn].Jij[1][2] = 0.0;//-1.0*atoms::t_exchange_list[nn].Jij[1][2];
		 	     atoms::t_exchange_list[nn].Jij[2][0] = 0.0;//-1.0*atoms::t_exchange_list[nn].Jij[2][0];
		 	     atoms::t_exchange_list[nn].Jij[2][1] = 0.0;//-1.0*atoms::t_exchange_list[nn].Jij[2][1];
				 atoms::t_exchange_list[nn].Jij[1][0] = 1.00*atoms::t_exchange_list[nn].Jij[0][0];
				 atoms::t_exchange_list[nn].Jij[0][0] = 1.0*atoms::t_exchange_list[nn].Jij[0][0];
				 atoms::t_exchange_list[nn].Jij[1][1] = 1.00*atoms::t_exchange_list[nn].Jij[1][1];
		 	     atoms::t_exchange_list[nn].Jij[2][2] = 1.0*atoms::t_exchange_list[nn].Jij[2][2];
				 	 if(atoms::x_coord_array[natom] < -0.01) {
				atoms::t_exchange_list[nn].Jij[1][0] *= -1.0;//*atoms::t_exchange_list[nn].Jij[1][1];
				atoms::t_exchange_list[nn].Jij[0][0] *= 1.0;//*atoms::t_exchange_list[nn].Jij[0][0];
				atoms::t_exchange_list[nn].Jij[2][2] *= 1.0;//*atoms::t_exchange_list[nn].Jij[2][2];
				atoms::t_exchange_list[nn].Jij[1][1] *= 1.0;//*atoms::t_exchange_list[nn].Jij[2][2];
				atoms::t_exchange_list[nn].Jij[0][1] *= -1;//*atoms::t_exchange_list[nn].Jij[0][0];
				}
				if(atoms::x_coord_array[atom] > -0.01) {
				atoms::t_exchange_list[nn].Jij[1][0] *= -1.0;//*atoms::t_exchange_list[nn].Jij[1][1];
				atoms::t_exchange_list[nn].Jij[0][0] *= 1.0;//*atoms::t_exchange_list[nn].Jij[0][0];
				atoms::t_exchange_list[nn].Jij[2][2] *= 1.0;//*atoms::t_exchange_list[nn].Jij[2][2];
				atoms::t_exchange_list[nn].Jij[1][1] *= 1.0;//*atoms::t_exchange_list[nn].Jij[2][2];
				atoms::t_exchange_list[nn].Jij[0][1] *= -1.0;//*atoms::t_exchange_list[nn].Jij[0][0];
				}
					 break;
			} 
			
		 	     }
		 	   }
		    }
		  }
		
		
		std::cout << vmpi::my_rank << ": " << count << std::endl;
		*/
		//works out which atoms are in which cells and sets cells based on wether the domain wall
		//is along x or y or z
		// if (sim::domain_wall_axis == 0){
		// 	for(int atom=0;atom<num_local_atoms;atom++){
		// 		int cat = atoms::category_array[atom];
		// 		int cell = atoms::x_coord_array[atom]/sim::domain_wall_discretisation;
		// 		//				std::cout << atom << '\t' << mat << '\t' << cell << "\t" << atom_to_cell_array.size() << "\t" << num_atoms_in_cell.size() << '\t' << num_dw_cells*mat + cell << std::endl;
		// 		atom_to_cell_array[atom] = cell;
		// 		num_atoms_in_cell[num_dw_cells*cat + cell] ++;
		// 		//	if (cell > num_dw_cells) std::cout << atoms::x_coord_array[atom] << '\t' << sim::domain_wall_discretisation <<std::endl;
		// 	}
		// }
		// if (sim::domain_wall_axis == 1){
		// 	for(int atom=0;atom<num_local_atoms;atom++){
		// 		int cat = atoms::sublayer_array[atom];
		// 		int cell = atoms::y_coord_array[atom]/sim::domain_wall_discretisation;
		// 		//				std::cout << atom << '\t' << mat << '\t' << cell << "\t" << atom_to_cell_array.size() << "\t" << num_atoms_in_cell.size() << '\t' << num_dw_cells*mat + cell << std::endl;
		// 		atom_to_cell_array[atom] = cell;
		// 		num_atoms_in_cell[num_dw_cells*cat + cell] ++;
		// 		// if(num_atoms_in_cell[num_dw_cells*cat + cell] < 0 || num_atoms_in_cell[num_dw_cells*cat + cell] > 20) std::cout << cell << ", " << atoms::x_coord_array[atom] << ", " << atoms::y_coord_array[atom] << ", " << atoms::z_coord_array[atom] << std::endl;
		// 		//	if (cell > num_dw_cells) std::cout << atoms::x_coord_array[atom] << '\t' << sim::domain_wall_discretisation <<std::endl;
		// 	}
		// }
		// if (sim::domain_wall_axis == 1){
		// 	for(int atom=0;atom<num_local_atoms;atom++){
		// 		int mat = atoms::type_array[atom];
		// 		int cell = atoms::y_coord_array[atom]/sim::domain_wall_discretisation;
		// 		atom_to_cell_array[atom] = cell;
		// 		num_atoms_in_cell[num_dw_cells*mat + cell] ++;
		// 	}
		// }
		// if (sim::domain_wall_axis == 2){
		// 	for(int atom=0;atom<num_local_atoms;atom++){
		// 		int mat = atoms::type_array[atom];
		// 		int cell = atoms::z_coord_array[atom]/sim::domain_wall_discretisation;
		// 		atom_to_cell_array[atom] = cell;
		// 		num_atoms_in_cell[num_dw_cells*mat + cell] ++;
		// 	}
		// }

		#ifdef MPICF
			MPI_Allreduce(MPI_IN_PLACE, &num_atoms_in_cell[0],  num_dw_cells,    MPI_INT,    MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(MPI_IN_PLACE, &cell_to_lattice_array[0],  3*num_dw_cells,    MPI_INT,    MPI_SUM, MPI_COMM_WORLD);
		#endif

		int avg_atoms_cells = 0;
		for(int i = 0; i < num_dw_cells; i++) {
			//if(num_atoms_in_cell[i] == 0) std::cout << i << std::endl;
			
			avg_atoms_cells += num_atoms_in_cell[i];
		}

			std::cout << "avg atoms per cell: " << avg_atoms_cells/num_dw_cells  << std::endl;
		// Set equilibration temperature only if continue checkpoint not loaded
		if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
		else{
			// Set equilibration temperature
			sim::temperature=sim::Teq;
		}
		// for (int cat = 1; cat < num_categories; cat ++) {
		// 		if( cat == 3) continue;
			for (int cell = 0; cell < 6*num_dw_cells; cell++) {
			
				// std::cout << mat << '\t' << cell << "\t" << mag_x[num_dw_cells*mat + cell] << "\t" << mag_y[num_dw_cells*mat + cell] << "\t" << mag_z[num_dw_cells*mat + cell] << std::endl;
					mag[cell] = 0.0;
					// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +1] = 0.0;
					// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +2] = 0.0;
					// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +3] = 0.0;
					// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +4] = 0.0;
					// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +5] = 0.0;

				// topological_charge[num_dw_cells*cat + cell] = 0.0;
			}

			// #ifdef MPICF
			// MPI_Allreduce(MPI_IN_PLACE, &mag_x[0],     num_dw_cells*mp::num_materials,    MPI_DOUBLE,    MPI_MIN, MPI_COMM_WORLD);
			// MPI_Allreduce(MPI_IN_PLACE, &mag_y[0],     num_dw_cells*mp::num_materials,    MPI_DOUBLE,    MPI_MIN, MPI_COMM_WORLD);
			// MPI_Allreduce(MPI_IN_PLACE, &mag_z[0],     num_dw_cells*mp::num_materials,    MPI_DOUBLE,    MPI_MIN, MPI_COMM_WORLD);
			// //MPI_Allreduce(MPI_IN_PLACE, &num_atoms_in_cell[0],     num_dw_cells*mp::num_materials,    MPI_INT,    MPI_MIN, MPI_COMM_WORLD);
			// #endif
		
		
			for(int atom=0;atom<num_local_atoms;atom++) {
				int cell = atom_to_cell_array[atom];
				//int cat = atoms::sublayer_array[atom];
				int mat = atoms::type_array[atom];
				double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
				mag[num_mag_cat*cell+ 0] += S[0];
				mag[num_mag_cat*cell+ 1] += S[1];
				mag[num_mag_cat*cell+ 2] += S[2];
				mag[num_mag_cat*cell+ 3] += exchange::single_spin_energy(atom, S[0], S[1], S[2]);
				mag[num_mag_cat*cell+ 4] += anisotropy::single_spin_energy(atom,mat , S[0], S[1], S[2], sim::temperature);
				mag[num_mag_cat*cell+ 5] += (cos(atan2(S[1],S[0])*2.0))*sim::internal::lot_lt_z[mat];
				// topological_charge[num_dw_cells*cat + cell] += M_PI*0.5 - atan2(atoms::y_spin_array[atom], atoms::x_spin_array[atom]);
			}


			#ifdef MPICF
			MPI_Allreduce(MPI_IN_PLACE, &mag[0], num_mag_cat*num_dw_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
		//	MPI_Allreduce(MPI_IN_PLACE, &mag_y[0],     num_dw_cells*num_categories,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
		//	MPI_Allreduce(MPI_IN_PLACE, &mag_z[0],     num_dw_cells*num_categories,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
			// MPI_Allreduce(MPI_IN_PLACE, &topological_charge[0],     num_dw_cells*num_categories,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
			#endif
			//	 std::cout << "a" <<std::endl;
		// if(vmpi::my_rank==0) { 
		// 	std::ofstream myfile;
		// 	string filename = "/dw/dw-0.txt";
		// 	myfile.open (string(directory) + filename);
		// 	for (int cell = 0; cell < num_dw_cells; cell++){
		// 			if (num_atoms_in_cell[cell] > 0){
		// 				int num = num_atoms_in_cell[cell];
		// 				myfile << cell <<"\t" <<  cell_to_lattice_array[3*cell + 0] << '\t' <<  cell_to_lattice_array[3*cell + 1] << '\t' <<  cell_to_lattice_array[3*cell + 2] << '\t' << mag[num_mag_cat*cell] / num << "\t" << mag[num_mag_cat*cell + 1] /num << "\t" << mag[num_mag_cat*cell +2]/num  << "\t" << //magnetisation data
		// 							mag[num_mag_cat*cell +3]/num << "\t" << //exchange energy
		// 							mag[num_mag_cat*cell +4]/num << "\t" << //anisotropy energy
		// 							mag[num_mag_cat*cell +5]/num << "\n"; //LOT energy
		// 			}
		// 		}
			
		// 	myfile.close();
		// }
			// Output data
			vout::data();
		program::fractional_electric_field_strength = 0.0;
		// Equilibrate system
		while(sim::time<sim::equilibration_time) {

			sim::integrate(sim::partial_time);

			// Calculate magnetisation statistics
			stats::update();
				
				// for (int cell = 0; cell < num_dw_cells; cell++){
				
				// // std::cout << mat << '\t' << cell << "\t" << mag_x[num_dw_cells*mat + cell] << "\t" << mag_y[num_dw_cells*mat + cell] << "\t" << mag_z[num_dw_cells*mat + cell] << std::endl;
				// 	mag[cell] = 0.0;
				// 	// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +1] = 0.0;
				// 	// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +2] = 0.0;
				// 	// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +3] = 0.0;
				// 	// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +4] = 0.0;
				// 	// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +5] = 0.0;

				// 	// topological_charge[num_dw_cells*cat + cell] = 0.0;
				// }
			

			// #ifdef MPICF
			// MPI_Allreduce(MPI_IN_PLACE, &mag_x[0],     num_dw_cells*mp::num_materials,    MPI_DOUBLE,    MPI_MIN, MPI_COMM_WORLD);
			// MPI_Allreduce(MPI_IN_PLACE, &mag_y[0],     num_dw_cells*mp::num_materials,    MPI_DOUBLE,    MPI_MIN, MPI_COMM_WORLD);
			// MPI_Allreduce(MPI_IN_PLACE, &mag_z[0],     num_dw_cells*mp::num_materials,    MPI_DOUBLE,    MPI_MIN, MPI_COMM_WORLD);
			// //MPI_Allreduce(MPI_IN_PLACE, &num_atoms_in_cell[0],     num_dw_cells*mp::num_materials,    MPI_INT,    MPI_MIN, MPI_COMM_WORLD);
			// #endif
			
			// Output data
		
			// Output data
			vout::data();
		}

		///	std::cout << "end equilbration" << std::endl;


		// Set temperature and reset stats only if continue checkpoint not loaded
		if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag) {}
		else{

			// set simulation temperature
			sim::temperature = temp;

			// Reset mean magnetisation counters
			stats::reset();

		}
		
		// Perform Time Series
		std::ofstream dw_pos;
	if(vmpi::my_rank==0){
		
			string filename = "/dw-pos.txt";
			dw_pos.open (string(directory) + filename);
			if(!dw_pos.is_open()) {
				std::cerr << "Fatal dw directory error for dw-" + std::to_string(sim::time) << std::endl;
			}
	}
		while(sim::time<sim::equilibration_time+sim::total_time) {

			double ftime = mp::dt_SI*double(sim::time-sim::equilibration_time);
			const double i_pump_time = 1.0/sim::pump_time;
  		 	double reduced_time = (ftime-1.5*sim::pump_time)*i_pump_time;
   			const double four_ln_2 = 2.77258872224; // 4 ln 2
   			double gaussian = exp(-four_ln_2*reduced_time*reduced_time);
			if(sim::enable_laser_torque_fields) {
				if(ftime < 1.5*sim::pump_time) sim::laser_torque_strength = gaussian;
				else if (ftime < 1.5*sim::pump_time + sim::double_pump_delay) sim::laser_torque_strength = 1.0;
				else {
					reduced_time = (ftime-1.5*sim::pump_time-sim::double_pump_delay)*i_pump_time;
					gaussian = exp(-four_ln_2*reduced_time*reduced_time);
					sim::laser_torque_strength = gaussian;
				}
			}
   			const double two_delta_sqrt_pi_ln_2 = 9394372.787;
	
   			const double pump= two_delta_sqrt_pi_ln_2*sim::pump_power*gaussian*i_pump_time;
   			const double Te = sim::TTTe;
   			const double Tp = sim::TTTp;
   			const double G  = sim::TTG;
   			const double Ce = sim::TTCe;
   			const double Cl = sim::TTCl;
   			const double dt = mp::dt_SI;

			// integrate two temperature model (floor in free elecron approximation (c prop to T) for low temperatures)
			if(Te>1.0) sim::TTTe = (-G*(Te-Tp)+pump)*dt/(Ce*Te) + Te;
			else sim::TTTe =       (-G*(Te-Tp)+pump)*dt/Ce + Te;
			sim::TTTp =            ( G*(Te-Tp)     )*dt/Cl + Tp - (Tp-sim::Teq)*sim::HeatSinkCouplingConstant*dt;

			  double time_from_start = mp::dt_SI * double(sim::time-sim::equilibration_time);

   		if(time_from_start < program::internal::electrical_pulse_rise_time ) {
      		program::fractional_electric_field_strength = time_from_start / program::internal::electrical_pulse_rise_time;
			// if(program::fractional_electric_field_strength > 0.0) std::cout << program::fractional_electric_field_strength << std::endl;
   		}
   // implement continuous current
   		else if(time_from_start < (program::internal::electrical_pulse_rise_time + program::internal::electrical_pulse_time) ){
     		 program::fractional_electric_field_strength = 1.0;
			//  std::cout << program::fractional_electric_field_strength << std::endl;
   		}
   // implement fall time
   		else if(time_from_start < (program::internal::electrical_pulse_rise_time + program::internal::electrical_pulse_time + program::internal::electrical_pulse_fall_time)) {
    	  const double fractional_fall_time = time_from_start - (program::internal::electrical_pulse_rise_time + program::internal::electrical_pulse_time);
    	  program::fractional_electric_field_strength = 1.0 - fractional_fall_time / program::internal::electrical_pulse_fall_time;
		//   if(program  ::fractional_electric_field_strength) std::cout << program::fractional_electric_field_strength << std::endl;
   		}
   // after pulse current = 0
   		else{
      		program::fractional_electric_field_strength = 0.0;
   		}
				//#pragma vector
				for (int cell = 0; cell < num_dw_cells; cell++) {
				
					// std::cout << mat << '\t' << cell << "\t" << mag_x[num_dw_cells*mat + cell] << "\t" << mag_y[num_dw_cells*mat + cell] << "\t" << mag_z[num_dw_cells*mat + cell] << std::endl;
					mag[num_mag_cat* cell + 0 ] = 0.0;
					mag[num_mag_cat* cell + 1 ] = 0.0;
					mag[num_mag_cat* cell + 2 ] = 0.0;
					mag[num_mag_cat* cell + 3 ] = 0.0;
					mag[num_mag_cat* cell + 4 ] = 0.0;
					mag[num_mag_cat* cell + 5 ] = 0.0;
					// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +1] = 0.0;
					// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +2] = 0.0;
					// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +3] = 0.0;
					// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +4] = 0.0;
					// mag[num_mag_cat*num_dw_cells*cat + num_mag_cat*cell +5] = 0.0;

					// topological_charge[num_dw_cells*cat + cell] = 0.0;
				
			}

			// #ifdef MPICF
			// MPI_Allreduce(MPI_IN_PLACE, &mag_x[0],     num_dw_cells*mp::num_materials,    MPI_DOUBLE,    MPI_MIN, MPI_COMM_WORLD);
			// MPI_Allreduce(MPI_IN_PLACE, &mag_y[0],     num_dw_cells*mp::num_materials,    MPI_DOUBLE,    MPI_MIN, MPI_COMM_WORLD);
			// MPI_Allreduce(MPI_IN_PLACE, &mag_z[0],     num_dw_cells*mp::num_materials,    MPI_DOUBLE,    MPI_MIN, MPI_COMM_WORLD);
			// // MPI_Allreduce(MPI_IN_PLACE, &num_atoms_in_cell[0],     num_dw_cells*mp::num_materials,    MPI_INT,    MPI_MIN, MPI_COMM_WORLD);
			// #endif

			// Integrate system
			sim::integrate(sim::partial_time);

			// Calculate magnetisation statistics
			stats::update();

			// std::ofstream myfile_1;
		
			for(int atom=0;atom<num_local_atoms;atom++) {
				int cell = atom_to_cell_array[atom];
			//	int cat = atoms::sublayer_array[atom];
				int mat = atoms::type_array[atom];
				double S[3] = {atoms::x_spin_array[atom], atoms::y_spin_array[atom],atoms::z_spin_array[atom]};
				mag[num_mag_cat*cell]   += S[0];
				mag[num_mag_cat*cell+1] += S[1];
				mag[num_mag_cat*cell+2] += S[2];
				mag[num_mag_cat*cell+3] += exchange::single_spin_energy(atom, S[0], S[1], S[2]);
				mag[num_mag_cat*cell+4] += anisotropy::single_spin_energy(atom,mat , S[0], S[1], S[2], sim::temperature);
				mag[num_mag_cat*cell+5] += sim::laser_torque_strength*(cos(atan2(S[1],S[0])*2.0))*sim::internal::lot_lt_z[mat];

				// topological_charge[num_dw_cells*cat + cell] += M_PI*1.5 - (M_PI+atan2(atoms::y_spin_array[atom], atoms::x_spin_array[atom]));
			}

			#ifdef MPICF
			MPI_Allreduce(MPI_IN_PLACE, &mag[0], num_mag_cat*num_dw_cells,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
		//	MPI_Allreduce(MPI_IN_PLACE, &mag_y[0],     num_dw_cells*num_categories,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
		//	MPI_Allreduce(MPI_IN_PLACE, &mag_z[0],     num_dw_cells*num_categories,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
			// MPI_Allreduce(MPI_IN_PLACE, &topological_charge[0],     num_dw_cells*num_categories,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
			#endif
		
			
			//sim::domain_wall_velocity = 1.0e-10* std::abs((minima_x[0]/minima_x[1]) - 10.0*sim::domain_wall_centre)/(sim::partial_time*mp::dt_SI); //A/0.1ft -> e16/ e10 -> m/s
			//sim::domain_wall_centre  = 0.1*minima_x[0]/minima_x[1];
		bool print = (sim::time % 500 == 0) ? true : false;
		
		if(vmpi::my_rank == 0) {
			double avg_topological_charge_acc = 0.0;
			std::vector< double> domain_tracks;
			int domain_counter = 0;
			std::ofstream dw_res;
			string filename = "/dw/dw-" + std::to_string(sim::time) + ".txt";
			if(print) dw_res.open (string(directory) + filename);
			if(print) if(!dw_res.is_open()) {
				std::cerr << "Fatal dw directory error for dw-" + std::to_string(sim::time) << std::endl;
			}
			
			
			dw_res.precision(10);
    		dw_res << std::scientific;
			int mat_type = 0;
			for (int z_cell = 0; z_cell < num_dw_cells_z-1; z_cell++) {
				for(int y_cell = 0; y_cell < num_dw_cells_y-1; y_cell++) {
					double mag_x_1 = 1.0/sqrt(2.0);
					double mag_y_1 = 1.0/sqrt(2.0);
					double mag_x_2 = 1.0/sqrt(2.0);
					double mag_y_2 = 1.0/sqrt(2.0);
					double d_topological_charge_acc = 0.0;
					double d_topological_charge = 0.0;

					if(print) {
					int cell = z_cell*num_dw_cells_x*num_dw_cells_y*num_mag_types + y_cell*num_dw_cells_x*num_mag_types + 1*num_mag_types + mat_type;	
					double num = num_atoms_in_cell[cell];
					dw_res << 0 << '\t' << y_cell << "\t" << z_cell << "\t" << mag_x_1 << "\t" << mag_y_1 << "\t" << "\n";
					dw_res << cell_to_lattice_array[3*cell+0]*sim::domain_wall_discretisation[0]  << '\t' << cell_to_lattice_array[3*cell + 1] << '\t' << cell_to_lattice_array[3*cell + 2]<< '\t' <<\
							        mag[num_mag_cat*cell] / num << "\t" << mag[num_mag_cat*cell + 1] /num << "\t" << //magnetisation data
									mag[num_mag_cat*cell +3]/num << "\t" << //exchange energy
									mag[num_mag_cat*cell +4]/num << "\t" << //anisotropy energy
									mag[num_mag_cat*cell +5]/num << "\t" << //LOT energy
									d_topological_charge << "\t"  << d_topological_charge_acc  <<  "\n";
					}
				for(int x_cell = 1; x_cell < num_dw_cells_x; x_cell++) {
						int cell = z_cell*num_dw_cells_x*num_dw_cells_y*num_mag_types + y_cell*num_dw_cells_x*num_mag_types + x_cell*num_mag_types + mat_type;	
					//std::cout << cell << std::endl;
					int num = num_atoms_in_cell[cell];
					if(num == 0) continue;
					//else std::cout << cell_to_lattice_array[3*cell+0]*sim::domain_wall_discretisation[0]  << '\t' << cell_to_lattice_array[3*cell + 1]*sim::domain_wall_discretisation[1] << '\t' << cell_to_lattice_array[3*cell + 2]*sim::domain_wall_discretisation[2] << '\t' <<\
							        // mag[num_mag_cat*cell] / num << "\t" << mag[num_mag_cat*cell + 1] /num << "\t" << //magnetisation data/
									// mag[num_mag_cat*cell +3]/num << "\t" << //exchange energy/
									// mag[num_mag_cat*cell +4]/num << "\t" << //anisotropy energy/
									// mag[num_mag_cat*cell +5]/num << "\t" << //LOT energy/
									// d_topological_charge << "\t"  << d_topological_charge_acc  <<  std::endl;
					mag_x_1 = mag[num_mag_cat*(cell-2) + 0] / double(num_atoms_in_cell[cell-2]);
					mag_y_1 = mag[num_mag_cat*(cell-2) + 1] /double(num_atoms_in_cell[cell-2]);
					mag_x_2 = mag[num_mag_cat*cell   +0] /double(num_atoms_in_cell[ cell]);
					mag_y_2 = mag[(num_mag_cat*cell) +1] / double(num_atoms_in_cell[cell]);
					
					if((mag_y_2 != 0.0 || mag_x_2 != 0.0) && (mag_y_1!= 0.0 || mag_x_1 != 0.0) ) {
						d_topological_charge =  atan2(mag_y_1, mag_x_1) -atan2(mag_y_2, mag_x_2);
						if(d_topological_charge > M_PI*0.5) d_topological_charge -= 2.0*M_PI;
						else if (d_topological_charge < -M_PI*0.5) d_topological_charge += 2.0*M_PI;
					}
					//if(std::abs(d_x) > 1e-4 || std::abs(d_y) > 1e-4) topological_charge_1[cat] = atan2(d_y,d_x);
					if( (mag_x_1*mag_x_2 < 0.0) != (mag_y_1*mag_y_2 < 0.0) ) {
						  //90 degree change only
						double x = cell_to_lattice_array[3*cell+0]*sim::domain_wall_discretisation[0];// + sl_offset[cell_to_lattice_array[3*cell+2]];
						double y = cell_to_lattice_array[3*cell+1]*sim::domain_wall_discretisation[1];// + sl_offset[cell_to_lattice_array[3*cell+2]];
						double z = cell_to_lattice_array[3*cell+2]*sim::domain_wall_discretisation[2];
							domain_tracks.push_back(x);
							domain_tracks.push_back(y);
							domain_tracks.push_back(z);
							domain_counter++;
					//	std::cout << cell << ", " << x << ", " << mag_x_1 << ", " << mag_x_2 << ", " << mag_y_1 << ", " << mag_y_2 << std::endl;
					}
					d_topological_charge_acc += d_topological_charge;
					avg_topological_charge_acc += d_topological_charge_acc;
					//if(d_topological_charge < 1e-6) std::cout << d_topological_charge << std::endl;
					if (num_atoms_in_cell[cell] > 0 ) {
						if(print && (d_topological_charge > 1e-7)) {
							double num = num_atoms_in_cell[cell];

							dw_res << cell_to_lattice_array[3*cell+0]*sim::domain_wall_discretisation[0]  << '\t' << cell_to_lattice_array[3*cell + 1] << '\t' << cell_to_lattice_array[3*cell + 2]<< '\t' <<\
							        mag[num_mag_cat*cell] / num << "\t" << mag[num_mag_cat*cell + 1] /num << "\t" << //magnetisation data
									mag[num_mag_cat*cell +3]/num << "\t" << //exchange energy
									mag[num_mag_cat*cell +4]/num << "\t" << //anisotropy energy
									mag[num_mag_cat*cell +5]/num << "\t" << //LOT energy
									d_topological_charge << "\t"  << d_topological_charge_acc  <<  "\n";
							// if(mag_x_2 != mag_x_2 || mag_y_2 != mag_y_2) {
							// 	std::cout << cell_to_lattice_array[3*cell+0]*sim::domain_wall_discretisation[0]  << '\t' << cell_to_lattice_array[3*cell + 1]*sim::domain_wall_discretisation[1] << '\t' << cell_to_lattice_array[3*cell + 2]*sim::domain_wall_discretisation[2] << '\t' <<\
							//         mag[num_mag_cat*cell] / num << "\t" << mag[num_mag_cat*cell + 1] /num << "\t" << //magnetisation data
							// 		mag[num_mag_cat*cell +3]/num << "\t" << //exchange energy
							// 		mag[num_mag_cat*cell +4]/num << "\t" << //anisotropy energy
							// 		mag[num_mag_cat*cell +5]/num << "\t" << //LOT energy
							// 		d_topological_charge << "\t"  << d_topological_charge_acc  <<  std::endl;
							// }
						}
					}
				}
				if(print) {
					int cell = (num_dw_cells_z-2)*num_dw_cells_x*num_dw_cells_y*num_mag_types + (num_dw_cells_y-2)*num_dw_cells_x*num_mag_types + (num_dw_cells_x-1)*num_mag_types + mat_type;
					double num = num_atoms_in_cell[cell];	
					dw_res << cell_to_lattice_array[3*cell+0]*sim::domain_wall_discretisation[0]  << '\t' << cell_to_lattice_array[3*cell + 1] << '\t' << cell_to_lattice_array[3*cell + 2]<< '\t' <<\
							        mag[num_mag_cat*cell] / num << "\t" << mag[num_mag_cat*cell + 1] /num << "\t" << //magnetisation data
									mag[num_mag_cat*cell +3]/num << "\t" << //exchange energy
									mag[num_mag_cat*cell +4]/num << "\t" << //anisotropy energy
									mag[num_mag_cat*cell +5]/num << "\t" << //LOT energy
									d_topological_charge << "\t"  << d_topological_charge_acc  <<  "\n";
				}
			}
			}
			/* second layer only
			// d_topological_charge = 0.0;
			// // double topological_charge_1[num_categories] = {0.0};
			// d_topological_charge_acc = 0.0;
			// for (int cell = (1*num_dw_cells_x*num_dw_cells_y + 2); cell != (2*num_dw_cells_x*num_dw_cells_y -2); cell++) {
			// 		// if( cell_to_lattice_array[3*cell + 2] != 2) {
			// 		// 	avg_topological_charge_acc += d_topological_charge_acc;
			// 		// 	d_topological_charge = 0.0;
			//  		// 	d_topological_charge_acc = 0.0;
			// 		// 	continue;
			// 		// }
			// 		if(num_atoms_in_cell[cell] != 2) continue;
			// 		mag_x_1 = mag[num_mag_cat*(cell-1) + 0] / double(num_atoms_in_cell[cell-1]);
			// 		mag_y_1 = mag[num_mag_cat*(cell-1) + 1] /double(num_atoms_in_cell[cell-1]);
			// 		mag_x_2 = mag[num_mag_cat*cell   +0] /double(num_atoms_in_cell[ cell]);
			// 		mag_y_2 = mag[(num_mag_cat*cell) +1] / double(num_atoms_in_cell[cell]);
					
			// 		if((mag_y_2 != 0.0 || mag_x_2 != 0.0) && (mag_y_1!= 0.0 || mag_x_1 != 0.0) ) {
			// 			d_topological_charge =  atan2(mag_y_1, mag_x_1) -atan2(mag_y_2, mag_x_2);
			// 			if(d_topological_charge > M_PI*0.5) d_topological_charge -= 2.0*M_PI;
			// 			else if (d_topological_charge < -M_PI*0.5) d_topological_charge += 2.0*M_PI;
			// 		}
			// 		//if(std::abs(d_x) > 1e-4 || std::abs(d_y) > 1e-4) topological_charge_1[cat] = atan2(d_y,d_x);
			// 		if( (mag_x_1*mag_x_2 < 0.0) != (mag_y_1*mag_y_2 < 0.0) ) {
			// 			  //90 degree change only
			// 			double x = cell_to_lattice_array[3*cell+0]*sim::domain_wall_discretisation[0] + sl_offset[cell_to_lattice_array[3*cell+2]];
			// 			double y = cell_to_lattice_array[3*cell+1]*sim::domain_wall_discretisation[1] + sl_offset[cell_to_lattice_array[3*cell+2]];
			// 			double z = cell_to_lattice_array[3*cell+2]*sim::domain_wall_discretisation[2];
			// 				domain_tracks.push_back(x);
			// 				domain_tracks.push_back(y);
			// 				domain_tracks.push_back(z);
			// 				domain_counter++;
			// 		//	std::cout << cell << ", " << x << ", " << mag_x_1 << ", " << mag_x_2 << ", " << mag_y_1 << ", " << mag_y_2 << std::endl;
			// 		}
			// 		d_topological_charge_acc += d_topological_charge;
			// 		//if(d_topological_charge < 1e-6) std::cout << d_topological_charge << std::endl;
			// 		if (num_atoms_in_cell[cell] > 0 && std::abs(d_topological_charge) > 1e-6 ) {
						
			// 			if(print) {
			// 				double num = num_atoms_in_cell[cell];
			// 				dw_res.precision(10);
    		// 				dw_res << std::scientific;
			// 				dw_res << cell_to_lattice_array[3*cell+0]  << '\t' <<\
			// 				        mag[num_mag_cat*cell] / num << "\t" << mag[num_mag_cat*cell + 1] /num << "\t" << //magnetisation data
			// 						mag[num_mag_cat*cell +3]/num << "\t" << //exchange energy
			// 						mag[num_mag_cat*cell +4]/num << "\t" << //anisotropy energy
			// 						mag[num_mag_cat*cell +5]/num << "\t" << //LOT energy
			// 						d_topological_charge << "\t"  << d_topological_charge_acc  <<  "\n";
			// 			}
			// 		}
			// } */
			
			if(print) dw_res.close();
			dw_pos << sim::time << "\t" << avg_topological_charge_acc/(num_dw_cells_y*num_dw_cells_z) << "\t" << domain_counter/(num_dw_cells_y*num_dw_cells_z) << "\t";
			//std::sort(domain_tracks.begin(), domain_tracks.end(), std::greater<double>());
			for(int d = 0; d < domain_tracks.size(); d++) {
				if(domain_counter == 0) break;
				dw_pos << domain_tracks[d] << "\t";
				
			} dw_pos << std::endl;
		}
			// Output data
			
		

		// if(domain_counter == 4 && print) {
		// 		std::ofstream dw_3d;
		// 		string filename;
		// 		if(vmpi::my_rank = 0) {filename = "/dw/dw-3d-" + std::to_string(sim::time) + ".txt";
		// 			if(print) dw_3d.open (string(directory) + filename);
		// 			if(print) if(!dw_3d.is_open()) {
		// 				std::cerr << "Fatal dw directory error for dw-" + std::to_string(sim::time) << std::endl;
		// 			}
		// 		}
		// 		// if(vmpi::my_rank)
		// 		int cell = (domain_tracks[2] + sim::unit_cell_x*0.5)/sim::unit_cell_x;
		// 		int start_cell = cell - horizontal_range/2;
		// 		int end_cell = cell + horizontal_range/2;
		// 		for(int cat = 1; cat < num_categories; cat++) {
		// 			if(cat == 3) continue;
		// 			for(cell = start_cell; cell < end_cell; cell++) {
		// 				for(int y_cell = 0; y_cell < depth_range; y_cell++) {
		// 					for(int z_cell = 0; z_cell < depth_range; z_cell++) {
		// 						int atom = cell_to_atom_array[num_dw_cells*cat + cell][vertical_range*y_cell + z_cell];
		// 						mag_x[num_dw_cells*cat + cell] += atoms::x_spin_array[atom];
		// 						mag_y[num_dw_cells*cat + cell] += atoms::y_spin_array[atom];
		// 						mag_z[num_dw_cells*cat + cell] += atoms::z_spin_array[atom];
		// 					}
		// 				}
		// 			}
		// 		}
		// 		#ifdef MPICF 
		// 		MPI_Allreduce(MPI_IN_PLACE, &mag_x[0], num_dw_cells*num_categories,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
		// 		MPI_Allreduce(MPI_IN_PLACE, &mag_y[0],     num_dw_cells*num_categories,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
		// 		MPI_Allreduce(MPI_IN_PLACE, &mag_z[0],     num_dw_cells*num_categories,    MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
		// 		#endif
		// 	if(vmpi::my_rank==0) {
		// 		for(int cat = 1; cat < num_categories; cat++) {
		// 			if(cat == 3) continue;
		// 			for(cell = start_cell; cell < end_cell; cell++) {
		// 				for(int y_cell = 0; y_cell < depth_range; y_cell++) {
		// 					for(int z_cell = 0; z_cell < depth_range; z_cell++) {
		// 						if(num_atoms_in_cell[num_dw_cells*cat +cell] > 0) {
		// 							dw_3d << cat << ", " << cell << ", " << y_cell << ", " << z_cell << ", " << mag_x[num_dw_cells*cat + cell] << ", " << mag_x[num_dw_cells*cat + cell] << ", " << mag_x[num_dw_cells*cat + cell] << "\n";
		// 						}
		// 					}
		// 				}
		// 			}
		// 		}
		// 	}	
		// 	dw_3d.close();
		// 	}
		// }
		
		vout::data();
		} 
		
		if(vmpi::my_rank==0) dw_pos.close();
	}

} //end of namespace program
