//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard Evans 2018. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <cmath>
#include <iostream>
#include <unordered_map>

// Vampire headers
#include "create.hpp"              // cs::pbc, cs::local_num_unit_cells, cs::total_num_unit_cells
#include "create_atoms_class.hpp"   // class definition for atoms in create module
#include "errors.hpp"
#include "neighbours.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

//-----------------------------------
// Fix for horrible windows compiler
//-----------------------------------
/*#ifdef WIN_COMPILE
   #define NOMINMAX
   #undef max
   #undef min
#endif*/

namespace neighbours{

//----------------------------------------------------------------------------------
// @brief Generate atomic neighbourlist for a generalised exchange template
//
// Assigns atoms to unit cells and then calculates all interactions between cells
//
// Partial cells can exist so ensure enough cells are generated
//
//    4    5    6    7    8
//    | ...|....|....|.   |
//
//  In this example offset=4, and max_cell = 8. Therefore 4 cells are needed.
//
//----------------------------------------------------------------------------------
void list_t::generate( std::vector<cs::catom_t>& atom_array,    // array of atoms (as reference for speed)
               unitcell::exchange_template_t& exchange, // exchange template to calculate neighbour list
               const uint64_t num_atoms_in_unit_cell,        // number of atoms in each cell to estimate interaction numbers
               const double ucdx,                       // unit cell size
               const double ucdy,
               const double ucdz
             ){

	// put number of atoms into temporary variable
	int64_t num_atoms = atom_array.size();

	// Reserve space for num_atoms
	list.reserve(num_atoms);

	// estimate number of interactions per atom
	const int64_t max_nn = int64_t( 1.1*( double(exchange.interaction_count) / double(num_atoms_in_unit_cell) ) );

	// Reserve space for each atom in neighbour list according to material type
	for(int64_t atom=0; atom < num_atoms; atom++){
		list.push_back(std::vector<neighbour_t>());
		//list[atom].reserve(max_nn);
	}

   // Calculate system dimensions and number of supercells
   const int64_t max_val=1000000000000;
   int64_t min[3] = {max_val,max_val,max_val}; // lowest cell id
   int64_t max[3] = {0,0,0}; // highest cell id
   int64_t offset[3];
   int64_t d[3];

   #ifdef MPICF
   // Use decomposition local cell range + halo. With PBC, halo atoms have
   // scx/scy/scz with offsets ±total_num_unit_cells, so min/max over atoms
   // would span the global grid and num_cells would explode (OOM).
   // Halo extent: when shared_memory_ucf use spatial cutoff (Å); otherwise use template interaction_range (cells).
   if(vmpi::mpi_mode==0){
      int64_t halo_cells[3];
      if(vmpi::shared_memory_ucf){
         const double cut = vmpi::shared_memory_halo_cutoff_angstrom;
         halo_cells[0] = static_cast<int64_t>(std::ceil(cut / ucdx));
         halo_cells[1] = static_cast<int64_t>(std::ceil(cut / ucdy));
         halo_cells[2] = static_cast<int64_t>(std::ceil(cut / ucdz));
      }
      else{
         const int64_t hr = static_cast<int64_t>(cs::unit_cell.interaction_range);
         halo_cells[0] = halo_cells[1] = halo_cells[2] = hr;
      }
      const int64_t min_cell[3] = {
         static_cast<int64_t>(vmpi::min_dimensions[0] / ucdx),
         static_cast<int64_t>(vmpi::min_dimensions[1] / ucdy),
         static_cast<int64_t>(vmpi::min_dimensions[2] / ucdz)
      };
      offset[0] = min_cell[0] - halo_cells[0];
      offset[1] = min_cell[1] - halo_cells[1];
      offset[2] = min_cell[2] - halo_cells[2];
      d[0] = static_cast<int64_t>(cs::local_num_unit_cells[0]) + 2*halo_cells[0];
      d[1] = static_cast<int64_t>(cs::local_num_unit_cells[1]) + 2*halo_cells[1];
      d[2] = static_cast<int64_t>(cs::local_num_unit_cells[2]) + 2*halo_cells[2];
   }
   else
   #endif
   {
      // find supercell range of atoms on this CPU (serial or replicated MPI)
      for(int64_t atom = 0; atom < num_atoms; atom++){
         int64_t c[3] = { atom_array[atom].scx, atom_array[atom].scy, atom_array[atom].scz };
         for(int i = 0; i < 3; i++){
            if( c[i] < min[i] ) min[i] = c[i];
            if( c[i] > max[i] ) max[i] = c[i];
         }
      }
      for(int i = 0; i < 3; i++){
         if(min[i] > max_val){
            std::cerr << "Programmer error! too many supercells in atom list" << std::endl;
         }
      }
      offset[0] = min[0]; offset[1] = min[1]; offset[2] = min[2];
      d[0] = max[0] - offset[0] + 1;
      d[1] = max[1] - offset[1] + 1;
      d[2] = max[2] - offset[2] + 1;
   }

	// Sparse supercell: one map per cell (uc_id -> atom index). Memory O(local atoms)
	// instead of O(local_cells * num_atoms_in_unit_cell), which caused OOM for large UCFs.
	const uint64_t num_cells = d[0]*d[1]*d[2];
	std::vector<std::unordered_map<uint64_t, int64_t> > supercell_map(num_cells);
	int64_t num_neighbours = 0;
	uint64_t total_neighbours = 0;

   #ifdef MPICF
      if(vmpi::master){
         zlog << zTs() << "Neighbourlist supercell uses sparse storage (O(local atoms) per rank)" << std::endl;
      }
   #else
      zlog << zTs() << "Neighbourlist supercell uses sparse storage" << std::endl;
   #endif

   zlog << zTs() << "Allocating memory for supercell (sparse) in neighbourlist calculation" << std::endl;
   zlog << zTs() << "\tAllocating memory done"<< std::endl;

	// declare 1D cell array to loop over
	std::vector< std::vector <int> > cell_coord_array;
	cell_coord_array.reserve(num_cells);
	for(uint64_t i=0; i < num_cells; i++){
		cell_coord_array.push_back(std::vector<int>());
		cell_coord_array[i].resize(3);
	}

	// Initialise cell_array
	uint64_t cell=0; // cell counter
	for(unsigned int x=0;x<d[0];x++){
		for(unsigned int y=0;y<d[1];y++){
			for(unsigned int z=0;z<d[2];z++){
				//save cell coordinates
				cell_coord_array[cell][0]=x;
				cell_coord_array[cell][1]=y;
				cell_coord_array[cell][2]=z;
				cell++;
			}
		}
	}

   // Inform user of time intensive process
   zlog << zTs() << "Populating supercell array for neighbourlist calculation..."<< std::endl;

	// Populate supercell array with atom numbers
	for(int64_t atom=0; atom < num_atoms; atom++){
		int64_t sc[3] = { atom_array[atom].scx, atom_array[atom].scy, atom_array[atom].scz };
		int64_t scc[3];
		#ifdef MPICF
		if(vmpi::mpi_mode==0){
			// With PBC, halo atoms have scx/scy/scz with ±total_num_unit_cells offset. Wrap into local grid.
			for(int i=0;i<3;i++) scc[i] = ((sc[i] - offset[i]) % d[i] + d[i]) % d[i];
		}
		else
		#endif
		{
			scc[0] = sc[0] - offset[0]; scc[1] = sc[1] - offset[1]; scc[2] = sc[2] - offset[2];
		}

		double c[3]={atom_array[atom].x, atom_array[atom].y, atom_array[atom].z};

		for(int i=0;i<3;i++){
			if( scc[i] < 0 || scc[i] >= d[i] ){
            terminaltextcolor(RED);
				std::cerr << "Error - atom out of supercell range in neighbourlist calculation!" << std::endl;
            terminaltextcolor(WHITE);
				#ifdef MPICF
				terminaltextcolor(RED);
				std::cerr << "\tCPU Rank: " << vmpi::my_rank << std::endl;
				terminaltextcolor(WHITE);
				#endif
				terminaltextcolor(RED);
				std::cerr << "\tAtom number:      " << atom << std::endl;
				std::cerr << "\tAtom coordinates: " << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << std::endl;
				std::cerr << "\tmin coordinates:  " << min[0] << "\t" << min[1] << "\t" << min[2] << "\t" << std::endl;
				std::cerr << "\tmax coordinates:  " << max[0] << "\t" << max[1] << "\t" << max[2] << "\t" << std::endl;
				std::cerr << "\tCell coordinates: " << scc[0] << "\t" << scc[1] << "\t" << scc[2] << "\t" << std::endl;
				std::cerr << "\tCell maxima:      " << d[0] << "\t" << d[1] << "\t" << d[2] << std::endl;
				std::cerr << "\tCell offset:      " << offset[0] << "\t" << offset[1] << "\t" << offset[2] << std::endl;
				std::cerr << "\tCell offest (dp): " << offset[0]*ucdx << "\t" << offset[1]*ucdy << "\t" << offset[2]*ucdz << std::endl;
				terminaltextcolor(WHITE);
				err::vexit();
			}
		}
		// Add atom to sparse supercell map (uc_id -> atom index)
		const uint64_t uc_id = atom_array[atom].uc_id;
		const uint64_t cell_linear = static_cast<uint64_t>(scc[0]) * static_cast<uint64_t>(d[1]) * static_cast<uint64_t>(d[2])
		                          + static_cast<uint64_t>(scc[1]) * static_cast<uint64_t>(d[2])
		                          + static_cast<uint64_t>(scc[2]);
		supercell_map[cell_linear][uc_id] = atom;
	}
   // Inform user of progress
   zlog << zTs() << "\tPopulating supercell array completed"<< std::endl;

   // Estimate neighbour list memory from this rank's atoms (list is filled by push_back below)
   const uint64_t est_neighbours_per_atom = (num_atoms_in_unit_cell > 0)
      ? (exchange.interaction_count / num_atoms_in_unit_cell) : 0;
   num_neighbours = static_cast<int64_t>(num_atoms) * static_cast<int64_t>(est_neighbours_per_atom);
   #ifdef MPICF
      total_neighbours = 0;
      if(vmpi::master){
         zlog << zTs() << "Memory estimate for neighbour list (per rank, from local atoms): " <<
         sizeof(neighbour_t)*num_neighbours/(1.0e6) << " MB (this rank: " << num_atoms << " atoms)" << std::endl;
      }
   #else
      zlog << zTs() << "Memory estimate for neighbour list: " <<
      sizeof(neighbour_t)*num_neighbours/1.0e6 << " MB" << std::endl;
   #endif

	// Generate neighbour list and inform user
	std::cout <<"Generating neighbour list"<< std::flush;
   zlog << zTs() << "Generating neighbour list..."<< std::endl;

   // temporary neighour to simplify memory management
	neighbour_t tmp_nt;

   // copy number of interactions to temporary constant
   const uint64_t num_interactions = exchange.interaction_count;

	// Loop over all cells
	for(uint64_t cell = 0; cell < num_cells; cell++){

      // Print out progress indicator to user
		if( cell % ( num_cells / 10 + 1 ) == 0 ){
			std::cout << "." << std::flush;
		}

      // get supercell coordinates of cell
		int scc[3]={ cell_coord_array[cell][0],
                   cell_coord_array[cell][1],
                   cell_coord_array[cell][2]};

		// Loop over all interactions in exchange template
		for(uint64_t i = 0; i < num_interactions; i++){

			const uint64_t atom = exchange.interaction_ptr[i].i;
			const uint64_t natom = exchange.interaction_ptr[i].j;

			int nx = exchange.interaction_ptr[i].dx + scc[0];
			int ny = exchange.interaction_ptr[i].dy + scc[1];
			int nz = exchange.interaction_ptr[i].dz + scc[2];

         // vector from i->j
         double vx=0.0;
         double vy=0.0;
         double vz=0.0;

         #ifdef MPICF
           // Parallel periodic boundaries are handled explicitly during the
           // halo region setup
         #else
         // Wrap around for periodic boundaries
         // Consider virtual atom position for position vector
         if(cs::pbc[0]==true){
            if(nx>=int(d[0])){
               nx=nx-d[0];
               vx=vx+d[0]*ucdx;
            }
            else if(nx<0){
               nx=nx+d[0];
               vx=vx-d[0]*ucdx;
            }
         }
         if(cs::pbc[1]==true){
            if(ny>=int(d[1])){
               ny=ny-d[1];
               vy=vy+d[1]*ucdy;
            }
            else if(ny<0){
               ny=ny+d[1];
               vy=vy-d[1]*ucdy;
            }
         }
         if(cs::pbc[2]==true){
            if(nz>=int(d[2])){
               nz=nz-d[2];
               vz=vz+d[2]*ucdz;
            }
            else if(nz<0){
               nz=nz+d[2];
               vz=vz-d[2]*ucdz;
            }
         }
         #endif
         #ifdef MPICF
         // With local grid + PBC, wrap neighbour cell into [0, d[]) for lookup
         if(vmpi::mpi_mode==0){
            nx = (nx % static_cast<int>(d[0]) + d[0]) % static_cast<int>(d[0]);
            ny = (ny % static_cast<int>(d[1]) + d[1]) % static_cast<int>(d[1]);
            nz = (nz % static_cast<int>(d[2]) + d[2]) % static_cast<int>(d[2]);
         }
         #endif
         // check for out-of-bounds access
         if( (nx >= 0 && static_cast<int64_t>(nx) < d[0] ) &&
             (ny >= 0 && static_cast<int64_t>(ny) < d[1] ) &&
             (nz >= 0 && static_cast<int64_t>(nz) < d[2] ) ){
            const uint64_t cell_i = static_cast<uint64_t>(scc[0])*static_cast<uint64_t>(d[1])*static_cast<uint64_t>(d[2])
                                  + static_cast<uint64_t>(scc[1])*static_cast<uint64_t>(d[2])
                                  + static_cast<uint64_t>(scc[2]);
            const uint64_t cell_j = static_cast<uint64_t>(nx)*static_cast<uint64_t>(d[1])*static_cast<uint64_t>(d[2])
                                  + static_cast<uint64_t>(ny)*static_cast<uint64_t>(d[2])
                                  + static_cast<uint64_t>(nz);
            auto it_i = supercell_map[cell_i].find(atom);
            auto it_j = supercell_map[cell_j].find(natom);
            if(it_i != supercell_map[cell_i].end() && it_j != supercell_map[cell_j].end()) {

               int64_t atomi = it_i->second;
               int64_t atomj = it_j->second;

               //std::cout << "int_id: " << i << "\tatom i: " << atomi << "\tatom j: " << atomj << "\tuc_i: " << atom << "\tuc_j: " << natom << std::endl;

               // Load atom positions
               double ix = atom_array[atomi].x; // Already in A
               double iy = atom_array[atomi].y;
               double iz = atom_array[atomi].z;
               double jx = atom_array[atomj].x;
               double jy = atom_array[atomj].y;
               double jz = atom_array[atomj].z;

               //std::cout << "\tpi:    " << ix << "\t" << iy << "\t" << iz << std::endl;
               //std::cout << "\tpj:    " << jx << "\t" << jy << "\t" << jz << std::endl;
               //std::cout << "\tv_uc:  " << vx << "\t" << vy << "\t" << vz << std::endl;

               vx += jx-ix;
               vy += jy-iy;
               vz += jz-iz;

               //std::cout << "\tv:     " << jx-ix << "\t" << jy-iy << "\t" << jz-iz << std::endl;
               //std::cout << "\tv_eff: " << vx << "\t" << vy << "\t" << vz << std::endl;

               // set neighbour data
               tmp_nt.nn = atomj;
               tmp_nt.i = i;
               tmp_nt.vx = vx;
               tmp_nt.vy = vy;
               tmp_nt.vz = vz;

               list[atomi].push_back(tmp_nt);

            }
			}
		}
	}

   // Inform user neighbour list calculation is complete
	if(vmpi::my_rank == 0){
		terminaltextcolor(GREEN);
		std::cout << "done!" << std::endl;
		terminaltextcolor(WHITE);
	}
   zlog << zTs() << "\tNeighbour list calculation complete"<< std::endl;

	// Deallocate sparse supercell map
   zlog << zTs() << "Deallocating supercell map for neighbour list calculation" << std::endl;
	for(uint64_t c = 0; c < num_cells; c++){
		supercell_map[c].clear();
	}
	supercell_map.clear();
   zlog << zTs() << "\tSupercell map deallocated" << std::endl;

	return;
}

} // End of namespace neighbours
