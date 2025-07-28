//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <fstream>
#include <sstream>

// Vampire headers
#include "create.hpp"
#include "material.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

// Internal create header
#include "internal.hpp"

namespace create{
namespace internal{

struct seed_point_t{
	double x; // x-position
	double y; // y-position
	double r; // point radius
};

// Function forward declaration
std::vector < seed_point_t > generate_random_seed_points(double sizex, double sizey, double scale, double randomness, double radius);
std::vector < std::vector <float> > generate_host_alloy_distribution(std::vector < seed_point_t >& seeds, const int x_cells, const int y_cells, const double smoothness, const double resolution);

//-----------------------------------------------------------------------------
//
// Function to generate a random or partially random alloy
//
// (c) R F L Evans 2016. All rights reserved.
//
//-----------------------------------------------------------------------------
void alloy(std::vector<cs::catom_t> & catom_array){

	// Alloy properties not guaranteed to exist
	// return here if unused to avoid segmentation fault
	if(create::internal::mp.size() == 0) return;

	// Print informative message to screen
	zlog << zTs() << "Calculating alloy properties of system" << std::endl;

	// Constants for distribution calculation
	const int num_alloy_materials = mp::num_materials; // local constant for number of materials
	double resolution = 2.5; // spatial reolution of concentration map (Angstroms)
	if(local_alloy) {
		resolution = local_alloy_radius;
		std::cout << resolution << ", " << int(cs::system_dimensions[0]/resolution)+1 << ", "<< int(cs::system_dimensions[1]/resolution)+1 << std::endl;
	}
	const double sizex = cs::system_dimensions[0];
	const double sizey = cs::system_dimensions[1];
	const int xcells = int(sizex/resolution)+1;
	const int ycells = int(sizey/resolution)+1;

	// array for alloy distributions for each material
	std::vector < std::vector < std::vector <float> > > distributions(num_alloy_materials);

   // generate host alloy distributions for all relevant materials
	for(int hm=0; hm<num_alloy_materials; ++hm){
		if(create::internal::mp[hm].alloy_master && (create::internal::mp[hm].host_alloy_distribution == random || create::internal::mp[hm].host_alloy_distribution == granular)){

			const double smoothness = create::internal::mp[hm].host_alloy_smoothness;
			const double scale = create::internal::mp[hm].host_alloy_scale;
			const double randomness = 0.2; // size distribution of seed points
			double radius = create::internal::mp[hm].alloy_radius;
			if(local_alloy) radius = create::internal::local_alloy_radius;
			// get seed points
			std::vector < seed_point_t > seed_points = generate_random_seed_points(sizex, sizey, scale, randomness, radius);
			// calculate interpolated probability distribution
			distributions[hm] = generate_host_alloy_distribution(seed_points, xcells, ycells, smoothness, resolution);
		}
	}

	
	// save distributions to file if required
	if(vmpi::my_rank == 0){
		for(int hm=0; hm<num_alloy_materials; ++hm){
			if(create::internal::mp[hm].save_host_alloy_profile){

				// determine filename
				std::string filename = "";
				std::stringstream fname_stream;
				// check for blank name
				if(create::internal::mp[hm].save_file_name == filename){
					fname_stream << "alloy-distribution-" << hm+1 << ".txt";
				  	filename = fname_stream.str();
				}
				else{
					filename = create::internal::mp[hm].save_file_name;
				}

				// print informative message to log file
				zlog << zTs() << "Saving alloy distribution for material " << hm+1 << " to " << filename << std::endl;

				// open file (with annnoying c_str format for pre C++11 compatibility. Maybe we can change this in 2020...)
				std::ofstream ofile;
				ofile.open(filename.c_str());

				// write data to file in gnuplot format
				for(int i=0; i<xcells; ++i){
					for(int j=0; j<ycells; ++j){
						ofile << double(i)*resolution << "\t" << double(j)*resolution << "\t" << distributions[hm][i][j] << std::endl;
					}
					ofile << std::endl;
				}
			}
		}
	}

	// Wait for all processors just in case anyone else times out
	vmpi::barrier();

	// re-seed random number generator on each CPU
	create::internal::grnd.seed(vmpi::parallel_rng_seed(create::internal::alloy_seed));

	// determine local probability
   for(unsigned int atom=0;atom<catom_array.size();atom++){

		// determine material of atom
      const int host_material=catom_array[atom].material;

		// if atom material is alloy master then reassign according to % chance
      if(create::internal::mp[host_material].alloy_master==true){

         //loop over all potential alloy materials for host
         for(int sm=0;sm<num_alloy_materials; sm++){
				if(create::internal::mp[host_material].slave_material[sm].fraction >= 0.0 || local_alloy) {
	            int slave_material = sm;
	            const double fraction = create::internal::mp[host_material].slave_material[slave_material].fraction;

					// if distribution is homogenoues calculate direct probability
					if(create::internal::mp[host_material].host_alloy_distribution == homogeneous){
						if(create::internal::grnd() < fraction) catom_array[atom].material=slave_material;
					}
					// otherwise determine probability from distribution
					else{
						double probability;
						const int i = catom_array[atom].x/resolution;
						const int j = catom_array[atom].y/resolution;
						const double variance = create::internal::mp[host_material].slave_material[slave_material].variance;
						//std::cout << i << ", " << j << ", " << distributions[host_material][i][j] << std::endl;
						switch(create::internal::mp[host_material].slave_material[slave_material].slave_alloy_distribution){
							
							case native:
								probability = fraction + variance * distributions[host_material][i][j];
								break;
						   case reciprocal:
								probability = fraction - variance * distributions[host_material][i][j];
						      break;
						   case uniform:
								probability = fraction;
						      break;
						   default: // native
								probability = fraction + variance * distributions[host_material][i][j];
						      break;
						}
							//if(probability > 0.0) std::cout << host_material << ", " << catom_array[atom].x << ", " << catom_array[atom].y << ", " << catom_array[atom].z << ", " << distributions[host_material][i][j] << \
															", " << probability << ", " << i << ", " << j << std::endl;
						// check if atom is to be replaced
						if(create::internal::grnd() < probability) catom_array[atom].material=slave_material;

					}
				} // if sm
			}
      }
   }

	return;

}

//--------------------------------------------------------------------------------------------------------
// Function to generate a list of seed points and radii within range 0-size with a characteristic scale
//
// Generate points within box extending 20% beyond sytem dimensions
//
//                  |----------------------------| 1.2
//                  |       .             .      |
//                  |  .  |----------------|  .  | 1
//                  |     |  .    .        |     |
//                  |     |        .     . |     |
//                  |     |----------------|.    | 0
//                  |    .     .         .       |
//                  |----------------------------| -0.2
//                -0.2    0                1    1.2
//
//	Note for MPI: custom generator is seeded identically on each processor, generating the same points
//
//--------------------------------------------------------------------------------------------------------
std::vector < seed_point_t > generate_random_seed_points(double size_x, double size_y, double scale, double randomness, double radius){

	if(scale == 0) {
		// empty vector to store non-touching seed points
		std::vector< seed_point_t > seeds(0);
		// determine number of trial points
		double resolution = 2.5;
		const int xcells = int(size_x/resolution)+1;
		const int ycells = int(size_y/resolution)+1;
		int num_points = xcells*ycells;
		double grain_randomness = 0.75;
		// re-seed generator on all processors
		create::internal::grnd.seed(create::internal::grain_seed);
		// generate random x,y,r trial point
		seed_point_t grain;
		grain.x = radius*4;
		grain.y = radius*4;
		grain.r = (1.0+randomness*(2.0*create::internal::grnd()-1.0))*radius;
		seeds.push_back(grain);
		for(int i=1; i<num_points*2; i++){

			// generate random x,y,r trial point
			seed_point_t grain;
			grain.x = (2*create::internal::grnd()-1.0)*resolution*grain_randomness*10.0 + (floor(i / xcells))*resolution;
			grain.y = (2*create::internal::grnd()-1.0)*resolution*grain_randomness*10.0 + (i % ycells)*resolution;
			grain.r = (1.0+randomness*(2.0*create::internal::grnd()-1.0))*radius;

			// flag to see if grains are touching
			bool touching=false;
			double minr = radius*4.5+radius*4.5*(1.0+grain_randomness*(2.0*create::internal::grnd()-1.0));
			// loop over all previous grains and check if point is not touching other grains
			for(unsigned int g=0;g<seeds.size(); g++){
				double dx = grain.x-seeds[g].x;
				double dy = grain.y-seeds[g].y;
				double rij = sqrt(dx*dx + dy*dy);
				// double minr = grain.r+seeds[g].r;
				if(rij<minr){
					touching = true;
					break;
				}
			}

			// save non-touching grains
			if(touching == false) seeds.push_back(grain);
			
		}
		return seeds;		
	}

	// empty vector to store non-touching seed points
	std::vector< seed_point_t > seeds(0);

	// determine number of trial points
	int num_points = int(25.0*size_x*size_y/(scale*scale));
	if(local_alloy) {
		num_points = 1;
		std::cout << "local alloy trials " << std::endl;
	}
	// re-seed generator on all processors
	create::internal::grnd.seed(create::internal::grain_seed);

	for(int i=0; i<num_points; i++){

		// generate random x,y,r trial point
		seed_point_t grain;
	   grain.x = (create::internal::grnd()*1.4-0.2)*size_x;
	   grain.y = (create::internal::grnd()*1.4-0.2)*size_y;
	   grain.r = (1.0+randomness*(2.0*create::internal::grnd()-1.0))*radius;
		if(local_alloy) {
			grain.x = local_alloy_x*size_x;
			grain.y = local_alloy_y*size_y;
			//grain.r = local_alloy_radius;
			std::cout << "local alloy: " << grain.x << ", " << grain.y << ", " << grain.r << std::endl;
		}
		// flag to see if grains are touching
		bool touching=false;
		// double minr = 400.0;//radius*4.0;
		// loop over all previous grains and check if point is not touching other grains
		for(unsigned int g=0;g<seeds.size(); g++){
			double dx = grain.x-seeds[g].x;
			double dy = grain.y-seeds[g].y;
			double rij = sqrt(dx*dx + dy*dy);
			double minr = grain.r+seeds[g].r;
			if(rij<minr){
				touching = true;
				break;
			}
		}

		// save non-touching grains
		if(touching == false) seeds.push_back(grain);

	}

	return seeds;

}

//--------------------------------------------------------------------------------------------------------
// Function to generate a probability distribution map using Gaussian fitting to an array of seed points
//
//--------------------------------------------------------------------------------------------------------
std::vector < std::vector <float> > generate_host_alloy_distribution(std::vector < seed_point_t >& seeds, const int x_cells, const int y_cells, const double smoothness, const double resolution){

	// resize array for stroing distribution
	std::vector < std::vector <float> > distribution;
	distribution.resize(x_cells);
	for(int i=0;i<x_cells; ++i) distribution[i].resize(y_cells);

	// running totals for normalization
	double max = -100000.0;
	double min = 100000.0;
	double sum = 0.0;

	// Calculate density at each cell
	for(int i=0; i<x_cells; ++i){
		for(int j=0; j<y_cells; ++j){
			// save cell position
			const double x = double(i)*resolution;
			const double y = double(j)*resolution;


         double density = 0.0;
			// loop over all seed points and calculate cumulative density
			for(unsigned int g=0;g<seeds.size(); g++){
            double gx = seeds[g].x;
            double gy = seeds[g].y;
            double gr = seeds[g].r;
            double dx = x-gx;
            double dy = y-gy;
            double rij2 = dx*dx + dy*dy;
            double factor = exp(-rij2/(smoothness*gr*gr));
			if(rij2 > gr*gr) factor = 0.0;
			else factor = 1.0;
            density +=factor;
         }

			// calculate new minima/maxima
			if(density > max) max = density;
			if(density < min) min = density;

			// add to running total for average calculation
			sum += density;

			// save total density to array
			distribution[i][j] = float(density);

		}
	}

	// Float constants
	const float mean = float(sum/(double(x_cells)*double(y_cells)));
	const float minf = float(min)-mean;
	const float maxf = float(max)-mean;

	// calculate largest scaling factor
	float lmax = 0.0;
	if(fabs(maxf) > lmax) lmax = fabs(maxf);
	if(fabs(minf) > lmax) lmax = fabs(minf);
	const float ifactor = 1.0f/lmax;

	// normalise distribution between -1 and +1
	for(int i=0; i<x_cells; ++i){
		for(int j=0; j<y_cells; ++j){
			distribution[i][j] = (distribution[i][j]-mean)*ifactor;
		}
	}

	return distribution;

}

} // end of internal namespace
} // end of create namespace
