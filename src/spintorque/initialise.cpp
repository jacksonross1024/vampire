//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and P Chureemart 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <complex>
#include <iostream>
#include <vector>

// Vampire headers
#include "errors.hpp"
#include "spintorque.hpp"
#include "vio.hpp"
#include "atoms.hpp"

// Spin Torque headers
#include "internal.hpp"

//-----------------------------------------------------------------------------
// Forward function declarations
//-----------------------------------------------------------------------------
namespace st{
   namespace internal{
      void set_microcell_properties(const std::vector<int>& atom_type_array, const int num_local_atoms);
   }
}

namespace st{

//-----------------------------------------------------------------------------
// Function for initialising spin torque data structure and variables
//-----------------------------------------------------------------------------
void initialise(const double system_dimensions_x,
                const double system_dimensions_y,
                const double system_dimensions_z,
                const std::vector<double>& atom_coords_x,
                const std::vector<double>& atom_coords_y,
                const std::vector<double>& atom_coords_z,
                const std::vector<int>& atom_type_array,
                const int num_local_atoms){

   //-------------------------------------------------------------------------------------
   // Check for spin torque calculation enabled, if not do nothing
   //-------------------------------------------------------------------------------------
   //if(!st::internal::enabled) return;
    if(st::internal::enabled==false) return;

   // output informative message
   zlog << zTs() << "Initialising data structures for spin torque calculation." << std::endl;

   //-------------------------------------------------------------------------------------
   // Determine transformation between x,y,z in actual and spin torque coordinate system
   //-------------------------------------------------------------------------------------
  
   if(st::internal::current_direction==0){
      st::internal::stx=2; // c[stx] = c[2] = atom_x // current direction
      st::internal::sty=1; // c[sty] = c[1] = atom_y
      st::internal::stz=0; // c[stz] = c[0] = atom_z
   }
   else if(st::internal::current_direction==1){
      st::internal::stx=0;// c[stx] = c[0] = atom_x
      st::internal::sty=2;// c[sty] = c[2] = atom_y // current direction
      st::internal::stz=1;// c[stz] = c[1] = atom_z
   }
   // st::internal::micro_cell_thickness = st::internal::micro_cell_size[stx];
   //-------------------------------------------------------------------------------------
   // Calculate number of stacks and microcells
   //-------------------------------------------------------------------------------------

   // Make a small array for system dimensions
   double system_dimensions[3]={system_dimensions_x,system_dimensions_y,system_dimensions_z};
   std::cout << "system dimensions " << system_dimensions[0] << ", " << system_dimensions[1] << ", " << system_dimensions[2] << std::endl;
   // determine number of cells in each stack  (global)
   st::internal::num_microcells_per_stack = 1+ceil((system_dimensions[st::internal::stz]+0.01)/st::internal::micro_cell_thickness);
   std::cout << "microcells per stack " << st::internal::num_microcells_per_stack << ", microcell volume (A^3): " <<  st::internal::micro_cell_size[st::internal::stx]*st::internal::micro_cell_size[st::internal::sty]*st::internal::micro_cell_thickness << std::endl;
   // determine number of stacks in x and y (global)
   st::internal::num_x_stacks = ceil((system_dimensions[st::internal::stx]+0.0)/st::internal::micro_cell_size[st::internal::stx]);
   std::cout << "x_stacks " << st::internal::num_x_stacks << ", " << system_dimensions[st::internal::stx]+0.0 << "/" << st::internal::micro_cell_size[st::internal::stx] << std::endl;
   st::internal::num_y_stacks = ceil((system_dimensions[st::internal::sty]+0.0)/st::internal::micro_cell_size[st::internal::sty]);
   std::cout << "y_stacks " << st::internal::num_y_stacks << ", " << system_dimensions[st::internal::sty]+0.0 << "/" << st::internal::micro_cell_size[st::internal::sty] << std::endl;
   // determine total number of stacks
   st::internal::num_stacks = st::internal::num_x_stacks*st::internal::num_y_stacks;
   std::cout << st::internal::num_stacks << std::endl;
   // allocate array to store index of first element of stack
   st::internal::stack_index.resize(st::internal::num_stacks);

   //-------------------------------------------------------------------------------------
   // allocate microcell data
   //-------------------------------------------------------------------------------------
   const int array_size = st::internal::num_stacks*st::internal::num_microcells_per_stack;

   st::internal::beta_cond.resize(array_size, 0.0); /// spin polarisation (conductivity)
   st::internal::beta_diff.resize(array_size, 0.0); /// spin polarisation (diffusion)
   st::internal::sa_infinity.resize(array_size, 0.0); /// intrinsic spin accumulation
   st::internal::lambda_sdl.resize(array_size, 0.0); /// spin diffusion length
   st::internal::diffusion.resize(array_size, 0.0); /// diffusion constant Do
   st::internal::sd_exchange.resize(array_size, 0.0); /// diffusion constant Do
   st::internal::a.resize(array_size, 0.0); // a parameter for spin accumulation
   st::internal::b.resize(array_size, 0.0); // b parameter for spin accumulation
   st::internal::coeff_ast.resize(array_size, 0.0);
   st::internal::coeff_nast.resize(array_size, 0.0);
   st::internal::cell_natom.resize(array_size, 0.0);


   const int three_vec_array_size = 3*array_size;
   st::internal::pos.resize(three_vec_array_size); /// microcell position
   st::internal::m.resize(three_vec_array_size); // magnetisation
   st::internal::j.resize(three_vec_array_size); // spin current
   st::internal::sa.resize(three_vec_array_size, 0.0); // spin accumulation
   st::internal::sa_sot.resize(three_vec_array_size); // spin accumulation
   st::internal::spin_torque.resize(three_vec_array_size); // spin torque
   st::internal::ast.resize(three_vec_array_size); // adiabatic spin torque
   st::internal::nast.resize(three_vec_array_size); // non-adiabatic spin torque
   st::internal::total_ST.resize(three_vec_array_size); // non-adiabatic spin torque
   st::internal::spin_acc_sign.resize(three_vec_array_size, 0.0);

   st::internal::sa_sum.resize(three_vec_array_size, 0.0);
         st::internal::sa_sot_sum.resize(three_vec_array_size, 0.0);
        // double m_sum[size] = {0.0};
         st::internal::j_sum.resize(three_vec_array_size, 0.0);
         st::internal::coeff_ast_sum.resize(three_vec_array_size, 0.0);
         st::internal::coeff_nast_sum.resize(three_vec_array_size, 0.0);
         st::internal::ast_sum.resize(three_vec_array_size, 0.0);
         st::internal::nast_sum.resize(three_vec_array_size, 0.0);
         st::internal::total_ST_sum.resize(three_vec_array_size, 0.0);
         st::internal::cell_natom_sum.resize(array_size, 0);

   //---------------------------------------------------
   // Noi Initialise j,sa, st, ast, nast here?
   //---------------------------------------------------
   for(int cell = 0; cell < array_size; ++cell){
      st::internal::sa[3*cell+0] = 0.0;
      st::internal::sa[3*cell+1] = 0.0;
      st::internal::sa[3*cell+2] = 0.0;
      st::internal::sa_sot[3*cell+0] = 0.0;
      st::internal::sa_sot[3*cell+1] = 0.0;
      st::internal::sa_sot[3*cell+2] = 0.0;
      st::internal::j [3*cell+0] = 0.0;
      st::internal::j [3*cell+1] = 0.0;
      st::internal::j [3*cell+2] = 0.0;
   }

   //---------------------------------------------------
   // Determine which atoms belong to which stacks
   //---------------------------------------------------
   {
   int ncx = st::internal::num_x_stacks; // temporary variables for readability
   int ncy = st::internal::num_y_stacks;
   int ncz = st::internal::num_microcells_per_stack;

   // Set cell and stack counters
   int cell=0;
   int stack=0;
      st::internal::cell_stack_index.resize(ncx*ncy*ncz);
   // Allocate space for 3D supercell array (ST coordinate system)
   std::vector<std::vector<std::vector<int> > > supercell_array;
   supercell_array.resize(ncx);

   for(int i=0;i<ncx;++i){
      supercell_array[i].resize(ncy);
      for(int j=0;j<ncy;++j){
         // std::cout << stack << std::endl;
         supercell_array[i][j].resize(ncz);
         // set starting cell for each stack
         st::internal::stack_index[stack]=cell;

         // increment stack counter
         stack++;
         // store cell coordinates
         for(int k=0; k<ncz; ++k){
            // associate cell with position i,j,k
            supercell_array[i][j][k]=cell;
            // save ijk coordinates as microcell positions
            st::internal::pos.at(3*cell+0)=i;
            st::internal::pos.at(3*cell+1)=j;
            st::internal::pos.at(3*cell+2)=k;
            st::internal::cell_stack_index[cell] = stack;
            // increment cell number
            cell++;
         }
      }
   }

   // define array to store atom-microcell associations
   st::internal::atom_st_index.resize(num_local_atoms);
 
   // Determine number of cells in x,y,z (ST coordinate system)
   const int d[3]={ncx,ncy,ncz};
   const double cs[3] = {st::internal::micro_cell_size[st::internal::stx], st::internal::micro_cell_size[st::internal::sty], st::internal::micro_cell_thickness}; // cell size

   // Assign atoms to cells
   for(int atom=0;atom<num_local_atoms;atom++){
      // if(st::internal::remove_nm) {
      //    if(atoms::type_array[atom] == st::internal::remove_nm) continue;
      // }
      // temporary for atom coordinates
      double c[3];
      // convert atom coordinates to st reference frame
      c[st::internal::stx]=atom_coords_x[atom]+0.0001;
      c[st::internal::sty]=atom_coords_y[atom]+0.0001;
      c[st::internal::stz]=atom_coords_z[atom]+0.0001;
      int scc[3]={0,0,0}; // super cell coordinates
      for(int i=0;i<3;i++){
         // Determine supercell coordinates for atom (rounding down)
         scc[i]=int(floor(c[i]/cs[i]));
       //  std::cout << i << ", " << scc[i] << ", " << c[i] << ", " << cs[i] << std::endl;
         // Always check cell in range
         if(scc[i]<0 || scc[i]>= d[i]){
            terminaltextcolor(RED);
            std::cerr << "Error - atom out of supercell range in spin torque microcell calculation!" << std::endl;
            terminaltextcolor(WHITE);
            #ifdef MPICF
            terminaltextcolor(RED);
            std::cerr << "\tCPU Rank: " << vmpi::my_rank << std::endl;
            terminaltextcolor(WHITE);
            #endif
            terminaltextcolor(RED);
            std::cerr << "\tAtom number:      " << atom << std::endl;
            std::cerr << "\tAtom coordinates: " << c[0] << "\t" << c[1] << "\t" << c[2] << "\t" << std::endl;
            std::cerr << "\tReal coordinates: " << atom_coords_x[atom] << "\t" << atom_coords_y[atom] << "\t" << atom_coords_z[atom] << "\t" << std::endl;
            std::cerr << "\tCell coordinates: " << scc[0] << "\t" << scc[1] << "\t" << scc[2] << "\t" << std::endl;
            std::cerr << "\tCell maxima:      " << d[0] << "\t" << d[1] << "\t" << d[2] << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
         }
      }
      // If no error for range then assign atom to cell.
      st::internal::atom_st_index[atom]= supercell_array[scc[0]][scc[1]][scc[2]+1]; // move cells up by one in z
      
   }
   } // end of supercell assignment of atoms

   //-------------------------------------------------------
   // Determine microcell properties from atomic properties
   //-------------------------------------------------------
   st::internal::set_microcell_properties(atom_type_array, num_local_atoms);

   //-------------------------------------------------------
   // Save value of local num atoms and resize field arrays
   //-------------------------------------------------------
   st::internal::num_local_atoms = num_local_atoms;
   st::internal::x_field_array.resize(num_local_atoms); // arrays to store atomic spin torque field
   st::internal::y_field_array.resize(num_local_atoms);
   st::internal::z_field_array.resize(num_local_atoms);

   // optionally output base microcell data
   st::internal::output_base_microcell_data();

       #ifdef MPICF 
     
   //determine mpi core stack list 
      int removed_stacks = st::internal::num_stacks;
      std::vector<int> mpi_stack_id(removed_stacks, 0);
   if(vmpi::num_processors > removed_stacks) {
      std::cout << "mpirun threads requested larger than spin-torque decomposition allows" << std::endl;
      err::vexit();
   }
   int residual = removed_stacks % vmpi::num_processors;
  // st::internal::mpi_stack_list.resize(int(floor(st::internal::num_stacks/vmpi::num_processors))+1);
   for(int s = 0; s < int(floor(removed_stacks/vmpi::num_processors)); s++) {
      st::internal::mpi_stack_list.push_back(s*vmpi::num_processors + vmpi::my_rank);
      mpi_stack_id.at(s*vmpi::num_processors + vmpi::my_rank) = 1;
   } 
   if(residual > 0 && vmpi::my_rank < residual){
      st::internal::mpi_stack_list.push_back(removed_stacks - vmpi::my_rank -1);
      mpi_stack_id.at(removed_stacks - vmpi::my_rank -1) = 1;
   }


   // if(err::check) {
      //check stt stack mpi setup 
      // int stacks[st::internal::num_stacks] = {0};
      // int mpi_stacks[st::internal::num_stacks] = {-1};
      int stack_sum = 0;
      int size = st::internal::mpi_stack_list.size();
      MPI_Reduce(&size,&stack_sum, 1,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &mpi_stack_id[0],  mpi_stack_id.size(),   MPI_INT,MPI_SUM, MPI_COMM_WORLD);
      bool error = false;
      if(vmpi::my_rank == 0 && stack_sum != removed_stacks) {
         std::cout << stack_sum << " != " << removed_stacks << std::endl;
         error = true;
      }
      for(int i = 0; i < mpi_stack_id.size(); i++) {
         if(mpi_stack_id[i] != 1) error = true;
      }
      MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);
      if(error) {

         for(int i = 0; i < mpi_stack_id.size(); i++) std::cout <<i << ", " << mpi_stack_id[i] << std::endl; 
          
         MPI_Barrier(MPI_COMM_WORLD);
         err::vexit();
      }
   // }
      std::cout << "MPI stack count " << removed_stacks << std::endl;
    #endif 
   return;
}

namespace internal{

   //--------------------------------------------------------------------------------
   // Function to determine spin torque properties from atomic material properties
   //--------------------------------------------------------------------------------
   void set_microcell_properties(const std::vector<int>& atom_type_array, const int num_local_atoms){
      //-------------------------------------------------------
      // Determine microcell properties from atomic properties
      //-------------------------------------------------------
      st::internal::default_properties.beta_cond = 0.0;//  0.23;
      st::internal::default_properties.beta_diff = 0.0;//0.56;
      st::internal::default_properties.sa_infinity = 0.0;// 1.48e7;
      st::internal::default_properties.lambda_sdl = 0.0;//16.0e-10; // m
      st::internal::default_properties.diffusion = 0.0;// 0.001; //m^2/s ? 
      st::internal::default_properties.sd_exchange = 0.0;//8.010883e-21; //Joule
      
     
      // Temporary array to hold number of atoms in each cell for averaging
      std::vector<double> count(st::internal::beta_cond.size(),0.0);

      // loop over all atoms
      for(int atom=0;atom<num_local_atoms;atom++) {

         // get material type
         int mat = atom_type_array[atom];

         // get microcell id
         int id = st::internal::atom_st_index[atom];
         // bool local = false;
         // #ifdef MPICH
         //    for(int i = 0; i < mpi_stack_list.size(); i++){
         //       if(mpi_stack_list[i] == )
         //    }
         //  std::cout << atoms::sublayer_array[atom] << ", " << st::internal::cell_stack_index[st::internal::atom_st_index[atom]] << std::endl;
         // determine atomic properties
         double beta_cond = st::internal::mp.at(mat).beta_cond; // beta
         double beta_diff = st::internal::mp.at(mat).beta_diff; // beta_prime
         double sa_infinity = st::internal::mp.at(mat).sa_infinity;
         double lambda_sdl = st::internal::mp.at(mat).lambda_sdl;
         double diffusion = st::internal::mp.at(mat).diffusion;
         double sd_exchange = st::internal::mp.at(mat).sd_exchange;

         //add atomic properties to microcells
         st::internal::beta_cond.at(id) += beta_cond;
         st::internal::beta_diff.at(id) += beta_diff;
         st::internal::sa_infinity.at(id) += sa_infinity;
         st::internal::lambda_sdl.at(id) += lambda_sdl;
         st::internal::diffusion.at(id) += diffusion;
         st::internal::sd_exchange.at(id) += sd_exchange;
         st::internal::spin_acc_sign.at(id) += (mat == 0) ? 0:((mat == 1) ? -1.0:1.0);
         count.at(id) += (mat == 0) ? 0:1;
      }

      // reduce microcell properties on all CPUs
      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::beta_cond[0],   st::internal::beta_cond.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::beta_diff[0],   st::internal::beta_diff.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::sa_infinity[0], st::internal::sa_infinity.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::lambda_sdl[0],  st::internal::lambda_sdl.size(),  MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::diffusion[0],   st::internal::diffusion.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::sd_exchange[0], st::internal::sd_exchange.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::spin_acc_sign[0], st::internal::spin_acc_sign.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &count[0],                     count.size(),                     MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
      #endif

      // Calculate average (mean) spin torque parameters
      for(size_t cell=0; cell<beta_cond.size(); ++cell){
         // const double nat = vmpi::num_processors*count.at(cell);
         const double nat = count.at(cell);
          st::internal::cell_natom[cell] = count.at(cell);
        
         // check for zero atoms in cell
         if(nat>0.0001){
            //  if(nat != 2)   std::cout << nat << std::endl;
            st::internal::beta_cond.at(cell)   /= nat;
            st::internal::beta_diff.at(cell)   /= nat;
            st::internal::sa_infinity.at(cell) /= nat;
            st::internal::lambda_sdl.at(cell)  /= nat;
            st::internal::diffusion.at(cell)   /= nat;
            st::internal::sd_exchange.at(cell) /= nat;
            st::internal::spin_acc_sign.at(cell) = std::round(st::internal::spin_acc_sign.at(cell)/count.at(cell));
            st::internal::default_properties.sa_infinity =  st::internal::sa_infinity.at(cell);
            if(st::internal::spin_acc_sign.at(cell) != 1.0 && st::internal::spin_acc_sign.at(cell) != -1.0) std::cout << st::internal::spin_acc_sign.at(cell) << std::endl;
         } else{
            st::internal::beta_cond.at(cell)   = st::internal::default_properties.beta_cond;
            st::internal::beta_diff.at(cell)   = st::internal::default_properties.beta_diff;
            st::internal::sa_infinity.at(cell) = st::internal::default_properties.sa_infinity;
            st::internal::lambda_sdl.at(cell)  = st::internal::default_properties.lambda_sdl;
            st::internal::diffusion.at(cell)   = st::internal::default_properties.diffusion;
            st::internal::sd_exchange.at(cell) = st::internal::default_properties.sd_exchange;
         }
      }



      // Determine a and b parameters
      const double hbar = 1.05457162e-34;
      for(size_t cell=0; cell<beta_cond.size(); ++cell){ 

         const double B  = st::internal::beta_cond[cell];
         const double Bp = st::internal::beta_diff[cell];
         const double lambda_sdl = st::internal::lambda_sdl[cell];
         const double Do = st::internal::diffusion[cell];
         const double Jsd = st::internal::sd_exchange[cell];

         const double BBp = 1.0/sqrt(1.0-B*Bp);
         const double lambda_sf = lambda_sdl*BBp;
         const double lambda_j = sqrt(2.0*hbar*Do/Jsd); // Angstroms
         const double lambda_sf2 = lambda_sf*lambda_sf;
         const double lambda_j2 = lambda_j*lambda_j;

         std::complex<double> inside (1.0/lambda_sf2, -1.0/lambda_j2);
         std::complex<double> inv_lplus = sqrt(inside);

         st::internal::a[cell] =  real(inv_lplus);
         st::internal::b[cell] = -imag(inv_lplus);

      }

      return;
   }
}

} // end of st namespace
