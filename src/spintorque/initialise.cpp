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
   if(st::internal::sot_sa ) {
       std::cout << "Spin Torque module set for SOT calculation." << std::endl;
       std::cout << "Make sure HM layer is set to 'keep'... " << std::endl;
   }
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

   if(st::internal::microcell_decomp_type == "A") {}
   else if (st::internal::microcell_decomp_type == "mpi") {
      st::internal::micro_cell_size[st::internal::stx] = system_dimensions[st::internal::stx] / st::internal::micro_cell_size[st::internal::stx];
      st::internal::micro_cell_size[st::internal::sty] = system_dimensions[st::internal::sty] / st::internal::micro_cell_size[st::internal::sty];
      st::internal::micro_cell_size[st::internal::stz] = system_dimensions[st::internal::stz] / st::internal::micro_cell_size[st::internal::stz];
   } else {
      terminaltextcolor(RED);
      std::cerr << "wrong microcell decomp type " << std::endl;
      terminaltextcolor(WHITE);
      err::vexit();
   }

   std::cout << "system dimensions " << system_dimensions[0] << ", " << system_dimensions[1] << ", " << system_dimensions[2] << std::endl;

   // determine number of cells in each stack  (global)
   st::internal::num_microcells_per_stack = ceil((system_dimensions[st::internal::stz]+0.0000)/st::internal::micro_cell_thickness);
   std::cout << "microcells per stack " << st::internal::num_microcells_per_stack << ", microcell volume (A^3): " <<  st::internal::micro_cell_size[st::internal::stx]*st::internal::micro_cell_size[st::internal::sty]*st::internal::micro_cell_thickness << std::endl;
   
   // determine number of stacks in x and y (global)
   st::internal::num_x_stacks = ceil((system_dimensions[st::internal::stx]+0.0000)/st::internal::micro_cell_size[st::internal::stx]);
   std::cout << "x_stacks " << st::internal::num_x_stacks << ", " << system_dimensions[st::internal::stx]+0.0001 << "/" << st::internal::micro_cell_size[st::internal::stx] << std::endl;
   st::internal::num_y_stacks = ceil((system_dimensions[st::internal::sty]+0.0000)/st::internal::micro_cell_size[st::internal::sty]);
   std::cout << "y_stacks " << st::internal::num_y_stacks << ", " << system_dimensions[st::internal::sty]+0.0001 << "/" << st::internal::micro_cell_size[st::internal::sty] << std::endl;
   
   // determine total number of stacks
   st::internal::num_stacks_y = st::internal::num_x_stacks*st::internal::num_y_stacks;
   st::internal::num_stacks_x = st::internal::num_y_stacks*st::internal::num_microcells_per_stack;
   std::cout << "Total stacks: " << st::internal::num_stacks_y << std::endl;
      if(st::internal::sot_sa) std::cout << "Total stacks x: " << st::internal::num_stacks_x << std::endl;
   
   // allocate array to store index of first element of stack
   st::internal::stack_index_y.resize(st::internal::num_stacks_y);
   st::internal::stack_index_x.resize(st::internal::num_stacks_x);
   //-------------------------------------------------------------------------------------
   // allocate microcell data
   //-------------------------------------------------------------------------------------
   const int array_size = st::internal::num_x_stacks*st::internal::num_y_stacks*st::internal::num_microcells_per_stack;

   //stt
   st::internal::beta_cond.resize(array_size, 0.0); /// spin polarisation (conductivity)
   st::internal::beta_diff.resize(array_size, 0.0); /// spin polarisation (diffusion)
   st::internal::sa_infinity.resize(array_size, 0.0); /// intrinsic spin accumulation
   st::internal::lambda_sdl.resize(array_size, 0.0); /// spin diffusion length
   st::internal::diffusion.resize(array_size, 0.0); /// diffusion constant Do
   st::internal::sd_exchange.resize(array_size, 0.0); /// diffusion constant Do
   st::internal::a.resize(array_size, 0.0); // a parameter for spin accumulation
   st::internal::b.resize(array_size, 0.0); // b parameter for spin accumulation

   //sot
   st::internal::sot_beta_cond.resize(array_size, 0.0); /// spin polarisation (conductivity)
   st::internal::sot_beta_diff.resize(array_size, 0.0); /// spin polarisation (diffusion)
   st::internal::sot_sa_infinity.resize(array_size, 0.0); /// intrinsic spin accumulation
   st::internal::sot_lambda_sdl.resize(array_size, 0.0); /// spin diffusion length
   st::internal::sot_diffusion.resize(array_size, 0.0); /// diffusion constant Do
   st::internal::sot_sd_exchange.resize(array_size, 0.0); /// diffusion constant Do
   st::internal::sot_a.resize(array_size, 0.0); // a parameter for spin accumulation
   st::internal::sot_b.resize(array_size, 0.0); // b parameter for spin accumulation
   st::internal::sot_sa_source.resize(array_size, false);
      
   st::internal::coeff_ast.resize(array_size, 0.0);
   st::internal::coeff_nast.resize(array_size, 0.0);
   st::internal::cell_natom.resize(array_size, 0.0);
   const int three_vec_array_size = 3*array_size;
   st::internal::spin_acc_sign.resize(three_vec_array_size, 0.0);

   st::internal::pos.resize(three_vec_array_size,0.0); /// microcell position
   st::internal::m.resize(three_vec_array_size,0.0); // magnetisation
   st::internal::j_final_up_y.resize(three_vec_array_size,0.0);
   st::internal::j_final_down_y.resize(three_vec_array_size,0.0);
   st::internal::j_int_up_y.resize(three_vec_array_size,0.0); // spin current
   st::internal::j_int_down_y.resize(three_vec_array_size,0.0); // spin current
   st::internal::j_init_up_y.resize(three_vec_array_size,0.0); // spin current
   st::internal::j_init_down_y.resize(three_vec_array_size,0.0); // spin current
   st::internal::j_final_up_x.resize(three_vec_array_size,0.0); // spin current
   st::internal::sa_final.resize(three_vec_array_size, 0.0); // spin accumulation
   st::internal::sa_int.resize(three_vec_array_size,0.0); // sot sa up stack
   st::internal::spin_torque.resize(three_vec_array_size,0.0); // spin torque
   st::internal::ast.resize(three_vec_array_size,0.0); // adiabatic spin torque
   st::internal::nast.resize(three_vec_array_size,0.0); // non-adiabatic spin torque
   st::internal::total_ST.resize(three_vec_array_size,0.0); // non-adiabatic spin torque

   st::internal::sa_sum.resize(three_vec_array_size, 0.0);
   st::internal::j_final_up_x_sum.resize(three_vec_array_size, 0.0);
   st::internal::j_final_up_y_sum.resize(three_vec_array_size, 0.0);
   st::internal::j_final_down_y_sum.resize(three_vec_array_size, 0.0);
   st::internal::coeff_ast_sum.resize(three_vec_array_size, 0.0);
   st::internal::coeff_nast_sum.resize(three_vec_array_size, 0.0);
   st::internal::ast_sum.resize(three_vec_array_size, 0.0);
   st::internal::nast_sum.resize(three_vec_array_size, 0.0);
   st::internal::total_ST_sum.resize(three_vec_array_size, 0.0);
   st::internal::cell_natom_sum.resize(array_size, 0);

   //---------------------------------------------------
   // Noi Initialise j,sa, st, ast, nast here?
   //---------------------------------------------------
   // for(int cell = 0; cell < array_size; ++cell){
   //    st::internal::sa[3*cell+0] = 0.0;
   //    st::internal::sa[3*cell+1] = 0.0;
   //    st::internal::sa[3*cell+2] = 0.0;
   //    st::internal::sa_sot[3*cell+0] = 0.0;
   //    st::internal::sa_sot[3*cell+1] = 0.0;
   //    st::internal::sa_sot[3*cell+2] = 0.0;
   //    st::internal::j [3*cell+0] = 0.0;
   //    st::internal::j [3*cell+1] = 0.0;
   //    st::internal::j [3*cell+2] = 0.0;
   // }

   //---------------------------------------------------
   // Determine which atoms belong to which stacks
   //---------------------------------------------------
   {
   int ncx = st::internal::num_x_stacks; // temporary variables for readability
   int ncy = st::internal::num_y_stacks;
   int ncz = st::internal::num_microcells_per_stack;

   // Set cell and stack counters
   int cell = 0;
   int stack_x = 0;
   int stack_y = 0;
      st::internal::cell_stack_index.resize(ncx*ncy*ncz);
   // Allocate space for 3D supercell array (ST coordinate system)
   std::vector<std::vector<std::vector<int> > > supercell_array;
   supercell_array.resize(ncx);
   st::internal::cell_index_x.resize(ncy*ncz);

   for(int i=0;i<ncx;++i){
      supercell_array[i].resize(ncy);
      stack_x = 0;

      for(int j=0;j<ncy;++j){
         
         supercell_array[i][j].resize(ncz);
         // set starting cell for y stack
         st::internal::stack_index_y.at(stack_y)=cell; 
        
         // increment y stack counter
         stack_y++;
         // store cell coordinates
         for(int k=0; k<ncz; ++k){
            if(i ==0) st::internal::stack_index_x.at(stack_x) = cell; 
            st::internal::cell_index_x.at(stack_x).push_back(cell); 
            stack_x++;
            // associate cell with position i,j,k
            supercell_array[i][j][k]=cell;
            // save ijk coordinates as microcell positions
            st::internal::pos.at(3*cell+0)=i;
            st::internal::pos.at(3*cell+1)=j;
            st::internal::pos.at(3*cell+2)=k;
            st::internal::cell_stack_index.at(cell) = stack_y;
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

      // temporary for atom coordinates
      double c[3];
      // convert atom coordinates to st reference frame
      c[st::internal::stx]=atom_coords_x[atom]+0.0000;
      c[st::internal::sty]=atom_coords_y[atom]+0.000;
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
      st::internal::atom_st_index[atom]= supercell_array[scc[0]][scc[1]][scc[2]]; 
      
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

      //added mpi decomposition for stacks
   #ifdef MPICF 
     
      //determine mpi core stack_y list 
      int removed_stacks_y = st::internal::num_stacks_y;
      std::vector<int> mpi_stack_id_y(removed_stacks_y, 0);
   
      if(vmpi::num_processors > removed_stacks_y ) {
         std::cout << "mpirun threads requested larger than spin-torque decomposition allows" << std::endl;
         err::vexit();
      }
      int residual = removed_stacks_y % vmpi::num_processors;
   // st::internal::mpi_stack_list.resize(int(floor(st::internal::num_stacks/vmpi::num_processors))+1);
      for(int s = 0; s < int(floor(removed_stacks_y/vmpi::num_processors)); s++) {
         st::internal::mpi_stack_list_y.push_back(s*vmpi::num_processors + vmpi::my_rank);
         mpi_stack_id_y.at(s*vmpi::num_processors + vmpi::my_rank) = 1;
      } 
      if(residual > 0 && vmpi::my_rank < residual) {
         st::internal::mpi_stack_list_y.push_back(removed_stacks_y - vmpi::my_rank -1);
         mpi_stack_id_y.at(removed_stacks_y - vmpi::my_rank -1) = 1;
      }   
   

      int stack_sum_y = 0;
      int size_y = st::internal::mpi_stack_list_y.size();
     
      MPI_Reduce(&size_y,&stack_sum_y, 1,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &mpi_stack_id_y[0],  mpi_stack_id_y.size(),   MPI_INT,MPI_SUM, MPI_COMM_WORLD);
     
      bool error = false;
      if(vmpi::my_rank == 0 && stack_sum_y != removed_stacks_y ) {
         std::cout << stack_sum_y << " != " << removed_stacks_y << std::endl;
         error = true;
      }  

      for(int i = 0; i < mpi_stack_id_y.size(); i++) {
         if(mpi_stack_id_y[i] != 1) error = true;
      }

      MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);
      if(error) {

         for(int i = 0; i < mpi_stack_id_y.size(); i++) std::cout << "mpi y stack: " << i << ", " << mpi_stack_id_y[i] << std::endl; 
         // for(int i = 0; i < mpi_stack_id_x.size(); i++) std::cout << "mpi x stack: " << i << ", " << mpi_stack_id_x[i] << std::endl; 
         MPI_Barrier(MPI_COMM_WORLD);
         err::vexit();
      }
      std::cout << "MPI y stack list " << st::internal::mpi_stack_list_y.size()  << std::endl;
      
      if(st::internal::sot_sa) {

         int removed_stacks_x = st::internal::num_stacks_x;
         std::vector<int> mpi_stack_id_x(removed_stacks_x, 0);

         if(vmpi::num_processors > removed_stacks_x ) {
            std::cout << "mpirun threads requested larger than spin-torque decomposition allows" << std::endl;
            err::vexit();
         }
         residual = removed_stacks_x % vmpi::num_processors;
         if(residual > 0) std::cout << "may have mpi decomp problem with residual " << residual << std::endl;
         // st::internal::mpi_stack_list.resize(int(floor(st::internal::num_stacks/vmpi::num_processors))+1);
         for(int s = 0; s < int(floor(removed_stacks_x/vmpi::num_processors)); s++) {
            st::internal::mpi_stack_list_x.push_back(s*vmpi::num_processors + vmpi::my_rank);
            mpi_stack_id_x.at(s*vmpi::num_processors + vmpi::my_rank) = 1;
         } 
         if(residual > 0 && vmpi::my_rank < residual){
            st::internal::mpi_stack_list_x.push_back(removed_stacks_x - vmpi::my_rank -1);
            mpi_stack_id_x.at(removed_stacks_x - vmpi::my_rank -1) = 1;
         }
         int stack_sum_x = 0;
         int size_x = st::internal::mpi_stack_list_x.size();
         MPI_Reduce(&size_x,&stack_sum_x, 1,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

         if(vmpi::my_rank == 0 && stack_sum_x != removed_stacks_x ) {
               std::cout << stack_sum_x << " != " << removed_stacks_x << std::endl;
               error = true;
            }
            MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);
            if(error) {

               // for(int i = 0; i < mpi_stack_id_y.size(); i++) std::cout << "mpi y stack: " << i << ", " << mpi_stack_id_y[i] << std::endl; 
               for(int i = 0; i < mpi_stack_id_x.size(); i++) std::cout << "mpi x stack: " << i << ", " << mpi_stack_id_x[i] << std::endl; 
               MPI_Barrier(MPI_COMM_WORLD);
               err::vexit();
            }
         MPI_Allreduce(MPI_IN_PLACE, &mpi_stack_id_x[0],  mpi_stack_id_x.size(),   MPI_INT,MPI_SUM, MPI_COMM_WORLD);
         

         for(int i = 0; i < mpi_stack_id_x.size(); i++) {
               if(mpi_stack_id_x[i] != 1) error = true;
            }
         MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);
            if(error) {

               // for(int i = 0; i < mpi_stack_id_y.size(); i++) std::cout << "mpi y stack: " << i << ", " << mpi_stack_id_y[i] << std::endl; 
               for(int i = 0; i < mpi_stack_id_x.size(); i++) std::cout << "mpi x stack: " << i << ", " << mpi_stack_id_x[i] << std::endl; 
               MPI_Barrier(MPI_COMM_WORLD);
               err::vexit();
            }
      // for(int i = 0; i < st::internal::mpi_stack_list_x.size(); i++ ) std::cout << "MPI " << vmpi::my_rank << ", x stack values: " << st::internal::mpi_stack_list_x[i]  << std::endl;
      std::cout << "MPI x stack list avg. " << st::internal::mpi_stack_list_x.size()  << std::endl;
         
      }
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
      st::internal::default_properties.beta_cond =  0.99;
      st::internal::default_properties.beta_diff = 0.99;
      st::internal::default_properties.sa_infinity =  1.48e7;
      st::internal::default_properties.lambda_sdl = 2000.0e-10; // m
      st::internal::default_properties.diffusion =  0.001; //m^2/s ? 
      st::internal::default_properties.sd_exchange = 8.010883e-25; //Joule
      
      st::internal::default_properties.sot_beta_cond =  0.99;
      st::internal::default_properties.sot_beta_diff = 0.99;
      st::internal::default_properties.sot_sa_infinity =  1.48e7;
      st::internal::default_properties.sot_lambda_sdl = 2000.0e-10; // m
      st::internal::default_properties.sot_diffusion =  0.001; //m^2/s ? 
      st::internal::default_properties.sot_sd_exchange = 8.010883e-25; //Joule
     
      // Temporary array to hold number of atoms in each cell for averaging
      std::vector<double> count(st::internal::beta_cond.size(),0.0);

      // loop over all atoms
      for(int atom=0;atom<num_local_atoms;atom++) {

         // get material type
         int mat = atom_type_array[atom];

         // get microcell id
         int id = st::internal::atom_st_index[atom];
   
         //STT determine atomic properties
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

         //SOT
         if(st::internal::sot_sa) {
            double sot_beta_cond = st::internal::mp.at(mat).sot_beta_cond; // beta
            double sot_beta_diff = st::internal::mp.at(mat).sot_beta_diff; // beta_prime
            double sot_sa_infinity = st::internal::mp.at(mat).sot_sa_infinity;
            double sot_lambda_sdl = st::internal::mp.at(mat).sot_lambda_sdl;
            double sot_diffusion = st::internal::mp.at(mat).sot_diffusion;
            double sot_sd_exchange = st::internal::mp.at(mat).sot_sd_exchange;

            //add atomic properties to microcells
            st::internal::sot_beta_cond.at(id) += sot_beta_cond;
            st::internal::sot_beta_diff.at(id) += sot_beta_diff;
            st::internal::sot_sa_infinity.at(id) += sot_sa_infinity;
            st::internal::sot_lambda_sdl.at(id) += sot_lambda_sdl;
            st::internal::sot_diffusion.at(id) += sot_diffusion;
            st::internal::sot_sd_exchange.at(id) += sot_sd_exchange;

            st::internal::spin_acc_sign.at(id) += (mat == 0) ? 0:((mat == 1) ? 1.0:-1.0); 
         }

         count.at(id) += (mat == 0) ? 1:1;
      }

      // reduce microcell properties on all CPUs
      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::beta_cond[0],   st::internal::beta_cond.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::beta_diff[0],   st::internal::beta_diff.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::sa_infinity[0], st::internal::sa_infinity.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::lambda_sdl[0],  st::internal::lambda_sdl.size(),  MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::diffusion[0],   st::internal::diffusion.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &st::internal::sd_exchange[0], st::internal::sd_exchange.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
          MPI_Allreduce(MPI_IN_PLACE, &count[0],                     count.size(),                     MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
      #endif

      if(st::internal::sot_sa) {
         #ifdef MPICF
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::spin_acc_sign[0], st::internal::spin_acc_sign.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::sot_beta_cond[0],   st::internal::sot_beta_cond.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::sot_beta_diff[0],   st::internal::sot_beta_diff.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::sot_sa_infinity[0], st::internal::sot_sa_infinity.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::sot_lambda_sdl[0],  st::internal::sot_lambda_sdl.size(),  MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::sot_diffusion[0],   st::internal::sot_diffusion.size(),   MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &st::internal::sot_sd_exchange[0], st::internal::sot_sd_exchange.size(), MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
         #endif
      }
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
            st::internal::default_properties.sa_infinity =  st::internal::sa_infinity.at(cell);
            // if(st::internal::spin_acc_sign.at(cell) != 1.0 && st::internal::spin_acc_sign.at(cell) != -1.0) std::cout << st::internal::spin_acc_sign.at(cell) << std::endl;
         } else{
            st::internal::beta_cond.at(cell)   = st::internal::default_properties.beta_cond;
            st::internal::beta_diff.at(cell)   = st::internal::default_properties.beta_diff;
            st::internal::sa_infinity.at(cell) = st::internal::default_properties.sa_infinity;
            st::internal::lambda_sdl.at(cell)  = st::internal::default_properties.lambda_sdl;
            st::internal::diffusion.at(cell)   = st::internal::default_properties.diffusion;
            st::internal::sd_exchange.at(cell) = st::internal::default_properties.sd_exchange;
         }
         if(st::internal::sot_sa) {
            if(nat>0.0001){
               //  if(nat != 2)   std::cout << nat << std::endl;
               st::internal::sot_beta_cond.at(cell)   /= nat;
               st::internal::sot_beta_diff.at(cell)   /= nat;
               st::internal::sot_sa_infinity.at(cell) /= nat;
               st::internal::sot_lambda_sdl.at(cell)  /= nat;
               st::internal::sot_diffusion.at(cell)   /= nat;
               st::internal::sot_sd_exchange.at(cell) /= nat;
               st::internal::spin_acc_sign.at(cell) = std::round(st::internal::spin_acc_sign.at(cell)/count.at(cell));
               st::internal::default_properties.sot_sa_infinity =  st::internal::sot_sa_infinity.at(cell);
               // if(st::internal::spin_acc_sign.at(cell) != 1.0 && st::internal::spin_acc_sign.at(cell) != -1.0) std::cout << st::internal::spin_acc_sign.at(cell) << std::endl;
            } else{
               st::internal::sot_beta_cond.at(cell)   = st::internal::default_properties.sot_beta_cond;
               st::internal::sot_beta_diff.at(cell)   = st::internal::default_properties.sot_beta_diff;
               st::internal::sot_sa_infinity.at(cell) = st::internal::default_properties.sot_sa_infinity;
               st::internal::sot_lambda_sdl.at(cell)  = st::internal::default_properties.sot_lambda_sdl;
               st::internal::sot_diffusion.at(cell)   = st::internal::default_properties.sot_diffusion;
               st::internal::sot_sd_exchange.at(cell) = st::internal::default_properties.sot_sd_exchange;
            }
            if(st::internal::spin_acc_sign.at(cell) == 0) st::internal::sot_sa_source.at(cell) = true;
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

         const double l_perp = 0;
         const double l_L = 0;
         const double lambda_phi2 = 2*lambda_j2;
         const double lambda_trans2 = lambda_phi2 + lambda_sf2;
         std::complex<double> inside (1.0/lambda_sf2, -1.0/lambda_j2);
         std::complex<double> inv_lplus = sqrt(inside);

         st::internal::a[cell] =  real(inv_lplus);
         st::internal::b[cell] = -imag(inv_lplus);

      //set sot-sa parameters
        if(st::internal::sot_sa){
         int mat = 0;

         const double sot_B  = st::internal::sot_beta_cond[cell];
         const double sot_Bp = st::internal::sot_beta_diff[cell];
         const double sot_lambda_sdl = st::internal::sot_lambda_sdl[cell];
         const double sot_Do = st::internal::sot_diffusion[cell];
         const double sot_Jsd = st::internal::sot_sd_exchange[cell];

         const double sot_BBp = 1.0/sqrt(1.0-sot_B*sot_Bp);
         const double sot_lambda_sf = sot_lambda_sdl*sot_BBp;
         const double sot_lambda_j = sqrt(2.0*hbar*sot_Do /sot_Jsd); // Angstroms
         const double sot_lambda_sf2 = sot_lambda_sf*sot_lambda_sf;
         const double sot_lambda_j2 = sot_lambda_j*sot_lambda_j;

         const double l_perp = 0;
         const double l_L = 0;
         const double sot_lambda_phi2 = 2*lambda_j2;
         const double sot_lambda_trans2 = lambda_phi2 + lambda_sf2;
         // std::complex<double> inside (1.0/lambda_trans2, -1.0/lambda_j2); 

         std::complex<double> sot_inside (1.0/sot_lambda_sf2, -1.0/sot_lambda_j2);
         std::complex<double> sot_inv_lplus = sqrt(sot_inside);

         st::internal::sot_a[cell] =  real(sot_inv_lplus);
         st::internal::sot_b[cell] = -imag(sot_inv_lplus);
         }
      }
      return;
   }
}

} // end of st namespace
