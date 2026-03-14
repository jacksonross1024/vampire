//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Jack B. Collings 2025. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//
// Shared-memory UCF initialisation: rank 0 reads the UCF locally (no
// MPI broadcast of the file text), broadcasts small metadata to all ranks,
// and sends the large interaction arrays to node leaders in chunks that
// respect MPI buffer limits.  Each node leader places the data into an
// MPI_Win_allocate_shared window; non-leaders obtain a pointer via
// MPI_Win_shared_query.

// C++ standard library headers
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

// Vampire headers
#include "errors.hpp"
#include "exchange.hpp"
#include "material.hpp"
#include "unitcell.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// unitcell module headers
#include "internal.hpp"

#ifdef MPICF

namespace unitcell{

//--------------------------------------------------------------------------
// Packed scalar metadata that every rank needs (tiny, safe to broadcast).
//--------------------------------------------------------------------------
struct ucf_metadata_t {
   double dimensions[3];
   double shape[3][3];
   unsigned int interaction_range;
   unsigned int lcsize;
   unsigned int hcsize;
   unsigned int surface_threshold;
   int bilinear_exchange_type;
   int biquadratic_exchange_type;
   bool bilinear_use_material;
   bool biquadratic_use_material;
   uint64_t num_atoms;
   uint64_t num_bilinear;
   uint64_t num_biquadratic;
   uint64_t num_bilinear_ni;
   uint64_t num_biquadratic_ni;
};

//--------------------------------------------------------------------------
// Chunked point-to-point send/recv for arbitrarily large byte buffers.
// Mirrors the pattern in vmpi::broadcast (wrapper.cpp) -- 1 GB chunks.
//--------------------------------------------------------------------------
static const uint64_t chunk_size = 1073741824ULL; // 1 GB

static void chunked_send(const void* buf, uint64_t total_bytes, int dest, int tag, MPI_Comm comm){
   const char* ptr = static_cast<const char*>(buf);
   uint64_t remaining = total_bytes;
   while(remaining > chunk_size){
      MPI_Send(ptr, static_cast<int>(chunk_size), MPI_BYTE, dest, tag, comm);
      ptr += chunk_size;
      remaining -= chunk_size;
   }
   if(remaining > 0){
      MPI_Send(ptr, static_cast<int>(remaining), MPI_BYTE, dest, tag, comm);
   }
}

static void chunked_recv(void* buf, uint64_t total_bytes, int source, int tag, MPI_Comm comm){
   char* ptr = static_cast<char*>(buf);
   uint64_t remaining = total_bytes;
   while(remaining > chunk_size){
      MPI_Recv(ptr, static_cast<int>(chunk_size), MPI_BYTE, source, tag, comm, MPI_STATUS_IGNORE);
      ptr += chunk_size;
      remaining -= chunk_size;
   }
   if(remaining > 0){
      MPI_Recv(ptr, static_cast<int>(remaining), MPI_BYTE, source, tag, comm, MPI_STATUS_IGNORE);
   }
}

//--------------------------------------------------------------------------
// Helper: set up an MPI shared-memory window on node_comm.
// Leader allocates total_bytes; non-leaders allocate 0.
// Returns the usable pointer (leader's allocation) on all ranks.
//
// If the allocation fails (e.g. /dev/shm too small), prints a
// diagnostic and aborts via err::zexit.
//--------------------------------------------------------------------------
static void setup_shared_window(
   MPI_Aint total_bytes,
   MPI_Win& win,
   void*& shared_ptr)
{
   MPI_Info info;
   MPI_Info_create(&info);
   MPI_Info_set(info, "alloc_shared_noncontig", "true");

   void* base_ptr = nullptr;
   MPI_Aint local_size = vmpi::node_leader ? total_bytes : 0;

   // Set a non-fatal error handler so we can catch failures ourselves
   MPI_Errhandler prev_handler;
   MPI_Comm_get_errhandler(vmpi::node_comm, &prev_handler);
   MPI_Comm_set_errhandler(vmpi::node_comm, MPI_ERRORS_RETURN);

   int mpi_err = MPI_Win_allocate_shared(
      local_size,
      1,
      info,
      vmpi::node_comm,
      &base_ptr,
      &win);

   // Restore previous error handler
   MPI_Comm_set_errhandler(vmpi::node_comm, prev_handler);

   MPI_Info_free(&info);

   if(mpi_err != MPI_SUCCESS){
      char err_string[MPI_MAX_ERROR_STRING];
      int err_len = 0;
      MPI_Error_string(mpi_err, err_string, &err_len);

      terminaltextcolor(RED);
      std::cerr << std::endl;
      std::cerr << "Error: MPI_Win_allocate_shared failed on rank " << vmpi::my_rank << std::endl;
      std::cerr << "  MPI error: " << std::string(err_string, err_len) << std::endl;
      std::cerr << "  Requested allocation: " << total_bytes / (1024ULL * 1024ULL) << " MB" << std::endl;
      std::cerr << "  This typically means /dev/shm is too small for the shared-memory UCF." << std::endl;
      std::cerr << "  Check with:  df -h /dev/shm" << std::endl;
      std::cerr << "  Fix with:    mount -o remount,size=<needed>G /dev/shm  (requires root)" << std::endl;
      std::cerr << "  Or disable:  remove 'sim:shared-memory-ucf' from the input file." << std::endl;
      terminaltextcolor(WHITE);

      zlog << zTs() << "Error: MPI_Win_allocate_shared failed on rank " << vmpi::my_rank << std::endl;
      zlog << zTs() << "  MPI error: " << std::string(err_string, err_len) << std::endl;
      zlog << zTs() << "  Requested: " << total_bytes / (1024ULL * 1024ULL) << " MB" << std::endl;
      zlog << zTs() << "  /dev/shm is likely too small. See stderr for remediation steps." << std::endl;

      err::zexit("MPI_Win_allocate_shared failed -- cannot allocate shared-memory UCF window");
   }

   if(vmpi::node_leader){
      shared_ptr = base_ptr;
   }
   else{
      MPI_Aint sz;
      int disp;
      MPI_Win_shared_query(win, 0, &sz, &disp, &shared_ptr);
   }
}

//--------------------------------------------------------------------------
// Read UCF file locally on calling process (no MPI calls).
// This replaces vin::get_string which would MPI_Bcast the entire file.
//--------------------------------------------------------------------------
static void read_unit_cell_local(unit_cell_t& unit_cell, std::string filename){

   std::cout << "Reading in unit cell data from disk (local)..." << std::flush;
   zlog << zTs() << "Reading in unit cell data from disk (local, no MPI broadcast)..." << std::endl;

   std::ifstream infile(filename.c_str());
   if(!infile.is_open()){
      terminaltextcolor(RED);
      std::cerr << "Error opening unit cell file \"" << filename << "\". Exiting." << std::endl;
      terminaltextcolor(WHITE);
      err::vexit();
   }

   // Read entire file into a stringstream (same approach as vin::get_string but local-only)
   std::string contents((std::istreambuf_iterator<char>(infile)), std::istreambuf_iterator<char>());
   infile.close();

   std::stringstream inputfile;
   inputfile.str(contents);

   std::cout << "done!\nProcessing unit cell data..." << std::flush;
   zlog << zTs() << "Reading data completed. Processing unit cell data..." << std::endl;

   uint64_t line_counter = 0;
   unsigned int line_id = 0;
   unsigned int interaction_range = 1;

   while(!inputfile.eof()){
      line_counter++;
      std::string line;
      getline(inputfile, line);

      std::string empty = "";
      if(line == empty) continue;

      const char* hash = "#";
      bool has_hash = false;
      for(unsigned int i = 0; i < line.length(); i++){
         if(line.at(i) == *hash){ has_hash = true; break; }
      }
      if(has_hash) continue;

      std::istringstream iss(line, std::istringstream::in);

      switch(line_id){
         case 0:
            iss >> unit_cell.dimensions[0] >> unit_cell.dimensions[1] >> unit_cell.dimensions[2];
            break;
         case 1:
            iss >> unit_cell.shape[0][0] >> unit_cell.shape[0][1] >> unit_cell.shape[0][2];
            break;
         case 2:
            iss >> unit_cell.shape[1][0] >> unit_cell.shape[1][1] >> unit_cell.shape[1][2];
            break;
         case 3:
            iss >> unit_cell.shape[2][0] >> unit_cell.shape[2][1] >> unit_cell.shape[2][2];
            break;
         case 4:{
            unsigned int num_uc_atoms;
            iss >> num_uc_atoms;
            if(num_uc_atoms > 0 && num_uc_atoms <= 100000000) unit_cell.atom.resize(num_uc_atoms);
            else{
               terminaltextcolor(RED);
               std::cerr << "Error! Requested number of atoms " << num_uc_atoms << " on line " << line_counter
                         << " of unit cell input file " << filename << " is outside of valid range 1-100,000,000. Exiting" << std::endl;
               terminaltextcolor(WHITE);
               err::vexit();
            }

            std::cout << "\nProcessing data for " << unit_cell.atom.size() << " atoms..." << std::flush;
            zlog << zTs() << "\tProcessing data for " << unit_cell.atom.size() << " unit cell atoms..." << std::endl;

            for(uint64_t i = 0; i < unit_cell.atom.size(); i++){
               line_counter++;
               uint64_t id = i;
               double cx = 2.0, cy = 2.0, cz = 2.0;
               int mat_id = 0, lcat_id = 0, hcat_id = 0;
               std::string atom_line;
               getline(inputfile, atom_line);
               std::istringstream atom_iss(atom_line, std::istringstream::in);
               atom_iss >> id >> cx >> cy >> cz >> mat_id >> lcat_id >> hcat_id;

               if(cx >= 0.0 && cx <= 1.0) unit_cell.atom[i].x = cx;
               else{ terminaltextcolor(RED); std::cerr << "Error! atom x-coordinate for atom " << id << " on line " << line_counter << " of " << filename << " is outside valid range 0.0-1.0. Exiting" << std::endl; terminaltextcolor(WHITE); err::vexit(); }
               if(cy >= 0.0 && cy <= 1.0) unit_cell.atom[i].y = cy;
               else{ terminaltextcolor(RED); std::cerr << "Error! atom y-coordinate for atom " << id << " on line " << line_counter << " of " << filename << " is outside valid range 0.0-1.0. Exiting" << std::endl; terminaltextcolor(WHITE); err::vexit(); }
               if(cz >= 0.0 && cz <= 1.0) unit_cell.atom[i].z = cz;
               else{ terminaltextcolor(RED); std::cerr << "Error! atom z-coordinate for atom " << id << " on line " << line_counter << " of " << filename << " is outside valid range 0.0-1.0. Exiting" << std::endl; terminaltextcolor(WHITE); err::vexit(); }
               if(mat_id >= 0 && mat_id < mp::num_materials) unit_cell.atom[i].mat = mat_id;
               else{ terminaltextcolor(RED); std::cerr << "Error! material id " << mat_id << " for atom " << id << " on line " << line_counter << " of " << filename << " >= num_materials (" << mp::num_materials << "). Exiting" << std::endl; terminaltextcolor(WHITE); err::vexit(); }
               unit_cell.atom[i].lc = lcat_id;
               unit_cell.atom[i].hc = hcat_id;
            }
            break;
         }
         case 5:{
            unit_cell.bilinear.read_interactions(unit_cell.atom.size(), inputfile, iss, filename, line_counter, interaction_range);
            break;
         }
         case 6:{
            unit_cell.biquadratic.read_interactions(unit_cell.atom.size(), inputfile, iss, filename, line_counter, interaction_range);
            break;
         }
         default:
            terminaltextcolor(RED);
            std::cerr << "Error! Unknown line type on line " << line_counter << " of unit cell input file " << filename << ". Exiting" << std::endl;
            terminaltextcolor(WHITE);
            err::vexit();
      }
      line_id++;
   }

   // Skip verify/normalise/shells (the shared-memory path treats the UCF as final)
   std::cout << "done!" << std::endl;
   zlog << zTs() << "\tProcessing unit cell interactions completed (verify/normalise/shells bypassed)" << std::endl;

   if(unit_cell.biquadratic.interaction.size() > 0){
      zlog << zTs() << "Enabling biquadratic interactions from unit cell file" << std::endl;
      exchange::biquadratic = true;
   }

   if(interaction_range > unit_cell.interaction_range) unit_cell.interaction_range = interaction_range;

   zlog << zTs() << "\tNumber of atoms read-in: " << unit_cell.atom.size() << std::endl;
   zlog << zTs() << "\tNumber of bilinear interactions read-in: " << unit_cell.bilinear.interaction.size() << std::endl;
   zlog << zTs() << "\tNumber of biquadratic interactions read-in: " << unit_cell.biquadratic.interaction.size() << std::endl;
   zlog << zTs() << "\tCalculated interaction range: " << unit_cell.interaction_range << " Unit Cells" << std::endl;
}

//--------------------------------------------------------------------------
// Main shared-memory UCF initialisation
//--------------------------------------------------------------------------
void initialise_shared(unit_cell_t& unit_cell){

   if(vmpi::my_rank == 0){
      std::cout << "Initialising shared-memory UCF mode" << std::endl;
   }
   zlog << zTs() << "Initialising shared-memory UCF mode" << std::endl;

   //----------------------------------------------------------------------
   // Log warnings about bypassed functionality
   //----------------------------------------------------------------------
   if(vmpi::my_rank == 0){
      std::cout << "   WARNING: shared-memory-ucf mode bypasses verify, normalise, shells, dmi, and kitaev" << std::endl;
      std::cout << "            processing of the UCF. The unit cell file must contain final interaction" << std::endl;
      std::cout << "            values; DMI and Kitaev exchange modifications will NOT be applied." << std::endl;
   }
   zlog << zTs() << "WARNING: shared-memory-ucf mode bypasses the following UCF post-processing steps:" << std::endl;
   zlog << zTs() << "   - verify (exchange symmetry checking)" << std::endl;
   zlog << zTs() << "   - normalise (exchange normalisation)" << std::endl;
   zlog << zTs() << "   - shells (crystal shell identification)" << std::endl;
   zlog << zTs() << "   - dmi (Dzyaloshinskii-Moriya interaction calculation)" << std::endl;
   zlog << zTs() << "   - kitaev (Kitaev exchange calculation)" << std::endl;
   zlog << zTs() << "The unit cell file must contain final, pre-computed interaction values." << std::endl;

   //----------------------------------------------------------------------
   // Log node leader determination
   //----------------------------------------------------------------------
   const int num_nodes_log = (vmpi::num_processors - 1) / vmpi::ppn + 1;
   zlog << zTs() << "Node communicator: ppn = " << vmpi::ppn
        << ", node_rank = " << vmpi::node_rank
        << ", node_size = " << vmpi::node_size
        << ", node_leader = " << (vmpi::node_leader ? "true" : "false")
        << " (global rank " << vmpi::my_rank << "); num_processors = " << vmpi::num_processors
        << " -> " << num_nodes_log << " nodes" << std::endl;
   if(vmpi::node_leader){
      zlog << zTs() << "   Rank " << vmpi::my_rank << " is node leader (node "
           << vmpi::my_rank / vmpi::ppn << " of " << num_nodes_log << ")" << std::endl;
   }

   //----------------------------------------------------------------------
   // 1. Rank 0 reads UCF locally (no MPI calls -- bypasses vin::get_string)
   //----------------------------------------------------------------------
   ucf_metadata_t meta;
   std::memset(&meta, 0, sizeof(meta));

   if(vmpi::my_rank == 0){

      std::string blank = "";
      if(uc::internal::unit_cell_filename.c_str() == blank){
         terminaltextcolor(RED);
         std::cerr << "Error: shared-memory-ucf mode requires a unit cell file. Exiting." << std::endl;
         terminaltextcolor(WHITE);
         err::vexit();
      }

      read_unit_cell_local(unit_cell, uc::internal::unit_cell_filename);

      // Pack metadata
      std::memcpy(meta.dimensions, unit_cell.dimensions, 3 * sizeof(double));
      std::memcpy(meta.shape, unit_cell.shape, 9 * sizeof(double));
      meta.interaction_range         = unit_cell.interaction_range;
      meta.lcsize                    = unit_cell.lcsize;
      meta.hcsize                    = unit_cell.hcsize;
      meta.surface_threshold         = unit_cell.surface_threshold;
      meta.bilinear_exchange_type    = static_cast<int>(unit_cell.bilinear.exchange_type);
      meta.biquadratic_exchange_type = static_cast<int>(unit_cell.biquadratic.exchange_type);
      meta.bilinear_use_material     = unit_cell.bilinear.use_material_exchange_constants;
      meta.biquadratic_use_material  = unit_cell.biquadratic.use_material_exchange_constants;
      meta.num_atoms                 = unit_cell.atom.size();
      meta.num_bilinear              = unit_cell.bilinear.interaction.size();
      meta.num_biquadratic           = unit_cell.biquadratic.interaction.size();
      meta.num_bilinear_ni           = unit_cell.bilinear.ni.size();
      meta.num_biquadratic_ni        = unit_cell.biquadratic.ni.size();

      std::cout << "   UCF parsed: " << meta.num_atoms << " atoms, "
                << meta.num_bilinear << " bilinear interactions ("
                << (meta.num_bilinear * sizeof(interaction_t)) / (1024ULL * 1024ULL) << " MB), "
                << meta.num_biquadratic << " biquadratic interactions" << std::endl;
      zlog << zTs() << "   UCF parsed: " << meta.num_atoms << " atoms, "
           << meta.num_bilinear << " bilinear interactions, "
           << meta.num_biquadratic << " biquadratic interactions" << std::endl;
   }

   //----------------------------------------------------------------------
   // 2. Broadcast small metadata to all ranks
   //----------------------------------------------------------------------
   MPI_Bcast(&meta, sizeof(ucf_metadata_t), MPI_BYTE, 0, MPI_COMM_WORLD);

   if(vmpi::my_rank != 0){
      std::memcpy(unit_cell.dimensions, meta.dimensions, 3 * sizeof(double));
      std::memcpy(unit_cell.shape, meta.shape, 9 * sizeof(double));
      unit_cell.interaction_range                          = meta.interaction_range;
      unit_cell.lcsize                                     = meta.lcsize;
      unit_cell.hcsize                                     = meta.hcsize;
      unit_cell.surface_threshold                          = meta.surface_threshold;
      unit_cell.bilinear.exchange_type                     = static_cast<exchange::exchange_t>(meta.bilinear_exchange_type);
      unit_cell.biquadratic.exchange_type                  = static_cast<exchange::exchange_t>(meta.biquadratic_exchange_type);
      unit_cell.bilinear.use_material_exchange_constants   = meta.bilinear_use_material;
      unit_cell.biquadratic.use_material_exchange_constants = meta.biquadratic_use_material;
   }

   //----------------------------------------------------------------------
   // 3. Broadcast atom array + ni vectors (small)
   //----------------------------------------------------------------------
   if(vmpi::my_rank != 0){
      unit_cell.atom.resize(meta.num_atoms);
   }
   MPI_Bcast(unit_cell.atom.data(),
             meta.num_atoms * sizeof(atom_t),
             MPI_BYTE, 0, MPI_COMM_WORLD);

   if(vmpi::my_rank != 0){
      unit_cell.bilinear.ni.resize(meta.num_bilinear_ni);
      unit_cell.biquadratic.ni.resize(meta.num_biquadratic_ni);
   }
   if(meta.num_bilinear_ni > 0){
      MPI_Bcast(unit_cell.bilinear.ni.data(),
                meta.num_bilinear_ni * sizeof(int),
                MPI_BYTE, 0, MPI_COMM_WORLD);
   }
   if(meta.num_biquadratic_ni > 0){
      MPI_Bcast(unit_cell.biquadratic.ni.data(),
                meta.num_biquadratic_ni * sizeof(int),
                MPI_BYTE, 0, MPI_COMM_WORLD);
   }

   if(meta.num_biquadratic > 0){
      exchange::biquadratic = true;
   }

   //----------------------------------------------------------------------
   // 4. Send interaction arrays from rank 0 to node leaders only.
   //    Uses chunked send/recv (1 GB chunks) to avoid overflowing the
   //    MPI int count limit, matching the vmpi::broadcast pattern.
   //----------------------------------------------------------------------
   const int num_nodes = (vmpi::num_processors - 1) / vmpi::ppn + 1;

   std::vector<interaction_t> recv_bilinear;
   std::vector<interaction_t> recv_biquadratic;

   const uint64_t bilinear_bytes   = meta.num_bilinear   * sizeof(interaction_t);
   const uint64_t biquadratic_bytes = meta.num_biquadratic * sizeof(interaction_t);

   if(meta.num_bilinear > 0){
      if(vmpi::my_rank == 0){
         for(int n = 1; n < num_nodes; n++){
            int leader_rank = n * vmpi::ppn;
            if(leader_rank < vmpi::num_processors){
               chunked_send(unit_cell.bilinear.interaction.data(),
                            bilinear_bytes, leader_rank, 100, MPI_COMM_WORLD);
            }
         }
      }
      else if(vmpi::node_leader){
         recv_bilinear.resize(meta.num_bilinear);
         chunked_recv(recv_bilinear.data(),
                      bilinear_bytes, 0, 100, MPI_COMM_WORLD);
      }
   }

   if(meta.num_biquadratic > 0){
      if(vmpi::my_rank == 0){
         for(int n = 1; n < num_nodes; n++){
            int leader_rank = n * vmpi::ppn;
            if(leader_rank < vmpi::num_processors){
               chunked_send(unit_cell.biquadratic.interaction.data(),
                            biquadratic_bytes, leader_rank, 101, MPI_COMM_WORLD);
            }
         }
      }
      else if(vmpi::node_leader){
         recv_biquadratic.resize(meta.num_biquadratic);
         chunked_recv(recv_biquadratic.data(),
                      biquadratic_bytes, 0, 101, MPI_COMM_WORLD);
      }
   }

   //----------------------------------------------------------------------
   // 5. Set up MPI shared-memory windows
   //----------------------------------------------------------------------
   void* bilinear_shared = nullptr;
   void* biquadratic_shared = nullptr;

   if(meta.num_bilinear > 0){
      setup_shared_window(
         static_cast<MPI_Aint>(bilinear_bytes),
         vmpi::bilinear_win,
         bilinear_shared);
   }

   if(meta.num_biquadratic > 0){
      setup_shared_window(
         static_cast<MPI_Aint>(biquadratic_bytes),
         vmpi::biquadratic_win,
         biquadratic_shared);
   }

   //----------------------------------------------------------------------
   // 6. Leader fills the shared window under lock, then barrier
   //----------------------------------------------------------------------
   // Cray MPICH (and some other implementations) require lock/unlock for
   // shared-memory window synchronization; MPI_Win_sync alone can trigger
   // "Wrong synchronization of RMA calls" when no RMA ops are used.
   // Leader holds exclusive lock on rank 0 (self) while writing; barrier
   // ensures all ranks see the data before non-leaders read.
   if(meta.num_bilinear > 0){
      if(vmpi::node_leader){
         MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, vmpi::bilinear_win);
         if(vmpi::my_rank == 0){
            std::memcpy(bilinear_shared,
                        unit_cell.bilinear.interaction.data(),
                        bilinear_bytes);
         }
         else{
            std::memcpy(bilinear_shared,
                        recv_bilinear.data(),
                        bilinear_bytes);
            std::vector<interaction_t>().swap(recv_bilinear);
         }
         MPI_Win_unlock(0, vmpi::bilinear_win);
      }
   }

   if(meta.num_biquadratic > 0){
      if(vmpi::node_leader){
         MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, vmpi::biquadratic_win);
         if(vmpi::my_rank == 0){
            std::memcpy(biquadratic_shared,
                        unit_cell.biquadratic.interaction.data(),
                        biquadratic_bytes);
         }
         else{
            std::memcpy(biquadratic_shared,
                        recv_biquadratic.data(),
                        biquadratic_bytes);
            std::vector<interaction_t>().swap(recv_biquadratic);
         }
         MPI_Win_unlock(0, vmpi::biquadratic_win);
      }
   }

   MPI_Barrier(vmpi::node_comm);

   //----------------------------------------------------------------------
   // 7. Set interaction_ptr on all ranks to point into the shared window.
   //    The exchange unrolling code (unroll.cpp etc.) still produces its
   //    own per-rank std::vector structures (i_exchange_list, etc.) for
   //    the actual simulation -- only the template lookup uses the pointer.
   //----------------------------------------------------------------------
   if(meta.num_bilinear > 0){
      unit_cell.bilinear.set_shared(
         static_cast<const interaction_t*>(bilinear_shared),
         meta.num_bilinear);
   }
   else{
      unit_cell.bilinear.set_shared(nullptr, 0);
   }

   if(meta.num_biquadratic > 0){
      unit_cell.biquadratic.set_shared(
         static_cast<const interaction_t*>(biquadratic_shared),
         meta.num_biquadratic);
   }
   else{
      unit_cell.biquadratic.set_shared(nullptr, 0);
   }

   // Free the local interaction vectors on rank 0 (data is now in the window)
   if(vmpi::my_rank == 0){
      std::vector<interaction_t>().swap(unit_cell.bilinear.interaction);
      std::vector<interaction_t>().swap(unit_cell.biquadratic.interaction);
   }

   if(vmpi::my_rank == 0){
      std::cout << "Shared-memory UCF initialisation complete" << std::endl;
   }
   zlog << zTs() << "Shared-memory UCF initialisation complete" << std::endl;

   return;
}

} // end of unitcell namespace

#endif // MPICF
