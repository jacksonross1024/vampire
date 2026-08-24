//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// program header
#include "vdc.hpp"

// openmp header
#ifdef _OPENMP
   #include <omp.h>
#else
   #define omp_get_thread_num() 0
   #define omp_get_max_threads() 1
   #define omp_set_num_threads(n) (void)(n)
#endif

namespace vdc{

namespace {

const char* binary_meta_filename = "vdc-binary.meta";

struct binary_output_record {
   binary_dump_spec spec;
   std::string last_file;
   std::string pattern;
   uint64_t n_rows;
   uint64_t n_cols;
};

std::vector<binary_output_record> g_binary_outputs;

//------------------------------------------------------------------------------
const char* host_endianness(){
   const uint16_t one = 1;
   return (*reinterpret_cast<const unsigned char*>(&one) != 0) ? "little" : "big";
}

//------------------------------------------------------------------------------
std::string numpy_uint64_dtype(){
   return (std::string(host_endianness()) == "little") ? "<u8" : ">u8";
}

//------------------------------------------------------------------------------
std::string numpy_float64_dtype(){
   return (std::string(host_endianness()) == "little") ? "<f8" : ">f8";
}

//------------------------------------------------------------------------------
// Replace a run of 8 digits (vdc frame id) with %08d for glob-style patterns.
//------------------------------------------------------------------------------
std::string filename_pattern(const std::string& name){

   std::string out = name;
   for(std::string::size_type i = 0; i + 8 <= out.size(); i++){
      bool eight = true;
      for(int k = 0; k < 8; k++){
         if(!std::isdigit(static_cast<unsigned char>(out[i + static_cast<std::string::size_type>(k)]))){
            eight = false;
            break;
         }
      }
      if(eight){
         out.replace(i, 8, "%08d");
         break;
      }
   }
   return out;

}

//------------------------------------------------------------------------------
std::string yaml_quote(const std::string& s){

   bool quote = s.empty();
   for(size_t i = 0; i < s.size(); i++){
      const char c = s[i];
      if(c == ':' || c == '#' || c == '"' || c == '\\' || c == '\n' ||
         c == '{' || c == '}' || c == '[' || c == ']' || c == '&' ||
         c == '*' || c == '!' || c == '|' || c == '>' || c == '\''){
         quote = true;
         break;
      }
   }
   if(!quote) return s;

   std::string o = "\"";
   for(size_t i = 0; i < s.size(); i++){
      const char c = s[i];
      if(c == '"' || c == '\\') o.push_back('\\');
      if(c == '\n'){ o += "\\n"; continue; }
      o.push_back(c);
   }
   o.push_back('"');
   return o;

}

} // end of anonymous namespace

//------------------------------------------------------------------------------
// Replace the last path extension with .bin (append .bin if none).
//------------------------------------------------------------------------------
std::string binary_filename(const std::string& text_filename){

   const std::string::size_type slash = text_filename.find_last_of("/\\");
   const std::string::size_type dot = text_filename.find_last_of('.');

   if(dot == std::string::npos || (slash != std::string::npos && dot < slash)){
      return text_filename + ".bin";
   }

   return text_filename.substr(0, dot) + ".bin";

}

//------------------------------------------------------------------------------
// --binary enables OpenMP only inside output_binary_file. All other regions
// (text writers, SSC compute, cell magnetisation) are forced to 1 thread.
//------------------------------------------------------------------------------
void configure_output_threads(){

#ifdef _OPENMP
   const int default_n = omp_get_max_threads();
   int write_n = 1;

   if(vdc::binary_dump){
      write_n = (vdc::omp_threads > 0) ? vdc::omp_threads : default_n;
   }

   vdc::omp_threads = write_n;
   omp_set_num_threads(1);

   if(vdc::verbose){
      if(vdc::binary_dump){
         std::cout << "Binary output path: " << write_n << " OpenMP threads on write helper (other work serial)" << std::endl;
      }
      else{
         std::cout << "Text output path: OpenMP threads set to 1" << std::endl;
      }
   }
#endif

}

//------------------------------------------------------------------------------
// Session metadata for --binary dumps. Aimed at script authors / AI agents.
//------------------------------------------------------------------------------
void write_binary_metadata(){

   if(!vdc::binary_dump) return;

   std::ofstream out(binary_meta_filename);
   if(!out.is_open()){
      std::cerr << "Error - cannot open \"" << binary_meta_filename << "\" for writing." << std::endl;
      std::exit(EXIT_FAILURE);
   }

   const std::string endian = host_endianness();
   const std::string dt_u = numpy_uint64_dtype();
   const std::string dt_f = numpy_float64_dtype();

   out << "# vdc --binary metadata\n";
   out << "# Written for humans and AI agents building loaders/scripts.\n";
   out << "# Each .bin file is self-describing in shape: read n_rows, n_cols from the 16-byte header.\n";
   out << "# Column names, units, and meaning live here (they are not stored inside the .bin).\n";
   out << "#\n";
   out << "# Do not parse PoVRAY .pov/.ini files for numeric data; those remain text scene files.\n";
   out << "\n";
   out << "format_name: VDCBIN\n";
   out << "format_version: 1\n";
   out << "writer: vdc\n";
   out << "endianness: " << endian << "\n";
   out << "header_nbytes: 16\n";
   out << "header:\n";
   out << "  - {name: n_rows, dtype: uint64, nbytes: 8, offset: 0}\n";
   out << "  - {name: n_cols, dtype: uint64, nbytes: 8, offset: 8}\n";
   out << "payload:\n";
   out << "  dtype: float64\n";
   out << "  ieee754: true\n";
   out << "  nbytes_per_value: 8\n";
   out << "  order: C\n";
   out << "  layout: data[row * n_cols + col]\n";
   out << "  offset_bytes: 16\n";
   out << "file_nbytes: 16 + n_rows * n_cols * 8\n";
   out << "padding: none\n";
   out << "checksum: none\n";
   out << "numpy_dtype_uint64: " << dt_u << "\n";
   out << "numpy_dtype_float64: " << dt_f << "\n";
   out << "\n";
   out << "python_load: |\n";
   out << "  import numpy as np\n";
   out << "  DT_U = np.dtype(\"" << dt_u << "\")\n";
   out << "  DT_F = np.dtype(\"" << dt_f << "\")\n";
   out << "  def load_vdc_bin(path):\n";
   out << "      with open(path, \"rb\") as f:\n";
   out << "          n_rows, n_cols = (int(x) for x in np.fromfile(f, dtype=DT_U, count=2))\n";
   out << "          data = np.fromfile(f, dtype=DT_F)\n";
   out << "      if data.size != n_rows * n_cols:\n";
   out << "          raise ValueError(\"payload size mismatch in %s\" % path)\n";
   out << "      return data.reshape(n_rows, n_cols)\n";
   out << "  # df = pandas.DataFrame(load_vdc_bin(path), columns=[c[\"name\"] for c in output[\"columns\"]])\n";
   out << "\n";
   out << "system:\n";
   out << "  size_Angstrom: [" << vdc::system_size[0] << ", " << vdc::system_size[1] << ", " << vdc::system_size[2] << "]\n";
   out << "  centre_Angstrom: [" << vdc::system_centre[0] << ", " << vdc::system_centre[1] << ", " << vdc::system_centre[2] << "]\n";
   out << "  n_atoms_unsliced: " << vdc::num_atoms << "\n";
   out << "  n_nm_atoms_unsliced: " << vdc::num_nm_atoms << "\n";
   out << "  n_sliced_magnetic_atoms: " << vdc::sliced_atoms_list.size() << "\n";
   out << "  n_sliced_nonmagnetic_atoms: " << vdc::sliced_nm_atoms_list.size() << "\n";
   out << "  n_materials: " << vdc::materials.size() << "\n";
   if(vdc::cells){
      out << "  cell_size_Angstrom: [" << vdc::cell_size[0] << ", " << vdc::cell_size[1] << ", " << vdc::cell_size[2] << "]\n";
      out << "  nx_cells: " << vdc::nx_cells << "\n";
      out << "  ny_cells: " << vdc::ny_cells << "\n";
      out << "  nz_cells: " << vdc::nz_cells << "\n";
      out << "  n_occupied_cells: " << vdc::total_cells << "\n";
      out << "  n_cell_material_slots: " << (1u + static_cast<unsigned int>(vdc::materials.size())) << "  # materials + total\n";
   }
   out << "\n";
   out << "materials:\n";
   if(vdc::materials.empty()){
      out << "  []\n";
   }
   else{
      for(size_t i = 0; i < vdc::materials.size(); i++){
         out << "  - index: " << i << "\n";
         out << "    id_1based: " << (vdc::materials[i].id + 1) << "\n";
         out << "    name: " << yaml_quote(vdc::materials[i].name) << "\n";
         out << "    element: " << yaml_quote(vdc::materials[i].element) << "\n";
         out << "    moment_muB: " << vdc::materials[i].moment << "\n";
      }
   }
   out << "\n";
   out << "frames:\n";
   out << "  requested_start: " << vdc::vdc_start_file_id << "\n";
   out << "  requested_final: " << vdc::vdc_final_file_id << "\n";
   out << "  processed_start: " << vdc::start_file_id << "\n";
   out << "  processed_final: " << vdc::final_file_id << "\n";
   out << "  filename_frame_digits: 8\n";
   out << "\n";
   out << "active_outputs:\n";
   out << "  xyz: " << (vdc::xyz ? "true" : "false") << "\n";
   out << "  text_or_spins: " << (vdc::txt ? "true" : "false") << "\n";
   out << "  cells: " << (vdc::cellsf ? "true" : "false") << "\n";
   out << "  vtk: " << (vdc::vtk ? "true" : "false") << "\n";
   out << "  ssc: " << (vdc::ssc ? "true" : "false") << "\n";
   out << "  povray: " << (vdc::povray ? "true" : "false") << "\n";
   out << "  track: " << (vdc::track ? "true" : "false") << "\n";
   out << "  povray_cells: " << (vdc::povcells ? "true" : "false") << "\n";
   out << "  spin_cells: " << (vdc::povspincells ? "true" : "false") << "\n";
   out << "  binary: true\n";
   out << "  omp_write_threads: " << vdc::omp_threads << "\n";
   out << "\n";
   out << "filename_collision_note: |\n";
   out << "  --text/--spins, --vtk, --povray, --track and --spin-cells all write spins-%08d.bin.\n";
   out << "  If more than one of those is enabled, the last writer in the snapshot loop overwrites the file.\n";
   out << "  Use last_file and kind in outputs[] to see which kind actually produced each pattern.\n";
   out << "\n";
   out << "outputs:\n";
   if(g_binary_outputs.empty()){
      out << "  []\n";
   }
   for(size_t o = 0; o < g_binary_outputs.size(); o++){
      const binary_output_record& rec = g_binary_outputs[o];
      out << "  - kind: " << yaml_quote(rec.spec.kind) << "\n";
      out << "    last_file: " << yaml_quote(rec.last_file) << "\n";
      out << "    filename_pattern: " << yaml_quote(rec.pattern) << "\n";
      out << "    n_rows: " << rec.n_rows << "\n";
      out << "    n_cols: " << rec.n_cols << "\n";
      if(rec.spec.columns.size() != rec.n_cols){
         out << "    column_count_mismatch: true  # spec has " << rec.spec.columns.size() << " names\n";
      }
      out << "    notes: " << yaml_quote(rec.spec.notes) << "\n";
      out << "    columns:\n";
      if(rec.spec.columns.empty()){
         out << "      []\n";
      }
      for(size_t c = 0; c < rec.spec.columns.size(); c++){
         const binary_column& col = rec.spec.columns[c];
         out << "      - index: " << c << "\n";
         out << "        name: " << yaml_quote(col.name) << "\n";
         out << "        unit: " << yaml_quote(col.unit) << "\n";
         out << "        description: " << yaml_quote(col.description) << "\n";
      }
   }

   out.close();

}

//------------------------------------------------------------------------------
// Parallel fill of a row-major double dump, then a serial write of the buffer.
//------------------------------------------------------------------------------
void output_binary_file(const std::string& filename, uint64_t n_rows, uint64_t n_cols,
                        binary_row_fn fill, void* ctx, const binary_dump_spec& spec){

   if(n_rows > 0 && fill == 0){
      std::cerr << "Error - binary dump of " << filename << " has no row filler." << std::endl;
      std::exit(EXIT_FAILURE);
   }

   std::ofstream ofile(filename.c_str(), std::ios::binary);
   if(!ofile.is_open()){
      std::cerr << "Error - cannot open binary file \"" << filename << "\" for writing." << std::endl;
      std::exit(EXIT_FAILURE);
   }

   ofile.write(reinterpret_cast<const char*>(&n_rows), sizeof(uint64_t));
   ofile.write(reinterpret_cast<const char*>(&n_cols), sizeof(uint64_t));

   if(n_rows > 0 && n_cols > 0){

      // Bound peak memory; OpenMP fills each chunk, then one write.
      const uint64_t chunk_bytes = 8ull * 1024ull * 1024ull;
      uint64_t rows_per_chunk = chunk_bytes / (n_cols * sizeof(double));
      if(rows_per_chunk < 1) rows_per_chunk = 1;

      std::vector<double> buf(static_cast<size_t>(rows_per_chunk * n_cols));

      for(uint64_t start = 0; start < n_rows; start += rows_per_chunk){

         const uint64_t n = std::min(rows_per_chunk, n_rows - start);
         const int nthr = (vdc::omp_threads > 0) ? vdc::omp_threads : 1;

         #pragma omp parallel num_threads(nthr)
         {
            #pragma omp for
            for(uint64_t i = 0; i < n; i++){
               fill(static_cast<size_t>(start + i), &buf[static_cast<size_t>(i * n_cols)], ctx);
            }
         }

         ofile.write(reinterpret_cast<const char*>(&buf[0]),
                     static_cast<std::streamsize>(n * n_cols * sizeof(double)));

         if(!ofile){
            std::cerr << "Error - failed writing binary file \"" << filename << "\"." << std::endl;
            std::exit(EXIT_FAILURE);
         }
      }
   }

   ofile.close();

   // upsert this kind and rewrite vdc-binary.meta
   bool found = false;
   for(size_t i = 0; i < g_binary_outputs.size(); i++){
      if(g_binary_outputs[i].spec.kind == spec.kind){
         g_binary_outputs[i].spec = spec;
         g_binary_outputs[i].last_file = filename;
         g_binary_outputs[i].pattern = filename_pattern(filename);
         g_binary_outputs[i].n_rows = n_rows;
         g_binary_outputs[i].n_cols = n_cols;
         found = true;
         break;
      }
   }
   if(!found){
      binary_output_record rec;
      rec.spec = spec;
      rec.last_file = filename;
      rec.pattern = filename_pattern(filename);
      rec.n_rows = n_rows;
      rec.n_cols = n_cols;
      g_binary_outputs.push_back(rec);
   }

   write_binary_metadata();

}

} // end of namespace vdc
