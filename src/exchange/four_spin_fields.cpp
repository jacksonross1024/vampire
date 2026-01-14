//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2023. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//------------------------------------------------------------------------------
//
// NOTE: Four-spin exchange internal data structures not available in current
//       implementation. This is a stub to allow compilation.

// C++ standard library headers
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "exchange.hpp"

// exchange module headers
#include "internal.hpp"

namespace exchange{
namespace internal{

void four_spin_exchange_fields(const int start_index,
                               const int end_index,
                               std::vector<double>& field_array_x,
                               std::vector<double>& field_array_y,
                               std::vector<double>& field_array_z){
   // Stub - four-spin exchange not implemented in this branch
   return;
}

} // end of internal namespace
} // end of exchange namespace
