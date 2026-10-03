# Plan: Migrate Thermal Gradient Calculation from ltmp to spin-currents Module

## Overview
This plan outlines the migration of laser-driven thermal gradient calculation from the `ltmp` module into the `spin-currents` module. The new implementation will use spin-currents' coarse-grained grid, maintain the link to sim-fields for thermal field updates, and route all interface options through the spin-torque module.

## Goals
1. Remove dependency on ltmp for thermal gradients when using spin-currents
2. Use spin-currents' existing coarse-grained microcell grid
3. Maintain compatibility with sim-fields thermal field calculation
4. Preserve all thermal physics (two-temperature model, Debye/Einstein heat capacity, etc.)
5. Keep all interface options in spin-torque module

## Current Architecture

### ltmp Module (to be partially replicated)
- **Files to reference:**
  - `src/ltmp/local_temperature_pulse.cpp` - Core two-temperature model solver
  - `src/ltmp/initialise.cpp` - Initialization with Debye model setup
  - `src/ltmp/field.cpp` - Thermal field calculation from temperatures
  - `src/ltmp/internal.hpp` - Data structures
  - `src/ltmp/data.cpp` - Material properties storage

### spin-currents Module (target location)
- **Current structure:**
  - `src/spintorque/spincurrents.cpp` - Main solver
  - `src/spintorque/initialise.cpp` - Initialization
  - `src/spintorque/internal.hpp` - Data structures
  - Uses coarse-grained grid: `num_microcells_per_stack` cells per stack

### sim-fields Integration
- **Current flow:**
  - `src/simulate/fields.cpp` line 173: `ltmp::get_localised_thermal_fields()` called
  - Thermal fields stored in `atoms::thermal_x_field`, `atoms::thermal_y_field`, `atoms::thermal_z_field`
  - Fields applied in LLG integrators (LLGMidpoint, LLGHeun, etc.)

## New File Structure

### New File: `src/spintorque/thermal_gradients.cpp`
**Purpose:** Core two-temperature model solver adapted for spin-currents coarse grid

**Key Functions:**
1. `calculate_thermal_gradients_1d_stack()` - Main solver for one stack
   - Input: stack index, time, dt
   - Output: Te and Tp arrays for coarse cells in stack
   - Based on `ltmp::internal::calculate_local_temperature_pulse()`

2. `einstein_model_phonon_heat_capacity()` - Debye/Einstein model
   - Copy from `ltmp/local_temperature_pulse.cpp` lines 33-36

3. `phonon_temperature_projector_step()` - Predictor-corrector for temperature-dependent Cp
   - Copy from `ltmp/local_temperature_pulse.cpp` lines 40-70
   - Two overloads: with and without substrate cooling

4. `update_thermal_gradients_all_stacks()` - Wrapper for all stacks
   - Loops over all local stacks
   - Calls `calculate_thermal_gradients_1d_stack()` for each

**Key Differences from ltmp:**
- Uses `num_microcells_per_stack` instead of `ltmp::internal::num_cells`
- Uses spin-currents' coarse grid z-positions (already calculated)
- Only vertical gradients (no lateral discretisation)
- Uses spin-currents' laser power density function (already integrated)

### New File: `src/spintorque/thermal_gradients_initialise.cpp`
**Purpose:** Initialize thermal gradient data structures using spin-currents grid

**Key Functions:**
1. `initialise_thermal_gradients()` - Main initialization
   - Called from `st::initialise()` after spin-currents grid is set up
   - Sets up per-coarse-cell material properties
   - Initializes temperature arrays
   - Sets up atom-to-cell mapping

2. `setup_thermal_material_properties()` - Extract material properties
   - Similar to `ltmp/initialise.cpp` lines 232-294
   - Maps material properties to coarse cells
   - Handles Debye constant lookup table initialization

3. `map_atoms_to_coarse_cells()` - Create atom-to-cell mapping
   - Maps each atom to its coarse microcell
   - Stores mapping for thermal field calculation
   - Similar to `ltmp/initialise.cpp` lines 180-242

**Data Structures (add to `internal.hpp`):**
```cpp
// Thermal gradient arrays (per coarse cell, per stack)
std::vector<double> sc1d_Te_coarse;  // Electron temperature (K)
std::vector<double> sc1d_Tp_coarse;  // Phonon temperature (K)
std::vector<double> sc1d_root_Te_coarse;  // sqrt(Te) for performance
std::vector<double> sc1d_root_Tp_coarse;  // sqrt(Tp) for performance

// Material properties (per coarse cell)
std::vector<double> sc1d_electron_heat_capacity;
std::vector<double> sc1d_phonon_heat_capacity;
std::vector<double> sc1d_electron_thermal_conductivity;
std::vector<double> sc1d_phonon_thermal_conductivity;
std::vector<double> sc1d_electron_phonon_coupling;
std::vector<double> sc1d_einstein_temperature;

// Debye lookup table (shared, initialized once)
std::vector<double> sc1d_debye_phonon_constant;

// Atom mapping (for thermal field calculation)
std::vector<int> sc1d_atom_to_coarse_cell;  // Maps atom index to coarse cell index
std::vector<int> sc1d_atom_temperature_type;  // 0=electron, 1=phonon (from material property)

// Thermal field arrays (per atom)
std::vector<double> sc1d_thermal_x_field;
std::vector<double> sc1d_thermal_y_field;
std::vector<double> sc1d_thermal_z_field;
```

### Modified File: `src/spintorque/spincurrents.cpp`
**Changes:**
1. Add call to `update_thermal_gradients_all_stacks()` in `calculate_spin_accumulation_1d()`
   - Call after charge transient update
   - Pass current time and dt

2. Modify `get_local_electron_temperature()` to use internal thermal gradients
   - When `sc1d_use_ltmp_temperatures = false` (new default)
   - Use `sc1d_Te_coarse` instead of `ltmp::get_electron_temperature_at_z()`

3. Remove dependency on ltmp for temperatures
   - Remove `#include "ltmp.hpp"` (or keep only for optional compatibility)
   - Remove calls to `ltmp::get_electron_temperature_at_z()`

### New File: `src/spintorque/thermal_fields.cpp`
**Purpose:** Calculate thermal fields from temperatures and apply to atoms

**Key Functions:**
1. `get_thermal_fields()` - Main interface (replaces `ltmp::get_localised_thermal_fields()`)
   - Input: start_index, end_index (atom range)
   - Output: Updates `atoms::thermal_x_field`, `atoms::thermal_y_field`, `atoms::thermal_z_field`
   - Based on `ltmp/field.cpp` lines 43-97

2. `calculate_thermal_field_for_atom()` - Per-atom calculation
   - Gets temperature from coarse cell
   - Applies temperature rescaling if enabled
   - Calculates `H_th_sigma * sqrt(T) * random_gaussian`

**Key Differences from ltmp:**
- Uses `sc1d_atom_to_coarse_cell` mapping instead of `ltmp::internal::atom_temperature_index`
- Uses `sc1d_root_Te_coarse` or `sc1d_root_Tp_coarse` based on material property
- Uses spin-currents' coarse cell temperatures

### Modified File: `src/spintorque/internal.hpp`
**Additions:**
- All data structures listed above
- Function declarations for thermal gradient functions
- Flags:
  - `sc1d_thermal_gradients_enable` - Enable thermal gradient calculation
  - `sc1d_thermal_gradients_initialised` - Initialization flag

### Modified File: `src/spintorque/interface.cpp`
**Changes:**
1. Add input parameter parsing for thermal gradient options:
   - `spin-torque:spin-currents-1d-thermal-gradients-enable`
   - Material properties (if not using ltmp defaults)
   - Substrate cooling constant

2. Add getter functions:
   - `get_electron_temperature(int coarse_cell)` - For external access
   - `get_phonon_temperature(int coarse_cell)` - For external access

### Modified File: `src/spintorque/initialise.cpp`
**Changes:**
1. Add call to `initialise_thermal_gradients()` in `st::initialise()`
   - After `initialise_spincurrents_1d()` is called
   - Only if `sc1d_thermal_gradients_enable == true`

2. Initialize Debye lookup table
   - Copy from `ltmp/initialise.cpp` lines 75-130
   - Store in `sc1d_debye_phonon_constant`

### Modified File: `src/simulate/fields.cpp`
**Changes:**
1. Replace `ltmp::get_localised_thermal_fields()` call (line 173)
   - With `st::get_thermal_fields()` when spin-currents thermal gradients enabled
   - Keep ltmp path for backward compatibility (when ltmp enabled but spin-currents thermal gradients disabled)

2. Logic:
```cpp
if(st::internal::sc1d_thermal_gradients_enable && st::internal::sc1d_thermal_gradients_initialised) {
   st::get_thermal_fields(atoms::thermal_x_field, atoms::thermal_y_field,
                          atoms::thermal_z_field, start_index, end_index);
} else if(ltmp::is_enabled()) {
   ltmp::get_localised_thermal_fields(atoms::thermal_x_field, atoms::thermal_y_field,
                                      atoms::thermal_z_field, start_index, end_index);
}
```

## Implementation Steps

### Phase 1: Core Thermal Gradient Solver
1. Create `thermal_gradients.cpp` with basic two-temperature model
   - Copy `calculate_local_temperature_pulse()` logic
   - Adapt to use coarse grid (only vertical, no lateral)
   - Use spin-currents' laser power density function
   - Test with simple case (no diffusion, only coupling)

2. Create `thermal_gradients_initialise.cpp`
   - Set up material properties per coarse cell
   - Initialize temperature arrays
   - Test initialization

### Phase 2: Integration with spin-currents
3. Modify `spincurrents.cpp`
   - Add call to update thermal gradients
   - Modify `get_local_electron_temperature()` to use internal temperatures
   - Test that temperatures are calculated correctly

4. Create `thermal_fields.cpp`
   - Implement thermal field calculation
   - Map atoms to coarse cells
   - Test field calculation

### Phase 3: sim-fields Integration
5. Modify `fields.cpp`
   - Add conditional call to `st::get_thermal_fields()`
   - Test that thermal fields are applied correctly in LLG

6. Modify `interface.cpp`
   - Add input parameter parsing
   - Add getter functions
   - Test parameter reading

### Phase 4: Material Properties and Debye Model
7. Add Debye/Einstein model support
   - Copy Debye constant lookup table initialization
   - Add `phonon_temperature_projector_step()` functions
   - Test temperature-dependent heat capacity

8. Add material property mapping
   - Extract properties from material file
   - Map to coarse cells
   - Handle empty cells (copy from neighbor)

### Phase 5: Testing and Validation
9. Compare with ltmp results
   - Run same simulation with both modules
   - Verify temperatures match (within numerical precision)
   - Verify thermal fields match

10. Remove ltmp dependency (optional)
    - Add flag to disable ltmp when using spin-currents thermal gradients
    - Test that simulation runs without ltmp

## Key Technical Details

### Grid Mapping
- **ltmp:** Uses 3D microcell grid (can be lateral + vertical)
- **spin-currents:** Uses 1D coarse grid per stack (vertical only)
- **Mapping:** Each atom maps to one coarse cell based on z-position

### Temperature Storage
- **ltmp:** `root_temperature_array[2*cell + 0/1]` stores sqrt(Te/Tp)
- **spin-currents:** `sc1d_root_Te_coarse[stack*num_cells + cell]` and `sc1d_root_Tp_coarse[...]`
- **Indexing:** Need stack index + cell index within stack

### Laser Source
- **ltmp:** Calculates laser power in `calculate_local_temperature_pulse()`
- **spin-currents:** Already has `laser_power_density()` function
- **Integration:** Use existing `laser_power_density()` in thermal gradient solver

### Diffusion Calculation
- **ltmp:** Uses neighbor list for 3D diffusion
- **spin-currents:** Only vertical diffusion (neighbors are cell k-1 and k+1)
- **Simplification:** Much simpler neighbor calculation

### Boundary Conditions
- **ltmp:** Handles 3D boundaries
- **spin-currents:** Only top/bottom boundaries (z=0 and z=z_max)
- **Top:** Laser source, zero flux boundary
- **Bottom:** Substrate cooling (sink term)

### MPI Considerations
- **ltmp:** Reduces material properties across MPI ranks
- **spin-currents:** Already handles MPI stack decomposition
- **Thermal gradients:** Each rank calculates for its local stacks
- **No MPI reduction needed:** Each stack is independent

## Output and Debugging

### Temperature Output
- Add option to output Te and Tp profiles (similar to spin-currents output)
- Use same output format as spin-currents (z-position, value)
- Output at same rate as spin-currents data

### Debugging Flags
- `sc1d_thermal_gradients_debug` - Print temperature values
- `sc1d_thermal_gradients_output_rate` - Control output frequency

## Backward Compatibility

### ltmp Module
- Keep ltmp module intact (no changes)
- Allow both modules to coexist
- Use flag to choose which module provides thermal gradients

### Input Parameters
- New parameters in spin-torque namespace
- Old ltmp parameters still work (if ltmp enabled)
- Clear documentation on which parameters to use

## Testing Checklist

- [ ] Thermal gradients initialize correctly
- [ ] Temperatures update correctly with laser pulse
- [ ] Diffusion works correctly (vertical only)
- [ ] Electron-phonon coupling works
- [ ] Substrate cooling works
- [ ] Thermal fields calculated correctly
- [ ] Thermal fields applied in LLG integrators
- [ ] Temperature-dependent heat capacity works
- [ ] Material property mapping works
- [ ] Empty cells handled correctly
- [ ] MPI parallelization works
- [ ] Output matches ltmp results (within tolerance)
- [ ] Performance is acceptable

## Files to Create/Modify Summary

### New Files:
1. `src/spintorque/thermal_gradients.cpp`
2. `src/spintorque/thermal_gradients_initialise.cpp`
3. `src/spintorque/thermal_fields.cpp`

### Modified Files:
1. `src/spintorque/spincurrents.cpp`
2. `src/spintorque/internal.hpp`
3. `src/spintorque/interface.cpp`
4. `src/spintorque/initialise.cpp`
5. `src/simulate/fields.cpp`
6. `makefile` (add new source files)

### Reference Files (no changes):
- `src/ltmp/*` - Use as reference for physics implementation

## Estimated Complexity

- **Core solver:** Medium (adapt existing code)
- **Initialization:** Medium (similar to ltmp but simpler grid)
- **Integration:** Low-Medium (well-defined interfaces)
- **Testing:** Medium (need to validate against ltmp)

## Notes

- This migration removes the timing mismatch issue between spin-currents and ltmp
- Simplifies the codebase by removing cross-module dependencies
- Uses spin-currents' existing grid infrastructure (no duplicate grid setup)
- Maintains all thermal physics from ltmp
- Keeps interface in spin-torque module (as requested)
