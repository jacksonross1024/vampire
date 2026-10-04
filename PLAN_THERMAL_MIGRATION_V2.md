# Migration Plan: Thermal Gradients from ltmp to spin-currents

## Objective
Move the thermal gradient calculation functionality from the `ltmp` module into the `spin-currents` module, using the existing spin-currents coarse-grained grid structure. This eliminates cross-module dependencies and timing synchronization issues.

## Current State Analysis

### ltmp Module Components
1. **Two-temperature model solver** (`local_temperature_pulse.cpp`)
   - Solves: Ce*dTe/dt = -G*(Te-Tp) + κe*∇²Te + S
   - Solves: Cp*dTp/dt = G*(Te-Tp) + κp*∇²Tp
   - Uses 3D microcell grid (can be lateral + vertical)
   - Stores temperatures as sqrt(T) for numerical stability

2. **Initialization** (`initialise.cpp`)
   - Sets up microcell grid
   - Maps atoms to cells
   - Extracts material properties per cell
   - Initializes Debye/Einstein heat capacity lookup tables

3. **Thermal field calculation** (`field.cpp`)
   - Converts temperatures to thermal fields
   - Applies to atoms via `get_localised_thermal_fields()`

### spin-currents Module Structure
- Uses 1D coarse grid: `num_microcells_per_stack` cells per stack
- Grid spacing: `micro_cell_thickness` (Angstroms)
- Already has laser power density calculation
- Already has temperature getter functions (currently calls ltmp)

### Integration Points
- `src/simulate/fields.cpp` line 173: Calls `ltmp::get_localised_thermal_fields()`
- LLG integrators use `atoms::thermal_x_field`, `atoms::thermal_y_field`, `atoms::thermal_z_field`

## Implementation Strategy

### File: `src/spintorque/thermal_solver.cpp`
**Purpose:** Core two-temperature model solver for spin-currents coarse grid

**Main Function:**
```cpp
void update_thermal_gradients_stack(int stack_idx, double time_s, double dt_si)
```
- Input: stack index, current time, time step
- Updates Te and Tp for all coarse cells in the stack
- Algorithm:
  1. Calculate diffusion terms (vertical only, neighbors are k-1 and k+1)
  2. Calculate coupling terms (electron-phonon exchange)
  3. Calculate laser source term (use existing `laser_power_density()`)
  4. Update temperatures: Te_new = Te_old + dt*(coupling + diffusion + source)/Ce
  5. Update temperatures: Tp_new = Tp_old + dt*(coupling + diffusion)/Cp(Tp)
  6. Store as sqrt(T) for stability

**Helper Functions:**
- `einstein_model_phonon_heat_capacity()` - Debye model calculation
- `phonon_temperature_projector_step()` - Predictor-corrector for temperature-dependent Cp
- `calculate_diffusion_terms()` - Vertical diffusion only (simpler than ltmp's 3D)

**Key Simplifications:**
- Only vertical diffusion (no lateral gradients)
- Neighbor calculation: cell k neighbors are k-1 and k+1
- Grid spacing: `micro_cell_thickness * 1e-10` (convert Angstroms to meters)

### File: `src/spintorque/thermal_init.cpp`
**Purpose:** Initialize thermal gradient data structures

**Main Function:**
```cpp
void initialise_thermal_gradients()
```
- Called from `st::initialise()` after spin-currents grid is ready
- Steps:
  1. Allocate temperature arrays (Te, Tp, sqrt(Te), sqrt(Tp)) per coarse cell
  2. Extract material properties from material file
  3. Map properties to coarse cells (average over atoms in cell)
  4. Initialize Debye constant lookup table
  5. Create atom-to-coarse-cell mapping
  6. Initialize temperatures to equilibrium value

**Data Structures (add to `internal.hpp`):**
```cpp
// Temperature arrays: [stack_idx * num_cells + cell_idx]
std::vector<double> sc1d_Te;           // Electron temperature (K)
std::vector<double> sc1d_Tp;          // Phonon temperature (K)  
std::vector<double> sc1d_sqrt_Te;     // sqrt(Te) for stability
std::vector<double> sc1d_sqrt_Tp;     // sqrt(Tp) for stability

// Material properties: [cell_idx] (same for all stacks, materials vary by z)
std::vector<double> sc1d_Ce;          // Electron heat capacity (J/m³/K)
std::vector<double> sc1d_Cp;          // Phonon heat capacity (J/m³/K)
std::vector<double> sc1d_kappa_e;     // Electron thermal conductivity (J/s/m/K)
std::vector<double> sc1d_kappa_p;     // Phonon thermal conductivity (J/s/m/K)
std::vector<double> sc1d_G;           // Electron-phonon coupling (J/s/m³/K)
std::vector<double> sc1d_T_Debye;     // Debye temperature (K)

// Debye lookup table (shared, size ~24000)
std::vector<double> sc1d_debye_table;

// Atom mapping
std::vector<int> sc1d_atom_cell_idx;      // Maps atom -> coarse cell index
std::vector<bool> sc1d_atom_use_phonon;   // true if atom couples to phonon temp
```

**Material Property Extraction:**
- Loop over all atoms
- For each atom, determine which coarse cell it belongs to (based on z-coordinate)
- Accumulate material properties for that cell
- After loop, average properties per cell (divide by number of atoms)
- Handle empty cells: copy from previous cell

### File: `src/spintorque/thermal_fields.cpp`
**Purpose:** Calculate thermal fields from temperatures and apply to atoms

**Main Function:**
```cpp
void get_thermal_fields(std::vector<double>& thermal_x,
                        std::vector<double>& thermal_y,
                        std::vector<double>& thermal_z,
                        int start_idx, int end_idx)
```
- Input: field arrays (to be filled), atom range
- Algorithm:
  1. Generate random Gaussian numbers for each atom
  2. For each atom:
     - Get coarse cell index from `sc1d_atom_cell_idx[atom]`
     - Get temperature: if `sc1d_atom_use_phonon[atom]` use Tp, else use Te
     - Get material: `atoms::type_array[atom]`
     - Get H_th_sigma from material properties
     - Calculate: `thermal_field = H_th_sigma * sqrt(T) * random_gaussian`
  3. Store in `atoms::thermal_x_field`, etc.

**Temperature Rescaling:**
- If material has `temperature_rescaling_Tc > 0`, apply rescaling
- Formula: `T_rescaled = Tc * (T/Tc)^alpha` if T < Tc, else T
- Use sqrt(T) for efficiency

### Modifications to Existing Files

#### `src/spintorque/spincurrents.cpp`
1. Add call to update thermal gradients:
   ```cpp
   void calculate_spin_accumulation_1d() {
      // ... existing code ...
      
      // Update thermal gradients if enabled
      if(sc1d_thermal_gradients_enable) {
         for(int stack : sc1d_local_stacks) {
            update_thermal_gradients_stack(stack, t_s, dt_si);
         }
      }
   }
   ```

2. Modify `get_local_electron_temperature()`:
   ```cpp
   double get_local_electron_temperature(double z_m) {
      if(sc1d_thermal_gradients_enable) {
         // Use internal thermal gradients
         int cell = get_coarse_cell_index(z_m);
         int stack = get_stack_index(z_m);  // or use current stack context
         int idx = stack * num_microcells_per_stack + cell;
         return sc1d_Te[idx];
      } else if(sc1d_use_ltmp_temperatures && ltmp::is_enabled()) {
         // Fallback to ltmp
         return ltmp::get_electron_temperature_at_z(z_m * 1e10);
      } else {
         return sc1d_reference_temperature;
      }
   }
   ```

#### `src/spintorque/initialise.cpp`
1. Add initialization call:
   ```cpp
   void st::initialise(...) {
      // ... existing initialization ...
      
      if(sc1d_thermal_gradients_enable) {
         initialise_thermal_gradients();
      }
   }
   ```

2. Initialize Debye table:
   ```cpp
   void initialise_debye_table() {
      // Copy from ltmp/initialise.cpp lines 75-130
      // Calculate Debye constant for T_D/T from 0 to 24.0
      sc1d_debye_table.resize(24001);
      for(int i = 0; i <= 24000; ++i) {
         double T_D_over_T = i / 1000.0;
         // ... Debye calculation ...
      }
   }
   ```

#### `src/spintorque/interface.cpp`
1. Add input parameter parsing:
   ```cpp
   bool match_input_parameter(...) {
      // Add:
      // "spin-torque:spin-currents-1d-thermal-gradients-enable"
      // Material property overrides (optional)
   }
   ```

#### `src/simulate/fields.cpp`
1. Modify thermal field calculation:
   ```cpp
   if(program::program == 13 || program::program == 55) {
      if(st::internal::sc1d_thermal_gradients_enable && 
         st::internal::sc1d_thermal_gradients_initialised) {
         // Use spin-currents thermal gradients
         st::get_thermal_fields(atoms::thermal_x_field,
                                atoms::thermal_y_field,
                                atoms::thermal_z_field,
                                start_index, end_index);
      } else if(ltmp::is_enabled()) {
         // Fallback to ltmp
         ltmp::get_localised_thermal_fields(atoms::thermal_x_field,
                                            atoms::thermal_y_field,
                                            atoms::thermal_z_field,
                                            start_index, end_index);
      }
   }
   ```

#### `src/spintorque/internal.hpp`
1. Add data structure declarations (see above)
2. Add function declarations:
   ```cpp
   void update_thermal_gradients_stack(int stack_idx, double time_s, double dt_si);
   void initialise_thermal_gradients();
   void get_thermal_fields(std::vector<double>& thermal_x,
                          std::vector<double>& thermal_y,
                          std::vector<double>& thermal_z,
                          int start_idx, int end_idx);
   ```

3. Add flags:
   ```cpp
   extern bool sc1d_thermal_gradients_enable;
   extern bool sc1d_thermal_gradients_initialised;
   ```

## Implementation Details

### Grid Mapping
- **Coarse cell index from z:** `cell = floor(z_angstrom / micro_cell_thickness)`
- **Stack index:** Already determined by spin-currents (y-position)
- **Array indexing:** `idx = stack_idx * num_microcells_per_stack + cell_idx`

### Diffusion Calculation
- **Vertical only:** Neighbors are cell k-1 and k+1
- **Boundary conditions:**
  - Top (k = num_cells-1): Zero flux (laser source handled separately)
  - Bottom (k = 0): Substrate cooling (sink term)
- **Harmonic mean for conductivity:** At interface between cells with different materials
- **Grid spacing:** `dz = micro_cell_thickness * 1e-10` (meters)

### Laser Source Integration
- Use existing `laser_power_density(z_m, t_s, step_counter, dt_si)` function
- Call for each coarse cell center: `z_cell = (cell + 0.5) * micro_cell_thickness * 1e-10`
- Source term: `S = laser_power_density(z_cell, ...)` (W/m³)
- Only applied to electron temperature equation

### Temperature-Dependent Heat Capacity
- **Electrons:** Linear: `Ce(T) = Ce0 * T` (free electron model)
- **Phonons:** Debye/Einstein model using lookup table
- **Predictor-corrector:** For phonons, use two-step method to handle T-dependent Cp

### Atom-to-Cell Mapping
- During initialization, loop over all atoms
- For each atom, calculate which coarse cell it belongs to
- Store mapping: `sc1d_atom_cell_idx[atom] = cell_idx`
- Also store whether atom couples to electron or phonon temperature

## Testing Strategy

1. **Unit Tests:**
   - Test diffusion calculation (compare with analytical solution)
   - Test coupling terms (verify energy conservation)
   - Test Debye model (compare with ltmp values)

2. **Integration Tests:**
   - Compare temperatures with ltmp for same input
   - Verify thermal fields match ltmp
   - Check that LLG integrators receive correct fields

3. **Performance Tests:**
   - Measure computation time vs ltmp
   - Verify MPI scaling

## Migration Steps

1. **Phase 1:** Create `thermal_solver.cpp` with basic two-temperature model (no diffusion)
2. **Phase 2:** Add diffusion terms (vertical only)
3. **Phase 3:** Add initialization code (`thermal_init.cpp`)
4. **Phase 4:** Add thermal field calculation (`thermal_fields.cpp`)
5. **Phase 5:** Integrate with spin-currents solver
6. **Phase 6:** Integrate with sim-fields
7. **Phase 7:** Add Debye/Einstein model
8. **Phase 8:** Testing and validation
9. **Phase 9:** Remove ltmp dependency (optional)

## Backward Compatibility

- Keep ltmp module unchanged
- Add flag to choose between ltmp and spin-currents thermal gradients
- Default: use spin-currents if enabled, otherwise fallback to ltmp
- Input parameters: new parameters in spin-torque namespace

## Notes

- This eliminates the timing synchronization issue between spin-currents and ltmp
- Uses existing spin-currents grid (no duplicate grid setup)
- Simpler than ltmp (no lateral gradients, only vertical)
- All physics preserved from ltmp
- Interface through spin-torque module only
