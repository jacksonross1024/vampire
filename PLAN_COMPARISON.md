# Comparison of Two Migration Plans

## Document 1: `PLAN_THERMAL_GRADIENTS_MIGRATION.md`
## Document 2: `PLAN_THERMAL_MIGRATION_V2.md`

## Major Structural Differences

### File Organization

**Document 1 proposes:**
- `thermal_gradients.cpp` - Core solver
- `thermal_gradients_initialise.cpp` - Initialization
- `thermal_fields.cpp` - Field calculation

**Document 2 proposes:**
- `thermal_solver.cpp` - Core solver
- `thermal_init.cpp` - Initialization
- `thermal_fields.cpp` - Field calculation

**Analysis:** Different naming conventions, but same functional breakdown. No contradiction.

### Data Structure Naming

**Document 1:**
- `sc1d_Te_coarse`, `sc1d_Tp_coarse`
- `sc1d_root_Te_coarse`, `sc1d_root_Tp_coarse`
- `sc1d_atom_to_coarse_cell`
- `sc1d_atom_temperature_type` (0=electron, 1=phonon)

**Document 2:**
- `sc1d_Te`, `sc1d_Tp`
- `sc1d_sqrt_Te`, `sc1d_sqrt_Tp`
- `sc1d_atom_cell_idx`
- `sc1d_atom_use_phonon` (bool)

**Analysis:** Different naming, but functionally equivalent. Document 2's bool is clearer than Document 1's int. Minor contradiction in naming style.

### Array Indexing Strategy

**Document 1:**
- `sc1d_Te_coarse[stack*num_cells + cell]`
- Explicitly mentions stack index in indexing

**Document 2:**
- `sc1d_Te[idx]` where `idx = stack_idx * num_microcells_per_stack + cell_idx`
- Same indexing, but Document 2 is more explicit about the calculation

**Analysis:** Same approach, Document 2 is more detailed. No contradiction.

### Material Property Storage

**Document 1:**
- Material properties stored per coarse cell
- Properties: `sc1d_electron_heat_capacity`, `sc1d_phonon_heat_capacity`, etc.
- No explicit mention of whether properties are per-stack or shared

**Document 2:**
- Material properties stored per coarse cell
- Properties: `sc1d_Ce`, `sc1d_Cp`, `sc1d_kappa_e`, `sc1d_kappa_p`, `sc1d_G`, `sc1d_T_Debye`
- Explicitly states: "same for all stacks, materials vary by z"

**Analysis:** Document 2 clarifies that properties are shared across stacks (since materials vary by z, not by x,y). Document 1 is ambiguous. Minor contradiction in clarity.

### Function Signatures

**Document 1:**
- `calculate_thermal_gradients_1d_stack(int start_cell, ...)` - mentions start_cell
- `update_thermal_gradients_all_stacks()` - wrapper function

**Document 2:**
- `update_thermal_gradients_stack(int stack_idx, double time_s, double dt_si)` - uses stack_idx
- No wrapper function mentioned

**Analysis:** Document 1's `start_cell` parameter is unclear (what does it mean?). Document 2's `stack_idx` is clearer. Document 1 has a wrapper function that Document 2 doesn't mention. Minor contradiction.

### Integration with spincurrents.cpp

**Document 1:**
- Calls `update_thermal_gradients_all_stacks()` in `calculate_spin_accumulation_1d()`
- Mentions calling after charge transient update

**Document 2:**
- Calls `update_thermal_gradients_stack()` in loop over `sc1d_local_stacks`
- Same location in `calculate_spin_accumulation_1d()`

**Analysis:** Document 1 uses wrapper function, Document 2 uses direct loop. Functionally equivalent, but Document 2 is more explicit. No contradiction.

### Temperature Getter Function

**Document 1:**
- `get_local_electron_temperature()` modification mentioned but not detailed

**Document 2:**
- Provides full implementation with:
  - Check for `sc1d_thermal_gradients_enable`
  - Fallback to ltmp if enabled
  - Fallback to reference temperature
  - Explicit cell and stack index calculation

**Analysis:** Document 2 is more complete. No contradiction, just different levels of detail.

### Debye Table Initialization

**Document 1:**
- Mentions copying from `ltmp/initialise.cpp` lines 75-130
- Stored in `sc1d_debye_phonon_constant`

**Document 2:**
- Mentions copying from `ltmp/initialise.cpp` lines 75-130
- Stored in `sc1d_debye_table`
- Explicitly states size: 24001 elements (for T_D/T from 0 to 24.0)

**Analysis:** Different naming (`sc1d_debye_phonon_constant` vs `sc1d_debye_table`), but Document 2 provides more detail. Minor contradiction in naming.

### Atom Mapping Details

**Document 1:**
- `sc1d_atom_to_coarse_cell` - maps atom to coarse cell
- `sc1d_atom_temperature_type` - 0=electron, 1=phonon

**Document 2:**
- `sc1d_atom_cell_idx` - maps atom to coarse cell
- `sc1d_atom_use_phonon` - bool, true if couples to phonon

**Analysis:** Same functionality, different naming. Document 2's bool is clearer. Minor contradiction.

### Grid Mapping Details

**Document 1:**
- Mentions using spin-currents' coarse grid
- No explicit formula for cell index calculation

**Document 2:**
- Provides explicit formula: `cell = floor(z_angstrom / micro_cell_thickness)`
- Explains stack index is already determined by spin-currents

**Analysis:** Document 2 is more explicit. No contradiction, just different levels of detail.

### Diffusion Calculation

**Document 1:**
- Mentions vertical diffusion only
- Mentions neighbors are k-1 and k+1
- Mentions harmonic mean for conductivity

**Document 2:**
- Same approach
- Adds explicit boundary conditions:
  - Top: Zero flux (laser source handled separately)
  - Bottom: Substrate cooling (sink term)
- Adds explicit grid spacing: `dz = micro_cell_thickness * 1e-10`

**Analysis:** Document 2 is more complete with boundary conditions. No contradiction.

### Laser Source Integration

**Document 1:**
- Mentions using existing `laser_power_density()` function
- No details on how to call it

**Document 2:**
- Explicitly states: call for each coarse cell center
- Provides formula: `z_cell = (cell + 0.5) * micro_cell_thickness * 1e-10`
- States source term is W/m³
- States only applied to electron temperature equation

**Analysis:** Document 2 is more detailed. No contradiction.

### Temperature-Dependent Heat Capacity

**Document 1:**
- Mentions Debye/Einstein model
- Mentions predictor-corrector for phonons

**Document 2:**
- Explicitly states: `Ce(T) = Ce0 * T` for electrons (free electron model)
- Explicitly states: Debye/Einstein model for phonons using lookup table
- Mentions predictor-corrector for phonons

**Analysis:** Document 2 provides more physics detail. No contradiction.

### sim-fields Integration

**Document 1:**
- Provides code snippet with conditional logic
- Checks `sc1d_thermal_gradients_enable && sc1d_thermal_gradients_initialised`
- Falls back to ltmp if not enabled

**Document 2:**
- Same approach
- Same conditional logic
- Same fallback

**Analysis:** Identical approach. No contradiction.

### Implementation Phases

**Document 1:**
- 9 phases, very detailed breakdown
- Phase 1: Core solver (no diffusion)
- Phase 2: Add diffusion
- Phase 3: Initialization
- Phase 4: Field calculation
- Phase 5: Integration with spin-currents
- Phase 6: Integration with sim-fields
- Phase 7: Material properties and Debye model
- Phase 8: Testing
- Phase 9: Remove ltmp dependency

**Document 2:**
- 9 phases, similar breakdown
- Same order and content

**Analysis:** Identical phases. No contradiction.

## Summary of Contradictions

### Minor Contradictions (Naming/Clarity):

1. **File naming:** `thermal_gradients.cpp` vs `thermal_solver.cpp` - Different names, same purpose
2. **Data structure naming:** Various differences (`_coarse` suffix, `sqrt_` vs `root_`, etc.)
3. **Function parameter:** `start_cell` vs `stack_idx` - Document 1's parameter is unclear
4. **Debye table name:** `sc1d_debye_phonon_constant` vs `sc1d_debye_table`
5. **Atom mapping type:** `int` (0/1) vs `bool` - Document 2's bool is clearer

### No Major Contradictions:

- Both plans agree on:
  - Core two-temperature model physics
  - Vertical-only diffusion
  - Use of spin-currents coarse grid
  - Integration points with sim-fields
  - Backward compatibility approach
  - Implementation phases
  - Testing strategy

## Recommendations

1. **Use Document 2's naming:** More consistent and clearer (bool vs int, explicit names)
2. **Use Document 2's detail level:** More explicit formulas and implementation details
3. **Keep Document 1's wrapper function:** `update_thermal_gradients_all_stacks()` is useful abstraction
4. **Combine best of both:** Document 2's clarity with Document 1's organizational structure

## Conclusion

The two documents are **fundamentally consistent** with no major contradictions. Differences are primarily:
- Naming conventions (minor)
- Level of implementation detail (Document 2 is more detailed)
- Some organizational choices (wrapper functions, etc.)

Both plans would lead to the same functional result, with Document 2 providing more implementation guidance.
