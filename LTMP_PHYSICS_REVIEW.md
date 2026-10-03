# Physics & Math Review: ltmp Module

## Standard Two-Temperature Model (TTM) Equations

The correct TTM equations are:
```
Ce * dTe/dt = -G*(Te - Tp) + κe*∇²Te + S
Cp * dTp/dt =  G*(Te - Tp) + κp*∇²Tp
```

Where:
- Ce, Cp = heat capacities (J/m³/K)
- G = electron-phonon coupling (J/s/m³/K)  
- κe, κp = thermal conductivities (J/s/m/K)
- S = laser source (W/m³ = J/s/m³)

---

## Issues Found

### ❌ **Issue 1: Laser Source Term Units (Line 78)**

**Current code:**
```cpp
const double pump = 1e10*ltmp::internal::pump_power*two_delta_sqrt_pi_ln_2*gaussian*i_pump_time/penetration_depth;
```

**Problem:**
- The factor `1e10` and division by `penetration_depth` are inconsistent with the reference implementation in `src/program/temperature_pulse.cpp` (line 148)
- Reference uses: `pump = two_delta_sqrt_pi_ln_2 * pump_power * gaussian * i_pump_time`
- The `penetration_depth` should be handled in the attenuation profile, not in the pump formula

**Expected formula (from reference):**
```cpp
const double pump = ltmp::internal::pump_power * two_delta_sqrt_pi_ln_2 * gaussian * i_pump_time;
```
Units: `(W/m³) * (dimensionless) * (dimensionless) * (1/s) = W/m³` ✓

**Fix:** Remove `1e10` and `/penetration_depth` factors.

---

### ❌ **Issue 2: Algorithm Structure - Diffusion Applied After Coupling**

**Current flow (lines 84-95, then 97-193):**
1. First loop: Update Te, Tp with coupling + laser (NO diffusion)
2. Second section: Calculate diffusion terms
3. Third loop: Apply diffusion deltas

**Problem:**
- The first loop (84-95) updates temperatures WITHOUT diffusion, then diffusion is added separately
- This creates a **split-step error** - the order matters and can cause instability
- Standard approach: calculate ALL terms (coupling + diffusion + source) together, then update

**Better approach:**
- Calculate diffusion terms FIRST (before any updates)
- Then update: `Te_new = Te + dt*(coupling + diffusion + source)/Ce`

---

### ❌ **Issue 3: Heat Diffusion Laplacian Implementation (Lines 113-114, 145-148)**

**Current code:**
```cpp
dTe_diff += (nTe - Te) * (κe[ncell] + κe[cell])*0.5 / (dr*1e-20);
```

**Problems:**

1. **Wrong Laplacian form**: The correct finite-difference Laplacian for 3D should be:
   ```
   ∇²T ≈ Σ_neighbors [ (T_neighbor - T_self) / dr² ]
   ```
   But the current code sums `κ*(T_neighbor - T_self)/dr²`, which gives `κ*∇²T` (correct for flux), BUT...

2. **Wrong κ averaging**: For flux continuity at interfaces, should use **harmonic mean**:
   ```
   κ_effective = 2*κ1*κ2/(κ1 + κ2)  // harmonic mean
   ```
   NOT arithmetic mean `(κ1 + κ2)/2`

3. **Missing cell volume normalization**: The Laplacian should account for the actual cell spacing. For a regular grid with spacing `Δx, Δy, Δz`, the 1D Laplacian is:
   ```
   d²T/dx² ≈ (T[i+1] - 2*T[i] + T[i-1]) / Δx²
   ```
   The current code uses `dr²` (3D distance squared), which is incorrect for a structured grid.

**Correct implementation should be:**
```cpp
// For each neighbor direction (x, y, z separately):
double dx = cell_spacing_x;  // actual grid spacing, not distance
dTe_diff += (nTe - Te) * κ_harmonic / (dx*dx*1e-20);
```

---

### ⚠️ **Issue 4: Electron Heat Capacity Temperature Dependence (Line 92)**

**Current code:**
```cpp
sqrt(Te + (G*(Tp-Te) + pump*attenuation)*dt/(Ce*Te))
```

**Analysis:**
- This assumes `Ce ∝ Te` (free electron model: Ce = γe * Te)
- The division by `Te` in the denominator is correct for this model
- However, the code uses `electron_heat_capacity[cell]` which is a constant (not temperature-dependent)
- **Inconsistency**: If using constant Ce, should be `dt/Ce`, not `dt/(Ce*Te)`

**Fix options:**
1. Use constant Ce: `Te + (G*(Tp-Te) + pump)*dt/Ce`
2. Use temperature-dependent: `Ce = Ce0 * Te`, then `dt/(Ce0*Te)` is correct

---

### ⚠️ **Issue 5: Phonon Heat Capacity - Debye Integral (Lines 336-351 in initialise.cpp)**

**Current implementation:**
```cpp
for(double T_D_over_T = 0.0; T_D_over_T < 12.0; T_D_over_T += 0.001) {
   int int_resolution = round((T_D_over_T) * 1000);
   double left = T_D_over_T + 0.00050;
   double right = T_D_over_T - 0.00050;
   double left_right_6 = (2*0.0005) / 8.0;
   integrand += left_right_6*(debeye_function(left) + 3.0*debeye_function(...) + ...);
   Debeye_phonon_constant[int_resolution] = 3.0*integrand/(T_D_over_T*T_D_over_T*T_D_over_T);
}
```

**Problems:**

1. **Integration bounds**: The integral should be from 0 to T_D/T, but the code accumulates from 0 to current T_D/T. This is correct for cumulative integration.

2. **Simpson's rule**: The `left_right_6` factor looks like Simpson's 3/8 rule, but the formula is wrong:
   - Correct Simpson's 3/8: `h/8 * (f0 + 3f1 + 3f2 + f3)` where `h = (b-a)/3`
   - Current: `(2*0.0005)/8.0 = 0.000125` but should be `h/8` where `h = 0.001`

3. **Debye function**: `debeye_function(x) = x⁴eˣ/(eˣ-1)²` is correct

4. **Final normalization**: `3.0*integrand/(T_D_over_T)³` should match the Debye formula, but the integral should be normalized differently.

**Correct Debye heat capacity:**
```
Cp(T) = 9*N*kB*(T/T_D)³ * ∫[0 to T_D/T] x⁴eˣ/(eˣ-1)² dx
```

The code's `3.0*integrand/(T_D_over_T)³` doesn't match this formula exactly.

---

### ❌ **Issue 6: Initialisation Averaging & Unused Normalization (Lines 232-237, 257, 265-270)**

**Current code:**
```cpp
// Line 257: Calculated but NEVER USED!
const double normalise_dz = (ltmp::internal::micro_cell_size[2]*ltmp::internal::micro_cell_size[2]*1e-20);

// Lines 232-237: Sum material properties
electron_heat_capacity[cell] += mp[mat].electron_heat_capacity;

// Lines 265-270: Average by atom count
electron_heat_capacity[cell] /= num_atoms_in_cell[cell];
```

**Problems:**

1. **Unused normalization factor**: `normalise_dz` is calculated (line 257) but never used. This suggests incomplete implementation.

2. **Averaging logic unclear**: 
   - Comments say properties are in **J/m³/K** (per-volume)
   - But code sums over atoms and divides by atom count
   - This gives: `(J/m³/K * N_atoms) / N_atoms = J/m³/K` ✓ (correct if material properties are per-volume)
   - OR: `(per-atom_value * N_atoms) / N_atoms = per-atom_value` (wrong units if we need per-volume)

3. **Missing volume normalization**: If material properties are per-atom, need to convert to per-volume:
   ```
   Ce_vol = Ce_atom * (atoms_per_cell / cell_volume)
   ```

**Fix needed:**
- Clarify if `mp[mat].electron_heat_capacity` is per-atom or per-volume
- If per-atom: multiply by atom density after averaging
- If per-volume: current averaging is wrong (should use volume-weighted average, not atom count)

---

### ⚠️ **Issue 7: Vertical Attenuation Coordinate (Line 368)**

**Current code:**
```cpp
double z = ltmp::internal::cell_position_array[3*cell+2];
vattn = exp(-z/ltmp::internal::penetration_depth);
```

**Analysis:**
- `z` is measured from the **bottom** (z=0 at bottom, z=system_dimensions_z at top)
- If laser enters from **top**, absorption should decrease with depth: `exp(-(system_dimensions_z - z)/penetration_depth)`
- If laser enters from **bottom**, current formula is correct: `exp(-z/penetration_depth)`
- **Need to verify convention**: Check if laser is assumed to enter from top or bottom

**Note:** The rutile branch uses the same formula, so this may be intentional if the convention is laser-from-bottom.

---

## Summary of Critical Issues

| Issue | Severity | Location | Impact |
|-------|----------|----------|--------|
| Laser source units | **HIGH** | Line 78 | Wrong heating magnitude (factor 1e10 error) |
| Diffusion Laplacian | **HIGH** | Lines 113-114, 145-148 | Incorrect heat flow (wrong κ averaging, wrong dr) |
| Algorithm order | **MEDIUM** | Lines 84-95 | Numerical instability (split-step error) |
| Ce temperature dep | **MEDIUM** | Line 92 | Inconsistent model (uses Ce*Te but Ce is constant) |
| Vertical attenuation | **LOW** | Line 368 | May be wrong if laser from top (needs convention check) |
| Debye integral | **LOW** | Lines 336-351 | Minor accuracy loss (Simpson's rule factor) |

## Additional Observations

### ✅ **Correct Implementations:**

1. **Einstein phonon model** (lines 30-33): Formula is correct for high T_D/T limit
2. **Phonon temperature projection** (lines 37-67): Uses predictor-corrector for temperature-dependent Cp - good approach
3. **Substrate cooling** (line 183): Correct form `-(Tp - Teq)*Tcool*dt`
4. **Root temperature storage**: Storing √T instead of T helps with numerical stability for low temperatures

### ⚠️ **Questionable but Possibly Intentional:**

1. **3D distance-based diffusion** (lines 112, 143): Using `dr² = dx²+dy²+dz²` instead of separate x,y,z Laplacians. This gives isotropic diffusion but may not match physical thermal conductivity anisotropy.
2. **Arithmetic mean κ** (lines 113, 145): Using `(κ1+κ2)/2` instead of harmonic mean. This is acceptable for small κ differences but wrong for large jumps (e.g., metal/insulator interfaces).

---

## Recommended Fixes Priority

1. **Fix laser source** (remove 1e10, remove /penetration_depth)
2. **Fix vertical attenuation** (use z_from_surface)
3. **Fix diffusion Laplacian** (use proper grid spacing, harmonic mean κ)
4. **Restructure algorithm** (calculate all terms before updating)
