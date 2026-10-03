# Analysis: Implementing Superdiffusive Transport and Seebeck Effect

## Executive Summary

The current implementation in `spincurrents.cpp` has a **well-structured foundation** that would make adding superdiffusive transport and Seebeck effect **moderately straightforward**. The code already has:
- ✅ Temperature infrastructure (`get_local_electron_temperature`, `Te_cell`)
- ✅ Charge current calculation framework (`Jc_edge`, `Jns_edge`)
- ✅ Spin current calculation framework (drift + diffusion)
- ✅ Material parameter infrastructure (`beta_cond`, `beta_diff`, `diffusion`)

**Estimated Implementation Difficulty**: 
- **Seebeck Effect**: ⭐⭐ (Easy-Moderate) - ~2-3 days
- **Superdiffusive Transport**: ⭐⭐⭐ (Moderate) - ~4-5 days

---

## Current Implementation Analysis

### 1. Charge Transient Solver (Lines 194-327)

**Current Structure:**
```cpp
// Equations solved:
// (1) d ns/dt = -d Jns/dz - ns/tau_s + laser_source
// (2) Jns = ve*ns - De*d(ns)/dz  (non-equilibrium current)
// (3) d ne/dt = -d Jc/dz
// (4) Jc = -De*d(ne)/dz + Jns  (total charge current)
```

**Key Components:**
- `ns`: Non-equilibrium charge density (hot electrons)
- `ne`: Excess charge density (equilibrium deviation)
- `Jns`: Non-equilibrium current = `ve*ns - De*d(ns)/dz`
- `Jc`: Total charge current = `-De*d(ne)/dz + Jns`
- `ve = De/dopt`: Simplified superdiffusive velocity (line 236)

**Current Limitations (Lines 201-203):**
```cpp
// Simplifications (per latest implementation plan):
//  - Seebeck term and explicit Te-gradient terms are not included yet.
//  - Temperature-gradient driven superdiffusive term is folded into ve = De/d.
```

### 2. Spin Accumulation Solver (Lines 827-906)

**Current Structure:**
```cpp
// Spin current components:
// Js = J_drift + J_diff
// J_drift = -beta_cond * Jc * mhat  (line 942-944)
// J_diff = -D * dS/dz  (line 990-992)
```

**Key Components:**
- Drift term uses `beta_cond` (spin polarization of conductivity)
- Diffusion term uses `D` (diffusion constant)
- Both are already material-dependent and interpolated to fine grid

---

## Implementation Requirements

### A. Seebeck Effect Implementation

#### Physics (from Serban AOS paper):
The Seebeck effect generates charge currents from temperature gradients. There are two cases:

**Case 1: Spin-independent Seebeck (S↑ = S↓)**
- **Charge current**: `J_seebeck = -S * dTe/dz` where `S` is Seebeck coefficient
- **Spin current**: Automatically included via drift term `J_drift = -beta * Jc * mhat`
  - When `Jc` includes Seebeck contribution, the spin current gets it through `beta * Jc * mhat`
  - **No separate spin Seebeck term needed** in this case

**Case 2: Spin-dependent Seebeck (S↑ ≠ S↓)**
- **Charge current**: `J_seebeck = -S_avg * dTe/dz` where `S_avg = (S↑ + S↓)/2`
- **Spin current**: Additional term `Js_seebeck = -S_s * dTe/dz * mhat` where `S_s = (S↑ - S↓)/2`
  - This creates a **pure spin current** (no net charge current) when S↑ ≠ S↓
  - This term is **independent** of the charge current contribution

**Key Insight**: If the Seebeck coefficient is the same for spin-up and spin-down electrons, adding Seebeck to the charge current is sufficient. The spin-dependent term is only needed when there's a difference between S↑ and S↓.

#### Code Changes Required:

**1. Add Material Parameters** (in `initialise.cpp` or material files):
```cpp
// New material properties needed:
extern std::vector<double> seebeck_coefficient;      // S (V/K)
extern std::vector<double> seebeck_spin_coefficient;  // S_s (V/K) - spin-dependent
```

**2. Modify Charge Current Calculation** (Lines 317-326):
```cpp
// Current: Jc_edge[e] = diff + Jns_edge[e];
// Add Seebeck term:
const double dTe_dz = (Te_cell[kR] - Te_cell[kL]) / dz;
const double S_avg = 0.5 * (seebeck_coefficient[cellL] + seebeck_coefficient[cellR]);
const double J_seebeck = -S_avg * dTe_dz;  // A/m^2
Jc_edge[e] = diff + Jns_edge[e] + J_seebeck;  // Add Seebeck contribution
```

**3. Modify Spin Current Calculation** (Lines 940-1001):
```cpp
// OPTIONAL: Only needed if S↑ ≠ S↓ (spin-dependent Seebeck)
// If S↑ = S↓, the Seebeck contribution to spin current comes automatically
// through the drift term: jdr = -beta * Jc * mhat (where Jc includes Seebeck)

// For spin-dependent Seebeck (S↑ ≠ S↓), add separate term:
if(sc1d_spin_dependent_seebeck_enable) {
   const double dTe_dz = (Te_fine[i+1] - Te_fine[i-1]) / (2.0*dz);  // central difference
   const double S_s = seebeck_spin_coefficient[i];  // S_s = (S↑ - S↓)/2
   const double js_seebeck_x = -S_s * dTe_dz * mhat[3*i+0];
   const double js_seebeck_y = -S_s * dTe_dz * mhat[3*i+1];
   const double js_seebeck_z = -S_s * dTe_dz * mhat[3*i+2];
   
   // Add to total spin current:
   jsx_sum += (jdr_x + jdiff_x + js_seebeck_x);
} else {
   // Standard case: Seebeck already included via Jc in drift term
   jsx_sum += (jdr_x + jdiff_x);
}
```

**4. Temperature Gradient Calculation**:
- ✅ Already have `Te_cell` array (line 222)
- ✅ Already have `get_local_electron_temperature()` function (line 91)
- Need to interpolate `Te` to fine grid (similar to material properties)

**Implementation Difficulty**: ⭐⭐ (Easy-Moderate)
- **Pros**: Temperature infrastructure exists, straightforward gradient calculation
- **Cons**: Need new material parameters, need to interpolate Te to fine grid
- **Lines to modify**: ~30-50 lines for charge current only, ~50-70 lines if including spin-dependent term

**Recommendation**: Start with charge current Seebeck only. The spin current will automatically get the Seebeck contribution through the existing drift term (`J_drift = -beta * Jc * mhat`). Only add the separate spin Seebeck term if materials have significantly different S↑ and S↓.

---

### B. Superdiffusive Transport Implementation

#### Physics (from Serban AOS paper):
Superdiffusive transport involves hot electrons with:
- **Time-dependent velocity**: `v_e(t) = v_e0 * exp(-t/tau_e)` where `tau_e` is energy relaxation time
- **Temperature-dependent velocity**: `v_e ~ sqrt(Te)` or `v_e ~ Te^alpha`
- **Energy-dependent mean free path**: Longer mean free path for hot electrons

#### Current Simplification (Line 236):
```cpp
ve_cell[k] = De / dopt;  // Simplified: constant velocity
```

#### Code Changes Required:

**1. Enhanced Velocity Model** (Lines 220-237):
```cpp
// Current: ve_cell[k] = De / dopt;
// Enhanced: Add time and temperature dependence

// Option A: Temperature-dependent velocity
const double Te_ratio = Te_cell[k] / sc1d_reference_temperature;
const double v_e0 = std::sqrt(De / dopt);  // Base velocity
const double v_e = v_e0 * std::sqrt(Te_ratio);  // v ~ sqrt(Te)

// Option B: Time-dependent decay (for hot electron relaxation)
const double t_hot = t_s - t_laser_peak;  // Time since laser peak
const double tau_e = 100e-15;  // Energy relaxation time (~100 fs)
const double v_e = v_e0 * std::exp(-t_hot / tau_e);

// Option C: Combined model
const double v_e_base = v_e0 * std::sqrt(Te_ratio);
const double v_e = v_e_base * std::exp(-t_hot / tau_e);
ve_cell[k] = v_e;
```

**2. Enhanced Non-Equilibrium Current** (Lines 291-301):
```cpp
// Current: Jns_edge[e] = adv + diff;
// Enhanced: Add energy-dependent mean free path

// Hot electron mean free path (longer for higher energy)
const double lambda_e = lambda_sdl[cell] * std::sqrt(Te_ratio);
const double De_hot = v_e * lambda_e;  // Hot electron diffusion

// Use hot electron properties for superdiffusive transport
const double adv = v_e * ns[kL];  // Use enhanced velocity
const double diff = -De_hot * (ns[kR] - ns[kL]) / dz;  // Use hot electron diffusion
Jns_edge[e] = adv + diff;
```

**3. Separate Hot/Cold Electron Populations** (Advanced):
If implementing full two-population model:
- Add `ns_hot` and `ns_cold` arrays
- Add coupling terms between populations
- Requires more extensive changes (~200+ lines)

**Implementation Difficulty**: ⭐⭐⭐ (Moderate)
- **Pros**: Velocity calculation already exists, just needs enhancement
- **Cons**: Need to decide on model complexity, may need new material parameters
- **Lines to modify**: ~30-50 lines for simple model, ~200+ for full two-population model

---

## Specific Code Locations for Changes

### Seebeck Effect:

1. **Material Parameter Initialization** (`initialise.cpp`, ~line 500):
   - Add `seebeck_coefficient` and `seebeck_spin_coefficient` to material loop
   - Add MPI reduction for these arrays

2. **Charge Current Calculation** (`spincurrents.cpp`, lines 317-326):
   - Calculate `dTe/dz` at edges
   - Add Seebeck term to `Jc_edge`

3. **Spin Current Calculation** (`spincurrents.cpp`, lines 940-1001):
   - Interpolate `Te` to fine grid (similar to material properties)
   - Calculate `dTe/dz` on fine grid
   - Add Seebeck spin current term

4. **Fine Grid Interpolation** (`spincurrents.cpp`, lines 416-463):
   - Add `Te_fine` array initialization
   - Interpolate `Te` similar to other material properties

### Superdiffusive Transport:

1. **Velocity Calculation** (`spincurrents.cpp`, lines 220-237):
   - Enhance `ve_cell` calculation with temperature/time dependence

2. **Non-Equilibrium Current** (`spincurrents.cpp`, lines 291-301):
   - Use enhanced velocity and hot electron diffusion

3. **Material Parameters** (if needed):
   - Add `tau_e` (energy relaxation time)
   - Add `v_e0` (base hot electron velocity)

---

## Recommended Implementation Order

1. **Phase 1: Seebeck Effect** (Easier, independent)
   - Add material parameters
   - Implement charge current Seebeck term
   - Implement spin current Seebeck term
   - Test and validate

2. **Phase 2: Enhanced Superdiffusive Transport** (Moderate)
   - Start with temperature-dependent velocity
   - Add time-dependent decay
   - Test and compare with simplified model

3. **Phase 3: Full Two-Population Model** (Advanced, optional)
   - Only if Phase 2 shows limitations
   - Requires more extensive refactoring

---

## Compatibility with Existing Code

✅ **Fully Compatible**: Both additions can be implemented as **optional enhancements**:
- Use flags: `sc1d_seebeck_enable`, `sc1d_superdiffusive_enhanced`
- Default to current behavior if flags are off
- No breaking changes to existing functionality

✅ **Infrastructure Ready**:
- Temperature arrays already exist
- Material parameter system is extensible
- Fine-grid interpolation framework is in place

---

## Estimated Code Changes

| Component | Lines Added | Lines Modified | Complexity |
|-----------|-------------|----------------|------------|
| Seebeck Material Params | ~15 | ~5 | Low |
| Seebeck Charge Current | ~15 | ~5 | Low |
| Seebeck Spin Current (optional) | ~20 | ~5 | Low |
| Superdiffusive Velocity | ~30 | ~10 | Medium |
| Superdiffusive Current | ~20 | ~10 | Medium |
| **Total (charge only)** | **~30** | **~10** | **Low** |
| **Total (with spin term)** | **~50** | **~15** | **Low-Medium** |

---

## Conclusion

The current codebase is **well-positioned** for adding both superdiffusive transport and Seebeck effect. The modular structure, existing temperature infrastructure, and material parameter system make these additions straightforward. The main challenges are:

1. **Seebeck**: Adding new material parameters and gradient calculations (easy)
2. **Superdiffusive**: Choosing appropriate model complexity (moderate)

Both can be implemented incrementally without disrupting existing functionality.
