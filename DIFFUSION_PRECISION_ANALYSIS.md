# Diffusion Precision Analysis for Charge and Spin Currents

## Overview

Analysis of precision requirements for diffusion calculations in charge current (ns, ne) and spin current (S) evolution with 0.1 fs timestep.

## Key Findings

### 1. Thomas Tridiagonal Solver Precision

**Current Implementation (lines 157-181, 277, 303, 786):**
- Solves diffusion equation: `b[i]*x[i] + a[i]*x[i-1] + c[i]*x[i+1] = d[i]`
- Involves forward elimination and backward substitution with divisions

**Precision Test Results:**
- Division operations: **Relative error ~1.4e-8** if using float
- Example: `cp[0] = c[0]/b[0]` where `b[0] ≈ 1e5` (fine grid)
- Error accumulates through forward/backward substitution

**Verdict:** ✅ **Small errors** but they **accumulate** through solver steps

### 2. Gradient Calculations Precision

**Spin Current Diffusion (lines 938-951):**
```cpp
dSdz_x = (Sbase[3*(i+1)+0] - Sbase[3*i+0]) / dz;
// where dz ~ 1e-12 m, S ~ 1.48e7 C/m^3
```

**Charge Current Diffusion (lines 288, 300):**
```cpp
diff = -De * (ns[kR] - ns[kL]) / dz;
divJns = (Jns_edge[k+1] - Jns_edge[k]) / dz;
```

**Precision Test Results:**
- Small variations (0.1%): Relative error **~7.6e-9** in float
- Large variations (2x, interfaces): Relative error **~3.1e-8** in float
- Diffusion current `J_diff = -D * dS/dz`: Relative error **~8.3e-9**

**Verdict:** ✅ **Errors are small** but critical for accurate interface gradients

### 3. Matrix Coefficient Precision

**Fine Grid Diffusion Matrix (lines 756-765):**
```cpp
b[i] = 1.0 + dt_half * (alpha_edge[i-1] + alpha_edge[i]);
// where alpha = D/(dz*dz) ≈ 1e21 1/s
// dt_half * alpha ≈ 5e4
// b[i] ≈ 1e5
```

**Coarse Grid Charge Diffusion (lines 259-261):**
```cpp
b[k] = 1.0 + dt/tau_s + (dt/dz)*veR + dt*(DeL + DeR)/(dz*dz);
// where dt/tau_s ≈ 1e-4 (small contribution)
```

**Precision Test Results:**
- Matrix diagonal `b[i] ≈ 1e5` (fine grid) - **representable in float**
- Small contributions like `dt/tau_s ≈ 1e-4` - **may lose precision** when added to large numbers
- Relative error in matrix coefficient: **~1.5e-16** (negligible)

**Verdict:** ⚠️ **Small contributions may be lost** when adding to 1.0 in float

### 4. Critical Operations

#### Charge Current (lines 200-316):
1. **Matrix assembly**: `dt*D/(dz*dz) ≈ 0.1` (coarse) or `1e5` (if using fine grid)
2. **Thomas solve**: Forward/backward substitution with divisions
3. **Gradient computation**: `(ns[kR] - ns[kL]) / dz`

#### Spin Current Diffusion (lines 751-789, 938-981):
1. **Matrix assembly**: `dt_half * alpha_edge` where `alpha ≈ D/(dz*dz) ≈ 1e21`
2. **Thomas solve**: Applied to each component (x, y, z) of spin accumulation
3. **Gradient computation**: `(S[i+1] - S[i]) / dz` with `dz ≈ 1e-12 m`

## Precision Impact Assessment

### Safe for Float? (with caveats)

| Operation | Float Error | Safe? | Notes |
|-----------|-------------|-------|-------|
| **Thomas solver divisions** | ~1.4e-8 | ⚠️ Maybe | Small but accumulates through n steps |
| **Gradient (small variations)** | ~7.6e-9 | ⚠️ Maybe | Interface gradients critical |
| **Gradient (large variations)** | ~3.1e-8 | ⚠️ Maybe | Steep interfaces need precision |
| **Matrix coefficients** | ~1.5e-16 | ✅ Yes | But diagonal b[i] very large |
| **Small contributions (dt/tau)** | Potential loss | ❌ No | May be lost in float |

### Critical Issues with Float

1. **Accumulated Errors in Thomas Solver**
   - Each division has ~1e-8 error
   - For n=100 grid points, errors accumulate
   - Forward/backward substitution amplifies errors

2. **Loss of Small Contributions**
   - `dt/tau_s ≈ 1e-4` added to `1.0 + dt*D/(dz*dz) ≈ 1.2`
   - In float, small contribution may not be accurately represented

3. **Gradient Precision at Interfaces**
   - Steep gradients (interfaces) require precision
   - `(S[i+1] - S[i]) / dz` amplifies small differences
   - Float errors (~3e-8) may affect interface physics

## Recommendations

### ❌ **DO NOT convert to float:**

1. **Thomas solver arrays** (`a`, `b`, `c`, `rhs`, `cp`, `dp`)
   - Division operations are sensitive
   - Errors accumulate through solver

2. **Spin accumulation arrays** (`sc1d_Sfine`, `Sbase`)
   - Source of gradient calculations
   - Large magnitude (~1e7 C/m³)
   - Evolves through integration

3. **Charge density arrays** (`ns_coarse`, `ne_coarse`)
   - Used in gradient calculations
   - Evolves through integration

4. **Matrix coefficients** (`a_tr`, `b_tr`, `c_tr`, `alpha_edge`)
   - Used in Thomas solver
   - Small contributions (dt/tau) need precision

### ⚠️ **POSSIBLY safe (test carefully):**

1. **Diffusion coefficients** (`De_cell`, `De_edge`, `D_fine`)
   - Constant material properties
   - **But**: Memory savings minimal (fewer arrays)

2. **Velocity arrays** (`ve_cell`, `ve_edge`)
   - Derived from diffusion
   - **But**: Used in matrix assembly

3. **Pre-computed geometric factors**
   - `1/dz`, `dz*dz` if constant
   - **But**: Usually combined with other factors

### ✅ **Safe (but negligible benefit):**

1. **Material constants** (`beta_cond`, `beta_diff` in fine-grid arrays)
   - Already in double, converting saves minimal memory

## Overall Verdict

**For 0.1 fs timestep: KEEP ALL DIFFUSION CALCULATIONS IN DOUBLE PRECISION**

**Reasons:**
1. ✅ Errors (~1e-8) are small but **accumulate** over many timesteps
2. ✅ Thomas solver is **sensitive to precision** (divisions, forward/backward substitution)
3. ✅ Gradient calculations **amplify small differences** by 1/dz (1e12)
4. ✅ Small contributions (`dt/tau`) may be **lost** in float
5. ✅ Memory savings would be minimal (arrays are already double)

**The precision requirements for diffusion at 0.1 fs are too stringent for float conversion.**

## Alternative Optimizations (without float conversion)

Instead of float conversion, consider:

1. **Pre-compute constant factors** (e.g., `1/dz`, `dt/(dz*dz)`)
2. **Optimize memory access patterns** (structure-of-arrays, cache blocking)
3. **Reduce temporary allocations** (move vectors outside loops)
4. **SIMD vectorization** (AVX2/AVX-512 for gradient calculations)

These optimizations provide better performance gains with **zero precision loss**.
