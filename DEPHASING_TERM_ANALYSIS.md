# Dephasing Term Analysis

## Problem
The spin torque fields blow up when the equilibration period ends, and the issue appears to originate from the dephasing term.

## Dephasing Term Calculation

The dephasing term is calculated as:
```cpp
deph = -sp / tau_phi
```

where:
- `sp` = perpendicular spin accumulation (S - (S·m̂)m̂)
- `tau_phi` = dephasing time = lambda_phi² / D

## Unit Conversion

From `interface.cpp` line 118:
- `lambda_phi` is read in Angstroms from the mat file
- Converted to meters: `lambda_phi = lphi * 1e-10`

## Values from CoCuCo-sv.mat

### Cobalt (Material 1 & 3) - UPDATED
- `spin-dephasing-length` = 32.0 Angstroms (updated from 21.0 to match PhysRevLett.109.117204)
- `diffusion-constant` = 0.0001 m²/s

**Calculated values (with updated lambda_phi):**
- `lambda_phi` = 32.0 × 1e-10 = 3.2e-9 m = 3.2 nm ✓ (from PhysRevLett.109.117204)
- `tau_phi` = (3.2e-9)² / 0.0001 = 1.024e-13 s = **102.4 fs** ✓ (improved from 44.1 fs)
- Dephasing rate = 1/tau_phi = 9.77e12 1/s = **9.77 1/ps** ✓ (improved from 22.68 1/ps)

**If |sp| = 1e7 C/m³ (typical spin accumulation):**
- |deph| = 1e7 / 1.024e-13 = **9.77e19 C/(m³·s)** ✓ (reduced by 2.3× from previous 2.27e20)

**Previous values (for reference):**
- Old `lambda_phi` = 21.0 Angstroms = 2.1 nm → `tau_phi` = 44.1 fs (VERY SHORT!)
- The update to 3.2 nm increases `tau_phi` by 2.3×, significantly improving numerical stability

### Copper (Material 2)
- `spin-dephasing-length` = 3500.0 Angstroms
- `diffusion-constant` = 0.001 m²/s

**Calculated values:**
- `lambda_phi` = 3500.0 × 1e-10 = 3.5e-7 m = 350 nm ✓ (reasonable)
- `tau_phi` = (3.5e-7)² / 0.001 = 1.225e-10 s = 122.5 ps ✓ (reasonable)
- Dephasing rate = 1/tau_phi = 8.16e9 1/s = 0.01 1/ps ✓ (reasonable)

**If |sp| = 1e7 C/m³:**
- |deph| = 1e7 / 1.225e-10 = 8.16e16 C/(m³·s) ✓ (reasonable)

## Root Cause

The problem is that for Cobalt:
1. **tau_phi = 44 fs is extremely short** - typical spin dephasing times are 100-500 fs or longer
2. **The dephasing rate is 22.68 1/ps** - this means any perpendicular spin accumulation will decay extremely rapidly
3. **If sp is not zero**, the dephasing term `deph = -sp / tau_phi` will be enormous and cause numerical instability

## Comparison with Literature and Other Material Files

### Serban_AOS_spinvalves.pdf
- Uses **De = 10⁻⁴ m²/s** for the spin valve system
- This **matches** the current Co value in `CoCuCo-sv.mat` (D = 0.0001 m²/s = 10⁻⁴ m²/s)

### IOP Science Paper (10.1088/0953-8984/10/48/005)
- Provides specific De values for Co and Cu
- **Note**: The actual values from this paper should be checked to verify if the current mat file values are physically reasonable

### Current Values in CoCuCo-sv.mat
- **Co**: D = 0.0001 m²/s = 10⁻⁴ m²/s ✓ (matches Serban)
- **Cu**: D = 0.001 m²/s = 10⁻³ m²/s (10× larger than Co)

### Physical Reasoning for Different De Values
- **Cu is a better conductor** than Co:
  - Cu resistivity ~ 1.7 μΩ·cm
  - Co resistivity ~ 6 μΩ·cm
  - Ratio ~ 3.5×
- **Better conductors have higher electron mobility** → higher diffusion constant
- **Expected ratio**: Cu should have D ~ 3-4× larger than Co
- **Current ratio**: 10× (might be too large, but direction is correct)

### Other Material Files in Codebase
- **FGaT**: D = 0.05714 m²/s (much larger than Co's 0.0001)
- **Mn2Au**: D = 0.01-0.002227 m²/s (much larger than Co's 0.0001)
- **Mn2Au-Jenkins**: D = 0.001 m²/s (same as current Cu value)

**Conclusion**: The Co value (10⁻⁴ m²/s) matches Serban's reference, so it may be physically correct. The issue is that this value, combined with `lambda_phi = 21 Angstroms`, creates a very short `tau_phi = 44 fs`, which is numerically challenging but may be physically accurate.

## Potential Issues

1. **Diffusion constant value**: D = 0.0001 m²/s for Co matches Serban's reference value
   - This suggests the value may be physically correct
   - However, it creates a very short tau_phi = 44 fs
   - The IOP Science paper (10.1088/0953-8984/10/48/005) should be consulted to verify if this is the correct value for Co
   - If the IOP paper gives a different value, that should be used instead

2. **Lambda_phi too small**: lambda_phi = 21 Angstroms = 2.1 nm
   - This is reasonable for transverse spin dephasing in ferromagnets
   - But combined with small D, it creates a very short tau_phi

3. **Numerical stability**: When the equilibration period ends:
   - Spin accumulation `S` may have built up
   - The perpendicular component `sp` may be non-zero
   - The dephasing term suddenly becomes active
   - If `sp` is large and `tau_phi` is very short, `deph` will be enormous and cause blow-up

## Recommendations

1. **Verify D values against IOP Science paper (10.1088/0953-8984/10/48/005)**
   - The paper provides specific De values for Co and Cu
   - Compare these to current values:
     - Co: D = 0.0001 m²/s (matches Serban)
     - Cu: D = 0.001 m²/s (10× larger than Co)
   - If the IOP paper gives different values, update the mat file accordingly
   - Note: The Co value matches Serban's reference, so it may be correct despite being numerically challenging

2. **lambda_phi updated to 32 Angstroms (3.2 nm)**
   - Updated from 21 Angstroms to match PhysRevLett.109.117204
   - This increases tau_phi from 44.1 fs to 102.4 fs (2.3× improvement)
   - Significantly improves numerical stability

3. **Add numerical stability check**:
   - If tau_phi < 1e-13 s (100 fs), issue a warning
   - Consider adding a minimum tau_phi floor (e.g., 100 fs) for numerical stability
   - Or add a check: if |deph| > threshold, reduce the time step or clamp the dephasing term

4. **Verify spin accumulation values**:
   - Check if `sp` (perpendicular component) is reasonable
   - If `sp` is unexpectedly large, investigate why

5. **Consider implicit treatment of dephasing**:
   - The dephasing term is currently treated explicitly in the Heun scheme
   - For very short tau_phi, an implicit treatment might be more stable

## Debug Output

The code now includes debug output that will print:
- `lambda_phi` (lphi[i]) in meters
- `D[i]` in m²/s
- `tau_phi` in seconds and femtoseconds
- `sp` (perpendicular spin accumulation) components
- `deph` (dephasing term) components
- Warnings if tau_phi is extremely short

This will help identify the exact values causing the blow-up.
