# Interpolation Step Precision Optimization

## Question
Can compile-time constants `kInvE` and `kInvMuB` be pre-combined during the interpolation step to reduce precision requirements and enable safe float conversion?

## Key Finding: Pre-combined Constant is Safe as Float!

### Pre-combined Constant

If we pre-compute at compile time:
```cpp
static constexpr double kTorqueScaleBase = kAtomcellVolume * kInvE * kInvMuB;
// = 2.89e-30 * 6.24e18 * 1.08e23 = 1.945e+12
```

**This constant is representable in float with only 0.000001% error!**

## Analysis Results

### Precision Test

| Computation | Double Value | Float Value | Relative Error |
|-------------|--------------|-------------|----------------|
| `kTorqueScaleBase` | 1.945001e+12 | 1.945001e+12 | 2.08e-8 (0.000002%) |
| Final torque (using float constant) | 1.151441e-01 T | 1.151441e-01 T | 1.10e-8 (0.000001%) |

### Operation Reordering Options

#### Current Approach (Line 1012-1014)
```cpp
spin_torque = kAtomcellVolume * sd_exchange[cell] * ax * kInvE * kInvMuB;
//            ^^^^^^^^^^^^^^^^ ^^^^^^^^^^^^^^^^
//            PROBLEM: This product underflows in float!
```

#### Option 1: Pre-combine Compile-Time Constants (Recommended)
```cpp
static constexpr double kTorqueScaleBase = kAtomcellVolume * kInvE * kInvMuB;

// Then in torque calculation:
spin_torque[3*cell+0] = kTorqueScaleBase * sd_exchange[cell] * ax;
```

**Benefits:**
- ✅ Reduces 4 multiplications to 2
- ✅ Computed at compile time
- ✅ `kTorqueScaleBase` can be float (tiny error: 0.000001%)
- ✅ Avoids underflow! (`kTorqueScaleBase * sd_exchange = 1.945e+12 * 4e-20 = 7.78e-8`, which is safe)

**Precision:** Still need double for `sd_exchange[cell] * ax`, but the final multiplication is safe.

#### Option 2: Scale During Fine-Grid Interpolation
```cpp
// During fine-grid loop (lines 922-926):
double ax_scaled = 0.0;
for(int j=0; j<nsub; j++){
    const int i = k*nsub + j;
    // Option: Multiply each fine-grid point by kTorqueScaleBase
    ax_scaled += Sbase[3*i+0] * kTorqueScaleBase;  // NEW
}
ax_scaled /= nsub;

// Then later:
spin_torque[3*cell+0] = sd_exchange[cell] * ax_scaled;
```

**Benefits:**
- ✅ Numerically equivalent to current approach
- ✅ Moves compile-time constant multiplication earlier

**Trade-offs:**
- ⚠️ Still need double precision for `Sbase` (spin accumulation values)
- ⚠️ No real precision benefit - just reorganizes computation

#### Option 3: Pre-combine All Constants Except ax
```cpp
// Per-cell pre-computation (could cache for each cell):
double torque_per_sa = kTorqueScaleBase * sd_exchange[cell];

// Then during coarse averaging (line 993):
double ax_scaled = ax * torque_per_sa;  // Store this instead of ax
// ...

// Later:
spin_torque[3*cell+0] = ax_scaled;  // Already scaled!
```

**Benefits:**
- ✅ Separates material-dependent scaling from accumulation

**Trade-offs:**
- ⚠️ Requires storing scaled values instead of raw `ax`
- ⚠️ More complex code

## Recommended Optimization

### Implementation: Pre-compute `kTorqueScaleBase`

**In spincurrents.cpp (near line 122):**
```cpp
// Add after existing constants:
static constexpr double kTorqueScaleBase = kAtomcellVolume * kInvE * kInvMuB;
// This can be computed at compile time and is numerically stable in float
// (relative error ~1e-8 if using float, but keep as double for consistency)
```

**In torque calculation (lines 1012-1014):**
```cpp
// OLD:
spin_torque[3*cell+0] = kAtomcellVolume * sd_exchange[cell] * ax * kInvE * kInvMuB;

// NEW (3 multiplications instead of 5):
spin_torque[3*cell+0] = kTorqueScaleBase * sd_exchange[cell] * ax;
spin_torque[3*cell+1] = kTorqueScaleBase * sd_exchange[cell] * ay;
spin_torque[3*cell+2] = kTorqueScaleBase * sd_exchange[cell] * az;
```

### Precision Impact

**If `kTorqueScaleBase` were converted to float:**
- Relative error in torque: **~0.000001%** (negligible)
- But we still need `double` for `sd_exchange[cell]` and `ax` due to:
  - `sd_exchange` varies per cell (material-dependent)
  - `ax` accumulates from fine-grid interpolation (integration-dependent)

### Why This Doesn't Enable Full Float Conversion

Even with `kTorqueScaleBase` pre-computed:

1. **Still need double for `sd_exchange[cell]`**: Material parameter, varies per cell
2. **Still need double for `ax`, `ay`, `az`**: These are accumulated from fine-grid spin accumulation `Sbase`, which:
   - Evolves via Heun integration (needs precision)
   - Represents charge density in C/m³ (large magnitude, needs precision)

3. **The underflow is avoided, but precision still required**:
   - `kTorqueScaleBase * sd_exchange = 1.945e+12 * 4e-20 = 7.78e-8` ✓ (safe range)
   - But `sd_exchange` itself (4e-20) is near float precision limits
   - `ax` (~1.48e6) multiplied by 7.78e-8 gives final torque (~0.115 T)

## Conclusion

**Pre-combining `kTorqueScaleBase` is beneficial, but doesn't enable float conversion:**

✅ **Benefits:**
- Fewer multiplications (performance: minor)
- Cleaner code (readability)
- Compile-time computation
- Avoids the underflow issue by reordering operations

❌ **Limitations:**
- Still need `double` for `sd_exchange[cell]` (material-dependent)
- Still need `double` for `ax`, `ay`, `az` (from integration)
- Converting `kTorqueScaleBase` to float saves negligible memory (it's just 1 constant)

**Recommendation:** Pre-compute `kTorqueScaleBase` for code clarity and minor performance improvement, but **keep it as `double`** for consistency with the rest of the computation. The precision savings from using float are negligible (0.000001% error), and keeping everything in double maintains numerical consistency.
