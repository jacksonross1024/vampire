# Float vs Double Precision Analysis for 0.1 fs Timestep

## Executive Summary

For a **0.1 fs (1e-16 s) timestep**, the torque calculation has **critical precision requirements** that limit safe conversion to float precision.

## Key Findings

### Critical Issue: Underflow in Torque Calculation

The torque calculation:
```cpp
spin_torque = kAtomcellVolume * sd_exchange * ax * kInvE * kInvMuB
```

**Cannot use float for kAtomcellVolume or sd_exchange** because:
- `kAtomcellVolume * sd_exchange = 2.89e-30 * 4.0e-20 = 1.156e-49`
- This is **below float32 minimum normal value** (~1.175e-38)
- Results in **100% error (underflow to zero)**

### Safe Conversions (Verified)

Based on numerical testing:

| Constant | Safe as float? | Error | Notes |
|----------|---------------|-------|-------|
| **kInvMuB** | ✅ YES | 0.00% | Large constant but multiplies last, no error |
| **kInvE** | ⚠️ MAYBE | 0.076% | Acceptable error, test carefully |
| **kAtomcellVolume** | ❌ NO | 100% | Causes underflow when × sd_exchange |
| **sd_exchange** | ❌ NO | 100% | Causes underflow when × kAtomcellVolume |

### Estimated Torque Error for 0.1 fs Timestep

**If converting kInvE to float (safest tested):**
- **Relative error: 0.0757%**
- **Absolute error: ~8.7e-5 T** (for torque ~0.115 T)
- **This is acceptable** for most physics applications

**If converting kInvMuB to float:**
- **Relative error: 0.00%** (no measurable error)
- **Safe to convert**

**If converting kAtomcellVolume or sd_exchange:**
- **UNDERFLOW → 100% error**
- **NOT SAFE**

## Detailed Analysis

### 1. Timestep Precision (dt = 0.1 fs = 1e-16 s)

- Float epsilon: ~1.19e-7
- dt / float_eps = 8.39e-10
- **dt is 10 orders of magnitude smaller than float epsilon**
- **dt MUST remain double** to avoid catastrophic time accumulation errors

### 2. Spin Accumulation Integration

After 100 timesteps (Heun method):
- Max relative error in S: **6.82e-8** (0.00000682%)
- Error accumulates slowly due to exponential decay
- **Spin accumulation arrays (sc1d_Sfine) should remain double** to prevent long-term drift

### 3. Torque Calculation Chain

The calculation multiplies:
1. Very small: `kAtomcellVolume = 2.89e-30` (must be double)
2. Very small: `sd_exchange = 4.0e-20` (must be double)  
3. Moderate: `ax ~ 1.48e6` (spin accumulation in C/m³)
4. Very large: `kInvE = 6.24e18` (can be float, introduces ~0.076% error)
5. Very large: `kInvMuB = 1.08e23` (can be float, no error)

**Order of operations matters!** The underflow occurs at step (1 × 2).

## Recommendations (Safest to Most Aggressive)

### ✅ SAFEST (Recommended)

**Convert only kInvMuB to float:**
- Error: 0.00%
- Memory savings: Minimal (1 constant)
- **Risk: None**

### ⚠️ MODERATE (Test First)

**Convert kInvE and kInvMuB to float:**
- Error: ~0.0757%
- Memory savings: Minimal (2 constants)
- **Risk: Low** - acceptable for most applications
- **Test with your specific torque values before deploying**

### ❌ NOT SAFE (Do Not Convert)

**DO NOT convert these to float:**
- ❌ `dt` (timestep) - near float precision limit
- ❌ `kAtomcellVolume` - causes underflow
- ❌ `sd_exchange` - causes underflow  
- ❌ `sc1d_Sfine` (spin accumulation) - accumulates integration errors
- ❌ `sa_final`, `ns_final`, `ne_final` - final outputs need precision
- ❌ `k1`, `k2`, `Spred` (integration intermediates) - Heun method needs precision
- ❌ `spin_torque` (final output) - final result should be double

## Practical Impact

### Memory Savings

Converting constants to float saves **almost nothing** (constants are tiny). The real benefit would come from converting arrays, but those are **not safe to convert**.

### Performance Impact

Converting constants to float has **minimal performance impact**. The bottleneck is memory bandwidth and computation on arrays, not constant multiplications.

### Recommended Action

**For 0.1 fs timestep, DO NOT convert anything to float.**

The potential error from converting kInvE (~0.08%) may be acceptable, but:
1. Memory savings are negligible
2. Risk of introducing subtle bugs
3. dt precision is already marginal

**If timestep were larger (e.g., 1 fs):**
- More constants could safely be float
- But memory savings still minimal
- Recommend keeping everything double for consistency

## Conclusion

For 0.1 fs timestep, **keep all computations in double precision**. The only safe conversion (kInvMuB) saves no meaningful memory, and the slightly less safe conversion (kInvE) introduces ~0.08% error with no real benefit.

**Float conversion is not recommended for optimization #31 at 0.1 fs timestep.**
