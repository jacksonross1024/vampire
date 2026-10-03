# Superdiffusive Transport Model Analysis
## Based on References 28, 34, 35 from Serban AOS Spin Valves Paper

## Current Implementation Analysis

### Current Charge Current Structure (Lines 194-327)

The current implementation has:
```cpp
// Eq. (1) equivalent: Jc = J_diffusive + J_superdiffusive (simplified)
// 
// Jns = ve*ns - De*d(ns)/dz  (non-equilibrium current, line 298-300)
// Jc = -De*d(ne)/dz + Jns    (total charge current, line 325)
//
// where ve = De/dopt (line 236) - SIMPLIFIED superdiffusive velocity
```

**Current Simplification** (Line 203):
- "Temperature-gradient driven superdiffusive term is folded into ve = De/d"

---

## What References 28, 34, 35 Tell Us About Superdiffusion

### Reference 28: Battiato, Carva, Oppeneer (Phys. Rev. B 86, 024404, 2012)

**Key Physics:**
- **Superdiffusive transport** involves hot electrons with **energy-dependent mean free paths**
- Hot electrons (excited above Fermi level) have **longer mean free paths** than thermal electrons
- Transport is **ballistic/quasi-ballistic** over short distances before scattering
- **Time-dependent velocity**: Hot electrons lose energy and slow down as they thermalize
- **Energy-dependent scattering**: Higher energy electrons scatter less frequently

**Model Characteristics:**
- Energy-resolved electron distribution: `f(E, z, t)`
- Energy-dependent velocity: `v(E) ~ sqrt(E - E_F)`
- Energy-dependent mean free path: `λ(E) = v(E) * τ(E)` where `τ(E)` is energy-dependent lifetime
- Time decay: Hot electrons thermalize with characteristic time `τ_e` (energy relaxation time)

**Equation Form:**
```
J_superdiffusive = ∫ v(E) * f_hot(E, z, t) * dE
```
where `f_hot(E, z, t)` is the non-equilibrium hot electron distribution.

---

### Reference 34: Najafi et al. (Nat. Commun. 8, 15177, 2017)

**Key Physics:**
- **Two-population model**: Hot electrons (non-thermal) vs. thermal electrons
- **Forward superdiffusive flow**: Initial ballistic transport of hot electrons
- **Backward diffusive flow**: Re-equilibration after thermalization
- **Asymmetric amplitudes**: Forward and backward flows have different magnitudes
- **Pulse-width dependence**: Different pulse widths lead to different forward/backward ratios

**Model Characteristics:**
- Separate populations: `n_hot` and `n_thermal`
- Hot electrons: High velocity, long mean free path, short lifetime
- Thermal electrons: Lower velocity, shorter mean free path, longer lifetime
- Coupling between populations: Hot electrons thermalize into thermal population

---

### Reference 35: Melnikov et al. (Phys. Rev. Lett. 107, 076601, 2011)

**Key Physics:**
- **Ultrafast demagnetization** driven by superdiffusive hot electron transport
- **Spin-dependent transport**: Different mean free paths for spin-up vs. spin-down
- **Energy window**: Only electrons in specific energy range contribute significantly
- **Interface effects**: Reflection/transmission at material interfaces

**Model Characteristics:**
- Spin-dependent mean free paths: `λ↑(E) ≠ λ↓(E)`
- Energy window: `E_F < E < E_F + ΔE` (typically ~1-2 eV above Fermi level)
- Interface transmission: `T(E)` depends on energy and material

---

## What Eq. 1 in Serban Paper Likely Contains

Based on the references and typical superdiffusion models, Eq. 1 probably has the form:

```
J_total = J_superdiffusive + J_diffusive
```

Where:

### Term 1: Superdiffusive Current
```
J_superdiffusive = ∫ v(E) * f_hot(E) * T(E) * dE
```
- `v(E)`: Energy-dependent velocity
- `f_hot(E)`: Hot electron distribution (non-equilibrium, laser-excited)
- `T(E)`: Transmission through layers/interfaces
- Integration over energy window above Fermi level

### Term 2: Diffusive Current  
```
J_diffusive = -σ * ∇μ - S * ∇T + D * ∇n
```
- Standard drift-diffusion terms
- Includes Seebeck effect (`-S * ∇T`)
- Charge/spin diffusion (`D * ∇n`)

---

## Boundary Conditions for Superdiffusive Transport

### Critical Boundary Condition Considerations (from Serban AOS Spin Valves, Page 5)

The first paragraph of page 5 in Serban's paper discusses important boundary condition considerations for superdiffusive transport in spin valves. Key points include:

#### 1. Interface Current Matching (Charge and Spin)

**Requirement**: At interfaces between layers (e.g., free layer/spacer/reference layer), both **charge current** and **spin current** must be continuous, accounting for both superdiffusive and diffusive components.

**Physical Basis**: 
- The forward superdiffusive flow of hot electrons carries both charge and spin
- The backward diffusive re-equilibration flow also carries charge and spin
- Interface transmission/reflection affects both components differently

**Implementation Implications**:
```cpp
// At interface between layers:
// Jc_superdiffusive + Jc_diffusive must be continuous
// Js_superdiffusive + Js_diffusive must be continuous
// 
// Current code (lines 239-249) handles edge-averaged properties,
// but may need enhancement for superdiffusive component
```

#### 2. Energy-Dependent Interface Transmission

**Requirement**: Hot electrons (superdiffusive) may have **different transmission probabilities** at interfaces compared to thermal electrons (diffusive).

**Physical Basis**:
- High-energy hot electrons can overcome interface barriers more easily
- Energy-dependent transmission: `T(E)` for hot electrons vs. `T_thermal` for diffusive
- Spin-dependent transmission: `T↑(E) ≠ T↓(E)` for spin-polarized interfaces

**Implementation Implications**:
- If implementing energy-resolved model (Option 3), need `T(E)` at each interface
- For simplified models (Option 1-2), may use effective transmission coefficient
- Current code uses `r_int_edge` for interface resistance, but this is for diffusive transport

#### 3. Surface Boundary Conditions

**Top Surface (Laser Side)**:
- **Source boundary**: Laser injects hot electrons, creating non-equilibrium distribution
- Boundary condition: `J_hot_in = laser_source(t)` for superdiffusive component
- May also have reflection: `J_hot_out = R * J_hot_in` where `R` is surface reflectivity

**Bottom Surface (Substrate Side)**:
- Typically **zero-flux** or **absorbing** boundary for spin current
- May allow charge current to flow into substrate
- Current code enforces `Jc_edge[0] = 0` and `Jc_edge[ncz] = 0` (lines 318-319)

**Implementation Implications**:
```cpp
// Current (lines 318-319):
Jc_edge[0] = 0.0;      // Bottom: zero charge current
Jc_edge[ncz] = 0.0;    // Top: zero charge current

// For superdiffusion, may need:
// Top: Jc_edge[ncz] = J_laser_source - R * J_reflected
// Bottom: Jc_edge[0] = -T_substrate * J_incident (if absorbing)
```

#### 4. Spin-Polarized Interface Reflections

**Requirement**: Interfaces must account for **spin-dependent reflection/transmission** of superdiffusive hot electrons.

**Physical Basis**:
- Minority vs. majority spin electrons have different scattering at interfaces
- This creates spin accumulation (key mechanism for switching)
- Forward superdiffusive flow: polarized by free layer
- Backward diffusive flow: repolarized by reference layer

**Implementation Implications**:
- Need spin-resolved transmission: `T↑` and `T↓` at each interface
- Current code tracks spin accumulation `S` but may need interface spin filtering
- The `r_int_edge` parameter (line 473) handles interface resistance but may need spin-resolved version

#### 5. Time-Dependent Boundary Behavior

**Requirement**: Boundary conditions may need to be **time-dependent** to capture the transition from superdiffusive to diffusive phases.

**Physical Basis**:
- Initial phase (0-100 fs): Superdiffusive hot electrons dominate, high transmission
- Later phase (100+ fs): Thermalized electrons, diffusive transport, different boundary behavior
- Boundary conditions should reflect this temporal evolution

**Implementation Implications**:
- May need time-dependent transmission coefficients
- Or separate boundary conditions for `ns` (hot) vs. `ne` (thermal) populations
- Current code uses constant `ve` and `De`, but time-dependent versions would help

#### 6. Nonlocal Boundary Effects

**Requirement**: Superdiffusive transport is **nonlocal** - carriers from deep in the material can reach boundaries, affecting boundary conditions.

**Physical Basis**:
- Hot electrons have long mean free paths (tens of nm)
- Can travel from one layer to another before scattering
- Boundary conditions must account for carriers arriving from non-adjacent regions

**Implementation Implications**:
- Current 1D model (z-direction only) should handle this, but need to ensure:
  - Integration over mean free path accounts for nonlocal contributions
  - Boundary flux includes contributions from entire penetration depth
- For energy-resolved model, integration over energy also includes nonlocal effects

#### 7. Charge/Spin Accumulation at Boundaries

**Requirement**: Boundaries must allow for **accumulation** of charge and spin, especially minority spins at the reference layer.

**Physical Basis**:
- Forward superdiffusive flow accumulates minority spins at reference layer
- This accumulation drives switching (P→AP transition)
- Boundary conditions that immediately sink spin accumulation would suppress this mechanism

**Implementation Implications**:
- Current code allows accumulation (no Dirichlet BC forcing `S=0`)
- But need to ensure superdiffusive component contributes to accumulation correctly
- The `ns` and `ne` arrays (lines 214-215) track charge, but spin accumulation `S` is separate

---

## Comparison with Current Implementation

### What's Currently Implemented:

**Lines 220-237: Velocity Calculation**
```cpp
ve_cell[k] = De / dopt;  // Simplified: constant velocity
```
- ❌ Missing: Energy dependence
- ❌ Missing: Time-dependent decay
- ❌ Missing: Temperature dependence

**Lines 291-301: Non-equilibrium Current**
```cpp
Jns_edge[e] = ve * ns[kL] - De * (ns[kR] - ns[kL]) / dz;
```
- ✅ Has: Advection term (`ve * ns`)
- ✅ Has: Diffusion term (`-De * d(ns)/dz`)
- ❌ Missing: Energy integration
- ❌ Missing: Separate hot/thermal populations

**Lines 317-326: Total Charge Current**
```cpp
Jc_edge[e] = -De * d(ne)/dz + Jns_edge[e];
```
- ✅ Has: Diffusive term for `ne`
- ✅ Has: Superdiffusive-like term via `Jns`
- ❌ Missing: Explicit energy-dependent superdiffusion
- ❌ Missing: Seebeck term (noted as TODO on line 202)

---

## Recommended Implementation Strategy

### Option 1: Enhanced Single-Population Model (Easier, ~50 lines)

**Enhance the existing `ve` calculation** with time and temperature dependence:

```cpp
// Enhanced velocity (lines 220-237)
const double Te_ratio = Te_cell[k] / sc1d_reference_temperature;
const double v_e0 = std::sqrt(De / dopt);  // Base velocity

// Temperature-dependent velocity (from energy distribution)
const double v_e_temp = v_e0 * std::sqrt(Te_ratio);

// Time-dependent decay (hot electrons thermalize)
const double t_hot = t_s - t_laser_peak;
const double tau_e = 100e-15;  // Energy relaxation time (~100 fs)
const double v_e_time = v_e0 * std::exp(-t_hot / tau_e);

// Combined model
ve_cell[k] = v_e_temp * (1.0 + (v_e_time - v_e0) / v_e0);
```

**Pros:**
- Minimal code changes
- Captures main physics (time decay, temperature dependence)
- Compatible with existing structure

**Cons:**
- Doesn't capture full energy dependence
- Single population (no separate hot/thermal)

---

### Option 2: Two-Population Model (More Complete, ~200+ lines)

**Add separate hot and thermal electron populations:**

```cpp
// New arrays needed:
std::vector<double> ns_hot;   // Hot electron density
std::vector<double> ns_thermal; // Thermal electron density
std::vector<double> ve_hot;    // Hot electron velocity (energy-dependent)
std::vector<double> ve_thermal; // Thermal electron velocity

// Separate evolution equations:
// d ns_hot/dt = laser_source - ns_hot/tau_e - coupling
// d ns_thermal/dt = coupling - ns_thermal/tau_s + diffusion
```

**Pros:**
- Captures full physics from references
- Can model forward/backward asymmetry
- Energy-dependent transport

**Cons:**
- Significant code restructuring
- More parameters needed
- More complex

---

### Option 3: Energy-Integrated Model (Most Complete, ~300+ lines)

**Full energy integration as in Battiato et al. (Ref 28):**

```cpp
// Energy bins
const int n_energy_bins = 20;
const double dE = 2.0 / n_energy_bins;  // 0 to 2 eV above Fermi

// For each energy bin:
for(int ie=0; ie<n_energy_bins; ie++){
   const double E = E_F + (ie + 0.5) * dE;
   const double v_E = std::sqrt(2.0 * E / m_e);  // Energy-dependent velocity
   const double lambda_E = v_E * tau_E;  // Energy-dependent mean free path
   const double f_hot_E = laser_excitation(E, t);  // Hot electron distribution
   
   J_superdiff += v_E * f_hot_E * transmission(E);
}
```

**Pros:**
- Most physically accurate
- Captures all energy-dependent effects

**Cons:**
- Very complex
- Requires energy-resolved distribution functions
- Computational overhead

---

## Recommended Approach

**Start with Option 1** (Enhanced Single-Population):
1. Add time-dependent velocity decay
2. Add temperature-dependent velocity scaling
3. Test and validate
4. If needed, upgrade to Option 2

**Key Parameters to Add:**
- `tau_e`: Energy relaxation time (~50-200 fs)
- `v_e0`: Base hot electron velocity
- `E_hot`: Characteristic hot electron energy
- `t_laser_peak`: Time of laser peak (for decay calculation)

**Code Locations:**
- Lines 220-237: Enhance `ve_cell` calculation
- Lines 291-301: Use enhanced velocity in `Jns`
- May need to track `t_laser_peak` from laser pulse timing

---

## Integration with Seebeck Effect

The **Seebeck effect** should be added **separately** to the diffusive term:

```cpp
// Current: Jc = -De*d(ne)/dz + Jns
// Enhanced: Jc = -De*d(ne)/dz + Jns + J_seebeck

const double dTe_dz = (Te_cell[kR] - Te_cell[kL]) / dz;
const double J_seebeck = -S * dTe_dz;
Jc_edge[e] = diff + Jns_edge[e] + J_seebeck;
```

The superdiffusive term (`Jns`) and Seebeck term (`J_seebeck`) are **independent**:
- **Superdiffusive**: Non-equilibrium hot electrons (fast, ~100 fs timescale)
- **Seebeck**: Thermal gradient-driven (slower, ~ps timescale)

---

## Summary

**Current State:**
- ✅ Basic structure in place (`Jns` term)
- ✅ Boundary conditions: Zero-flux at surfaces (`Jc_edge[0] = Jc_edge[ncz] = 0`)
- ✅ Interface handling: `r_int_edge` for interface resistance
- ❌ Simplified velocity model (`ve = De/dopt`)
- ❌ Missing time/temperature dependence
- ❌ Missing energy integration
- ❌ Missing energy-dependent interface transmission for superdiffusive component
- ❌ Missing spin-resolved boundary conditions

**Boundary Condition Gaps:**
- Current boundaries enforce zero charge current, but may need:
  - Laser source injection at top surface for superdiffusive hot electrons
  - Energy-dependent transmission at interfaces
  - Spin-resolved reflection/transmission
  - Time-dependent boundary behavior (superdiffusive → diffusive transition)

**Recommended Next Steps:**
1. Enhance `ve` with time decay and temperature scaling (Option 1)
2. Add Seebeck term to charge current
3. **Review and enhance boundary conditions** for superdiffusive transport:
   - Add laser source boundary condition at top surface
   - Consider energy-dependent interface transmission
   - Verify spin accumulation boundary conditions allow minority spin buildup
4. Test and validate against experimental data
5. Consider Option 2 if more accuracy needed

**Estimated Effort:**
- Option 1: ~50 lines, 1-2 days
- Boundary condition enhancements: ~30-50 lines, 1 day
- Option 2: ~200 lines, 4-5 days  
- Option 3: ~300 lines, 1-2 weeks
