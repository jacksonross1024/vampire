# Physics Equations: Spin Currents and Thermal Solvers

This document collects all relevant physics equations implemented in `spincurrents.cpp` and the thermal solvers (`thermal_solver.cpp` and `local_temperature_pulse.cpp`).

---

## Table of Contents

1. [Physical Constants](#physical-constants)
2. [Spin Accumulation Equations](#spin-accumulation-equations)
3. [Charge Current Equations](#charge-current-equations)
4. [Thermal Gradient Equations](#thermal-gradient-equations)
5. [Effective Field and Torque](#effective-field-and-torque)
6. [Numerical Methods](#numerical-methods)

---

## Physical Constants

```cpp
ℏ = 1.05457162e-34 J·s        // Reduced Planck constant
μ_B = 9.27400968e-24 J/T       // Bohr magneton
e = 1.60217662e-19 C           // Elementary charge
h = 6.62607015e-34 J·s         // Planck constant
c = 2.99792458e8 m/s           // Speed of light
V_atom = 2.89e-30 m³           // Atomic cell volume
```

---

## Spin Accumulation Equations

### Main Spin Accumulation Evolution (Serban Equation 5)

The spin accumulation **S** (charge density, C/m³) evolves according to:

```
∂S/∂t = -∇·J_drift - ∇·J_diff + ω × S + De·[(S - sa_inf·m̂)/λ_sf²] + De·[m̂ × (S × m̂)/λ_φ²] + S_demag + S_pump
```

Where:
- **S** = (S_x, S_y, S_z) is the spin accumulation vector (C/m³)
- **J_drift** = drift spin current (A/m²)
- **J_diff** = diffusive spin current (A/m²)
- **ω** = exchange precession frequency (rad/s)
- **De** = electron diffusion constant (m²/s)
- **sa_inf** = equilibrium spin accumulation (C/m³)
- **m̂** = normalized magnetization direction
- **λ_sf** = spin-flip length (m)
- **λ_φ** = dephasing length (m)
- **S_demag** = demagnetization-driven source term
- **S_pump** = spin pumping source term

### Component Breakdown

#### 1. Drift Current Divergence
```
∇·J_drift = -∇·(β_c · J_c · m̂)
```
Where:
- **β_c** = spin Hall angle (dimensionless)
- **J_c** = charge current density (A/m²)
- **m̂** = normalized magnetization direction

Finite difference form (1D, z-direction):
```
(∇·J_drift)_i = -(J_drift[i+1] - J_drift[i]) / dz
J_drift[i] = -β_c[i] · J_c[i] · m̂[i]
```

#### 2. Diffusive Current Divergence
```
∇·J_diff = -∇·(De · ∇S)
```
Finite difference form (1D):
```
(∇·J_diff)_i = -De · (S[i+1] - 2·S[i] + S[i-1]) / dz²
```

#### 3. Exchange Precession
```
ω × S = -ω · (S × m̂)
```
Where:
```
ω = J_sd / (2·ℏ)
```
- **J_sd** = sd-exchange constant (J)

Component form:
```
(ω × S)_x = -ω · (S_y·m̂_z - S_z·m̂_y)
(ω × S)_y = -ω · (S_z·m̂_x - S_x·m̂_z)
(ω × S)_z = -ω · (S_x·m̂_y - S_y·m̂_x)
```

#### 4. Spin-Flip Relaxation
```
De · [(S - sa_inf·m̂) / λ_sf²]
```

Where:
```
λ_sf = λ_sf0 · √(1 - β_c·β_d)
```
- **λ_sf0** = base spin-flip length (m)
- **β_d** = diffusion spin Hall angle (dimensionless)

#### 5. Dephasing
```
De · [m̂ × (S × m̂) / λ_φ²]
```

Where **S_perp** = S - (S·m̂)m̂ is the perpendicular component:
```
S_perp = S - (S·m̂)·m̂
```

#### 6. Demagnetization-Driven Source
```
S_demag = -χ · (d|m|/dt) · m̂_ref
```
Where:
- **χ** = demagnetization susceptibility (dimensionless)
- **d|m|/dt** = rate of change of magnetization magnitude
- **m̂_ref** = reference magnetization direction (frozen when |m| is small)

#### 7. Spin Pumping Source
```
S_pump = -χ_sp · (dm/dt)_perp
```
Where:
- **χ_sp** = spin pumping coefficient (C·s/m³)
- **(dm/dt)_perp** = perpendicular component of magnetization rate of change

Spin pumping coefficient:
```
χ_sp = (μ_B · σ · ℏ) / (e² · λ_J²)
```
Where:
- **σ** ≈ 1e6 S/m (typical metal conductivity)
- **λ_J** = exchange length = √(De · ℏ / J_sd)

---

## Charge Current Equations

### Non-Equilibrium Charge Density (ns) - Equation 2

```
∂ns/∂t = -∇·J_ns - ns/τ_s + S_laser
```

Where:
- **ns** = non-equilibrium hot electron density (C/m³)
- **J_ns** = superdiffusive current (A/m²)
- **τ_s** = population relaxation time (s)
- **S_laser** = laser source term (C/(m³·s))

#### Superdiffusive Current
```
J_ns = v_e · ns - De · ∇ns
```

Where:
- **v_e** = hot electron velocity (m/s)
- **De** = electron diffusion constant (m²/s)

#### Hot Electron Velocity (Superdiffusive Model)

Base velocity:
```
v_e0 = √(De / d_opt)
```
or user-specified:
```
v_e0 = v_e0_user
```

Temperature-dependent scaling:
```
v_e(T) = v_e0 · √(Te / T_ref)
```

Time-dependent decay (energy relaxation):
```
v_e(t) = v_e(T) · exp(-(t - t_peak) / τ_e)  for t > t_peak
```
Where:
- **τ_e** = energy relaxation time (~50-200 fs)
- **t_peak** = laser pulse peak time

#### Laser Source Term
```
S_laser = (e · η · Q · λ) / (h · c)
```
Where:
- **η** = quantum efficiency (dimensionless)
- **Q** = laser power density (W/m³)
- **λ** = laser wavelength (m)
- **h** = Planck constant
- **c** = speed of light

Laser power density (spatial attenuation):
```
Q(z,t) = Q_0 · exp(-(z_max - z) / d) · gaussian(t - t_0 - t_eq)
```
Where:
- **Q_0** = peak power density (W/m³)
- **d** = optical absorption length (m)
- **z_max** = top surface position (m)
- **gaussian(t)** = temporal envelope

Gaussian temporal envelope:
```
gaussian(t) = exp(-0.5 · ((t - t_0) / σ)²)
σ = FWHM / (2·√(2·ln(2)))
```

### Excess Charge Density (ne) - Equation 3

```
∂ne/∂t = -∇·J_c
```

Where:
- **ne** = excess charge density (C/m³)
- **J_c** = total charge current (A/m²)

### Total Charge Current (J_c) - Equation 1

```
J_c = -De · ∇ne + J_ns + J_seebeck
```

#### Seebeck Current
```
J_seebeck = -S · ∇Te
```
Where:
- **S** = Seebeck coefficient (V/K)
- **Te** = electron temperature (K)

Temperature-dependent Seebeck coefficient:
```
S(Te) = S_0 · (1 + α_S · (Te - T_ref) / T_ref)
```
Where:
- **α_S** ≈ 0.5 (temperature coefficient)
- **T_ref** = reference temperature (K)

---

## Thermal Gradient Equations

### Two-Temperature Model (TTM)

The thermal solver implements the standard two-temperature model:

```
Ce · dTe/dt = -G·(Te - Tp) + ∇·(κe·∇Te) + S
Cp(Tp) · dTp/dt = G·(Te - Tp) + ∇·(κp·∇Tp) - Q_sink
```

Where:
- **Te** = electron temperature (K)
- **Tp** = phonon temperature (K)
- **Ce** = electron heat capacity (J/(m³·K))
- **Cp(Tp)** = phonon heat capacity (temperature-dependent, J/(m³·K))
- **G** = electron-phonon coupling constant (J/(s·m³·K))
- **κe** = electron thermal conductivity (J/(s·m·K))
- **κp** = phonon thermal conductivity (J/(s·m·K))
- **S** = laser source term (W/m³)
- **Q_sink** = substrate cooling term (W/m³)

### Electron Heat Capacity

For free electron model:
```
Ce = Ce_0 · Te
```

However, the code uses constant **Ce** (per cell), so:
```
dTe/dt = [-G·(Te - Tp) + ∇·(κe·∇Te) + S] / (Ce · Te)
```

### Phonon Heat Capacity

#### Debye Model (Low Temperature, T_D/T ≥ 12)
```
Cp(T) = Cp_0 · (4π⁴ / 5) · (T / T_D)³
```

Where:
- **T_D** = Debye temperature (K)
- **Cp_0** = high-temperature phonon heat capacity (J/(m³·K))

#### Einstein Model (High Temperature, T_D/T < 12)
```
Cp(T) = Cp_0 · f_Debye(T_D/T)
```

Where **f_Debye** is computed from Debye integral:
```
f_Debye(x) = 3 · ∫[0 to x] (u⁴·e^u / (e^u - 1)²) du / x³
```

### Temperature-Dependent Electron Thermal Conductivity

Chen–Beraun (e–ph Drude limit):
```
κe(Te, Tp) = κe_0 · (Te / Tp)
```

Where:
- **κe_0** = equilibrium electron thermal conductivity (J/(s·m·K)), the material input
- **Tp** = local phonon temperature (K)

Continuum product rule (identity only):
```
∇·(κe ∇Te) = κe ∇²Te + ∇κe · ∇Te
```

with κe = κe_0 Te/Tp. This is **not** added as a second finite-difference term. The conservative stencil already is the discrete divergence of the face flux q = −κ ∇Te, so ∇κ·∇Te is included when left and right faces have different κ.

### Thermal Diffusion (Finite Difference)

```
∇·(κ∇T) ≈ Σ_neighbors [κ_harmonic · (T_neighbor − T_self) / Δz²]
```

Electron interface conductivity:
```
κe_cell = κe_0 · (Te / Tp)
κ_harmonic = 2·κ_self·κ_neighbor / (κ_self + κ_neighbor)
```

Phonons use κp_0 with no Te/Tp factor.

### Electron-Phonon Coupling

```
dTe/dt (coupling) = -G·(Te - Tp) / (Ce · Te)
dTp/dt (coupling) = G·(Te - Tp) / Cp(Tp)
```

### Laser Source Term

```
S(z,t) = Q_0 · exp(-(z_max - z) / d) · gaussian(t - t_0 - t_eq)
```

Where:
- **Q_0** = peak power density (W/m³)
- **d** = optical absorption length (m)
- **z_max** = top surface position (m)
- **gaussian(t)** = temporal envelope

### Substrate Cooling (Bottom Boundary)

```
Q_sink = T_cool · (Tp - T_eq)
```

Where:
- **T_cool** = heat sink coupling constant (1/s)
- **T_eq** = equilibration temperature (K)

---

## Effective Field and Torque

### Effective Field from Spin Accumulation

```
H_eff = (V_atom · J_sd · Sa) / (e · μ_B)
```

Where:
- **V_atom** = atomic cell volume (m³)
- **J_sd** = sd-exchange constant (J)
- **Sa** = spin accumulation (C/m³)
- **e** = elementary charge (C)
- **μ_B** = Bohr magneton (J/T)

Units: (m³ · J · C/m³) / (C · J/T) = T (Tesla)

### Spin Torque

```
τ = m × H_eff
```

Component form:
```
τ_x = m_y · H_eff_z - m_z · H_eff_y
τ_y = m_z · H_eff_x - m_x · H_eff_z
τ_z = m_x · H_eff_y - m_y · H_eff_x
```

---

## Numerical Methods

### Time Integration

#### Spin Accumulation: Strang Splitting

1. **First diffusion half-step** (implicit, backward Euler):
   ```
   S* = S_n + (dt/2) · ∇·(De·∇S*)
   ```

2. **Explicit step** (Heun for first step, AB2 for subsequent):
   - **Heun (RK2)** for first step:
     ```
     k1 = f(S*)
     S_pred = S* + dt · k1
     k2 = f(S_pred)
     S** = S* + 0.5·dt·(k1 + k2)
     ```
   - **Adams-Bashforth 2 (AB2)** for subsequent steps:
     ```
     S** = S* + dt · (1.5·f_n - 0.5·f_{n-1})
     ```

3. **Second diffusion half-step** (implicit, backward Euler):
   ```
   S_{n+1} = S** + (dt/2) · ∇·(De·∇S_{n+1})
   ```

#### Charge Density: Implicit Backward Euler

For **ns** (non-equilibrium density):
```
(I + dt·A) · ns_{n+1} = ns_n + dt·S_laser
```

Where **A** is the advection-diffusion-reaction operator:
```
A = (v_e/dz) · (upwind) + (De/dz²) · (Laplacian) + (1/τ_s) · I
```

For **ne** (excess density):
```
(I + dt·De·∇²) · ne_{n+1} = ne_n - dt·∇·J_ns
```

#### Thermal Solver: Split-Step

1. **Calculate diffusion terms** (using old temperatures):
   ```
   ΔT_diff = dt · ∇·(κ·∇T) / (C·T)
   ```

2. **Update with coupling + diffusion + source**:
   ```
   T_new = T_old + ΔT_diff + dt·(coupling + source) / (C·T)
   ```

### Spatial Discretization

#### Fine Grid (Spin Accumulation)
- Fine grid spacing: **dz_fine = dz_coarse / nsub**
- **nsub** = subdivisions per coarse microcell
- Interpolation between coarse cells (optional, disabled by default)

#### Coarse Grid (Charge Current)
- Coarse grid spacing: **dz_coarse = micro_cell_thickness** (Å)
- Edge-averaged properties for interfaces

### Boundary Conditions

#### Spin Accumulation
- **Neumann boundaries**: Total spin current **J_total = 0** at surfaces
- Drift and diffusive currents cancel at boundaries

#### Charge Current
- **Zero flux**: **J_c = 0** at top and bottom surfaces
- Open circuit condition

#### Thermal
- **Bottom**: Substrate cooling (phonon temperature)
- **Top**: Zero flux (insulated)

### Interface Treatment

#### Interface Resistance (Spin)
```
R_int = r_int_edge (Ω·m²)
```

Interface coupling:
```
α_edge = 1 / (dz/(2·D_L) + R_int + dz/(2·D_R))
```

#### Harmonic Mean (Thermal Conductivity)
```
κ_harmonic = 2·κ_L·κ_R / (κ_L + κ_R)
```

Ensures flux continuity at material interfaces.

---

## Key Physical Relationships

### Spin Diffusion Length
```
λ_sf = λ_sf0 · √(1 - β_c·β_d)
```

### Exchange Length
```
λ_J = √(De · ℏ / J_sd)
```

### Spin-Flip Time
```
τ_sf = λ_sf² / De
```

### Dephasing Time
```
τ_φ = λ_φ² / De
```

### Exchange Precession Frequency
```
ω = J_sd / (2·ℏ)
```

### Spin Pumping Coefficient
```
χ_sp = (μ_B · σ · ℏ) / (e² · λ_J²)
```

---

## Temperature Scaling

### Diffusion Constant
```
De(Te) = De_0 · f(Te/T_ref)
```
Currently: **f(Te/T_ref) = 1.0** (no scaling, but mechanism in place)

### Spin Diffusion Length
```
λ_sdl(Te) = λ_sdl_0 · f(Te/T_ref)
```

### Population Relaxation Time
```
τ_s(Te) = τ_s_0 · f(Te/T_ref)
```

### Seebeck Coefficient
```
S(Te) = S_0 · (1 + α_S · (Te - T_ref) / T_ref)
```
Where **α_S ≈ 0.5** for typical metals.

---

## Numerical Stability

### Root Temperature Storage
Temperatures stored as **√T** to improve numerical stability at low temperatures:
```
T_stored = √T
T_actual = T_stored²
```

### RHS Clamping
Right-hand side terms clamped to prevent overflow:
```
RHS_clamped = max(-1e25, min(1e25, RHS))
```

### Finite Checks
All computed values checked for **NaN** and **Inf**:
```
if(!isfinite(value)) value = 0.0;
```

### Gradient Limiting
Spin accumulation gradients limited to prevent numerical issues:
```
|∇S|_max = 1e20 C/m⁴
```

---

## References

1. **Serban et al., AOS Spin Valves** - Main reference for spin accumulation equations
2. **Lepadatu, Unified treatment of spin torques** - Spin pumping and effective field
3. **Two-Temperature Model (TTM)** - Standard electron-phonon thermal model
4. **Debye Model** - Phonon heat capacity at low temperatures
5. **Einstein Model** - Phonon heat capacity at high temperatures

---

*Document generated from code analysis of `spincurrents.cpp`, `thermal_solver.cpp`, and `local_temperature_pulse.cpp`*
