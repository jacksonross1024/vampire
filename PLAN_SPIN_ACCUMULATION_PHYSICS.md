# Plan: Physics Upgrades for 1D Spin Accumulation Solver

This document outlines the plan to introduce spin pumping terms into the 1D transient spin-current solver (`spincurrents.cpp`). These changes align the solver with the drift-diffusion framework of Lepadatu (2017, 2021) and the Spin Valves paper (Serban et al.).

## Literature Audit Summary

### Lepadatu 2017 ("Unified treatment of spin torques...")
Key equations for spin accumulation dynamics (Eq. 3, 8):
$$
\frac{\partial \mathbf{S}}{\partial t} = -\nabla \cdot \mathbf{J}_S - D_e \left[ \frac{\mathbf{S}}{\lambda_{sf}^2} + \frac{\mathbf{S} \times \mathbf{m}}{\lambda_J^2} + \frac{\mathbf{m} \times (\mathbf{S} \times \mathbf{m})}{\lambda_\phi^2} \right]
$$

Spin pumping (Eq. 14) as an **interfacial flux**:
$$
\mathbf{J}_S^{pump} = \frac{\mu_B}{2\pi} \left[ \text{Re}\{g^{\uparrow\downarrow}\} \mathbf{m} \times \frac{\partial \mathbf{m}}{\partial t} + \text{Im}\{g^{\uparrow\downarrow}\} \frac{\partial \mathbf{m}}{\partial t} \right]
$$

### Spin Valves Paper (Serban et al.)
Modified spin accumulation equation (Eq. 5) with **bulk spin pumping source**:
$$
\frac{\partial \mathbf{S}}{\partial t} = -\nabla \cdot \mathbf{J}_S - \chi_{sp} \frac{\partial \mathbf{m}}{\partial t} - D_e \left[ \frac{\mathbf{S}}{\lambda_{sf}^2} + \frac{\mathbf{S} \times \mathbf{m}}{\lambda_J^2} + \frac{\mathbf{m} \times (\mathbf{S} \times \mathbf{m})}{\lambda_\phi^2} \right]
$$

where the spin pumping coefficient is:
$$
\chi_{sp} = \frac{\mu_B \sigma}{e^2 \lambda_J^2}
$$

**Note**: The pumping term is $-\chi_{sp} \frac{\partial \mathbf{m}}{\partial t}$ (proportional to the vector time derivative), NOT $\mathbf{m} \times \frac{\partial \mathbf{m}}{\partial t}$.

---

## Audit Conclusions

### 1. Equilibrium Source Term (`sa_infinity`)

**Original Proposal (INCORRECT)**: Add source term `+ S_inf * m / tau_sf` to relax spin accumulation to an equilibrium value.

**Literature Finding**: In Lepadatu's framework, `S` represents *non-equilibrium* spin accumulation. The relaxation terms drive `S → 0`, not to some finite `S_infinity`.

**Decision**: Do NOT implement equilibrium source term. The existing relaxation (`-S/tau_sf`) is correct.

**Note on `sa_infinity`**: This parameter may be relevant for steady-state DC spin transport (Valet-Fert model), but is NOT used in the transient ultrafast regime modeled by Lepadatu.

### 2. Spin Pumping Source Term

**Original Proposal (PARTIALLY INCORRECT)**: Add `-div(J_pump)` where `J_pump ~ m x dm/dt`.

**Literature Finding**: The Spin Valves paper uses a simpler **bulk source term**:
$$
\frac{\partial \mathbf{S}}{\partial t}\bigg|_{pump} = -\chi_{sp} \frac{\partial \mathbf{m}}{\partial t}
$$

This is NOT `m x dm/dt`. It is simply `-dm/dt` scaled by `chi_sp`.

**Decision**: Implement spin pumping as a bulk source term proportional to $-\frac{\partial \mathbf{m}}{\partial t}$.

---

## Implementation Plan

### Step 1: Track Previous Magnetization Vector
We need the vector time derivative $\frac{d\mathbf{m}}{dt}$. Currently, only the magnitude rate is tracked.

*   **Files**: `src/spintorque/data.cpp`, `internal.hpp`
*   **Action**: Add `std::vector<double> sc1d_m_prev_vec` (size $3 \times N_{microcells}$).
*   **File**: `src/spintorque/spincurrents.cpp` (`calculate_spin_accumulation_1d`)
*   **Action**: Initialize and update `sc1d_m_prev_vec` at each coarse timestep.

### Step 2: Compute Vector Time Derivative
Inside the coarse loop in `calculate_spin_accumulation_1d`:

```cpp
// Vector dm/dt for spin pumping
double dmdt_vec_x = (mx - sc1d_m_prev_vec[3*cell+0]) / dt;
double dmdt_vec_y = (my - sc1d_m_prev_vec[3*cell+1]) / dt;
double dmdt_vec_z = (mz - sc1d_m_prev_vec[3*cell+2]) / dt;

// Update previous magnetization
sc1d_m_prev_vec[3*cell+0] = mx;
sc1d_m_prev_vec[3*cell+1] = my;
sc1d_m_prev_vec[3*cell+2] = mz;
```

### Step 3: Calculate Spin Pumping Coefficient
The coefficient from the Spin Valves paper:
$$
\chi_{sp} = \frac{\mu_B \sigma}{e^2 \lambda_J^2}
$$

In the code:
```cpp
const double chi_sp = (kMuB * sigma) / (kE * kE * lambda_J * lambda_J);
```

where:
*   `sigma` = electrical conductivity (material parameter, needs to be added or derived)
*   `lambda_J` = exchange rotation length (related to `sd_exchange` via $\lambda_J = \sqrt{D_e / J_{sd}}$)

### Step 4: Add Spin Pumping Source to Heun Solver
In the `compute_rhs` lambda function:

```cpp
// Spin pumping source: -chi_sp * dm/dt
// chi_sp interpolated to fine grid like other parameters
const double src_pump_x = -chi_sp[i] * dmdt_vec_x_interp[i];
const double src_pump_y = -chi_sp[i] * dmdt_vec_y_interp[i];
const double src_pump_z = -chi_sp[i] * dmdt_vec_z_interp[i];

kout[3*i+0] += src_pump_x;
kout[3*i+1] += src_pump_y;
kout[3*i+2] += src_pump_z;
```

### Step 5: Add Enable Flag and Material Parameter
*   **Flag**: `sc1d_spin_pumping_enable` (default: false)
*   **Material Parameter**: May need to add `sigma` (conductivity) to material file, or compute from existing parameters.

---

## Verification Tests

1.  **Static Test**: With $\frac{d\mathbf{m}}{dt} = 0$ and no current, $\mathbf{S}$ should relax to zero (confirming no spurious equilibrium source).

2.  **Demagnetization Test**: During laser-induced demagnetization ($\frac{d|\mathbf{m}|}{dt} < 0$), spin accumulation should be generated with direction anti-parallel to $\frac{d\mathbf{m}}{dt}$.

3.  **Precession Test**: During precessional dynamics, the spin pumping source should inject transverse spin accumulation that diffuses away.

---

## Notes

*   **Independence from Demag Source**: The spin pumping term is independent of the existing demagnetization source (`-chi_demag * d|m|/dt * m`). The demag source depends on magnitude change; spin pumping depends on vector change.

*   **Typical Values** (from Spin Valves paper):
    *   Low pumping regime: $\chi_{sp} = 2.4$ kA/m ($\sigma = 10^6$ S/m)
    *   High pumping regime: $\chi_{sp} = 24$ kA/m ($\sigma = 10^7$ S/m)
    *   $\lambda_J = 4$ nm, $D_e = 10^{-4}$ m²/s

*   **Alternative Approach**: Lepadatu 2017 treats spin pumping as an interfacial flux boundary condition. This could be implemented as an alternative mode if needed for comparison.

