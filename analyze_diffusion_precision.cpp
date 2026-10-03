// Analysis: Precision requirements for charge and spin current diffusion
// Including Thomas tridiagonal solver and gradient calculations

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <vector>

// Typical values for 0.1 fs timestep
const double dt = 0.1e-15;          // 1e-16 s (0.1 fs)
const double dz_coarse = 10e-10;    // 1e-9 m (coarse microcell ~10 Angstrom)
const double dz_fine = 1e-12;       // 1e-12 m (fine grid ~0.01 Angstrom)
const double D = 0.001;             // 1e-3 m^2/s (typical diffusion coefficient)
const double ns_typical = 1e10;     // C/m^3 (typical charge density)
const double S_typical = 1.48e7;    // C/m^3 (typical spin accumulation)

// Thomas solver test
void test_thomas_precision() {
    std::cout << "\n=== THOMAS TRIDIAGONAL SOLVER PRECISION ===\n\n";
    
    const int n = 100;  // Typical fine grid size
    std::vector<double> a_d(n, 0.0), b_d(n, 1.0), c_d(n, 0.0), d_d(n, 1.0);
    std::vector<float> a_f(n, 0.0f), b_f(n, 1.0f), c_f(n, 0.0f), d_f(n, 1.0f);
    
    // Build typical diffusion matrix
    // b[k] = 1.0 + dt*alpha, where alpha = D/(dz*dz) or similar
    double alpha = D / (dz_fine * dz_fine);  // 1e-3 / 1e-24 = 1e21 (!)
    double dt_alpha = dt * alpha;            // 1e-16 * 1e21 = 1e5
    
    std::cout << "Typical diffusion matrix coefficients:\n";
    std::cout << "  alpha = D/(dz*dz) = " << D << " / (" << dz_fine << ")^2 = " << alpha << "\n";
    std::cout << "  dt*alpha = " << dt << " * " << alpha << " = " << dt_alpha << "\n";
    
    // Diagonal dominance: b[k] = 1.0 + dt*alpha
    // For stability: dt*alpha < 1 (implicit method is unconditionally stable, but...)
    double b_diag = 1.0 + dt_alpha;
    std::cout << "  b[k] (diagonal) = 1.0 + dt*alpha = " << b_diag << "\n";
    std::cout << "  Off-diagonal: a[k] = c[k] = -dt*alpha = " << -dt_alpha << "\n\n";
    
    // Test division precision in Thomas algorithm
    // Step 1: cp[0] = c[0] / b[0]
    double cp0_d = -dt_alpha / b_diag;
    float cp0_f = static_cast<float>(-dt_alpha) / static_cast<float>(b_diag);
    
    std::cout << "Thomas algorithm division precision:\n";
    std::cout << "  cp[0] = c[0]/b[0] (double): " << cp0_d << "\n";
    std::cout << "  cp[0] = c[0]/b[0] (float):  " << cp0_f << "\n";
    std::cout << "  Relative error: " << std::abs((cp0_d - cp0_f) / cp0_d) << "\n\n";
    
    // For fine grid with dt=0.1fs, dt*alpha can be very large!
    // Actually, wait - let me recalculate
    // alpha_edge[i] = C/dz where C = 1/Rseries, and Rseries ~ dz/(2*D)
    // So alpha ~ 2*D/dz^2? Let me use a more realistic value
    double alpha_realistic = 2.0 * D / (dz_fine * dz_fine);  // Still large but more realistic
    double dt_alpha_real = dt * alpha_realistic;
    std::cout << "More realistic alpha (considering interface resistance):\n";
    std::cout << "  alpha ~ 2*D/dz^2 = " << alpha_realistic << "\n";
    std::cout << "  dt*alpha = " << dt_alpha_real << "\n";
    
    if (dt_alpha_real > 1e6) {
        std::cout << "  WARNING: dt*alpha is very large - may cause precision issues!\n";
    }
}

void test_gradient_precision() {
    std::cout << "\n=== GRADIENT CALCULATION PRECISION ===\n\n";
    
    // Spin accumulation gradient: dS/dz
    // S ~ 1.48e7 C/m^3, dz ~ 1e-12 m
    // dS/dz = (S[i+1] - S[i]) / dz
    
    // Small variation case (typical)
    double S1 = S_typical;
    double S2 = S_typical * 1.001;  // 0.1% difference
    double dS_dz_double = (S2 - S1) / dz_fine;
    float dS_dz_float = (static_cast<float>(S2) - static_cast<float>(S1)) / static_cast<float>(dz_fine);
    
    std::cout << "Spin accumulation gradient (small variation):\n";
    std::cout << "  S[i] = " << S1 << " C/m^3\n";
    std::cout << "  S[i+1] = " << S2 << " C/m^3 (0.1% larger)\n";
    std::cout << "  dS/dz (double): " << dS_dz_double << " C/m^4\n";
    std::cout << "  dS/dz (float):  " << dS_dz_float << " C/m^4\n";
    std::cout << "  Relative error: " << std::abs((dS_dz_double - dS_dz_float) / dS_dz_double) << "\n\n";
    
    // Large gradient case (interface, shock)
    double S1_large = S_typical;
    double S2_large = S_typical * 2.0;  // 2x difference (steep interface)
    double dS_dz_large_d = (S2_large - S1_large) / dz_fine;
    float dS_dz_large_f = (static_cast<float>(S2_large) - static_cast<float>(S1_large)) / static_cast<float>(dz_fine);
    
    std::cout << "Spin accumulation gradient (large variation - interface):\n";
    std::cout << "  S[i] = " << S1_large << " C/m^3\n";
    std::cout << "  S[i+1] = " << S2_large << " C/m^3 (2x larger)\n";
    std::cout << "  dS/dz (double): " << dS_dz_large_d << " C/m^4\n";
    std::cout << "  dS/dz (float):  " << dS_dz_large_f << " C/m^4\n";
    std::cout << "  Relative error: " << std::abs((dS_dz_large_d - dS_dz_large_f) / dS_dz_large_d) << "\n\n";
    
    // Diffusion current: J_diff = -D * dS/dz
    double J_diff_d = -D * dS_dz_double;
    float J_diff_f = static_cast<float>(-D) * dS_dz_float;
    
    std::cout << "Diffusion current calculation:\n";
    std::cout << "  J_diff = -D * dS/dz\n";
    std::cout << "  J_diff (double): " << J_diff_d << " A/m^2\n";
    std::cout << "  J_diff (float):  " << J_diff_f << " A/m^2\n";
    std::cout << "  Relative error: " << std::abs((J_diff_d - J_diff_f) / J_diff_d) << "\n\n";
}

void test_charge_current_diffusion() {
    std::cout << "\n=== CHARGE CURRENT DIFFUSION PRECISION ===\n\n";
    
    // Charge density ns: C/m^3
    // Diffusion equation: d(ns)/dt = D * d^2(ns)/dz^2 - ns/tau_s + source
    // Implicit: (ns_new - ns_old)/dt = D * d^2(ns_new)/dz^2 - ...
    
    // Matrix coefficient: dt*D/(dz*dz)
    double dt_D_dz2_coarse = dt * D / (dz_coarse * dz_coarse);
    double dt_D_dz2_fine = dt * D / (dz_fine * dz_fine);
    
    std::cout << "Charge diffusion matrix coefficients:\n";
    std::cout << "  Coarse grid (dz = " << dz_coarse << " m):\n";
    std::cout << "    dt*D/(dz*dz) = " << dt_D_dz2_coarse << "\n";
    std::cout << "  Fine grid (dz = " << dz_fine << " m):\n";
    std::cout << "    dt*D/(dz*dz) = " << dt_D_dz2_fine << "\n\n";
    
    // Diagonal: b[k] = 1.0 + dt/tau_s + dt*D/(dz*dz) + ...
    double tau_s = 1e-12;  // 1 ps (typical)
    double dt_tau = dt / tau_s;  // 1e-16 / 1e-12 = 1e-4
    double b_diag_coarse = 1.0 + dt_tau + 2.0 * dt_D_dz2_coarse;
    
    std::cout << "Diagonal coefficient b[k] (coarse):\n";
    std::cout << "  b[k] = 1.0 + dt/tau_s + 2*dt*D/(dz*dz)\n";
    std::cout << "       = 1.0 + " << dt_tau << " + " << (2.0*dt_D_dz2_coarse) << "\n";
    std::cout << "       = " << b_diag_coarse << "\n";
    
    // Check if small contributions are lost
    double small_contrib = dt_tau;
    float small_contrib_f = static_cast<float>(dt_tau);
    
    std::cout << "\nSmall contribution precision:\n";
    std::cout << "  dt/tau_s (double): " << small_contrib << "\n";
    std::cout << "  dt/tau_s (float):  " << small_contrib_f << "\n";
    std::cout << "  When added to 1.0: 1.0 + dt/tau_s (double) = " << (1.0 + small_contrib) << "\n";
    std::cout << "  When added to 1.0: 1.0 + dt/tau_s (float)  = " << (1.0f + small_contrib_f) << "\n";
    if (small_contrib < 1e-7) {
        std::cout << "  WARNING: dt/tau_s may be lost in float precision when added to 1.0!\n";
    }
}

void test_spin_diffusion_matrix() {
    std::cout << "\n=== SPIN DIFFUSION MATRIX PRECISION ===\n\n";
    
    // Spin diffusion: alpha_edge[i] = C/dz where C = 1/Rseries
    // Rseries = dz/(2*Di) + Rint + dz/(2*Dj)
    // For typical: Di = Dj = D, Rint = 0
    // C = 1 / (dz/D) = D/dz
    // alpha = C/dz = D/(dz*dz)
    
    double alpha = D / (dz_fine * dz_fine);  // Same as charge diffusion
    double dt_alpha_half = 0.5 * dt * alpha;  // dt_half * alpha
    
    std::cout << "Spin diffusion (fine grid):\n";
    std::cout << "  alpha = D/(dz*dz) = " << alpha << " 1/s\n";
    std::cout << "  dt_half * alpha = " << dt_alpha_half << "\n";
    
    // Matrix: b[i] = 1.0 + dt_half*(alpha[i-1] + alpha[i])
    double b_diag = 1.0 + dt_alpha_half * 2.0;
    
    std::cout << "  b[i] = 1.0 + dt_half*(alpha[i-1] + alpha[i]) = " << b_diag << "\n";
    
    // Off-diagonal: a[i] = -dt_half*alpha[i-1]
    double a_off = -dt_alpha_half;
    std::cout << "  a[i] = c[i] = -dt_half*alpha = " << a_off << "\n\n";
    
    // Test precision
    float b_diag_f = 1.0f + static_cast<float>(dt_alpha_half * 2.0);
    float a_off_f = static_cast<float>(-dt_alpha_half);
    
    std::cout << "Precision test:\n";
    std::cout << "  b[i] (double): " << b_diag << "\n";
    std::cout << "  b[i] (float):  " << b_diag_f << "\n";
    std::cout << "  Relative error: " << std::abs((b_diag - b_diag_f) / b_diag) << "\n";
    
    if (b_diag > 1e6) {
        std::cout << "  WARNING: Matrix diagonal is very large - precision may be lost in float!\n";
    }
}

int main() {
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "=== DIFFUSION PRECISION ANALYSIS FOR 0.1 fs TIMESTEP ===\n";
    
    test_charge_current_diffusion();
    test_spin_diffusion_matrix();
    test_thomas_precision();
    test_gradient_precision();
    
    std::cout << "\n=== SUMMARY & RECOMMENDATIONS ===\n\n";
    std::cout << "CRITICAL PRECISION REQUIREMENTS:\n\n";
    
    std::cout << "1. THOMAS TRIDIAGONAL SOLVER:\n";
    std::cout << "   - Requires high precision for division operations\n";
    std::cout << "   - Diagonal dominance important: b[i] >> |a[i]|, |c[i]|\n";
    std::cout << "   - For dt*alpha >> 1, b[i] becomes very large\n";
    std::cout << "   - RECOMMENDATION: Keep double precision for all solver arrays\n\n";
    
    std::cout << "2. GRADIENT CALCULATIONS (dS/dz, d(ns)/dz):\n";
    std::cout << "   - Involves subtraction of similar large numbers\n";
    std::cout << "   - Small relative differences amplified by 1/dz (very small)\n";
    std::cout << "   - Example: (S[i+1] - S[i]) / 1e-12 can lose precision\n";
    std::cout << "   - RECOMMENDATION: Keep double for S, ns, ne arrays\n\n";
    
    std::cout << "3. MATRIX COEFFICIENTS:\n";
    std::cout << "   - dt*D/(dz*dz) can be very large for fine grids\n";
    std::cout << "   - Small contributions (dt/tau) may be lost when added to 1.0\n";
    std::cout << "   - RECOMMENDATION: Keep double for all matrix coefficients\n\n";
    
    std::cout << "4. SAFE FOR FLOAT:\n";
    std::cout << "   - Diffusion coefficient D (if constant per cell)\n";
    std::cout << "   - Geometric factors (dz, 1/dz) if pre-computed\n";
    std::cout << "   - Material constants (beta_cond, beta_diff) in fine-grid arrays\n";
    std::cout << "   - BUT: Memory savings minimal, risk significant\n\n";
    
    std::cout << "OVERALL RECOMMENDATION:\n";
    std::cout << "  Keep ALL diffusion-related computations in DOUBLE precision.\n";
    std::cout << "  The precision requirements for 0.1 fs timestep are too stringent\n";
    std::cout << "  for float, especially in Thomas solver and gradient calculations.\n";
    
    return 0;
}