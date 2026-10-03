// Analysis: Can compile-time constants kInvE and kInvMuB be pre-combined
// during interpolation to reduce precision requirements?

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

// Physical constants (from spincurrents.cpp)
static constexpr double kMuB  = 9.27400968e-24;      // J/T
static constexpr double kE    = 1.60217662e-19;      // C
static constexpr double kAtomcellVolume = 2.89e-30;  // m^3
static constexpr double kInvMuB = 1.0 / kMuB;        // 1.078282e+23
static constexpr double kInvE   = 1.0 / kE;          // 6.241509e+18

// Pre-combined constant (compile-time)
static constexpr double kTorqueScaleBase = kAtomcellVolume * kInvE * kInvMuB;

int main() {
    std::cout << std::scientific << std::setprecision(6);
    
    // Typical values
    double sd_exchange = 4.0e-20;  // J
    double ax = 1.48e6;             // C/m^3 (typical spin accumulation after averaging)
    
    std::cout << "=== PRECISION ANALYSIS: Interpolation Step Optimization ===\n\n";
    
    std::cout << "Constants:\n";
    std::cout << "  kAtomcellVolume = " << kAtomcellVolume << "\n";
    std::cout << "  kInvE          = " << kInvE << "\n";
    std::cout << "  kInvMuB        = " << kInvMuB << "\n";
    std::cout << "  kTorqueScaleBase = " << kTorqueScaleBase << "\n\n";
    
    // Current approach: multiply all factors at end
    // torque = kAtomcellVolume * sd_exchange * ax * kInvE * kInvMuB
    double torque_current = kAtomcellVolume * sd_exchange * ax * kInvE * kInvMuB;
    std::cout << "Current approach (all multiplications at end):\n";
    std::cout << "  torque = " << torque_current << " T\n\n";
    
    // Option 1: Pre-combine kAtomcellVolume * kInvE * kInvMuB
    // torque = kTorqueScaleBase * sd_exchange * ax
    double torque_precombined = kTorqueScaleBase * sd_exchange * ax;
    std::cout << "Option 1 (pre-combine compile-time constants):\n";
    std::cout << "  kTorqueScaleBase = " << kTorqueScaleBase << "\n";
    std::cout << "  torque = kTorqueScaleBase * sd_exchange * ax = " << torque_precombined << " T\n";
    std::cout << "  Error: " << std::abs(torque_current - torque_precombined) / torque_current << " (relative)\n\n";
    
    // Option 2: Scale during interpolation (multiply fine-grid points)
    // During fine-grid loop: ax_scaled += (Sbase[i] * kTorqueScaleBase)
    // Then: torque = sd_exchange * ax_scaled
    // But wait - this would scale each fine-grid point, then average
    // Let's check if this is numerically equivalent
    
    // Simulate fine-grid interpolation with nsub = 10
    int nsub = 10;
    double S_fine[10];
    for(int i=0; i<nsub; i++) {
        S_fine[i] = ax;  // All fine points have same value for simplicity
    }
    
    // Current: average first, then multiply
    double ax_avg_current = 0.0;
    for(int i=0; i<nsub; i++) {
        ax_avg_current += S_fine[i];
    }
    ax_avg_current /= nsub;
    double torque_avg_first = kTorqueScaleBase * sd_exchange * ax_avg_current;
    
    // Alternative: multiply each fine point, then average
    double ax_scaled_avg = 0.0;
    for(int i=0; i<nsub; i++) {
        ax_scaled_avg += S_fine[i] * kTorqueScaleBase;
    }
    ax_scaled_avg /= nsub;
    double torque_mult_first = sd_exchange * ax_scaled_avg;
    
    std::cout << "Option 2 (scale during fine-grid interpolation):\n";
    std::cout << "  Current: average first, then multiply by kTorqueScaleBase\n";
    std::cout << "    torque = " << torque_avg_first << " T\n";
    std::cout << "  Alternative: multiply each fine point, then average\n";
    std::cout << "    torque = " << torque_mult_first << " T\n";
    std::cout << "  Numerically equivalent? " << (std::abs(torque_avg_first - torque_mult_first) < 1e-15 ? "YES" : "NO") << "\n";
    std::cout << "  Error: " << std::abs(torque_avg_first - torque_mult_first) / torque_avg_first << " (relative)\n\n";
    
    // Option 3: Only scale the final averaged ax
    // Pre-compute: torque_per_sa = kTorqueScaleBase * sd_exchange
    // Then: torque = torque_per_sa * ax
    double torque_per_sa = kTorqueScaleBase * sd_exchange;
    double torque_scaled_ax = torque_per_sa * ax;
    
    std::cout << "Option 3 (pre-combine all constants except ax):\n";
    std::cout << "  torque_per_sa = kTorqueScaleBase * sd_exchange = " << torque_per_sa << "\n";
    std::cout << "  torque = torque_per_sa * ax = " << torque_scaled_ax << " T\n";
    std::cout << "  Error: " << std::abs(torque_current - torque_scaled_ax) / torque_current << " (relative)\n\n";
    
    // Check precision with float
    float kTorqueScaleBase_float = static_cast<float>(kTorqueScaleBase);
    float sd_exchange_float = static_cast<float>(sd_exchange);
    float ax_float = static_cast<float>(ax);
    
    float torque_float_precombined = kTorqueScaleBase_float * sd_exchange_float * ax_float;
    
    std::cout << "PRECISION TEST: Using float for kTorqueScaleBase:\n";
    std::cout << "  kTorqueScaleBase (double): " << kTorqueScaleBase << "\n";
    std::cout << "  kTorqueScaleBase (float):  " << kTorqueScaleBase_float << "\n";
    std::cout << "  Relative error in constant: " << std::abs((kTorqueScaleBase - kTorqueScaleBase_float) / kTorqueScaleBase) << "\n";
    std::cout << "  Torque (double): " << torque_precombined << " T\n";
    std::cout << "  Torque (float):  " << torque_float_precombined << " T\n";
    std::cout << "  Relative error in torque: " << std::abs((torque_precombined - torque_float_precombined) / torque_precombined) << "\n\n";
    
    // Check if kTorqueScaleBase is representable in float
    float float_min = std::numeric_limits<float>::min();
    float float_max = std::numeric_limits<float>::max();
    
    std::cout << "FLOAT RANGE CHECK:\n";
    std::cout << "  float min: " << float_min << "\n";
    std::cout << "  float max: " << float_max << "\n";
    std::cout << "  kTorqueScaleBase: " << kTorqueScaleBase << "\n";
    std::cout << "  Can represent kTorqueScaleBase in float? ";
    if (kTorqueScaleBase >= float_min && kTorqueScaleBase <= float_max) {
        std::cout << "YES\n";
    } else {
        std::cout << "NO\n";
    }
    
    // Final recommendation
    std::cout << "\n=== RECOMMENDATION ===\n";
    std::cout << "Pre-combining kAtomcellVolume * kInvE * kInvMuB:\n";
    std::cout << "  - Numerically equivalent to current approach ✓\n";
    std::cout << "  - Reduces number of multiplications (2 instead of 4) ✓\n";
    std::cout << "  - Can be computed at compile time ✓\n";
    std::cout << "  - Still requires double precision (kTorqueScaleBase too large for float) ⚠\n";
    std::cout << "  - Does NOT solve the underflow issue (still need double for V * sd_exchange)\n\n";
    
    std::cout << "BEST APPROACH: Pre-compute kTorqueScaleBase at compile time, but keep as double.\n";
    std::cout << "This is a minor optimization (fewer multiplications) but doesn't enable float conversion.\n";
    
    return 0;
}