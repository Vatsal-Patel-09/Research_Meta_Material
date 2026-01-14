"""
Main simulation script for Meta-dissipation analysis.

This script reproduces the key results from the paper:
- Θ_PC ≈ 7.49 (Phononic Crystal)
- Θ_AM ≈ 14.52 (Acoustic Metamaterial)

Reference: Meta-dissipation paper (Banerjee, Bera, Adhikari, 2025)
"""

import numpy as np
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.models import (
    get_benchmark_pc_parameters,
    get_benchmark_am_parameters,
    validate_equal_sound_speed,
    build_phononic_crystal_matrices,
    build_acoustic_metamaterial_matrices,
    phononic_crystal_dispersion,
    acoustic_metamaterial_dispersion,
)

from src.core import (
    compute_meta_dissipation,
    impulse_response_modal,
    compute_energy_from_displacement,
)

from src.utils import (
    plot_dispersion_comparison,
    plot_response_comparison,
    plot_energy_decay,
    plot_theta_comparison,
    plot_modal_parameters,
    create_summary_figure,
)


def print_section(title: str):
    """Print formatted section header."""
    print("\n" + "=" * 60)
    print(f" {title}")
    print("=" * 60)


def run_validation():
    """
    Run complete Meta-dissipation validation.
    
    Reproduces:
    - Θ_PC ≈ 7.49
    - Θ_AM ≈ 14.52
    """
    print_section("META-DISSIPATION FRAMEWORK VALIDATION")
    print("Reference: Banerjee, Bera, Adhikari (2025)")
    
    # =========================================================================
    # 1. Load Benchmark Parameters
    # =========================================================================
    print_section("1. BENCHMARK PARAMETERS")
    
    pc_params = get_benchmark_pc_parameters()
    am_params = get_benchmark_am_parameters()
    
    print("\nPhononic Crystal (PC) Parameters:")
    print(f"  m₁ = {pc_params.m1} kg")
    print(f"  m₂ = {pc_params.m2} kg")
    print(f"  k₁ = {pc_params.k1} N/m")
    print(f"  k₂ = {pc_params.k2} N/m")
    print(f"  c₁ = {pc_params.c1} Ns/m")
    print(f"  c₂ = {pc_params.c2} Ns/m")
    print(f"  C_stat = {pc_params.static_sound_speed:.2f} m/s")
    
    print("\nAcoustic Metamaterial (AM) Parameters:")
    print(f"  m₁ = {am_params.m1} kg")
    print(f"  m₂ = {am_params.m2} kg")
    print(f"  k₁ = {am_params.k1} N/m")
    print(f"  k₂ = {am_params.k2} N/m")
    print(f"  c₁ = {am_params.c1} Ns/m")
    print(f"  c₂ = {am_params.c2} Ns/m")
    print(f"  C_stat = {am_params.static_sound_speed:.2f} m/s")
    
    # Validate equal sound speeds
    if validate_equal_sound_speed(pc_params, am_params):
        print("\n✓ Sound speeds match - fair comparison ensured")
    else:
        print("\n⚠ Warning: Sound speeds differ")
    
    # =========================================================================
    # 2. Build Consistent Unit Cell Matrices
    # =========================================================================
    print_section("2. CONSISTENT UNIT CELL MATRICES")
    
    matrices_pc = build_phononic_crystal_matrices(pc_params)
    matrices_am = build_acoustic_metamaterial_matrices(am_params)
    
    print("\nPC Matrices (Eq S15):")
    print(f"  M = {matrices_pc.M}")
    print(f"  C = {matrices_pc.C}")
    print(f"  K = {matrices_pc.K}")
    
    print("\nAM Matrices (Eq S16):")
    print(f"  M = {matrices_am.M}")
    print(f"  C = {matrices_am.C}")
    print(f"  K = {matrices_am.K}")
    
    # =========================================================================
    # 3. Compute Meta-dissipation
    # =========================================================================
    print_section("3. META-DISSIPATION ANALYSIS")
    
    # Prepare parameter dictionaries for alpha weight calculation
    pc_param_dict = {
        'm1': pc_params.m1,
        'm2': pc_params.m2,
        'k1': pc_params.k1,
        'k2': pc_params.k2,
        'c1': pc_params.c1,
        'c2': pc_params.c2,
    }
    
    am_param_dict = {
        'm1': am_params.m1,
        'm2': am_params.m2,
        'k1': am_params.k1,
        'k2': am_params.k2,
        'c1': am_params.c1,
        'c2': am_params.c2,
    }
    
    result_pc = compute_meta_dissipation(matrices_pc, system_type='PC', params=pc_param_dict)
    result_am = compute_meta_dissipation(matrices_am, system_type='AM', params=am_param_dict)
    
    print("\nPhononic Crystal Results:")
    print(f"  Natural frequencies: ω₁={result_pc.modal_params.omega[0]:.2f}, "
          f"ω₂={result_pc.modal_params.omega[1]:.2f} rad/s")
    print(f"  Damping ratios: ξ₁={result_pc.modal_params.xi[0]:.4f}, "
          f"ξ₂={result_pc.modal_params.xi[1]:.4f}")
    print(f"  Decay coefficients: θᵤ={result_pc.theta_u:.2f}, "
          f"θᵥ={result_pc.theta_v:.2f}")
    print(f"  → Unified Θ_PC = {result_pc.Theta:.2f}")
    
    print("\nAcoustic Metamaterial Results:")
    print(f"  Natural frequencies: ω₁={result_am.modal_params.omega[0]:.2f}, "
          f"ω₂={result_am.modal_params.omega[1]:.2f} rad/s")
    print(f"  Damping ratios: ξ₁={result_am.modal_params.xi[0]:.4f}, "
          f"ξ₂={result_am.modal_params.xi[1]:.4f}")
    print(f"  Decay coefficients: θᵤ={result_am.theta_u:.2f}, "
          f"θᵥ={result_am.theta_v:.2f}")
    print(f"  → Unified Θ_AM = {result_am.Theta:.2f}")
    
    # =========================================================================
    # 4. Comparison with Paper Values
    # =========================================================================
    print_section("4. VALIDATION AGAINST PAPER")
    
    Theta_PC_paper = 7.49
    Theta_AM_paper = 14.52
    
    error_pc = abs(result_pc.Theta - Theta_PC_paper) / Theta_PC_paper * 100
    error_am = abs(result_am.Theta - Theta_AM_paper) / Theta_AM_paper * 100
    
    print(f"\nPaper values:")
    print(f"  Θ_PC (paper) = {Theta_PC_paper}")
    print(f"  Θ_AM (paper) = {Theta_AM_paper}")
    
    print(f"\nComputed values:")
    print(f"  Θ_PC (computed) = {result_pc.Theta:.2f}  (error: {error_pc:.1f}%)")
    print(f"  Θ_AM (computed) = {result_am.Theta:.2f}  (error: {error_am:.1f}%)")
    
    improvement = (result_am.Theta / result_pc.Theta - 1) * 100
    print(f"\n→ AM shows {improvement:.0f}% higher dissipation than PC")
    
    if error_pc < 10 and error_am < 10:
        print("\n✓ VALIDATION SUCCESSFUL: Results match paper within 10%")
    else:
        print("\n⚠ Results deviate from paper - check implementation")
    
    return result_pc, result_am, pc_params, am_params, matrices_pc, matrices_am


def run_time_domain_simulation(result_pc, result_am, matrices_pc, matrices_am):
    """Run time-domain simulations and create plots."""
    
    print_section("5. TIME-DOMAIN SIMULATION")
    
    # Time vector
    t = np.linspace(0, 0.5, 5000)
    
    # Impulse response
    u_pc, v_pc = impulse_response_modal(matrices_pc, t)
    u_am, v_am = impulse_response_modal(matrices_am, t)
    
    print(f"\nPC max displacement: {np.max(np.abs(u_pc)):.4f}")
    print(f"AM max displacement: {np.max(np.abs(u_am)):.4f}")
    
    # Energy calculation
    _, _, E_pc = compute_energy_from_displacement(matrices_pc, t, u_pc, v_pc)
    _, _, E_am = compute_energy_from_displacement(matrices_am, t, u_am, v_am)
    
    return t, u_pc, v_pc, u_am, v_am, E_pc, E_am


def create_all_figures(t, u_pc, v_pc, u_am, v_am, E_pc, E_am, 
                       result_pc, result_am, pc_params, am_params):
    """Create and save all figures."""
    
    print_section("6. GENERATING FIGURES")
    
    import matplotlib.pyplot as plt
    
    # Figure 1: Dispersion comparison
    mu = np.linspace(0, np.pi, 200)
    omega_pc = phononic_crystal_dispersion(pc_params, mu)
    omega_am = acoustic_metamaterial_dispersion(am_params, mu)
    
    fig1, _ = plot_dispersion_comparison(mu, omega_pc, omega_am)
    fig1.savefig('dispersion_comparison.png', dpi=150, bbox_inches='tight')
    print("  Saved: dispersion_comparison.png")
    
    # Figure 2: Response comparison
    fig2, _ = plot_response_comparison(t, u_pc, u_am)
    fig2.savefig('response_comparison.png', dpi=150, bbox_inches='tight')
    print("  Saved: response_comparison.png")
    
    # Figure 3: Energy decay
    fig3, _ = plot_energy_decay(t, E_pc, E_am, result_pc.Theta, result_am.Theta)
    fig3.savefig('energy_decay.png', dpi=150, bbox_inches='tight')
    print("  Saved: energy_decay.png")
    
    # Figure 4: Theta comparison
    fig4, _ = plot_theta_comparison(result_pc.Theta, result_am.Theta)
    fig4.savefig('theta_comparison.png', dpi=150, bbox_inches='tight')
    print("  Saved: theta_comparison.png")
    
    # Figure 5: Modal parameters
    fig5, _ = plot_modal_parameters(
        result_pc.modal_params.omega, result_pc.modal_params.xi,
        result_am.modal_params.omega, result_am.modal_params.xi
    )
    fig5.savefig('modal_parameters.png', dpi=150, bbox_inches='tight')
    print("  Saved: modal_parameters.png")
    
    # Figure 6: Summary
    fig6, _ = create_summary_figure(
        t, u_pc, u_am, E_pc, E_am,
        result_pc.Theta, result_am.Theta,
        title="Meta-dissipation Analysis: PC vs AM"
    )
    fig6.savefig('summary_figure.png', dpi=150, bbox_inches='tight')
    print("  Saved: summary_figure.png")
    
    plt.close('all')
    print("\n✓ All figures generated successfully")


def main():
    """Main entry point."""
    # Run validation
    result_pc, result_am, pc_params, am_params, matrices_pc, matrices_am = run_validation()
    
    # Time domain simulation
    t, u_pc, v_pc, u_am, v_am, E_pc, E_am = run_time_domain_simulation(
        result_pc, result_am, matrices_pc, matrices_am
    )
    
    # Generate figures
    try:
        create_all_figures(
            t, u_pc, v_pc, u_am, v_am, E_pc, E_am,
            result_pc, result_am, pc_params, am_params
        )
    except Exception as e:
        print(f"\n⚠ Figure generation failed: {e}")
        print("  (This may be expected in non-GUI environments)")
    
    print_section("SIMULATION COMPLETE")
    print("\nKey Results:")
    print(f"  Θ_PC = {result_pc.Theta:.2f} (expected: ~7.49)")
    print(f"  Θ_AM = {result_am.Theta:.2f} (expected: ~14.52)")
    
    return result_pc, result_am


if __name__ == "__main__":
    main()
