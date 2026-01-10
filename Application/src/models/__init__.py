"""
Mathematical models for meta-material damping.

Modules:
- parameters: System parameter dataclasses (PC, AM)
- consistent_unit_cell: Matrix builders for consistent unit cells
- dispersion: Dispersion relation solvers
"""

from .parameters import (
    MonoatomicParameters,
    PhononicCrystalParameters,
    AcousticMetamaterialParameters,
    get_benchmark_pc_parameters,
    get_benchmark_am_parameters,
    validate_equal_sound_speed,
)

from .consistent_unit_cell import (
    SystemMatrices,
    build_monoatomic_matrices,
    build_phononic_crystal_matrices,
    build_acoustic_metamaterial_matrices,
    compute_undamped_frequencies,
    compute_modal_matrices,
    compute_modal_damping_ratios,
)

from .dispersion import (
    monoatomic_dispersion,
    phononic_crystal_dispersion,
    acoustic_metamaterial_dispersion,
    compute_bandgap,
    compute_group_velocity,
)

__all__ = [
    # Parameters
    "MonoatomicParameters",
    "PhononicCrystalParameters", 
    "AcousticMetamaterialParameters",
    "get_benchmark_pc_parameters",
    "get_benchmark_am_parameters",
    "validate_equal_sound_speed",
    # Matrices
    "SystemMatrices",
    "build_monoatomic_matrices",
    "build_phononic_crystal_matrices",
    "build_acoustic_metamaterial_matrices",
    "compute_undamped_frequencies",
    "compute_modal_matrices",
    "compute_modal_damping_ratios",
    # Dispersion
    "monoatomic_dispersion",
    "phononic_crystal_dispersion",
    "acoustic_metamaterial_dispersion",
    "compute_bandgap",
    "compute_group_velocity",
]
