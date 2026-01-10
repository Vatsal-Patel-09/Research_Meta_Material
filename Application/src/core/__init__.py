"""
Core functionality for meta-material damping analysis.

Modules:
- meta_dissipation: Main algorithm for extracting Θ
- impulse_response: Time-domain simulation
"""

from .meta_dissipation import (
    ModalParameters,
    EnvelopeCoefficients,
    DecayParameters,
    MetaDissipationResult,
    extract_modal_parameters,
    compute_envelope_coefficients,
    solve_cubic_cardano,
    compute_cubic_coefficients,
    extract_decay_coefficient,
    compute_unified_decay_coefficient,
    compute_meta_dissipation,
)

from .impulse_response import (
    impulse_response_modal,
    impulse_response_numerical,
    compute_envelope,
    biexponential_envelope,
    single_exponential_fit,
    compute_energy,
    compute_energy_from_displacement,
)

__all__ = [
    # Meta-dissipation
    "ModalParameters",
    "EnvelopeCoefficients",
    "DecayParameters",
    "MetaDissipationResult",
    "extract_modal_parameters",
    "compute_envelope_coefficients",
    "solve_cubic_cardano",
    "compute_cubic_coefficients",
    "extract_decay_coefficient",
    "compute_unified_decay_coefficient",
    "compute_meta_dissipation",
    # Impulse response
    "impulse_response_modal",
    "impulse_response_numerical",
    "compute_envelope",
    "biexponential_envelope",
    "single_exponential_fit",
    "compute_energy",
    "compute_energy_from_displacement",
]
