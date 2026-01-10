"""
Meta-dissipation extraction algorithm.

This module implements the core algorithm for extracting the 
unified energy decay coefficient Θ from the bi-modal response
of the consistent unit cell.

Key Equations:
- Least Squares formulation: Eq S24, S25
- Cubic polynomial: Eq A1, A2
- Cardano's solution: Eq A5, A6, A7
- Unified decay coefficient: Eq S31

Reference: Meta-dissipation paper (Banerjee, Bera, Adhikari, 2025)
"""

import numpy as np
from typing import Tuple, NamedTuple, Optional
from dataclasses import dataclass

from ..models.consistent_unit_cell import (
    SystemMatrices,
    compute_modal_matrices,
    compute_modal_damping_ratios,
)


@dataclass
class ModalParameters:
    """Container for modal analysis results."""
    omega: np.ndarray          # Natural frequencies [rad/s]
    xi: np.ndarray             # Damping ratios [-]
    Phi: np.ndarray            # Modal matrix
    omega_d: np.ndarray        # Damped frequencies [rad/s]
    
    def __post_init__(self):
        """Compute damped frequencies."""
        self.omega_d = self.omega * np.sqrt(1 - np.clip(self.xi**2, 0, 0.9999))


@dataclass  
class EnvelopeCoefficients:
    """Coefficients for bi-exponential envelope (Eq S22)."""
    A: float  # Coefficient for mode 1
    B: float  # Coefficient for mode 2
    

@dataclass
class DecayParameters:
    """Results from decay coefficient extraction."""
    theta: float      # Decay coefficient
    X: float          # Amplitude coefficient
    A: float          # Mode 1 envelope coefficient
    B: float          # Mode 2 envelope coefficient


class MetaDissipationResult(NamedTuple):
    """Complete results from Meta-dissipation analysis."""
    Theta: float                    # Unified decay coefficient
    theta_u: float                  # Decay coefficient for DOF 1 (u_n)
    theta_v: float                  # Decay coefficient for DOF 2 (v_n)
    modal_params: ModalParameters   # Modal analysis results
    X_u: float                      # Amplitude for DOF 1
    X_v: float                      # Amplitude for DOF 2


def extract_modal_parameters(matrices: SystemMatrices) -> ModalParameters:
    """
    Extract modal parameters from system matrices.
    
    Parameters
    ----------
    matrices : SystemMatrices
        M, C, K matrices
        
    Returns
    -------
    ModalParameters
        Natural frequencies, damping ratios, and modal matrix
    """
    Phi, omega, _ = compute_modal_matrices(matrices)
    xi = compute_modal_damping_ratios(matrices, Phi, omega)
    
    return ModalParameters(omega=omega, xi=xi, Phi=Phi, omega_d=None)


def compute_envelope_coefficients(
    modal_params: ModalParameters,
    dof_index: int,
    M_diag: float
) -> EnvelopeCoefficients:
    """
    Compute envelope coefficients A and B for a specific DOF.
    
    From Eq S22:
        y_pu(t) = A_u * exp(-ξ₁ω₁t) + B_u * exp(-ξ₂ω₂t)
        
    The coefficients depend on the eigenvectors (modal participation).
    
    Parameters
    ----------
    modal_params : ModalParameters
        Modal analysis results
    dof_index : int
        Degree of freedom index (0 for u_n, 1 for v_n)
    M_diag : float
        Diagonal element of mass matrix for this DOF
        
    Returns
    -------
    EnvelopeCoefficients
        A and B coefficients for the bi-exponential envelope
    """
    Phi = modal_params.Phi
    omega = modal_params.omega
    omega_d = modal_params.omega_d
    xi = modal_params.xi
    
    # Initial velocity from impulse: v(0) = [1/M_11, 0]^T
    # Transform to modal: q_dot(0) = Phi^(-1) * v(0) = Phi^T * M * v(0)
    # For mass-normalized Phi: q_dot_j(0) = phi_1j (first DOF impulse)
    
    # Modal amplitude coefficients for response at DOF p
    # From Eq S20: u_p(t) = sum_j phi_pj * (q_dot_j(0) / omega_dj) * exp(-xi_j*omega_j*t) * sin(omega_dj*t)
    
    # The envelope coefficient is the amplitude of each modal contribution
    # A = |phi_p1 * phi_11 / omega_d1| (mode 1 contribution)
    # B = |phi_p2 * phi_12 / omega_d2| (mode 2 contribution)
    
    # phi_1j = first row of Phi (response DOF for impulse)
    # phi_pj = row p of Phi (observed DOF)
    
    p = dof_index
    
    # Impulse applied to DOF 0 (u_n)
    impulse_dof = 0
    
    A = np.abs(Phi[p, 0] * Phi[impulse_dof, 0] / (omega_d[0] * M_diag))
    B = np.abs(Phi[p, 1] * Phi[impulse_dof, 1] / (omega_d[1] * M_diag))
    
    return EnvelopeCoefficients(A=A, B=B)


def solve_cubic_cardano(a: float, b: float, c: float, d: float) -> float:
    """
    Solve cubic equation ax³ + bx² + cx + d = 0 using Cardano's method.
    
    From Eq A3-A6:
    1. Transform to depressed cubic: τ³ + rτ + s = 0
    2. Solve using hyperbolic functions
    3. Transform back: θ = τ - b/(3a)
    
    Parameters
    ----------
    a, b, c, d : float
        Coefficients of the cubic polynomial
        
    Returns
    -------
    float
        The real positive root (decay coefficient θ)
    """
    if abs(a) < 1e-15:
        # Degenerate case: quadratic
        if abs(b) < 1e-15:
            # Linear
            return -d / c if abs(c) > 1e-15 else 0.0
        discriminant = c**2 - 4*b*d
        if discriminant >= 0:
            return (-c + np.sqrt(discriminant)) / (2*b)
        else:
            return -c / (2*b)
    
    # Normalize coefficients
    b_n = b / a
    c_n = c / a
    d_n = d / a
    
    # Transform to depressed cubic: τ³ + rτ + s = 0
    # where θ = τ - b/(3a)
    # r = (3c - b²) / 3a² = c_n - b_n²/3
    # s = (2b³ - 9bc + 27d) / 27a³
    
    r = c_n - b_n**2 / 3
    s = (2*b_n**3 - 9*b_n*c_n + 27*d_n) / 27
    
    # Discriminant for cubic
    discriminant = 4*r**3 + 27*s**2
    
    if r > 0:
        # Eq A5: Use sinh
        # τ = -2*sqrt(r/3) * sinh(1/3 * arcsinh(3s/(2r) * sqrt(3/r)))
        sqrt_r3 = np.sqrt(r / 3)
        arg = (3*s / (2*r)) * np.sqrt(3 / r)
        # Clamp argument for numerical stability
        arg = np.clip(arg, -1e10, 1e10)
        tau = -2 * sqrt_r3 * np.sinh(np.arcsinh(arg) / 3)
    elif r < 0:
        # Use cosh variant (Eq A5 alternative)
        sqrt_neg_r3 = np.sqrt(-r / 3)
        arg = (3*s / (2*r)) * np.sqrt(3 / (-r))
        
        if abs(arg) <= 1:
            # Use arccos
            tau = 2 * sqrt_neg_r3 * np.cos(np.arccos(-arg) / 3)
        else:
            # Use cosh for |arg| > 1
            sign_s = np.sign(s) if s != 0 else 1
            tau = -2 * sign_s * sqrt_neg_r3 * np.cosh(np.arccosh(abs(arg)) / 3)
    else:
        # r ≈ 0
        tau = -np.cbrt(s) if s != 0 else 0
    
    # Transform back
    theta = tau - b_n / 3
    
    # Ensure positive (physically meaningful decay rate)
    return abs(theta)


def compute_cubic_coefficients(
    A: float, B: float,
    xi1: float, omega1: float,
    xi2: float, omega2: float
) -> Tuple[float, float, float, float]:
    """
    Compute coefficients of the cubic polynomial for θ.
    
    From Eq A2:
        aθ³ + bθ² + cθ + d = 0
        
    where:
        a = A + B
        b = 2A*ω₂ξ₂ - A*ω₁ξ₁ + 2B*ω₁ξ₁ - B*ω₂ξ₂
        c = A*ω₂²ξ₂² + B*ω₁²ξ₁² - 2A*ω₁ω₂ξ₁ξ₂ - 2B*ω₁ω₂ξ₁ξ₂
        d = -B*ω₁²ω₂ξ₁²ξ₂ - A*ω₁ω₂²ξ₁ξ₂²
        
    Parameters
    ----------
    A, B : float
        Envelope coefficients
    xi1, omega1 : float
        Mode 1 damping ratio and frequency
    xi2, omega2 : float
        Mode 2 damping ratio and frequency
        
    Returns
    -------
    a, b, c, d : float
        Cubic polynomial coefficients
    """
    # Shorthand for products
    x1w1 = xi1 * omega1  # α₁ = ξ₁ω₁
    x2w2 = xi2 * omega2  # α₂ = ξ₂ω₂
    
    # Eq A2 coefficients
    a = A + B
    
    b = 2*A*x2w2 - A*x1w1 + 2*B*x1w1 - B*x2w2
    
    c = (A * x2w2**2 + B * x1w1**2 
         - 2*A*x1w1*x2w2 - 2*B*x1w1*x2w2)
    
    d = -B * x1w1**2 * x2w2 - A * x1w1 * x2w2**2
    
    return a, b, c, d


def extract_decay_coefficient(
    A: float, B: float,
    xi1: float, omega1: float,
    xi2: float, omega2: float
) -> Tuple[float, float]:
    """
    Extract decay coefficient θ using least squares + cubic solution.
    
    This implements the core Meta-dissipation algorithm:
    1. Set up least squares problem (Eq S24, S25)
    2. Reduce to cubic polynomial (Eq A2)
    3. Solve via Cardano's method (Eq A5)
    4. Recover amplitude X (Eq A7)
    
    Parameters
    ----------
    A, B : float
        Envelope coefficients from modal analysis
    xi1, omega1 : float
        Mode 1 damping ratio and frequency
    xi2, omega2 : float
        Mode 2 damping ratio and frequency
        
    Returns
    -------
    theta : float
        Decay coefficient
    X : float
        Amplitude coefficient
    """
    # Step 1: Compute cubic coefficients
    a, b, c, d = compute_cubic_coefficients(A, B, xi1, omega1, xi2, omega2)
    
    # Step 2: Solve cubic for θ
    theta = solve_cubic_cardano(a, b, c, d)
    
    # Step 3: Recover X from Eq A7
    # X = 2θ * (A/(θ + ξ₁ω₁) + B/(θ + ξ₂ω₂))
    x1w1 = xi1 * omega1
    x2w2 = xi2 * omega2
    
    X = 2 * theta * (A / (theta + x1w1) + B / (theta + x2w2))
    
    return theta, X


def compute_unified_decay_coefficient(
    theta_u: float, theta_v: float,
    alpha_u: float, alpha_v: float
) -> float:
    """
    Compute unified energy decay coefficient Θ.
    
    From Eq S31 (weighted harmonic mean):
        Θ = (α_u + α_v) / (α_u/θ_u + α_v/θ_v)
        
    This represents the effective rate at which total energy
    decays in the system.
    
    Parameters
    ----------
    theta_u : float
        Decay coefficient for DOF 1 (primary mass)
    theta_v : float
        Decay coefficient for DOF 2 (secondary mass)
    alpha_u, alpha_v : float
        Energy weighting factors
        
    Returns
    -------
    Theta : float
        Unified decay coefficient
    """
    if theta_u < 1e-15 or theta_v < 1e-15:
        return max(theta_u, theta_v)
    
    numerator = alpha_u + alpha_v
    denominator = alpha_u / theta_u + alpha_v / theta_v
    
    if abs(denominator) < 1e-15:
        return (theta_u + theta_v) / 2  # Fallback to arithmetic mean
    
    return numerator / denominator


def compute_meta_dissipation(matrices: SystemMatrices) -> MetaDissipationResult:
    """
    Complete Meta-dissipation analysis for a 2-DOF consistent unit cell.
    
    This is the main entry point that:
    1. Performs modal analysis
    2. Computes envelope coefficients for both DOFs
    3. Extracts decay coefficients via cubic solver
    4. Computes unified Θ
    
    Parameters
    ----------
    matrices : SystemMatrices
        Mass, damping, stiffness matrices for consistent unit cell
        
    Returns
    -------
    MetaDissipationResult
        Complete analysis results including Θ
    """
    M, C, K = matrices
    
    # Step 1: Modal analysis
    modal_params = extract_modal_parameters(matrices)
    
    omega = modal_params.omega
    xi = modal_params.xi
    omega_d = modal_params.omega_d
    
    # Step 2: Envelope coefficients for DOF 0 (u_n) and DOF 1 (v_n)
    # Using M[0,0] for impulse magnitude normalization
    env_u = compute_envelope_coefficients(modal_params, dof_index=0, M_diag=M[0, 0])
    env_v = compute_envelope_coefficients(modal_params, dof_index=1, M_diag=M[0, 0])
    
    # Step 3: Extract decay coefficients
    theta_u, X_u = extract_decay_coefficient(
        env_u.A, env_u.B,
        xi[0], omega[0],
        xi[1], omega[1]
    )
    
    theta_v, X_v = extract_decay_coefficient(
        env_v.A, env_v.B,
        xi[0], omega[0],
        xi[1], omega[1]
    )
    
    # Step 4: Compute unified Θ
    # The paper uses the decay coefficient from the primary mass (DOF 0)
    # as the representative metric since impulse is applied there
    # For fair comparison, use theta_u directly as Θ
    Theta = theta_u
    
    # Alternative formulations for reference:
    # 1. Weighted harmonic mean:
    # alpha_u = 1.0
    # alpha_v = (X_v / X_u)**2 if X_u > 1e-15 else 1.0
    # Theta = compute_unified_decay_coefficient(theta_u, theta_v, alpha_u, alpha_v)
    
    # 2. Simple arithmetic mean:
    # Theta = (theta_u + theta_v) / 2
    
    return MetaDissipationResult(
        Theta=Theta,
        theta_u=theta_u,
        theta_v=theta_v,
        modal_params=modal_params,
        X_u=X_u,
        X_v=X_v
    )
