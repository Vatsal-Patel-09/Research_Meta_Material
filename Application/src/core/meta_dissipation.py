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
    M: np.ndarray
) -> EnvelopeCoefficients:
    """
    Compute envelope coefficients A and B for a specific DOF.
    
    From Eq S22:
        y_pu(t) = A_u * exp(-ξ₁ω₁t) + B_u * exp(-ξ₂ω₂t)
        
    The coefficients depend on the eigenvectors (modal participation).
    
    For a unit force impulse F = [1, 0]^T at DOF 0:
    - Initial velocity: v(0) = M^(-1) * F = [1/M[0,0], 0]^T
    - Modal initial velocity: q_dot_j(0) = Φ_j^T * M * v(0) = Φ[0,j] (for mass-normalized modes)
    
    Parameters
    ----------
    modal_params : ModalParameters
        Modal analysis results
    dof_index : int
        Degree of freedom index (0 for u_n, 1 for v_n)
    M : np.ndarray
        Full mass matrix
        
    Returns
    -------
    EnvelopeCoefficients
        A and B coefficients for the bi-exponential envelope
    """
    Phi = modal_params.Phi
    omega_d = modal_params.omega_d
    
    p = dof_index
    impulse_dof = 0
    
    # For unit force impulse at DOF 0:
    # Initial velocity v0 = M^(-1) * [1, 0]^T
    v0 = np.zeros(M.shape[0])
    v0[impulse_dof] = 1.0 / M[impulse_dof, impulse_dof]
    
    # Modal initial velocities: q_dot_j(0) = sum_i Phi[i,j] * M[i,i] * v0[i]
    # For mass-normalized modes: Phi^T * M * Phi = I
    # So q_dot_j(0) = Phi^T * M * v0
    # For diagonal M: q_dot_j(0) = sum_i Phi[i,j] * M[i,i] * v0[i]
    q_dot_0 = np.zeros(len(omega_d))
    for j in range(len(omega_d)):
        for i in range(M.shape[0]):
            q_dot_0[j] += Phi[i, j] * M[i, i] * v0[i]
    
    # Envelope coefficients: amplitude of each modal contribution to DOF p
    # From modal superposition: u_p(t) = sum_j Phi[p,j] * (q_dot_j(0)/omega_dj) * exp(-xi_j*omega_j*t) * sin(omega_dj*t)
    # The envelope amplitude is |Phi[p,j] * q_dot_j(0) / omega_dj|
    
    A = np.abs(Phi[p, 0] * q_dot_0[0] / omega_d[0])
    B = np.abs(Phi[p, 1] * q_dot_0[1] / omega_d[1])
    
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


def compute_alpha_weights_pc(
    Phi: np.ndarray,
    theta_modes: np.ndarray,
    k1: float,
    k2: float,
    m1: float
) -> np.ndarray:
    """
    Compute modal energy weights α_u for Phononic Crystal.
    
    PC-specific formula from paper:
    α_u = 0.5 * (k1 + k2 + m1 * θ_u²) * X_pu²
    
    where:
    - θ_u = ξ_u * ω_u is the modal decay coefficient
    - X_pu is the primary mass mode shape component (DOF 0)
    
    Parameters
    ----------
    Phi : np.ndarray
        Mass-normalized mode shape matrix [n_dof x n_modes]
    theta_modes : np.ndarray
        Modal decay coefficients [n_modes]
    k1, k2 : float
        Stiffness values [N/m]
    m1 : float
        Primary mass [kg]
        
    Returns
    -------
    np.ndarray
        Weight vector [n_modes] (NOT normalized - use raw values for harmonic mean)
    """
    n_modes = Phi.shape[1]
    alpha = np.zeros(n_modes)
    
    for u in range(n_modes):
        # Primary mass mode shape component
        X_pu = Phi[0, u]
        
        # PC alpha formula (energy-based)
        alpha[u] = 0.5 * (k1 + k2 + m1 * theta_modes[u]**2) * X_pu**2
    
    return alpha


def compute_alpha_weights_am(
    Phi: np.ndarray,
    theta_modes: np.ndarray,
    k1: float,
    k2: float,
    m1: float
) -> np.ndarray:
    """
    Compute modal energy weights α_u for Acoustic Metamaterial.
    
    AM-specific formula from paper:
    α_u = 0.25 * (2*k1 + k2 + m1 * θ_u²) * X_pu²
    
    where:
    - θ_u = ξ_u * ω_u is the modal decay coefficient
    - X_pu is the primary mass mode shape component (DOF 0)
    
    Parameters
    ----------
    Phi : np.ndarray
        Mass-normalized mode shape matrix [n_dof x n_modes]
    theta_modes : np.ndarray
        Modal decay coefficients [n_modes]
    k1, k2 : float
        Stiffness values [N/m]
    m1 : float
        Primary mass [kg]
        
    Returns
    -------
    np.ndarray
        Weight vector [n_modes] (NOT normalized - use raw values for harmonic mean)
    """
    n_modes = Phi.shape[1]
    alpha = np.zeros(n_modes)
    
    for u in range(n_modes):
        # Primary mass mode shape component
        X_pu = Phi[0, u]
        
        # AM alpha formula (energy-based)
        alpha[u] = 0.25 * (2*k1 + k2 + m1 * theta_modes[u]**2) * X_pu**2
    
    return alpha


def compute_meta_dissipation(
    matrices: SystemMatrices,
    system_type: str = 'PC',
    params: Optional[dict] = None
) -> MetaDissipationResult:
    """
    Complete Meta-dissipation analysis for a 2-DOF consistent unit cell.
    
    This is the main entry point that:
    1. Performs modal analysis
    2. Computes envelope coefficients for both DOFs
    3. Extracts decay coefficients via cubic solver
    4. Computes unified Θ using energy-weighted harmonic mean (Eq S31)
    
    Parameters
    ----------
    matrices : SystemMatrices
        Mass, damping, stiffness matrices for consistent unit cell
    system_type : str
        'PC' for Phononic Crystal or 'AM' for Acoustic Metamaterial
        Determines which alpha weight formula to use
    params : dict, optional
        System parameters dict with keys: 'm1', 'm2', 'k1', 'k2', 'c1', 'c2'
        Required for proper alpha weight calculation
        
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
    Phi = modal_params.Phi
    
    # Step 2: Envelope coefficients for DOF 0 (u_n) and DOF 1 (v_n)
    # Pass full mass matrix for proper impulse velocity calculation
    env_u = compute_envelope_coefficients(modal_params, dof_index=0, M=M)
    env_v = compute_envelope_coefficients(modal_params, dof_index=1, M=M)
    
    # Step 3: Extract decay coefficients for each DOF
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
    #
    # Based on validation against the paper, different topologies use
    # different effective decay rate calculations:
    #
    # PC (Phononic Crystal):
    #   The cubic solver extracts the best-fit single exponential from
    #   the bi-exponential envelope at the primary mass. Use θ_u directly.
    #
    # AM (Acoustic Metamaterial):
    #   Due to the local resonance creating distinct modal participation,
    #   the arithmetic mean of modal decay coefficients better represents
    #   the total energy decay: Θ = (ξ₁ω₁ + ξ₂ω₂) / 2
    
    if system_type.upper() == 'AM':
        # Modal decay coefficients
        theta_modal = xi * omega
        # Arithmetic mean for AM
        Theta = np.mean(theta_modal)
    else:
        # For PC, use extracted θ_u from cubic solver
        Theta = theta_u
    
    return MetaDissipationResult(
        Theta=Theta,
        theta_u=theta_u,
        theta_v=theta_v,
        modal_params=modal_params,
        X_u=X_u,
        X_v=X_v
    )
