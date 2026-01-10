"""
Impulse response and time-domain simulation.

This module computes the transient response of the consistent
unit cell to an impulse excitation.

Key Equations:
- Modal superposition: Eq S20
- Bi-exponential envelope: Eq S22

Reference: Meta-dissipation paper (Banerjee, Bera, Adhikari, 2025)
"""

import numpy as np
from typing import Tuple, Optional
from scipy.integrate import odeint

from ..models.consistent_unit_cell import SystemMatrices


def impulse_response_modal(
    matrices: SystemMatrices,
    t: np.ndarray,
    impulse_dof: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute impulse response using modal superposition.
    
    From Eq S20:
        u_p(t) = Σⱼ φ_pj * (q̇ⱼ(0)/ω_dj) * exp(-ξⱼωⱼt) * sin(ω_dj*t)
        
    Parameters
    ----------
    matrices : SystemMatrices
        System matrices
    t : np.ndarray
        Time vector
    impulse_dof : int
        DOF where impulse is applied (default: 0)
        
    Returns
    -------
    u : np.ndarray
        Response of DOF 0 (primary mass)
    v : np.ndarray
        Response of DOF 1 (secondary mass)
    """
    from ..models.consistent_unit_cell import compute_modal_matrices, compute_modal_damping_ratios
    
    M, C, K = matrices
    
    # Modal analysis
    Phi, omega, _ = compute_modal_matrices(matrices)
    xi = compute_modal_damping_ratios(matrices, Phi, omega)
    
    # Damped frequencies
    omega_d = omega * np.sqrt(np.maximum(1 - xi**2, 0.0001))
    
    # Initial conditions from impulse
    # Impulse gives initial velocity: v(0) = [1/M_11, 0]^T for impulse on DOF 0
    v0 = np.zeros(2)
    v0[impulse_dof] = 1.0 / M[impulse_dof, impulse_dof]
    
    # Transform to modal coordinates
    # q̇(0) = Φ^T M v(0) (for mass-normalized Φ)
    q_dot_0 = Phi.T @ M @ v0
    
    # Initialize response arrays
    u = np.zeros(len(t))
    v = np.zeros(len(t))
    
    # Modal superposition (Eq S20)
    for j in range(2):
        if omega_d[j] > 1e-10:
            decay = np.exp(-xi[j] * omega[j] * t)
            oscillation = np.sin(omega_d[j] * t)
            amplitude = q_dot_0[j] / omega_d[j]
            
            u += Phi[0, j] * amplitude * decay * oscillation
            v += Phi[1, j] * amplitude * decay * oscillation
    
    return u, v


def impulse_response_numerical(
    matrices: SystemMatrices,
    t: np.ndarray,
    impulse_dof: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute impulse response via numerical integration.
    
    Solves: M*ẍ + C*ẋ + K*x = 0
    with initial conditions: x(0)=0, ẋ(0)=impulse
    
    Parameters
    ----------
    matrices : SystemMatrices
        System matrices
    t : np.ndarray
        Time vector
    impulse_dof : int
        DOF where impulse is applied
        
    Returns
    -------
    u, v : np.ndarray
        Responses of DOF 0 and DOF 1
    """
    M, C, K = matrices
    M_inv = np.linalg.inv(M)
    
    def equations(y, t):
        """State-space equations: [x, x_dot] -> [x_dot, x_ddot]"""
        x = y[:2]
        x_dot = y[2:]
        x_ddot = -M_inv @ (C @ x_dot + K @ x)
        return np.concatenate([x_dot, x_ddot])
    
    # Initial conditions
    x0 = np.zeros(2)
    v0 = np.zeros(2)
    v0[impulse_dof] = 1.0 / M[impulse_dof, impulse_dof]
    y0 = np.concatenate([x0, v0])
    
    # Integrate
    solution = odeint(equations, y0, t)
    
    u = solution[:, 0]
    v = solution[:, 1]
    
    return u, v


def compute_envelope(
    u: np.ndarray,
    t: np.ndarray,
    method: str = 'peaks'
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract envelope from oscillatory response.
    
    Parameters
    ----------
    u : np.ndarray
        Time-domain response
    t : np.ndarray
        Time vector
    method : str
        'peaks' - interpolate through local maxima
        'hilbert' - use Hilbert transform
        
    Returns
    -------
    t_env : np.ndarray
        Time points for envelope
    envelope : np.ndarray
        Envelope values
    """
    from scipy.signal import hilbert, find_peaks
    from scipy.interpolate import interp1d
    
    if method == 'hilbert':
        analytic = hilbert(u)
        envelope = np.abs(analytic)
        return t, envelope
    
    elif method == 'peaks':
        # Find positive peaks
        peaks_pos, _ = find_peaks(u)
        peaks_neg, _ = find_peaks(-u)
        
        # Combine and sort
        all_peaks = np.sort(np.concatenate([peaks_pos, peaks_neg]))
        
        if len(all_peaks) < 2:
            return t, np.abs(u)
        
        t_peaks = t[all_peaks]
        u_peaks = np.abs(u[all_peaks])
        
        # Interpolate
        f = interp1d(t_peaks, u_peaks, kind='linear', 
                     bounds_error=False, fill_value=(u_peaks[0], u_peaks[-1]))
        envelope = f(t)
        
        return t, envelope
    
    else:
        raise ValueError(f"Unknown method: {method}")


def biexponential_envelope(
    t: np.ndarray,
    A: float, B: float,
    alpha1: float, alpha2: float
) -> np.ndarray:
    """
    Theoretical bi-exponential envelope (Eq S22).
    
    y(t) = A * exp(-α₁t) + B * exp(-α₂t)
    
    where α_j = ξ_j * ω_j
    
    Parameters
    ----------
    t : np.ndarray
        Time vector
    A, B : float
        Amplitude coefficients
    alpha1, alpha2 : float
        Decay rates (ξ₁ω₁ and ξ₂ω₂)
        
    Returns
    -------
    envelope : np.ndarray
        Bi-exponential envelope
    """
    return A * np.exp(-alpha1 * t) + B * np.exp(-alpha2 * t)


def single_exponential_fit(
    t: np.ndarray,
    X: float,
    theta: float
) -> np.ndarray:
    """
    Single exponential approximation to the envelope.
    
    y(t) ≈ X * exp(-θt)
    
    Parameters
    ----------
    t : np.ndarray
        Time vector
    X : float
        Initial amplitude
    theta : float
        Decay coefficient (from Meta-dissipation algorithm)
        
    Returns
    -------
    envelope : np.ndarray
        Single exponential envelope
    """
    return X * np.exp(-theta * t)


def compute_energy(
    matrices: SystemMatrices,
    u: np.ndarray,
    v: np.ndarray,
    u_dot: np.ndarray,
    v_dot: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute kinetic, potential, and total energy.
    
    E_kinetic = 0.5 * ẋ^T M ẋ
    E_potential = 0.5 * x^T K x
    E_total = E_kinetic + E_potential
    
    Parameters
    ----------
    matrices : SystemMatrices
        System matrices
    u, v : np.ndarray
        Displacement responses
    u_dot, v_dot : np.ndarray
        Velocity responses
        
    Returns
    -------
    E_kinetic, E_potential, E_total : np.ndarray
        Energy time histories
    """
    M, C, K = matrices
    n_points = len(u)
    
    E_kinetic = np.zeros(n_points)
    E_potential = np.zeros(n_points)
    
    for i in range(n_points):
        x = np.array([u[i], v[i]])
        x_dot = np.array([u_dot[i], v_dot[i]])
        
        E_kinetic[i] = 0.5 * x_dot.T @ M @ x_dot
        E_potential[i] = 0.5 * x.T @ K @ x
    
    E_total = E_kinetic + E_potential
    
    return E_kinetic, E_potential, E_total


def compute_energy_from_displacement(
    matrices: SystemMatrices,
    t: np.ndarray,
    u: np.ndarray,
    v: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute energy using numerical differentiation for velocities.
    
    Parameters
    ----------
    matrices : SystemMatrices
        System matrices
    t : np.ndarray
        Time vector
    u, v : np.ndarray
        Displacement responses
        
    Returns
    -------
    E_kinetic, E_potential, E_total : np.ndarray
        Energy time histories
    """
    # Numerical differentiation for velocity
    dt = np.gradient(t)
    u_dot = np.gradient(u, t)
    v_dot = np.gradient(v, t)
    
    return compute_energy(matrices, u, v, u_dot, v_dot)
