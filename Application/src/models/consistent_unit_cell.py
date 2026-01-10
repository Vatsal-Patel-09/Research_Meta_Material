"""
Consistent Unit Cell matrices for Meta-dissipation framework.

This module builds the Mass, Damping, and Stiffness matrices
for the Consistent Unit Cells derived via Brillouin Zone integration.

Key Equations:
- Phononic Crystal: Eq S15
- Acoustic Metamaterial: Eq S16

Reference: Meta-dissipation paper (Banerjee, Bera, Adhikari, 2025)
"""

import numpy as np
from typing import Tuple, NamedTuple

from .parameters import (
    MonoatomicParameters,
    PhononicCrystalParameters,
    AcousticMetamaterialParameters,
)


class SystemMatrices(NamedTuple):
    """Container for M, C, K matrices."""
    M: np.ndarray  # Mass matrix
    C: np.ndarray  # Damping matrix
    K: np.ndarray  # Stiffness matrix


def build_monoatomic_matrices(params: MonoatomicParameters) -> SystemMatrices:
    """
    Build consistent unit cell matrices for monoatomic chain.
    
    From Eq S10:
    (M/2)ü + Cú + Ku = 0
    
    This is a 1-DOF system with effective mass M/2.
    
    Parameters
    ----------
    params : MonoatomicParameters
        System parameters
        
    Returns
    -------
    SystemMatrices
        1x1 matrices for the SDOF system
    """
    M = np.array([[params.M / 2]])
    C = np.array([[params.C]])
    K = np.array([[params.K]])
    
    return SystemMatrices(M=M, C=C, K=K)


def build_phononic_crystal_matrices(
    params: PhononicCrystalParameters
) -> SystemMatrices:
    """
    Build consistent unit cell matrices for Phononic Crystal.
    
    From Eq S15 (Table 1):
    
    Mass Matrix:
        M = diag(m1, m2/2)
        
    Damping Matrix:
        C = [[c1+c2, -c2],
             [-c2,    c2]]
             
    Stiffness Matrix:
        K = [[k1+k2, -k2],
             [-k2,    k2]]
    
    Physical interpretation:
    - m1 retains full mass (connected to ground via k1, c1)
    - m2 is halved (boundary sharing in infinite chain)
    - m2 connected to m1 via k2, c2
    
    Parameters
    ----------
    params : PhononicCrystalParameters
        System parameters
        
    Returns
    -------
    SystemMatrices
        2x2 matrices for the 2-DOF system
    """
    # Mass matrix: diag(m1, m2/2)
    M = np.array([
        [params.m1, 0.0],
        [0.0, params.m2 / 2]
    ])
    
    # Damping matrix
    C = np.array([
        [params.c1 + params.c2, -params.c2],
        [-params.c2, params.c2]
    ])
    
    # Stiffness matrix
    K = np.array([
        [params.k1 + params.k2, -params.k2],
        [-params.k2, params.k2]
    ])
    
    return SystemMatrices(M=M, C=C, K=K)


def build_acoustic_metamaterial_matrices(
    params: AcousticMetamaterialParameters
) -> SystemMatrices:
    """
    Build consistent unit cell matrices for Acoustic Metamaterial.
    
    From Eq S16 (Table 1):
    
    Mass Matrix:
        M = diag(m1/2, m2/2)
        
    Damping Matrix:
        C = [[c1 + c2/2, -c2/2],
             [-c2/2,      c2/2]]
             
    Stiffness Matrix:
        K = [[k1 + k2/2, -k2/2],
             [-k2/2,      k2/2]]
    
    Physical interpretation:
    - Both masses halved (symmetric boundary)
    - Internal resonator (k2, c2) contributions halved
    - This topology enables enhanced dissipation
    
    Parameters
    ----------
    params : AcousticMetamaterialParameters
        System parameters
        
    Returns
    -------
    SystemMatrices
        2x2 matrices for the 2-DOF system
    """
    # Mass matrix: diag(m1/2, m2/2)
    M = np.array([
        [params.m1 / 2, 0.0],
        [0.0, params.m2 / 2]
    ])
    
    # Damping matrix (note: c2 terms are halved)
    C = np.array([
        [params.c1 + params.c2 / 2, -params.c2 / 2],
        [-params.c2 / 2, params.c2 / 2]
    ])
    
    # Stiffness matrix (note: k2 terms are halved)
    K = np.array([
        [params.k1 + params.k2 / 2, -params.k2 / 2],
        [-params.k2 / 2, params.k2 / 2]
    ])
    
    return SystemMatrices(M=M, C=C, K=K)


def compute_undamped_frequencies(matrices: SystemMatrices) -> np.ndarray:
    """
    Compute undamped natural frequencies from M and K matrices.
    
    Solves the eigenvalue problem:
        det(K - ω²M) = 0
        
    Parameters
    ----------
    matrices : SystemMatrices
        System matrices
        
    Returns
    -------
    np.ndarray
        Natural frequencies in ascending order [rad/s]
    """
    # Solve generalized eigenvalue problem: K*phi = omega^2 * M * phi
    eigenvalues, _ = np.linalg.eig(np.linalg.solve(matrices.M, matrices.K))
    
    # Take real part and square root for frequencies
    omega_squared = np.real(eigenvalues)
    omega_squared = np.maximum(omega_squared, 0)  # Ensure non-negative
    
    frequencies = np.sqrt(omega_squared)
    return np.sort(frequencies)


def compute_modal_matrices(
    matrices: SystemMatrices
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute mass-normalized modal matrix and modal parameters.
    
    Parameters
    ----------
    matrices : SystemMatrices
        System matrices
        
    Returns
    -------
    Phi : np.ndarray
        Mass-normalized modal matrix (eigenvectors as columns)
    omega : np.ndarray
        Undamped natural frequencies
    modal_M : np.ndarray
        Modal mass matrix (should be identity if normalized)
    """
    M, C, K = matrices
    
    # Solve generalized eigenvalue problem
    eigenvalues, eigenvectors = np.linalg.eig(np.linalg.solve(M, K))
    
    # Sort by eigenvalue magnitude
    idx = np.argsort(np.real(eigenvalues))
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    # Compute natural frequencies
    omega = np.sqrt(np.maximum(np.real(eigenvalues), 0))
    
    # Mass-normalize the eigenvectors
    Phi = np.zeros_like(eigenvectors)
    for j in range(eigenvectors.shape[1]):
        phi_j = eigenvectors[:, j]
        modal_mass = phi_j.T @ M @ phi_j
        Phi[:, j] = phi_j / np.sqrt(np.abs(modal_mass))
    
    # Verify: modal_M should be identity
    modal_M = Phi.T @ M @ Phi
    
    return Phi, omega, modal_M


def compute_modal_damping_ratios(
    matrices: SystemMatrices,
    Phi: np.ndarray,
    omega: np.ndarray
) -> np.ndarray:
    """
    Compute modal damping ratios assuming proportional damping.
    
    For proportional (Rayleigh) damping:
        ξ_j = (Φ_j^T C Φ_j) / (2 ω_j)
        
    Parameters
    ----------
    matrices : SystemMatrices
        System matrices
    Phi : np.ndarray
        Mass-normalized modal matrix
    omega : np.ndarray
        Natural frequencies
        
    Returns
    -------
    np.ndarray
        Modal damping ratios
    """
    modal_C = Phi.T @ matrices.C @ Phi
    
    # For proportional damping, modal_C is diagonal
    # ξ_j = modal_C[j,j] / (2 * omega_j)
    xi = np.zeros(len(omega))
    for j in range(len(omega)):
        if omega[j] > 1e-10:
            xi[j] = modal_C[j, j] / (2 * omega[j])
        else:
            xi[j] = 0.0
            
    return xi
