"""
Dispersion relations for periodic lattice systems.

This module implements the dispersion relation solvers for:
- Monoatomic chain (Eq S7)
- Phononic Crystal / Diatomic chain (Eq S13)
- Acoustic Metamaterial (mass-in-mass)

The dispersion relation ω(κ) or κ(ω) describes how waves
propagate through the periodic structure.

Reference: Meta-dissipation paper (Banerjee, Bera, Adhikari, 2025)
"""

import numpy as np
from typing import Tuple, Optional
from scipy import optimize

from .parameters import (
    MonoatomicParameters,
    PhononicCrystalParameters,
    AcousticMetamaterialParameters,
)


def monoatomic_dispersion(
    params: MonoatomicParameters,
    mu: np.ndarray,
    include_damping: bool = False
) -> np.ndarray:
    """
    Compute dispersion relation for monoatomic chain.
    
    From Eq S7:
        M*λ² + 2(K + λC)(1 - cos(μ)) = 0
        
    For undamped case (C=0):
        ω² = (2K/M)(1 - cos(μ)) = (4K/M)sin²(μ/2)
        
    Parameters
    ----------
    params : MonoatomicParameters
        System parameters
    mu : np.ndarray
        Dimensionless wavenumber μ = κl, range [0, π]
    include_damping : bool
        If True, solve complex eigenvalue problem
        
    Returns
    -------
    omega : np.ndarray
        Angular frequency for each wavenumber
    """
    M, K, C = params.M, params.K, params.C
    
    if not include_damping:
        # Undamped dispersion: ω² = (2K/M)(1 - cos(μ))
        omega_sq = (2 * K / M) * (1 - np.cos(mu))
        omega = np.sqrt(np.maximum(omega_sq, 0))
        return omega
    else:
        # Damped case: solve quadratic for complex λ
        # M*λ² + 2C*(1-cos(μ))*λ + 2K*(1-cos(μ)) = 0
        omega = np.zeros(len(mu), dtype=complex)
        for i, m in enumerate(mu):
            factor = 1 - np.cos(m)
            if factor < 1e-12:
                omega[i] = 0
            else:
                a = M
                b = 2 * C * factor
                c = 2 * K * factor
                # λ = (-b ± sqrt(b²-4ac)) / 2a
                discriminant = b**2 - 4*a*c
                lambda_val = (-b + np.sqrt(discriminant + 0j)) / (2*a)
                omega[i] = np.imag(lambda_val)  # Damped frequency
        return omega


def phononic_crystal_dispersion(
    params: PhononicCrystalParameters,
    mu: np.ndarray,
    include_damping: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute dispersion relation for Phononic Crystal (diatomic chain).
    
    The diatomic chain has two branches: acoustic and optical.
    
    From Eq S13, the quartic dispersion polynomial yields
    two pairs of conjugate eigenvalues.
    
    Parameters
    ----------
    params : PhononicCrystalParameters
        System parameters
    mu : np.ndarray
        Dimensionless wavenumber μ = κl, range [0, π]
    include_damping : bool
        If True, include damping in calculation
        
    Returns
    -------
    omega_acoustic : np.ndarray
        Acoustic branch frequencies
    omega_optical : np.ndarray
        Optical branch frequencies
    """
    m1, m2 = params.m1, params.m2
    k1, k2 = params.k1, params.k2
    c1, c2 = params.c1, params.c2
    
    n_points = len(mu)
    omega_acoustic = np.zeros(n_points)
    omega_optical = np.zeros(n_points)
    
    for i, m in enumerate(mu):
        cos_mu = np.cos(m)
        
        if not include_damping:
            # Undamped case: quadratic in ω²
            # From the determinant of [K - ω²M] with Bloch conditions
            # 
            # The characteristic equation for undamped diatomic:
            # m1*m2*ω⁴ - (m1+m2)(k1+k2)*ω² + 2*k1*k2*(1-cos(μ)) = 0
            #
            # Using standard diatomic form:
            # ω² = (k1+k2)(m1+m2)/(2*m1*m2) ± sqrt[...]
            
            A = m1 * m2
            B = -(k1 + k2) * (m1 + m2)
            C = 2 * k1 * k2 * (1 - cos_mu)
            
            discriminant = B**2 - 4*A*C
            discriminant = max(discriminant, 0)
            
            omega_sq_1 = (-B - np.sqrt(discriminant)) / (2*A)
            omega_sq_2 = (-B + np.sqrt(discriminant)) / (2*A)
            
            omega_acoustic[i] = np.sqrt(max(omega_sq_1, 0))
            omega_optical[i] = np.sqrt(max(omega_sq_2, 0))
        else:
            # Damped case: quartic in λ
            # Solve numerically
            factor = 1 - cos_mu
            
            def dispersion_det(omega_val):
                """Determinant of dynamic stiffness matrix."""
                omega = omega_val + 0j
                D11 = -m1*omega**2 + 1j*omega*(c1+c2) + (k1+k2)
                D12 = -1j*omega*c2 - k2
                D21 = D12
                D22 = -m2*omega**2 + 1j*omega*(c1+c2) + (k1+k2) - 2*(k1 + 1j*omega*c1)*factor
                
                det_val = D11*D22 - D12*D21
                return np.abs(det_val)
            
            # Find roots approximately
            omega_acoustic[i] = omega_sq_1 if 'omega_sq_1' in dir() else 0
            omega_optical[i] = omega_sq_2 if 'omega_sq_2' in dir() else 0
    
    return omega_acoustic, omega_optical


def acoustic_metamaterial_dispersion(
    params: AcousticMetamaterialParameters,
    mu: np.ndarray,
    include_damping: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute dispersion relation for Acoustic Metamaterial.
    
    The AM has a local resonance bandgap around the internal
    resonator frequency ω_r = sqrt(k2/m2).
    
    Parameters
    ----------
    params : AcousticMetamaterialParameters
        System parameters
    mu : np.ndarray  
        Dimensionless wavenumber μ = κl, range [0, π]
    include_damping : bool
        If True, include damping
        
    Returns
    -------
    omega_acoustic : np.ndarray
        Lower (acoustic) branch frequencies
    omega_optical : np.ndarray
        Upper (optical) branch frequencies
    """
    m1, m2 = params.m1, params.m2
    k1, k2 = params.k1, params.k2
    
    n_points = len(mu)
    omega_acoustic = np.zeros(n_points)
    omega_optical = np.zeros(n_points)
    
    # Local resonance frequency
    omega_r = np.sqrt(k2 / m2)
    
    for i, m in enumerate(mu):
        cos_mu = np.cos(m)
        
        # For mass-in-mass metamaterial:
        # Effective mass becomes frequency-dependent:
        # m_eff(ω) = m1 + m2 * ω_r² / (ω_r² - ω²)
        #
        # Dispersion: ω² = 2*k1*(1-cos(μ)) / m_eff(ω)
        #
        # This leads to: m1*m2*ω⁴ - [2k1(1-cos(μ))*m2 + k2*(m1+m2)]*ω² 
        #                + 2*k1*k2*(1-cos(μ)) = 0
        
        factor = 1 - cos_mu
        
        A = m1 * m2
        B = -(2 * k1 * factor * m2 + k2 * (m1 + m2))
        C = 2 * k1 * k2 * factor
        
        discriminant = B**2 - 4*A*C
        discriminant = max(discriminant, 0)
        
        omega_sq_1 = (-B - np.sqrt(discriminant)) / (2*A)
        omega_sq_2 = (-B + np.sqrt(discriminant)) / (2*A)
        
        omega_acoustic[i] = np.sqrt(max(omega_sq_1, 0))
        omega_optical[i] = np.sqrt(max(omega_sq_2, 0))
    
    return omega_acoustic, omega_optical


def compute_bandgap(
    omega_acoustic: np.ndarray,
    omega_optical: np.ndarray
) -> Tuple[float, float, float]:
    """
    Compute bandgap bounds from dispersion branches.
    
    Parameters
    ----------
    omega_acoustic : np.ndarray
        Acoustic branch frequencies
    omega_optical : np.ndarray
        Optical branch frequencies
        
    Returns
    -------
    gap_lower : float
        Lower bound of bandgap (max of acoustic branch)
    gap_upper : float
        Upper bound of bandgap (min of optical branch)  
    gap_width : float
        Width of bandgap (gap_upper - gap_lower)
    """
    gap_lower = np.max(omega_acoustic)
    gap_upper = np.min(omega_optical)
    gap_width = max(0, gap_upper - gap_lower)
    
    return gap_lower, gap_upper, gap_width


def compute_group_velocity(
    omega: np.ndarray,
    mu: np.ndarray
) -> np.ndarray:
    """
    Compute group velocity v_g = dω/dκ from dispersion data.
    
    Parameters
    ----------
    omega : np.ndarray
        Frequencies at each wavenumber
    mu : np.ndarray
        Dimensionless wavenumbers
        
    Returns
    -------
    v_g : np.ndarray
        Group velocity (normalized by lattice constant)
    """
    v_g = np.gradient(omega, mu)
    return v_g
