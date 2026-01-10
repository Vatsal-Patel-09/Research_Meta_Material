"""
Parameter data classes for Meta-dissipation framework.

This module defines the system parameters for:
- Monoatomic chains
- Phononic Crystals (PC) - diatomic chains  
- Acoustic Metamaterials (AM) - mass-in-mass topology

Reference: Meta-dissipation paper (Banerjee, Bera, Adhikari, 2025)
"""

from dataclasses import dataclass, field
from typing import Optional
import numpy as np


@dataclass
class MonoatomicParameters:
    """
    Parameters for a monoatomic periodic chain.
    
    Attributes
    ----------
    M : float
        Total mass of the unit cell [kg]
    K : float  
        Spring stiffness [N/m]
    C : float
        Damping coefficient [Ns/m]
    l : float
        Lattice constant (unit cell length) [m]
    """
    M: float
    K: float
    C: float
    l: float = 1.0
    
    @property
    def omega_n(self) -> float:
        """Undamped natural frequency of consistent unit cell."""
        return np.sqrt(2 * self.K / self.M)
    
    @property
    def zeta(self) -> float:
        """Damping ratio of consistent unit cell."""
        return self.C / (2 * np.sqrt(self.K * self.M / 2))


@dataclass
class DiatomicParameters:
    """
    Base parameters for two-mass periodic systems.
    
    Attributes
    ----------
    m1 : float
        Primary mass [kg]
    m2 : float
        Secondary mass [kg]  
    k1 : float
        Primary stiffness [N/m]
    k2 : float
        Secondary/coupling stiffness [N/m]
    c1 : float
        Primary damping coefficient [Ns/m]
    c2 : float
        Secondary/coupling damping coefficient [Ns/m]
    l : float
        Lattice constant [m]
    """
    m1: float
    m2: float
    k1: float
    k2: float
    c1: float
    c2: float
    l: float = 1.0
    
    @property
    def total_mass(self) -> float:
        """Total mass of the unit cell."""
        return self.m1 + self.m2
    
    @property
    def mass_ratio(self) -> float:
        """Ratio of secondary to primary mass."""
        return self.m2 / self.m1


@dataclass
class PhononicCrystalParameters(DiatomicParameters):
    """
    Parameters for Phononic Crystal (diatomic chain).
    
    The PC topology consists of two distinct masses connected
    in series by alternating springs k1 and k2.
    
    Static sound speed (Eq S33):
    C_stat = l * sqrt(k1*k2 / ((m1+m2)*(k1+k2)))
    """
    
    @property
    def static_sound_speed(self) -> float:
        """
        Long-wavelength (static) sound speed.
        Eq S33 for PC topology.
        """
        numerator = self.k1 * self.k2
        denominator = self.total_mass * (self.k1 + self.k2)
        return self.l * np.sqrt(numerator / denominator)
    
    def __post_init__(self):
        """Validate parameters after initialization."""
        if any(x <= 0 for x in [self.m1, self.m2, self.k1, self.k2]):
            raise ValueError("Mass and stiffness values must be positive")


@dataclass  
class AcousticMetamaterialParameters(DiatomicParameters):
    """
    Parameters for Acoustic Metamaterial (mass-in-mass topology).
    
    The AM topology has an internal resonator (m2, k2, c2) 
    suspended inside the primary mass (m1), which connects
    to adjacent cells via k1, c1.
    
    Static sound speed (Eq S33):
    C_stat = l * sqrt(k1 / (m1 + m2))
    """
    
    @property
    def static_sound_speed(self) -> float:
        """
        Long-wavelength (static) sound speed.
        Eq S33 for AM topology.
        """
        return self.l * np.sqrt(self.k1 / self.total_mass)
    
    @property
    def local_resonance_frequency(self) -> float:
        """Natural frequency of the internal resonator."""
        return np.sqrt(self.k2 / self.m2)
    
    def __post_init__(self):
        """Validate parameters after initialization."""
        if any(x <= 0 for x in [self.m1, self.m2, self.k1, self.k2]):
            raise ValueError("Mass and stiffness values must be positive")


# ============================================================================
# Benchmark Parameters from Paper (Section 7.1)
# ============================================================================

def get_benchmark_pc_parameters() -> PhononicCrystalParameters:
    """
    Get benchmark Phononic Crystal parameters from the paper.
    
    Returns parameters that yield C_stat ≈ 83.33
    """
    return PhononicCrystalParameters(
        m1=1.0,
        m2=0.8,
        k1=40906.0,
        k2=18000.0,
        c1=20.0,
        c2=8.8,
        l=1.0
    )


def get_benchmark_am_parameters() -> AcousticMetamaterialParameters:
    """
    Get benchmark Acoustic Metamaterial parameters from the paper.
    
    Returns parameters that yield C_stat ≈ 83.33
    """
    return AcousticMetamaterialParameters(
        m1=1.0,
        m2=0.8,
        k1=12500.0,
        k2=5500.0,
        c1=20.0,
        c2=8.8,
        l=1.0
    )


def validate_equal_sound_speed(
    pc: PhononicCrystalParameters,
    am: AcousticMetamaterialParameters,
    tol: float = 0.1
) -> bool:
    """
    Validate that PC and AM have equal static sound speeds.
    
    This is required for fair comparison of dissipation performance.
    
    Parameters
    ----------
    pc : PhononicCrystalParameters
        Phononic crystal parameters
    am : AcousticMetamaterialParameters  
        Acoustic metamaterial parameters
    tol : float
        Tolerance for equality check
        
    Returns
    -------
    bool
        True if sound speeds match within tolerance
    """
    diff = abs(pc.static_sound_speed - am.static_sound_speed)
    return diff < tol
