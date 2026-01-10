"""
Visualization module for Meta-dissipation framework.

This module provides plotting functions for:
- Dispersion curves
- Time-domain responses
- Energy decay
- Comparison charts

Reference: Meta-dissipation paper (Banerjee, Bera, Adhikari, 2025)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from typing import Tuple, Optional, List

# Set publication-quality defaults
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 11,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'figure.figsize': (10, 6),
    'lines.linewidth': 2,
    'axes.grid': True,
    'grid.alpha': 0.3,
})


def plot_dispersion_comparison(
    mu: np.ndarray,
    omega_pc: Tuple[np.ndarray, np.ndarray],
    omega_am: Tuple[np.ndarray, np.ndarray],
    title: str = "Dispersion Relations: PC vs AM"
) -> Tuple[Figure, Axes]:
    """
    Plot dispersion curves comparing PC and AM.
    
    Parameters
    ----------
    mu : np.ndarray
        Dimensionless wavenumber
    omega_pc : tuple
        (acoustic, optical) branches for PC
    omega_am : tuple
        (acoustic, optical) branches for AM
    title : str
        Plot title
        
    Returns
    -------
    fig, ax : Figure, Axes
        Matplotlib figure and axes
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # PC curves
    ax.plot(mu/np.pi, omega_pc[0], 'b-', label='PC Acoustic', linewidth=2)
    ax.plot(mu/np.pi, omega_pc[1], 'b--', label='PC Optical', linewidth=2)
    
    # AM curves  
    ax.plot(mu/np.pi, omega_am[0], 'r-', label='AM Acoustic', linewidth=2)
    ax.plot(mu/np.pi, omega_am[1], 'r--', label='AM Optical', linewidth=2)
    
    ax.set_xlabel(r'Normalized Wavenumber $\mu/\pi$')
    ax.set_ylabel(r'Angular Frequency $\omega$ [rad/s]')
    ax.set_title(title)
    ax.legend()
    ax.set_xlim([0, 1])
    
    plt.tight_layout()
    return fig, ax


def plot_impulse_response(
    t: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    title: str = "Impulse Response",
    label_u: str = r"$u_n$ (Primary mass)",
    label_v: str = r"$v_n$ (Secondary mass)"
) -> Tuple[Figure, Axes]:
    """
    Plot time-domain impulse response.
    
    Parameters
    ----------
    t : np.ndarray
        Time vector
    u, v : np.ndarray
        Displacement responses
    title : str
        Plot title
        
    Returns
    -------
    fig, ax : Figure, Axes
    """
    fig, ax = plt.subplots(figsize=(12, 5))
    
    ax.plot(t, u, 'b-', label=label_u, linewidth=1.5)
    ax.plot(t, v, 'r-', label=label_v, linewidth=1.5, alpha=0.8)
    
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Displacement')
    ax.set_title(title)
    ax.legend()
    ax.axhline(y=0, color='k', linewidth=0.5)
    
    plt.tight_layout()
    return fig, ax


def plot_response_comparison(
    t: np.ndarray,
    u_pc: np.ndarray,
    u_am: np.ndarray,
    title: str = "Response Comparison: PC vs AM"
) -> Tuple[Figure, Axes]:
    """
    Compare impulse responses of PC and AM.
    
    Parameters
    ----------
    t : np.ndarray
        Time vector
    u_pc, u_am : np.ndarray
        Primary mass responses for PC and AM
    title : str
        Plot title
        
    Returns
    -------
    fig, ax : Figure, Axes
    """
    fig, ax = plt.subplots(figsize=(12, 5))
    
    ax.plot(t, u_pc, 'b-', label='Phononic Crystal', linewidth=1.5)
    ax.plot(t, u_am, 'r-', label='Acoustic Metamaterial', linewidth=1.5)
    
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Displacement')
    ax.set_title(title)
    ax.legend()
    ax.axhline(y=0, color='k', linewidth=0.5)
    
    plt.tight_layout()
    return fig, ax


def plot_envelope_analysis(
    t: np.ndarray,
    response: np.ndarray,
    envelope_biexp: np.ndarray,
    envelope_single: np.ndarray,
    theta: float,
    title: str = "Envelope Analysis"
) -> Tuple[Figure, Axes]:
    """
    Plot response with bi-exponential and single exponential envelopes.
    
    Parameters
    ----------
    t : np.ndarray
        Time vector
    response : np.ndarray
        Time-domain response
    envelope_biexp : np.ndarray
        Bi-exponential envelope (Eq S22)
    envelope_single : np.ndarray  
        Single exponential fit (Meta-dissipation)
    theta : float
        Decay coefficient
    title : str
        Plot title
        
    Returns
    -------
    fig, ax : Figure, Axes
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.plot(t, response, 'b-', alpha=0.5, linewidth=1, label='Response')
    ax.plot(t, envelope_biexp, 'g--', linewidth=2, 
            label='Bi-exponential envelope')
    ax.plot(t, envelope_single, 'r-', linewidth=2.5,
            label=f'Single exp. fit (θ={theta:.2f})')
    ax.plot(t, -envelope_biexp, 'g--', linewidth=2)
    ax.plot(t, -envelope_single, 'r-', linewidth=2.5)
    
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Displacement')
    ax.set_title(title)
    ax.legend()
    ax.axhline(y=0, color='k', linewidth=0.5)
    
    plt.tight_layout()
    return fig, ax


def plot_energy_decay(
    t: np.ndarray,
    E_pc: np.ndarray,
    E_am: np.ndarray,
    Theta_pc: float,
    Theta_am: float,
    normalize: bool = True,
    title: str = "Energy Decay Comparison"
) -> Tuple[Figure, Axes]:
    """
    Plot energy decay for PC and AM.
    
    Parameters
    ----------
    t : np.ndarray
        Time vector
    E_pc, E_am : np.ndarray
        Total energy for PC and AM
    Theta_pc, Theta_am : float
        Unified decay coefficients
    normalize : bool
        If True, normalize to initial energy
    title : str
        Plot title
        
    Returns
    -------
    fig, ax : Figure, Axes
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if normalize:
        E_pc = E_pc / E_pc[0] if E_pc[0] > 0 else E_pc
        E_am = E_am / E_am[0] if E_am[0] > 0 else E_am
        ylabel = 'Normalized Energy'
    else:
        ylabel = 'Energy'
    
    ax.semilogy(t, E_pc, 'b-', linewidth=2, label=f'PC (Θ={Theta_pc:.2f})')
    ax.semilogy(t, E_am, 'r-', linewidth=2, label=f'AM (Θ={Theta_am:.2f})')
    
    # Theoretical decay lines
    E0_pc = E_pc[0] if E_pc[0] > 0 else 1
    E0_am = E_am[0] if E_am[0] > 0 else 1
    ax.semilogy(t, E0_pc * np.exp(-2*Theta_pc*t), 'b--', linewidth=1.5, 
                alpha=0.7, label=f'PC theory: exp(-2Θt)')
    ax.semilogy(t, E0_am * np.exp(-2*Theta_am*t), 'r--', linewidth=1.5,
                alpha=0.7, label=f'AM theory: exp(-2Θt)')
    
    ax.set_xlabel('Time [s]')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend()
    
    plt.tight_layout()
    return fig, ax


def plot_theta_comparison(
    Theta_pc: float,
    Theta_am: float,
    title: str = "Meta-dissipation Coefficient Comparison"
) -> Tuple[Figure, Axes]:
    """
    Bar chart comparing Θ values.
    
    Parameters
    ----------
    Theta_pc : float
        Unified decay coefficient for PC
    Theta_am : float
        Unified decay coefficient for AM
    title : str
        Plot title
        
    Returns
    -------
    fig, ax : Figure, Axes
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    systems = ['Phononic Crystal\n(PC)', 'Acoustic Metamaterial\n(AM)']
    values = [Theta_pc, Theta_am]
    colors = ['#3498db', '#e74c3c']
    
    bars = ax.bar(systems, values, color=colors, edgecolor='black', linewidth=1.5)
    
    # Add value labels on bars
    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax.annotate(f'Θ = {val:.2f}',
                   xy=(bar.get_x() + bar.get_width()/2, height),
                   xytext=(0, 3),
                   textcoords="offset points",
                   ha='center', va='bottom',
                   fontsize=14, fontweight='bold')
    
    ax.set_ylabel(r'Unified Decay Coefficient $\Theta$ [1/s]')
    ax.set_title(title)
    
    # Add improvement annotation
    improvement = (Theta_am / Theta_pc - 1) * 100
    ax.annotate(f'AM shows {improvement:.0f}% higher dissipation',
               xy=(0.5, 0.95), xycoords='axes fraction',
               ha='center', fontsize=12,
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    return fig, ax


def plot_modal_parameters(
    omega_pc: np.ndarray, xi_pc: np.ndarray,
    omega_am: np.ndarray, xi_am: np.ndarray,
    title: str = "Modal Parameters"
) -> Tuple[Figure, List[Axes]]:
    """
    Plot modal frequencies and damping ratios.
    
    Parameters
    ----------
    omega_pc, xi_pc : np.ndarray
        PC natural frequencies and damping ratios
    omega_am, xi_am : np.ndarray
        AM natural frequencies and damping ratios
    title : str
        Plot title
        
    Returns
    -------
    fig, axes : Figure, List[Axes]
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Frequencies
    ax1 = axes[0]
    x = np.arange(2)
    width = 0.35
    
    ax1.bar(x - width/2, omega_pc, width, label='PC', color='#3498db')
    ax1.bar(x + width/2, omega_am, width, label='AM', color='#e74c3c')
    
    ax1.set_ylabel(r'Natural Frequency $\omega$ [rad/s]')
    ax1.set_xlabel('Mode')
    ax1.set_xticks(x)
    ax1.set_xticklabels(['Mode 1', 'Mode 2'])
    ax1.legend()
    ax1.set_title('Natural Frequencies')
    
    # Damping ratios
    ax2 = axes[1]
    ax2.bar(x - width/2, xi_pc * 100, width, label='PC', color='#3498db')
    ax2.bar(x + width/2, xi_am * 100, width, label='AM', color='#e74c3c')
    
    ax2.set_ylabel(r'Damping Ratio $\xi$ [%]')
    ax2.set_xlabel('Mode')
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Mode 1', 'Mode 2'])
    ax2.legend()
    ax2.set_title('Modal Damping Ratios')
    
    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    return fig, axes


def create_summary_figure(
    t: np.ndarray,
    u_pc: np.ndarray, u_am: np.ndarray,
    E_pc: np.ndarray, E_am: np.ndarray,
    Theta_pc: float, Theta_am: float,
    title: str = "Meta-dissipation Analysis Summary"
) -> Tuple[Figure, np.ndarray]:
    """
    Create comprehensive summary figure with multiple subplots.
    
    Parameters
    ----------
    t : np.ndarray
        Time vector
    u_pc, u_am : np.ndarray
        Displacement responses
    E_pc, E_am : np.ndarray
        Energy histories
    Theta_pc, Theta_am : float
        Unified decay coefficients
    title : str
        Overall title
        
    Returns
    -------
    fig, axes : Figure, ndarray of Axes
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # (a) Time response comparison
    ax1 = axes[0, 0]
    ax1.plot(t, u_pc, 'b-', linewidth=1.5, label='PC')
    ax1.plot(t, u_am, 'r-', linewidth=1.5, label='AM')
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('Displacement')
    ax1.set_title('(a) Impulse Response Comparison')
    ax1.legend()
    ax1.axhline(y=0, color='k', linewidth=0.5)
    
    # (b) Energy decay (log scale)
    ax2 = axes[0, 1]
    E_pc_norm = E_pc / E_pc[0] if E_pc[0] > 0 else E_pc
    E_am_norm = E_am / E_am[0] if E_am[0] > 0 else E_am
    ax2.semilogy(t, E_pc_norm, 'b-', linewidth=2, label=f'PC')
    ax2.semilogy(t, E_am_norm, 'r-', linewidth=2, label=f'AM')
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Normalized Energy')
    ax2.set_title('(b) Energy Decay')
    ax2.legend()
    
    # (c) Zoomed initial response
    ax3 = axes[1, 0]
    t_zoom = t[t < 0.1] if len(t[t < 0.1]) > 10 else t[:100]
    n_zoom = len(t_zoom)
    ax3.plot(t_zoom, u_pc[:n_zoom], 'b-', linewidth=1.5, label='PC')
    ax3.plot(t_zoom, u_am[:n_zoom], 'r-', linewidth=1.5, label='AM')
    ax3.set_xlabel('Time [s]')
    ax3.set_ylabel('Displacement')
    ax3.set_title('(c) Initial Response (Zoomed)')
    ax3.legend()
    ax3.axhline(y=0, color='k', linewidth=0.5)
    
    # (d) Theta comparison bar chart
    ax4 = axes[1, 1]
    systems = ['PC', 'AM']
    values = [Theta_pc, Theta_am]
    colors = ['#3498db', '#e74c3c']
    bars = ax4.bar(systems, values, color=colors, edgecolor='black')
    for bar, val in zip(bars, values):
        ax4.annotate(f'Θ={val:.2f}',
                    xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', fontsize=12, fontweight='bold')
    ax4.set_ylabel(r'$\Theta$ [1/s]')
    ax4.set_title('(d) Meta-dissipation Coefficient')
    
    fig.suptitle(title, fontsize=16, fontweight='bold')
    plt.tight_layout()
    return fig, axes
