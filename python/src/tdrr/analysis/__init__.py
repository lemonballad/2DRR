"""Spectroscopy analysis modules.

Contains functions for:
- Mode mapping
- Trajectory generation
- Correlation analysis
- Force calculations
"""

from tdrr.analysis.mode_mapping import build_mode_map, interpolate_mode_map
from tdrr.analysis.trajectory import (
    generate_trajectory,
    compute_correlation,
    compute_spectral_density,
    compute_2d_correlation,
    fit_exponential_decay,
    compute_vibrational_coupling,
    compute_2d_spectral_density,
    process_multiple_trajectories,
)
from tdrr.analysis.force import (
    compute_coulomb_energy,
    compute_lennard_jones_energy,
    compute_total_interaction_energy,
    compute_energy_along_normal,
    compute_plane_normal,
    get_residue_range,
    compute_force_vectors,
)

__all__ = [
    # Mode mapping
    "build_mode_map",
    "interpolate_mode_map",
    # Trajectory
    "generate_trajectory",
    "compute_correlation",
    "compute_spectral_density",
    "compute_2d_correlation",
    "fit_exponential_decay",
    "compute_vibrational_coupling",
    "compute_2d_spectral_density",
    "process_multiple_trajectories",
    # Force
    "compute_coulomb_energy",
    "compute_lennard_jones_energy",
    "compute_total_interaction_energy",
    "compute_energy_along_normal",
    "compute_plane_normal",
    "get_residue_range",
    "compute_force_vectors",
]
