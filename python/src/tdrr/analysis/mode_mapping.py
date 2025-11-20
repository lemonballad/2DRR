"""Mode mapping functions for vibrational analysis.

Functions for building 2D maps of vibrational modes across dihedral angle space.
"""

import numpy as np
from numpy.typing import NDArray


def build_mode_map(
    frequencies: NDArray[np.float64],
    l_angles: NDArray[np.float64],
    r_angles: NDArray[np.float64],
    mode_index: int,
    angle_resolution: int = 361,
) -> NDArray[np.float64]:
    """Build 2D frequency map across dihedral angle space.

    Args:
        frequencies: Array of frequencies for each angle combination.
            Shape: (n_l_angles, n_r_angles) or (n_points,)
        l_angles: Left dihedral angles in degrees.
        r_angles: Right dihedral angles in degrees.
        mode_index: Index of the vibrational mode to map.
        angle_resolution: Resolution of output map (default 361 for 1 degree).

    Returns:
        2D array of frequencies indexed by [l_angle + 180, r_angle + 180].
    """
    freq_map = np.zeros((angle_resolution, angle_resolution), dtype=np.float64)

    # Normalize angles to [0, 360] range for indexing
    l_idx = (np.round(l_angles).astype(int) + 180) % angle_resolution
    r_idx = (np.round(r_angles).astype(int) + 180) % angle_resolution

    if frequencies.ndim == 1:
        for i, (li, ri) in enumerate(zip(l_idx, r_idx)):
            freq_map[li, ri] = frequencies[i]
    else:
        freq_map = frequencies

    return freq_map


def interpolate_mode_map(
    freq_map: NDArray[np.float64],
    l_angles: NDArray[np.float64],
    r_angles: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Interpolate frequencies from mode map for given angles.

    Args:
        freq_map: 2D frequency map from build_mode_map.
        l_angles: Query left dihedral angles in degrees.
        r_angles: Query right dihedral angles in degrees.

    Returns:
        Array of interpolated frequencies.
    """
    from scipy.interpolate import RectBivariateSpline

    # Create coordinate grids
    n = freq_map.shape[0]
    angles = np.arange(n) - 180

    # Create interpolator
    interp = RectBivariateSpline(angles, angles, freq_map)

    # Normalize query angles
    l_norm = np.mod(l_angles + 180, 360) - 180
    r_norm = np.mod(r_angles + 180, 360) - 180

    return interp(l_norm, r_norm, grid=False)


def compute_normal_mode_overlap(
    mode_coords: NDArray[np.float64],
    reference_coords: NDArray[np.float64],
    atom_indices: NDArray[np.int64] | None = None,
) -> float:
    """Compute overlap between normal mode vectors.

    Args:
        mode_coords: Normal mode displacement coordinates (N_atoms x 3).
        reference_coords: Reference mode coordinates (N_atoms x 3).
        atom_indices: Optional subset of atom indices to use.

    Returns:
        Overlap value (0 to 1).
    """
    if atom_indices is not None:
        mode_coords = mode_coords[atom_indices]
        reference_coords = reference_coords[atom_indices]

    # Normalize vectors
    mode_norm = np.sqrt(np.trace(mode_coords.T @ mode_coords))
    ref_norm = np.sqrt(np.trace(reference_coords.T @ reference_coords))

    mode_coords = mode_coords / mode_norm
    reference_coords = reference_coords / ref_norm

    # Compute overlap
    overlap = np.abs(np.sqrt(np.trace(mode_coords.T @ reference_coords)))

    return float(overlap)


def find_mode_correspondence(
    coordinates: NDArray[np.float64],
    reference_mode: int,
    l_angle: int,
    r_angle: int,
    atom_indices: NDArray[np.int64],
    threshold: float = 0.7,
) -> int:
    """Find corresponding mode at different dihedral angles.

    Args:
        coordinates: 5D array of normal mode coordinates
            (n_modes, n_atoms, 3, n_l_angles, n_r_angles).
        reference_mode: Index of reference mode at (0, 0) angles.
        l_angle: Target left dihedral angle index.
        r_angle: Target right dihedral angle index.
        atom_indices: Atom indices to use for comparison.
        threshold: Minimum overlap threshold for correspondence.

    Returns:
        Index of corresponding mode at target angles, or -1 if not found.
    """
    n_modes = coordinates.shape[0]
    ref_center = coordinates.shape[3] // 2

    # Get reference mode coordinates
    ref_coords = coordinates[reference_mode, :, :, ref_center, ref_center]
    ref_subset = ref_coords[atom_indices]
    ref_norm = np.sqrt(np.trace(ref_subset.T @ ref_subset))
    ref_subset = ref_subset / ref_norm

    best_overlap = 0.0
    best_mode = -1

    for mode_idx in range(n_modes):
        mode_coords = coordinates[mode_idx, :, :, l_angle, r_angle]
        mode_subset = mode_coords[atom_indices]
        mode_norm = np.sqrt(np.trace(mode_subset.T @ mode_subset))
        mode_subset = mode_subset / mode_norm

        overlap = np.abs(np.sqrt(np.trace(mode_subset.T @ ref_subset)))

        if overlap > best_overlap:
            best_overlap = overlap
            best_mode = mode_idx

    if best_overlap >= threshold:
        return best_mode

    return -1
