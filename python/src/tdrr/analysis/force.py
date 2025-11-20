"""Force and energy analysis functions.

Functions for computing electrostatic and van der Waals interactions
between molecular components.
"""

import numpy as np
from numpy.typing import NDArray


# Physical constants
EPSILON_0 = 8.85e-12  # Vacuum permittivity (F/m)
ANG_PER_M = 1e10  # Angstroms per meter
C_PER_E = 1.6e-19  # Coulombs per electron charge
WVN_PER_J = 5.034e22  # Wavenumbers per Joule
WVN_MOL_PER_KJ = 83.593  # Wavenumbers per kJ/mol


def compute_coulomb_energy(
    coords1: NDArray[np.float64],
    charges1: NDArray[np.float64],
    coords2: NDArray[np.float64],
    charges2: NDArray[np.float64],
) -> float:
    """Compute Coulomb interaction energy between two sets of atoms.

    Args:
        coords1: Coordinates of first atom set (N1 x 3) in Angstroms.
        charges1: Charges of first atom set (N1,) in electron charges.
        coords2: Coordinates of second atom set (N2 x 3) in Angstroms.
        charges2: Charges of second atom set (N2,) in electron charges.

    Returns:
        Coulomb energy in wavenumbers (cm^-1).
    """
    prefactor = ANG_PER_M * C_PER_E ** 2 * WVN_PER_J / (4 * np.pi * EPSILON_0)

    # Compute pairwise distances
    diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff ** 2, axis=2))

    # Compute charge products
    charge_products = charges1[:, np.newaxis] * charges2[np.newaxis, :]

    # Sum Coulomb interactions
    energy = prefactor * np.sum(charge_products / distances)

    return float(energy)


def compute_lennard_jones_energy(
    coords1: NDArray[np.float64],
    c6_1: NDArray[np.float64],
    c12_1: NDArray[np.float64],
    coords2: NDArray[np.float64],
    c6_2: NDArray[np.float64],
    c12_2: NDArray[np.float64],
) -> float:
    """Compute Lennard-Jones interaction energy between two sets of atoms.

    Args:
        coords1: Coordinates of first atom set (N1 x 3) in Angstroms.
        c6_1: C6 coefficients for first set.
        c12_1: C12 coefficients for first set.
        coords2: Coordinates of second atom set (N2 x 3) in Angstroms.
        c6_2: C6 coefficients for second set.
        c12_2: C12 coefficients for second set.

    Returns:
        Lennard-Jones energy in wavenumbers (cm^-1).
    """
    # Compute pairwise distances
    diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff ** 2, axis=2))

    # Combine LJ coefficients (geometric mean)
    c6 = np.sqrt(c6_1[:, np.newaxis] * c6_2[np.newaxis, :])
    c12 = np.sqrt(c12_1[:, np.newaxis] * c12_2[np.newaxis, :])

    # Compute LJ energy
    r6 = distances ** 6
    r12 = distances ** 12
    energy = WVN_MOL_PER_KJ * np.sum(c12 / r12 - c6 / r6)

    return float(energy)


def compute_total_interaction_energy(
    coords1: NDArray[np.float64],
    charges1: NDArray[np.float64],
    c6_1: NDArray[np.float64],
    c12_1: NDArray[np.float64],
    coords2: NDArray[np.float64],
    charges2: NDArray[np.float64],
    c6_2: NDArray[np.float64],
    c12_2: NDArray[np.float64],
) -> tuple[float, float, float]:
    """Compute total interaction energy (Coulomb + LJ).

    Args:
        coords1: Coordinates of first atom set (N1 x 3).
        charges1: Charges of first atom set.
        c6_1: C6 LJ coefficients for first set.
        c12_1: C12 LJ coefficients for first set.
        coords2: Coordinates of second atom set (N2 x 3).
        charges2: Charges of second atom set.
        c6_2: C6 LJ coefficients for second set.
        c12_2: C12 LJ coefficients for second set.

    Returns:
        Tuple of (coulomb_energy, lj_energy, total_energy) in cm^-1.
    """
    ec = compute_coulomb_energy(coords1, charges1, coords2, charges2)
    lj = compute_lennard_jones_energy(coords1, c6_1, c12_1, coords2, c6_2, c12_2)

    return ec, lj, ec + lj


def compute_energy_along_normal(
    static_coords: NDArray[np.float64],
    static_charges: NDArray[np.float64],
    static_c6: NDArray[np.float64],
    static_c12: NDArray[np.float64],
    dynamic_coords: NDArray[np.float64],
    dynamic_charges: NDArray[np.float64],
    dynamic_c6: NDArray[np.float64],
    dynamic_c12: NDArray[np.float64],
    normal_vector: NDArray[np.float64],
    displacements: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Compute energy profile along a normal mode direction.

    Args:
        static_coords: Coordinates of static atoms (N_static x 3).
        static_charges: Charges of static atoms.
        static_c6: C6 coefficients for static atoms.
        static_c12: C12 coefficients for static atoms.
        dynamic_coords: Coordinates of dynamic atoms (N_dynamic x 3).
        dynamic_charges: Charges of dynamic atoms.
        dynamic_c6: C6 coefficients for dynamic atoms.
        dynamic_c12: C12 coefficients for dynamic atoms.
        normal_vector: Unit normal vector for displacement (3,).
        displacements: Array of displacement values in Angstroms.

    Returns:
        Tuple of (coulomb_energies, lj_energies, total_energies) arrays.
    """
    n_disp = len(displacements)
    ec = np.zeros(n_disp)
    lj = np.zeros(n_disp)

    for i, dr in enumerate(displacements):
        # Displace dynamic atoms along normal
        displaced_coords = dynamic_coords + dr * normal_vector

        ec[i], lj[i], _ = compute_total_interaction_energy(
            displaced_coords,
            dynamic_charges,
            dynamic_c6,
            dynamic_c12,
            static_coords,
            static_charges,
            static_c6,
            static_c12,
        )

    return ec, lj, ec + lj


def compute_plane_normal(
    p1: NDArray[np.float64],
    p2: NDArray[np.float64],
    p3: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Compute unit normal vector to plane defined by three points.

    Args:
        p1: First point (3,).
        p2: Second point (3,).
        p3: Third point (3,).

    Returns:
        Unit normal vector (3,).
    """
    v1 = p1 - p2
    v2 = p1 - p3
    normal = np.cross(v1, v2)
    return normal / np.linalg.norm(normal)


def get_residue_range(
    res_ids: NDArray[np.int64],
    target_residues: list[int],
) -> NDArray[np.bool_]:
    """Get boolean mask for atoms in specified residues.

    Args:
        res_ids: Array of residue IDs for each atom.
        target_residues: List of residue IDs to select.

    Returns:
        Boolean mask array.
    """
    mask = np.zeros(len(res_ids), dtype=bool)
    for res_id in target_residues:
        mask |= res_ids == res_id
    return mask


def compute_force_vectors(
    coords: NDArray[np.float64],
    forces: NDArray[np.float64],
    res_ids: NDArray[np.int64],
    residue_pairs: list[tuple[int, int]],
    scale: float = 7e-3,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Compute force vector origins and endpoints for visualization.

    Args:
        coords: Atomic coordinates (N x 3).
        forces: Force vectors (3 x N_residue_pairs).
        res_ids: Residue IDs for each atom.
        residue_pairs: List of (res1, res2) tuples to average over.
        scale: Scale factor for force visualization.

    Returns:
        Tuple of (origins, endpoints) arrays, each (N_pairs x 3).
    """
    n_pairs = len(residue_pairs)
    origins = np.zeros((n_pairs, 3))
    endpoints = np.zeros((n_pairs, 3))

    for i, (res1, res2) in enumerate(residue_pairs):
        mask1 = res_ids == res1
        mask2 = res_ids == res2

        # Compute center of residue pair
        center = np.concatenate([coords[mask1], coords[mask2]], axis=0).mean(axis=0)
        origins[i] = center

        # Compute endpoint
        endpoints[i] = center + scale * forces[:, i]

    return origins, endpoints
