"""Myoglobin energy and force constant analysis.

Functions for reading Gaussian log files, extracting energy data,
and converting force constants to frequencies for myoglobin analysis.
"""

import numpy as np
from numpy.typing import NDArray
from pathlib import Path


# Physical constants for unit conversion
HARTREE_TO_CM = 219474.6313705  # Hartree to cm^-1
U_KG = 1.66e-27  # Atomic mass unit in kg
C_CMS = 3e10  # Speed of light in cm/s
MDYN_FACTOR = 1e2  # mDyn/Angstrom to SI conversion


def read_gaussian_energies(
    list_file: str | Path,
    energy_pattern: str,
    dihedral_pattern: str,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Read energies and dihedral angles from Gaussian log files.

    Reads a list file containing paths to Gaussian .log files,
    extracts SCF energies and dihedral angle values.

    Args:
        list_file: Path to file containing list of Gaussian .log paths.
        energy_pattern: Pattern to match energy line (e.g., ' SCF Done:  E(UB3LYP) =  ').
        dihedral_pattern: Pattern to match dihedral angle (e.g., ' Input=MBN_').

    Returns:
        Tuple of (energies in Hartree, dihedral angles in degrees).
    """
    list_path = Path(list_file)
    fail_pattern = " Convergence failure -- run terminated."

    with open(list_path, "r") as f:
        lines = f.readlines()

    # First line contains number of files
    num_files = int(lines[0].strip())
    energies = np.zeros(num_files, dtype=np.float64)
    dihedrals = np.zeros(num_files, dtype=np.float64)

    # Process each log file
    for idx, log_path in enumerate(lines[1 : num_files + 1]):
        log_path = log_path.strip()
        if not log_path:
            continue

        with open(log_path, "r") as f:
            for line in f:
                # Check for convergence failure
                if fail_pattern in line:
                    energies[idx] = 0.0
                    dihedrals[idx] = 0.0
                    break

                # Extract energy
                if energy_pattern in line:
                    # Energy value follows the pattern
                    after_pattern = line.split(energy_pattern)[1]
                    energy_str = after_pattern.split()[0]
                    energies[idx] = float(energy_str)

                # Extract dihedral angle
                if dihedral_pattern in line:
                    after_pattern = line.split(dihedral_pattern)[1]
                    dihedral_str = after_pattern.split()[0]
                    # Remove any trailing characters like '.log'
                    dihedral_str = "".join(c for c in dihedral_str if c.isdigit() or c in ".-")
                    if dihedral_str:
                        dihedrals[idx] = float(dihedral_str)

    return energies, dihedrals


def process_myoglobin_energies(
    energies: NDArray[np.float64],
    dihedrals: NDArray[np.float64],
    reference_energy: float = 0.0,
    convert_to_cm: bool = False,
    ensure_ascending: bool = True,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Process raw myoglobin energy data.

    Filters zero values, applies reference offset, optionally converts
    units, and ensures ascending dihedral order.

    Args:
        energies: Raw energies in Hartree.
        dihedrals: Dihedral angles in degrees.
        reference_energy: Reference energy to subtract (Hartree).
        convert_to_cm: If True, convert energies to cm^-1.
        ensure_ascending: If True, sort by ascending dihedral angle.

    Returns:
        Tuple of (processed energies, filtered dihedrals).
    """
    # Filter out zero values (failed calculations)
    mask = energies != 0
    e_filtered = energies[mask].copy()
    d_filtered = dihedrals[mask].copy()

    # Apply reference energy offset
    if reference_energy != 0:
        e_filtered = e_filtered + reference_energy

    # Convert to wavenumbers if requested
    if convert_to_cm:
        e_filtered = e_filtered * HARTREE_TO_CM

    # Ensure ascending dihedral order
    if ensure_ascending and len(d_filtered) > 1:
        if d_filtered[0] > d_filtered[-1]:
            d_filtered = d_filtered[::-1]
            e_filtered = e_filtered[::-1]

    return e_filtered, d_filtered


def force_constant_to_frequency(
    force_constant: NDArray[np.float64],
    reduced_mass: float,
) -> NDArray[np.float64]:
    """Convert force constants to vibrational frequencies.

    Uses the harmonic oscillator relation: omega = sqrt(k/m) / (2*pi*c)

    Args:
        force_constant: Force constants in mDyn/Angstrom.
        reduced_mass: Reduced mass in atomic mass units.

    Returns:
        Frequencies in cm^-1.
    """
    # omega = (1/c) * sqrt(k * conversion / (m * u_kg))
    frequencies = (
        (1 / C_CMS)
        * np.sqrt(force_constant * MDYN_FACTOR / (reduced_mass * U_KG))
    )
    return frequencies


def extract_mode_frequencies(
    force_constants: NDArray[np.float64],
    dihedral_angles: NDArray[np.float64],
    dihedral_index: int,
    bond_indices: tuple[int, int],
    angle_index: int,
    reduced_mass_stretch: float = 6.0,
    reduced_mass_bend: float = 2.0,
) -> NDArray[np.float64]:
    """Extract vibrational mode frequencies from force constant data.

    Converts force constants for stretching and bending modes to frequencies
    at specified dihedral angle positions.

    Args:
        force_constants: 3D array of force constants (dihedral, mode, angle).
        dihedral_angles: Array of dihedral angle values.
        dihedral_index: Index into first dimension for dihedral selection.
        bond_indices: Tuple of (bond1, bond2) mode indices.
        angle_index: Index for bending mode.
        reduced_mass_stretch: Reduced mass for stretching modes (amu).
        reduced_mass_bend: Reduced mass for bending modes (amu).

    Returns:
        Array of shape (n_angles, 4) with columns:
        [dihedral_angle, freq_bond1, freq_bond2, freq_angle].
    """
    n_angles = len(dihedral_angles)
    data = np.zeros((n_angles, 4), dtype=np.float64)

    data[:, 0] = dihedral_angles

    # Extract force constants for each mode
    fc_bond1 = force_constants[dihedral_index, bond_indices[0], :]
    fc_bond2 = force_constants[dihedral_index, bond_indices[1], :]
    fc_angle = force_constants[dihedral_index, angle_index, :]

    # Convert to frequencies
    data[:, 1] = force_constant_to_frequency(fc_bond1, reduced_mass_stretch)
    data[:, 2] = force_constant_to_frequency(fc_bond2, reduced_mass_stretch)
    data[:, 3] = force_constant_to_frequency(fc_angle, reduced_mass_bend)

    return data


def read_mbo_mbn_energies(
    list_file_1: str | Path,
    list_file_2: str | Path,
    state: str = "MBN",
    reference_energies: tuple[float, float] = (0.0, 0.0),
    convert_to_cm: bool = True,
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
]:
    """Read and process MBO/MBN myoglobin energy data.

    Convenience function for reading paired energy datasets from
    different SCF methods (e.g., UB3LYP and ROB3LYP).

    Args:
        list_file_1: Path to first list file.
        list_file_2: Path to second list file.
        state: Binding state ('MBN' or 'MBO').
        reference_energies: Tuple of (ref1, ref2) reference energies.
        convert_to_cm: If True, convert to cm^-1.

    Returns:
        Tuple of (energies_1, dihedrals_1, energies_2, dihedrals_2).
    """
    # Define patterns based on state
    if state == "MBN":
        pattern_1 = (" SCF Done:  E(UB3LYP) =  ", f" Input={state}_")
        pattern_2 = (" SCF Done:  E(ROB3LYP) =  ", f" Input=MBO_")
    else:  # MBO
        pattern_1 = (" SCF Done:  E(RB3LYP) =  ", f" Input=1{state}_")
        pattern_2 = (" SCF Done:  E(UB3LYP) =  ", f" Input=1MBO_")

    # Read raw data
    e1_raw, d1_raw = read_gaussian_energies(
        list_file_1, pattern_1[0], pattern_1[1]
    )
    e2_raw, d2_raw = read_gaussian_energies(
        list_file_2, pattern_2[0], pattern_2[1]
    )

    # Process data
    e1, d1 = process_myoglobin_energies(
        e1_raw, d1_raw, reference_energies[0], convert_to_cm
    )
    e2, d2 = process_myoglobin_energies(
        e2_raw, d2_raw, reference_energies[1], convert_to_cm
    )

    return e1, d1, e2, d2
