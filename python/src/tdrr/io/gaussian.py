"""Gaussian quantum chemistry file I/O.

Functions for reading and writing Gaussian input/output files,
extracting spectroscopic properties like frequencies, force constants,
and normal modes.
"""

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray


@dataclass
class GaussianData:
    """Container for Gaussian spectroscopic output data.

    Attributes:
        frequencies: Vibrational frequencies in cm^-1.
        reduced_masses: Reduced masses in amu.
        force_constants: Force constants in mdyne/A.
        ir_intensities: IR intensities in KM/Mole.
        raman_activities: Raman activities in A^4/amu.
        depol_p: Depolarization ratios (perpendicular).
        depol_u: Depolarization ratios (unpolarized).
    """

    frequencies: NDArray[np.float64]
    reduced_masses: NDArray[np.float64]
    force_constants: NDArray[np.float64]
    ir_intensities: NDArray[np.float64]
    raman_activities: NDArray[np.float64]
    depol_p: NDArray[np.float64]
    depol_u: NDArray[np.float64]


def extract_frequencies(log_path: str | Path) -> GaussianData:
    """Extract vibrational data from Gaussian log file.

    Args:
        log_path: Path to Gaussian .log output file.

    Returns:
        GaussianData containing all extracted spectroscopic properties.
    """
    log_path = Path(log_path)
    content = log_path.read_text()

    # Patterns for extracting data
    freq_pattern = r"Frequencies --\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
    red_mass_pattern = r"Red\. masses --\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
    frc_const_pattern = r"Frc consts  --\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
    ir_pattern = r"IR Inten    --\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
    raman_pattern = r"Raman Activ --\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
    depol_p_pattern = r"Depolar \(P\) --\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
    depol_u_pattern = r"Depolar \(U\) --\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"

    def extract_values(pattern: str) -> list[float]:
        matches = re.findall(pattern, content)
        values = []
        for match in matches:
            values.extend([float(v) for v in match])
        return values

    frequencies = np.array(extract_values(freq_pattern), dtype=np.float64)
    reduced_masses = np.array(extract_values(red_mass_pattern), dtype=np.float64)
    force_constants = np.array(extract_values(frc_const_pattern), dtype=np.float64)
    ir_intensities = np.array(extract_values(ir_pattern), dtype=np.float64)
    raman_activities = np.array(extract_values(raman_pattern), dtype=np.float64)
    depol_p = np.array(extract_values(depol_p_pattern), dtype=np.float64)
    depol_u = np.array(extract_values(depol_u_pattern), dtype=np.float64)

    return GaussianData(
        frequencies=frequencies,
        reduced_masses=reduced_masses,
        force_constants=force_constants,
        ir_intensities=ir_intensities,
        raman_activities=raman_activities,
        depol_p=depol_p,
        depol_u=depol_u,
    )


def read_energy(
    list_path: str | Path,
    energy_pattern: str,
    distance_pattern: str,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Read energies and distances from Gaussian log files.

    Args:
        list_path: Path to file containing list of log file paths.
        energy_pattern: Pattern to match energy line.
        distance_pattern: Pattern to match distance line.

    Returns:
        Tuple of (energies, distances) arrays.
    """
    list_path = Path(list_path)
    lines = list_path.read_text().strip().split("\n")

    num_files = int(lines[0])
    energies = np.zeros(num_files, dtype=np.float64)
    distances = np.zeros(num_files, dtype=np.float64)

    for i, log_file in enumerate(lines[1 : num_files + 1]):
        log_content = Path(log_file.strip()).read_text()

        # Extract energy
        energy_match = re.search(f"{energy_pattern}\\s*([\\d.E+-]+)", log_content)
        if energy_match:
            energies[i] = float(energy_match.group(1))

        # Extract distance
        dist_match = re.search(f"{distance_pattern}\\s*([\\d.]+)", log_content)
        if dist_match:
            distances[i] = float(dist_match.group(1))

        # Check for convergence failure
        if "Convergence failure -- run terminated" in log_content:
            energies[i] = 0.0
            distances[i] = 0.0

    return energies, distances


def combine_data(*datasets: GaussianData) -> GaussianData:
    """Combine multiple GaussianData objects.

    Args:
        *datasets: Variable number of GaussianData objects.

    Returns:
        Combined GaussianData with concatenated arrays.
    """
    return GaussianData(
        frequencies=np.concatenate([d.frequencies for d in datasets]),
        reduced_masses=np.concatenate([d.reduced_masses for d in datasets]),
        force_constants=np.concatenate([d.force_constants for d in datasets]),
        ir_intensities=np.concatenate([d.ir_intensities for d in datasets]),
        raman_activities=np.concatenate([d.raman_activities for d in datasets]),
        depol_p=np.concatenate([d.depol_p for d in datasets]),
        depol_u=np.concatenate([d.depol_u for d in datasets]),
    )
