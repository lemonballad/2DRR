"""Input/Output functions for various file formats.

Supports reading and writing:
- Gaussian quantum chemistry output files
- GROMACS molecular dynamics files
- PDB protein structure files
- Myoglobin energy and force constant data
"""

from tdrr.io.gaussian import GaussianData, extract_frequencies, read_energy
from tdrr.io.gromacs import (
    extract_xvg,
    read_gro,
    write_gro,
    read_write_angles,
    get_residue_ranges,
    displace_atoms_along_normal,
    generate_displaced_trajectory,
    read_force_xvg,
)
from tdrr.io.pdb import (
    PDBData,
    parse_pdb,
    write_b_factors,
    write_pdb_trajectory,
    generate_displaced_pdb_trajectory,
)
from tdrr.io.myoglobin import (
    read_gaussian_energies,
    process_myoglobin_energies,
    force_constant_to_frequency,
    extract_mode_frequencies,
    read_mbo_mbn_energies,
)

__all__ = [
    # Gaussian
    "GaussianData",
    "extract_frequencies",
    "read_energy",
    # GROMACS
    "extract_xvg",
    "read_gro",
    "write_gro",
    "read_write_angles",
    "get_residue_ranges",
    "displace_atoms_along_normal",
    "generate_displaced_trajectory",
    "read_force_xvg",
    # PDB
    "PDBData",
    "parse_pdb",
    "write_b_factors",
    "write_pdb_trajectory",
    "generate_displaced_pdb_trajectory",
    # Myoglobin
    "read_gaussian_energies",
    "process_myoglobin_energies",
    "force_constant_to_frequency",
    "extract_mode_frequencies",
    "read_mbo_mbn_energies",
]
