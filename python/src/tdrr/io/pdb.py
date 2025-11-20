"""PDB protein structure file I/O.

Functions for reading and writing Protein Data Bank format files.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import NDArray


@dataclass
class PDBData:
    """Container for PDB structure data.

    Attributes:
        record_name: Record type (ATOM, HETATM, etc.).
        atom_num: Serial numbers of atoms.
        atom_name: Element names.
        alt_loc: Alternate location indicators.
        res_name: Residue names.
        chain_id: Chain identifiers.
        res_num: Residue sequence numbers.
        x: X coordinates in Angstroms.
        y: Y coordinates in Angstroms.
        z: Z coordinates in Angstroms.
        occupancy: Occupancy values.
        beta_factor: Temperature factors.
        element: Element symbols.
        charge: Formal charges.
    """

    record_name: list[str] = field(default_factory=list)
    atom_num: NDArray[np.int64] = field(default_factory=lambda: np.array([], dtype=np.int64))
    atom_name: list[str] = field(default_factory=list)
    alt_loc: list[str] = field(default_factory=list)
    res_name: list[str] = field(default_factory=list)
    chain_id: list[str] = field(default_factory=list)
    res_num: NDArray[np.int64] = field(default_factory=lambda: np.array([], dtype=np.int64))
    x: NDArray[np.float64] = field(default_factory=lambda: np.array([], dtype=np.float64))
    y: NDArray[np.float64] = field(default_factory=lambda: np.array([], dtype=np.float64))
    z: NDArray[np.float64] = field(default_factory=lambda: np.array([], dtype=np.float64))
    occupancy: NDArray[np.float64] = field(default_factory=lambda: np.array([], dtype=np.float64))
    beta_factor: NDArray[np.float64] = field(default_factory=lambda: np.array([], dtype=np.float64))
    element: list[str] = field(default_factory=list)
    charge: list[str] = field(default_factory=list)


def parse_pdb(pdb_path: str | Path) -> PDBData:
    """Parse PDB file and extract atom coordinates.

    Args:
        pdb_path: Path to PDB file.

    Returns:
        PDBData containing all atom information.

    Examples:
        >>> data = parse_pdb("protein.pdb")
        >>> coords = np.column_stack([data.x, data.y, data.z])
    """
    pdb_path = Path(pdb_path)

    record_names: list[str] = []
    atom_nums: list[int] = []
    atom_names: list[str] = []
    alt_locs: list[str] = []
    res_names: list[str] = []
    chain_ids: list[str] = []
    res_nums: list[int] = []
    x_coords: list[float] = []
    y_coords: list[float] = []
    z_coords: list[float] = []
    occupancies: list[float] = []
    beta_factors: list[float] = []
    elements: list[str] = []
    charges: list[str] = []

    with open(pdb_path) as f:
        for line in f:
            # Only process ATOM and HETATM records
            record = line[0:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue

            # Check for valid coordinate line
            if len(line) < 54:
                continue

            # Check that coordinate fields are numeric
            coord_section = line[22:54]
            if any(c.isalpha() for c in coord_section):
                continue

            record_names.append(record)
            atom_nums.append(int(line[6:11]))
            atom_names.append(line[12:16].strip())
            alt_locs.append(line[16:17])
            res_names.append(line[17:20].strip())
            chain_ids.append(line[21:22])
            res_nums.append(int(line[22:26]))
            x_coords.append(float(line[30:38]))
            y_coords.append(float(line[38:46]))
            z_coords.append(float(line[46:54]))

            # Optional fields
            if len(line) >= 60:
                occupancies.append(float(line[54:60]))
            else:
                occupancies.append(1.0)

            if len(line) >= 66:
                beta_factors.append(float(line[60:66]))
            else:
                beta_factors.append(0.0)

            if len(line) >= 78:
                elements.append(line[76:78].strip())
            else:
                elements.append("")

            if len(line) >= 80:
                charges.append(line[78:80].strip())
            else:
                charges.append("")

    return PDBData(
        record_name=record_names,
        atom_num=np.array(atom_nums, dtype=np.int64),
        atom_name=atom_names,
        alt_loc=alt_locs,
        res_name=res_names,
        chain_id=chain_ids,
        res_num=np.array(res_nums, dtype=np.int64),
        x=np.array(x_coords, dtype=np.float64),
        y=np.array(y_coords, dtype=np.float64),
        z=np.array(z_coords, dtype=np.float64),
        occupancy=np.array(occupancies, dtype=np.float64),
        beta_factor=np.array(beta_factors, dtype=np.float64),
        element=elements,
        charge=charges,
    )


def write_b_factors(
    pdb_data: PDBData,
    b_factors: NDArray[np.float64],
    output_path: str | Path,
) -> None:
    """Write PDB file with updated B-factors.

    Args:
        pdb_data: Original PDB data.
        b_factors: New B-factor values for each atom.
        output_path: Path for output PDB file.
    """
    output_path = Path(output_path)

    lines = []
    for i in range(len(pdb_data.record_name)):
        line = (
            f"{pdb_data.record_name[i]:<6}"
            f"{pdb_data.atom_num[i]:5d} "
            f"{pdb_data.atom_name[i]:4s}"
            f"{pdb_data.alt_loc[i]:1s}"
            f"{pdb_data.res_name[i]:3s} "
            f"{pdb_data.chain_id[i]:1s}"
            f"{pdb_data.res_num[i]:4d}    "
            f"{pdb_data.x[i]:8.3f}"
            f"{pdb_data.y[i]:8.3f}"
            f"{pdb_data.z[i]:8.3f}"
            f"{pdb_data.occupancy[i]:6.2f}"
            f"{b_factors[i]:6.2f}"
        )
        lines.append(line)

    output_path.write_text("\n".join(lines) + "\n")
