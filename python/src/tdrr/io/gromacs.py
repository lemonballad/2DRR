"""GROMACS molecular dynamics file I/O.

Functions for reading and writing GROMACS file formats
including XVG data files and GRO coordinate files.
"""

import re
from pathlib import Path

import numpy as np
from numpy.typing import NDArray


def extract_xvg(xvg_path: str | Path) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Extract data from GROMACS XVG file.

    Args:
        xvg_path: Path to XVG file.

    Returns:
        Tuple of (x_data, y_data) arrays.
    """
    xvg_path = Path(xvg_path)

    x_values: list[float] = []
    y_values: list[float] = []

    with open(xvg_path) as f:
        for line in f:
            # Skip comments and metadata
            if line.startswith(("#", "@")):
                continue

            parts = line.split()
            if len(parts) >= 2:
                x_values.append(float(parts[0]))
                y_values.append(float(parts[1]))

    return np.array(x_values), np.array(y_values)


def write_gro(
    output_path: str | Path,
    title: str,
    atom_names: list[str],
    res_names: list[str],
    res_nums: NDArray[np.int64],
    coords: NDArray[np.float64],
    box: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> None:
    """Write GROMACS GRO coordinate file.

    Args:
        output_path: Path for output GRO file.
        title: Title line for the file.
        atom_names: List of atom names.
        res_names: List of residue names.
        res_nums: Array of residue numbers.
        coords: Nx3 array of coordinates in nm.
        box: Box dimensions (x, y, z) in nm.
    """
    output_path = Path(output_path)
    num_atoms = len(atom_names)

    lines = [title, str(num_atoms)]

    for i in range(num_atoms):
        line = (
            f"{res_nums[i]:5d}"
            f"{res_names[i]:5s}"
            f"{atom_names[i]:5s}"
            f"{i + 1:5d}"
            f"{coords[i, 0]:8.3f}"
            f"{coords[i, 1]:8.3f}"
            f"{coords[i, 2]:8.3f}"
        )
        lines.append(line)

    lines.append(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}")

    output_path.write_text("\n".join(lines) + "\n")


def read_gro(gro_path: str | Path) -> tuple[
    list[str],
    list[str],
    NDArray[np.int64],
    NDArray[np.float64],
    tuple[float, float, float],
]:
    """Read GROMACS GRO coordinate file.

    Args:
        gro_path: Path to GRO file.

    Returns:
        Tuple of (atom_names, res_names, res_nums, coords, box).
    """
    gro_path = Path(gro_path)
    lines = gro_path.read_text().strip().split("\n")

    num_atoms = int(lines[1])

    atom_names: list[str] = []
    res_names: list[str] = []
    res_nums: list[int] = []
    coords: list[list[float]] = []

    for i in range(2, 2 + num_atoms):
        line = lines[i]
        res_nums.append(int(line[0:5]))
        res_names.append(line[5:10].strip())
        atom_names.append(line[10:15].strip())
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        coords.append([x, y, z])

    # Parse box line
    box_line = lines[2 + num_atoms].split()
    box = (float(box_line[0]), float(box_line[1]), float(box_line[2]))

    return (
        atom_names,
        res_names,
        np.array(res_nums, dtype=np.int64),
        np.array(coords, dtype=np.float64),
        box,
    )


def read_write_angles(
    input_path: str | Path,
    output_path: str | Path,
) -> NDArray[np.float64]:
    """Read dihedral angles from GROMACS output and write processed data.

    Args:
        input_path: Path to input XVG file with angle data.
        output_path: Path for output processed angles.

    Returns:
        Array of dihedral angles in degrees.
    """
    time, angles = extract_xvg(input_path)

    # Convert to degrees and normalize to [-180, 180]
    angles_deg = np.degrees(angles)
    angles_deg = np.mod(angles_deg + 180, 360) - 180

    # Write to output
    output_path = Path(output_path)
    np.savetxt(output_path, np.column_stack([time, angles_deg]), fmt="%.6f")

    return angles_deg
