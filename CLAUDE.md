# CLAUDE.md

## Project Overview

**2DRR** (Two-Dimensional Resonant Raman Spectroscopy) is an archived computational chemistry research project that analyzes molecular vibrational modes and spectral properties of heme-containing proteins (myoglobin in MBO/MBN states).

**Status**: Archived research project - repository serves as an archive.

## Technology Stack

- **Primary Language**: MATLAB 2016
- **Dependencies**: MATLAB Standard Library only (no external toolboxes required)
- **Build System**: None - scripts run directly in MATLAB

## Code Structure

```
matlab/
├── sub_task_1/    # Data reading and preprocessing
│                  # Reads Gaussian output files, extracts frequencies, force constants
├── sub_task_2/    # Force analysis and visualization
│                  # GROMACS MD trajectories, force vectors, energy plots
└── sub_task_3/    # Main spectroscopic analysis
                   # 2D spectral plots, mode mapping, trajectory generation
```

## Key File Formats

- **Gaussian**: `.log` files (quantum chemistry calculations)
- **GROMACS**: `.gro`, `.xvg` (molecular dynamics)
- **Protein structures**: `.pdb` (Protein Data Bank)
- **Data**: `.txt`, `.mat` (MATLAB native)

## Coding Conventions

### Naming
- Scripts: UPPER_SNAKE_CASE or MixedCase (`Make_map.m`, `PROCESSOR.m`)
- Utility functions: lowercase (`apodfun.m`, `lagwind.m`)
- Classes: PascalCase (`GauData.m`)

### Configuration Pattern
Scripts start with hardcoded parameters:
```matlab
clear all
path_source='C:\...';
molecule_name='MBO';
% Define parameters before main logic
```

### Common Variables
- `freq`, `redM`, `frcConsts` - spectroscopic properties
- `MBN`, `MBO` - myoglobin binding states
- `L_indices`, `R_indices` - left/right dihedral angles

### Data Structures
- 3D arrays: `(modes × left_angle × right_angle)` for spectral maps
- 2D matrices: Time-domain correlation functions
- Vectors: Frequency data, molecular coordinates

## Key Algorithms

- **Mode mapping**: 2D maps of vibrational modes across dihedral angle space
- **FFT-based spectral analysis**: Frequency domain transformations
- **Trajectory generation**: Synthetic molecular trajectories from MD data
- **Higher-order spectral analysis**: Bispectrum calculations

## Running Code

1. Open MATLAB 2016+
2. Navigate to appropriate sub_task directory
3. Update hardcoded paths at script top to match your system
4. Run scripts directly (e.g., `PROCESSOR.m`, `SD2DPlots.m`)

## Important Notes

- **Windows paths**: Scripts contain hardcoded Windows paths that need updating
- **No tests**: No unit test framework exists
- **Research code**: Contains debug code, iterations, and commented alternatives
- **Data-heavy**: sub_task_2 and sub_task_3 contain large data files (~77MB and ~72MB)

## Key Files

- `/matlab/sub_task_1/GauData.m` - Core data class definition
- `/matlab/sub_task_3/PROCESSOR.m` - Main processing pipeline
- `/matlab/sub_task_3/SD2DPlots.m` - 2D spectral plot generation
- `/matlab/sub_task_3/export_fig.m` - External utility for figure export

## Domain Knowledge

- **Vibrational spectroscopy**: Analysis of molecular vibrations
- **Dihedral angles**: Rotation angles (72 discrete values for left/right)
- **Normal modes**: Quantum chemistry vibrational mode analysis
- **Resonant Raman**: 2D correlation of vibrational frequencies
