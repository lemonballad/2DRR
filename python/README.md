# TDRR - Two-Dimensional Resonant Raman Spectroscopy

Python package for analyzing molecular vibrational modes and spectral properties.

## Installation

```bash
pip install -e .
```

For development:

```bash
pip install -e ".[dev]"
```

## Usage

### Core Spectral Functions

```python
import numpy as np
from tdrr.core.spectral import apodfun, lag_window, bispectrum

# Create apodization window
x = np.linspace(0, 1000, 1001)
window = apodfun(x, 500.0, 50.0, 50.0, 100.0)

# Generate lag window
w = lag_window(100, "hamming")

# Compute bispectrum
signal = np.random.randn(1000, 1)
bisp, freq, cum, lag = bispectrum(signal, sample_rate=1000, max_lag=50)
```

### Reading Gaussian Output

```python
from tdrr.io.gaussian import extract_frequencies

data = extract_frequencies("output.log")
print(f"Frequencies: {data.frequencies}")
print(f"Force constants: {data.force_constants}")
```

### Reading PDB Files

```python
from tdrr.io.pdb import parse_pdb

pdb_data = parse_pdb("protein.pdb")
coords = np.column_stack([pdb_data.x, pdb_data.y, pdb_data.z])
```

## Testing

```bash
pytest
```

## Project Structure

```
src/tdrr/
├── core/           # Core spectral analysis functions
│   └── spectral.py # apodfun, lag_window, bispectrum
├── io/             # File I/O
│   ├── gaussian.py # Gaussian file parsing
│   └── pdb.py      # PDB file parsing
├── analysis/       # Analysis modules
└── plotting/       # Visualization
```

## License

MIT
