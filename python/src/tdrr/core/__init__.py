"""Core spectral analysis functions.

This module contains fundamental signal processing and spectral analysis
functions including apodization, windowing, and bispectrum calculations.
"""

from tdrr.core.spectral import apodfun, lag_window, bispectrum, toeplitz_matrix
from tdrr.core.fitting import fit_exponential, fit_spectral_exponential
from tdrr.core.filters import lowpass_filter, bandpass_filter, wigner_distribution
from tdrr.core.matrix import circulant_matrix, symmetric_toeplitz

__all__ = [
    # Spectral
    "apodfun",
    "lag_window",
    "bispectrum",
    "toeplitz_matrix",
    # Fitting
    "fit_exponential",
    "fit_spectral_exponential",
    # Filters
    "lowpass_filter",
    "bandpass_filter",
    "wigner_distribution",
    # Matrix
    "circulant_matrix",
    "symmetric_toeplitz",
]
