"""Core spectral analysis functions.

This module contains fundamental signal processing and spectral analysis
functions including apodization, windowing, and bispectrum calculations.
"""

from tdrr.core.spectral import apodfun, lag_window, bispectrum

__all__ = ["apodfun", "lag_window", "bispectrum"]
