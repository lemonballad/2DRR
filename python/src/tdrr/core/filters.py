"""Signal filtering and processing functions.

Low-pass filters, frequency domain operations, and related utilities.
"""

import numpy as np
from numpy.typing import NDArray

from tdrr.core.spectral import apodfun


def lowpass_filter(
    signal: NDArray[np.float64],
    dt: float,
    cutoff_freq: float,
    filter_type: str = "gaussian",
) -> NDArray[np.float64]:
    """Apply low-pass filter to signal in frequency domain.

    Args:
        signal: Input time-domain signal.
        dt: Time step in picoseconds.
        cutoff_freq: Cutoff frequency in cm^-1.
        filter_type: Filter type ("gaussian" or "ideal").

    Returns:
        Filtered signal in time domain.
    """
    n = len(signal)
    c = 0.00003  # Speed of light in cm/ps

    # Compute frequency axis
    freq = np.fft.fftshift(np.fft.fftfreq(n, dt)) * 2 * np.pi / c

    # Create filter
    if filter_type == "gaussian":
        # Gaussian low-pass
        filt = 1 - np.exp(-freq ** 2 / (2 * cutoff_freq ** 2))
    elif filter_type == "ideal":
        # Ideal low-pass
        filt = np.abs(freq) <= cutoff_freq
    else:
        raise ValueError(f"Unknown filter type: {filter_type}")

    # Apply filter in frequency domain
    signal_fft = np.fft.fft(signal)
    filtered_fft = np.fft.ifftshift(np.fft.fftshift(signal_fft) * filt)
    filtered = np.fft.ifft(filtered_fft)

    return np.real(filtered)


def bandpass_filter(
    signal: NDArray[np.float64],
    dt: float,
    low_freq: float,
    high_freq: float,
    width: float = 10.0,
) -> NDArray[np.float64]:
    """Apply bandpass filter to signal.

    Args:
        signal: Input time-domain signal.
        dt: Time step in picoseconds.
        low_freq: Low cutoff frequency in cm^-1.
        high_freq: High cutoff frequency in cm^-1.
        width: Filter transition width.

    Returns:
        Filtered signal in time domain.
    """
    n = len(signal)
    c = 0.00003  # Speed of light in cm/ps

    # Compute frequency axis
    freq = np.fft.fftshift(np.fft.fftfreq(n, dt)) * 2 * np.pi / c

    # Create bandpass filter using apodfun
    center = (low_freq + high_freq) / 2
    band_width = high_freq - low_freq
    filt = apodfun(np.abs(freq), center, width, width, band_width, invert=False)

    # Apply filter
    signal_fft = np.fft.fft(signal)
    filtered_fft = np.fft.ifftshift(np.fft.fftshift(signal_fft) * filt)
    filtered = np.fft.ifft(filtered_fft)

    return np.real(filtered)


def compute_power_spectrum(
    signal: NDArray[np.float64],
    dt: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Compute power spectrum of signal.

    Args:
        signal: Input time-domain signal.
        dt: Time step in picoseconds.

    Returns:
        Tuple of (frequencies in cm^-1, power spectrum).
    """
    n = len(signal)
    c = 0.00003  # Speed of light in cm/ps

    # Compute FFT
    signal_fft = np.fft.fftshift(np.fft.fft(signal))
    power = np.abs(signal_fft) ** 2

    # Compute frequency axis
    freq = np.fft.fftshift(np.fft.fftfreq(n, dt)) * 2 * np.pi / c

    return freq, power


def wigner_distribution(
    x: NDArray[np.float64],
    p: NDArray[np.float64],
    mass: float,
    omega: float,
    n: int = 0,
    hbar: float = 5300.0,
) -> NDArray[np.float64]:
    """Compute Wigner quasi-probability distribution.

    Args:
        x: Position grid.
        p: Momentum grid.
        mass: Particle mass.
        omega: Angular frequency in cm^-1.
        n: Quantum number.
        hbar: Reduced Planck constant in cm^-1 * ps.

    Returns:
        Wigner distribution on (x, p) grid.
    """
    from scipy.special import eval_genlaguerre

    c = 0.00003  # Speed of light in cm/ps

    # Create meshgrid if 1D arrays provided
    if x.ndim == 1 and p.ndim == 1:
        X, P = np.meshgrid(x, p, indexing='ij')
    else:
        X, P = x, p

    # Compute classical energy
    u = mass * (2 * np.pi * c * omega) ** 2 * X ** 2 / 2 + P ** 2 / (2 * mass)

    # Compute Wigner function for nth state
    # W_n = (-1)^n / (pi * hbar) * L_n(4u/hbar/omega) * exp(-2u/hbar/omega)
    prefactor = (-1) ** n / (np.pi * hbar)
    laguerre = eval_genlaguerre(n, 0, 4 * u / (hbar * omega))
    exponential = np.exp(-2 * u / (hbar * omega))

    return prefactor * laguerre * exponential
