"""Trajectory generation and analysis functions.

Functions for generating synthetic molecular trajectories and analyzing
time-dependent spectroscopic data.
"""

import numpy as np
from numpy.typing import NDArray

from tdrr.core.spectral import apodfun


def generate_trajectory(
    freq_map: NDArray[np.float64],
    l_angles: NDArray[np.float64],
    r_angles: NDArray[np.float64],
    dt: float = 0.004,
) -> NDArray[np.float64]:
    """Generate frequency trajectory from dihedral angle time series.

    Args:
        freq_map: 2D frequency map (361 x 361 for 1-degree resolution).
        l_angles: Time series of left dihedral angles in degrees.
        r_angles: Time series of right dihedral angles in degrees.
        dt: Time step in picoseconds.

    Returns:
        Array of frequencies at each time point.
    """
    n_points = len(l_angles)
    frequencies = np.zeros(n_points, dtype=np.float64)

    # Convert angles to map indices
    l_idx = (np.round(l_angles).astype(int) + 180) % freq_map.shape[0]
    r_idx = (np.round(r_angles).astype(int) + 180) % freq_map.shape[1]

    for i in range(n_points):
        frequencies[i] = freq_map[l_idx[i], r_idx[i]]

    return frequencies


def compute_correlation(
    signal: NDArray[np.float64],
    normalize: bool = True,
) -> NDArray[np.float64]:
    """Compute autocorrelation function of signal.

    Args:
        signal: Input time series.
        normalize: If True, normalize by variance.

    Returns:
        Autocorrelation function (positive lags only).
    """
    n = len(signal)
    mean_signal = np.mean(signal)
    centered = signal - mean_signal

    # Compute full correlation
    corr = np.correlate(centered, centered, mode="full")

    # Take positive lags only
    corr = corr[n - 1 :]

    if normalize:
        corr = corr / corr[0]

    return corr


def compute_spectral_density(
    correlation: NDArray[np.float64],
    dt: float,
    apod_decay: float = 60.0,
    n_fft: int | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Compute spectral density from correlation function.

    Args:
        correlation: Autocorrelation function.
        dt: Time step in picoseconds.
        apod_decay: Apodization decay parameter.
        n_fft: FFT size (default: next power of 2).

    Returns:
        Tuple of (frequencies in cm^-1, spectral density).
    """
    n = len(correlation)
    time = np.arange(n) * dt

    # Apply apodization
    apod = apodfun(time, 0, 0, apod_decay, apod_decay, invert=False)
    apodized = correlation * apod

    # Compute FFT
    if n_fft is None:
        n_fft = 2 ** (int(np.ceil(np.log2(n))) + 2)

    spectrum = np.real(np.fft.fftshift(np.fft.fft(apodized, n_fft)))

    # Compute frequency axis in cm^-1
    # Convert from ps^-1 to cm^-1 using c = 3e10 cm/s
    c_cm_ps = 0.00003  # speed of light in cm/ps
    freq = np.fft.fftshift(np.fft.fftfreq(n_fft, dt)) * 2 * np.pi / c_cm_ps

    return freq, spectrum


def compute_2d_correlation(
    signal: NDArray[np.float64],
    dt: float,
    max_lag: int,
) -> NDArray[np.float64]:
    """Compute 2D correlation function for third-order spectroscopy.

    Args:
        signal: Input time series.
        dt: Time step.
        max_lag: Maximum lag in samples.

    Returns:
        2D correlation array of shape (max_lag, max_lag).
    """
    n = len(signal)
    mean_signal = np.mean(signal)
    centered = signal - mean_signal

    corr_2d = np.zeros((max_lag, max_lag), dtype=np.float64)

    for i1 in range(max_lag):
        for i2 in range(i1, max_lag):
            total = 0.0
            count = 0
            for t in range(n):
                t1 = t - i1
                t2 = t - i2
                if t1 >= 0 and t2 >= 0:
                    total += centered[t] * centered[t1] * centered[t2]
                    count += 1
            if count > 0:
                corr_2d[i1, i2] = total / count
                corr_2d[i2, i1] = corr_2d[i1, i2]

    return corr_2d


def fit_exponential_decay(
    time: NDArray[np.float64],
    correlation: NDArray[np.float64],
    n_exp: int = 2,
) -> dict[str, float]:
    """Fit exponential decay to correlation function.

    Args:
        time: Time array.
        correlation: Correlation values.
        n_exp: Number of exponential terms (1 or 2).

    Returns:
        Dictionary with fit parameters (a, b for single; a, b, c, d for double).
    """
    from scipy.optimize import curve_fit

    if n_exp == 1:

        def model(t: NDArray[np.float64], a: float, b: float) -> NDArray[np.float64]:
            return a * np.exp(b * t)

        popt, _ = curve_fit(model, time, correlation, p0=[1.0, -1.0])
        return {"a": popt[0], "b": popt[1]}

    else:

        def model(
            t: NDArray[np.float64], a: float, b: float, c: float, d: float
        ) -> NDArray[np.float64]:
            return a * np.exp(b * t) + c * np.exp(d * t)

        popt, _ = curve_fit(
            model, time, correlation, p0=[0.5, -10.0, 0.5, -1.0], maxfev=5000
        )
        return {"a": popt[0], "b": popt[1], "c": popt[2], "d": popt[3]}
