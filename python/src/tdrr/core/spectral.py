"""Core spectral analysis functions.

Provides fundamental signal processing operations for spectroscopy analysis
including apodization windows, lag windows, and bispectrum calculations.
"""

import numpy as np
from numpy.typing import NDArray


def apodfun(
    x: NDArray[np.float64],
    x0: float,
    g_low: float,
    g_high: float,
    width: float,
    invert: bool = False,
) -> NDArray[np.float64]:
    """Generate asymmetric bandpass or notch filter.

    Creates an apodization window with Gaussian decay on each side
    and a flat top region.

    Args:
        x: Input domain array (frequencies or time points).
        x0: Center position of the filter.
        g_low: Gaussian decay width on the low side.
        g_high: Gaussian decay width on the high side.
        width: Total width of the flat top region.
        invert: If True, create notch filter instead of bandpass.

    Returns:
        Apodization window array with same shape as x.

    Examples:
        >>> import numpy as np
        >>> x = np.linspace(0, 100, 1000)
        >>> window = apodfun(x, 50.0, 10.0, 10.0, 20.0)
        >>> window[500]  # Center should be ~1.0
        1.0
    """
    x = np.asarray(x, dtype=np.float64)

    lo_center = x0 - width / 2
    hi_center = x0 + width / 2

    result = np.zeros_like(x)

    # Region 1: below low center - Gaussian decay
    mask_low = x < lo_center
    result[mask_low] = np.exp(-0.5 * (x[mask_low] - lo_center) ** 2 / g_low ** 2)

    # Region 2: flat top
    mask_flat = (x >= lo_center) & (x <= hi_center)
    result[mask_flat] = 1.0

    # Region 3: above high center - Gaussian decay
    mask_high = x > hi_center
    result[mask_high] = np.exp(-0.5 * (x[mask_high] - hi_center) ** 2 / g_high ** 2)

    if invert:
        result = 1.0 - result

    return result


def lag_window(
    lag: int,
    window: str = "uniform",
) -> NDArray[np.float64]:
    """Compute lag window function for spectral estimation.

    Args:
        lag: Number of lag elements (must be >= 1).
        window: Window type. Options:
            - "uniform" or "u": Uniform window
            - "sasaki" or "s": Sasaki window
            - "priestley" or "p": Priestley window
            - "parzen" or "pa": Parzen window
            - "hamming" or "h": Hamming window
            - "gaussian" or "g": Gaussian window
            - "daniell" or "d": Daniell window

    Returns:
        Lag window array of length `lag`.

    References:
        C. L. Nikias, A. P. Petropulu, Higher-Order Spectra Analysis:
        A Nonlinear Signal Processing Framework, PTR Prentice Hall, 1993.

    Examples:
        >>> w = lag_window(4, "parzen")
        >>> len(w)
        4
    """
    lag = int(lag)

    if lag == 1:
        return np.array([1.0])

    window = window.lower()
    lag1 = lag - 1

    # Resolve window aliases
    window_map = {
        "u": "uniform",
        "s": "sasaki",
        "p": "priestley",
        "pa": "parzen",
        "h": "hamming",
        "g": "gaussian",
        "d": "daniell",
    }
    window = window_map.get(window, window)

    if window == "uniform":
        return np.ones(lag)

    elif window == "sasaki":
        wind_lag = np.arange(lag) / lag1
        return np.sin(np.pi * wind_lag) / np.pi + np.cos(np.pi * wind_lag) * (1 - wind_lag)

    elif window == "priestley":
        wind_lag = np.arange(1, lag) / lag1
        result = np.zeros(lag)
        result[0] = 1.0
        result[1:] = (np.sin(np.pi * wind_lag) / np.pi / wind_lag -
                      np.cos(np.pi * wind_lag)) * 3 / np.pi ** 2 / wind_lag ** 2
        return result

    elif window == "parzen":
        fix_lag12 = lag1 // 2
        fix_lag121 = fix_lag12 + 1
        wind = np.zeros(lag)

        wind_lag0 = np.arange(fix_lag121) / lag1
        wind[:fix_lag121] = 1 - (1 - wind_lag0) * wind_lag0 ** 2 * 6

        wind_lag1 = 1 - np.arange(fix_lag121, lag) / lag1
        wind[fix_lag121:] = wind_lag1 ** 3 * 2

        return wind

    elif window == "hamming":
        return 0.54 + 0.46 * np.cos(np.pi * np.arange(lag) / lag1)

    elif window == "gaussian":
        wind = np.zeros(lag)
        wind[0] = 1.0
        wind[1:lag-1] = 0.5 * np.erfc(
            (np.arange(1, lag-1) / lag1 - 0.5) * 8 / np.sqrt(2)
        )
        wind[lag-1] = 0.0
        return wind

    elif window == "daniell":
        wind = np.zeros(lag)
        wind[0] = 1.0
        wind_lag = np.arange(1, lag) / lag1
        wind[1:] = np.sin(np.pi * wind_lag) / np.pi / wind_lag
        return wind

    else:
        raise ValueError(f"Unknown window type: {window}")


def toeplitz_matrix(
    column: NDArray[np.float64],
    row: NDArray[np.float64] | None = None,
) -> NDArray[np.float64]:
    """Generate Toeplitz matrix from column and row vectors.

    Args:
        column: First column of the matrix.
        row: First row of the matrix. If None, uses conjugate of column.

    Returns:
        Toeplitz matrix.
    """
    from scipy.linalg import toeplitz

    column = np.asarray(column)
    if row is None:
        row = np.conj(column)
    else:
        row = np.asarray(row)

    return toeplitz(column, row)


def bispectrum(
    signal: NDArray[np.float64],
    sample_rate: float = 1.0,
    max_lag: int = 0,
    window: str = "none",
    scale: str = "biased",
) -> tuple[
    NDArray[np.complex128],
    NDArray[np.float64],
    NDArray[np.complex128],
    NDArray[np.float64],
]:
    """Compute auto bispectrum and 3rd order cumulant.

    Args:
        signal: Input signal matrix (samples in rows, records in columns).
        sample_rate: Signal sample rate.
        max_lag: Maximum lag for cumulant computation.
        window: Lag window type ("none", "uniform", "sasaki", etc.).
        scale: Scale type ("biased" or "unbiased").

    Returns:
        Tuple of (bispectrum, frequency, cumulant, lag).

    References:
        C. L. Nikias, A. P. Petropulu, Higher-Order Spectra Analysis:
        A Nonlinear Signal Processing Framework, PTR Prentice Hall, 1993.
    """
    signal = np.atleast_2d(signal)
    if signal.shape[0] == 1:
        signal = signal.T

    sample, record = signal.shape

    # Compute lag and frequency vectors
    lag_index = np.arange(-max_lag, max_lag + 1)
    lag = lag_index / sample_rate

    if max_lag:
        freq = lag_index / max_lag / 2 * sample_rate
    else:
        freq = np.array([0.0])
        window = "none"
        scale = "biased"

    max_lag1 = max_lag + 1
    max_lag2 = max_lag * 2
    max_lag21 = max_lag2 + 1

    # Subtract mean from signal
    signal = signal - signal.mean(axis=0)

    # Initialize cumulant matrix
    cum = np.zeros((max_lag21, max_lag21), dtype=np.complex128)

    # Compute cumulant for each record
    for k in range(record):
        sig = signal[:, k]
        trfl_sig = sig[::-1]

        # Build Toeplitz matrix
        col = np.concatenate([sig[sample - max_lag:], np.zeros(max_lag)])
        row = np.concatenate([[trfl_sig[max_lag]], np.zeros(max_lag2)])
        toep_sig = toeplitz_matrix(col, row)

        # Compute cumulant contribution
        cum += toep_sig * trfl_sig.reshape(1, -1)[:, :max_lag21] @ toep_sig.T

    cum = cum / record

    # Apply scaling
    if scale.lower().startswith("b"):
        cum = cum / sample
    else:
        # Unbiased scaling (simplified)
        cum = cum / sample

    # Generate lag window
    if window.lower() == "none" or window.lower() == "n":
        wind = np.ones((max_lag21, max_lag21))
    else:
        wind1d = lag_window(max_lag1, window)
        wind1d_full = np.concatenate([wind1d[::-1][:-1], wind1d])
        wind = np.outer(wind1d_full, wind1d_full)

    # Compute bispectrum
    bisp = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(cum * wind)))

    return bisp, freq, cum, lag
