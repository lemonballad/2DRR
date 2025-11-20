"""Spectral plotting functions.

Matplotlib-based visualization for spectroscopy analysis.
"""

import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes


def plot_2d_spectrum(
    freq1: NDArray[np.float64],
    freq2: NDArray[np.float64],
    spectrum: NDArray[np.float64],
    n_contours: int = 50,
    colormap: str = "jet",
    freq_range: tuple[float, float] | None = None,
    title: str = "",
    xlabel: str = r"$\omega_1$ (cm$^{-1}$)",
    ylabel: str = r"$\omega_2$ (cm$^{-1}$)",
) -> tuple[Figure, Axes]:
    """Plot 2D spectral contour map.

    Args:
        freq1: First frequency axis.
        freq2: Second frequency axis.
        spectrum: 2D spectral data.
        n_contours: Number of contour levels.
        colormap: Matplotlib colormap name.
        freq_range: Optional (min, max) frequency limits.
        title: Plot title.
        xlabel: X-axis label.
        ylabel: Y-axis label.

    Returns:
        Tuple of (figure, axes).
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    contour = ax.contourf(
        freq1, freq2, spectrum.T, n_contours, cmap=colormap
    )
    fig.colorbar(contour, ax=ax)

    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_title(title, fontsize=18)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=14)

    if freq_range is not None:
        ax.set_xlim(freq_range)
        ax.set_ylim(freq_range)

    return fig, ax


def plot_correlation(
    time: NDArray[np.float64],
    correlation: NDArray[np.float64],
    fit_curve: NDArray[np.float64] | None = None,
    title: str = "Correlation Function",
    xlabel: str = "t (ps)",
    ylabel: str = "C(t)",
    log_scale: bool = False,
) -> tuple[Figure, Axes]:
    """Plot correlation function with optional fit.

    Args:
        time: Time array.
        correlation: Correlation values.
        fit_curve: Optional fitted curve.
        title: Plot title.
        xlabel: X-axis label.
        ylabel: Y-axis label.
        log_scale: Use logarithmic x-axis.

    Returns:
        Tuple of (figure, axes).
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    if log_scale:
        ax.semilogx(time, correlation, "b-", linewidth=2, label="Data")
        if fit_curve is not None:
            ax.semilogx(time, fit_curve, "r--", linewidth=2, label="Fit")
    else:
        ax.plot(time, correlation, "b-", linewidth=2, label="Data")
        if fit_curve is not None:
            ax.plot(time, fit_curve, "r--", linewidth=2, label="Fit")

    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_title(title, fontsize=18)
    ax.tick_params(labelsize=14)
    ax.set_aspect("auto")

    if fit_curve is not None:
        ax.legend(fontsize=12)

    return fig, ax


def plot_spectral_density(
    freq: NDArray[np.float64],
    spectrum: NDArray[np.float64],
    freq_range: tuple[float, float] = (0, 100),
    title: str = "Spectral Density",
    xlabel: str = r"$\omega$ (cm$^{-1}$)",
    ylabel: str = r"C($\omega$)",
    normalize: bool = True,
) -> tuple[Figure, Axes]:
    """Plot spectral density.

    Args:
        freq: Frequency array.
        spectrum: Spectral density values.
        freq_range: Frequency range to display.
        title: Plot title.
        xlabel: X-axis label.
        ylabel: Y-axis label.
        normalize: Normalize to maximum value.

    Returns:
        Tuple of (figure, axes).
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    if normalize:
        spectrum = spectrum / np.max(np.abs(spectrum))

    ax.plot(freq, spectrum, "b-", linewidth=2)

    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_title(title, fontsize=18)
    ax.set_xlim(freq_range)
    ax.tick_params(labelsize=14)

    return fig, ax


def plot_energy_profile(
    displacement: NDArray[np.float64],
    energy: NDArray[np.float64],
    energy_type: str = "Total",
    reference_lines: list[float] | None = None,
    xlabel: str = "dr (Å)",
    ylabel: str = r"E (cm$^{-1}$)",
) -> tuple[Figure, Axes]:
    """Plot energy profile along displacement coordinate.

    Args:
        displacement: Displacement values.
        energy: Energy values.
        energy_type: Type of energy for title.
        reference_lines: Optional x-positions for vertical reference lines.
        xlabel: X-axis label.
        ylabel: Y-axis label.

    Returns:
        Tuple of (figure, axes).
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.plot(displacement, energy, "b-", linewidth=4)

    if reference_lines is not None:
        for x_pos in reference_lines:
            ax.axvline(x=x_pos, color="k", linestyle=":", linewidth=2)

    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_title(f"{energy_type} Energy", fontsize=18)
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(labelsize=14)

    return fig, ax
