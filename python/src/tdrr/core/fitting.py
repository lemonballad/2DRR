"""Curve fitting functions for spectral analysis.

Functions for fitting exponential decay models to correlation functions
and spectral data.
"""

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import curve_fit


def exponential_model(
    t: NDArray[np.float64],
    *params: float,
) -> NDArray[np.float64]:
    """Evaluate sum of exponentials model.

    Args:
        t: Time array.
        *params: Model parameters:
            - 3 params: [a1, tau1, offset]
            - 5 params: [a1, tau1, a2, tau2, offset]
            - 7 params: [a1, tau1, a2, tau2, a3, tau3, offset]

    Returns:
        Model values at each time point.
    """
    t = np.asarray(t).flatten()
    n_params = len(params)

    if n_params == 3:
        a1, tau1, offset = params
        return a1 * np.exp(-t / tau1) + offset

    elif n_params == 5:
        a1, tau1, a2, tau2, offset = params
        return a1 * np.exp(-t / tau1) + a2 * np.exp(-t / tau2) + offset

    elif n_params == 7:
        a1, tau1, a2, tau2, a3, tau3, offset = params
        return (
            a1 * np.exp(-t / tau1)
            + a2 * np.exp(-t / tau2)
            + a3 * np.exp(-t / tau3)
            + offset
        )

    else:
        raise ValueError(f"Expected 3, 5, or 7 parameters, got {n_params}")


def fit_exponential(
    t: NDArray[np.float64],
    y: NDArray[np.float64],
    n_exp: int = 1,
    initial_guess: NDArray[np.float64] | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], float]:
    """Fit exponential decay model to data.

    Args:
        t: Time array.
        y: Data values to fit.
        n_exp: Number of exponential terms (1, 2, or 3).
        initial_guess: Initial parameter guess. If None, uses defaults.

    Returns:
        Tuple of (parameters, fitted_curve, r_squared).

    Examples:
        >>> t = np.linspace(0, 10, 100)
        >>> y = 2.0 * np.exp(-t / 3.0) + 0.5
        >>> params, fit, rsq = fit_exponential(t, y, n_exp=1)
    """
    t = np.asarray(t).flatten()
    y = np.asarray(y).flatten()

    # Set up initial guess and bounds based on number of exponentials
    if n_exp == 1:
        if initial_guess is None:
            initial_guess = np.array([y[0], t[-1] / 5, y[-1]])
        bounds = (
            [0, 1e-10, -np.inf],
            [np.inf, np.inf, np.inf],
        )
    elif n_exp == 2:
        if initial_guess is None:
            initial_guess = np.array([y[0] * 0.5, t[-1] / 10, y[0] * 0.5, t[-1] / 2, y[-1]])
        bounds = (
            [0, 1e-10, 0, 1e-10, -np.inf],
            [np.inf, np.inf, np.inf, np.inf, np.inf],
        )
    elif n_exp == 3:
        if initial_guess is None:
            initial_guess = np.array([
                y[0] * 0.3, t[-1] / 20,
                y[0] * 0.3, t[-1] / 5,
                y[0] * 0.3, t[-1] / 2,
                y[-1]
            ])
        bounds = (
            [0, 1e-10, 0, 1e-10, 0, 1e-10, -np.inf],
            [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf],
        )
    else:
        raise ValueError(f"n_exp must be 1, 2, or 3, got {n_exp}")

    # Perform fit
    popt, _ = curve_fit(
        exponential_model,
        t,
        y,
        p0=initial_guess,
        bounds=bounds,
        maxfev=10000,
    )

    # Generate fitted curve
    y_fit = exponential_model(t, *popt)

    # Calculate R-squared
    ss_res = np.sum((y - y_fit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    return popt, y_fit, r_squared


def fit_spectral_exponential(
    t: NDArray[np.float64],
    y: NDArray[np.float64],
    initial_params: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], float]:
    """Fit exponential model to spectral data with residual analysis.

    Args:
        t: Time/frequency array.
        y: Data values to fit.
        initial_params: Initial parameter guess (3, 5, or 7 elements).

    Returns:
        Tuple of (parameters, fitted_curve, residuals, r_squared).
    """
    t = np.asarray(t).flatten()
    y = np.asarray(y).flatten()
    n_params = len(initial_params)

    # Determine number of exponentials
    if n_params == 3:
        n_exp = 1
    elif n_params == 5:
        n_exp = 2
    elif n_params == 7:
        n_exp = 3
    else:
        raise ValueError(f"Expected 3, 5, or 7 parameters, got {n_params}")

    # Fit
    params, y_fit, r_squared = fit_exponential(t, y, n_exp, initial_params)

    # Calculate relative residuals
    residuals = (y_fit - y) / np.mean(y)

    return params, y_fit, residuals, r_squared
