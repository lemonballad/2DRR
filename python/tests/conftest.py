"""Pytest configuration and fixtures."""

import numpy as np
import pytest


@pytest.fixture
def sample_signal() -> np.ndarray:
    """Generate a simple test signal."""
    t = np.linspace(0, 10, 1000)
    return np.sin(2 * np.pi * t) + 0.5 * np.sin(4 * np.pi * t)


@pytest.fixture
def frequency_axis() -> np.ndarray:
    """Generate a standard frequency axis."""
    return np.linspace(0, 1000, 1001)
