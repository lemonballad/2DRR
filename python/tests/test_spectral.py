"""Tests for core spectral analysis functions."""

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from tdrr.core.spectral import apodfun, lag_window, bispectrum


class TestApodfun:
    """Tests for apodfun function."""

    def test_flat_top_region(self) -> None:
        """Flat top region should be 1.0."""
        x = np.linspace(0, 100, 1001)
        result = apodfun(x, 50.0, 10.0, 10.0, 20.0)

        # Check flat top region (40-60)
        flat_mask = (x >= 40) & (x <= 60)
        assert_array_almost_equal(result[flat_mask], 1.0)

    def test_gaussian_decay(self) -> None:
        """Edges should decay as Gaussian."""
        x = np.linspace(0, 100, 1001)
        result = apodfun(x, 50.0, 10.0, 10.0, 20.0)

        # Check decay at edges
        assert result[0] < 0.01
        assert result[-1] < 0.01

    def test_asymmetric_decay(self) -> None:
        """Different g_low and g_high should create asymmetry."""
        x = np.linspace(0, 100, 1001)
        result = apodfun(x, 50.0, 5.0, 20.0, 10.0)

        # Low side should decay faster
        idx_35 = 350  # x = 35
        idx_65 = 650  # x = 65
        assert result[idx_35] < result[idx_65]

    def test_invert(self) -> None:
        """Inverted filter should be 1 - normal."""
        x = np.linspace(0, 100, 1001)
        normal = apodfun(x, 50.0, 10.0, 10.0, 20.0, invert=False)
        inverted = apodfun(x, 50.0, 10.0, 10.0, 20.0, invert=True)

        assert_array_almost_equal(normal + inverted, 1.0)

    def test_output_shape(self) -> None:
        """Output should have same shape as input."""
        x = np.linspace(0, 100, 500)
        result = apodfun(x, 50.0, 10.0, 10.0, 20.0)

        assert result.shape == x.shape


class TestLagWindow:
    """Tests for lag_window function."""

    def test_uniform_window(self) -> None:
        """Uniform window should be all ones."""
        w = lag_window(10, "uniform")
        assert_array_almost_equal(w, np.ones(10))

    def test_hamming_window(self) -> None:
        """Hamming window should start at 1.0."""
        w = lag_window(10, "hamming")
        assert w[0] == pytest.approx(1.0)
        assert w[-1] < w[0]

    def test_parzen_window(self) -> None:
        """Parzen window should be normalized."""
        w = lag_window(10, "parzen")
        assert w[0] == pytest.approx(1.0)
        assert all(w >= 0)
        assert all(w <= 1)

    def test_window_aliases(self) -> None:
        """Short aliases should work."""
        w_full = lag_window(10, "hamming")
        w_alias = lag_window(10, "h")
        assert_array_almost_equal(w_full, w_alias)

    def test_single_lag(self) -> None:
        """Lag of 1 should return [1.0]."""
        w = lag_window(1, "parzen")
        assert_array_almost_equal(w, [1.0])

    def test_output_length(self) -> None:
        """Output length should equal lag parameter."""
        for lag in [5, 10, 20, 100]:
            w = lag_window(lag, "uniform")
            assert len(w) == lag

    def test_unknown_window_raises(self) -> None:
        """Unknown window type should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown window"):
            lag_window(10, "invalid_window")

    def test_all_window_types(self) -> None:
        """All window types should produce valid output."""
        windows = ["uniform", "sasaki", "priestley", "parzen",
                   "hamming", "gaussian", "daniell"]
        for window in windows:
            w = lag_window(20, window)
            assert len(w) == 20
            assert not np.any(np.isnan(w))


class TestBispectrum:
    """Tests for bispectrum function."""

    def test_output_shapes(self) -> None:
        """Output shapes should match expected dimensions."""
        signal = np.random.randn(100, 1)
        bisp, freq, cum, lag = bispectrum(signal, max_lag=10)

        assert bisp.shape == (21, 21)
        assert freq.shape == (21,)
        assert cum.shape == (21, 21)
        assert lag.shape == (21,)

    def test_zero_lag(self) -> None:
        """Zero max_lag should return scalar-like outputs."""
        signal = np.random.randn(100, 1)
        bisp, freq, cum, lag = bispectrum(signal, max_lag=0)

        assert bisp.shape == (1, 1)
        assert freq.shape == (1,)

    def test_symmetric_signal(self) -> None:
        """Symmetric signal should produce real bispectrum."""
        # Create symmetric test signal
        t = np.linspace(0, 2 * np.pi, 100)
        signal = np.cos(t).reshape(-1, 1)

        bisp, _, _, _ = bispectrum(signal, max_lag=5)

        # Imaginary part should be small for symmetric signal
        assert np.max(np.abs(bisp.imag)) < np.max(np.abs(bisp.real)) * 0.1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
