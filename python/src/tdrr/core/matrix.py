"""Matrix utility functions for spectral analysis.

Includes Toeplitz matrix generation and other linear algebra utilities
used in higher-order spectral analysis.
"""

import numpy as np
from numpy.typing import NDArray


def toeplitz_matrix(
    column: NDArray[np.float64],
    row: NDArray[np.float64] | None = None,
) -> NDArray[np.float64]:
    """Generate a Toeplitz matrix from column and row vectors.

    A Toeplitz matrix has constant diagonals. This is used in
    bispectrum and cumulant calculations.

    Args:
        column: First column of the matrix.
        row: First row of the matrix. If None, uses column.

    Returns:
        Toeplitz matrix of shape (len(column), len(row)).
    """
    column = np.atleast_1d(column).flatten()

    if row is None:
        row = column.copy()
    else:
        row = np.atleast_1d(row).flatten()

    n_col = len(column)
    n_row = len(row)

    # Build vector [row(end:-1:2), column]
    vec = np.concatenate([row[-1:0:-1], column])

    # Generate Toeplitz matrix
    mat = np.zeros((n_col, n_row), dtype=column.dtype)
    for k in range(n_row):
        mat[:, k] = vec[n_row - 1 - k : n_row - 1 - k + n_col]

    return mat


def symmetric_toeplitz(vector: NDArray[np.float64]) -> NDArray[np.float64]:
    """Generate a symmetric Toeplitz matrix.

    Args:
        vector: First row/column of the symmetric matrix.

    Returns:
        Symmetric Toeplitz matrix.
    """
    return toeplitz_matrix(vector, vector)


def circulant_matrix(vector: NDArray[np.float64]) -> NDArray[np.float64]:
    """Generate a circulant matrix from a vector.

    A circulant matrix is a special Toeplitz matrix where each row
    is a cyclic shift of the previous row.

    Args:
        vector: First row of the circulant matrix.

    Returns:
        Circulant matrix of shape (n, n) where n = len(vector).
    """
    n = len(vector)
    mat = np.zeros((n, n), dtype=vector.dtype)

    for i in range(n):
        mat[i, :] = np.roll(vector, i)

    return mat
