"""Tests for grafica.cluster_matrices."""

import numpy as np
import pytest
from grafica.cluster_matrices import cluster_matrices


class TestClusterMatrices:
    def _make_data(self, rng):
        """Helper: create a simple 2-cluster problem."""
        labels = np.array([1, 1, 1, 2, 2, 2])
        n = len(labels)
        T = 3
        K = 2
        # Simple adjacency: block diagonal
        A = np.zeros((n, n))
        A[:3, :3] = 1.0
        A[3:, 3:] = 1.0
        np.fill_diagonal(A, 0)
        # Random spectral projections
        zn = rng.standard_normal((n, T))
        return labels, K, T, A, zn

    def test_output_shapes(self, rng):
        """B and C should both be (T, T)."""
        labels, K, T, A, zn = self._make_data(rng)
        B, C = cluster_matrices(labels, K, T, A, zn)
        assert B.shape == (T, T)
        assert C.shape == (T, T)

    def test_symmetric(self, rng):
        """B and C should be symmetric (sums of outer products)."""
        labels, K, T, A, zn = self._make_data(rng)
        B, C = cluster_matrices(labels, K, T, A, zn)
        np.testing.assert_array_almost_equal(B, B.T)
        np.testing.assert_array_almost_equal(C, C.T)

    def test_positive_semidefinite(self, rng):
        """B and C should be positive semi-definite."""
        labels, K, T, A, zn = self._make_data(rng)
        B, C = cluster_matrices(labels, K, T, A, zn)
        eigvals_B = np.linalg.eigvalsh(B)
        eigvals_C = np.linalg.eigvalsh(C)
        assert np.all(eigvals_B >= -1e-10)
        assert np.all(eigvals_C >= -1e-10)

    def test_identical_signals_zero_B(self):
        """When all nodes in a cluster have the same zn, B should be zero."""
        labels = np.array([1, 1, 2, 2])
        T = 2
        K = 2
        A = np.ones((4, 4)) - np.eye(4)
        # Same spectral projection within each cluster
        zn = np.array([[1.0, 2.0], [1.0, 2.0], [3.0, 4.0], [3.0, 4.0]])
        B, C = cluster_matrices(labels, K, T, A, zn)
        np.testing.assert_array_almost_equal(B, np.zeros((T, T)))

    def test_validation_non_square_A(self, rng):
        """Should raise ValueError for non-square A."""
        labels = np.array([1, 1, 2])
        zn = rng.standard_normal((3, 2))
        with pytest.raises(ValueError, match="square"):
            cluster_matrices(labels, 2, 2, np.ones((3, 4)), zn)

    def test_validation_zn_shape_mismatch(self, rng):
        """Should raise ValueError if zn shape doesn't match."""
        labels = np.array([1, 1, 2])
        A = np.eye(3)
        zn = rng.standard_normal((3, 5))  # T=2 but zn has 5 columns
        with pytest.raises(ValueError, match="zn must have shape"):
            cluster_matrices(labels, 2, 2, A, zn)
