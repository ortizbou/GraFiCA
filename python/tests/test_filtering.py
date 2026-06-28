"""Tests for grafica.filtering."""

import numpy as np
import pytest
from grafica.filtering import GraFiCAResult, grafica


class TestGraFiCA:
    def _make_attributed_graph(self, rng):
        """Helper: create attributed graph with clear cluster structure."""
        n = 30
        K = 3
        labels = np.repeat([1, 2, 3], 10)

        # Adjacency: planted partition
        A = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                p = 0.5 if labels[i] == labels[j] else 0.03
                if rng.random() < p:
                    A[i, j] = 1.0
                    A[j, i] = 1.0

        # Attributes correlated with clusters
        f = np.zeros((n, 2))
        for k in range(1, K + 1):
            mask = labels == k
            f[mask, :] = rng.normal(loc=k * 2.0, scale=0.5, size=(mask.sum(), 2))

        return A, f, K, labels

    def test_returns_result_object(self, rng):
        """Should return a GraFiCAResult."""
        A, f, K, gt = self._make_attributed_graph(rng)
        result = grafica(A, f, K, T=3, max_iter=3, kmeans_replicates=5, rng=rng)
        assert isinstance(result, GraFiCAResult)

    def test_result_shapes(self, rng):
        """Result arrays should have expected shapes."""
        A, f, K, gt = self._make_attributed_graph(rng)
        n = A.shape[0]
        T = 3
        result = grafica(A, f, K, T=T, max_iter=3, kmeans_replicates=5, rng=rng)
        assert result.h.shape == (T,)
        assert result.H.shape == (n, n)
        assert result.f_tilde.shape == f.shape
        assert result.labels.shape == (n,)
        assert result.W.shape == (n, n)

    def test_labels_range(self, rng):
        """Labels should be in {1, ..., K}."""
        A, f, K, gt = self._make_attributed_graph(rng)
        result = grafica(A, f, K, T=3, max_iter=3, kmeans_replicates=5, rng=rng)
        assert set(result.labels).issubset(set(range(1, K + 1)))

    def test_n_iter_within_bounds(self, rng):
        """n_iter should be between 1 and max_iter."""
        A, f, K, gt = self._make_attributed_graph(rng)
        max_iter = 5
        result = grafica(A, f, K, T=3, max_iter=max_iter, kmeans_replicates=5, rng=rng)
        assert 1 <= result.n_iter <= max_iter

    def test_with_ground_truth(self, rng):
        """With ground truth, NMI and ARI should be finite."""
        A, f, K, gt = self._make_attributed_graph(rng)
        result = grafica(
            A, f, K, T=3, ground_truth=gt,
            max_iter=3, kmeans_replicates=5, rng=rng,
        )
        assert np.isfinite(result.nmi)
        assert np.isfinite(result.ari)

    def test_without_ground_truth(self, rng):
        """Without ground truth, NMI and ARI should be NaN."""
        A, f, K, _ = self._make_attributed_graph(rng)
        result = grafica(A, f, K, T=3, max_iter=3, kmeans_replicates=5, rng=rng)
        assert np.isnan(result.nmi)
        assert np.isnan(result.ari)

    def test_reproducible(self, rng):
        """Same rng seed should produce same result."""
        A, f, K, _ = self._make_attributed_graph(rng)
        r1 = grafica(A, f, K, T=3, max_iter=3, kmeans_replicates=5,
                      rng=np.random.default_rng(99))
        r2 = grafica(A, f, K, T=3, max_iter=3, kmeans_replicates=5,
                      rng=np.random.default_rng(99))
        np.testing.assert_array_equal(r1.h, r2.h)
        np.testing.assert_array_equal(r1.labels, r2.labels)

    def test_validation_non_square_A(self):
        """Should raise ValueError for non-square A."""
        with pytest.raises(ValueError, match="square"):
            grafica(np.ones((3, 4)), np.ones((3, 2)), K=2, T=2)

    def test_validation_f_mismatch(self):
        """Should raise ValueError when f rows don't match A."""
        with pytest.raises(ValueError, match="N rows"):
            grafica(np.eye(5), np.ones((3, 2)), K=2, T=2)

    def test_validation_bad_K(self):
        """Should raise ValueError for K < 1."""
        with pytest.raises(ValueError, match="positive"):
            grafica(np.eye(5), np.ones((5, 2)), K=0, T=2)
