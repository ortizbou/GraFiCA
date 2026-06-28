"""Tests for grafica.community."""

import numpy as np
import pytest
from grafica.community import community_detection


class TestCommunityDetection:
    def _make_block_matrix(self, rng):
        """Helper: create a similarity matrix with clear block structure."""
        n = 30
        K = 3
        labels = np.repeat([1, 2, 3], 10)
        # Build a similarity matrix where within-cluster entries are low
        # (since we take K smallest eigenvectors, we want clusters to
        # correspond to small eigenvalues)
        W = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                if labels[i] == labels[j]:
                    W[i, j] = -1.0 + 0.1 * rng.standard_normal()
                else:
                    W[i, j] = 1.0 + 0.1 * rng.standard_normal()
                W[j, i] = W[i, j]
        return W, K, labels

    def test_returns_four_values(self, rng):
        """Should return (labels, nmi, ari, cost)."""
        W, K, gt = self._make_block_matrix(rng)
        result = community_detection(W, K, ground_truth=gt, kmeans_replicates=10, rng=rng)
        assert len(result) == 4

    def test_labels_shape(self, rng):
        """Labels should have length N."""
        W, K, gt = self._make_block_matrix(rng)
        labels, _, _, _ = community_detection(W, K, kmeans_replicates=10, rng=rng)
        assert len(labels) == W.shape[0]

    def test_labels_range(self, rng):
        """Labels should be in {1, ..., K}."""
        W, K, gt = self._make_block_matrix(rng)
        labels, _, _, _ = community_detection(W, K, kmeans_replicates=10, rng=rng)
        assert set(labels).issubset(set(range(1, K + 1)))

    def test_nmi_ari_nan_without_gt(self, rng):
        """NMI and ARI should be NaN when no ground truth is given."""
        W, K, _ = self._make_block_matrix(rng)
        _, nmi, ari, _ = community_detection(W, K, kmeans_replicates=10, rng=rng)
        assert np.isnan(nmi)
        assert np.isnan(ari)

    def test_nmi_ari_computed_with_gt(self, rng):
        """NMI and ARI should be finite floats in [0, 1] with ground truth."""
        W, K, gt = self._make_block_matrix(rng)
        _, nmi, ari, _ = community_detection(
            W, K, ground_truth=gt, kmeans_replicates=10, rng=rng
        )
        assert np.isfinite(nmi)
        assert 0.0 <= nmi <= 1.0
        assert np.isfinite(ari)

    def test_cost_is_scalar(self, rng):
        """Cost should be a scalar float."""
        W, K, _ = self._make_block_matrix(rng)
        _, _, _, cost = community_detection(W, K, kmeans_replicates=10, rng=rng)
        assert isinstance(cost, float)

    def test_validation_non_square(self):
        """Should raise ValueError for non-square W."""
        with pytest.raises(ValueError, match="square"):
            community_detection(np.ones((3, 4)), 2)

    def test_validation_bad_K(self):
        """Should raise ValueError for K < 1."""
        with pytest.raises(ValueError, match="positive"):
            community_detection(np.eye(3), 0)
