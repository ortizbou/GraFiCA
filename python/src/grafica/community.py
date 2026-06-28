"""Spectral clustering for community detection.

Eigendecomposition of a similarity matrix followed by k-means clustering.

Reference:
    M. Ortiz-Bouza and S. Aviyente, "Learning Optimal Graph Filters for
    Clustering of Attributed Graphs," IEEE Trans. Signal Inform. Process.
    over Networks, vol. 11, pp. 520-534, 2025.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import eigh
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def community_detection(
    W: np.ndarray,
    K: int,
    ground_truth: np.ndarray | None = None,
    kmeans_replicates: int = 500,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, float, float, float]:
    """Spectral clustering on a similarity matrix.

    Computes the K smallest eigenvectors of W, then applies k-means.

    Parameters
    ----------
    W : np.ndarray
        Similarity matrix (N x N), symmetric.
    K : int
        Number of clusters.
    ground_truth : np.ndarray, optional
        Ground truth labels (length N). When provided, NMI and ARI are
        computed against the detected communities.
    kmeans_replicates : int
        Number of k-means replicates. Default 500.
    rng : np.random.Generator, optional
        Random number generator for k-means reproducibility.

    Returns
    -------
    tuple[np.ndarray, float, float, float]
        (labels, nmi, ari, cost) where labels are 1-indexed cluster
        assignments (length N), and nmi/ari are NaN if ground_truth
        is None.
    """
    W = np.asarray(W, dtype=float)
    if W.ndim != 2 or W.shape[0] != W.shape[1]:
        raise ValueError("W must be a square matrix.")
    if K < 1:
        raise ValueError("K must be a positive integer.")

    N = W.shape[0]

    # K smallest eigenvectors
    eigenvalues, eigenvectors = eigh(W, subset_by_index=[0, K - 1])
    V = eigenvectors  # (N, K)

    cost = float(np.trace(V.T @ W @ V))

    # k-means clustering
    seed = int(rng.integers(2**31)) if rng is not None else None
    km = KMeans(n_clusters=K, n_init=kmeans_replicates, random_state=seed)
    labels_0indexed = km.fit_predict(V)
    labels = labels_0indexed + 1  # 1-indexed to match MATLAB convention

    # Evaluation metrics
    if ground_truth is not None:
        gt = np.asarray(ground_truth).ravel()
        nmi = float(normalized_mutual_info_score(gt, labels))
        ari = float(adjusted_rand_score(gt, labels))
    else:
        nmi = float("nan")
        ari = float("nan")

    return labels, nmi, ari, cost
