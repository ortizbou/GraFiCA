"""Core GraFiCA algorithm: graph filter learning for attributed graph clustering.

Iterates between:
  1. Filter learning — minimize within-cluster signal variation (eigenvalue problem)
  2. Community detection — spectral clustering on filtered attribute distances

Reference:
    M. Ortiz-Bouza and S. Aviyente, "Learning Optimal Graph Filters for
    Clustering of Attributed Graphs," IEEE Trans. Signal Inform. Process.
    over Networks, vol. 11, pp. 520-534, 2025.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial.distance import cdist

from grafica.cluster_matrices import cluster_matrices
from grafica.community import community_detection
from grafica.utils import normalized_adjacency


@dataclass
class GraFiCAResult:
    """Container for GraFiCA results."""

    h: np.ndarray = field(default_factory=lambda: np.array([]))
    H: np.ndarray = field(default_factory=lambda: np.array([]))
    f_tilde: np.ndarray = field(default_factory=lambda: np.array([]))
    labels: np.ndarray = field(default_factory=lambda: np.array([]))
    W: np.ndarray = field(default_factory=lambda: np.array([]))
    n_iter: int = 0
    nmi: float = float("nan")
    ari: float = float("nan")


def grafica(
    A: np.ndarray,
    f: np.ndarray,
    K: int,
    T: int,
    alpha: float = 0.5,
    gamma: float = 1.0,
    ground_truth: np.ndarray | None = None,
    max_iter: int = 20,
    conv_thresh: float = 1e-3,
    kmeans_replicates: int = 500,
    rng: np.random.Generator | None = None,
) -> GraFiCAResult:
    """Learn optimal graph filters for clustering of attributed graphs.

    Parameters
    ----------
    A : np.ndarray
        Adjacency matrix (N x N).
    f : np.ndarray
        Node attribute matrix (N x d).
    K : int
        Number of clusters.
    T : int
        Filter order.
    alpha : float
        Weight for topology vs attributes. Default 0.5.
    gamma : float
        Weight for between-cluster term. Default 1.0.
    ground_truth : np.ndarray, optional
        Ground truth labels (length N). If provided, NMI/ARI are computed.
    max_iter : int
        Maximum iterations. Default 20.
    conv_thresh : float
        Convergence threshold on filter coefficients. Default 1e-3.
    kmeans_replicates : int
        Number of k-means replicates. Default 500.
    rng : np.random.Generator, optional
        Random number generator for reproducibility.

    Returns
    -------
    GraFiCAResult
        Result object with filter coefficients, filtered signal, labels, etc.
    """
    A = np.asarray(A, dtype=float)
    f = np.asarray(f, dtype=float)

    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("A must be a square matrix.")
    N = A.shape[0]
    if f.shape[0] != N:
        raise ValueError("f must have N rows matching A.")
    if K < 1:
        raise ValueError("K must be a positive integer.")
    if T < 1:
        raise ValueError("T must be a positive integer.")

    if f.ndim == 1:
        f = f.reshape(-1, 1)

    # Preprocessing: eigendecomposition and spectral projections
    An = normalized_adjacency(A)
    eigenvalues, U = np.linalg.eigh(A)

    # Spectral projections: zn[i, t] = U[i, :] @ (eigenvalues ** (t-1))
    # Shape: (N, T)
    zn = np.column_stack([U @ (eigenvalues ** t) for t in range(T)])

    # Initial clustering to bootstrap B and C
    labels_init, _, _, _ = community_detection(
        An, K, ground_truth=ground_truth,
        kmeans_replicates=kmeans_replicates, rng=rng,
    )
    B, C = cluster_matrices(labels_init, K, T, An, zn)

    # Iterative optimization
    h_old = np.zeros(T)
    result = GraFiCAResult()

    for n_iter in range(1, max_iter + 1):
        # Filter learning step
        S = B - gamma * C
        eigvals_S, eigvecs_S = np.linalg.eigh(S)
        h_new = eigvecs_S[:, 0]  # eigenvector for smallest eigenvalue

        # Construct filter matrix: H = sum_t h(t) * diag(eigenvalues^(t-1))
        H_diag = np.zeros(N)
        for t in range(T):
            H_diag += h_new[t] * (eigenvalues ** t)
        H = np.diag(H_diag)

        # Apply filter to signal
        f_tilde = U @ H @ U.T @ f

        # Community detection step
        # Pairwise attribute distances
        Af1 = cdist(f_tilde, f_tilde, metric='euclidean')
        Af = normalized_adjacency(Af1)

        # Weighted adjacency
        W = Af - alpha * An

        # Spectral clustering
        labels_n, nmi_val, ari_val, _ = community_detection(
            W, K, ground_truth=ground_truth,
            kmeans_replicates=kmeans_replicates, rng=rng,
        )

        # Convergence check
        if np.linalg.norm(h_new - h_old) < conv_thresh:
            result.h = h_new
            result.H = U @ H @ U.T
            result.f_tilde = f_tilde
            result.labels = labels_n
            result.W = W
            result.n_iter = n_iter
            result.nmi = nmi_val
            result.ari = ari_val
            break

        # Update B and C for next iteration
        B, C = cluster_matrices(labels_n, K, T, An, zn)
        h_old = h_new.copy()
    else:
        # Did not converge, return last iteration results
        result.h = h_new
        result.H = U @ H @ U.T
        result.f_tilde = f_tilde
        result.labels = labels_n
        result.W = W
        result.n_iter = max_iter
        result.nmi = nmi_val
        result.ari = ari_val

    return result
