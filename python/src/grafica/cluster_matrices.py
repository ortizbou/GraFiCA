"""Within-cluster and between-cluster matrix computation.

Builds the T x T matrices B and C used in the filter optimization
step of GraFiCA.

Reference:
    M. Ortiz-Bouza and S. Aviyente, "Learning Optimal Graph Filters for
    Clustering of Attributed Graphs," IEEE Trans. Signal Inform. Process.
    over Networks, vol. 11, pp. 520-534, 2025.
"""

from __future__ import annotations

import numpy as np


def cluster_matrices(
    labels: np.ndarray,
    K: int,
    T: int,
    A: np.ndarray,
    zn: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute within-cluster (B) and between-cluster (C) matrices.

    B captures within-cluster signal differences and C captures
    between-cluster signal differences, both weighted by cluster volume.

    Parameters
    ----------
    labels : np.ndarray
        Cluster labels (length N), integers in {1, ..., K}.
    K : int
        Number of clusters.
    T : int
        Filter order.
    A : np.ndarray
        Adjacency matrix (N x N), for computing cluster volumes.
    zn : np.ndarray
        Spectral projections, shape (N, T).

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (B, C) each of shape (T, T).
    """
    labels = np.asarray(labels).ravel()
    A = np.asarray(A, dtype=float)
    zn = np.asarray(zn, dtype=float)

    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("A must be a square matrix.")
    if zn.shape != (len(labels), T):
        raise ValueError(f"zn must have shape ({len(labels)}, {T}), got {zn.shape}.")

    B = np.zeros((T, T))
    C = np.zeros((T, T))

    for k in range(1, K + 1):
        members = np.where(labels == k)[0]
        non_members = np.where(labels != k)[0]
        vol_k = A[members, :].sum()

        # Within-cluster differences
        Bsum = np.zeros((T, T))
        for i in range(len(members)):
            for j in range(i, len(members)):
                diff = zn[members[i]] - zn[members[j]]
                Bsum += np.outer(diff, diff)
        B += Bsum / vol_k

        # Between-cluster differences
        Csum = np.zeros((T, T))
        for i in range(len(members)):
            for j in range(len(non_members)):
                diff = zn[members[i]] - zn[non_members[j]]
                Csum += np.outer(diff, diff)
        C += Csum / vol_k

    return B, C
