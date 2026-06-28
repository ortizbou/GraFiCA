"""GraFiCA Tutorial: Python implementation.

This script demonstrates the full GraFiCA pipeline on a synthetic
attributed graph with known community structure.

Usage:
    cd python/
    pip install -e .
    python examples/tutorial.py
"""

import numpy as np
from grafica import grafica


def generate_planted_partition(labels, p_in=0.4, p_out=0.05, rng=None):
    """Generate a symmetric adjacency matrix from a planted partition model."""
    if rng is None:
        rng = np.random.default_rng()
    n = len(labels)
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            p = p_in if labels[i] == labels[j] else p_out
            if rng.random() < p:
                A[i, j] = 1.0
                A[j, i] = 1.0
    return A


def main():
    rng = np.random.default_rng(42)

    # -------------------------------------------------------------------------
    # Step 0: Generate synthetic attributed graph
    # -------------------------------------------------------------------------
    print("=" * 60)
    print("GraFiCA Tutorial")
    print("=" * 60)

    n = 60
    K = 3
    T = 4
    labels_true = np.repeat([1, 2, 3], 20)

    A = generate_planted_partition(labels_true, p_in=0.4, p_out=0.05, rng=rng)

    # Attributes correlated with clusters
    f = np.zeros((n, 2))
    for k in range(1, K + 1):
        mask = labels_true == k
        f[mask, :] = rng.normal(loc=k * 2.0, scale=0.5, size=(mask.sum(), 2))

    print(f"\nGenerated attributed graph:")
    print(f"  Nodes: {n}, Attributes: {f.shape[1]}")
    print(f"  Clusters: {K}, Filter order: {T}")

    # -------------------------------------------------------------------------
    # Step 1: Run GraFiCA with ground truth
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("Step 1: Running GraFiCA with ground truth")
    print("-" * 60)

    result = grafica(
        A, f, K, T,
        alpha=0.5, gamma=1.0,
        ground_truth=labels_true,
        max_iter=20,
        kmeans_replicates=50,
        rng=np.random.default_rng(42),
    )

    print(f"\n  Converged in {result.n_iter} iterations")
    print(f"  NMI: {result.nmi:.4f}")
    print(f"  ARI: {result.ari:.4f}")
    print(f"  Filter coefficients: {np.array2string(result.h, precision=4)}")

    # -------------------------------------------------------------------------
    # Step 2: Run without ground truth
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("Step 2: Running GraFiCA without ground truth")
    print("-" * 60)

    result_real = grafica(
        A, f, K, T,
        alpha=0.5, gamma=1.0,
        max_iter=20,
        kmeans_replicates=50,
        rng=np.random.default_rng(42),
    )

    print(f"\n  Converged in {result_real.n_iter} iterations")
    print(f"  Cluster sizes: ", end="")
    for k in range(1, K + 1):
        print(f"{np.sum(result_real.labels == k)} ", end="")
    print()

    # -------------------------------------------------------------------------
    # Step 3: Experiment with parameters
    # -------------------------------------------------------------------------
    print("\n" + "-" * 60)
    print("Step 3: Parameter sensitivity (alpha)")
    print("-" * 60)

    for alpha in [0.1, 0.5, 1.0, 2.0]:
        r = grafica(
            A, f, K, T,
            alpha=alpha, gamma=1.0,
            ground_truth=labels_true,
            max_iter=20,
            kmeans_replicates=50,
            rng=np.random.default_rng(42),
        )
        print(f"  alpha={alpha:.1f}  NMI={r.nmi:.4f}  ARI={r.ari:.4f}  iters={r.n_iter}")

    print("\n" + "=" * 60)
    print("Tutorial complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
