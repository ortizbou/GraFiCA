"""GraFiCA: Learning optimal graph filters for clustering of attributed graphs.

Learns FIR graph filter coefficients optimized for community detection on
attributed graphs, combining graph structure with node attribute similarity.

Reference:
    M. Ortiz-Bouza and S. Aviyente, "Learning Optimal Graph Filters for
    Clustering of Attributed Graphs," IEEE Trans. Signal Inform. Process.
    over Networks, vol. 11, pp. 520-534, 2025.
"""

from grafica.cluster_matrices import cluster_matrices
from grafica.community import community_detection
from grafica.filtering import GraFiCAResult, grafica
from grafica.utils import normalized_adjacency, row_normalize

__all__ = [
    "cluster_matrices",
    "community_detection",
    "grafica",
    "GraFiCAResult",
    "normalized_adjacency",
    "row_normalize",
]

__version__ = "0.1.0"
