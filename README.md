# GraFiCA

Learning optimal graph filters for clustering of attributed graphs.

**Authors:** Meiby Ortiz-Bouza and Selin Aviyente\
Department of Electrical and Computer Engineering, Michigan State University, MI\
https://ieeexplore.ieee.org/abstract/document/11021005

## Overview

GraFiCA is a graph signal processing approach that learns the parameters of FIR graph filters optimized for clustering attributed graphs. The method formulates a two-step iterative optimization: (1) learn interpretable graph filters that maximize cluster separation, and (2) apply spectral clustering on the filtered attribute distances.

![framewokB2](https://github.com/user-attachments/assets/f86b9441-21ba-4250-a23b-92f983a23bce)

## Repository Structure

```
GraFiCA/
├── matlab/                        # MATLAB implementation
│   ├── GraFiCA.m                    # Main algorithm (function)
│   ├── CommDetec.m                  # Spectral clustering
│   ├── commMatricesTxT.m           # Within/between-cluster matrices
│   ├── examples/
│   │   ├── demo_synthetic.m         # Example with synthetic data + ground truth
│   │   └── demo_real.m              # Example without ground truth
│   └── helpers/
│       ├── getNMI.m                 # Normalized Mutual Information
│       ├── normadj.m                # Normalized adjacency matrix
│       ├── rand_index.m             # Adjusted Rand Index
│       ├── rnorm.m                  # Row normalization
│       └── vec.m                    # Vectorization utility
├── python/                        # Python implementation
│   ├── pyproject.toml               # Package metadata and dependencies
│   ├── src/grafica/
│   │   ├── filtering.py             # Core GraFiCA algorithm
│   │   ├── community.py             # Spectral clustering
│   │   ├── cluster_matrices.py      # Within/between-cluster matrices
│   │   └── utils.py                 # Graph utilities
│   ├── tests/                       # pytest test suite
│   └── examples/
│       └── tutorial.py              # End-to-end tutorial script
└── data/                          # Shared example data
```

## Installation

### Python

```bash
cd python
pip install -e .          # basic install
pip install -e ".[dev]"   # with dev tools (pytest, ruff)
```

Requires Python >= 3.9, NumPy >= 1.22, SciPy >= 1.8, scikit-learn >= 1.0.

### MATLAB

No installation needed. Add `matlab/` and `matlab/helpers/` to your MATLAB path.

## Quick Start

### Python

```python
import numpy as np
from grafica import grafica

# A: adjacency matrix (N x N)
# f: node attributes (N x d)
# K: number of clusters
# T: filter order

result = grafica(A, f, K=3, T=4, alpha=0.5, gamma=1.0)

result.labels     # cluster assignments
result.h          # learned filter coefficients
result.f_tilde    # filtered signal
result.n_iter     # iterations to convergence

# With ground truth for evaluation
result = grafica(A, f, K=3, T=4, ground_truth=labels)
result.nmi        # Normalized Mutual Information
result.ari        # Adjusted Rand Index
```

See `python/examples/tutorial.py` for a complete walkthrough.

### MATLAB

```matlab
% Run GraFiCA
results = GraFiCA(A, f, K, T, 'alpha', 0.5, 'gamma', 1.0);

% With ground truth
results = GraFiCA(A, f, K, T, 'ground_truth', GT);

% Access results
results.idx       % cluster labels
results.h         % filter coefficients
results.NMI       % Normalized Mutual Information
results.ARI       % Adjusted Rand Index
```

See `matlab/examples/demo_synthetic.m` and `matlab/examples/demo_real.m` for complete walkthroughs.

## Testing

```bash
cd python
pytest -v    # 31 tests
```

## Citation

```bibtex
@article{ortiz2025learning,
  title={Learning Optimal Graph Filters for Clustering of Attributed Graphs},
  author={Ortiz-Bouza, Meiby and Aviyente, Selin},
  journal={IEEE Transactions on Signal and Information Processing over Networks},
  volume={11},
  pages={520--534},
  year={2025},
  publisher={IEEE}
}
```
