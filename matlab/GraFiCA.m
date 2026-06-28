function results = GraFiCA(A, f, K, T, opts)
%  GraFiCA: Learning optimal graph filters for clustering of attributed graphs.
%
%  Learns FIR graph filter coefficients that jointly optimize graph structure
%  and node attributes for community detection. Iterates between:
%    1. Filter learning: minimize within-cluster signal variation
%    2. Community detection: spectral clustering on filtered attributes
%
%  Input:
%         A:    adjacency matrix (N x N)
%         f:    node attribute matrix (N x d), d attributes per node
%         K:    number of clusters
%         T:    filter order
%
%  Optional name-value arguments:
%         alpha:             weight for topology vs attributes, default 0.5
%         gamma:             weight for between-cluster term, default 1.0
%         ground_truth:      ground truth labels (N x 1) or [] if unavailable
%         max_iter:          maximum iterations, default 20
%         conv_thresh:       convergence threshold on filter coefficients, default 1e-3
%         kmeans_replicates: number of k-means replicates, default 500
%
%  Output: struct with fields:
%         h:          learned filter coefficients (T x 1)
%         H:          filter matrix (N x N)
%         f_tilde:    filtered signal (N x d)
%         idx:        cluster labels (N x 1)
%         W:          weighted adjacency matrix (N x N)
%         n_iter:     number of iterations performed
%         NMI:        Normalized Mutual Information (NaN if no ground truth)
%         ARI:        Adjusted Rand Index (NaN if no ground truth)
%
%  Author: Meiby Ortiz-Bouza
%  Address: Michigan State University, ECE
%  email: ortizbou@msu.edu
%
%  Reference:
%  [1] M. Ortiz-Bouza and S. Aviyente, "Learning Optimal Graph Filters for
%      Clustering of Attributed Graphs," IEEE Trans. Signal Inform. Process.
%      over Networks, vol. 11, pp. 520-534, 2025.

arguments
    A
    f
    K    (1,1) double
    T    (1,1) double
    opts.alpha             (1,1) double = 0.5
    opts.gamma             (1,1) double = 1.0
    opts.ground_truth                   = []
    opts.max_iter          (1,1) double = 20
    opts.conv_thresh       (1,1) double = 1e-3
    opts.kmeans_replicates (1,1) double = 500
end

%% Input validation
assert(ismatrix(A) && size(A,1) == size(A,2), ...
    'GraFiCA:invalidInput', 'A must be a square matrix.');
N = size(A, 1);
assert(size(f, 1) == N, ...
    'GraFiCA:invalidInput', 'f must have N rows matching A.');
assert(K > 0 && K == round(K), ...
    'GraFiCA:invalidInput', 'K must be a positive integer.');
assert(T > 0 && T == round(T), ...
    'GraFiCA:invalidInput', 'T must be a positive integer.');

%% Preprocessing: eigendecomposition and spectral projections
An = normadj(A);
[U, D_mat] = eig(full(A));
D = diag(D_mat);  % eigenvalues as vector

% Compute spectral projections zn (cell array of T x 1 vectors)
for i = 1:N
    zn{i} = zeros(T, 1);
    for t = 1:T
        zn{i}(t) = U(i, :) * (D .^ (t - 1));
    end
end

%% Initial clustering to get B and C
[idx_init, ~, ~, ~] = CommDetec(An, opts.ground_truth, K, ...
    'kmeans_replicates', opts.kmeans_replicates);
[B, C] = commMatricesTxT(idx_init, K, T, An, zn);

%% Iterative optimization
h_old = zeros(T, 1);

for n = 1:opts.max_iter
    %% Filter learning step
    S = B - opts.gamma * C;
    [u, v] = eig(S);
    [~, ind] = sort(diag(v), 'ascend');
    u = u(:, ind);
    hnew = u(:, 1);

    %% Construct filter matrix
    H = zeros(N);
    for t = 1:T
        H = H + hnew(t) * diag(D .^ (t - 1));
    end

    %% Apply filter to signal
    f_tilde = U * H * U' * f;

    %% Community detection step
    % Pairwise attribute distances
    Af1 = zeros(N);
    for i = 1:N
        Af1(:, i) = vecnorm((f_tilde - f_tilde(i, :)), 2, 2);
    end
    Af = normadj(Af1);

    % Weighted adjacency
    W = Af - opts.alpha * An;

    % Spectral clustering
    [idxn, NMI_val, ARI_val, ~] = CommDetec(W, opts.ground_truth, K, ...
        'kmeans_replicates', opts.kmeans_replicates);

    %% Convergence check
    if norm(hnew - h_old) < opts.conv_thresh
        break
    end

    %% Update B and C for next iteration
    [B, C] = commMatricesTxT(idxn, K, T, An, zn);
    h_old = hnew;
end

%% Pack results
results.h = hnew;
results.H = U * diag(diag(H)) * U';
results.f_tilde = f_tilde;
results.idx = idxn;
results.W = W;
results.n_iter = n;
results.NMI = NMI_val;
results.ARI = ARI_val;

end
