function [idx, nmi, ari, cost] = CommDetec(W, GT, K, opts)
%  Spectral clustering on a similarity matrix.
%
%  Computes the K smallest eigenvectors of W, then applies k-means
%  clustering. When ground truth labels are provided, NMI and ARI
%  are computed; otherwise they are returned as NaN.
%
%  Input:
%         W:    similarity matrix (N x N)
%         GT:   ground truth labels (N x 1), or [] if unavailable
%         K:    number of clusters
%
%  Optional name-value arguments:
%         kmeans_replicates:  number of k-means replicates, default 500
%
%  Output:
%         idx:   cluster labels (N x 1)
%         nmi:   Normalized Mutual Information vs GT (NaN if GT is empty)
%         ari:   Adjusted Rand Index vs GT (NaN if GT is empty)
%         cost:  trace(V' * W * V) objective value
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
    W
    GT
    K                      (1,1) double
    opts.kmeans_replicates (1,1) double = 500
end

%% Input validation
assert(ismatrix(W) && size(W,1) == size(W,2), ...
    'CommDetec:invalidInput', 'W must be a square matrix.');
assert(isscalar(K) && K > 0 && K == round(K), ...
    'CommDetec:invalidInput', 'K must be a positive integer.');

%% Spectral decomposition
[Vnew, ~] = eigs(W, K, 'sa');
cost = trace(Vnew' * W * Vnew);

%% k-means clustering
[idx, ~] = kmeans(Vnew, K, 'emptyaction', 'singleton', ...
    'replicates', opts.kmeans_replicates);

%% Evaluation metrics
if ~isempty(GT)
    nmi = getNMI(idx, GT);
    ari = rand_index(idx, GT, 'adjusted');
else
    nmi = NaN;
    ari = NaN;
end

end
