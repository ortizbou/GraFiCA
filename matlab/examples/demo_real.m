%% GraFiCA Demo: Real Attributed Graph
% This script demonstrates how to run GraFiCA on real data where
% ground truth labels are not available.
%
% Reference:
%   M. Ortiz-Bouza and S. Aviyente, "Learning Optimal Graph Filters for
%   Clustering of Attributed Graphs," IEEE Trans. Signal Inform. Process.
%   over Networks, vol. 11, pp. 520-534, 2025.

%% Add paths
addpath('..');
addpath('../helpers');

%% Step 0: Load your data
% Replace this section with your own data loading.
% A: adjacency matrix (N x N)
% f: node attributes (N x d)
% K: number of clusters
% T: filter order
%
% Example:
%   load('my_graph.mat', 'A', 'f');
%   K = 5;
%   T = 4;

% --- Placeholder: generate example data ---
n = 60;
K = 3;
T = 4;
GT = [ones(20,1); 2*ones(20,1); 3*ones(20,1)];
A = generate_planted_partition(GT, 0.4, 0.05);
rng(42);
f = zeros(n, 2);
for k = 1:K
    mask = (GT == k);
    f(mask, :) = k * 2 + 0.5 * randn(sum(mask), 2);
end
% --- End placeholder ---

%% Step 1: Run GraFiCA (no ground truth)
fprintf('Running GraFiCA...\n');

results = GraFiCA(A, f, K, T, ...
    'alpha', 0.5, ...
    'gamma', 1.0, ...
    'max_iter', 20, ...
    'kmeans_replicates', 100);

%% Step 2: Examine results
fprintf('\nResults:\n');
fprintf('  Converged in %d iterations\n', results.n_iter);
fprintf('  Filter coefficients: [%s]\n', num2str(results.h', '%.4f '));

fprintf('\n  Cluster sizes:\n');
for k = 1:K
    fprintf('    Cluster %d: %d nodes\n', k, sum(results.idx == k));
end

%% Helper function
function A = generate_planted_partition(labels, p_in, p_out)
    n = length(labels);
    A = zeros(n);
    for i = 1:n
        for j = i+1:n
            if labels(i) == labels(j)
                A(i,j) = rand() < p_in;
            else
                A(i,j) = rand() < p_out;
            end
        end
    end
    A = A + A';
end
