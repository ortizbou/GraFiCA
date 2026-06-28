%% GraFiCA Demo: Synthetic Attributed Graph
% This script demonstrates how to run GraFiCA on a synthetic attributed
% graph with known community structure.
%
% Reference:
%   M. Ortiz-Bouza and S. Aviyente, "Learning Optimal Graph Filters for
%   Clustering of Attributed Graphs," IEEE Trans. Signal Inform. Process.
%   over Networks, vol. 11, pp. 520-534, 2025.

%% Add paths
addpath('..');
addpath('../helpers');

%% Step 0: Generate a synthetic attributed graph
n = 60;           % number of nodes
K = 3;            % number of clusters
T = 4;            % filter order

% Ground truth labels
GT = [ones(20,1); 2*ones(20,1); 3*ones(20,1)];

% Generate adjacency matrix (planted partition model)
A = generate_planted_partition(GT, 0.4, 0.05);

% Generate node attributes correlated with cluster membership
rng(42);
f = zeros(n, 2);
for k = 1:K
    mask = (GT == k);
    f(mask, :) = k * 2 + 0.5 * randn(sum(mask), 2);
end

fprintf('Generated attributed graph:\n');
fprintf('  Nodes: %d, Attributes: %d\n', n, size(f, 2));
fprintf('  Clusters: %d, Filter order: %d\n', K, T);

%% Step 1: Run GraFiCA with ground truth
fprintf('\nRunning GraFiCA...\n');

results = GraFiCA(A, f, K, T, ...
    'alpha', 0.5, ...
    'gamma', 1.0, ...
    'ground_truth', GT, ...
    'kmeans_replicates', 100);

fprintf('\nResults:\n');
fprintf('  Converged in %d iterations\n', results.n_iter);
fprintf('  NMI: %.4f\n', results.NMI);
fprintf('  ARI: %.4f\n', results.ARI);
fprintf('  Filter coefficients: [%s]\n', num2str(results.h', '%.4f '));

%% Step 2: Run without ground truth
fprintf('\nRunning GraFiCA without ground truth...\n');

results_real = GraFiCA(A, f, K, T, ...
    'alpha', 0.5, ...
    'gamma', 1.0, ...
    'kmeans_replicates', 100);

fprintf('  Converged in %d iterations\n', results_real.n_iter);
fprintf('  Cluster sizes: ');
for k = 1:K
    fprintf('%d ', sum(results_real.idx == k));
end
fprintf('\n');

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
