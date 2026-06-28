function [B, C] = commMatricesTxT(idx, K, T, A, zn)
%  Compute within-cluster (B) and between-cluster (C) matrices.
%
%  Builds the T x T matrices B and C used in the filter optimization
%  step of GraFiCA. B captures within-cluster signal differences and
%  C captures between-cluster signal differences, both weighted by
%  cluster volume.
%
%  Input:
%         idx:  cluster labels (N x 1)
%         K:    number of clusters
%         T:    filter order
%         A:    adjacency matrix (N x N), for computing volumes
%         zn:   cell array of spectral projections (each T x 1)
%
%  Output:
%         B:    within-cluster matrix (T x T)
%         C:    between-cluster matrix (T x T)
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
    idx
    K    (1,1) double
    T    (1,1) double
    A
    zn   cell
end

%% Input validation
assert(isvector(idx), ...
    'commMatricesTxT:invalidInput', 'idx must be a vector.');
assert(ismatrix(A) && size(A,1) == size(A,2), ...
    'commMatricesTxT:invalidInput', 'A must be a square matrix.');

%% Compute B and C
B = zeros(T);
C = zeros(T);

for k = 1:K
    a = find(idx == k);
    volC = sum(sum(A(a, :)));
    b = find(idx ~= k);

    %% Within-cluster differences
    Bsum = zeros(T);
    for i = 1:length(a)
        for j = i:length(a)
            diff_ij = zn{a(i)} - zn{a(j)};
            Bsum = Bsum + diff_ij * diff_ij';
        end
    end
    B = B + (1 / volC) * Bsum;

    %% Between-cluster differences
    Csum = zeros(T);
    for i = 1:length(a)
        for j = 1:length(b)
            diff_ij = zn{a(i)} - zn{b(j)};
            Csum = Csum + diff_ij * diff_ij';
        end
    end
    C = C + (1 / volC) * Csum;
end

end
