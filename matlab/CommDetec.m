function [idx,nmi,ari,cost]=CommDetec(W,GT,K)

%% Initial Communities
[Vnew,~] = eigs(W,K,'sa');  %% eigendecomposition, getting the K-smallest eigenvalues and their corresponding eigenvecctors
cost=trace(Vnew'*W*Vnew);
for m=1:1   % change if want to repeat k-means more times and select the best results
[idx2(:,m),ctrs] = kmeans(Vnew,K,'emptyaction','singleton','replicates',500);
nmi2(m)=getNMI(idx2(:,m),GT);  % compute NMI
ari(m)=rand_index(idx2(:,m),GT,'adjusted'); % compute ARI
end
[nmi,m]=max(nmi2);
idx=idx2(:,m);   % community labels
ari=max(ari);
end