%%%% Input: A adjacency, f graph signal, K number of clusters, T, alpha
% Code for GraFiCA "Learning Optimal Graph Filters for Clustering of
% Attributed Graphs", Meiby Ortiz-Bouza and Selin Aviyente
% Michigan State University, email: ortizbou@msu.edu


%% Learn filter Step
h_old=zeros(T,1);
for n=1:20
    S=B-gamma*C;
    [u,v]=eig(S);
    [~,ind] = sort(diag(v),'ascend');
    u = u(:,ind);
    hnew=u(:,1);
    H=0;
    for t=1:T
      Hnew=H+hnew(t)*diag(D.^(t-1));
      H=Hnew;  
    end
    f_tilde=U*H*U'*f;
    %%%community detection
    Af1=zeros(N);
    parfor i=1:N
        Af1(:,i)=vecnorm((f-f(i,:)),2,2); %matrix created from the norm of the difference of the attributes, called W in the paper
    end
    Af=normadj(Af1); %Wn in the paper
    W = Af-alpha*An;
    [idxn,NMI,ARI]=CommDetec(W,GT,K); %%% community detection step
    if norm(hnew-h_old)<1e-3
       break
    end    
    %%% create B and C for the next iteration
    [B,C]=commMatricesTxT(idxn,K,T,An,zn);
    h_old=hnew;
end



