function [B,C]=commMatricesTxT(idx,K,T,A,zn)

B=zeros(T);
C=zeros(T);

B0=zeros(T);
C0=zeros(T);
for k=1:K
a=find(idx==k);
volC=vol(a);
b=find(idx~=k);
m=0;
for i=1:length(a)
    for j=i:length(a)
        m=m+1;
        Btemp(:,:,m)=(zn{a(i)}-zn{a(j)})*(zn{a(i)}-zn{a(j)})';
    end
end
B=B0+(1/volC)*sum(Btemp,3);
clear Btemp
B0=B;

m=0;
for i=1:length(a)
    for j=1:length(b)
        m=m+1;
        Ctemp(:,:,m)=(zn{a(i)}-zn{b(j)})*(zn{a(i)}-zn{b(j)})';
    end
end
C=C0+(1/volC)*sum(Ctemp,3);
clear Ctemp
C0=C;
end

    function volC= vol(a)
         volC= sum(sum(A(a,:)));
    end

end