function [comp,Spca]=compat(results,q,inic)  
% This program computes the compatibility matrix.

% The similarity of PPCA models based on the index by Krzanowski
for i=1:size(results.S,3)
    [M,sk]=eig(results.S(:,:,i));
    M=M(:,end-q+1:end);
    for j=1:size(results.S,3);
        [L,sk]=eig(results.S(:,:,j));
        L=L(:,end-q+1:end);
        Spca(i,j)=trace(L'*M*M'*L)/q; 
    end
end

% The euclidian distance of cluster centers in the featuring space
dd=[];
a=results.mu';
for i=1:size(a,2)
    d=a-repmat(a(:,i),1,size(a,2));
    dd=[dd; (diag(d'*d))'];
end

% Compatible cluster merging by Kaymak & Babuska
nu1=1/(inic*(inic-1))*sum(sum(Spca - diag(diag(Spca))));
nu2=1/(inic*(inic-1))*sum(sum(sqrt(dd)));
    
c1=1/(1-nu1).*Spca-nu1/(1-nu1);
index1=find(c1<0);
c1(index1)=0;
    
c2=-1/nu2.*(dd)+1;
index2=find(c2<0);
c2(index2)=0;
    
p=2;
comp=(1/2*(c1.^p + c2.^p)).^(1/p);
