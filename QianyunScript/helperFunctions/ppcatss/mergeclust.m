function results2=mergeclust(results,q,index)
% Cluster merging based on the method of P. M. Kelly.

[nsamples dim]=size(results.x);
x1=ones(nsamples,1);
nclusters=size(results.pc,2);

% a priori probabilites of the cluster
pc=sum(results.pc(index:index+1));

% cluster center and covariance matrix
m(1,:)=results.mu(index,:);
m(2,:)=results.mu(index+1,:);

    n1=nsamples*results.pc(index);
    n2=nsamples*results.pc(index+1);
    sumn=n1+n2;
    mu=n1/sumn*m(1,:)+n2/sumn*m(2,:);
    S=(n1-1)/(sumn-1)*results.S(:,:,index)+(n2-1)/(sumn-1)*results.S(:,:,index+1) ...
      + (n1*n2)/(sumn*(sumn-1))*( (m(1,:)-m(2,:))'*(m(1,:)-m(2,:)) );
        

% variance and weight matrix
[V,D]=eig(S);
[Y,I]=sort(-diag(D));
D=diag(D);
D=D(I);
V=V(:,I);

szoras=max(1/(dim-q)*sum(D(q+1:end)),1e-3);
W=V(:,1:q)*sqrt(diag(D(1:q)-szoras));

% p(x_k | \eta_i)
xv=results.x-x1*mu;
C=szoras*eye(dim,dim)+W*W';
E=sum((xv*pinv(C).*xv)')';

px=(2*pi)^(-dim/2)*det(C)^(-1/2)*exp(-1/2*E);



% cluster center and variance in time
mt=results.mut(index:index+1);

    mut=n1/sumn*mt(1)+n2/sumn*mt(2);
    szorast=(n1-1)/(sumn-1)*results.szorast(index)+(n2-1)/(sumn-1)*results.szorast(index+1) ...
      + (n1*n2)/(sumn*(sumn-1))*( (mt(1)-mt(2))'*(mt(1)-mt(2)) );


% p(t_k | \eta_i)
pt = prod((exp(-(results.t-kron(x1,mut)).^2 ...
  ./ (2*kron(x1,szorast)))./ kron(x1,sqrt(2*pi*szorast))),2);




% Refresh the results and the indeces
results2.W(:,:,1:index-1)=results.W(:,:,1:index-1);
results2.W(:,:,index)=W;
results2.W(:,:,index+1:nclusters-1)=results.W(:,:,index+2:nclusters);

results2.mu(1:index-1,:)=results.mu(1:index-1,:);
results2.mu(index,:)=mu;
results2.mu(index+1:nclusters-1,:)=results.mu(index+2:nclusters,:);

results2.szoras(1:index-1)=results.szoras(1:index-1);
results2.szoras(index)=szoras;
results2.szoras(index+1:nclusters-1)=results.szoras(index+2:nclusters);

results2.pc(1:index-1)=results.pc(1:index-1);
results2.pc(index)=pc;
results2.pc(index+1:nclusters-1)=results.pc(index+2:nclusters);

results2.px(:,1:index-1)=results.px(:,1:index-1);
results2.px(:,index)=px;
results2.px(:,index+1:nclusters-1)=results.px(:,index+2:nclusters);

results2.pt(:,1:index-1)=results.pt(:,1:index-1);
results2.pt(:,index)=pt;
results2.pt(:,index+1:nclusters-1)=results.pt(:,index+2:nclusters);

