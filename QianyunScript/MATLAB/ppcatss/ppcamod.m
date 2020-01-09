function [results]=ppcamod(t,x,inic,q,niterations,timew,m);
% Clustering algorithm for segmentation of 
% multivariate time-series.

min_var = 1e-4;

[nsamples dim]=size(x);
[nsamples dimt]=size(t);
x1=ones(nsamples,1);

x=(x-kron(mean(x),ones(size(x,1),1)));
x=x./(kron(std(x),ones(size(x,1),1)));

% Initializing
if isfield(inic,'W')
    nclusters=size(inic.pc,2);

    W=inic.W;
    mu=inic.mu;
    szoras=inic.szoras;
    pc=inic.pc;
    px=inic.px;
    pt=inic.pt;
else
    nclusters=inic;
    
    [PC, SCORE, LATENT, TSQUARE] = princomp(x);
    for cluster=1:nclusters
        W(:,:,cluster)=PC(:,1:q);
    end
    mu=ones(nclusters,dim);
    szoras=ones(1,nclusters);
    pc=ones(1,nclusters)/nclusters;

    px=zeros(nsamples,nclusters);
    dumm=ceil(nsamples/nclusters);
    for i=1:nclusters
        px((i-1)*dumm+1:min(i*dumm,length(x)),i)=1;
    end
    pt=px;
end    


% Iteration
step = 0;
ppold=zeros(nsamples,nclusters);
while (step <= niterations)
    
   if m==1 
     if timew
        p = sum((pt .*px .* kron(x1,pc)),2);
        pp = (pt .* px .* kron(x1,pc)) ./ (kron(p,ones(1,nclusters))+eps);
    else
        p = sum((px .* kron(x1,pc)),2);
        pp = (px .* kron(x1,pc)) ./ (kron(p,ones(1,nclusters))+eps);
    end   
  else      
     if timew
        p = sum((pt .*px .* kron(x1,pc)),2).^(m-1);
        pp = (pt .* px .* kron(x1,pc)).^(m-1) ./ (kron(p,ones(1,nclusters))+eps);
    else
        p = sum((px .* kron(x1,pc)),2).^(m-1);
        pp = (px .* kron(x1,pc)).^(m-1) ./ (kron(p,ones(1,nclusters))+eps);
    end   
  end      

         pc = sum(pp) / nsamples;

    if max(max(abs(pp-ppold))) < min_var
        break
    end
       
        
        ppold=pp;
        Wold=W;

        
        pp=pp.^m;
        sumpp=sum(pp);
        
        
    for cluster = 1:nclusters
        

        
        
        M(:,:,cluster)=szoras(cluster)*eye(q,q)+W(:,:,cluster)'*W(:,:,cluster);
        varhatx=[pinv(M(:,:,cluster))*W(:,:,cluster)'*(x-x1*mu(cluster,:))']';
        mu(cluster,:)=sum((x-(W(:,:,cluster)*varhatx')' ) ...  
                      .* kron(pp(:,cluster),ones(1,dim))) / (sumpp(cluster));
        for d=1:dim
            S(d,:,cluster) =  sum( (x - kron(x1,mu(cluster,:))) .* ...
                 (kron(x(:,d),ones(1,dim)) - ones(nsamples,dim)*mu(cluster,d)) ...
                 .* kron(pp(:,cluster),ones(1,dim))) / (sumpp(cluster));
        end
        S(:,:,cluster)=S(:,:,cluster)+eye(size(S(:,:,cluster),1))*1e-5;
        
        W(:,:,cluster)= S(:,:,cluster)*W(:,:,cluster)*pinv(szoras(:,cluster)*eye(q,q)+pinv(M(:,:,cluster))*W(:,:,cluster)'*S(:,:,cluster)*W(:,:,cluster));
        szoras(:,cluster)=max(1/dim*trace(S(:,:,cluster)-S(:,:,cluster)*Wold(:,:,cluster)*pinv(M(:,:,cluster))*W(:,:,cluster)'),1e-3);
    end %for cluster = 1:nclusters
        
    
    for cluster = 1:nclusters
        xv=x-x1*mu(cluster,:);
        C(:,:,cluster)=szoras(cluster)*eye(dim,dim)+W(:,:,cluster)*W(:,:,cluster)';
        E(:,:,cluster)=sum((xv*pinv(C(:,:,cluster)).*xv)')';
    end %for cluster = 1:nclusters
    
    for cluster = 1:nclusters
        px(:,cluster)=(2*pi)^(-dim/2)*det(C(:,:,cluster))^(-1/2)*exp(-1/2*E(:,:,cluster));
    end %for cluster = 1:nclusters

    
    
    for cluster = 1:nclusters
         mut(cluster,:) = sum(t .* kron(pp(:,cluster),ones(1,dimt))) / (sumpp(cluster));
         szorast(cluster,:) = 1e-4 + sum((t - kron(x1,mut(cluster,:))).^2 .* ...
                      kron(pp(:,cluster),ones(1,dimt))) / (sumpp(cluster));
     
         pt(:,cluster) = prod((exp(-(t-kron(x1,mut(cluster,:))).^2 ...
            ./ (2*kron(x1,szorast(cluster,:))))./ kron(x1,sqrt(2*pi*szorast(cluster,:)))),2);
    end %for cluster = 1:nclusters
    
    % Figure
    if 1
        a=pt./repmat(max(pt),nsamples,1);
        suma=sum(a,2);
        figure(5)
        subplot(3,1,1)
        plot(t,a./repmat(suma,1,size(a,2)))
        axis([t(1) t(end) 0 1])
        ylabel('\beta_i(t_k)')
        subplot(3,1,2)
        plot(t,a)
        axis([t(1) t(end) 0 1])
        ylabel('A_i(t_k)')
        subplot(3,1,3)
        plot(t,pp)
        axis([t(1) t(end) 0 1])
        ylabel('p({\bf z}_k | \eta_i)')
        xlabel('time')
        drawnow
    end
    
step = step + 1
pause(0.01)
end % Iteration

[y ind]=sort(mut);
results.W(:,:,:)=W(:,:,ind);
results.S(:,:,:)=S(:,:,ind);
results.C(:,:,:)=C(:,:,ind);
results.mu=mu(ind,:);
results.mut=mut(ind);
results.szoras=szoras(:,ind);
results.szorast=szorast(ind,:);
results.pc=pc(ind);
results.px=px(:,ind);
results.pt=pt(:,ind);
results.pp=pp(:,ind);
results.t=t;
results.x=x;


