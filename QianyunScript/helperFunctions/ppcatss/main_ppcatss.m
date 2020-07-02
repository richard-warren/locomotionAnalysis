clear all
close all
% This program is able to segment multivariate
% time-series based on modified Gath-Geva fuzzy
% clustering method (mixture of Probabilistic
% Principal Component models).
%
% Janos Abonyi (abonyij@fmt.vein.hu) and
% Balazs Feil (feilb@fmt.vein.hu) 2004
%

rand('seed',564)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .: Initial parameters :. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inic=10;     % Initial number of segments
m=2;         % Fuzzyness parameter
thres=0.75;  % Threshold for compatible cluster merging
stop1=100;   % After every "stop1" number of iterations 
             % the cluster merging criterion will be analyzed 
stop2=100;   % Number of iterations after there are
             % no mergeable cluster

inicBU=7;      % Number of segments for Bottom-Up algorithm
inicfromBU=0;  % Initialize the PPCA-TSS based on the results of Bottom-Up
inicfromT2=0;  % If inicfromBU==1, then based on Hotelling T^2 (1) or Q reconstruction error (0)

%%%%%%%%%%%%%%%%%%%%%%%%%
% .: Data generation :. %
%%%%%%%%%%%%%%%%%%%%%%%%%
N=2000;  % Number of samples
nl=2;    % Number of latent variables
nx=6;    % Number of observable variables

[t,x,latgen]=datagen(N,nl,nx);
    
% Another dataset can be loaded
%         samples in the rows
%         time in the first column
% load data
% t=data(:,1);
% x=data(:,2:end);
% nx=size(x,2);

% Plot the observable variables
figure(2)
for j=1:nx
    subplot(nx,1,j)
    plot(t,x(:,j))
    xlabel('time')
    dumm=['x_' num2str(j)];
    ylabel(dumm)
    ymin=min(x(:,j));
    ymax=max(x(:,j));
    dumm=0.05*(ymax-ymin);
    axis([0 N ymin-dumm ymax+dumm]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .: Number of principal components :. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncov=inic;
dumm=[1 ceil([1:ncov].*N./ncov)];
    
for i=1:ncov
    eigvect=[]; eigval=[];
    [eigvect eigval]=eig(cov(x(dumm(i):dumm(i+1),:)));
    eigval=-1*(sort(-1*diag(eigval)));
    figure(3)
    subplot(2,1,1)
    plot(eigval,'-o')
    ylabel('eigenvalues')
    hold on
    title('Screeplot')
    subplot(2,1,2)
    plot(cumsum(eigval)/sum(eigval),'-o')
    xlabel('number of eigenvalues')
    ylabel('rate of cumulative sum')
    hold on
end
drawnow

disp('See the screeplot.')
q = input('Number of principal components will be ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .: Bottom-Up algorithm :. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag=[0 1];    % 0:T2, 1:Q  for Bottom-Up algorithm
for i=1:length(flag)
    segment=[]; tc=[]; index=[];
    [segment,tc] = pcaseg(x,inicBU,q,flag(i));
    for k=1:size(segment,2);
        for l=1:size(segment,2);
            L=segment(k).pc(:,1:q);
            M=segment(l).pc(:,1:q);
            Spca=trace(L'*M*M'*L)/q;
            index=[index; [k l Spca]]; 
        end
    end

    if flag(i)==0
        segmentT2=segment;
        indexT2=index;
    elseif flag(i)==1
        segmentQ=segment;
        indexQ=index;
    end
end

    % Initialize the PPCA-TSS based on the results of Bottom-Up
    if inicfromBU
        segment=[];
        if inicfromT2
            segment=segmentT2;
        else
            segment=segmentQ;
        end
        for i=1:size([segment.lx],2)
            inic.W(:,:,i)=segment(i).pc(:,1:q);
            inic.mu(i,:)=segment(i).avg;
            inic.szoras(i)=sum(segment(i).latent(q+1:end));
            inic.pc(i)=(segment(i).rx-segment(i).lx)/size(x,1);
            inic.px(:,i)=zeros(size(x,1),1);
            inic.px(segment(i).lx:segment(i).rx-1,i)=1;
            inic.pt=inic.px;
        end
    end

    % Plot the results of Bottom-Up algorithm
    figure(4)
    subplot(2,1,1)
    hold on
    for j=1:size(segmentT2,2)
       line([segmentT2(j).lx  segmentT2(j).lx], [0 1]);
    end
    axis([t(1) t(end) 0 1])
    ylabel('Hotelling T^2')
    hold off
    subplot(2,1,2)
    hold on
    for j=1:size(segmentQ,2)
       line([segmentQ(j).lx  segmentQ(j).lx], [0 1]);
    end
    axis([t(1) t(end) 0 1])
    ylabel('Q Reconstruction error')
    hold off
    xlabel('time')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .: PPCA-TSS algorithm :. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indcomp=0; % index of comparable clusters
while ~isempty(indcomp)
    %[results]=ppcamod(t,x,inic,q,niterations,timew,m);
     [results]=ppcamod(t,x,inic,q,stop1,1,m);

    % The comapibility matrix
    if isfield(inic,'W')
        nclusters=size(inic.pc,2);
    else
        nclusters=inic;
    end
    [comp, Spca]=compat(results,q,nclusters);
    PPCAModelsSimilarity=diag(Spca,1)
    CompatibilityIndeces=diag(comp,1)
    aa=diag(comp,1);
    indcomp=find(aa>thres);
    if ~isempty(indcomp)
        indcomp=find( aa==max(aa) )
        inic=mergeclust(results,q,indcomp);
 %       inic=inic-1;  % Reinitializing
    end
end
[results]=ppcamod(t,x,inic,q,stop2,1,m);
PPCAModelsSimilarity=diag(Spca,1)
CompatibilityIndeces=diag(comp,1)
    
