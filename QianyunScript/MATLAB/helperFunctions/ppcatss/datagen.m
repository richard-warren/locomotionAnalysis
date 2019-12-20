function [t,x,latgen]=datagen(N,nl,nx)
% This program generates data
% that have been used for illustration.
% There is a change in the relationship between
% the latent variables at the quarter and
% a change in the mean at the half.
%
%     input: N - number of samples
%            nl - number of latent variables
%            nx - number of observable variables
%     output: t - time (columnt vector)
%             x - observable variables (samples in the rows)
%             latgen - latent variables (samples in the rows)


x=[]; 
latgen=[]; 

% Sample time interval of the latent variables
% (i.e. the difference between two changes)
swt=[120 220 70 300];
    
% Generate the latent variables
for j=1:nl 
    n=ceil(N/swt(j));
    dumm=[];
    for i=1:n 
        dumm=[dumm; (ones(swt(j),1)*(rand-0.5)*7)];
    end    
    latgen=[latgen dumm(1:N)] ;
end   

% Change in the mean at the half
latgen(round(length(latgen)/2):end,:)=latgen(round(length(latgen)/2):end,:)+3;
    
% Additional white "noise"
latgen=latgen+rand(size(latgen))*2-1;

% Change in the relationship between
% the latent variables at the quarter
inter=[0 ceil(N/4) N];
relationship=rand(N*2,nx);
[pc1, score1, latent1, tsquare1] = princomp(relationship(1:N,:));
Worig{1}=pc1(:,1:nl)*sqrt(diag(latent1(1:nl)));
[pc2, score2, latent2, tsquare2] = princomp(relationship(N+1:2*N,:));
Worig{2}=pc2(:,1:nl)*sqrt(diag(latent1(1:nl)));
    
% Generate the observable variables
for i=1:length(inter)-1
    x=[x;(Worig{i}*latgen(inter(i)+1:inter(i+1),:)')'; ]; 
end  

% Additional white noise
x=x+rand(size(x))*0.1-0.05;
    
t=[1:N]';

% Plot the latent variables
figure(1)
for j=1:nl
    subplot(nl,1,j)
    plot(latgen(:,j))
    xlabel('time')
    dumm=['y_' num2str(j)];
    ylabel(dumm)
end
