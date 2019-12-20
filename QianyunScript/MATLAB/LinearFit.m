function linearModel = LinearFit(probeCoords)

% 3D Orthogonal Distance Regression (ODR) line, using SVD
N = size(probeCoords, 1);
avg = mean(probeCoords, 1);
linearModel.avg = avg;

dX = bsxfun(@minus, probeCoords, avg);
% [coeff,score,roots] = pca(probeCoords);
% dirVect = coeff(:,1);
% 
% t = [min(score(:,1))-.2, max(score(:,1))+.2];
% endpts = [meanX + t(1)*dirVect'; meanX + t(2)*dirVect'];
% plot3(endpts(:,1),endpts(:,2),endpts(:,3),'k-');


C=(dX'*dX)/(N-1);
[R,D]=svd(C,0);
linearModel.dirVect = R(:, 1)';
% L(t)= avg + t*R(:,1)', where t is a real number



% Get the R^2 value
D=diag(D);
R2=D(1)/sum(D);
linearModel.R2 = R2;


x=dX*R(:,1);    % project residuals on R(:,1) 
x_min=min(x);
x_max=max(x);
dx=x_max-x_min;
Xa=(x_min-0.05*dx)*R(:,1)' + avg;
Xb=(x_max+0.05*dx)*R(:,1)' + avg;
X_end=[Xa;Xb];

hold on
plot3(X_end(:,1),X_end(:,3),X_end(:,2),'-g','LineWidth',3) % best fit line


end