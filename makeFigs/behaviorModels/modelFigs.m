

% get data for model
experiment = 'baseline';

load([getenv('OBSDATADIR') 'matlabData\' experiment 'KinematicData.mat'], 'data', 'touchClassNames');
kinData = data; clear data;
modelData = getModelData(kinData, touchClassNames);
disp([experiment ' kinematic data loaded!'])


%% simple logistic regression

dv = 'isSuccess';

% set up predictor and reponse matrices
X = cat(2, [modelData.obsPos]', ...
           [modelData.swingStartDistance]', ...
           [modelData.vel]', ...
           [modelData.angle]', ...
           [modelData.obsHeightsVid]', ...
           [modelData.isLightOn]');
Y = [modelData.(dv)]'; % predict whether mouse takes one big step and whether trial is successful

[B,dev,stats] = mnrfit(X, categorical(Y));

% plot results
xScat = cat(2,X,ones(size(X,1),1)) * B;
yScat = Y;
x = min(xScat):.001:max(xScat);
y = (1 + exp(x)).^(-1);

disp(mean(((1 + exp(xScat)).^(-1)>1) == yScat))

%%
figure('Color', 'white', 'Position', [900 400 800 300]);
plot(x,y); hold on;
scatter(xScat, yScat);
set(gca, 'YLim', [-.1 1.1], 'Box', 'off');
title(dv)






