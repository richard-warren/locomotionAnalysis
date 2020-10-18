% function fitNeuronGlm(session, neuron)

% fit (several) glms for a neuron: full model, single model for each
% predictor group, and models with each predictor group removed

% temp
session = '181020_001'; neuron = 69;  % 37    54    65    66    69    83

% settings
s.lambdas = logspace(-8, -1, 50);    % ridge regression coefficients
s.timeDegrees = 3;                   % fit timeDegrees order polynomial to account for drift over time
s.folds = 5;                         % cross-validation folds
s.parallel = false;
s.plot = true;


% load (or make) design matrix
[dmat, t, reward_all] = makeDesignMatrix(session, 'timeDegrees', s.timeDegrees);

% load neuron
neuralData = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [session '_neuralData.mat']));
spkRate = interp1(neuralData.timeStamps, neuralData.spkRates(neuralData.unit_ids==neuron,:), t);
inds = find(~isnan(spkRate),1,'first') : find(~isnan(spkRate),1,'last');  % valid inds for neuron
t = t(inds);
dt = t(2) - t(1);
spkTimes = neuralData.spkTimes{neuralData.unit_ids==neuron};

% prep data for model (restrict to valid times for neuron)
options = struct('alpha', 0, 'lambda', s.lambdas, 'standardize', true);
X = table2array(dmat(inds,:));
y = histcounts(spkTimes, t(1)-dt/2 : dt : t(end)+dt/2)';

% define cross-validation folds
r = reward_all(reward_all>t(1) & reward_all<t(end))';  % epochs spanning beginning of one reward to the end of the next
r = [t(1) r t(end)];
fold_id = nan(size(t));
[~,~,f] = histcounts(randperm(length(r)-1), linspace(1,length(r)-1,s.folds+1));  % evey ind assigned an int on [1,k]
for i = 1:s.folds
    for j = find(f==i)                     % for each reward epoch in fold
        fold_id(t>=r(j) & t<=r(j+1)) = i;  % assign each time point within epoch to fold i
    end
end

%% fit model
tic; fprintf('\nfitting model (%i folds, %i lambdas)... ', s.folds, length(s.lambdas));
fit = cvglmnet(X, y, 'poisson', options, [], s.folds, fold_id, s.parallel);
fprintf('finished in %.1f minutes\n', toc/60);
minInd = find(fit.lambda==fit.lambda_min);
lambdaMinInd = find(fit.lambda==fit.lambda_min);

% plot predicted vs true firing rate
if s.plot
    yhat = cvglmnetPredict(fit, X, [], 'response') / dt;
    close all; figure('position', [2.00 1024.00 1278.00 332.00], 'color', 'white'); hold on
    plot(t, spkRateGaus)
    plot(t, yhat, 'linewidth', 2)
    xlabel('time (s)'); ylabel('firing rate (hz)')
    legend('actual', 'predicted')
    fprintf('r squared:          %.2f\n', corr(spkRateGaus', yhat)^2)
    fprintf('deviance explained: %.2f\n', fit.glmnet_fit.dev(lambdaMinInd))

    % plot regularization
    figure('position', [214.00 155.00 609.00 430.00], 'color', 'white'); hold on
    plot(fit.lambda, fit.cvm, 'linewidth', 2)
    scatter(fit.lambda(minInd), fit.cvm(minInd), 20, 'red', 'filled')
    set(gca, 'xscale', 'log')
    xlabel('lambda'); ylabel('cross validated error')
end
