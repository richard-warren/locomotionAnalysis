%% play around with GLMs :)


%% inits

session = '181020_001';
neuron = 69;

load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'designMatrices', [session '_designMatrix.mat']), 'dmat', 't')
neuralData = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [session '_neuralData.mat']));
spkRate = interp1(neuralData.timeStamps, neuralData.spkRates(neuralData.unit_ids==neuron,:), t);
inds = find(~isnan(spkRate),1,'first') : find(~isnan(spkRate),1,'last');  % valid inds for neuron
spkRate = spkRate(inds);
t = t(inds);
dt = t(2) - t(1);

spkTimes = neuralData.spkTimes{neuralData.unit_ids==neuron};

X = table2array(dmat(inds,:));
y = histcounts(spkTimes, [t-dt/2 t(end)+dt/2])';


%% compute deviance

%% cv model

lambdas = logspace(-8, -1, 20);
parallel = true;
options = struct('alpha', 0, 'lambda', lambdas, 'standardize', true);
fit = cvglmnet(X, y, 'poisson', options, [], 5, [], parallel, true);

%%
deviance = cvdeviance(fit, 'holdout', false);

close all; figure('color', 'white', 'position', [2.00 722.00 1278.00 634.00]);
subplot(121); hold on
scatter(fit.glmnet_fit.dev, deviance);
plot([0 1], [0 1])
set(gca, 'xlim', [0 1], 'ylim', [0 1])
subplot(122)
histogram(fit.glmnet_fit.dev - deviance)
xlabel('deviance differences')













