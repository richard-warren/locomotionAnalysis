%% play around with GLMs :)


%% inits

session = '181020_001'; neuron = 66;  % 37    54    65    66    69    83
timeDegrees = 3;  % 0: .56, 1: .72, 2: .72, 3: .72

% make design matrix
[dmat, t] = makeDesignMatrix(session, 'timeDegrees', timeDegrees);

% load cells
neuralData = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'rewardTimes');
spkRate = interp1(neuralData.timeStamps, neuralData.spkRates(neuralData.unit_ids==neuron,:), t);
inds = find(~isnan(spkRate),1,'first') : find(~isnan(spkRate),1,'last');  % valid inds for cell
t = t(inds);
dt = t(2) - t(1);

spkTimes = neuralData.spkTimes{neuralData.unit_ids==neuron};
[spkRateGaus, tTemp] = getFiringRate(spkTimes, 'fs', 1000, 'kernel', 'gauss');
spkRateGaus = interp1(tTemp, spkRateGaus, t);

% cross validation splits
% (each epoch extends from one reward to the next // assign each epoch to
% one of k partitions randomly)

% settings
k = 5;

% epochs spanning beginning of one reward to the end of the next
r = rewardTimes(rewardTimes>t(1) & rewardTimes<t(end))';
r = [t(1) r t(end)];

% assign each epoch to a parition p
partitions = true(length(t), k);  % logical masks for each partition
[~,~,p] = histcounts(randperm(length(r)-1), linspace(1,length(r)-1,k+1));  % evey ind assigned an int on [1,k]
for i = 1:k
    for j = find(p==i)
        partitions(t>=r(j) & t<=r(j+1), i) = false;
    end
end

disp('init complete!')

%% linear glm (manual)

lambdas = [0 2.^(0:15)];


X = zscore(table2array(dmat(inds,:)));
X(:,end-timeDegrees) = 1;  % restore the constant column to ones
Imat = eye(size(X,2));
y = spkRate(inds)';

r2 = nan(1,length(lambdas));
for i = 1:length(lambdas)
    disp(i)
    reg = Imat*lambdas(i);
    lr2 = nan(1,k);
    for j = 1:k
        b = partitions(:,j);
        w = (X(b,:)' * X(b,:) + reg) \ X(b,:)' * y(b,:);
        yhat = X*w;  % prediction for held out data
        lr2(j) = corr(yhat(~b), y(~b))^2;
    end
    r2(i) = mean(lr2);
end

% no cross validation
% w = (X'*X) \ X'*y;
% yhat = X*w;  % prediction for held out data
% r2 = corr(yhat, y);


close all; figure('position', [2.00 1024.00 1278.00 332.00], 'color', 'white'); hold on
plot(t, y)

yhatMasked = yhat; yhatMasked(~b) = nan;
plot(t, yhatMasked, 'linewidth', 2)
yhatMasked = yhat; yhatMasked(b) = nan;
plot(t, yhatMasked, 'linewidth', 2)

xlabel('time (s)'); ylabel('firing rate (hz)')
legend('actual', 'train', 'test')
fprintf('r squared: %.2f\n', mean(r2))

%% poisson glm (glmfit)

X = zscore(table2array(dmat(inds,:)));
X(:,end-timeDegrees) = 1;  % restore the constant column to ones
spkTimes = neuralData.spkTimes{neuralData.unit_ids==neuron};
y = histcounts(spkTimes, t(1)-dt/2 : dt : t(end)+dt/2)';

tic; [w, dev, stats] = glmfit(X, y, 'poisson', 'constant', 'off'); toc
yhat = exp(X*w) / dt;

close all; figure('position', [2.00 1024.00 1278.00 332.00], 'color', 'white'); hold on
plot(t, spkRateGaus)
plot(t, yhat, 'linewidth', 2)
xlabel('time (s)'); ylabel('firing rate (hz)')
legend('actual', 'predicted')
fprintf('r squared: %.2f\n', corr(spkRateGaus', yhat)^2)


%% poisson glm (lassoglm)

X = zscore(table2array(dmat(inds,:)));
X(:,end-timeDegrees) = 1;  % restore the constant column to ones
spkTimes = neuralData.spkTimes{neuralData.unit_ids==neuron};
y = histcounts(spkTimes, t(1)-dt/2 : dt : t(end)+dt/2)';

tic; [B, fitInfo] = lassoglm(X, y, 'poisson', ...
    'alpha', .5, 'lambda', [0 1], 'standardize', true, 'CV', 2); toc

minInd = fitInfo.IndexMinDeviance;
minInd = 2;
B0 = fitInfo.Intercept(minInd);
w = [B0; B(:,minInd)];
yhat = glmval(w, X, 'log') / dt;

close all; figure('position', [2.00 1024.00 1278.00 332.00], 'color', 'white'); hold on
plot(t, spkRateGaus)
plot(t, yhat, 'linewidth', 2)
xlabel('time (s)'); ylabel('firing rate (hz)')
legend('actual', 'predicted')
fprintf('r squared: %.2f\n', corr(spkRateGaus', yhat)^2)

%% glmnet package (no cross validation, no regularizaton)

X = zscore(table2array(dmat(inds,:)));
X(:,end-timeDegrees) = 1;  % restore the constant column to ones
y = histcounts(spkTimes, t(1)-dt/2 : dt : t(end)+dt/2)';

% git 
options = struct('alpha', 0, 'lambda', 0, 'standardize', false);
fit = glmnet(X, y, 'poisson', options);

% plot predicted vs true firing rate
yhat = glmnetPredict(fit, X, [], 'response') / dt;
close all; figure('position', [2.00 1024.00 1278.00 332.00], 'color', 'white'); hold on
plot(t, spkRateGaus)
plot(t, yhat, 'linewidth', 2)
xlabel('time (s)'); ylabel('firing rate (hz)')
legend('actual', 'predicted')
fprintf('r squared:          %.2f\n', corr(spkRateGaus', yhat)^2)
fprintf('deviance explained: %.2f\n', fit.dev)


%% glmnet package (cross validation, regularizaton)

% settings
lambdas = logspace(-8, 1, 50);
% options = struct('alpha', 0, 'nlambda', 5, 'lambda_min', 1e-20, 'standardize', false);
options = struct('alpha', 0, 'lambda', lambdas, 'standardize', false);

% inits
X = zscore(table2array(dmat(inds,:)));
X(:,end-timeDegrees) = 1;  % restore the constant column to ones
y = histcounts(spkTimes, t(1)-dt/2 : dt : t(end)+dt/2)';

fold_id = nan(size(t));
for i = 1:size(partitions,1)
    fold_id(i) = find(~partitions(i,:));
end

% fit model
tic; fit = cvglmnet(X, y, 'poisson', options, [], k, fold_id);
fprintf('\nfit model in %.1f minutes\n', toc/60);
lambdaMinInd = find(fit.lambda==fit.lambda_min);

% plot predicted vs true firing rate
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
set(gca, 'xscale', 'log')
xlabel('lambda'); ylabel('cross validated error')














