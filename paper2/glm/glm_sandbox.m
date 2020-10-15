%% play around with GLMs

%% inits

session = '181020_001'; cell = 83;  % 37    54    65    66    69    83
timeDegrees = 2;  % 0: .56, 1: .72, 2: .72, 3: .72

% make design matrix
[dmat, t] = makeDesignMatrix(session, 'timeDegrees', timeDegrees);

% load cells
neuralData = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'unit_ids', 'spkRates', 'timeStamps');
spkRate = interp1(neuralData.timeStamps, neuralData.spkRates(neuralData.unit_ids==cell,:), t);
inds = find(~isnan(spkRate),1,'first') : find(~isnan(spkRate),1,'last');  % valid inds for cell
tsub = t(inds);

% X = zscore(table2array(dmat(inds,:)));
X = zscore(table2array(dmat(inds,:)));
X(:,end-timeDegrees) = 1;  % restore the constant column to ones
y = spkRate(inds)';
w = (X'*X) \ X'*y;
yhat = X * w;


close all; figure('position', [2.00 1024.00 1278.00 332.00], 'color', 'white'); hold on
plot(tsub, y)
plot(tsub, yhat, 'linewidth', 2)
% set(gca, 'xlim', [626.9596  639.5402])
xlabel('time (s)'); ylabel('firing rate (hz)')
legend('actual', 'predicted')
fprintf('r squared: %.2f\n', corr(y, yhat)^2)

% figure; plot(X(:,end-timeDegrees:end))