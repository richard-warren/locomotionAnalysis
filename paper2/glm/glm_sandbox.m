%% play around with GLMs

%% inits

% session = '181020_001'; cell = 66;
session = '181020_001'; cell = 66;

% load raw predictors
load(['C:\Users\richa\Desktop\lab_files\paper2\modelling\predictors\' session '_predictors.mat'], 'predictors')
t = predictors.t{1};  % assumes first predictor is continuous
dt = t(2) - t(1);

% load cells
neuralData = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'unit_ids', 'spkRates', 'timeStamps');
spkRate = interp1(neuralData.timeStamps, neuralData.spkRates(neuralData.unit_ids==cell,:), t);
inds = find(~isnan(spkRate),1,'first') : find(~isnan(spkRate),1,'last');  % valid inds for cell
tsub = t(inds);

%% prepare predictors
nkernels = 10;
wiskTimes = predictors{'whiskerContact', 'data'}{1};
wisk = histcounts(wiskTimes, tsub(1)-dt/2 :dt : tsub(end)+dt/2);
[basis, tBasis] = makeCosBasis(-0, .25, nkernels, 'dt', dt);
wiskBases = nan(nkernels, length(tsub))';
for i = 1:nkernels; wiskBases(:,i) = conv(wisk, basis(i,:), 'same'); end

close all; figure('position', [2.00 722.00 1278.00 634.00], 'color', 'white')
subplot(311)
xlims = [tBasis(1) tBasis(end)] + predictors{'whiskerContact', 'data'}{1}(3);
plot(tsub, wiskBases)
set(gca, 'xlim', xlims)

% wisk only model

% design matrix
X = [wiskBases ones(length(tsub),1)];
y = spkRate(inds)';
w = (X'*X) \ X'*y;
yhat = X * w;


subplot(312); hold on
plot(tsub, yhat)
plot(tsub, y)
plot(repmat(wiskTimes,1,2)', ylim, 'color', 'black')
set(gca, 'xlim', xlims)

subplot(313); hold on
plot(tBasis, (basis .* w(1:end-1))')
set(gca, 'xlim', [tBasis(1) tBasis(end)])





