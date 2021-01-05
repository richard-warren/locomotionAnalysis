function [models, fitdata] = fitFinekinGlm(session, neuron, varargin)



% settings
s.grossKin = {'velocity', 'bodyAngle', 'buttHeight', 'bodyAngle_vel', 'buttHeight_vel', 'velocity_vel', 'paw4RH_stride'};
% s.grossKin = {'velocity', 'bodyAngle', 'buttHeight', 'bodyAngle_vel', 'buttHeight_vel', 'velocity_vel'};
s.fineKin = {'paw1LH_x', 'paw2LF_x', 'paw3RF_x', 'paw4RH_x', 'paw1LH_y', 'paw2LF_y', 'paw3RF_y', 'paw4RH_y', 'paw1LH_z', 'paw2LF_z', 'paw3RF_z', 'paw4RH_z', 'paw1LH_x_vel', 'paw2LF_x_vel', 'paw3RF_x_vel', 'paw4RH_x_vel', 'paw1LH_y_vel', 'paw2LF_y_vel', 'paw3RF_y_vel', 'paw4RH_y_vel', 'paw1LH_z_vel', 'paw2LF_z_vel', 'paw3RF_z_vel', 'paw4RH_z_vel'};
s.lambdas = logspace(-8, -1, 40);    % ridge regression coefficients
s.folds = 5;                         % cross-validation folds
s.parallel = false;                  % whether crossval analyses are parallelized
s.verbose = true;
s.save = true;
s.outputFileName = fullfile(...
    getenv('SSD'), 'paper2', 'modelling', 'glms', 'finekin_glms', ...
    [session '_cell_' num2str(neuron) '_glm.mat']);



% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
if s.verbose; fprintf('%s: fitting models for neuron %i... ', session, neuron); end

% load design matrix
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'designMatrices', [session '_designMatrix.mat']), ...
    'dmat', 't', 'reward_all')

% load neuron
neuralData = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [session '_neuralData.mat']));
spkRate = interp1(neuralData.timeStamps, neuralData.spkRates(neuralData.unit_ids==neuron,:), t);
inds = find(~isnan(spkRate),1,'first') : find(~isnan(spkRate),1,'last');  % valid inds for neuron
spkRate = spkRate(inds);
t = t(inds);
dt = t(2) - t(1);
spkTimes = neuralData.spkTimes{neuralData.unit_ids==neuron};
if isempty(spkTimes)
    fprintf('%s: WARNING! No spikes for neuron %i! Something is wrong. Aborting...\n', session, neuron);
    return
end

% prep data for model
dmat = dmat(inds,:);
y = histcounts(spkTimes, [t-dt/2 t(end)+dt/2])';


% define cross-validation folds
r = reward_all(reward_all>t(1) & reward_all<t(end))';  % epochs spanning beginning of one reward to the end of the next
r = [t(1) r t(end)];
fold_id = nan(size(t));
[~,~,f] = histcounts(randperm(length(r)-1), linspace(1,length(r)-1,s.folds+1));  % every ind assigned an int on [1,k]
for i = 1:s.folds
    for j = find(f==i)                     % for each reward epoch in fold
        fold_id(t>=r(j) & t<=r(j+1)) = i;  % assign each time point within epoch to fold i
    end
end


% initialize table
models = table(cell(4,1), nan(4,1), ...
    'RowNames', {'gross', 'fine', 'grossfine', 'full'}, 'VariableNames', {'model', 'dev'});

% get names of columns in X
colNames = cell(1,length(dmat.Properties.VariableNames));
for i = 1:length(dmat.Properties.VariableNames)
    colNames{i} = repelem({dmat.Properties.VariableNames{i}}, length(dmat{1,i}));  % number of 
end
colNames = cat(2, colNames{:});


% full
X_full = table2array(dmat);
model = fitModel(X_full, y, s.lambdas);
models{'full', 'model'}{1} = model;
models{'full', 'dev'} = cvdeviance(X_full, y, model, 'holdout', true, 'bestLambdaOnly', true);

% gross
colBins = ismember(colNames, s.grossKin);
X = X_full(:,colBins);
model = fitModel(X, y, s.lambdas);
models{'gross', 'model'}{1} = model;
models{'gross', 'dev'} = cvdeviance(X, y, model, 'holdout', true, 'bestLambdaOnly', true);

% fine
colBins = ismember(colNames, s.fineKin);
X = X_full(:,colBins);
model = fitModel(X, y, s.lambdas);
models{'fine', 'model'}{1} = model;
models{'fine', 'dev'} = cvdeviance(X, y, model, 'holdout', true, 'bestLambdaOnly', true);

% gross + fine
colBins = ismember(colNames, [s.grossKin s.fineKin]);
X = X_full(:,colBins);
model = fitModel(X, y, s.lambdas);
models{'grossfine', 'model'}{1} = model;
models{'grossfine', 'dev'} = cvdeviance(X, y, model, 'holdout', true, 'bestLambdaOnly', true);

if s.verbose; disp('all done!'); end


% save some objects for convenience...
fitdata.t = t;
fitdata.y = y;
fitdata.yRate = spkRate;
fitdata.session = session;
fitdata.neuron = neuron;

% save
if s.save; save(s.outputFileName, 'models', 'fitdata'); end


function fit = fitModel(X, y, lambdas, bins)
    % fit model with k fold cross validations
    if nargin<4; bins = true(1,size(X,1)); end
    if length(lambdas)==1; lambdas = [0 lambdas]; end  % a hack, because cvglmnet requires multiple lambdas
    
    % temporary hack to ensure all folds are represented
    % (should really re-assign temporally discontinuous chunks to different folds)
    if length(unique(fold_id(bins)))==s.folds
        ids = fold_id(bins);
    else
        ids = [];
    end
    
    if any(all(X==0,1)); disp('WARNING! Column with all zero predictors!!!'); end  % check for all-zero predictors
    
    options = struct('alpha', 0, 'lambda', lambdas, 'standardize', true);
    fit = cvglmnet(X(bins,:), y(bins), 'poisson', options, [], s.folds, ids, s.parallel, true);
    fit.fit_preval = fit.fit_preval(:,fit.lambda_min_id);  % only keep predictions for best lambda (save disk space)
end


end



