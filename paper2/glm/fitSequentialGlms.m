function [models, fitdata] = fitSequentialGlms(session, neuron, varargin)

% train null GLMs, and gradually add predictor groups one by one until we
% get to the full model.

% settings
s.lambdas = logspace(-8, -1, 40);    % ridge regression coefficients
s.folds = 5;                         % cross-validation folds
s.parallel = false;                  % whether crossval analyses are parallelized
s.verbose = true;
s.save = true;
s.outputFileName = fullfile(...
    getenv('SSD'), 'paper2', 'modelling', 'glms', 'sequential_glms', ...
    [session '_cell_' num2str(neuron) '_glm.mat']);



% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
if s.verbose; fprintf('%s: fitting models for neuron %i... ', session, neuron); end

% load design matrix
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'designMatrices', 'residual', [session '_designMatrix.mat']), ...
    'dmat', 't', 'reward_all')

% load predictor info
% TODO: update settings spreadsheet with groups and order of interest
filename = fullfile(getenv('GITDIR'), 'locomotionAnalysis', 'paper2', 'glm', 'settings', 'sequential_predictorSettings.xlsx');
predictorInfo = readtable(filename, 'sheet', 'predictors', 'ReadRowNames', true);
predictorInfo = predictorInfo(dmat.Properties.VariableNames(1:end-1), :);  % end-1 because last predictor is time in dmat
groupInfo = readtable(filename, 'sheet', 'groups', 'ReadRowNames', true);
groups = groupInfo.group;

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
[~,~,f] = histcounts(randperm(length(r)-1), linspace(1,length(r)-1,s.folds+1));  % evey ind assigned an int on [1,k]
for i = 1:s.folds
    for j = find(f==i)                     % for each reward epoch in fold
        fold_id(t>=r(j) & t<=r(j+1)) = i;  % assign each time point within epoch to fold i
    end
end


% initialize table
% TODO: update col and row names (what needs to be stored?)
nrows = length(groups) + 1;
rowNames = ['null', groups'];
models = table(cell(nrows,1), nan(nrows,1), ...
    'RowNames', rowNames, 'VariableNames', {'model', 'dev'});

% get names of columns in X
colNames = cell(1,length(dmat.Properties.VariableNames));
for i = 1:length(dmat.Properties.VariableNames)
    colNames{i} = repelem({dmat.Properties.VariableNames{i}}, length(dmat{1,i}));
end
colNames = cat(2, colNames{:});

% train full model
% (the full model corresponds to the addition of the final group)
X_full = table2array(dmat);
keyboard
modelFull = fitModel(X_full, y, s.lambdas);
models{groups{end}, 'model'}{1} = modelFull;
models{groups{end}, 'dev'} = cvdeviance(X_full, y, modelFull, 'holdout', true, 'bestLambdaOnly', true);
lambda_min = modelFull.lambda_min;

% train time only model
colBins = strcmp(colNames, 'time');
X = X_full(:, colBins);
model = fitModel(X, y, lambda_min);
models{'null', 'model'}{1} = model;
models{'null', 'dev'} = cvdeviance(X, y, model, 'holdout', true, 'bestLambdaOnly', true);

% check that rick's deviance matches glmnet deviance for full model evaluated on training data
glmnet_dev = modelFull.glmnet_fit.dev(modelFull.lambda_min_id);
rick_dev = cvdeviance(X_full, y, modelFull, 'holdout', true, 'bestLambdaOnly', true);  % don't use heldout data to match glmnet_fit deviance (not sure if this is correct...)
if abs(glmnet_dev - rick_dev)>.02
    fprintf('WARNING (%s_unit_%i)! Deviance diff of %.3f between glmnet (%.3f) and rick (%.3f)!\n', ...
        session, neuron, (glmnet_dev - rick_dev), glmnet_dev, rick_dev)
end


% add groups one by one! (note we already computed the final, full model)
for i = 1:(length(groups)-1)
    % TODO: compute models for reward and non-reward times separately?
    
    % find predictors belonging to groups 1:i
    members = predictorInfo.name(ismember(predictorInfo.group, groups(1:i)));
    
    % train model
    colBins = ismember(colNames, [members; 'time']);
    X = X_full(:, colBins);
    model = fitModel(X, y, lambda_min);
    models{groups{i}, 'model'}{1} = model;
    models{groups{i}, 'dev'} = cvdeviance(X, y, model, ...
        'holdout', true, 'bestLambdaOnly', true);  % should i subtract the null deviance?
end

if s.verbose; disp('all done!'); end


% save some objects for convenience...
yhat = exp(models{groups{end}, 'model'}{1}.fit_preval) / dt;
fitdata.t = t;
fitdata.y = y;
fitdata.yRate = spkRate;
fitdata.yhat = yhat;
fitdata.groups = groups;
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
    
%     if any(all(X==0,1)); disp('WARNING! Column with all zero predictors!!!'); end  % check for all-zero predictors
    
    options = struct('alpha', 0, 'lambda', lambdas, 'standardize', true);
    fit = cvglmnet(X(bins,:), y(bins), 'poisson', options, [], s.folds, ids, s.parallel, true);
    fit.fit_preval = fit.fit_preval(:,fit.lambda_min_id);  % only keep predictions for best lambda (save disk space)
end


end



