function [models, fitdata] = fitResidualGlm(session, neuron, varargin)

% for each group of predictors, trained GLM on RESIDUALS after regressing
% away user-specified groups

% temp
% session = '181020_001';
% neuron = 69;

% settings
s.lambdas = logspace(-8, -1, 40);    % ridge regression coefficients
s.folds = 5;                         % cross-validation folds
s.parallel = false;                  % whether crossval analyses are parallelized
s.verbose = true;
s.save = true;
s.outputFileName = fullfile(...
    getenv('SSD'), 'paper2', 'modelling', 'glms', 'residual_glms', ...
    [session '_cell_' num2str(neuron) '_glm.mat']);



% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
if s.verbose; fprintf('%s: fitting models for neuron %i... ', session, neuron); end

% load design matrix
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'designMatrices', [session '_designMatrix.mat']), ...
    'dmat', 't', 'reward_all')

% load predictor info
filename = fullfile(getenv('GITDIR'), 'locomotionAnalysis', 'paper2', 'glm', 'predictorSettings.xlsx');
predictorInfo = readtable(filename, 'sheet', 'predictors', 'ReadRowNames', true);
predictorInfo = predictorInfo(dmat.Properties.VariableNames(1:end-1), :);  % end-1 because last predictor is time in dmat
groups = unique(predictorInfo.group);
groupInfo = readtable(filename, 'sheet', 'groups', 'ReadRowNames', true);

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
% y = y + randi(2, length(y), 1);


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
nrows = length(groups) + 1;
rowNames = ['full', groups'];
models = table(cell(nrows,1), cell(nrows,1), cell(nrows,1), nan(nrows, 1), nan(nrows,1), true(nrows, length(t)), ...
    'RowNames', rowNames, 'VariableNames', {'groups_pre', 'model_pre', 'model', 'dev_pre', 'dev', 'bins'});

% get names of columns in X
colNames = cell(1,length(dmat.Properties.VariableNames));
for i = 1:length(dmat.Properties.VariableNames)
    colNames{i} = repelem({dmat.Properties.VariableNames{i}}, length(dmat{1,i}));  % number of 
end
colNames = cat(2, colNames{:});


% train full model
X_full = table2array(dmat);
modelFull = fitModel(X_full, y, s.lambdas);
models{'full', 'model'}{1} = modelFull;
models{'full', 'dev'} = cvdeviance(X_full, y, modelFull, 'holdout', true, 'bestLambdaOnly', true);
lambda_min = modelFull.lambda_min;

% train time only model
colBins = strcmp(colNames, 'time');
X = X_full(:, colBins);
model = fitModel(X, y, lambda_min);
models{'full', 'model_pre'}{1} = model;
models{'full', 'dev_pre'} = cvdeviance(X, y, model, 'holdout', true, 'bestLambdaOnly', true);

% check that rick's deviance matches glmnet deviance for full model evaluated on training data
glmnet_dev = modelFull.glmnet_fit.dev(modelFull.lambda_min_id);
rick_dev = cvdeviance(X_full, y, modelFull, 'holdout', false, 'bestLambdaOnly', true);  % don't use heldout data to match glmnet_fit deviance
if abs(glmnet_dev - rick_dev)>.02
    fprintf('WARNING! Deviance computed by Glmnet %.3f off from rick deviance!\n', ...
        glmnet_dev - rick_dev)
end




% assess importance of groups of features
for i = 1:length(groups)
%     disp(groups{i})
    
    % determine time bins relevant to single group models
    % deviance will only be computed at these times
    restrictTime = ~isnan(groupInfo{groups{i}, 'restrict_time'}{1});
    
    if restrictTime
        colBins = ismember(colNames, groupInfo{groups{i}, 'restrict_time'});  % only restrict time based on predictor listed in restrict_time column
        timeBins = any(X_full(:,colBins)>0,2);
        models{groups{i}, 'bins'} = timeBins';
    else
        timeBins = true(size(X_full,1),1);
    end
    
    if any(timeBins)
    
        % find predictors for 'pre' group
        members = {'time'};  % members of pre group
        if ~isnan(groupInfo{groups{i}, 'residuals'}{1})
            preGroups = split(groupInfo{groups{i}, 'residuals'}{1}, ' ');
            models{groups{i}, 'groups_pre'}{1} = preGroups;
            members = [members; predictorInfo.Properties.RowNames(ismember(predictorInfo.group, preGroups))];
        end

        % train pre model
        preColBins = ismember(colNames, members);
        X = X_full(:, preColBins);
        model = fitModel(X, y, lambda_min, timeBins);
        models{groups{i}, 'model_pre'}{1} = model;
        models{groups{i}, 'dev_pre'} = cvdeviance(X(timeBins,:), y(timeBins), model, ...
            'holdout', true, 'bestLambdaOnly', true);

        % train model
        members = predictorInfo.Properties.RowNames(strcmp(predictorInfo.group, groups{i}));  % members of group
        colBins = ismember(colNames, members) | preColBins;
        X = X_full(:, colBins);
%         if i==5; keyboard; end
        model = fitModel(X, y, lambda_min, timeBins);
        models{groups{i}, 'model'}{1} = model;
        dev = cvdeviance(X(timeBins,:), y(timeBins), model, ...
            'holdout', true, 'bestLambdaOnly', true);
        models{groups{i}, 'dev'} = dev - models{groups{i}, 'dev_pre'};
    end
end

if s.verbose; disp('all done!'); end


% save some objects for convenience...
yhat = exp(models{'full', 'model'}{1}.fit_preval) / dt;
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
    
    if any(all(X==0,1)); disp('WARNING! Column with all zero predictors!!!'); end  % check for all-zero predictors
    
    options = struct('alpha', 0, 'lambda', lambdas, 'standardize', true);
    fit = cvglmnet(X(bins,:), y(bins), 'poisson', options, [], s.folds, ids, s.parallel, true);
    fit.fit_preval = fit.fit_preval(:,fit.lambda_min_id);  % only keep predictions for best lambda (save disk space)
end


end



