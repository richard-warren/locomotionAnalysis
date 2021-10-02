function [models, fitdata] = fitUpperLowerGlm(session, neuron, varargin)

% fit (several) glms for a neuron: full model, single model for each
% predictor group, and models with each predictor group removed // either
% fits single model on all data (method 'none'), or additionally assesses
% importance of each feature group by making models with just that group
% included, and just that group excluded // either does this by refitting
% models for with and without each group (method 'refit') or by shuffling
% each group and all predictors other than that group in time (method
% 'shuffle') // all models use ridge regression, with lambda automatically
% chosen to maximize test set performance // lambda is chosen on the full
% model and applied to nested models

% for refitting, best lambda determined by full model and applied to other
% models, although 0 reg is used if it's generalization is better... // 
% for shuffling, best lambda used for each nested model...


% settings
s.lambdas = logspace(-8, -1, 40);    % ridge regression coefficients
s.folds = 5;                         % cross-validation folds
s.parallel = false;                  % whether crossval analyses are parallelized
s.method = 'refit';                  % 'shuffle', or 'refit', or 'mask' (see description above...)
s.verbose = true;
s.save = true;
s.outputFileName = fullfile(...
    getenv('SSD'), 'paper2', 'modelling', 'glms', 'upperlower_glms', ...
    [session '_cell_' num2str(neuron) '_glm.mat']);




% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
if s.verbose; fprintf('%s: fitting models for neuron %i... ', session, neuron); end

% load design matrix
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'designMatrices', 'upperLower', [session '_designMatrix.mat']), ...
    'dmat', 't', 'reward_all')

% load predictor info
filename = fullfile(getenv('GITDIR'), 'locomotionAnalysis', 'paper2', 'glm', 'settings', 'upperlower_predictorSettings.xlsx');
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
nrows = length(groups) + 1;
rowNames = ['full', groups'];
models = table(cell(nrows,1), cell(nrows,1), nan(nrows,1), nan(nrows,1), 'RowNames', rowNames, ...
    'VariableNames', {'model_in', 'model_out', 'dev_in', 'dev_out'});

% train full model
X_full = table2array(dmat);
model = fitModel(X_full, y, s.lambdas);
models{'full', 'model_in'}{1} = model;
dev_full = cvdeviance(X_full, y, model, 'holdout', true, 'bestLambdaOnly', true);
models{'full', 'dev_in'} = dev_full;
models{'full', 'dev_out'} = 0;
lambda_min = model.lambda_min;

% check that rick's deviance matches glmnet deviance for full model evaluated on training data
glmnet_dev = model.glmnet_fit.dev(model.lambda_min_id);
rick_dev = cvdeviance(X_full, y, model, 'holdout', false, 'bestLambdaOnly', true);  % don't use heldout data to match glmnet_fit deviance
if abs(glmnet_dev - rick_dev)>.02
    fprintf('WARNING! Deviance computed by Glmnet %.3f off from rick deviance!\n', ...
        glmnet_dev - rick_dev)
end

% get names of columns in X
colNames = cell(1,length(dmat.Properties.VariableNames));
for i = 1:length(dmat.Properties.VariableNames)
    colNames{i} = repelem({dmat.Properties.VariableNames{i}}, length(dmat{1,i}));  % number of 
end
colNames = cat(2, colNames{:});


% assess importance of groups of features
for i = 1:length(groups)
    members = predictorInfo.Properties.RowNames(strcmp(predictorInfo.group, groups{i}));  % members of group
    
    % determine time bins relevant to single group models
    % deviance will only be computed at these times
    restrictTime = ~isnan(groupInfo{groups{i}, 'restrict_time'}{1});
    
    if restrictTime
        colBins = ismember(colNames, groupInfo{groups{i}, 'restrict_time'});  % only restrict time based on predictor listed in restrict_time column
        timeBins = any(X_full(:,colBins)>0,2);
        dev_full_group = cvdeviance(X_full, y, models{'full', 'model_in'}{1}, ...
            'holdout', true, 'bestLambdaOnly', true, 'timeBins', timeBins);
    else
        timeBins = true(size(X_full,1),1);
        dev_full_group = dev_full;
    end
    
    
    % single group model
    switch s.method
        case 'shuffle'
            X = X_full;
            colBins = ~ismember(colNames, [members; 'time']);
            for j = find(colBins)'; X(:,j) = X(randperm(size(X,1)), j); end  % shuffle
            models{groups{i}, 'dev_in'} = cvdeviance(X, y, model, ...
                'holdout', true, 'bestLambdaOnly', true, 'timeBins', timeBins);
        
        case 'mask'
            colBins = ismember(colNames, [members; 'time']);
            X = X_full .* colBins;
            models{groups{i}, 'dev_in'} = cvdeviance(X, y, model, ...
                'holdout', true, 'bestLambdaOnly', true, 'timeBins', timeBins);
            
        case 'refit'
            colBins = ismember(colNames, [members; 'time']);
            X = X_full(:, colBins);
            model = fitModel(X, y, lambda_min);
            models{groups{i}, 'model_in'}{1} = model;
            models{groups{i}, 'dev_in'} = cvdeviance(X, y, model, ...
                'holdout', true, 'bestLambdaOnly', true, 'timeBins', timeBins);
    end
    
    
    % single group removed
    switch s.method
        case 'shuffle'
            X = X_full;
            colBins = ismember(colNames, members);
            for j = find(colBins)'; X(:,j) = X(randperm(size(X,1)), j); end  % shuffle
            models{groups{i}, 'dev_out'} = dev_full_group - ...
                cvdeviance(X, y, model, ...
                'holdout', true, 'bestLambdaOnly', true, 'timeBins', timeBins);
        
        case 'mask'
            colBins = ~ismember(colNames, members);
            X = X_full .* colBins;
            models{groups{i}, 'dev_out'} = dev_full_group - ...
                cvdeviance(X, y, model, ...
                'holdout', true, 'bestLambdaOnly', true, 'timeBins', timeBins);
            
        case 'refit'
            colBins = ~ismember(colNames, members);
            X = X_full(:, colBins);
            model = fitModel(X, y, lambda_min);
            models{groups{i}, 'model_out'}{1} = model;
            models{groups{i}, 'dev_out'} = dev_full_group - ...
                cvdeviance(X, y, model, ...
                'holdout', true, 'bestLambdaOnly', true, 'timeBins', timeBins);
    end
end

if s.verbose; disp('all done!'); end

% save some objects for convenience...
yhat = exp(models{'full', 'model_in'}{1}.fit_preval) / dt;
fitdata.t = t;
fitdata.y = y;
fitdata.yRate = spkRate;
fitdata.yhat = yhat;
fitdata.groups = groups;
fitdata.session = session;
fitdata.neuron = neuron;

% save
if s.save; save(s.outputFileName, 'models', 'fitdata'); end


function fit = fitModel(X, y, lambdas)
    % fit model with k fold cross validations
    if length(lambdas)==1; lambdas = [0 lambdas]; end  % a hack, because cvglmnet requires multiple lambdas
    options = struct('alpha', 0, 'lambda', lambdas, 'standardize', true);
    fit = cvglmnet(X, y, 'poisson', options, [], s.folds, fold_id, s.parallel, true);
    fit.fit_preval = fit.fit_preval(:,fit.lambda_min_id);  % only keep predictions for best lambda (save disk space)
end


end



