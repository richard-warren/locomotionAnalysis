function [models, fitdata] = fitNeuronGlm(session, neuron, varargin)

% fit (several) glms for a neuron: full model, single model for each
% predictor group, and models with each predictor group removed // either
% fits single model on all data (method 'none'), or additional assesses
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
s.outputFileName = '';               % whether to save output to file
s.verbose = true;




% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
if s.verbose; fprintf('%s: fitting models for neuron %i... ', session, neuron); end

% load design matrix
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'designMatrices', [session '_designMatrix.mat']), ...
    'dmat', 't', 'reward_all')

% load predictor info
predictorInfo = readtable(fullfile(getenv('GITDIR'), 'locomotionAnalysis', 'paper2', 'glm', 'predictorSettings.csv'));
groups = unique(predictorInfo.group);

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
X = table2array(dmat);
model = fitModel(X, y, s.lambdas);
models{'full', 'model_in'}{1} = model;
dev_full = cvdeviance(X, y, model, 'holdout', true, 'bestLambdaOnly', true);
models{'full', 'dev_in'} = dev_full;
models{'full', 'dev_out'} = 0;
lambda_min = model.lambda_min;

% check that rick's deviance matches glmnet deviance for full model evaluated on training data
glmnet_dev = model.glmnet_fit.dev(model.lambda_min_id);
rick_dev = cvdeviance(X, y, model, 'holdout', true, 'bestLambdaOnly', true);
if abs(glmnet_dev - rick_dev)>.01
    disp('WARNING! Deviance computed by Glmnet disagrees with rick deviance!')
end


% assess importance of groups of features
for i = 1:length(groups)
    members = predictorInfo.name(strcmp(predictorInfo.group, groups{i}));  % members of group
    
    % single group model
    switch s.method
        case 'shuffle'
            colBins = ~ismember(dmat.Properties.VariableNames, [members; 'time']);
            X = dmat;
            for j = find(colBins)'; X(:,j) = X(randperm(height(dmat)), j); end  % shuffle
            X = table2array(X);
            models{groups{i}, 'dev_in'} = cvdeviance(X, y, model, 'holdout', true, 'bestLambdaOnly', true);
        
        case 'mask'
            colBins = ~ismember(dmat.Properties.VariableNames, [members; 'time']);
            X = dmat;
            for j = find(colBins)'; X{:,j} = zeros(size(X{:,j})); end  % mask with zeros
            X = table2array(X);
            models{groups{i}, 'dev_in'} = cvdeviance(X, y, model, 'holdout', true, 'bestLambdaOnly', true);
            
        case 'refit'
            colBins = ismember(dmat.Properties.VariableNames, members);
            X = dmat(:, colBins);
            X = table2array([X dmat(:,'time')]);
            model = fitModel(X, y, lambda_min);
            models{groups{i}, 'model_in'}{1} = model;
            models{groups{i}, 'dev_in'} = cvdeviance(X, y, model, 'holdout', true, 'bestLambdaOnly', true);
    end
    
    % single group removed
    switch s.method
        case 'shuffle'
            colBins = ismember(dmat.Properties.VariableNames, members);
            X = dmat;
            for j = find(colBins)'; X(:,j) = X(randperm(height(dmat)), j); end  % shuffle
            X = table2array(X);
            models{groups{i}, 'dev_out'} = dev_full - cvdeviance(X, y, model, 'holdout', true, 'bestLambdaOnly', true);
        
        case 'mask'
            colBins = ismember(dmat.Properties.VariableNames, members);
            X = dmat;
            for j = find(colBins)'; X{:,j} = zeros(size(X{:,j})); end  % mask with zeros
            X = table2array(X);
            models{groups{i}, 'dev_out'} = dev_full - cvdeviance(X, y, model, 'holdout', true, 'bestLambdaOnly', true);
            
        case 'refit'
            X = dmat(:, ~ismember(dmat.Properties.VariableNames, members));
            X = table2array(X);
            model = fitModel(X, y, lambda_min);
            models{groups{i}, 'model_out'}{1} = model;
            models{groups{i}, 'dev_out'} = dev_full - cvdeviance(X, y, model, 'holdout', true, 'bestLambdaOnly', true);
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
if ~isempty(s.outputFileName); save(s.outputFileName, 'models', 'fitdata'); end


function fit = fitModel(X, y, lambdas)
    % fit model with k fold cross validations
    if any(isnan(X(:))); keyboard; end
    if length(lambdas)==1; lambdas = [0 lambdas]; end  % a hack, because cvglmnet requires multiple lambdas
    options = struct('alpha', 0, 'lambda', lambdas, 'standardize', true);
    fit = cvglmnet(X, y, 'poisson', options, [], s.folds, fold_id, s.parallel, true);
    fit.fit_preval = fit.fit_preval(:,fit.lambda_min_id);  % only keep predictions for best lambda (save disk space)
end


end



