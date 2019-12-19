function [X, y, predictorNames, isCategorical] = ...
    prepareDecisionModelData(flat, predictors, outcome, useAllPaws, referenceModPaw, normalizeData, varargin)

% given data struct flat containing predictors and outcome, creates X and y
% matrices that can be used to train models // useAllPaws is boolean
% determining whether to use all versions of variables that apply to each
% paw, or only the version corresponding to mod paw // if predictors is set
% to 'all', uses all available predcitors // otherwise, only uses those
% specified in cell array

% settings
allPredictors = {'velAtWiskContact', 'angleAtWiskContact', 'obsHgt', 'wiskContactPosition', ...
    'modPawStartingDistance', 'modStepStart', 'isStance', 'x', 'z', 'xVel', 'zVel', 'modPawPredictedDistanceToObs'};
s.removeNans = true; % whether to remove rows in X and y containing nans
s.balanceClasses = true; % whether to remove entries of modal y s.t. groups are perfectly balanced



% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
if strcmp(predictors, 'all'); predictors = allPredictors; end
cols = flat.Properties.VariableNames;
predictorNames = {};
for i = 1:length(predictors)

    % add to list if predictor appears in column as is
    if ismember(predictors{i}, cols)
        predictorNames{end+1} = predictors{i};

    % otherwise check if it is a paw specific variable
    elseif ismember([predictors{i} '_paw' num2str(referenceModPaw)], cols)

        if useAllPaws % add variable for each paw
            for j = 1:4; predictorNames{end+1} = [predictors{i} '_paw' num2str(j)]; end        
        else % add only for mod paw
            predictorNames{end+1} = [predictors{i} '_paw' num2str(referenceModPaw)];
        end
    end
    
end
isCategorical = contains(predictorNames, 'isStance');


% set up data
[~, colInds] = ismember(predictorNames, flat.Properties.VariableNames);
X = table2array(flat(:,colInds));
y = flat.(outcome);

% remove nans
if s.removeNans
    validBins = ~any(isnan([X y]),2);
    X = X(validBins,:);
    y = y(validBins);
end

% subsample to even out class sizes
if s.balanceClasses
    modalClass = mode(y);
    inds = find(y==modalClass);
    inds = inds(randperm(length(inds), sum(y~=modalClass)));
    inds = sort([inds; find(y~=modalClass)]);
    X = X(inds,:);
    y = y(inds);
end

if normalizeData
    Xtemp = X;
    X = normalize(X,1); % !!! should not normalize logical vars...
    X(:,isCategorical) = Xtemp(:,isCategorical); % make sure you don't z score logical / categorical vars!
end





