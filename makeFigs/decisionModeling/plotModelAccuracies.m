function plotModelAccuracies(data, varargin)

% train models to predict big vs. little step for different experimental
% conditions to see whether behavioral determinants are affected by
% manipulations // models are trained per mouse


% settings
s.condition = '';  % name of field in 'data' that contains different experimental conditions
s.levels = {''};  % levels of s.condition to plot
s.colors = [];

s.kFolds = 4;  % for k folds cross validation

s.successOnly = false;  % whether to only include successful trials
s.modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
s.lightOffOnly = false;  % whether to restrict to light on trials

s.saveLocation = '';  % if provided, save figure automatically to this location


% initialization
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if isempty(s.colors); s.colors = jet(length(s.levels)+1); end

flat = flattenData(data, {'mouse', 'session', 'trial', 'modPawOnlySwing', 'isTrialSuccess', ...
    'isBigStep', 'obsHgt', 'velAtWiskContact', 'wiskContactPosition', 'modPawKin', 'contactInd', 'isLightOn', s.condition});

% restrict to desired trials
if s.successOnly; flat = flat([flat.isTrialSuccess]==1); end
if s.lightOffOnly; flat = flat([flat.isLightOn]==0); end
if s.modPawOnlySwing; flat = flat([flat.modPawOnlySwing]==1); end


% compute x position of first mod paw at moment of whisker contact
for i = 1:length(flat)
    ind = max(flat(i).contactInd,1);
    ind = min(ind, size(flat(i).modPawKin,2));
    flat(i).modPawX = flat(i).modPawKin(1,ind);
end


% prepare predictor and target data
X = [[flat.modPawX]; [flat.obsHgt]; [flat.velAtWiskContact]; [flat.wiskContactPosition]]';
y = [flat.isBigStep];


% remove bins with NaNs
validBins = all(~isnan([X,y']), 2);
flat = flat(validBins);
X = X(validBins,:);
y = y(validBins);
y = logical(y);


mice = unique({flat.mouse});
[~, condition] = ismember({flat.(s.condition)}, s.levels);  % turn the 'condition' into numbers
models = cell(length(s.levels)+1, length(mice));  % one model per mouse per experimental condition, plus one per mouse for shuffled models
[accuracies, f1Scores] = deal(nan(length(s.levels)+1, length(mice)));


% loop over mice
for i = 1:length(mice)
    mouseBins = strcmp({flat.mouse}, mice{i});
    
    % models per condition
    for j = 1:length(s.levels)
        conditionBins = condition==j;
        X_sub = X(mouseBins & conditionBins,:);
        y_sub = y(mouseBins & conditionBins);
        if isempty(y_sub); keyboard; end
        [models{j,i}, accuracies(j,i), f1Scores(j,i)] = ...
            computeModel(X_sub, y_sub, s.kFolds);  % mouse model for this condition
    end
    
    % shuffled
    X_sub = X(mouseBins,:);
    y_sub = y(mouseBins);
    y_sub = y_sub(randperm(length(y_sub)));
    [models{end,i}, accuracies(end,i), f1Scores(end,i)] = ...
            computeModel(X_sub, y_sub, s.kFolds);
end




function [model, accuracy, f1] = computeModel(X, y, kFolds)
    % compute model accuracies and f1 score // accuracy and f1 score are
    % average of kFold partitions // model is created across all trials
    
    partitions = cvpartition(length(y), 'kfold', kFolds);  % cross validation splits
    [acc, f1s] = deal(nan(1, kFolds));
    
    for k = 1:kFolds
        model = fitglm(X(partitions.training(k),:), y(partitions.training(k)), 'Distribution', 'binomial');
        yhat = predict(model, X(partitions.test(k),:)) > .5;
        acc(k) = mean(y(partitions.test(k))==yhat');

        confusion = confusionmat(y(partitions.test(k))', yhat, 'Order', [false true]);
        precision = confusion(2,2)/sum(confusion(:,2));
        recall = confusion(2,2)/sum(confusion(2,:));
        f1s(k) = harmmean([precision, recall]);
%         if isnan(f1s(k)); keyboard; end
    end
    
    model = fitglm(X, y, 'Distribution', 'binomial');  % fit model on all data
    accuracy = nanmean(acc);
    f1 = nanmean(f1s);
end



% plot everything
figure('position', [2040.00 703.00 600 255.00], 'color', 'white', 'menubar', 'none')

% accuracies
subplot(1,2,1)
barFancy(accuracies, 'ylabel', 'model accuracy', 'levelNames', {[s.levels, 'shuffled']}, 'colors', s.colors)

% f1 scores
subplot(1,2,2)
barFancy(f1Scores, 'ylabel', 'f1 score', 'levelNames', {[s.levels, 'shuffled']}, 'colors', s.colors)


% save
if ~isempty(s.saveLocation); saveas(gcf, s.saveLocation, 'svg'); end

end

