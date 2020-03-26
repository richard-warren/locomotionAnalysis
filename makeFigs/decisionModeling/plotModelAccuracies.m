function [accuracies, f1Scores, flat] = plotModelAccuracies(flat, predictors, target, varargin)

% train models to predict big vs. little step for different experimental
% conditions to see whether behavioral determinants are affected by
% manipulations // models are trained per mouse // 'flat' is data struct with
% 'mouse' fields, and optional field for 'condition' // shold have fields
% for all 'predictors' listed, which will be used to construct predictors
% matrix X // flat_sub is a version of flat with rows removed that violate
% trial restrictions, and with a [s.outcome '_predicted'] column containing
% the predictions for each row

% todo: add class balancing option?


% settings
s.condition = '';  % name of field in 'data' that contains different experimental conditions
s.levels = {''};  % levels of s.condition to plot
s.colors = [];
s.model = 'glm';  % 'glm' or 'ann'
s.modelTransfers = [];  % nX2 matrix describing which models trained under one condition (left column) to test on another condition (right column)

s.kFolds = 10;  % for k folds cross validation
s.balanceClasses = false;  % whether to balance classes by subsampling
s.weightClasses = false;  % whether to balance classes by applying weights (only applies when model is 'glm')

s.successOnly = false;  % whether to only include successful trials
s.modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
s.lightOffOnly = false;  % whether to restrict to light on trials
s.deltaMin = 0;  % exclude little step trials where modPawDeltaLength is less than deltaLim standard deviations

s.plot = true;  % whether to generate plot
s.barProps = {};  % properties to pass to barFancy
s.saveLocation = '';  % if provided, save figure automatically to this location

s.hiddenUnits = 100;  % if ann is used, this defines number of hidden units in 3 layers perceptron


% initialization
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if isempty(s.colors); s.colors = jet(length(s.levels)); end
if isstruct(flat); flat = struct2table(flat); end
cNum = length(s.levels) + size(s.modelTransfers,1);  % total number of conditions


% restrict to desired trials
if s.successOnly; flat = flat(flat.isTrialSuccess==1, :); end
if s.lightOffOnly; flat = flat(flat.isLightOn==0, :); end
if s.modPawOnlySwing; flat = flat(flat.modPawOnlySwing==1, :); end
if s.deltaMin; flat = flat(~(abs(zscore(flat.modPawDeltaLength))<s.deltaMin & flat.isBigStep==0), :); end


% prepare predictor and target
[~, predictorInds] = ismember(predictors, flat.Properties.VariableNames);
X = table2array(flat(:, predictorInds));
y = flat.(target);


% remove bins with NaNs
validBins = all(~isnan([X,y]), 2);
flat = flat(validBins,:);
X = X(validBins,:);
y = y(validBins);
y = logical(y);


mice = unique(flat.mouse);
if ~isempty(s.condition)
    [~, condition] = ismember(flat.(s.condition), s.levels);  % turn the 'condition' into numbers
else
    condition = ones(height(flat), 1);
end
models = cell(length(s.levels), length(mice));  % (condition) X (mice)
[accuracies, f1Scores] = deal(nan(cNum, 2, length(mice)));  % (condition) X (shuffled vs. non-shuffled) X (mice)
flat.([target '_predicted']) = nan(height(flat),1);


% loop over mice
for i = 1:length(mice)
    mouseBins = strcmp(flat.mouse, mice{i});
    
    % models per condition
    for j = 1:cNum
        
        % for normal conditions
        if j<=length(s.levels)  
            conditionBins = condition==j;
            model = {};
            
        % when transfering model from one condition to another
        else  
            conditionBins = condition==s.modelTransfers(j-length(s.levels),2);
            model = {models{s.modelTransfers(j-length(s.levels),1), i}};
        end
        
        X_sub = X(mouseBins & conditionBins,:);
        y_sub = y(mouseBins & conditionBins);
        
        if s.balanceClasses
            n = min(sum(y_sub), sum(~y_sub));
            
            inds_t = find(y_sub);
            inds_t = inds_t(randperm(length(inds_t), n));
            inds_f = find(~y_sub);
            inds_f = inds_f(randperm(length(inds_f), n));
            
            inds = sort([inds_t; inds_f]);
            X_sub = X_sub(inds,:);
            y_sub = y_sub(inds);
        end
        
        if ~isempty(y_sub)
            
            % train model
            [models{j,i}, accuracies(j,2,i), f1Scores(j,2,i), predictions] = ...
                computeModel(X_sub, y_sub, s.kFolds, model{:});  % mouse model for this condition
            
            % store model predictions
            flat.([target '_predicted'])(mouseBins & conditionBins) = predictions;
            
            % train model on shuffled data
            [~, accuracies(j,1,i), f1Scores(j,1,i)] = ...
                computeModel(X_sub, y_sub(randperm(length(y_sub))), s.kFolds, model{:});  % mouse model for this condition
        end
    end
end




function [model, accuracy, f1, predictions] = computeModel(X, y, kFolds, prevModel)
    % compute model accuracies and f1 score // accuracy and f1 score are
    % average of kFold partitions // model is created across all trials
    
%     disp(size(X,1));
    partitions = cvpartition(length(y), 'kfold', kFolds);  % cross validation splits
    [acc, f1s] = deal(nan(1, kFolds));
    
    for k = 1:kFolds
        
        % train model
        switch s.model
            case 'glm'
                
                if s.weightClasses
                    w_f = sum(y) / length(y);
                    weights = ones(length(y), 1) * w_f;
                    weights(y) = 1 - w_f;
                else
                    weights = ones(length(y), 1);
                end
                
                if exist('prevModel', 'var')
                    model = prevModel;
                else
                    try
                    model = fitglm(X(partitions.training(k),:), y(partitions.training(k)), ...
                        'Distribution', 'binomial', 'Weights', weights(partitions.training(k)));
                    catch; keyboard; end
                end
                yhat = predict(model, X(partitions.test(k),:)) > .5;
                
            case 'ann'
                if exist('prevModel', 'var')
                    model = prevModel;
                else 
                    model = patternnet(s.hiddenUnits);
                    model = train(model, X(partitions.training(k),:)', y(partitions.training(k))');
                end
                yhat = model(X(partitions.test(k),:)')' > .5;
        end
        
        % accuracy
        acc(k) = mean(y(partitions.test(k))==yhat);

        % f1 scores
        confusion = confusionmat(y(partitions.test(k))', yhat, 'Order', [false true]);
        precision = confusion(2,2)/sum(confusion(:,2));
        recall = confusion(2,2)/sum(confusion(2,:));
        f1s(k) = harmmean([precision, recall]);
    end
    
    % compuate accuracy and f1 across partitions
    accuracy = nanmean(acc);
    f1 = nanmean(f1s);
    
    % train model on on samples
    model = fitglm(X, y, 'Distribution', 'binomial');  % fit model on all data
    predictions = predict(model, X);
end


% add model transfers to condition names and colors
for i = 1:size(s.modelTransfers,1)
    name = [s.levels{s.modelTransfers(i,1)} ' -> ' s.levels{s.modelTransfers(i,2)}];
    s.levels{end+1} = name;
    s.colors(end+1,:) = s.colors(s.modelTransfers(i,1),:);
end

if s.plot
    % plot everything
    figure('position', [2040.00 703.00 600 255.00], 'color', 'white', 'menubar', 'none')
    colors = repelem(s.colors, 2, 1);
    colors = colors .* repmat([.5;1],length(s.levels),1);

    % accuracies
    subplot(1,2,1)
    barFancy(accuracies, 'ylabel', 'model accuracy', 'levelNames', {s.levels}, 'colors', colors, s.barProps{:})

    % f1 scores
    subplot(1,2,2)
    barFancy(f1Scores, 'ylabel', 'f1 score', 'levelNames', {[s.levels, 'shuffled']}, 'colors', colors, s.barProps{:})


    % save
    if ~isempty(s.saveLocation); saveas(gcf, s.saveLocation, 'svg'); end
end

flat = table2struct(flat)';

end
