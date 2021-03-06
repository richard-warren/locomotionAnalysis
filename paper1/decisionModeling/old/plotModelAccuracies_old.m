function [accuracies, f1Scores] = plotModelAccuracies_old(flat, predictors, target, varargin)

% train models to predict big vs. little step for different experimental
% conditions to see whether behavioral determinants are affected by
% manipulations // models are trained per mouse // 'flat' is data struct with
% 'mouse' fields, and optional field for 'condition' // shold have fields
% for all 'predictors' listed, which will be used to construct predictors
% matrix X

% todo: add class balancing option?


% settings
s.condition = '';  % name of field in 'data' that contains different experimental conditions
s.levels = {''};  % levels of s.condition to plot
s.colors = [];
s.model = 'glm';  % 'glm' or 'ann'

s.kFolds = 4;  % for k folds cross validation
s.balanceClasses = false;  % whether to balance classes by subsampling
s.weightClasses = false;  % whether to balance classes by applying weights (only applies when model is 'glm')

s.successOnly = false;  % whether to only include successful trials
s.modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
s.lightOffOnly = false;  % whether to restrict to light on trials
s.deltaMin = 0;  % exclude little step trials where modPawDeltaLength is less than deltaLim standard deviations

s.barProps = {};  % properties to pass to barFancy
s.saveLocation = '';  % if provided, save figure automatically to this location

s.hiddenUnits = 100;  % if ann is used, this defines number of hidden units in 3 layers perceptron


% initialization
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if isempty(s.colors); s.colors = jet(length(s.levels)+1); end
if isstruct(flat); flat = struct2table(flat); end
if s.deltaMin; flat = flat(~(abs(zscore(flat.modPawDeltaLength))<s.deltaMin & flat.isBigStep==0), :); end


% restrict to desired trials
if s.successOnly; flat = flat(flat.isTrialSuccess==1, :); end
if s.lightOffOnly; flat = flat(flat.isLightOn==0, :); end
if s.modPawOnlySwing; flat = flat(flat.modPawOnlySwing==1, :); end


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
[accuracies, f1Scores] = deal(nan(length(s.levels), 2, length(mice)));  % (condition) X (shuffled vs. non-shuffled) X (mice)


% loop over mice
for i = 1:length(mice)
    mouseBins = strcmp(flat.mouse, mice{i});
    
    % models per condition
    for j = 1:length(s.levels)
        conditionBins = condition==j;
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
            [models{j,i}, accuracies(j,2,i), f1Scores(j,2,i)] = ...
                computeModel(X_sub, y_sub, s.kFolds);  % mouse model for this condition
            
            % train model on shuffled data
            [~, accuracies(j,1,i), f1Scores(j,1,i)] = ...
                computeModel(X_sub, y_sub(randperm(length(y_sub))), s.kFolds);  % mouse model for this condition
        end
    end
end




function [model, accuracy, f1] = computeModel(X, y, kFolds)
    % compute model accuracies and f1 score // accuracy and f1 score are
    % average of kFold partitions // model is created across all trials
    
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
                
                model = fitglm(X(partitions.training(k),:), y(partitions.training(k)), ...
                    'Distribution', 'binomial', 'Weights', weights(partitions.training(k)));
                yhat = predict(model, X(partitions.test(k),:)) > .5;
                
            case 'ann'
                net = patternnet(s.hiddenUnits);
                net = train(net, X(partitions.training(k),:)', y(partitions.training(k))');
                yhat = net(X(partitions.test(k),:)')' > .5;
        end
        
        % accuracy
        acc(k) = mean(y(partitions.test(k))==yhat);

        % f1 scores
        confusion = confusionmat(y(partitions.test(k))', yhat, 'Order', [false true]);
        precision = confusion(2,2)/sum(confusion(:,2));
        recall = confusion(2,2)/sum(confusion(2,:));
        f1s(k) = harmmean([precision, recall]);
    end
    
    model = fitglm(X, y, 'Distribution', 'binomial');  % fit model on all data
    accuracy = nanmean(acc);
    f1 = nanmean(f1s);
end


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
