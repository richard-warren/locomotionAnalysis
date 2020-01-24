function plotModelAccuracies(data, varargin)

% train models to predict big vs. little step for different experimental
% conditions to see whether behavioral determinants are affected by
% manipulations // models are trained per mouse


% settings
s.kFolds = 5;  % for k folds cross validation

s.successOnly = false;  % whether to only include successful trials
s.modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
s.lightOffOnly = false;  % whether to restrict to light on trials


% initialization
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin

flat = flattenData(data, {'mouse', 'session', 'trial', 'conditionNum', 'modPawOnlySwing', 'condition', 'isTrialSuccess', ...
    'isBigStep', 'obsHgt', 'velAtWiskContact', 'wiskContactPosition', 'modPawKin', 'contactInd', 'isLightOn', 'modPawPredictedDistanceToObs'});

% restrict to desired trials
if s.lightOffOnly; flat = flat([flat.isLightOn]==0); end
if s.modPawOnlySwing; flat = flat([flat.modPawOnlySwing]==1); end
if s.modPawOnlySwing; flat = flat([flat.isTrialSuccess]==1); end

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


mice = unique({flat.mouse});

for i = 1:length(mice)
    
    mouseBins = strcmp({flat.mouse}, mice{i});
    models = cell(1,2);
    
    for j = 1:2
        X_sub = X(mouseBins & allBins{j},:);
        y_sub = y(mouseBins & allBins{j});
        crossVals = cvpartition(length(y_sub), 'kfold', kFolds);  % cross validation splits
        [mouseAccuracies, mouseF1s] = deal(nan(1, kFolds));

        for k = 1:kFolds
            glm = fitglm(X_sub(crossVals.training(k),:), y_sub(crossVals.training(k)), 'Distribution', 'binomial');
            predictions = round(predict(glm, X_sub(crossVals.test(k),:)));
            mouseAccuracies(k) = mean(y_sub(crossVals.test(k))==predictions');
            
            confusion = confusionmat(y_sub(crossVals.test(k)), logical(predictions));
            precision = confusion(2,2)/sum(confusion(:,2));
            recall = confusion(2,2)/sum(confusion(2,:));
            mouseF1s(k) = harmmean([precision, recall]);
        end
        models{j} = fitglm(X_sub, y_sub, 'Distribution', 'binomial');  % mouse model for this condition
        accuracies(j,i) = mean(mouseAccuracies);
        f1Scores(j,i) = mean(mouseF1s);
    end
    
    
    % pre model evaluated on post data
    X_sub = X(mouseBins & postBins,:);
    y_sub = y(mouseBins & postBins);
    predictions = round(predict(models{1}, X_sub));
    accuracies(3,i) = mean(predictions == y_sub'); 
    
    confusion = confusionmat(y_sub, logical(predictions));
    precision = confusion(2,2)/sum(confusion(:,2));
    recall = confusion(2,2)/sum(confusion(2,:));
    f1Scores(3,i) = harmmean([precision, recall]);
    
    
    % shuffled
    X_sub = X(mouseBins,:);
    y_sub = y(mouseBins);
    crossVals = cvpartition(length(y_sub), 'kfold', kFolds);  % corss validation splits
    mouseAccuracies = nan(1, kFolds);
    
    for j = 1:kFolds
        glm = fitglm(X_sub(crossVals.training(j),:), y_sub(crossVals.training(j)), 'Distribution', 'binomial');
        predictions = round(predict(glm, X_sub(crossVals.test(j),:)));
        predictions = predictions(randperm(length(predictions)));  % shuffle
        mouseAccuracies(j) = mean(y_sub(crossVals.test(j))==predictions');
        
        confusion = confusionmat(y_sub(crossVals.test(j)), logical(predictions));
        precision = confusion(2,2)/sum(confusion(:,2));
        recall = confusion(2,2)/sum(confusion(2,:));
        f1Scores(4,i) = harmmean([precision, recall]);
    end
    accuracies(4,i) = mean(mouseAccuracies);
end

if strcmp(dataset, 'mtc_muscimol')
    colors_temp = [ctlStepColor; muscimolColor; mean([ctlStepColor; muscimolColor],1); ctlStepColor*.5];
else
    colors_temp = [ctlStepColor; lesionColor; mean([ctlStepColor; lesionColor],1); ctlStepColor*.5];
end

% accuracies
figure('position', [2018.00 100.00 252.00 239.00], 'color', 'white', 'menubar', 'none')
barFancy(accuracies, 'ylabel', 'model accuracy', 'levelNames', {{'pre', 'post', 'pre->post', 'shuffled'}}, 'colors', colors_temp, barProperties{:})
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_decisionModels' suffix1 suffix2]), 'svg');

% f1 scores
figure('position', [2018.00 100.00 252.00 239.00], 'color', 'white', 'menubar', 'none')
barFancy(f1Scores, 'ylabel', 'f1 score', 'levelNames', {{'pre', 'post', 'pre->post', 'shuffled'}}, 'colors', colors_temp, barProperties{:})
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_decisionModelsF1' suffix1 suffix2]), 'svg');



% % compare decision thresholds
% thresholds = nan(2, length(mice));  % (pre/post) X (mouse)
% for i = 1:length(mice)
%     mouseBins = strcmp({flat.mouse}, mice{i});
%     preBins = strcmp({flat.condition}, matchConditions(1));
%     postBins = strcmp({flat.condition}, matchConditions(2));
%     allBins = {preBins, postBins};
%     
%     for j = 1:2
%         x = [flat(mouseBins & allBins{j}).modPawPredictedDistanceToObs] * 1000;
%         y = [flat(mouseBins & allBins{j}).isBigStep];
%         
%         glm = fitglm(x', y', 'Distribution', 'binomial');
%         coeffs = glm.Coefficients.Estimate;
%         thresholds(j,i) = (-coeffs(1)) / coeffs(2); % solve for prediction = 0
%     end
% end

figure('position', [2018.00 100.00 252.00 239.00], 'color', 'white', 'menubar', 'none')
barFancy(thresholds, 'ylabel', 'decision threshold', 'levelNames', {vars.condition.levelNames}, ...
    'colors', colors, barProperties{:}, 'textRotation', 0)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'manipulations', [dataset '_decisionThresholds' suffix1 suffix2]), 'svg');

% show histograms to sanity check decision thresholds
% figure;
% bins = -50:1:50;
% x = [flat.modPawPredictedDistanceToObs] * 1000;
% y = [flat.isBigStep];
% subplot(2,1,1);
% histogram(x(y & preBins), bins); hold on; histogram(x(~y & preBins), bins)
% subplot(2,1,2);
% histogram(x(y & postBins), bins); hold on; histogram(x(~y & postBins), bins)
