%% load experiment data
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

% settings
modPawOnlySwing = true;  % if true, only include trials where the modified paw is the only one in swing
lightOffOnly = true;  % whether to restrict to light on trials
referenceModPaw = 2;  % flip trials around the y axis if first modified paw is not referenceModPaw
velTime = .05;  % how many samples to compute velocity over
useReferencePaw = true;  % if true, flip everything with respect to referenceModPaw
frameRate = 250;  % frame rate for videos


% initializations
mice = {data.data.mouse};
global_config;
velSmps = round(velTime * frameRate);
plotcol



%% compute predictors for behavioral model


% note: this code is a little strange, in that i flip predictors relative
% to referenceModPaw that are subsequently used by
% prepareDecisionModelData, but i don't flip other predictors
% (modPawKinInterp for example), which should be fine because i don't use
% these in the model, but still is a little strange... //  would probably
% be best to flip everything in one loop right off the bat, so everything
% downstream can inherit correctly flipped data...

flat = flattenData(data, {'mouse', 'session', 'trial', 'isTrialSuccess', ...
    'firstModPaw', 'velAtWiskContact', 'angleAtWiskContact', 'obsHgt', 'wiskContactPosition', 'isBigStep', ...
    'modPawKinInterp', 'preModPawKinInterp', 'isLightOn', 'modPawPredictedDistanceToObs', 'modPawDistanceToObs', 'modPawOnlySwing'});
flat = struct2table(flat);

% flip predictors in flat so everything is relative to firstModPaw
if useReferencePaw
    flipBins = [flat.firstModPaw]==referenceModPaw;
    flat.angleAtWiskContact(flipBins) = -flat.angleAtWiskContact(flipBins);
end


% add predictors not in flat already
sessions = unique(flat.session);

for i = 1:length(sessions)
    
    fprintf('%s: computing extra predictor variables\n', sessions{i})
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'kinData.mat'), 'kinData')
    sessionBins = strcmp(flat.session, sessions{i});
    if sum(sessionBins)~=length(kinData); disp('WTF! mismatch in number of trials! fucking shit!!!'); end
    
    % loop over paws
    for j = find([kinData.isTrialAnalyzed])
        
        % flip everything relative to first modified paw
        pawSequence = [1 2 3 4];
        if useReferencePaw && kinData(j).firstModPaw~=referenceModPaw; pawSequence = [4 3 2 1]; end
        flatBin = sessionBins & [flat.trial]==j;
        flat.firstModPaw(flatBin) = pawSequence(flat.firstModPaw(flatBin));
        
        for k = 1:4
            
            % starting paw position
            flat.(['modStepStart_paw' num2str(k)])(flatBin) = kinData(j).modifiedLocations{pawSequence(k)}(1,1,1);
            
            % is in stance
            flat.(['isStance_paw' num2str(k)])(flatBin) = kinData(j).stanceBins(kinData(j).contactInd,pawSequence(k));
            
            % x and z positions
            flat.(['x_paw' num2str(k)])(flatBin) = kinData(j).locations(kinData(j).contactInd, 1, pawSequence(k));
            flat.(['z_paw' num2str(k)])(flatBin) = kinData(j).locations(kinData(j).contactInd, 3, pawSequence(k));
            
            % x and z velocity
            kin = kinData(j).locations(kinData(j).contactInd-velSmps+1:kinData(j).contactInd, :, pawSequence(k));
            flat.(['xVel_paw' num2str(k)])(flatBin) = (kin(end,1)-kin(1,1)) / (velSmps/frameRate);
            flat.(['zVel_paw' num2str(k)])(flatBin) = (kin(end,3)-kin(1,3)) / (velSmps/frameRate);
        end
        
        % ind in modPawKinInterp at which whiskers contact obs
        flat.contactInd(flatBin) = kinData(j).pawObsPosIndInterp(kinData(j).firstModPaw);
    end
    
    % remove unanalyzed trials
    validBins = ~(sessionBins & ismember(flat.trial, find(~[kinData.isTrialAnalyzed]))); % all trial not in ssession OR in session but analyzed
    flat = flat(validBins,:);
end
disp('all done!')


% restrict to desired trials
flat_all = flat;
if lightOffOnly; flat = flat(logical(~[flat.isLightOn]),:); end
if modPawOnlySwing; flat = flat([flat.modPawOnlySwing]==1,:); end


%% prepare training data
[X, y, predictorNames, isCategorical] = ...  % less inclusive model
    prepareDecisionModelData(flat, ...
    {'velAtWiskContact', 'angleAtWiskContact', 'obsHgt', 'wiskContactPosition', 'x', 'xVel', 'z', 'zVel'}, ...
    'isBigStep', 'balanceClasses', true, 'useAllPaws', false, 'normalizeData', true, 'referencePaw', referenceModPaw);

% [X, y, predictorNames, isCategorical] = ...  % more inclusive model
%     prepareDecisionModelData(flat, ...
%     {'velAtWiskContact', 'angleAtWiskContact', 'obsHgt', 'wiskContactPosition', 'x', 'xVel', 'z', 'zVel'}, ...
%     'isBigStep', 'balanceClasses', true, 'useAllPaws', true, 'normalizeData', true, 'referencePaw', referenceModPaw);

% [X, y, predictorNames, isCategorical] = ...  % most inclusive model
%     prepareDecisionModelData(flat, ...
%     {'velAtWiskContact', 'angleAtWiskContact', 'obsHgt', 'wiskContactPosition', 'modStepStart', 'x', 'z', 'xVel', 'zVel'}, ...
%     'isBigStep', 'balanceClasses', true, 'useAllPaws', true, 'normalizeData', true, 'referencePaw', referenceModPaw);


% %% stepwise glm
% 
% model = stepwiseglm(X, y, 'linear', 'Distribution', 'binomial', ...
%     'PredictorVars', predictorNames, 'CategoricalVars', isCategorical);
% [~, sortInds] = sort(abs(glm.Coefficients.tStat(2:end)));


% %% lasso glm
% 
% % settings
% kFolds = 4;  % number of folds for k-folds cross validation
% 
% [B, fitInfo] = lassoglm(X, logical(y), 'binomial', ...
%     'CV', kFolds, 'PredictorNames', predictorNames, 'Standardize', false);
% % lassoPlot(B, fitInfo, 'PlotType', 'CV');  legend('show')
% % lassoPlot(B, fitInfo, 'PlotType', 'Lambda');  legend('show')
% 
% lambdaInd = fitInfo.IndexMinDeviance;
% % lambdaInd = fitInfo.Index1SE;
% [vals, sortInds] = sort(abs(B(:,lambdaInd)));
% fitInfo.PredictorNames(sortInds(vals~=0))
% yhat = X*B(:,lambdaInd)>0;
% accuracy = mean(yhat==y);
% fprintf('accuracy: %.2f\n', accuracy)


%% forward feature selection

% settings
kFolds = 10;

% function that determines misclassification rate of the model
% classf = @(xtrain, ytrain, xtest, ytest) ...
%     sum(ytest~=classify(xtest, xtrain, ytrain, 'linear'));
classf = @(xtrain, ytrain, xtest, ytest) ...
    sum(ytest ~= round(predict(fitglm(xtrain,ytrain,'linear','Distribution','binomial'), xtest)));



crossVals = cvpartition(size(X,1), 'kfold', kFolds);
[fs, fsHistory] = sequentialfs(classf, X, y, 'cv', crossVals, 'Nf', length(predictorNames));
% [fs, fsHistory] = sequentialfs(classf, X, y, 'cv', crossVals);

% figure out order in which features were added
numFeatures = size(fsHistory.In,1);
predictorColors = lines(numFeatures);
sortInds = nan(1,numFeatures);
sortInds(1) = find(fsHistory.In(1,:));
for i = 2:numFeatures; sortInds(i) = find(fsHistory.In(i,:) - fsHistory.In(i-1,:)); end

figure('Color', 'white', 'MenuBar', 'none', 'Position', [2375.00 747.00 297.00 232.00]); hold on
plot(1:numFeatures, 1-fsHistory.Crit, ...
    'LineWidth', 1, 'Color', [.2 .2 .2]);
scatter(1:numFeatures, 1-fsHistory.Crit, 60, predictorColors, 'filled');
xlabel('number of features')
ylabel('cross-validation accuracy')
set(gca, 'YLim', [.7 .9], 'YTick', .7:.1:.9, 'XLim', [1, numFeatures])
fprintf('max accuracy: %.2f\n', max(1-fsHistory.Crit))
fprintf('PREDICTORS: '); fprintf('%s ', predictorNames{sortInds}); fprintf('\n')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'baselineModelAccuracy');
saveas(gcf, file, 'svg');


% %% glm
% 
% glm = fitglm(X, y, 'linear', 'Distribution', 'binomial', ...
%     'PredictorVars', predictorNames, 'CategoricalVars', isCategorical);
% accuracy = mean(round(predict(glm,X))==y);
% fprintf('accuracy: %.2f\n', accuracy);
% [~, sortInds] = sort(abs(glm.Coefficients.tStat(2:end)));
% sortInds = flipud(sortInds);  % inds of best predictors from top to bottom
% glm

%% predictor scatters and histograms

% close all
figure('Color', 'white', 'MenuBar', 'none', 'Position', [1939.00 431.00 778.00 587.00]);
scatSz = 8;
scatAlpha = .4;
maxScatters = 1000;  % only plot this many per condition to avoid large vector images
percentileLims = [5 99];
binNum = 15;
maxPlots = 4;


sz = min(maxPlots, numFeatures);
randInds = randsample(size(X,1), min(maxScatters, size(X,1)));
X_sub = X(randInds,:);
y_sub = y(randInds,:);
lims = prctile(X, percentileLims, 1)';

for r = 1:sz
    for c = 1:sz
        
        subplot(sz, sz, (r-1)*sz + c); hold on
        r_X = sortInds(r);
        c_X = sortInds(c);
        
        % if along diagonal, plot histogram
        if r==c
            edges = linspace(lims(r_X,1), lims(r_X,2), binNum+1);
            histogram(X(~logical(y),r_X), edges, ...
                'FaceColor', decisionColors(1,:));
            histogram(X(logical(y),r_X), edges, ...
                'FaceColor', decisionColors(2,:))
        
        % otherwise scatter
        elseif r>c
            scatter(X_sub(~logical(y_sub),c_X), X_sub(~logical(y_sub),r_X), ...
                scatSz, decisionColors(1,:), 'filled', 'MarkerFaceAlpha', scatAlpha);
            scatter(X_sub(logical(y_sub),c_X), X_sub(logical(y_sub),r_X), ...
                scatSz, decisionColors(2,:), 'filled', 'MarkerFaceAlpha', scatAlpha);
            set(gca, 'XLim', lims(c_X,:), 'YLim', lims(r_X,:))
        else
            set(gca, 'Visible', 'off')
        end
        
        % pimp fig
        set(gca, 'XTick', [], 'YTick', [])
        if c==1
            ylabel(predictorNames{r_X}, 'Interpreter', 'none', 'FontWeight', 'bold', 'Color', predictorColors(r,:));
        end
        if r==sz
            xlabel(predictorNames{c_X}, 'Interpreter', 'none', 'FontWeight', 'bold', 'Color', predictorColors(c,:));
        end
    end
    pause(.001)
end

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'baselinePredictors');
saveas(gcf, file, 'svg');




% %% neural network vs. glm models
% 
% % settings
% iterations = 20;
% barColors = [0 .24 .49; .6 .75 .23; .2 .2 .2];
%     
%     
% accuracies = nan(3,iterations);
% for i = 1:iterations
%     
%     % NEURAL NET
%     net = patternnet(100);
%     net.divideParam.trainRatio = .7;
%     net.divideParam.valRatio = .15;
%     net.divideParam.testRatio = .15;
%     [net, tr] = train(net, X(:,~predDistBin)', y');
%     outputs = net(X(:,~predDistBin)');
%     accuracies(1,i) = mean(round(outputs(tr.testInd))==y(tr.testInd)');
%     
%     % NN SHUFFLED
%     shuffleInds = randperm(length(tr.testInd));
%     accuracies(3,i) = mean(round(outputs(tr.testInd))==y(tr.testInd(shuffleInds))'); % NN shuffled
%     
%     % GLM
% %     glmPredictorBins = true(1,size(X,2));  % bins of predictors to include in GLM (columns of X)
%     glmPredictorBins = predDistBin;
%     glm = fitglm(X([tr.trainInd tr.valInd], glmPredictorBins), y([tr.trainInd tr.valInd]), ...
%         'Distribution', 'binomial', 'CategoricalVars', isCategorical(glmPredictorBins));
%     accuracies(2,i) = mean(round(predict(glm, X(tr.testInd,glmPredictorBins)))==y(tr.testInd));
%     
% end
% 
% %
% figure('Color', 'white', 'MenuBar', 'none', 'Position', [1964 646 339 300]);
% barPlotRick(accuracies, {'conditionNames', {{'NN', 'GLM', 'shuffled'}}, ...
%     'lineThickness', 2, 'addBars', true, 'scatColors', 'hsv', 'scatAlpha', .5, 'showStats', false, ...
%     'ylim', [0 1], 'ytick', [0 .5 1], 'ylabel', 'accuracy', 'conditionColors', barColors})
% set(gca, 'Position', [0.15 0.2 0.75 0.79])
% 
% % save
% file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
%         'baseline_decision_modelAccuracies');
% saveas(gcf, file, 'svg');


%% schematic imgs for step type decision

session = '180630_001';
yMax = [];  % set to 140 to trim the bottom view out

close all
showDecisionFrames(session, 'stepColors', decisionColors, 'yMax', yMax, ...
    'contrastLims', contrast, 'leftPadding', 0, 'rightPadding', 70, 'drawKinematics', false);



%% predictor scatters


figure('Color', 'white', 'Position', [2000 400 450 350], 'MenuBar', 'none');
scatterHistoRick([flat.x_paw2]*1000, [flat.trialVel]*1000, ...
    {'groupId', [flat.isBigStep]+1, 'colors', decisionColors, ...
    'xlabel', 'paw position (mm)', 'ylabel', 'obstalce height (mm)', ...
    'showCrossHairs', false, 'scatSize', 4, 'scatAlpha', .4, ...
    'xLims', [-70 -15], 'yLims', [3 11], 'groupHistoLineWidth', 3, ...
    'groupHistoAlpha', .8, 'histoAlpha', 0});

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baseline_decision_PositionObsHgtScatter');
saveas(gcf, file, 'svg');


%% binned kinematics

% settings
binNum = 5;
pctileBins = true;


% initializations
kin = permute(cat(3,flat.modPawKinInterp{:}), [3,1,2]);
kinCtl = permute(cat(3,flat.preModPawKinInterp{:}), [3,1,2]);
[X, y, predictorNames, isCategorical] = ...
        prepareDecisionModelData(flat, 'all', 'isBigStep', true, referenceModPaw, normalizeData, ...
        {'balanceClasses', false, 'removeNans', false});  % must recom 

% prepareDecisionModelData(flat, 'all', 'isBigStep', true, referenceModPaw, normalizeData, {'balanceClasses', true});


% choose binning variable
binVar = net(X(:,~predDistBin)'); % neural network output
% binVar = -flat.x_paw2; % position of mod paw at moment of contact
% binVar = flat.velAtWiskContact; % vel
% binVar = flat.modStepStart_paw2; % starting position of mod paw
% binVar = flat.modPawPredictedDistanceToObs;


% perctentile bins    
if pctileBins
    binEdges = prctile(binVar, linspace(0, 100, binNum+1));
% evenly spaced bins
else
    binEdges = linspace(min(binVar), max(binVar), binNum+1);
end
condition = discretize(binVar, binEdges);


figure('color', 'white', 'menubar', 'none', 'position', [2000 100 750 100*binNum], 'InvertHardcopy', 'off');
plotBigStepKin(kin(:,[1,3],:), kinCtl(:,[1,3],:), flat.obsHgt, condition, flat.isBigStep, ...
    {'colors', stepTypeColors, 'xLims', [-.09 .04], 'addHistos', false, 'lineWid', 3, ...
    'contactInds', flat.contactInd, 'histoHgt', .5, 'showSmpNum', false})

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baseline_decision_kinematics');
saveas(gcf, file, 'svg');


%% heatmap

% settings
xLims = [-.03 .015];
yLims = [-.03 .03];

figure('Color', 'white', 'Position', [2006 540 384 395], 'MenuBar', 'none');
heatmapRick(flat.modPawPredictedDistanceToObs, flat.modPawDistanceToObs, ...
    {'xLims', xLims, 'yLims', yLims, ...
    'xlabel', 'predicted distance to obstalce (m)', 'ylabel', 'distance to obstalce (m)'})
set(gca, 'DataAspectRatio', [1 1 1])
line(xLims, xLims, 'color', [.5 .5 1 .8], 'lineWidth', 3)

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baseline_decision_heatmaps');
saveas(gcf, file, 'svg');



