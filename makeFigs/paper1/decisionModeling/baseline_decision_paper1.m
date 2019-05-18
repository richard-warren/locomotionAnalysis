%% COMPUTE PREDICTORS

% global settings
referenceModPaw = 2;
velSmps = 10; % how many samples to compute velocity over
frameRate = 250;
stepTypeColors = [0.850 0.325 0.098; 0 0.447 0.741]; % first entry is small step, second is big step
useReferencePaw = true; % if true, flip everything with respect to referenceModPaw

% load experiment data
fprintf('loading...'); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

% flatten data and compute predictors
flat = flattenData(data, {'mouse', 'session', 'trial', 'isTrialSuccess', ...
    'firstModPaw', 'velAtWiskContact', 'angleAtWiskContact', 'obsHgt', 'wiskContactPosition', 'isBigStep', ...
    'modPawKinInterp', 'preModPawKinInterp', 'isLightOn', 'modPawPredictedDistanceToObs', 'modPawDistanceToObs'});
flat = struct2table(flat);

% flip predictors in flat so everything is relative to firstModPaw
flipBins = [flat.firstModPaw]==referenceModPaw;
flat.angleAtWiskContact(flipBins) = -flat.angleAtWiskContact(flipBins);


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
            flat.(['xVel_paw' num2str(k)])(flatBin) = (kin(end,1)-kin(1,1)) / velSmps / frameRate;
            flat.(['zVel_paw' num2str(k)])(flatBin) = (kin(end,3)-kin(1,3)) / velSmps / frameRate;
        end
        
        % ind in modPawKinInterp at which whiskers contact obs
        flat.contactInd(flatBin) = kinData(j).pawObsPosIndInterp(kinData(j).firstModPaw);
    end
    
    % remove unanalyzed trials
    validBins = ~(sessionBins & ismember(flat.trial, find(~[kinData.isTrialAnalyzed]))); % all trial not in ssession OR in session but analyzed
    flat = flat(validBins,:);
end
disp('all done!')

% keep copy of data prior to applying restrictions below
flatAll = flat;




% %% BUILD MODELS
% 
% % NOTES: important considerations: do we restrict to light off // do we
% % restrict to trials where only mod paw is in the air // in theory neural
% % net should not need things flipped wrt mod paw, bc mod paw can be
% % inferred from stance positions... // should see if this is true
% 
% % apply restrictions to dataset
% flat = flatAll;
% flat = flat(~flat.isLightOn,:);
% flat = flat((flat.isStance_paw2+flat.isStance_paw3)==1, :); % only trials where only mod paw is in air at contact!
% 
% 
% 
% % settings
% useAllPaws = [true false false];
% predictors = {{'all'}, {'all'}, {'x_paw2', 'obsHgt', 'velAtWiskContact'}};
% outcome = 'isBigStep';
% iterations = 20;
% normalizeData = true;
% 
% 
% accuracies = nan(length(useAllPaws), 2, iterations);
% [nnets, glms] = deal(cell(1,length(useAllPaws)));
% 
% for i = 1:length(useAllPaws)
%     
%     [X, y, ~, isCategorical] = ...
%         prepareDecisionModelData(flat, predictors{i}, outcome, useAllPaws(i), referenceModPaw, normalizeData, ...
%         {'balanceClasses', true});
%     
%     [accuracies(i,:,:), nnets{i}, glms{i}] = trainDecisionModels(X, y, iterations, isCategorical);
%     
% end
% 
% 
% % PLOT MODEL ACCURACY
% 
% % figure('Color', 'white', 'MenuBar', 'none', 'Position', [1964 646 600 300]);
% % barPlotRick(accuracies, {'lineThickness', 2, 'addBars', true, 'scatColors', 'hsv', 'scatAlpha', .5, 'showStats', false, ...
% %     'ylim', [.5 1], 'ytick', [.5 .75 1], 'ylabel', 'accuracy'})
% % set(gca, 'Position', [0.1300 0.2 0.7750 0.79])
% 
% 
% % get 4 out of all 6 conditions
% accuraciesSub = nan(4,iterations);
% accuraciesSub(1,:) = accuracies(1,1,:); % all features, NN
% accuraciesSub(2,:) = accuracies(2,1,:); % one paw features, NN
% accuraciesSub(3,:) = accuracies(2,2,:); % one paw features, GLM
% accuraciesSub(4,:) = accuracies(3,2,:); % minimal features, GLM
% colors = repelem([0 .24 .49; .6 .75 .23],2,1) .* [1 .5 1 .5]';
% 
% 
% figure('Color', 'white', 'MenuBar', 'none', 'Position', [1964 646 339 300]);
% barPlotRick(accuraciesSub, {'conditionNames', {{'NN full', 'NN 1paw', 'GLM 1paw', 'GLM pos+vel+hgt'}}, ...
%     'lineThickness', 2, 'addBars', true, 'scatColors', 'hsv', 'scatAlpha', .5, 'showStats', false, ...
%     'ylim', [.5 1], 'ytick', [.5 .75 1], 'ylabel', 'accuracy', 'conditionColors', colors})
% set(gca, 'Position', [0.15 0.2 0.75 0.79])
% 
% % save
% file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
%         'baseline_decision_modelAccuracies');
% saveas(gcf, file, 'svg');


%% BUILD MODELS

% settings
iterations = 20;
normalizeData = true;
barColors = [0 .24 .49; .6 .75 .23; .2 .2 .2];

accuracies = nan(3,iterations);

for i = 1:iterations
    
    predDistBin = ismember(predictorNames, 'modPawPredictedDistanceToObs'); % this var will ONLY be included in the GLM, not the neural net
    
    % NEURAL NET
    [X, y, predictorNames, isCategorical] = ...
        prepareDecisionModelData(flat, 'all', 'isBigStep', true, referenceModPaw, normalizeData, {'balanceClasses', true});
    net = patternnet(100);
    net.divideParam.trainRatio = .7;
    net.divideParam.valRatio = .15;
    net.divideParam.testRatio = .15;
    [net, tr] = train(net, X(:,~predDistBin)', y');
    outputs = net(X(:,~predDistBin)');
    accuracies(1,i) = mean(round(outputs(tr.testInd))==y(tr.testInd)');
    
    % NN SHUFFLED
    shuffleInds = randperm(length(tr.testInd));
    accuracies(3,i) = mean(round(outputs(tr.testInd))==y(tr.testInd(shuffleInds))'); % NN shuffled
    
    % GLM
    glm = fitglm(X([tr.trainInd tr.valInd], predDistBin), y([tr.trainInd tr.valInd]), ...
        'Distribution', 'binomial', 'CategoricalVars', isCategorical(predDistBin));
    accuracies(2,i) = mean(round(predict(glm, X(tr.testInd,colBins)))==y(tr.testInd));
    
end

%%
figure('Color', 'white', 'MenuBar', 'none', 'Position', [1964 646 339 300]);
barPlotRick(accuracies, {'conditionNames', {{'NN', 'GLM', 'shuffled'}}, ...
    'lineThickness', 2, 'addBars', true, 'scatColors', 'hsv', 'scatAlpha', .5, 'showStats', false, ...
    'ylim', [0 1], 'ytick', [0 .5 1], 'ylabel', 'accuracy', 'conditionColors', barColors})
set(gca, 'Position', [0.15 0.2 0.75 0.79])

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baseline_decision_modelAccuracies');
saveas(gcf, file, 'svg');


%% PREDICTOR SCATTERS


figure('Color', 'white', 'Position', [2000 400 450 350], 'MenuBar', 'none');
scatterHistoRick([flat.x_paw2]*1000, [flat.obsHgt]*1000, ...
    {'groupId', [flat.isBigStep]+1, 'colors', stepTypeColors, ...
    'xlabel', 'paw position (mm)', 'ylabel', 'obstalce height (mm)', ...
    'showCrossHairs', false, 'scatSize', 4, 'scatAlpha', .4, ...
    'xLims', [-70 -15], 'yLims', [3 11], 'groupHistoLineWidth', 3, ...
    'groupHistoAlpha', .8, 'histoAlpha', 0});

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baseline_decision_PositionObsHgtScatter');
saveas(gcf, file, 'svg');


%% BINNED KINEMATICS AND HEATMAP

% settings
netNum = 1;
binNum = 5;
pctileBins = true;


% initializations
kin = permute(cat(3,flat.modPawKinInterp{:}), [3,1,2]);
kinCtl = permute(cat(3,flat.preModPawKinInterp{:}), [3,1,2]);
[X, y, predictorNames, isCategorical] = ...
        prepareDecisionModelData(flat, predictors{netNum}, outcome, useAllPaws(netNum), referenceModPaw, normalizeData, ...
        {'balanceClasses', false, 'removeNans', false});

% choose binning variable
binVar = nnets{netNum}(X'); % neural network output
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



close all;
figure('color', 'white', 'menubar', 'none', 'position', [2000 100 750 100*binNum], 'InvertHardcopy', 'off');
plotBigStepKin(kin(:,[1,3],:), kinCtl(:,[1,3],:), flat.obsHgt, condition, flat.isBigStep, ...
    {'colors', stepTypeColors, 'xLims', [-.09 .04], 'addHistos', false, 'lineWid', 3, ...
    'contactInds', flat.contactInd, 'histoHgt', .5, 'showSmpNum', false})

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baseline_decision_kinematics');
saveas(gcf, file, 'svg');


%% HEATMAP

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

%% SANDBOX



% LASSO
[glmLasso, lassoFit] = lassoglm(xNorm([tr.trainInd tr.valInd],:), y([tr.trainInd tr.valInd]), ...
    'binomial', 'PredictorNames', predictors); % !!! categofical vars not explicitly defined???
% fprintf('GLM_lasso train: %.2f, test: %.3f\n', mean(round(predict(glmLasso, xNorm(tr.trainInd,:)))==y(tr.trainInd)), ...
%     mean(round(predict(glmLasso, xNorm(tr.testInd,:)))==y(tr.testInd)))

zeroInds = nan(1, length(predictors));
for i = 1:size(glmLasso,1); zeroInds(i) = find(abs(glmLasso(i,:))>0,1,'last'); end
[~, sortInds] = sort(zeroInds, 'descend');
predictorsSorted = predictors(sortInds)

% lassoPlot(glmLasso, lassoFit, 'PlotType', 'Lambda', 'PredictorNames', predictors); legend(predictors)

% !!! check accuracy with different numbers of predictors here...




% STEPWISE REGRESSION

glmStepwise = stepwiseglm(xNorm([tr.trainInd tr.valInd],:), y([tr.trainInd tr.valInd]), 'constant', ...
    'Upper', 'linear', 'Distribution', 'binomial', 'VarNames', [predictors outcome], 'PEnter', .001, 'CategoricalVars', isCategorical);
fprintf('GLM_stepwise train: %.2f, test: %.3f\n', mean(round(predict(glmStepwise, xNorm(tr.trainInd,:)))==y(tr.trainInd)), ...
    mean(round(predict(glmStepwise, xNorm(tr.testInd,:)))==y(tr.testInd)))













