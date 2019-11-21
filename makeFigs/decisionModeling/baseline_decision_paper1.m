%% COMPUTE PREDICTORS

% !!! should i restrict to light off trials???

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



%% BUILD MODELS

% settings
iterations = 20;
normalizeData = true;
barColors = [0 .24 .49; .6 .75 .23; .2 .2 .2];

% prepare training data
[X, y, predictorNames, isCategorical] = ...
    prepareDecisionModelData(flat, 'all', 'isBigStep', true, referenceModPaw, normalizeData, {'balanceClasses', true});
predDistBin = ismember(predictorNames, 'modPawPredictedDistanceToObs'); % this var will ONLY be included in the GLM, not the neural net
    
    
accuracies = nan(3,iterations);
for i = 1:iterations
    
    % NEURAL NET
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
    glmPredictorBins = true(1,size(X,2));  % bins of predictors to include in GLM (columns of X)
    glm = fitglm(X([tr.trainInd tr.valInd],glmPredictorBins), y([tr.trainInd tr.valInd]), ...
        'Distribution', 'binomial', 'CategoricalVars', isCategorical(glmPredictorBins));
    accuracies(2,i) = mean(round(predict(glm, X(tr.testInd,glmPredictorBins)))==y(tr.testInd));
    
end

%
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



