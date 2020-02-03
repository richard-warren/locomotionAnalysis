%% load experiment data
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

% settings
outcome = 'isModPawLengthened';  % isBigStep or isModPawLengthened
modPawOnlySwing = false;  % if true, only include trials where the modified paw is the only one in swing
lightOffOnly = true;  % whether to restrict to light on trials
successOnly  = false;  % whether to restrict to successful trials
deltaMin = .5;  % exclude little step trials where modPawDeltaLength is less than deltaLim standard deviations
predictors = {'velAtWiskContact', 'angleAtWiskContact', 'obsHgt', 'wiskContactPosition', 'modPawX', 'modPawXVel', 'modPawZ', 'modPawZVel'};

global_config;


%% schematic imgs for step type decision

session = '180630_001';
yMax = [];  % set to 140 to trim the bottom view out

close all
showDecisionFrames(session, 'stepColors', decisionColors, 'yMax', yMax, ...
    'contrastLims', contrast, 'leftPadding', 0, 'rightPadding', 70, 'drawKinematics', false);


%% overlay images

% settings
sessions = {'180628_004', '180703_001', '180711_000', '180803_003'};
trialsToOverlay = 10;

for s = 1:length(sessions)
    vidTop = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{s}, 'runTop.mp4'));
    vidBot = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{s}, 'runBot.mp4'));
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{s}, 'runAnalyzed.mat'), ...
        'frameTimeStamps', 'obsOnTimes', 'wiskContactTimes');
    wiskContactTimes = wiskContactTimes(~isnan(wiskContactTimes));
    frameTimes = sort(wiskContactTimes(randperm(length(wiskContactTimes), trialsToOverlay)));

    imgsTop = uint8(zeros(vidTop.Height, vidTop.Width, trialsToOverlay));
    imgsBot = uint8(zeros(vidBot.Height, vidBot.Width, trialsToOverlay));

    for i = 1:trialsToOverlay
        imgsTop(:,:,i) = rgb2gray(read(vidTop, knnsearch(frameTimeStamps, frameTimes(i))));
        imgsBot(:,:,i) = rgb2gray(read(vidBot, knnsearch(frameTimeStamps, frameTimes(i))));
    end

    overlayTop = overlayImgs(imgsTop, 'colors', 'jet', 'contrastLims', [.3 .75], 'cutoff', 100, 'projection', 'mean');
    overlayBot = overlayImgs(imgsBot, 'colors', 'jet', 'contrastLims', [.3 .5], 'cutoff', 100, 'projection', 'mean');
    overlay = cat(1,overlayTop, overlayBot);
    figure(); imshow(overlay);


    % write image to desk
    file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'imgs', ['variability' num2str(s) '.png']);
    fprintf('writing %s to disk...\n', file)
    imwrite(overlay, file);
end


%% whisker contact frame image for making schematic of predictor variables

% settings
session = '180630_001';
trial = 19;


vidRun = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4'));
vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'frameTimeStampsWisk', 'wiskContactFrames')
runContactFrame = knnsearch(frameTimeStamps, frameTimeStampsWisk(wiskContactFrames(trial)));
frame = getFrameWithWisk(vidRun, vidWisk, frameTimeStamps, frameTimeStampsWisk, runContactFrame, ...
    'runContrast', contrast);

figure('position', [1937.00 518.00 size(frame,2) size(frame,1)], 'MenuBar', 'none');
colormap gray
image(frame, 'cdatamapping', 'scaled')
set(gca, 'Position', [0 0 1 1], 'Visible', 'off')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'imgs', 'contactFrame.png'));



%% decision making

flat = flattenData(data, ...
    [m.predictors, {'mouse', 'isModPawLengthened', 'modPawDeltaLength', 'isBigStep', 'isLightOn', 'modPawOnlySwing', 'isTrialSuccess', 'modPawPredictedDistanceToObs', 'modPawDistanceToObs', 'modPawKinInterp', 'preModPawKinInterp', 'firstModPaw'}]);

%% heatmaps
plotDecisionHeatmaps(flat, ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'avgMice', true, 'plotMice', false, 'colors', sensColors, 'outcome', 'isModPawLengthened', ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceHeatmaps'));

%% trials scatters
plotDecisionTrials(flat, 'outcome', 'isBigStep', ...
    'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', decisionColors, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceDecisionKin'));

%% model accuracies
plotModelAccuracies(flat, m.predictors, 'isModPawLengthened', ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', [sensColors; .2 .2 .2], 'barProps', barProperties, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceModels'));

%% decision threshold
plotDecisionThresholds(flat, 'outcome', 'isModPawLengthened', ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', sensColors, 'barProps', barProperties, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceThresholds'));



%% !!! forward feature selection



%% !!! predictor scatters and histograms

% f1 = figure('Color', 'white', 'MenuBar', 'none', 'Position', [1939.00 431.00 778.00 587.00]);
% f2 = figure('Color', 'white', 'MenuBar', 'none', 'Position', [2750 554.00 523.00 464.00]);
% set(0, 'CurrentFigure', f1)
% scatSz = 8;
% scatAlpha = .4;
% maxScatters = 1000;  % only plot this many per condition to avoid large vector images
% percentileLims = [1 99];
% binNum = 15;
% maxPlots = 8;
% 
% 
% sz = min(maxPlots, length(predictors));
% randInds = randsample(size(X,1), min(maxScatters, size(X,1)));
% X_sub = X(randInds,:);
% y_sub = y(randInds,:);
% lims = prctile(X, percentileLims, 1)';
% 
% for r = 1:sz
%     for c = 1:sz
%         
%         subplot(sz, sz, (r-1)*sz + c); hold on
%         r_X = sortInds(r);
%         c_X = sortInds(c);
%         
%         % if along diagonal, plot histogram
%         if r==c
%             edges = linspace(lims(r_X,1), lims(r_X,2), binNum+1);
%             histogram(X(~logical(y),r_X), edges, ...
%                 'FaceColor', decisionColors(1,:));
%             histogram(X(logical(y),r_X), edges, ...
%                 'FaceColor', decisionColors(2,:))
%             
%             % add to second figure that has only histograms, no scatters
%             set(0, 'CurrentFigure', f2)
%             subplot(2,2,r); hold on
%             histogram(X(~logical(y),r_X), edges, ...
%                 'FaceColor', decisionColors(1,:));
%             histogram(X(logical(y),r_X), edges, ...
%                 'FaceColor', decisionColors(2,:))
%             title(predictorNames{r_X}, 'Color', predictorColors(r,:), 'Interpreter', 'none');
%             set(gca, 'XLim', lims(r_X,:), 'Box', 'off', 'XTick', [], 'YTick', [])
%             set(0, 'CurrentFigure', f1)
% 
% 
%         
%         % otherwise scatter
%         elseif r>c
%             scatter(X_sub(~logical(y_sub),c_X), X_sub(~logical(y_sub),r_X), ...
%                 scatSz, decisionColors(1,:), 'filled', 'MarkerFaceAlpha', scatAlpha);
%             scatter(X_sub(logical(y_sub),c_X), X_sub(logical(y_sub),r_X), ...
%                 scatSz, decisionColors(2,:), 'filled', 'MarkerFaceAlpha', scatAlpha);
%             set(gca, 'XLim', lims(c_X,:), 'YLim', lims(r_X,:))
%         else
%             set(gca, 'Visible', 'off')
%         end
%         
%         % pimp fig
%         set(gca, 'XTick', [], 'YTick', [])
%         if c==1
%             ylabel(predictorNames{r_X}, 'Interpreter', 'none', 'FontWeight', 'bold', 'Color', predictorColors(r,:));
%         end
%         if r==sz
%             xlabel(predictorNames{c_X}, 'Interpreter', 'none', 'FontWeight', 'bold', 'Color', predictorColors(c,:));
%         end
%     end
%     pause(.001)
% end
% 
% % save
% file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselinePredictors');
% saveas(f1, file, 'svg');
% file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselinePredictorsHistosOnly');
% saveas(f2, file, 'svg');




% %% !!! binned kinematics
% 
% % settings
% binNum = 5;
% pctileBins = false;
% xLims = [-.08 .03];
% 
% 
% % initializations
% kin = permute(cat(3,flat.modPawKinInterp{:}), [3,1,2]);
% kinCtl = permute(cat(3,flat.preModPawKinInterp{:}), [3,1,2]);
% [X, y, predictorNames_full, isCategorical] = ...  % most inclusive model
%     prepareDecisionModelData(flat, ...
%     {'velAtWiskContact', 'angleAtWiskContact', 'obsHgt', 'wiskContactPosition', 'modStepStart', 'x', 'z', 'xVel', 'zVel'}, ...
%     outcome, 'balanceClasses', false, 'useAllPaws', true, 'normalizeData', true, 'referencePaw', referenceModPaw);  % note that balanceClasses is FALSE, which means rows of 'X' correspond to rows of 'flat'
% glmPredictorBins = ismember(predictorNames_full, predictorNames(sortInds(1:maxFeatures)));  % bins of predictors to be used in restricted dataset
% 
% % train model with ALL data
% glm = fitglm(X(:, glmPredictorBins), y(:), ...
%     'Distribution', 'binomial', 'CategoricalVars', isCategorical(glmPredictorBins));
% 
% 
% 
% % choose binning variable
% binVar = predict(glm, X(:,glmPredictorBins));
% % binVar = net(X(:,~predDistBin)'); % neural network output
% % binVar = -flat.x_paw2; % position of mod paw at moment of contact
% % binVar = flat.velAtWiskContact; % vel
% % binVar = flat.modStepStart_paw2; % starting position of mod paw
% % binVar = flat.modPawPredictedDistanceToObs;
% 
% 
% % perctentile bins    
% if pctileBins
%     binEdges = prctile(binVar, linspace(0, 100, binNum+1));
% % evenly spaced bins
% else
%     binEdges = linspace(min(binVar), max(binVar), binNum+1);
% end
% condition = discretize(binVar, binEdges);
% 
% 
% figure('color', 'white', 'menubar', 'none', 'position', [2000.00 100.00 750.00 707.00], 'InvertHardcopy', 'off');
% plotBigStepKin(kin(:,[1,3],:), kinCtl(:,[1,3],:), flat.obsHgt, condition, flat.(outcome), ...
%     'colors', decisionColors, 'xLims', xLims, 'addHistos', false, 'lineWid', 3, ...
%     'contactInds', flat.contactInd, 'histoHgt', .5, 'showSmpNum', false, 'obsColor', obsColor)
% 
% % save
% file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', ...
%         'baseline_decision_kinematics');
% saveas(gcf, file, 'svg');



% %% !!! success rate by model prediction adherence
% 
% % settings
% confidenceLims = [.4 .6];  % only include trials where model is very confident
% 
% % compute overall success rate
% glm = fitglm(X(:, glmPredictorBins), y, ...
%     'Distribution', 'binomial', 'CategoricalVars', isCategorical(glmPredictorBins));
% predictions = predict(glm, X(:,glmPredictorBins));
% confidentBins = predictions<confidenceLims(1) | predictions>confidenceLims(2);
% isModelCorrect = round(predictions) == [flat.(outcome)];  % bins where decision agrees with model prediction
% fprintf('POOLED MODEL\n')
% fprintf('correct prediction success:   %.3f (n=%i)\n', ...
%     mean(flat.isTrialSuccess(isModelCorrect & confidentBins)), sum(isModelCorrect & confidentBins))
% fprintf('incorrect prediction success: %.3f (n=%i)\n\n', ...
%     mean(flat.isTrialSuccess(~isModelCorrect & confidentBins)), sum(~isModelCorrect & confidentBins))
% 
% % build per-mouse models
% mice = unique(flat.mouse);
% successRates = nan(2, length(mice));  % (correct, incorrect model predictions) X mouse
% accuracies = nan(1,length(mice));
% 
% for i = 1:length(mice)
%     mouseBins = strcmp(flat.mouse, mice{i});
%     glm_temp = fitglm(X(mouseBins, glmPredictorBins), y(mouseBins), ...
%         'Distribution', 'binomial', 'CategoricalVars', isCategorical(glmPredictorBins));
%     
%     predictions = predict(glm, X(mouseBins, glmPredictorBins));
%     confidentBins = predictions<confidenceLims(1) | predictions>confidenceLims(2);
%     isModelCorrect = round(predictions) == flat.(outcome)(mouseBins);
%     mouseSuccess = flat.isTrialSuccess(mouseBins);
%     
%     successRates(1,i) = mean(mouseSuccess(confidentBins & isModelCorrect));
%     successRates(2,i) = mean(mouseSuccess(confidentBins & ~isModelCorrect));
%     
%     accuracies(i) = mean(isModelCorrect);
% end
% 
% fprintf('INDIVIDUAL MOUSE MODELS\n')
% fprintf('correct prediction success:   %.3f\n', mean(successRates(1,:)))
% fprintf('incorrect prediction success: %.3f\n', mean(successRates(2,:)))
% fprintf('overall accuracy: %.3f\n\n', mean(accuracies))
% [h, p] = ttest(successRates(1,:), successRates(2,:));
% fprintf('significance: %.3f\n\n', p)
% 
% 
% figure('Color', 'white', 'Position', [1987.00 448.00 300.00 291.00], 'MenuBar', 'none');
% barFancy(successRates, barProperties{:}, ...
%     'levelNames', {{'correct', 'incorrect'}}, 'textRotation', 0, 'colors', [modelColor; modelColor*.2])
% set(gca, 'YTick', 0:.5:1)
% 
% file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineModelIncorrectAccuracy');
% fprintf('writing %s to disk...\n', file)
% saveas(gcf, file, 'svg');


%% distance and time to contact


% settings
xLims = [15 50];
yLims = [0 150];
scatSize = 4;
scatAlpha = .3;
mouseColors = true;
scatPlotSize = .7;
border = .15;


% initializations
d = flat.distances;
t = flat.times;
[~,~,mouseIds] = unique(flat.mouse);

% plot that shit
figure('Color', 'white', 'Position', [2000 400 450 350], 'MenuBar', 'none');
scatterHistoRick(d,t, ...
    'groupId', mouseIds, 'colors', 'jet', 'groupFcn', @nanmedian, ...
    'xlabel', 'distance to contact (mm)', 'ylabel', 'time to contact (ms)', ...
    'xLims', xLims, 'yLims', yLims, 'showCrossHairs', true, 'scatSize', scatSize, 'scatAlpha', scatAlpha);

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'baselineDistanceTimeToContact');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

% get median distances and times for all mice
[distances, times] = deal(nan(1,length(mice)));
for i = 1:length(mice)
    mouseBins = strcmp(flat.mouse, mice{i});
    distances(i) = median(flat.distances(mouseBins));
    times(i) = median(flat.times(mouseBins));
end

fprintf('\ndistance to contact:   %.1f +- %.1f SEM\n', mean(distances), std(distances)/sqrt(length(distances)))
fprintf('time to contact:       %.1f +- %.1f SEM\n', mean(times), std(times)/sqrt(length(times)))

