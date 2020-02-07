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

flat = flattenData(data, ...
    [m.predictors, {'mouse', 'isModPawLengthened', 'modPawDeltaLength', 'isBigStep', 'isLightOn', 'modPawOnlySwing', 'isTrialSuccess', 'modPawPredictedDistanceToObs', 'modPawDistanceToObs', 'modPawKinInterp', 'preModPawKinInterp', 'firstModPaw', 'contactIndInterp'}]);

% get version of flat with trial restrictions applied, and with predictions of model stored as a colum
[~, ~, flat_sub] = plotModelAccuracies(flat, m.predictors, outcome, ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'weightClasses', true, 'plot', false, 'kFolds', 10);



%% schematic imgs for step type decision

session = '180630_001';
yMax = [];  % set to 140 to trim the bottom view out

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
    'weightClasses', true, 'barProps', barProperties, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceModels'));

%% decision threshold
plotDecisionThresholds(flat, 'outcome', 'isModPawLengthened', ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', sensColors, 'barProps', barProperties, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceThresholds'));

%% forward feature selection

% settings
kFolds = 10;


% initializations
bestPredictors = {};
remainingPredictors = m.predictors;
accuracies = nan(1, length(m.predictors));
predictorColors = lines(length(m.predictors));

for i = 1:length(m.predictors)
    
    fprintf('predictor %i: ', i)
    a = nan(1, length(remainingPredictors));
    
    for j = 1:length(remainingPredictors)
        fprintf('%i/%i ', j, length(remainingPredictors))
        a_sub = plotModelAccuracies(flat, [bestPredictors, remainingPredictors{j}], 'isModPawLengthened', ...
            'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
            'weightClasses', true, 'plot', false, 'kFolds', kFolds);
        a(j) = mean(a_sub(1,2,:));
    end
    fprintf('\n')
    
    % store best predictor
    [accuracies(i), ind] = max(a);
    bestPredictors{end+1} = remainingPredictors{ind};
    remainingPredictors = remainingPredictors(~strcmp(remainingPredictors, remainingPredictors{ind}));
end

% plot
figure('Color', 'white', 'MenuBar', 'none', 'Position', [2375.00 747.00 297.00 232.00]); hold on
plot(1:length(m.predictors), accuracies, ...
    'LineWidth', 1, 'Color', [.2 .2 .2]);
scatter(1:length(m.predictors), accuracies, 60, predictorColors, 'filled');
xlabel('number of features')
ylabel('cross-validation accuracy')
set(gca, 'XLim', [1, length(m.predictors)])
fprintf('max accuracy: %.2f\n', max(accuracies))
fprintf('PREDICTORS: '); fprintf('%s ', bestPredictors{:}); fprintf('\n')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineModelAccuracy');
saveas(gcf, file, 'svg');


%% predictor scatters and histograms

% settings
% bestPredictors = {'modPawX', 'obsHgt', 'velAtWiskContact', 'wiskContactPosition', 'modPawZ', 'modPawZVel', 'modPawXVel', 'angleAtWiskContact'}
binNum = 15;
maxPlots = 5;
scatAlpha = .4;
maxScatters = 2000;  % only plot this many per condition to avoid large vector images
percentileLims = [1 99];
scatSz = 8;


f1 = figure('Color', 'white', 'MenuBar', 'none', 'Position', [1939.00 431.00 778.00 587.00]);
f2 = figure('Color', 'white', 'MenuBar', 'none', 'Position', [2750.00 900.00 160*maxPlots 118.00]);
set(0, 'CurrentFigure', f1)

sz = min(maxPlots, length(predictors));
randInds = randsample(length(flat_sub), min(maxScatters, length(flat_sub)));
randBins = false(1, length(flat_sub));
randBins(randInds) = true;

for r = 1:sz
    x_row = [flat_sub.(bestPredictors{r})];
    lims_row = prctile(x_row, percentileLims);
    
    for c = 1:sz
        x_col = [flat_sub.(bestPredictors{c})];
        lims_col = prctile(x_col, percentileLims);
        
        subplot(sz, sz, (r-1)*sz + c); hold on
        
        % if along diagonal, plot histogram
        if r==c
            edges = linspace(lims_row(1), lims_row(2), binNum+1);
            histogram(x_row([flat_sub.(outcome)]~=1), edges, ...
                'FaceColor', decisionColors(1,:));
            histogram(x_row([flat_sub.(outcome)]==1), edges, ...
                'FaceColor', decisionColors(2,:));
            
            % add to second figure that has only histograms, no scatters
            set(0, 'CurrentFigure', f2)
            subplot(1,maxPlots,r); hold on
            histogram(x_row([flat_sub.(outcome)]~=1), edges, ...
                'FaceColor', decisionColors(1,:));
            histogram(x_row([flat_sub.(outcome)]==1), edges, ...
                'FaceColor', decisionColors(2,:));
            title(bestPredictors{r}, 'Color', predictorColors(r,:), 'Interpreter', 'none');
            set(gca, 'XLim', lims_col, 'Box', 'off', 'XTick', [], 'YTick', [])
            set(0, 'CurrentFigure', f1)

        
        % otherwise scatter
        elseif r>c
            scatter(x_col([flat_sub.(outcome)]~=1 & randBins), x_row([flat_sub.(outcome)]~=1 & randBins), ...
                scatSz, decisionColors(1,:), 'filled', 'MarkerFaceAlpha', scatAlpha);
            scatter(x_col([flat_sub.(outcome)]==1 & randBins), x_row([flat_sub.(outcome)]==1 & randBins), ...
                scatSz, decisionColors(2,:), 'filled', 'MarkerFaceAlpha', scatAlpha);
            
            set(gca, 'XLim', lims_col, 'YLim', lims_row)
        else
            set(gca, 'Visible', 'off')
        end
        
        % pimp fig
        set(gca, 'XTick', [], 'YTick', [])
        if c==1
            ylabel(bestPredictors{r}, 'Interpreter', 'none', 'FontWeight', 'bold', 'Color', predictorColors(r,:));
        end
        if r==sz
            xlabel(bestPredictors{c}, 'Interpreter', 'none', 'FontWeight', 'bold', 'Color', predictorColors(c,:));
        end
    end
    pause(.001)
end

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselinePredictors');
saveas(f1, file, 'svg');
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselinePredictorsHistosOnly');
saveas(f2, file, 'svg');




%% binned kinematics

% settings
binNum = 5;
pctileBins = false;
xLims = [-.08 .03];


% initializations
kin = permute(cat(3,flat_sub.modPawKinInterp), [3,1,2]);
kinCtl = permute(cat(3,flat_sub.preModPawKinInterp), [3,1,2]);


% choose binning variable
binVar = [flat_sub.([outcome '_predicted'])];

% perctentile bins    
if pctileBins
    binEdges = prctile(binVar, linspace(0, 100, binNum+1));
% evenly spaced bins
else
    binEdges = linspace(min(binVar), max(binVar), binNum+1);
end
condition = discretize(binVar, binEdges);


figure('color', 'white', 'menubar', 'none', 'position', [2000.00 100.00 750.00 707.00], 'InvertHardcopy', 'off');
plotBigStepKin(kin(:,[1,3],:), kinCtl(:,[1,3],:), [flat_sub.obsHgt], condition, [flat_sub.(outcome)], ...
    'colors', decisionColors, 'xLims', xLims, 'addHistos', false, 'lineWid', 3, ...
    'contactInds', [flat_sub.contactIndInterp], 'histoHgt', .5, 'showSmpNum', false, 'obsColor', obsColor)

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', ...
        'baseline_decision_kinematics');
saveas(gcf, file, 'svg');



%% success rate by model prediction adherence (!!! need to make a new flat_sub trained without unsuccessful trials removed)

% settings
confidenceLims = [.4 .6];  % only include trials where model is very confident

% compute overall success rate
predictions = [flat_sub.isModPawLengthened_predicted];
confidentBins = predictions<confidenceLims(1) | predictions>confidenceLims(2);
isModelCorrect = round(predictions) == [flat_sub.(outcome)];  % bins where decision agrees with model prediction
fprintf('POOLED MODEL\n')
fprintf('correct prediction success:   %.3f (n=%i)\n', ...
    mean([flat_sub(isModelCorrect & confidentBins).isTrialSuccess]), sum(isModelCorrect & confidentBins))
fprintf('incorrect prediction success: %.3f (n=%i)\n\n', ...
    mean([flat_sub(~isModelCorrect & confidentBins).isTrialSuccess]), sum(~isModelCorrect & confidentBins))


% individual mouse rates
mice = unique({flat_sub.mouse});
successRates = nan(2, length(mice));  % (correct, incorrect model predictions) X mouse
accuracies = nan(1,length(mice));

for i = 1:length(mice)
    mouseBins = strcmp({flat_sub.mouse}, mice{i});
    
    predictions = flat_sub(mouseBins).([outcome '_predicted']);
    confidentBins = predictions<confidenceLims(1) | predictions>confidenceLims(2);
    isModelCorrect = round(predictions) == [flat_sub(mouseBins).(outcome)];
    mouseSuccess = [flat_sub(mouseBins).isTrialSuccess];
    
    successRates(1,i) = mean(mouseSuccess(confidentBins & isModelCorrect));
    successRates(2,i) = mean(mouseSuccess(confidentBins & ~isModelCorrect));
    
    accuracies(i) = mean(isModelCorrect);
end

fprintf('INDIVIDUAL MOUSE MODELS\n')
fprintf('correct prediction success:   %.3f\n', mean(successRates(1,:)))
fprintf('incorrect prediction success: %.3f\n', mean(successRates(2,:)))
fprintf('overall accuracy: %.3f\n\n', mean(accuracies))
[h, p] = ttest(successRates(1,:), successRates(2,:));
fprintf('significance: %.3f\n\n', p)


figure('Color', 'white', 'Position', [1987.00 448.00 300.00 291.00], 'MenuBar', 'none');
barFancy(successRates, barProperties{:}, ...
    'levelNames', {{'correct', 'incorrect'}}, 'textRotation', 0, 'colors', [modelColor; modelColor*.2])
set(gca, 'YTick', 0:.5:1)

file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineModelIncorrectAccuracy');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% see how model accuracy suffers without whisker info

excludeVars = {'obsHgt', 'wiskContactPosition'};
predictorsNoWisk = m.predictors(~ismember(m.predictors, excludeVars));
[acc_full] = plotModelAccuracies(flat, m.predictors, outcome, ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'weightClasses', true, 'plot', false, 'kFolds', 10);
[acc_noWisk] = plotModelAccuracies(flat, predictorsNoWisk, outcome, ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'weightClasses', true, 'plot', false, 'kFolds', 10);

mat = [squeeze(acc_full(1,2,:))'; squeeze(acc_noWisk(1,2,:))'];  % (full vs. noWisk) X (mouse)
figure('Color', 'white', 'Position', [1987.00 448.00 300.00 291.00], 'MenuBar', 'none');
barFancy(mat, barProperties{:}, ...
    'levelNames', {{'full model', 'noWisk'}});

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

