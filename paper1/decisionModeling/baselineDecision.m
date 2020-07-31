%% load experiment data
% tic; fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!'); toc
fprintf('loading... '); load(fullfile('C:\Users\richa\Desktop\', 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!');  % temp, use to load from local drive because engram is slow over ethernet

% settings
outcome = 'isModPawLengthened';  % isBigStep or isModPawLengthened

global_config;

flat = flattenData(data, ...
    [m.predictorsAll, {'mouse', 'session', 'trial', 'isModPawLengthened', 'modPawDeltaLength', 'isBigStep', 'isLightOn', ...
    'modPawOnlySwing', 'isTrialSuccess', 'modPawPredictedDistanceToObs', 'modPawDistanceToObs', 'modPawKinInterp', ...
    'preModPawKinInterp', 'firstModPaw', 'contactIndInterp', 'preModPawDeltaLength', 'contactInd', 'modSwingContacts'}]);

predictorColors = lines(length(m.predictorsAll));
mice = unique({flat.mouse});

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
plotDecisionHeatmaps(flat, 'normalize', 'col', 'outcome', 'isModPawLengthened', ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'avgMice', false, 'plotMice', false, 'colors', decisionColors(1,:), 'xLims', [-20 15], 'plotProbs', false);
set(gcf, 'position', [200 379 295 362])
xlabel('predicted landing distance (mm)')
ylabel('landing distance (mm)')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineHeatmaps'), 'svg');

%% (temp, for exploring whether heatmaps are blurrier with no light) heatmaps
plotDecisionHeatmaps(flat, 'condition', 'isLightOn', 'levels', {0,1}, 'normalize', 'col', 'outcome', 'isModPawLengthened', ...
    'deltaMin', 0, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', false, ...
    'avgMice', true, 'plotMice', true, 'colors', decisionColors(1,:), 'xLims', [-20 15], 'plotProbs', false);

%% (temp, for exploring whether 'better' mice have sharper decision boundaries) explore whether 'better' mice have sharper decisions

% get heatmaps
plotDecisionHeatmaps(flat, 'normalize', '', 'outcome', 'isModPawLengthened', ...
    'deltaMin', 0, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'avgMice', true, 'plotMice', true, 'colors', decisionColors(1,:), 'xLims', [-20 15], 'plotProbs', false, ...
    'subplotDims', [2, ceil((length(mice)+1)/2)]);
close(gcf);

% get entropy
entropies = plotEntropies(flat, ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'barProps', barProperties);
close(gcf);

% get success and velocity
[success, velocity] = deal(nan(1,length(mice)));
for i = 1:length(mice)
    bins = strcmp({flat.mouse}, mice{i});
    success(i) = nanmean([flat(bins).isTrialSuccess]);
    velocity(i) = nanmean([flat(bins).velAtWiskContact]);
end

%%
close all
figure('Position', [36 1038 806 305], 'color', 'white', 'MenuBar', 'none');

subplot(1,2,1)
scatter(success, entropies, 'filled')
xlabel('success rate')
ylabel('entropy')

subplot(1,2,2)
scatter(velocity, entropies, 'filled')
xlabel('velocity (m/s)')
ylabel('entropy')

% subplot(1,3,3)
% scatter3(success, velocity, entropies, 'filled')
% xlabel('success')
% ylabel('velocity (m/s)')
% zlabel('entropy')


% plot heatmaps sorted by entropy
figure('Position', [36 522 1282 821], 'color', 'white', 'MenuBar', 'none');

[~, inds] = sort(entropies);

entropiesSorted = entropies(inds);
heatmapsSorted = heatmaps(inds);
velocitySorted = velocity(inds);
successSorted = success(inds);

colormap('hot')
for i = 1:length(heatmapsSorted)
    subplot(3,ceil(length(heatmapsSorted)/3),i)
    imagesc(heatmapsSorted{i})
    set(gca, 'YDir', 'normal', 'TickDir', 'out', 'Box', 'off')
    
    title({sprintf('entropy: %.2f', entropiesSorted(i)), ...
           sprintf('velocity: %.2f', velocitySorted(i)), ...
           sprintf('success: %.2f', successSorted(i))})
end



%% trials scatters
rng(1)
plotDecisionTrials(flat, 'outcome', 'isModPawLengthened', ...
    'modSwingContactsMax', 0, 'deltaMin', m.deltaMin, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...  % don't limit to successful trials for this plot
    'colors', decisionColors, 'view', 'top', 'xLims', [-.11 .06], 'obsColor', obsColor, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineDecisionKin'));

%% model accuracies
[~,~,flat_restricted] = plotModelAccuracies(flat, m.predictors, 'isModPawLengthened', 'model', 'glm', ...
    'modSwingContactsMax', m.modSwingContactsMax, 'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'weightClasses', true, 'barProps', barProperties, 'kFolds', 15, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineModels'));

%% model predictors
plotPredictors(flat, m.predictors, 'isModPawLengthened', 'avgMice', true, 'colors', predictorColors,...
    'modSwingContactsMax', m.modSwingContactsMax, 'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'mouseAlpha', .2, 'subplotDims', [2 4], 'names', m.predictorsNamed, 'figPos', [100 400 671 300], ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselinePredictors'));

%% model accuracies (predicted distance only)
acc = plotModelAccuracies(flat, {'modPawPredictedDistanceToObs'}, 'isModPawLengthened', 'model', 'glm', ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'weightClasses', true, 'barProps', barProperties, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineModelsPredDistOnly'));
fprintf('predicted distance only model accuracy: %.3f\n', mean(acc(1,2,:)));

%% decision threshold
plotDecisionThresholds(flat, 'outcome', 'isModPawLengthened', ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', sensColors, 'barProps', barProperties, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineThresholds'));

%% distribution of step modifications (this demonstrates the effect of different deltaMin strategies)
close all
xLims = prctile([flat_sub.modPawDeltaLength flat_sub.preModPawDeltaLength], [1 99]);
binEdges = linspace(xLims(1), xLims(2), 75);

figure('position', [63.00 403.00 836.00 782.00], 'color', 'white', 'menubar', 'none')

extraBins = [flat.isBigStep]==0;
% extraBins = true(1, length(flat));

subplot(4,1,1)
bins = extraBins;
histogram([flat(bins).modPawDeltaLength], binEdges); hold on
histogram([flat(bins).preModPawDeltaLength], binEdges); hold on
set(gca, 'box', 'off', 'YColor', 'none')

subplot(4,1,2)
minDif = std([flat.preModPawDeltaLength]) * 1;
bins = abs([flat.modPawDeltaLength])>minDif & extraBins;
histogram([flat(bins).modPawDeltaLength], binEdges); hold on
histogram([flat(bins).preModPawDeltaLength], binEdges); hold on
set(gca, 'box', 'off', 'YColor', 'none')

subplot(4,1,3)
bins = ~(abs(zscore([flat.modPawDeltaLength]))<m.deltaMin) & extraBins;
histogram([flat(bins).modPawDeltaLength], binEdges); hold on
histogram([flat(bins).preModPawDeltaLength], binEdges); hold on
set(gca, 'box', 'off', 'YColor', 'none')

subplot(4,1,4)
bins = ~(abs(zscore([flat.modPawDeltaLength]))<m.deltaMin & [flat.isBigStep]==0) & extraBins;
histogram([flat(bins).modPawDeltaLength], binEdges); hold on
histogram([flat(bins).preModPawDeltaLength], binEdges); hold on
set(gca, 'box', 'off', 'YColor', 'none')

xlabel('\delta length')
print -clipboard -dbitmap

%% forward feature selection

% settings
kFolds = 15;
rng(0);  % random seed initialization

% initializations
bestPredictors = {};
remainingPredictors = m.predictorsAll;
mice = unique({flat.mouse});
accuracies = nan(length(mice), length(m.predictorsAll));

for i = 1:length(m.predictorsAll)
    
    fprintf('predictor %i: ', i)
    a = nan(length(mice), length(remainingPredictors));
    
    for j = 1:length(remainingPredictors)
        fprintf('%i/%i ', j, length(remainingPredictors))
        temp = plotModelAccuracies(flat, [bestPredictors, remainingPredictors{j}], 'isModPawLengthened', ...
            'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
            'weightClasses', true, 'plot', false, 'kFolds', kFolds, 'verbose', false);
        a(:,j) = temp(1,2,:);
    end
    fprintf('\n')
    
    % store best predictor
    [~, ind] = max(nanmean(a,1));
    accuracies(:,i) = a(:, ind);
    bestPredictors{end+1} = remainingPredictors{ind};
    remainingPredictors = remainingPredictors(~strcmp(remainingPredictors, remainingPredictors{ind}));
end

% plot forward feature selection
figure('Color', 'white', 'MenuBar', 'none', 'Position', [200 747.00 297.00 232.00]); hold on
plot(0:length(m.predictorsAll), [.5 nanmean(accuracies,1)], ...
    'LineWidth', 1, 'Color', [.2 .2 .2]);
scatter(repelem(1:length(m.predictorsAll), length(mice)), ...
    accuracies(:), 20, [1 1 1]*.8, 'filled', 'MarkerFaceAlpha', .4);
% plot(repmat(1:length(m.predictorsAll),length(mice),1)', ...
%     accuracies', 20, [1 1 1]*.8, 'color', [0 0 0 .1]);
scatter(1:length(m.predictorsAll), nanmean(accuracies,1), 60, predictorColors, 'filled');



xlabel('number of features')
ylabel('cross-validation accuracy')
set(gca, 'XLim', [0, length(m.predictorsAll)])
fprintf('max accuracy: %.2f\n', max(mean(accuracies,1)))
fprintf('PREDICTORS: {');fprintf('''%s'', ', bestPredictors{:}); fprintf('}\n')

% save
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineModelAccuracy'), 'svg');


%% predictor scatters and histograms

[~, ~, flat_sub] = plotModelAccuracies(flat, m.predictors, outcome, ...
    'modSwingContactsMax', m.modSwingContactsMax, 'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'weightClasses', true, 'plot', false, 'kFolds', 2);

% settings
binNum = 15;  % number of histogram bins
maxPlots = length(m.predictors);
scatAlpha = .4;
maxScatters = 500;  % only plot this many per condition to avoid large vector images
percentileLims = [1 99];
scatSz = 8;


f1 = figure('Color', 'white', 'MenuBar', 'none', 'Position', [169 169 1037 843]);
set(0, 'CurrentFigure', f1)

sz = min(maxPlots, length(m.predictors));
randInds = randsample(length(flat_sub), min(maxScatters, length(flat_sub)));
randBins = false(1, length(flat_sub));
randBins(randInds) = true;

for r = 1:sz
    x_row = [flat_sub.(m.predictors{r})];
    lims_row = prctile(x_row, percentileLims);
    
    for c = 1:sz
        x_col = [flat_sub.(m.predictors{c})];
        lims_col = prctile(x_col, percentileLims);
        
        subplot(sz, sz, (r-1)*sz + c); hold on
        
        % if along diagonal, plot histogram
        if r==c
            edges = linspace(lims_row(1), lims_row(2), binNum+1);
            histogram(x_row([flat_sub.(outcome)]~=1), edges, ...
                'FaceColor', decisionColors(1,:));
            histogram(x_row([flat_sub.(outcome)]==1), edges, ...
                'FaceColor', decisionColors(2,:));
        
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
            ylabel(m.predictorsNamed{r}, 'Interpreter', 'none', 'FontWeight', 'bold', 'Color', predictorColors(r,:));
        end
        if r==sz
            xlabel(m.predictorsNamed{c}, 'Interpreter', 'none', 'FontWeight', 'bold', 'Color', predictorColors(c,:));
        end
    end
    pause(.001)
end

% save
saveas(f1, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselinePredictorScatters'), 'svg');


%% binned kinematics

% settings
binNum = 5;
pctileBins = true;
xLims = [-.08 .03];

% initializations
kin = permute(cat(3,flat_sub.modPawKinInterp), [3,1,2]);
kinCtl = permute(cat(3,flat_sub.preModPawKinInterp), [3,1,2]);

% choose binning variable
binVar = [flat_sub.([outcome '_predicted'])];
% binVar = [flat_sub.modPawX];

% perctentile bins    
if pctileBins
    binEdges = prctile(binVar, linspace(0, 100, binNum+1));
% evenly spaced bins
else
    binEdges = linspace(min(binVar), max(binVar), binNum+1);
end
condition = discretize(binVar, binEdges);

figure('color', 'white', 'menubar', 'none', 'position', [200 100.00 750.00 707.00], 'InvertHardcopy', 'off');
plotBigStepKin(kin(:,[1,3],:), kinCtl(:,[1,3],:), [flat_sub.obsHgt], condition, [flat_sub.isBigStep]==1, ...
    'colors', decisionColors, 'xLims', xLims, 'addHistos', false, 'lineWid', 3, ...
    'contactInds', [flat_sub.contactIndInterp], 'histoHgt', .5, 'showSmpNum', false, 'obsColor', obsColor)
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baseline_decision_kinematics'), 'svg');



%% success rate by model prediction adherence (!!! need to make a new flat_sub trained without unsuccessful trials removed)

% settings
confidenceLims = [.5 .5];  % only include trials where model is very confident

% get model predictors (stored in flat_temp), when including all failure trials
[~,~,flat_temp] = plotModelAccuracies(flat, m.predictors, 'isModPawLengthened', 'model', 'glm', ...
    'modSwingContactsMax', m.modSwingContactsMax, 'deltaMin', m.deltaMin, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'weightClasses', true, 'barProps', barProperties, 'kFolds', 15, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineModels'));



%% compute overall success rate
predictions = [flat_temp.isModPawLengthened_predicted];
confidentBins = predictions<=confidenceLims(1) | predictions>confidenceLims(2);
isModelCorrect = round(predictions) == [flat_temp.(outcome)];  % bins where decision agrees with model prediction
fprintf('POOLED MODEL\n')
fprintf('correct prediction success:   %.3f (n=%i)\n', ...
    mean([flat_temp(isModelCorrect & confidentBins).isTrialSuccess]), sum(isModelCorrect & confidentBins))
fprintf('incorrect prediction success: %.3f (n=%i)\n\n', ...
    mean([flat_temp(~isModelCorrect & confidentBins).isTrialSuccess]), sum(~isModelCorrect & confidentBins))


% individual mouse rates
mice = unique({flat_temp.mouse});
successRates = nan(2, length(mice));  % (correct, incorrect model predictions) X mouse
accuracies = nan(1,length(mice));

for i = 1:length(mice)
    mouseBins = strcmp({flat_temp.mouse}, mice{i});
    
    predictions = [flat_temp(mouseBins).([outcome '_predicted'])];
    confidentBins = predictions<confidenceLims(1) | predictions>confidenceLims(2);
    isModelCorrect = round(predictions) == [flat_temp(mouseBins).(outcome)];
    mouseSuccess = [flat_temp(mouseBins).isTrialSuccess];
    
    successRates(1,i) = mean(mouseSuccess(confidentBins & isModelCorrect));
    successRates(2,i) = mean(mouseSuccess(confidentBins & ~isModelCorrect));
    
    accuracies(i) = mean(isModelCorrect);
end

fprintf('INDIVIDUAL MOUSE MODELS\n')
fprintf('correct prediction success:   %.3f\n', nanmean(successRates(1,:)))
fprintf('incorrect prediction success: %.3f\n', nanmean(successRates(2,:)))
fprintf('overall accuracy: %.3f\n\n', mean(accuracies))


figure('Color', 'white', 'Position', [1864 322 218 291], 'MenuBar', 'none');
barFancy(successRates, barProperties{:}, ...
    'levelNames', {{'correct', 'incorrect'}}, 'textRotation', 0, 'colors', [modelColor; modelColor*.2], ...
    'comparisons', [1 2], 'test', 'ttest', 'ylabel', 'success rate')
set(gca, 'YTick', 0:.5:1)

file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineModelIncorrectAccuracy');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% see how model accuracy suffers without whisker info

rng(0)
excludeVars = {'obsHgt', 'wiskContactPosition'};
predictorsNoWisk = m.predictors(~ismember(m.predictors, excludeVars));
[acc_full] = plotModelAccuracies(flat, m.predictors, outcome, ...
    'modSwingContactsMax', m.modSwingContactsMax, 'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'weightClasses', true, 'plot', false, 'kFolds', 10);
[acc_noWisk] = plotModelAccuracies(flat, predictorsNoWisk, outcome, ...
    'modSwingContactsMax', m.modSwingContactsMax, 'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'weightClasses', true, 'plot', false, 'kFolds', 10);

mat = [squeeze(acc_full(1,2,:))'; squeeze(acc_noWisk(1,2,:))'];  % (full vs. noWisk) X (mouse)

figure('Color', 'white', 'Position', [1864 322 218 291], 'MenuBar', 'none');
barFancy(mat, barProperties{:}, ...
    'levelNames', {{'full model', 'no whisker predictors'}}, barProperties{:}, ...
    'comparisons', [1 2], 'test', 'ttest', 'scatterColors', 'lines', 'scatterCondColor', false, ...
    'ylabel', 'cross-validation accuracy', 'scatterAlpha', .6);
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineModelNoWiskPredModels'), 'svg');


%% distance and time to contact


% compute distance and time to contact
trialSmps = 100;  % fit linear model of leading forepaw position over this many samples
fps = 250;
sessions = unique({flat.session});
% sessions = sessions(1:6); % !!! temp

% initialize new fields
temp = num2cell(nan(1,length(flat)));
[flat.distances] = temp{:};
[flat.times] = temp{:};

for i = 1:length(sessions)
    
    fprintf('%s: computing distance and time to contact\n', sessions{i})
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'kinData.mat'), 'kinData');
    sessionBins = strcmp({flat.session}, sessions{i});
    if sum(sessionBins)~=length(kinData); disp('WTF! mismatch in number of trials! fucking shit!!!'); end
    
    % loop over paws
    for j = find([kinData.isTrialAnalyzed])
        
        flatBin = sessionBins & [flat.trial]==j;
        
        % find distance and time to contact
        flat(flatBin).distances = abs(max(kinData(j).locations(kinData(j).contactInd,1,:)))*1000;  % distance of leading paw at contact
        trialX = max(kinData(j).locations(kinData(j).contactInd-trialSmps+1:kinData(j).contactInd,1,:), [], 3);  % x position of the leading most forepaw across entire trial
        predictedAtObsInd = polyval(polyfit(trialX', 1:trialSmps, 1), 0);  % predicted the ind at which paw would intercept obstacle
        flat(flatBin).times = abs((predictedAtObsInd-trialSmps) / fps)*1000; % frame until contact / (frames/second)
    end
end
disp('all done!')





%% settings
binNum = 100;
mouseAlpha = .25;


% initializations
d = [flat.distances];
t = [flat.times];
[mice,~,mouseIds] = unique({flat.mouse});

% plot that shit
close all
figure('Color', 'white', 'Position', [200 400 450 350], 'MenuBar', 'none');
dvs = {d, t};
labels = {'distance to contact (mm)', 'time to contact (ms)'};
xLims = [10 60; 0 150];
colors = lines(length(mice));

for i = 1:2
    subplot(2,1,i); hold on
    
    xGrid = linspace(xLims(i,1), xLims(i,2), binNum);
    
    mouseMedians = nan(1,length(mice));  % mouse means for two conditions
    pdfs = nan(length(mice), binNum);
    
    for j = 1:length(mice)
        bins = mouseIds==j;
        pdfs(j,:) = ksdensity(dvs{i}(bins), xGrid);
        plot(xGrid, pdfs(j,:), 'Color', [colors(j,:) mouseAlpha], 'LineWidth', 1)
        mouseMedians(j) = nanmedian(dvs{i}(bins));
    end
    
    plot(xGrid, nanmean(pdfs,1), 'Color', [1 1 1]*.15, 'LineWidth', 3)
    yLims = get(gca, 'ylim');
    plot([1 1]*mean(mouseMedians), yLims, 'Color', [.15 .15 .15])  % add line at the mean across mice :)
    set(gca, 'YColor', 'none', 'Box', 'off', 'xlim', xLims(i,:), 'xtick', linspace(xLims(i,1), xLims(i,2), 4), 'TickDir', 'out')
    xlabel(labels{i})
    text(mean(mouseMedians)+.05*range(xLims(i,:)), yLims(2), sprintf('%.1f', mean(mouseMedians)), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
end


% fprintf('\ndistance to contact:   %.1f +- %.1f SEM\n', mean(distances), std(distances)/sqrt(length(distances)))
% fprintf('time to contact:       %.1f +- %.1f SEM\n', mean(times), std(times)/sqrt(length(times)))

saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineDistanceTimeToContact'), 'svg');
