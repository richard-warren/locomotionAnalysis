%% compute experiment data from scratch

sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'sensoryDependenceNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);  % remove not included sessions and empty lines
mice = unique(sessionInfo.mouse);

% sessionInfo = sessionInfo(1:6,:);  % !!! temp

data = cell(1,length(mice));
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data{1}.data = cellfun(@(x) x.data, data); data = data{1};
fprintf('saving...'); save(fullfile(getenv('SSD'), 'paper1', 'sensoryDependence_data.mat'), 'data'); disp('data saved!')



%% load experiment data
fprintf('loading... '); load(fullfile(getenv('SSD'), 'paper1', 'sensoryDependence_data.mat'), 'data'); disp('sensoryDependence data loaded!');  % temp, use to load from local drive because engram is slow over ethernet

% global settings
paper1_config;
varsToAvg = {'mouse'};
miceToExclude = {'sen1'};


% initializations
data.data = data.data(~ismember({data.data.mouse}, miceToExclude));

vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.sensoryCondition = struct('name', 'sensoryCondition', 'levels', {{'WL', 'W', 'L', '-'}}, 'levelNames', {{'W+V', 'W', 'V', '-'}});
vars.whiskers = struct('name', 'whiskers', 'levels', {{'full', 'none'}}, 'levelNames', {{'whiskers', 'no whiskers'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore paw', 'hind paw'}});
vars.isBigStep = struct('name', 'isBigStep', 'levels', [0 1], 'levelNames', {{'little step', 'big step'}});

conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
conditionals.none = struct('name', '', 'condition', @(x) x);
conditionals.modPawOnlySwing = struct('name', 'modPawOnlySwing', 'condition', @(x) x==1);


%% ----------
% PLOT THINGS
%  ----------


%% bar plots

% success
figure('position', [200 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
conditions = [vars.sensoryCondition];
mat = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg);
barFancy(mat, 'levelNames', {{'whiskers + vision', 'whiskers', 'vision', 'neither'}}, 'ylabel', 'success rate', ...
    'colors', sensColors, barProperties{:}, 'YLim', [0 1], 'YTick', [0 .5 1], ...
    'comparisons', [1 2; 1 3; 1 4], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceSuccessBars'), 'svg');


%% velocity at whisker contact
figure('position', [200 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
conditions = [vars.sensoryCondition];
mat = getDvMatrix(data, 'velAtWiskContact', conditions, varsToAvg);
barFancy(mat, 'levelNames', {conditions.levelNames}, 'ylabel', 'velocity at whisker contact (m/s)', ...
    'colors', sensColors, barProperties{:}, 'YTick', [.3 .5 .7], ...
    'comparisons', [1 2; 1 3; 1 4], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceVelBars'), 'svg');


%% step over height for all paws (measured before paw reaches obstacle)

figure('position', [200 472.00 600 328.00], 'color', 'white', 'menubar', 'none');
conditions = [vars.isFore; vars.isLeading; vars.sensoryCondition];
mat = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg) * 1000;
barFancy(mat, 'levelNames', {conditions.levelNames}, 'ylabel', 'step over height (mm)', ...
    'colors', repmat(sensColors,4,1), barProperties{:}, 'YTick', [4 10 16], ...
    'summaryFunction', @nanmean, 'constantEdgeColor', [.15 .15 .15], ...
    'comparisons', repmat([1 2; 1 3; 1 4],4,1) + repelem([0,4,8,12],3)', 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceStepHgts'), 'svg');


%% success of each paw
figure('position', [200 472.00 600 328.00], 'color', 'white', 'menubar', 'none');
conditions = [vars.isFore; vars.isLeading; vars.sensoryCondition];
mat = 1-getDvMatrix(data, 'anyTouchFrames', conditions, varsToAvg);
barFancy(mat, 'levelNames', {conditions.levelNames}, 'ylabel', 'success rate', 'YTick', [0 .5 1], ...
    'colors', repmat(sensColors,4,1), barProperties{:}, 'summaryFunction', @nanmean, 'constantEdgeColor', [.15 .15 .15], ...
    'comparisons', repmat([1 2; 1 3; 1 4],4,1) + repelem([0,4,8,12],3)', 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependencePawSuccess'), 'svg');


%% baseline step heights
figure('position', [2800.00 100 700 328.00], 'color', 'white', 'menubar', 'none');
figVars = [vars.isFore; vars.sensoryCondition];
dv = getDvMatrix(data, 'controlStepHgt', figVars, varsToAvg)*1000;
colorRepeats = prod(cellfun(@length, {figVars(1:end-1).levels}));
barFancy(dv, 'ylabel', 'control step height (mm)', 'levelNames', {figVars().levelNames}, 'colors', repmat(sensColors,colorRepeats,1), barProperties{:})
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependencePawSuccess'), 'svg');


%% speed vs. position

% settings
yLims = [.3 .7];
plotSequence = [4 3 2 1]; % determine which lines on plotted on top of which

% initializations
flat = flattenData(data, {'mouse', 'session', 'trial', 'sensoryCondition', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'velContinuousAtContact', 'velContinuousAtContactX', 'isWheelBreak', 'wiskContactPosition'});
colorsTemp = [sensColors(1:end-1,:); .6 .6 .6]; % the no vision no whisker condition can be a little lighter here

% speed vs position
figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [200 50 600 300], 'inverthardcopy', 'off')
x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
rectangle('Position', [x(1) yLims(1) diff(x) diff(yLims)], ...
    'FaceColor', [obsColor obsAlpha], 'EdgeColor', 'none');
line([0 0], yLims, 'linewidth', 2, 'color', get(gca, 'XColor'))
velData = plotDvPsth(flat, 'velVsPosition', 'sensoryCondition', ...
    {'showLegend', false, 'conditionColors', colorsTemp(plotSequence,:), 'xlim', [-.5 .2], ... 
     'plotConditions', vars.sensoryCondition.levels(plotSequence), 'errorAlpha', .1, 'lineWidth', 4, 'flipx', true});
set(gca, 'YLim', yLims, 'YTick', linspace(yLims(1),yLims(2),3));
xlabel('position relative to nose (m)')
ylabel('velocity (m/s)')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceVel');
saveas(gcf, file, 'svg');


% get vel at moment obstacle is under nose
atNoseInd = find(flat(1).velVsPositionX>=0,1,'first');
noseVels = squeeze(velData(:,:,atNoseInd));

for i = 1:length(vars.sensoryCondition.levels)
    fprintf('\n%s -> obs at nose vel: %.2f +- %.2f SEM', ...
        vars.sensoryCondition.levels{plotSequence(i)}, mean(noseVels(i,:)), std(noseVels(i,:))/sqrt(size(noseVels,2)))
end
fprintf('\n')


%% height shaping

% settings
xLims = [3 10];
isHgtPreObs = true; % do you measure the height at peak (false) or before the paw reaches the obstacle (true)

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isLeading', 'stepOverMaxHgt', 'sensoryCondition'});
flat = flat(~isnan([flat.obsHgt]) & ...
            ~isnan([flat.stepOverMaxHgt]) & ...
            ~isnan([flat.preObsHgt]));  % add conditionals here
if isHgtPreObs; pawHgts = [flat.preObsHgt]*1000; else; pawHgts = [flat.stepOverMaxHgt]*1000; end
obsHgts = [flat.obsHgt]*1000;
[~, conditions] = ismember({flat.sensoryCondition}, vars.sensoryCondition.levels);


% collect data
mice = unique({flat.mouse});
[corrs, slopes] = deal(nan(2, 2, 4, length(mice))); % isFore(10) X isLeading(10) X sensory condition X mouse
foreSequence = [true false];
leadingSequence = [true false];

for i = 1:2  % isFore
    for j = 1:2  % isLeading
        for c = 1:4  % sensory condition
            for k = 1:length(mice)
                bins = [flat.isFore]==foreSequence(i) & ...
                       [flat.isLeading]==leadingSequence(j) & ...
                       conditions==c & ...
                       strcmp({flat.mouse}, mice{k});
                x = obsHgts(bins);
                y = pawHgts(bins);
                fit = polyfit(x, y, 1);
                corrs(i,j,c,k) = corr(x', y');
                slopes(i,j,c,k) = fit(1);
            end
        end
    end
end



%% plot things

% all paws
figure('position', [200.00 472.00 600 328.00], 'color', 'white', 'menubar', 'none');
temp = [vars.isFore; vars.isLeading; vars.sensoryCondition];
barFancy(corrs, 'levelNames', {temp.levelNames}, 'ylabel', 'paw:hurdle correlation', 'YTick', [0 .4 .8], ...
    'colors', repmat(sensColors,4,1), 'YLim', [], barProperties{:}, 'constantEdgeColor', [.15 .15 .15], ...
    'comparisons', repmat([1 2; 1 3; 1 4],4,1) + repelem([0,4,8,12],3)', 'test', 'ttest')

file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceCorrsAllPaws');
saveas(gcf, file, 'svg');


%% leading forepaw only
figure('position', [200.00 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
barFancy(squeeze(corrs(1,1,:,:)), 'levelNames', {vars.sensoryCondition.levelNames}, 'ylabel', 'paw:hurdle correlation', ...
    'colors', repmat(sensColors,4,1), 'YLim', [-.2 .6], barProperties{:}, 'YTick', -.2:.2:.6, ...
    'comparisons', [1 2; 1 3; 1 4], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceCorrs'), 'svg');


%% moving averages
figure('Color', 'white', 'Position', [200 400 500 400], 'MenuBar', 'none');
plot([0 xLims(2)], [0 xLims(2)], 'Color', [obsColor .4], 'LineWidth', 3) % add unity line
lfBins = [flat.isLeading] & [flat.isFore];% leading forepaw bins
logPlotRick(obsHgts(lfBins), pawHgts(lfBins), ...
    'colors', sensColors, 'conditions', conditions(lfBins), 'xlabel', 'hurdle height (mm)', 'ylabel', 'paw height (mm)', 'plotMice', false, ...
    'xlim', [4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat(lfBins).mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)))
set(gca, 'xlim', [4 10], 'YTick', 4:2:12, 'YLim', [4 12])

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceHeightShapingMovingAvgs');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% kinematics

% settings
obsHgtBins = 4; % discretize obstacle heights into this many bins
fading = .25; % within a plot, each paw's color fades from fading*color -> color
xLims = [-.06 .0];
yLims = [0 .016];

% initializations
flat = flattenData(data, {'mouse', 'session', 'trial', 'sensoryCondition', 'isWheelBreak', ...
    'obsHgt', 'isFore', 'isLeading', 'preObsKin', 'controlStepKinInterp', 'stepOverKinInterp'});
obsHgtDiscretized = num2cell(discretize([flat.obsHgt], linspace(3.4, 10, obsHgtBins+1)/1000));
[flat.obsHgtDiscretized] = obsHgtDiscretized{:};
flat = flat(~isnan([flat.obsHgtDiscretized]) & ...
            ~[flat.isWheelBreak] & ...
            [flat.isLeading] & ...
            [flat.isFore]);
[~, conditions] = ismember({flat.sensoryCondition}, vars.sensoryCondition.levels);
kinData = permute(cat(3, flat.preObsKin), [3,1,2]);
% kinData = permute(cat(3, flat.stepOverKinInterp), [3,1,2]);
kinDataCtl = permute(cat(3, flat.controlStepKinInterp), [3,1,2]);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1); % change the x starting x position of ctl steps to match steps over


figure('position', [2000 200 400 700], 'color', 'white', 'menubar', 'none'); hold on;
colorsTemp = [sensColors(1:end-1,:); .8 .8 .8];
obsColors = repmat(obsColor, obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)';

for i = 1:4
    subplot(5,1,i)
    bins = conditions==i;
    plotColor = repmat(colorsTemp(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    
    plotKinematics(kinDataCtl(bins,[1,3],:), [flat(bins).obsHgt], ones(1,sum(bins)), ...
        'colors', ctlStepColor, 'lineAlpha', .8, 'plotObs', false) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    plotKinematics(kinData(bins,[1,3],:), [flat(bins).obsHgt], [flat(bins).obsHgtDiscretized], ...
        'colors', plotColor, 'obsAlpha', 1, 'lineAlpha', .8, 'mouseNames', {flat(bins).mouse}, 'obsColors', obsColors) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    set(gca, 'XLim', xLims, 'YLim', yLims)
end

subplot(5,1,5)
plotKinematics(kinData(:,[1,3],:), [flat.obsHgt], conditions, ...
    'colors', sensColors, 'obsAlpha', 1, 'lineAlpha', 1, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @nanstd, 'lineWidth', 3, 'obsColors', repmat(obsColor, obsHgtBins, 1)) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
set(gca, 'XLim', xLims, 'YLim', yLims)

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceKinematics');
saveas(gcf, file, 'svg');

%% decision making initializations
flat = flattenData(data, ...
    [m.predictors, {'mouse', 'modSwingContacts', 'isModPawLengthened', 'modPawDeltaLength', 'isBigStep', 'isLightOn', ...
    'modPawOnlySwing', 'isTrialSuccess', 'sensoryCondition', 'modPawPredictedDistanceToObs', 'modPawDistanceToObs', ...
    'modPawKinInterp', 'preModPawKinInterp', 'firstModPaw', 'preModPawDeltaLength', 'modPawXVel', 'modPawZVel', 'trialVel'}]);

%% heatmaps (figures f3h)
plotDecisionHeatmaps(flat, 'condition', 'sensoryCondition', 'levels', vars.sensoryCondition.levels, 'outcome', 'isModPawLengthened', ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'avgMice', false, 'plotMice', false, 'colors', sensColors, 'xLims', [-20 15], 'normalize', 'col', ...
    'plotProbs', false, 'subplotDims', [4 1], ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceHeatmaps'));

%% trial scatters (figures f3g)
% figures
close all
rng(2)
plotDecisionTrials(flat, 'condition', 'sensoryCondition', 'levels', vars.sensoryCondition.levels, 'outcome', 'isModPawLengthened', ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', decisionColors, 'showTitles', false, 'rowColors', sensColors, 'obsColor', obsColor, ...
    'showHistos', true, 'poolHistos', true, 'evenlySampleTrials', true, 'histosOnBottom', true, 'histoFillAlpha', 0, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceDecisionKin'));

%% model accuracies (figures s3f)
rng(1)
[~,~,temp] = plotModelAccuracies(flat, m.predictors, 'isModPawLengthened', ...
    'condition', 'sensoryCondition', 'levels', vars.sensoryCondition.levels, 'weightClasses', true, ...
    'modSwingContactsMax', m.modSwingContactsMax, 'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', sensColors, 'barProps', [barProperties, 'YLim', [.2 1], 'comparisons', [2 4; 2 6; 2 8; 4 10], 'test', 'ttest'], 'kFolds', 10, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceModels'), ...
    'modelTransfers', [2 3]);

%% landing position entropy (and landing position histogram)
plotEntropies(flat, 'condition', 'sensoryCondition', 'levels', vars.sensoryCondition.levels, ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'xLimsAbs', [-.11 .08], 'colors', sensColors, ...
    'barProps', [barProperties, {'comparisons', [1 2; 1 3; 1 4; 5 6; 5 7; 5 8]}], ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceEntropy'));
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceDistributions'), 'svg');

%% plot predictors
plotPredictors(flat, [m.predictors {'modPawPredictedDistanceToObs'}], 'isModPawLengthened', 'colors', sensColors, 'avgMice', true, ...
    'condition', 'sensoryCondition', 'levels', vars.sensoryCondition.levels, ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'mouseAlpha', .4);

%% bimodality
% rng(2)
close all
plotBimodalities(flat, 'bic', true, 'condition', 'sensoryCondition', 'levels', vars.sensoryCondition.levels, 'outcome', 'isModPawLengthened', ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', sensColors, 'barProps', [barProperties, 'comparisons', [2 4; 2 6; 2 8], 'test', 'ttest'], ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceBimodality'));

%% lagging forepaw planting distance variability

flat = flattenData(data, {'mouse', 'session', 'trial', 'sensoryCondition', 'laggingPenultKin'});
mice = unique({flat.mouse});
landingPositions = cellfun(@(x) x(1,end), {flat.laggingPenultKin});

distVar = nan(4, length(mice));  % variability in planting paw distance to hurdle // (sensory condition) X (mouse)
for i = 1:4
    for j = 1:length(mice)
        bins = strcmp({flat.sensoryCondition}, vars.sensoryCondition.levels{i}) & ...
            strcmp({flat.mouse}, mice{j});
        distVar(i,j) = nanstd(landingPositions(bins));
    end
end

figure('position', [200 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
barFancy(distVar*1000, 'levelNames', {vars.sensoryCondition.levelNames}, 'ylabel', 'planting distance variability (mm)', ...
    'colors', sensColors, barProperties{:}, 'YTick', 6:4:22, ...
    'comparisons', [2 3], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependencePlantingVariability'), 'svg');

%% distribution of planting landing positions

flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOff', 'obsHgt', ...
    'sensoryCondition', 'laggingPenultKin'});
plotLaggingKinematics(flat, 'condition', 'sensoryCondition', 'levels', vars.sensoryCondition.levels, 'colors', sensColors, ...
    'trialsToShow', 50, 'xLims', [-.15 .01], 'randSeed', 5, 'obsColor', obsColor, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependencePlantingKinematics'));


%% landing distance bar plots

% landing distance
figure('position', [200 472.00 300 255], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'modPawDistanceToObsAbs', [vars.isBigStep; vars.sensoryCondition], {'mouse'}, [conditionals.modPawOnlySwing])*1000;
barFancy(dv, 'ylabel', 'landing distance (mm)', 'levelNames', {vars.isBigStep.levelNames}, 'colors', repmat(sensColors,2,1), barProperties{:}, ...
    'comparisons', [1 2; 1 3; 1 4; 5 6; 5 7; 5 8; 2 3; 6 7], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependence_landingDistance'), 'svg');

%% mod paw deltas
close all

binEdges = linspace(-.05, .15, 50);
stepTypeBins = [flat.isBigStep]==1;

figure('color', 'white', 'position', [91.00 627.00 1099.00 402.00]); hold on
for i = vars.sensoryCondition.levels
    bins = strcmp({flat.sensoryCondition}, i{1});
    x = [flat(bins & stepTypeBins).modPawDeltaLength];
    histogram(x, binEdges)
    fprintf('%s: %.2f\n', i{1}, nanmean((x)*1000))
end

legend(vars.sensoryCondition.levelNames)








