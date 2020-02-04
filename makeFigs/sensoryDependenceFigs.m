%% compute experiment data from scratch

sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'sensoryDependenceNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);  % remove not included sessions and empty lines
mice = unique(sessionInfo.mouse);

data = cell(1,length(mice));
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data{1}.data = cellfun(@(x) x.data, data); data = data{1};
fprintf('saving...'); save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'sensoryDependence_data.mat'), 'data'); disp('data saved!')



%% load experiment data
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'sensoryDependence_data.mat'), 'data'); disp('sensoryDependence data loaded!')

% global settings
global_config;
varsToAvg = {'mouse'};
miceToExclude = {'sen1'};


% initializations
data.data = data.data(~ismember({data.data.mouse}, miceToExclude));

vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.sensoryCondition = struct('name', 'sensoryCondition', 'levels', {{'WL', 'W', 'L', '-'}}, 'levelNames', {{'W+V', 'W', 'V', '-'}});
vars.whiskers = struct('name', 'whiskers', 'levels', {{'full', 'none'}}, 'levelNames', {{'whiskers', 'no whiskers'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore paw', 'hind paw'}});

conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
conditionals.none = struct('name', '', 'condition', @(x) x);


%% ----------
% PLOT THINGS
%  ----------


%% bar plots

% success
figure('position', [2000 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
conditions = [vars.sensoryCondition];
mat = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg);
barFancy(mat, 'levelNames', {conditions.levelNames}, 'ylabel', 'success rate', ...
    'colors', sensColors, barProperties{:}, 'YLim', [0 1])
set(gca, 'YTick', 0:.5:1, 'TickDir', 'out')

file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'sensoryDependenceSuccessBars');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% velocity at whisker contact
figure('position', [2200 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
conditions = [vars.sensoryCondition];
mat = getDvMatrix(data, 'velAtWiskContact', conditions, varsToAvg);
barFancy(mat, 'levelNames', {conditions.levelNames}, 'ylabel', 'velocity at whisker contact (m/s)', ...
    'colors', sensColors, barProperties{:})
set(gca, 'YTick', 0:.2:.8, 'TickDir', 'out')

file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'sensoryDependenceVelBars');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% step over height for all paws (measured before paw reaches obstacle)
figure('position', [2400 472.00 600 328.00], 'color', 'white', 'menubar', 'none');
conditions = [vars.isFore; vars.isLeading; vars.sensoryCondition];
mat = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg) * 1000;
barFancy(mat, 'levelNames', {conditions.levelNames}, 'ylabel', 'step over height (mm)', ...
    'colors', repmat(sensColors,4,1), barProperties{:})
set(gca, 'YTick', 0:4:16, 'TickDir', 'out')

file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'sensoryDependenceStepHgts');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% success of each paw
figure('position', [2400 472.00 600 328.00], 'color', 'white', 'menubar', 'none');
conditions = [vars.isFore; vars.isLeading; vars.sensoryCondition];
mat = 1-getDvMatrix(data, 'anyTouchFrames', conditions, varsToAvg);
barFancy(mat, 'levelNames', {conditions.levelNames}, 'ylabel', 'success rate', ...
    'colors', repmat(sensColors,4,1), barProperties{:})
set(gca, 'YTick', 0:.5:1, 'TickDir', 'out')

file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'sensoryDependencePawSuccess');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% baseline step heights
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
figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 300], 'inverthardcopy', 'off')
x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
rectangle('Position', [x(1) yLims(1) diff(x) diff(yLims)], ...
    'FaceColor', [obsColor obsAlpha], 'EdgeColor', 'none');
line([0 0], yLims, 'linewidth', 2, 'color', get(gca, 'XColor'))
velData = plotDvPsth(flat, 'velVsPosition', 'sensoryCondition', ...
    {'showLegend', false, 'conditionColors', colorsTemp(plotSequence,:), 'xlim', [-.5 .2], ... 
     'plotConditions', vars.sensoryCondition.levels(plotSequence), 'errorAlpha', .1, 'lineWidth', 4});
set(gca, 'YLim', yLims, 'YTick', linspace(yLims(1),yLims(2),3));
xlabel('position relative to nose (m)')
ylabel('velocity (m/s)')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'sensoryDependenceVel');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% get vel at moment obstacle is under nose
atNoseInd = find(flat(1).velVsPositionX>=0,1,'first');
noseVels = squeeze(velData(1,:,atNoseInd));

fprintf('\nobs at nose: %.2f +- %.2f SEM\n', mean(noseVels), std(noseVels)/sqrt(length(noseVels)))


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
figure('position', [2000.00 472.00 600 328.00], 'color', 'white', 'menubar', 'none');
temp = [vars.isFore; vars.isLeading; vars.sensoryCondition];
barFancy(corrs, 'levelNames', {temp.levelNames}, 'ylabel', 'paw obstacle correlation', ...
    'colors', repmat(sensColors,4,1), 'YLim', [], barProperties{:})
set(gca, 'YTick', -.2:.2:.8, 'TickDir', 'out', 'position', [.15 .0 .8 .9])

file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceCorrsAllPaws');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% leading forepaw only
figure('position', [2600.00 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
barFancy(squeeze(corrs(1,1,:,:)), 'levelNames', {vars.sensoryCondition.levelNames}, 'ylabel', 'leading forepaw obstacle correlation', ...
    'colors', repmat(sensColors,4,1), 'YLim', [-.2 .6], barProperties{:})
set(gca, 'YTick', -.2:.2:1, 'TickDir', 'out', 'position', [.15 .0 .8 .9])

file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceCorrs');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% moving averages
figure('Color', 'white', 'Position', [2800 400 500 400], 'MenuBar', 'none');
plot([0 xLims(2)], [0 xLims(2)], 'Color', obsColor, 'LineWidth', 3) % add unity line
lfBins = [flat.isLeading] & [flat.isFore];% leading forepaw bins
logPlotRick(obsHgts(lfBins), pawHgts(lfBins), ...
    {'colors', sensColors, 'conditions', conditions(lfBins), 'xlabel', 'obstacle height', 'ylabel', 'paw height (mm)', 'plotMice', false, ...
    'xlim', [3.4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat(lfBins).mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1))})
set(gca, 'xlim', [4 10], 'YTick', 4:2:12, 'TickDir', 'out')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceHeightShapingMovingAvgs');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% kinematics

% settings
close all
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
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'sensoryDependenceKinematics');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

%% decision making

flat = flattenData(data, ...
    [m.predictors, {'mouse', 'isModPawLengthened', 'modPawDeltaLength', 'isBigStep', 'isLightOn', 'modPawOnlySwing', 'isTrialSuccess', 'sensoryCondition', 'modPawPredictedDistanceToObs', 'modPawDistanceToObs', 'modPawKinInterp', 'preModPawKinInterp', 'firstModPaw'}]);

%% heatmaps
plotDecisionHeatmaps(flat, 'condition', 'sensoryCondition', 'levels', vars.sensoryCondition.levels, ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'avgMice', true, 'plotMice', true, 'colors', sensColors, 'outcome', 'isModPawLengthened', ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceHeatmaps'));

%% trials scatters
plotDecisionTrials(flat, 'condition', 'sensoryCondition', 'levels', vars.sensoryCondition.levels, 'outcome', 'isModPawLengthened', ...
    'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', decisionColors, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceDecisionKin'));

%% model accuracies
plotModelAccuracies_new(flat, m.predictors, 'isModPawLengthened', 'modelTransfers', [1 3], ...
    'condition', 'sensoryCondition', 'levels', vars.sensoryCondition.levels, ...
    'balanceClasses', false, 'weightClasses', true, 'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'shuffledConditions', [1 2], 'colors', sensColors, 'barProps', barProperties, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceModels'));

%% decision threshold
plotDecisionThresholds(flat, 'condition', 'sensoryCondition', 'levels', vars.sensoryCondition.levels, 'outcome', 'isModPawLengthened', ...
    'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', sensColors, 'barProps', barProperties, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'sensoryDependenceThresholds'));






