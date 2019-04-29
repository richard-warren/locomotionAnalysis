% load experiment data
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'sensoryDependence_data.mat'), 'data'); disp('baseline data loaded!')

%% global settings
colorWisk = [51 204 255]; %[255 204 51];
colorVision = [255 221 21]; %[51 204 255];
colorNone = [.2 .2 .2];
ctlColor = [.5 .5 .5];
varsToAvg = {'mouse'};
miceToExclude = {'sen1'};


% initializations
data.data = data.data(~ismember({data.data.mouse}, miceToExclude));
colors = [mean([colorWisk;colorVision],1); colorWisk; colorVision; colorNone]/255;

vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.sensoryCondition = struct('name', 'sensoryCondition', 'levels', {{'WL', 'W', 'L', '-'}}, 'levelNames', {{'W+V', 'W', 'V', '-'}});
vars.whiskers = struct('name', 'whiskers', 'levels', {{'full', 'none'}}, 'levelNames', {{'whiskers', 'no whiskers'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore paw', 'hind paw'}});

conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);




%% ----------
% PLOT THINGS
%  ----------


%% BARS

% initializations
rows = 2;
cols = 4;
figure('name', 'sensoryDependence', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1200 400])

% success
subplot(rows, cols, 1);
conditions = [vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg);
barPlotRick(dvMatrix, {'conditionNames', {conditions.levelNames}, 'ylabel', 'success rate', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colors, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [0 1], 'ytick', 0:.5:1, ...
    'compareToFirstOnly', true})

% velocity
subplot(rows, cols, 2);
conditions = [vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'trialVel', conditions, varsToAvg);
barPlotRick(dvMatrix, {'conditionNames', {conditions.levelNames}, 'ylabel', 'wheel velocity (m/s)', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colors, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [0 .75], 'ytick', 0:.25:.75, ...
    'compareToFirstOnly', false})

% paw error rate
subplot(rows, cols, 3:4);
conditions = [vars.isFore; vars.isLeading; vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'isPawSuccess', conditions, varsToAvg);
barPlotRick(dvMatrix, {'conditionNames', {conditions.levelNames}, 'ylabel', 'paw error rate', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', repmat(colors,4,1), 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [0 1], 'ytick', 0:.5:1, ...
    'compareToFirstOnly', false, 'lineWidth', 1.5})

% step over height for all paws
subplot(rows, cols, 5:6);
conditions = [vars.isFore; vars.isLeading; vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg) * 1000;
barPlotRick(dvMatrix, {'conditionNames', {conditions.levelNames}, 'ylabel', 'step over height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', repmat(colors,4,1), 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, ...
    'compareToFirstOnly', false})

% step over height for leading forelimbs
subplot(rows, cols, 7)
conditions = [vars.sensoryCondition];
figConditionals = [conditionals.isFore; conditionals.isLeading];
matMod = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg, figConditionals) * 1000;
matBl = getDvMatrix(data, 'baselineStepHgt', conditions, varsToAvg, figConditionals) * 1000;
dvMatrix = permute(cat(3,matBl,matMod), [1 3 2]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(ctlColor,8,1);
colorsWithBl(2:2:8,:) = colors;
barPlotRick(dvMatrix, {'conditionNames', {conditions.levelNames}, 'ylabel', {'leading fore paw', 'step over height (mm)'}, ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colorsWithBl, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, ...
    'compareToFirstOnly', false, 'lineWidth', 1.5})

% wisk contact position (light, manip)
subplot(rows, cols, 8);
conditions = [vars.sensoryCondition];
dvMatrix = abs(getDvMatrix(data, 'wiskContactPosition', conditions, varsToAvg)) * 1000;
barPlotRick(dvMatrix, {'conditionNames', {conditions.levelNames}, 'ylabel', {'distance to nose', 'at whisker contact (mm)'}, ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colors, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, ...
    'compareToFirstOnly', false})

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'sensoryDependenceBars');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');



%% log plots

flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsHgt', ...
    'modPawPredictedDistanceToObs', 'isBigStep', 'sensoryCondition'});
plotVar = vars.sensoryCondition;

% big step prob by predicted distance to obs (manip)
conditions = cellfun(@(x) find(ismember(plotVar.levels,x)), {flat.(plotVar.name)});
mice = unique({flat.mouse});
figure('name', 'sensoryDependence', 'color', 'white', 'menubar', 'none', 'position', [2000 200 300*length(mice) 250])

for i = 1:length(mice)
    subplot(1,length(mice),i)
    bins = strcmp({flat.mouse}, mice{i});
    logPlotRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).isBigStep], ...
        {'predicted distance to obstacle (m)', 'big step probability'}, conditions(bins))
    title(mice{i})
end
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'sensoryDependence', 'sensoryDependence_bigStepProbability_mice.fig'))

figure('name', 'sensoryDependence', 'color', 'white', 'menubar', 'none', 'position', [2000 200 500 300])
logPlotRick([flat.modPawPredictedDistanceToObs], [flat.isBigStep], ...
    {'predicted distance to obstacle (m)', 'big step probability'}, conditions, plotVar.levelNames)
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'sensoryDependence', 'sensoryDependence_bigStepProbability.fig'))


%% speed vs. position / time plots

yLims = [.25 .6];
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'sensoryCondition', 'isLightOn', 'obsOnPositions', ...
    'velContinuousAtContact', 'velVsPosition', 'isWheelBreak', 'wiskContactPosition', 'isBigStep'});
flat = flat(~[flat.isWheelBreak]);

% speed vs position, 
figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 700 600], 'inverthardcopy', 'off')
subplot(2,1,1)
plotDvPsth(flat, 'velVsPosition', [-.5 .2], 'sensoryCondition')
line(repmat(nanmean([flat.obsOnPositions]),1,2), yLims, 'color', [.5 .5 .5])
line([0 0], yLims, 'color', [.5 .5 .5])
set(gca, 'YLim', yLims);
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')

% speed vs time centered around whisker contract
subplot(2,1,2)
plotDvPsth(flat, 'velContinuousAtContact', [-.5 .5], 'sensoryCondition')
line([0 0], yLims, 'color', [.5 .5 .5])
set(gca, 'YLim', yLims);
xlabel('time relative to whisker contact (s)')
ylabel('velocity (m/s)')
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'sensoryDependence', 'sensoryDependence_speed.fig'))



%% height shaping scatters

figure('name', 'sensoryDependence', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1000 900])

% settings
rowVar = vars.isFore;
colVar = vars.isLeading;
scatVar = vars.sensoryCondition;
xLims = [3 10];
yLims = [0 20];

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'sensoryCondition'});
% flat = flat(strcmp({flat.mouse})); % add conditionals here
obsHgts = [flat.obsHgt]*1000;
pawHgts = [flat.preObsHgt]*1000;

% initializations
conditions = cellfun(@(x) find(ismember(scatVar.levels, x)), {flat.(scatVar.name)});
rows = length(rowVar.levels);
cols = length(colVar.levels);
        
plotInd = 1;
for i = 1:rows
    for j = 1:cols
        subplot(rows, cols, plotInd)
        bins = cellfun(@(x) isequal(x, rowVar.levels(i)), {flat.(rowVar.name)}) & ...
               cellfun(@(x) isequal(x, colVar.levels(j)), {flat.(colVar.name)});
        scatterPlotRick(cat(1,obsHgts(bins),pawHgts(bins)), {'obstacle height', 'paw height'}, conditions(bins), scatVar.levelNames)
        plotInd = plotInd+1;
        title(sprintf('%s, %s', rowVar.levelNames{i}, colVar.levelNames{j}))
        set(gca, 'XLim', xLims, 'YLim', yLims)
    end
end

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'sensoryDependence', 'sensoryDependence_heightShaping.fig'))

%% heat maps

% settings
rowVar = vars.isLightOn;
colVar = vars.whiskers;
xLims = [-.03 .015];
yLims = [-.03 .03];

% predicted vs. actual mod paw distance
figure('name', 'sensoryDependence', 'color', 'white', 'menubar', 'none', 'position', [2000 100 700 900])
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', 'sensoryCondition', 'whiskers', ...
    'modPawDistanceToObs', 'modPawPredictedDistanceToObs', 'isTrialSuccess'});
% flat = flat(~[flat.isTrialSuccess]); % set conditionals here
mice = unique({flat.mouse});
rows = length(rowVar.levels);
cols = length(colVar.levels);

plotInd = 1;
for i = 1:rows
    for j = 1:cols
        subplot(rows, cols, plotInd)
        bins = cellfun(@(x) isequal(x, rowVar.levels(i)), {flat.(rowVar.name)}) & ...
               cellfun(@(x) isequal(x, colVar.levels{j}), {flat.(colVar.name)});
        heatmapRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).modPawDistanceToObs], ...
            {'predicted distance to obs', 'actual distance'}, xLims, yLims); hold on
        plot(xLims, xLims, 'color', [.6 .6 1], 'LineWidth', 2)
        title(sprintf('%s, %s', rowVar.levelNames{i}, colVar.levelNames{j}))
        plotInd = plotInd + 1;
    end
end
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'sensoryDependence', 'sensoryDependence_predictedDistanceHeatmaps.fig'))

% !!! predicted vs actual planting distance, one map per paw



%% kinematics

figure('name', 'sensoryDependence', 'color', 'white', 'menubar', 'none', 'position', [2000 566 1250 384])

% settings
rowVar = vars.isFore;
colVar = vars.isLeading;
plotVar = vars.sensoryCondition;

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'isTrialSuccess', 'trial', 'isLightOn', 'isWheelBreak', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'condition', 'stepOverKinInterp', 'preObsKin', 'isBigStep', 'sensoryCondition'});
% flat = flat(~[flat.isBigStep]); % add conditionals here

% initializations
conditions = cellfun(@(x) find(ismember(plotVar.levels, x)), {flat.(plotVar.name)});
rows = length(rowVar.levels);
cols = length(colVar.levels);
% kinData = permute(cat(3, flat.stepOverKinInterp), [3,1,2]);
kinData = permute(cat(3, flat.preObsKin), [3,1,2]);
kinData = kinData(:,[1,3],:); % keep only x and z dimensions
        
plotInd = 1;
for i = 1:rows
    for j = 1:cols
        subplot(rows, cols, plotInd)
        bins = cellfun(@(x) isequal(x, rowVar.levels(i)), {flat.(rowVar.name)}) & ...
               cellfun(@(x) isequal(x, colVar.levels(j)), {flat.(colVar.name)});        
        plotKinematics(kinData(bins,:,:), [flat(bins).obsHgt], conditions(bins), plotVar.levelNames)
        plotInd = plotInd+1;
        title(sprintf('%s, %s', rowVar.levelNames{i}, colVar.levelNames{j}))
    end
end

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'sensoryDependence', 'sensoryDependence_kinematics.fig'))




