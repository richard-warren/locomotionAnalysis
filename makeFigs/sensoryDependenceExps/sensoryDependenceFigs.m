%% -------------- 
% INITIALIZATIONS
% ---------------

% settings
varsToAvg = {'mouse'};

% load session metadata
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'sensoryDependenceNotes');
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session), :); % remove empty rows

% set categorical vars
vars.paw = struct('name', 'paw', 'levels', 1:4, 'levelNames', {{'LH', 'LF', 'RF', 'RH'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'lead', 'lag'}});
vars.isFore = struct('name', 'isFore', 'levels', [0 1], 'levelNames', {{'hind', 'fore'}});
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.sensoryCondition = struct('name', 'sensoryCondition', 'levels', {{'WL', 'W', 'L', '-'}}, 'levelNames', {{'WL', 'W', 'L', '-'}});
vars.whiskers = struct('name', 'whiskers', 'levels', {{'full', 'none'}}, 'levelNames', {{'whiskers', 'no whiskers'}});

% set conditionals
conditionals.lightOff = struct('name', 'isLightOn', 'condition', @(x) x==0);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);

figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals


%% compute kinData for all sessions (only need to do once)
sessions = unique(sessionInfo(logical([sessionInfo.include]),:).session);
parfor i = 1:length(sessions); getKinematicData5(sessions{i}); end

%% compute experiment data
data = getExperimentData(sessionInfo, 'all');
save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'sensoryDependence_data.mat'), 'data'); disp('done saving')

%% load experiment data
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'sensoryDependence_data.mat'), 'data');
disp('sensoryDependence data loaded!')


%% ----------
% PLOT THINGS
%  ----------


%% bar plots


% settings
rowVar = 4;
colVar = 3;


% initializations
figure('name', 'sensoryDependence', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1000 900])
plotInd = 0;

% success (light, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'success rate', true)

% velocity
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'trialVel', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'velocity (m/s)', true)

% baseline step height
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'baselineStepHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dvMatrix, {conditions.levelNames}, 'baseline step height (mm)', true)

% penult step length (light, fore/hind, ipsi/contra, manip) - hgt, vel?
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd:plotInd);
conditions = [vars.isFore; vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'penultStepLength', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'penultimate step length', true)

% paw error rate
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isFore; vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'isPawSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'success rate', true)

% planting step distance (light, fore/hind, ipsi/contra, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isFore; vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'stepOverStartingDistance', conditions, varsToAvg, [figConditionals; conditionals.isLagging])*-1000; % only take lagging paws
barPlotRick(dvMatrix, {conditions.levelNames}, 'planting foot distance (mm)', true)

% step over height (light, fore/hind, ipsi/contra, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isFore; vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'step over anticipatory height', true)

% height shaping
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isFore; vars.sensoryCondition];
dvMatrix = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, conditions, {'mouse'}, {'session'}, figConditionals); % perform regression for each session, then average slopes across sessions for each mouse
barPlotRick(dvMatrix, {conditions.levelNames}, 'height shaping', true)

% big step probability (light, modPawContra, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'isBigStep', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'big step probability', true)

% wisk contact position (light, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'wiskContactPosition', conditions, varsToAvg, figConditionals)*-1000;
barPlotRick(dvMatrix, {conditions.levelNames}, {'obstacle distance', 'to nose at contact (mm)'}, true)

% tail height
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'tailHgt', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'tail height (m)', true)

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'sensoryDependence', 'sensoryDependennce_bars.fig'))


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
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'condition', 'stepOverKinInterp', 'isBigStep', 'sensoryCondition'});
% flat = flat(~[flat.isBigStep]); % add conditionals here

% initializations
conditions = cellfun(@(x) find(ismember(plotVar.levels, x)), {flat.(plotVar.name)});
rows = length(rowVar.levels);
cols = length(colVar.levels);
kinData = permute(cat(3, flat.stepOverKinInterp), [3,1,2]);
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




