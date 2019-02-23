%% -------------- 
% INITIALIZATIONS
% ---------------

% settings
varsToAvg = {'mouse'};

% load session metadata
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'whiskerTrimNotes');
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session),:); % remove empty rows, not included sessions, and those without correct brain region

% set categorical vars
vars.paw = struct('name', 'paw', 'levels', 1:4, 'levelNames', {{'LH', 'LF', 'RF', 'RH'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'lead', 'lag'}});
vars.isContra = struct('name', 'isContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
vars.isFore = struct('name', 'isFore', 'levels', [0 1], 'levelNames', {{'hind', 'fore'}});
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.isModPawContra = struct('name', 'isModPawContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
vars.condition = struct('name', 'condition', ...
    'levels', {{'bilateral full', 'unilateral full', 'unilateral int1', 'unilateral int2', 'unilateral int3', 'unilateral deltaOnly'}}, ...
    'levelNames', {{'biFull', 'unFull', 'un1', 'un2', 'un3', 'delta'}});

% set conditionals
conditionals.lightOff = struct('name', 'isLightOn', 'condition', @(x) x==0);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);

figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals

%% compute kinData for all sessions (only need to do once)
sessions = unique(sessionInfo(logical([sessionInfo.include]),:).session);
parfor i = 1:length(sessions); getKinematicData5(sessions{i}); end

%% compute experiment data
data = getExperimentData(sessionInfo, 'all');
save(fullfile(getenv('OBSDATADIR'), 'matlabData', [brainRegion '_' manipulation '_data.mat']), 'data');

%% load experiment data
load(fullfile(getenv('OBSDATADIR'), 'matlabData', [brainRegion '_' manipulation '_data.mat']), 'data');
disp([manipulation ' data loaded!'])


%% ----------
% PLOT THINGS
%  ----------


%% bar plots


% settings
close all
rowVar = 4;
colVar = 5;


% initializations
figure('name', 'whiskerTrim', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1800 900])
plotInd = 0;

% success
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = vars.condition;
dvMatrix = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'success rate', true)

% velocity
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'trialVel', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'velocity (m/s)', true)

% baseline step height
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd:plotInd+2); plotInd = plotInd+2;
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'baselineStepHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dvMatrix, {conditions.levelNames}, 'baseline step height (mm)', true)

% body angle
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'trialAngleContra', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'body angle (towards contra)', true)

% contra first rate
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'isContraFirst', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'contra paw first rate', true)

% paw success rate
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd:plotInd+2); plotInd = plotInd+2;
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'isPawSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'success rate', true)

% step over height (light, fore/hind, ipsi/contra, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd:plotInd+2); plotInd = plotInd+2;
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'step over anticipatory height', true)

% big step probability (light, modPawContra, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'isBigStep', conditions, varsToAvg, conditionals.isLagging);
barPlotRick(dvMatrix, {conditions.levelNames}, 'big step probability', true)

% wisk contact position (light, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'wiskContactPosition', conditions, varsToAvg, figConditionals)*-1000;
barPlotRick(dvMatrix, {conditions.levelNames}, {'obstacle distance', 'to nose at contact (mm)'}, true)

% height shaping !!! why is this different than before!
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd:plotInd+2); plotInd = plotInd+2;
%%
close all; figure
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, conditions, {'mouse'}, {'session'}, figConditionals); % perform regression for each session, then average slopes across sessions for each mouse
barPlotRick(dvMatrix, {conditions.levelNames}, 'height shaping', true)

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', 'whiskerTrim_bars.fig'))


%% log plots

flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', 'obsHgt', ...
    'modPawPredictedDistanceToObs', 'isBigStep', 'condition'});
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>3 & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions

% big step prob by predicted distance to obs (manip)
conditions = cellfun(@(x) find(ismember(vars.condition.levels,x)), {flat.condition});
mice = unique({flat.mouse});
figure('name', [brainRegion '_' manipulation], 'color', 'white', 'menubar', 'none', 'position', [2000 200 300*length(mice) 600])

for i = 1:length(mice)
    subplot(2,length(mice),i)
    bins = strcmp({flat.mouse}, mice{i});
    logPlotRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).isBigStep], ...
        {'predicted distance to obstacle (m)', 'big step probability'}, conditions(bins))
    title(mice{i})
end

subplot(2,length(mice),length(mice)+1:length(mice)*2)
logPlotRick([flat.modPawPredictedDistanceToObs], [flat.isBigStep], ...
    {'predicted distance to obstacle (m)', 'big step probability'}, conditions, vars.condition.levelNames)
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', [brainRegion '_' manipulation '_bigStepProbability.fig']))

%% big step kinematics

% settings
rowVar = 'modPawPredictedDistanceToObs';
numRows = 5;

flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'conditionNum', 'modPawKinInterp', 'preModPawKinInterp', 'isBigStep', 'isLightOn', ...
    rowVar, 'preModPawDeltaLength', 'modPawDeltaLength', 'obsHgt', 'condition', 'isTrialSuccess', 'isWheelBreak'});
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>3 & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions
flat = flat(~[flat.isLightOn]); % add conditionals here
controlBins = strcmp({flat.condition}, manipConditions{1});
manipBins = strcmp({flat.condition}, manipConditions{2});
lims = prctile([flat.(rowVar)], [5 95]);
rowInds = discretize([flat.(rowVar)], linspace(lims(1), lims(2), numRows+1));

plotBigStepKin(flat(controlBins), rowInds(controlBins));
set(gcf, 'Name', [brainRegion '_' manipulation '_control.fig'])
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation '_bigStepKinControl.fig']))

plotBigStepKin(flat(manipBins), rowInds(manipBins));
set(gcf, 'Name', [brainRegion '_' manipulation '_manip.fig'])
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', [brainRegion '_' manipulation '_bigStepKinManip.fig']))


%% sessions over time

% success, vel, body angle, baseline step height, 
dvs = {'isTrialSuccess', 'trialVel', 'trialAngleContra', 'isContraFirst', 'isBigStep', 'tailHgt'};
flat = getNestedStructFields(data, cat(2, {'mouse', 'session', 'trial', 'condition', 'sessionNum', 'conditionNum'}, dvs));
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>3 & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions
plotAcrossSessions2(flat, dvs);
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', [brainRegion '_' manipulation '_sessionsOverTime.fig']))


%% speed vs. position / time plots

flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'conditionNum', 'isLightOn', 'condition', 'obsOnPositions', ...
    'velContinuousAtContact', 'velVsPosition', 'isWheelBreak', 'wiskContactPosition'});
flat = flat(~[flat.isWheelBreak]);
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>3 & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions

% speed vs position, control vs manip
figure('name', [brainRegion '_' manipulation], 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 500], 'inverthardcopy', 'off')
plotDvPsth(flat, 'velVsPosition', [-.5 .2], 'condition', 'isLightOn')
conditions = unique({flat.condition});
plotTitles = {'light off', 'light on'};
for i = 1:length(conditions)
    subplot(length(conditions),1,i)
    line(repmat(nanmean([flat.obsOnPositions]),1,2), get(gca, 'YLim'), 'color', [.5 .5 .5])
    title(plotTitles{i})
end
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation 'speedVsPosition.fig']))

% speed vs time centered around whisker contract
figure('name', [brainRegion '_' manipulation], 'Color', 'white', 'MenuBar', 'none', 'Position', [2600 50 600 500], 'inverthardcopy', 'off')
plotDvPsth(flat, 'velContinuousAtContact', [-.5 .5], 'condition', 'isLightOn')
for i = 1:length(conditions)
    subplot(length(conditions),1,i)
    title(plotTitles{i})
end
xlabel('time relative to whisker contact (s)')
ylabel('velocity (m/s)')
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', [brainRegion '_' manipulation 'speedAroundContact.fig']))


%% height shaping scatters

figure('name', [brainRegion '_' manipulation], 'color', 'white', 'menubar', 'none', 'position', [2000 50 1000 900])

% settings
rowVar = vars.isFore;
colVar = vars.isLeading;
scatVar = vars.condition;
xLims = [2 10];
yLims = [0 20];

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'condition'});
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>3 & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions
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

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', [brainRegion '_' manipulation '_heightShaping.fig']))

%% heat maps

% settings
xLims = [-.03 .015];
yLims = [-.03 .03];

% predicted vs. actual mod paw distance
figure('name', [brainRegion '_' manipulation], 'color', 'white', 'menubar', 'none', 'position', [2000 100 700 900])
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'condition', 'trial', 'isLightOn', ...
    'modPawDistanceToObs', 'modPawPredictedDistanceToObs', 'isTrialSuccess', 'isModPawContra'});
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>3 & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions
validBins = true(size(flat));
mice = unique({flat.mouse});

for i = 1:length(vars.condition.levels)
    for j = 1:length(vars.isModPawContra.levels)
        subplot(2,2,i+(j-1)*2)
        bins = strcmp(vars.condition.levels{i}, {flat.condition}) & validBins & [flat.isModPawContra]==vars.isModPawContra.levels(j);

        heatmapRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).modPawDistanceToObs], ...
            {'predicted distance to ', 'actual distance'}, xLims, yLims); hold on
        plot(xLims, xLims, 'color', [.6 .6 1], 'LineWidth', 2)
        title(vars.condition.levelNames{i})
        if i==1; ylabel(vars.isModPawContra.levelNames{j}); end
    end
end
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation '_predictedDistanceHeatmaps.fig']))

% plot for individual mice
figure('name', [brainRegion '_' manipulation], 'color', 'white', 'menubar', 'none', 'position', [2000 100 300*length(mice) 600])
for mouse = 1:length(mice)
    for i = 1:length(vars.condition.levels)
        subplot(2,length(mice), mouse + (i-1)*length(mice))
        bins = strcmp(vars.condition.levels{i}, {flat.condition}) & validBins & strcmp({flat.mouse}, mice{mouse});

        heatmapRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).modPawDistanceToObs], ...
            {'predicted distance to ', 'actual distance'}, xLims, yLims); hold on
        plot(xLims, xLims, 'color', [.6 .6 1], 'LineWidth', 2)
        title(vars.condition.levelNames{i})
    end
    xlabel(mice{mouse})
end
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', [brainRegion '_' manipulation '_predictedDistanceHeatmapsMice.fig']))


% !!! predicted vs actual planting distance, one map per paw



%% kinematics

figure('name', [brainRegion '_' manipulation], 'color', 'white', 'menubar', 'none', 'position', [2000 566 1250 384])

% settings
rowVar = vars.isFore;
colVar = vars.isLeading;
plotVar = vars.condition;

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'isTrialSuccess', 'trial', 'isLightOn', 'isWheelBreak', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'condition', 'stepOverKinInterp', 'isBigStep'});
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>3 & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions
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

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', [brainRegion '_' manipulation 'kinematics.fig']))








