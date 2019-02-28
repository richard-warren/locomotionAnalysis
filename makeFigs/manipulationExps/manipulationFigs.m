% -------------- 
% INITIALIZATIONS
% ---------------

% settings
varsToAvg = {'mouse'};
manipulation = 'muscimol';
brainRegion = 'mtc';
maxEarlySession = 3;

% load session metadata
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', [manipulation 'Notes']);
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session) & ...
                          strcmp(sessionInfo.brainRegion, brainRegion),:); % remove empty rows, not included sessions, and those without correct brain region
mice = unique(sessionInfo.mouse);

% set categorical vars
vars.paw = struct('name', 'paw', 'levels', 1:4, 'levelNames', {{'LH', 'LF', 'RF', 'RH'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'lead', 'lag'}});
vars.isContra = struct('name', 'isContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
vars.isFore = struct('name', 'isFore', 'levels', [0 1], 'levelNames', {{'hind', 'fore'}});
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.isModPawContra = struct('name', 'isModPawContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
if strcmp(manipulation, 'muscimol')
    vars.condition = struct('name', 'condition', 'levels', {{'saline', 'muscimol'}}, 'levelNames', {{'sal', 'mus'}});
elseif strcmp(manipulation, 'lesion')
    vars.condition = struct('name', 'condition', 'levels', {{'pre', 'post'}}, 'levelNames', {{'pre', 'post'}});
end
manipConditions = vars.condition.levels;

% set conditionals
conditionals.lightOff = struct('name', 'isLightOn', 'condition', @(x) x==0);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
conditionals.isPre = struct('name', 'condition', 'condition', @(x) strcmp(x, 'pre'));
conditionals.isPost = struct('name', 'condition', 'condition', @(x) strcmp(x, 'post'));
conditionals.isEarly = struct('name', 'conditionNum', 'condition', @(x) x<=maxEarlySession);
conditionals.isLate = struct('name', 'conditionNum', 'condition', @(x) x>=5 & x<=8);

if strcmp(manipulation, 'lesion'); figConditionals = [conditionals.isEarly]; else; figConditionals = struct('name', '', 'condition', @(x) x); end

%% compute kinData for all sessions (only need to do once)
sessions = unique(sessionInfo(logical([sessionInfo.include]),:).session);
parfor i = 1:length(sessions); getKinematicData5(sessions{i}); end

%% compute experiment data

data = cell(1,length(mice));
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data = cat(2,data{:});
save(fullfile(getenv('OBSDATADIR'), 'matlabData', [brainRegion '_' manipulation '_data.mat']), 'data'); disp('data saved')

%% load experiment data
load(fullfile(getenv('OBSDATADIR'), 'matlabData', [brainRegion '_' manipulation '_data.mat']), 'data');
disp([manipulation ' data loaded!'])


%% ----------
% PLOT THINGS
%  ----------


%% bar plots


% settings
rowVar = 4;
colVar = 4;


% initializations
figure('name', [brainRegion '_' manipulation], 'color', 'white', 'menubar', 'none', 'position', [2000 50 1800 900])
plotInd = 0;

% success (light, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isLightOn; vars.condition];
dvMatrix = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'success rate', true)

% velocity
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isLightOn; vars.condition];
dvMatrix = getDvMatrix(data, 'trialVel', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'velocity (m/s)', true)

% body angle
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isLightOn; vars.condition];
dvMatrix = getDvMatrix(data, 'trialAngleContra', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'body angle (towards contra)', true)

% baseline step height
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'baselineStepHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dvMatrix, {conditions.levelNames}, 'baseline step height (mm)', true)

% contra first rate (light, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isLightOn; vars.condition];
dvMatrix = getDvMatrix(data, 'isContraFirst', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'contra paw first rate', true)

% penult step length (light, fore/hind, ipsi/contra, manip) - hgt, vel?
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd:plotInd);
conditions = [vars.condition; vars.isFore; vars.isLeading];
dvMatrix = getDvMatrix(data, 'penultStepLength', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'penultimate step length', true)

% paw error rate
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'isPawSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'success rate', true)

% planting step distance (light, fore/hind, ipsi/contra, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'stepOverStartingDistance', conditions, varsToAvg, [figConditionals; conditionals.isLagging])*-1000; % only take lagging paws
barPlotRick(dvMatrix, {conditions.levelNames}, 'planting foot distance (mm)', true)

% step over height (light, fore/hind, ipsi/contra, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'step over anticipatory height', true)

% height shaping
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isFore; vars.isLeading; vars.condition];
dvMatrix = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, conditions, {'mouse'}, {'session'}, figConditionals); % perform regression for each session, then average slopes across sessions for each mouse
barPlotRick(dvMatrix, {conditions.levelNames}, 'height shaping', true)

% big step probability (light, modPawContra, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isLightOn; vars.isModPawContra; vars.condition];
dvMatrix = getDvMatrix(data, 'isBigStep', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'big step probability', true)

% wisk contact position (light, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isLightOn; vars.condition];
dvMatrix = getDvMatrix(data, 'wiskContactPosition', conditions, varsToAvg, figConditionals)*-1000;
barPlotRick(dvMatrix, {conditions.levelNames}, {'obstacle distance', 'to nose at contact (mm)'}, true)

% tail height
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isLightOn; vars.condition];
dvMatrix = getDvMatrix(data, 'tailHgt', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'tail height (m)', true)


% compare early and late success
if strcmp(manipulation, 'lesion')
    dvs = {'isTrialSuccess', 'trialVel', 'tailHgt'};
    for dv = dvs
        plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
        conditions = [vars.isLightOn];
        pre = getDvMatrix(data, dv{1}, conditions, varsToAvg, conditionals.isPre);
        early = getDvMatrix(data, dv{1}, conditions, varsToAvg, [conditionals.isPost; conditionals.isEarly]);
        late = getDvMatrix(data, dv{1}, conditions, varsToAvg, [conditionals.isPost; conditionals.isLate]);

        dvMatrix = cat(length(conditions)+2, pre, early, late);
        dvMatrix = permute(dvMatrix, [1:length(conditions), length(conditions)+2, length(conditions)+1]); % swap last two dimensions
        barPlotRick(dvMatrix, {conditions.levelNames, {'pre', 'early', 'late'}}, dv{1}, true)
    end
end

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation '_bars.fig']))


%% log plots

% settings
rows = 2;

% initializations
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', 'obsHgt', ...
    'modPawPredictedDistanceToObs', 'isBigStep', 'condition'});
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>maxEarlySession & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions
conditions = cellfun(@(x) find(ismember(vars.condition.levels,x)), {flat.condition});
cols = ceil(length(mice)/rows);
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [100 100 300*cols 200*rows])

for i = 1:length(mice)
    subplot(rows, cols,i)
    bins = strcmp({flat.mouse}, mice{i});
    logPlotRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).isBigStep], ...
        {'predicted distance to obstacle (m)', 'big step probability'}, conditions(bins))
    title(mice{i})
end
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation '_bigStepProbability_mice.fig']))

figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 200 600 400])
logPlotRick([flat.modPawPredictedDistanceToObs], [flat.isBigStep], ...
    {'predicted distance to obstacle (m)', 'big step probability'}, conditions, vars.condition.levelNames)
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation '_bigStepProbability.fig']))

%% sessions over time

% success, vel, body angle, baseline step height, 
dvs = {'isTrialSuccess', 'trialVel', 'trialAngleContra', 'isContraFirst', 'isBigStep', 'tailHgt'};
flat = getNestedStructFields(data, cat(2, {'mouse', 'session', 'trial', 'condition', 'sessionNum', 'conditionNum'}, dvs));
plotAcrossSessions2(flat, dvs);
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation '_sessionsOverTime.fig']))


%% speed vs. position / time plots

flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'conditionNum', 'isLightOn', 'condition', 'obsOnPositions', ...
    'velContinuousAtContact', 'velVsPosition', 'isWheelBreak', 'wiskContactPosition'});
flat = flat(~[flat.isWheelBreak]);
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>maxEarlySession & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions
yLims = [.2 .6];

% speed vs position, control vs manip
figure('name', [brainRegion '_' manipulation], 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 500], 'inverthardcopy', 'off')
plotDvPsth(flat, 'velVsPosition', [-.5 .2], 'condition', 'isLightOn')
conditions = unique({flat.condition});
plotTitles = {'light off', 'light on'};
for i = 1:length(conditions)
    subplot(length(conditions),1,i)
    title(plotTitles{i})
    line(repmat(nanmean([flat.obsOnPositions]),1,2), yLims, 'color', [.5 .5 .5])
    line([0 0], yLims, 'color', [.5 .5 .5])
    set(gca, 'YLim', yLims)
end
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation '_speedVsPosition.fig']))

% speed vs time centered around whisker contract
figure('name', [brainRegion '_' manipulation], 'Color', 'white', 'MenuBar', 'none', 'Position', [2600 50 600 500], 'inverthardcopy', 'off')
plotDvPsth(flat, 'velContinuousAtContact', [-.5 .5], 'condition', 'isLightOn')
for i = 1:length(conditions)
    subplot(length(conditions),1,i)
    title(plotTitles{i})
    line([0 0], yLims, 'color', [.5 .5 .5])
    set(gca, 'YLim', yLims)
end
xlabel('time relative to whisker contact (s)')
ylabel('velocity (m/s)')
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation '_speedAroundContact.fig']))


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
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>maxEarlySession & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions
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

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation '_heightShaping.fig']))

%% heat maps

% settings
rowVar = vars.isModPawContra;
colVar = vars.condition;
xLims = [-.03 .015];
yLims = [-.03 .03];

% predicted vs. actual mod paw distance
figure('name', [brainRegion '_' manipulation], 'color', 'white', 'menubar', 'none', 'position', [2000 100 800 800])
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', 'sensoryCondition', 'condition', ...
    'modPawDistanceToObs', 'modPawPredictedDistanceToObs', 'isTrialSuccess', 'isModPawContra', 'conditionNum'});
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>maxEarlySession & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions
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
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation '_predictedDistanceHeatmaps.fig']))

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
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>maxEarlySession & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions
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

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation '_kinematics.fig']))








