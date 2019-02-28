%% -------------- 
% INITIALIZATIONS
% ---------------

% settings
varsToAvg = {'mouse'};
miceToExclude = {'den12'};

% load session metadata
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'whiskerTrimNotes');
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session),:); % remove empty rows, not included sessions, and those without correct brain region
mice = unique(sessionInfo.mouse); mice = mice(~ismember(mice, miceToExclude));

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
vars.conditionSub = struct('name', 'condition', ...
    'levels', {{'bilateral full', 'unilateral full', 'unilateral int3', 'unilateral deltaOnly'}}, ...
    'levelNames', {{'biFull', 'unFull', 'un3', 'delta'}}); % only a subset of all whisker trim conditions
vars.mouse = struct('name', 'mouse', 'levels', {mice}, 'levelNames', {mice});

% set conditionals
conditionals.lightOff = struct('name', 'isLightOn', 'condition', @(x) x==0);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
conditionals.isHind = struct('name', 'isFore', 'condition', @(x) x==0);
conditionals.goodMiceOnly = struct('name', 'mouse', 'condition', @(x) ~ismember(x, miceToExclude));


figConditionals = [conditionals.lightOff; conditionals.goodMiceOnly]; % no conditionals

%% compute kinData for all sessions (only need to do once)
sessions = unique(sessionInfo(logical([sessionInfo.include]),:).session);
parfor i = 1:length(sessions); getKinematicData5(sessions{i}); end

%% compute experiment data
data = cell(1,length(mice));
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data = cat(2,data{:});
save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'whiskerTrim_data.mat'), 'data'); disp('done saving data!')

%% load experiment data
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'whiskerTrim_data.mat'), 'data');
disp('whiskerTrim data loaded!')


%% ----------
% PLOT THINGS
%  ----------

%% bar plots


% settings
rows = 4;
cols = 5;


% initializations
figure('name', 'whiskerTrim', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1500 900])
plotInd = 0;

% success
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = vars.condition;
dvMatrix = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'success rate', true)

% velocity
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'trialVel', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'velocity (m/s)', true)

% baseline step height
plotInd = plotInd+1; subplot(rows, cols, plotInd:plotInd+2); plotInd = plotInd+2;
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'baselineStepHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dvMatrix, {conditions.levelNames}, 'baseline step height (mm)', true)

% body angle
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'trialAngleContra', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'body angle (towards contra)', true)

% contra first rate
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'isContraFirst', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'contra paw first rate', true)

% paw success rate
plotInd = plotInd+1; subplot(rows, cols, plotInd:plotInd+2); plotInd = plotInd+2;
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'isPawSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'success rate', true)

% step over height (light, fore/hind, ipsi/contra, manip)
plotInd = plotInd+1; subplot(rows, cols, plotInd:plotInd+2); plotInd = plotInd+2;
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'step over anticipatory height', true)

% big step probability (light, modPawContra, manip)
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'isBigStep', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'big step probability', true)

% wisk contact position (light, manip)
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'wiskContactPosition', conditions, varsToAvg, figConditionals)*-1000;
barPlotRick(dvMatrix, {conditions.levelNames}, {'obstacle distance', 'to nose at contact (mm)'}, true)

% height shaping !!! why is this different than before?! should check with a more straightforward method
plotInd = plotInd+1; subplot(rows, cols, plotInd:plotInd+2);
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, conditions, {'mouse'}, {'session'}, figConditionals); % perform regression for each session, then average slopes across sessions for each mouse
barPlotRick(dvMatrix, {conditions.levelNames}, 'height shaping', true)

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', 'whiskerTrim_bars.fig'))


%% log plots

% settings
rows = 2;

% initializations
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', 'obsHgt', ...
    'modPawPredictedDistanceToObs', 'isBigStep', 'condition'});
flat = flat(~[flat.isLightOn] & ~ismember({flat.mouse}, miceToExclude));
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
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', 'whiskerTrim_bigStepProbability_mice.fig'))

figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 200 600 400])
logPlotRick([flat.modPawPredictedDistanceToObs], [flat.isBigStep], ...
    {'predicted distance to obstacle (m)', 'big step probability'}, conditions, vars.condition.levelNames)
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', 'whiskerTrim_bigStepProbability.fig'))


%% sessions over time

% success, vel, body angle, baseline step height, 
dvs = {'isTrialSuccess', 'trialVel', 'trialAngleContra', 'isContraFirst', 'isBigStep', 'wiskContactPosition'};
flat = getNestedStructFields(data, cat(2, {'mouse', 'session', 'trial', 'condition', 'sessionNum', 'conditionNum', 'isLightOn'}, dvs));
flat = flat(~[flat.isLightOn] & ~ismember({flat.mouse}, miceToExclude));
plotAcrossSessions2(flat, dvs);
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', 'whiskerTrim_sessionsOverTime.fig'))


%% speed vs. position / time plots

yLims = [.25 .6];
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'condition', 'isLightOn', 'obsOnPositions', ...
    'velContinuousAtContact', 'velVsPosition', 'isWheelBreak', 'wiskContactPosition'});
flat = flat(~[flat.isWheelBreak] & ~[flat.isLightOn] & ...
    ismember({flat.condition}, vars.conditionSub.levels) & ~ismember({flat.mouse}, miceToExclude));


% speed vs position, 
figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 700 600], 'inverthardcopy', 'off')
subplot(2,1,1)
plotDvPsth(flat, 'velVsPosition', [-.5 .2], 'condition')
line(repmat(nanmean([flat.obsOnPositions]),1,2), yLims, 'color', [.5 .5 .5])
line([0 0], yLims, 'color', [.5 .5 .5])
set(gca, 'YLim', yLims);
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')

% speed vs time centered around whisker contract
subplot(2,1,2)
plotDvPsth(flat, 'velContinuousAtContact', [-.5 .5], 'condition')
line([0 0], yLims, 'color', [.5 .5 .5])
set(gca, 'YLim', yLims);
xlabel('time relative to whisker contact (s)')
ylabel('velocity (m/s)')
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', 'whiskerTrim_speed.fig'))


%% height shaping scatters

figure('name', 'whiskerTrim', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1000 900])

% settings
rowVar = vars.isFore;
colVar = vars.isContra;
scatVar = vars.conditionSub;
xLims = [2 10];
yLims = [0 20];

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'condition'});
flat = flat(~[flat.isLightOn] & ismember({flat.condition}, vars.conditionSub.levels) & ~ismember({flat.mouse}, miceToExclude));

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

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', 'whiskerTrim_heightShaping.fig'))

%% heat maps

% settings
rowVar = vars.isModPawContra;
colVar = vars.condition;
xLims = [-.03 .015];
yLims = [-.03 .03];

% predicted vs. actual mod paw distance
figure('name', 'sensoryDependence', 'color', 'white', 'menubar', 'none', 'position', [2000 100 1600 800])
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', 'sensoryCondition', 'condition', ...
    'modPawDistanceToObs', 'modPawPredictedDistanceToObs', 'isTrialSuccess', 'isModPawContra'});
% flat = flat(~[flat.isLightOn] & ismember({flat.condition}, vars.conditionSub.levels) & ~ismember({flat.mouse}, miceToExclude));
flat = flat(~[flat.isLightOn] & ~ismember({flat.mouse}, miceToExclude));
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
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', 'whiskerTrim_predictedDistanceHeatmaps.fig'))

% !!! predicted vs actual planting distance, one map per paw



%% kinematics

figure('name', 'sensoryDependence', 'color', 'white', 'menubar', 'none', 'position', [2000 566 1250 384])

% settings
rowVar = vars.isContra;
colVar = vars.isLeading;
plotVar = vars.conditionSub;

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'isTrialSuccess', 'trial', 'isLightOn', 'isWheelBreak', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'condition', 'stepOverKinInterp', 'isBigStep', 'sensoryCondition'});
flat = flat(~[flat.isLightOn] & ismember({flat.condition}, vars.conditionSub.levels) & ...
    ~ismember({flat.mouse}, miceToExclude) & [flat.isFore]);

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

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'whiskerTrim', 'whiskerTrim_kinematics.fig'))




