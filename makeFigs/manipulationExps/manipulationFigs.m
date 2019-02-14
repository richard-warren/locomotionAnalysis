%% -------------- 
% INITIALIZATIONS
% ---------------

% settings
manipulation = 'muscimol';
brainRegion = 'mtc';
varsToAvg = {'mouse', 'session'};

% load session metadata
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', [manipulation 'Notes']);
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session) & ...
                          [sessionInfo.include]==1 & ...
                          strcmp(sessionInfo.brainRegion, brainRegion),:); % remove empty rows, not included sessions, and those without correct brain region

% set categorical vars
vars.paw = struct('name', 'paw', 'levels', 1:4, 'levelNames', {{'LH', 'LF', 'RF', 'RH'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'lead', 'lag'}});
vars.isContra = struct('name', 'isContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
vars.isFore = struct('name', 'isFore', 'levels', [0 1], 'levelNames', {{'hind', 'fore'}});
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.isModPawContra = struct('name', 'isModPawContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
if strcmp(manipulation, 'muscimol')
    vars.condition = struct('name', 'condition', 'levels', {{'saline', 'muscimol'}}, 'levelNames', {{'sal', 'mus'}});
else
    vars.condition = struct('name', 'condition', 'levels', {{'pre', 'post'}}, 'levelNames', {{'pre', 'post'}});
end

% set conditionals
conditionals.lightOff = struct('name', 'isLightOn', 'condition', @(x) x==0);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
conditionals.isPre = struct('name', 'condition', 'condition', @(x) strcmp(x, 'pre'));
conditionals.isPost = struct('name', 'condition', 'condition', @(x) strcmp(x, 'post'));
conditionals.isEarly = struct('name', 'conditionNum', 'condition', @(x) x<=3);
conditionals.isLate = struct('name', 'conditionNum', 'condition', @(x) (x>=5 & x<=8));

if strcmp(manipulation, 'lesion'); figConditionals = [conditionals.isEarly]; else; figConditionals = struct('name', '', 'condition', @(x) x); end

%% compute kinData for all sessions (only need to do once)
sessions = unique(sessionInfo.session);
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

flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', 'obsHgt', ...
    'modPawPredictedDistanceToObs', 'isBigStep', 'condition'});
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>3 & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions

% big step prob by predicted distance to obs (manip)
conditions = cellfun(@(x) find(ismember(vars.condition.levels,x)), {flat.condition});
figure('name', [brainRegion '_' manipulation], 'color', 'white', 'menubar', 'none', 'position', [2000 550 750 350])
logPlotRick([flat.modPawPredictedDistanceToObs], [flat.isBigStep], ...
    {'predicted distance to obstacle (m)', 'big step probability'}, conditions, vars.condition.levelNames)
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation '_bigStepProbability.fig']))

%% big step kinematics

% !!! how to show this for different conditions?


%% sessions over time

% success, vel, body angle, baseline step height, 

%% speed vs. position / time plots

% speed vs position, control vs manip


% speed vs time centered around whisker contract


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
xLims = [-.03 .015];
yLims = [-.03 .03];

% predicted vs actual planting distance, one map per paw
figure('name', [brainRegion '_' manipulation], 'color', 'white', 'menubar', 'none', 'position', [1962 459 682 386])
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'condition', 'trial', 'isLightOn', ...
    'modPawDistanceToObs', 'modPawPredictedDistanceToObs', 'isTrialSuccess'});
if strcmp(manipulation, 'lesion'); flat = flat(~([flat.conditionNum]>3 & strcmp({flat.condition}, 'post'))); end % only use first 3 lesion sessions
validBins = true(size(flat));

for i = 1:length(vars.condition.levels)
    subplot(1,2,i)
    bins = strcmp(vars.condition.levels{i}, {flat.condition}) & validBins;
    
    heatmapRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).modPawDistanceToObs], ...
        {'predicted distance to ', 'actual distance'}, xLims, yLims); hold on
    plot(xLims, xLims, 'color', [.6 .6 1], 'LineWidth', 2)
    title(vars.condition.levelNames{i})
    
end
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'manipulations', [brainRegion '_' manipulation '_predictedDistanceHeatmaps.fig']))

% delta mod step length vs predicted distance to obs (contra, manip)


%% kinematics

% fore/hind, contra, manip (leading only, light off only)


% fore/hind, leading, manip (contra only, light off only)













