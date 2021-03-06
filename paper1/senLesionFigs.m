%% -------------- 
% INITIALIZATIONS
% ---------------

% settings
varsToAvg = {'mouse'};
maxEarlySession = 3;
minObsHgt = .008;

% load session metadata
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'senLesionNotes');
sessionInfo = sessionInfo(179:end,:); % !!! this temporary hack that includes only the new senLesion mice
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session),:); % remove empty rows, not included sessions, and those without correct brain region
mice = unique(sessionInfo.mouse);

% set categorical vars
vars.paw = struct('name', 'paw', 'levels', 1:4, 'levelNames', {{'LH', 'LF', 'RF', 'RH'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'lead', 'lag'}});
vars.isContra = struct('name', 'isContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
vars.isFore = struct('name', 'isFore', 'levels', [0 1], 'levelNames', {{'hind', 'fore'}});
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.isModPawContra = struct('name', 'isModPawContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
vars.condition = struct('name', 'condition', 'levels', {{'preTrim', 'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}}, 'levelNames', {{'preTrim', 'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}});
vars.conditionSub = struct('name', 'condition', 'levels', {{'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}}, 'levelNames', {{'pre', 'postIpsi', 'postContra', 'postBi', 'noWisk'}});
vars.sessionNum = struct('name', 'sessionNum', 'levels', 1:100, 'levelNames', {cellfun(@num2str, num2cell(1:100), 'UniformOutput', false)});
manipConditions = vars.condition.levels;
vars.mouse = struct('name', 'mouse', 'levels', {mice}, 'levelNames', {mice});

% set conditionals
conditionals.lightOff = struct('name', 'isLightOn', 'condition', @(x) x==0);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
conditionals.isPre = struct('name', 'condition', 'condition', @(x) strcmp(x, 'pre'));
conditionals.isEarly = struct('name', 'conditionNum', 'condition', @(x) x<=maxEarlySession);
conditionals.isLate = struct('name', 'conditionNum', 'condition', @(x) x>=5 & x<=8);
conditionals.isObsHigh = struct('name', 'obsHgt', 'condition', @(x) x>minObsHgt);

% figConditionals = [conditionals.isObsHigh];
figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals


%% load experiment data
disp('loading...'); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'senLesion_data.mat'), 'data'); disp('senLesion data loaded!')

%% add new data to loaded data
data = getExperimentData(sessionInfo, 'all', data);
disp('saving data...'); save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'senLesion_data.mat'), 'data', '-v7.3'); disp('data saved');

%% compute experiment from scratch, in parallel
data = cell(1,length(mice));    
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data{1}.data = cellfun(@(x) x.data, data); data = data{1};
disp('saving data...'); save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'senLesion_data.mat'), 'data', '-v7.3'); disp('data saved');

%% ----------
% PLOT THINGS
%  ----------

%% sessions over time

% success, vel, body angle, baseline step height, 
dvs = {'isTrialSuccess', 'trialVel', 'trialAngleContra', 'isContraFirst', 'isBigStep', 'tailHgt'};
flat = flattenData(data, cat(2, {'mouse', 'session', 'trial', 'condition', 'sessionNum', 'conditionNum'}, dvs));
plotAcrossSessions2(flat, dvs);
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'senLesion', 'senLesion_sessionsOverTime.fig'))


%% bar plots


% settings
rowVar = 4;
colVar = 4;


% initializations
figure('name', 'senLesion', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1800 900])
plotInd = 0;

% success
% plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
figure
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg, figConditionals);
barFancy(dvMatrix, 'levelNames', {conditions.levelNames}, 'ylabel', 'success rate')

%% velocity
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'trialVel', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'velocity (m/s)', true)

% body angle
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'trialAngleContra', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'body angle (towards contra)', true)

% baseline step height
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'baselineStepHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dvMatrix, {conditions.levelNames}, 'baseline step height (mm)', true)

% contra first rate (light, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.condition];
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
conditions = [vars.isModPawContra; vars.condition];
dvMatrix = getDvMatrix(data, 'isBigStep', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'big step probability', true)

% wisk contact position (light, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'wiskContactPosition', conditions, varsToAvg, figConditionals)*-1000;
barPlotRick(dvMatrix, {conditions.levelNames}, {'obstacle distance', 'to nose at contact (mm)'}, true)

% tail height
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.condition];
dvMatrix = getDvMatrix(data, 'tailHgt', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'tail height (m)', true, mice)


% compare early and late success
% dvs = {'isTrialSuccess', 'trialVel', 'tailHgt'};
% for dv = dvs
%     plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
%     conditions = [vars.isLightOn];
%     pre = getDvMatrix(data, dv{1}, conditions, varsToAvg, conditionals.isPre);
%     early = getDvMatrix(data, dv{1}, conditions, varsToAvg, [conditionals.isPost; conditionals.isEarly]);
%     late = getDvMatrix(data, dv{1}, conditions, varsToAvg, [conditionals.isPost; conditionals.isLate]);
% 
%     dvMatrix = cat(length(conditions)+2, pre, early, late);
%     dvMatrix = permute(dvMatrix, [1:length(conditions), length(conditions)+2, length(conditions)+1]); % swap last two dimensions
%     barPlotRick(dvMatrix, {conditions.levelNames, {'pre', 'early', 'late'}}, dv{1}, true)
% end

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'senLesion', 'senLesion_bars.fig'))


%% log plots

% settings
rows = 2;

% initializations
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', 'obsHgt', ...
    'modPawPredictedDistanceToObs', 'isBigStep', 'condition'});
flat = flat(~([flat.conditionNum]>maxEarlySession & contains({flat.condition}, 'post'))); % only use first 3 lesion sessions
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
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'senLesion', 'senLesion_bigStepProbability_mice.fig'))

figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 200 600 400])
logPlotRick([flat.modPawPredictedDistanceToObs], [flat.isBigStep], ...
    {'predicted distance to obstacle (m)', 'big step probability'}, conditions, vars.condition.levelNames)
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'senLesion', 'senLesion_bigStepProbability.fig'))


%% speed vs. position / time plots

flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'conditionNum', 'isLightOn', 'condition', 'obsOnPositions', ...
    'velContinuousAtContact', 'velVsPosition', 'isWheelBreak', 'wiskContactPosition'});
flat = flat(~[flat.isWheelBreak]);
flat = flat(~([flat.conditionNum]>maxEarlySession & contains({flat.condition}, 'post'))); % only use first 3 lesion sessions
yLims = [.2 .8];

% speed vs position, control vs manip
figure('name', 'senLesion', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 500], 'inverthardcopy', 'off')

subplot(2,1,1)
plotDvPsth(flat, 'velVsPosition', [-.5 .2], 'condition')
line([0 0], yLims, 'color', [.5 .5 .5])
line(repmat(nanmean([flat.obsOnPositions]),1,2), yLims, 'color', [.5 .5 .5])
set(gca, 'YLim', yLims)
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')

subplot(2,1,2)
plotDvPsth(flat, 'velContinuousAtContact', [-.5 .5], 'condition')
line([0 0], yLims, 'color', [.5 .5 .5])
set(gca, 'YLim', yLims)
xlabel('time relative to whisker contact (s)')
ylabel('velocity (m/s)')    

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'senLesion', 'senLesion_speed.fig'))


%% height shaping scatters

figure('name', 'senLesion', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1000 900])

% settings
rowVar = vars.isFore;
colVar = vars.isLeading;
scatVar = vars.condition;
xLims = [2 10];
yLims = [0 20];

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'condition'});
flat = flat(~([flat.conditionNum]>maxEarlySession & contains({flat.condition}, 'post'))); % only use first 3 lesion sessions
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

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'senLesion', 'senLesion_heightShaping.fig'))

%% heat maps

% settings
rowVar = vars.isModPawContra;
colVar = vars.condition;
xLims = [-.03 .015];
yLims = [-.03 .03];

% predicted vs. actual mod paw distance
figure('name', 'senLesion', 'color', 'white', 'menubar', 'none', 'position', [2000 100 400*length(vars.condition.levels) 800])
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', 'sensoryCondition', 'condition', 'sessionNum', ...
    'modPawDistanceToObs', 'modPawPredictedDistanceToObs', 'isTrialSuccess', 'isModPawContra', 'conditionNum'});
% colVar.levels = colVar.levels(1:max([flat.sessionNum])); colVar.levelNames = colVar.levelNames(1:max([flat.sessionNum]));
flat = flat(~([flat.conditionNum]>maxEarlySession & contains({flat.condition}, 'post'))); % only use first 3 lesion sessions
rows = length(rowVar.levels);
cols = length(colVar.levels);

plotInd = 1;
for i = 1:rows
    for j = 1:cols
        subplot(rows, cols, plotInd)
        bins = cellfun(@(x) isequal(x, rowVar.levels(i)), {flat.(rowVar.name)}) & ...
               cellfun(@(x) isequal(x, colVar.levels{j}), {flat.(colVar.name)});
        if any(bins)
            heatmapRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).modPawDistanceToObs], ...
                {'predicted distance to obs', 'actual distance'}, xLims, yLims); hold on
            plot(xLims, xLims, 'color', [.6 .6 1], 'LineWidth', 2)
        end
        title(sprintf('%s, %s', rowVar.levelNames{i}, colVar.levelNames{j}))
        plotInd = plotInd + 1;
    end
end
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'senLesion', 'senLesion_predictedDistanceHeatmaps'))

%% mouse heat maps over time

% settings
% mouse = 'sen7';
rowVar = vars.mouse;
colVar = vars.sessionNum;
xLims = [-.03 .015];
yLims = [-.03 .03];

% predicted vs. actual mod paw distance
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', 'sensoryCondition', 'condition', 'sessionNum', ...
    'modPawDistanceToObs', 'modPawPredictedDistanceToObs', 'isTrialSuccess', 'isModPawContra', 'conditionNum'});
colVar.levels = colVar.levels(1:max([flat.sessionNum])); colVar.levelNames = colVar.levelNames(1:max([flat.sessionNum]));
% flat = flat(~([flat.conditionNum]>maxEarlySession & contains({flat.condition}, 'post'))); % only use first 3 lesion sessions
rows = length(rowVar.levels);
cols = length(colVar.levels);
figure('name', 'senLesion', 'color', 'white', 'menubar', 'none', 'position', [100 50 200*cols 400*rows])

plotInd = 1;
for i = 1:rows
    for j = 1:cols
        subplot(rows, cols, plotInd)
        bins = cellfun(@(x) isequal(x, rowVar.levels{i}), {flat.(rowVar.name)}) & ...
               cellfun(@(x) isequal(x, colVar.levels(j)), {flat.(colVar.name)});
        if any(bins)
            heatmapRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).modPawDistanceToObs], ...
                {'predicted distance to obs', 'actual distance'}, xLims, yLims); hold on
            plot(xLims, xLims, 'color', [.6 .6 1], 'LineWidth', 2)
            condition = unique({flat(bins).condition});
            title(sprintf('%s, %s', rowVar.levelNames{i}, condition{1}));
        end
        plotInd = plotInd + 1;
    end
end
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'senLesion', 'senLesion_predictedDistanceHeatmaps_mice'))




%% kinematics

figure('name', 'senLesion', 'color', 'white', 'menubar', 'none', 'position', [2000 566 1250 384])

% settings
rowVar = vars.isFore;
colVar = vars.isLeading;
plotVar = vars.condition;

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'isTrialSuccess', 'trial', 'isLightOn', 'isWheelBreak', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'condition', 'stepOverKinInterp', 'isBigStep', 'preObsKin'});
flat = flat(~([flat.conditionNum]>maxEarlySession & contains({flat.condition}, 'post'))); % only use first 3 lesion sessions
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

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'senLesion', 'senLesion_kinematics.fig'))


%% render experiment vids

trialsPerSession = 15;
flat = flattenData(data, {'mouse', 'session', 'condition', 'trial', 'velAtWiskContact'});
overWriteVids = false;

for i = 1:length(mice)

    mouseConditions = unique({flat(strcmp({flat.mouse}, mice{i})).condition});
    for j = 1:length(mouseConditions)
        conditionSessions = unique({flat(strcmp({flat.mouse}, mice{i}) & ...
                                         strcmp({flat.condition}, mouseConditions{j})).session});                                     
        session = conditionSessions{end};
        trials = sort(randperm(max([flat(strcmp({flat.session}, session)).trial]), trialsPerSession));
        
        file = fullfile(getenv('OBSDATADIR'), 'editedVid', 'senLesion', [mice{i} '_' session '_' mouseConditions{j}]);
        if ~exist(file, 'file') || overWriteVids; makeVidWisk(file, session, [-.05 .1], .15, trials); end        
    end    
end





