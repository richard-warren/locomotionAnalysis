%% -------------- 
% INITIALIZATIONS
% ---------------

% settings
varsToAvg = {'mouse'}; % {'mouse', 'session'}
obsHgtBins = 5; % discretize obstacle heights into this many bins

% load session metadata
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'baselineNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);
mice = unique(sessionInfo.mouse);

% set categorical vars
vars.isTrialSuccess = struct('name', 'isTrialSuccess', 'levels', [1 0], 'levelNames', {{'success', 'not success'}});
vars.paw = struct('name', 'paw', 'levels', 1:4, 'levelNames', {{'LH', 'LF', 'RF', 'RH'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'lead', 'lag'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore', 'hind'}});
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.obsHgtDiscretized = struct('name', 'obsHgtDiscretized', 'levels', 1:obsHgtBins, 'levelNames', {num2cell(1:obsHgtBins)});

% set conditionals
conditionals.lightOff = struct('name', 'isLightOn', 'condition', @(x) x==0);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals


%% recompute kinData

for i = 1:length(mice)
    sessions = sessionInfo.session(strcmp(sessionInfo.mouse, mice{i}));
    for j = 1:length(sessions); getKinematicData5(sessions{j}); end
end


%% compute experiment data
data = cell(1,length(mice));
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data{1}.data = cellfun(@(x) x.data, data); data = data{1};
fprintf('saving...'); save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('data saved!')

%% load experiment data
fprintf('loading...'); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')


%% ----------
% PLOT THINGS
%  ----------

%% bar plots

% settings
rows = 2;
cols = 2;


% initializations
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 50 700 500])
plotInd = 0;

% penult step length
plotInd = plotInd+1; subplot(rows, cols, plotInd:plotInd);
conditions = [vars.isFore; vars.isLeading];
dvMatrix = getDvMatrix(data, 'penultStepLength', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'penultimate step length', true)

% step over starting distance
plotInd = plotInd+1; subplot(rows, cols, plotInd:plotInd);
conditions = [vars.isFore; vars.isLeading];
dvMatrix = getDvMatrix(data, 'stepOverStartingDistance', conditions, varsToAvg)*-1000;
barPlotRick(dvMatrix, {conditions.levelNames}, 'planting foot distance (mm)', true)

% step over height
plotInd = plotInd+1; subplot(rows, cols, plotInd:plotInd);
conditions = [vars.isFore];
dvMatrix = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'step over anticipatory height', true)

% height shaping
plotInd = plotInd+1; subplot(rows, cols, plotInd:plotInd);
conditions = [vars.isFore; vars.isLeading];
dvMatrix = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, conditions, {'mouse'}, {'session'}, figConditionals); % perform regression for each session, then average slopes across sessions for each mouse
barPlotRick(dvMatrix, {conditions.levelNames}, 'height shaping', true)

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'baseline', 'baseline_bars.fig'))


%% big step prob by predicted distance to obs

% settings
rows = 3;

% initializations
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', 'obsHgt', ...
    'modPawPredictedDistanceToObs', 'isBigStep', 'condition'});
% flat = flat(~strcmp({flat.mouse}, 'mtc1')); % add conditionals here
conditions = discretize([flat.obsHgt], linspace(3.4, 10, 4)/1000);
mice = unique({flat.mouse});
cols = ceil(length(mice)/rows);
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [100 100 300*cols 200*rows])

for i = 1:length(mice)
    subplot(rows, cols,i)
    bins = strcmp({flat.mouse}, mice{i});
    logPlotRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).isBigStep], ...
        {'predicted distance to obstacle (m)', 'big step probability'}, conditions(bins))
    title(mice{i})
end
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'baseline', 'baseline_bigStepProbability_mice.fig'))

figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 200 600 400])
logPlotRick([flat.modPawPredictedDistanceToObs], [flat.isBigStep], ...
    {'predicted distance to obstacle (m)', 'big step probability'}, conditions, {'low', 'medium', 'high'})
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'baseline', 'baseline_bigStepProbability.fig'))


%% speed vs. position / time plots

yLims = [.25 .55];
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', ...
    'velContinuousAtContact', 'velContinuousAtContactX', 'velVsPosition', 'velVsPositionX', ...
    'isWheelBreak', 'wiskContactPosition', 'isBigStep'});
flat = flat(~[flat.isWheelBreak]);

% speed vs position, 
figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 500], 'inverthardcopy', 'off')
subplot(2,1,1)
plotDvPsth(flat, 'velVsPosition', [-.5 .2], 'isLightOn')
line(repmat(nanmean([flat.obsOnPositions]),1,2), yLims, 'color', [.5 .5 .5])
line([0 0], yLims, 'color', [.5 .5 .5])
set(gca, 'YLim', yLims);
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')

% speed vs time centered around whisker contract
subplot(2,1,2)
plotDvPsth(flat, 'velContinuousAtContact', [-.5 .5], 'isLightOn')
line([0 0], yLims, 'color', [.5 .5 .5])
set(gca, 'YLim', yLims);
xlabel('time relative to whisker contact (s)')
ylabel('velocity (m/s)')
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'baseline', 'baseline_speed.fig'))

%% height shaping scatters

figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1000 900])

% settings
rowVar = vars.isFore;
colVar = vars.isLeading;
scatVar = vars.isLightOn;
xLims = [2 10];
yLims = [0 20];

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading'});
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

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'baseline', 'baseline_heightShaping.fig'))

%% heat maps

% settings
rowVar = vars.isLightOn;
colVar = vars.isTrialSuccess;
xLims = [-.03 .015];
yLims = [-.03 .03];

% predicted vs. actual mod paw distance
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 100 700 900])
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsHgt', ...
    'modPawDistanceToObs', 'modPawPredictedDistanceToObs', 'isTrialSuccess', 'wiskContactPosition'});
% flat = flat(~[flat.isTrialSuccess]); % set conditionals here
mice = unique({flat.mouse});
rows = length(rowVar.levels);
cols = length(colVar.levels);

plotInd = 1;
for i = 1:rows
    for j = 1:cols
        subplot(rows, cols, plotInd)
        bins = cellfun(@(x) isequal(x, rowVar.levels(i)), {flat.(rowVar.name)}) & ...
               cellfun(@(x) isequal(x, colVar.levels(j)), {flat.(colVar.name)});
        heatmapRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).modPawDistanceToObs], ...
            {'xLims', xLims}); hold on
        plot(xLims, xLims, 'color', [.6 .6 1], 'LineWidth', 2)
        title(sprintf('%s, %s', rowVar.levelNames{i}, colVar.levelNames{j}))
        plotInd = plotInd + 1;
    end
end
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'baseline', 'baseline_predictedDistanceHeatmaps.fig'))

% plot for individual mice
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 100 300*cols 900])
rows = 4;
cols = ceil(length(mice)/rows);

for i = 1:length(mice)
    subplot(rows, cols, i)
    bins = strcmp({flat.mouse}, mice{i});
    heatmapRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).modPawDistanceToObs], ...
        {'predicted distance to obs', 'actual distance'}, xLims, yLims); hold on
    plot(xLims, xLims, 'color', [.6 .6 1], 'LineWidth', 2)
    title(mice{i})
end
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'baseline', 'baseline_predictedDistanceHeatmaps_mice.fig'))


% whisker contact position vs obstacle height
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 100 300 400])
heatmapRick([flat.obsHgt]*1000, [flat.wiskContactPosition]*1000, ...
        {'obstacle height (mm)', 'distance from nose at conact (mm)'}); hold on
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'baseline', 'baseline_obsHgtContactPosHeatmaps_mice.fig'))



%% kinematics

figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 566 1250 384])

% settings
rowVar = vars.isFore;w
colVar = vars.isLeading;
plotVar = vars.obsHgtDiscretized;

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'isTrialSuccess', 'trial', 'isLightOn', 'isWheelBreak', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'stepOverKinInterp', 'isBigStep'});
% flat = flat(~isnan([flat.obsHgt]));
obsHgtDiscretized = num2cell(discretize([flat.obsHgt], linspace(3.4, 10, obsHgtBins+1)/1000));
[flat.obsHgtDiscretized] = obsHgtDiscretized{:};
if isequal(plotVar, vars.obsHgtDiscretized); flat = flat(~isnan([flat.obsHgtDiscretized])); end


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
        plotKinematics(kinData(bins,:,:), [flat(bins).obsHgt], conditions(bins))
        plotInd = plotInd+1;
        title(sprintf('%s, %s', rowVar.levelNames{i}, colVar.levelNames{j}))
    end
end

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'baseline', 'baseline_kinematics.fig'))






