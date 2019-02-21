%% -------------- 
% INITIALIZATIONS
% ---------------

% settings
varsToAvg = {'mouse', 'session'};

% load session metadata
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'baselineNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);

% set categorical vars
vars.paw = struct('name', 'paw', 'levels', 1:4, 'levelNames', {{'LH', 'LF', 'RF', 'RH'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'lead', 'lag'}});
vars.isFore = struct('name', 'isFore', 'levels', [0 1], 'levelNames', {{'hind', 'fore'}});
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});

% set conditionals
conditionals.lightOff = struct('name', 'isLightOn', 'condition', @(x) x==0);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals


%% compute kinData for all sessions (only need to do once)
sessions = unique(sessionInfo.session);
parfor i = 1:length(sessions); getKinematicData5(sessions{i}); end


%% compute experiment data
data = getExperimentData(sessionInfo, 'all');
save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data');

%% load experiment data
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data');
disp('baseline data loaded!')


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
close all
rows = 3;

% initializations
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', 'obsHgt', ...
    'modPawPredictedDistanceToObs', 'isBigStep', 'condition'});
% flat = flat(~strcmp({flat.mouse}, 'mtc1')); % add conditionals here
conditions = discretize([flat.obsHgt], linspace(3.4, 10, 4)/1000);
mice = unique({flat.mouse});
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [100 100 300*cols 200*rows])
cols = ceil(length(mice)/rows);

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


%% big step kinematics

% settings
rowVar = 'modPawPredictedDistanceToObs';
numRows = 4;

flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'modPawKinInterp', 'preModPawKinInterp', 'isBigStep', 'isLightOn', ...
    rowVar, 'preModPawDeltaLength', 'modPawDeltaLength', 'obsHgt', 'isTrialSuccess', 'isWheelBreak'});
% flat = flat(~[flat.isLightOn]); % add conditionals here
lims = prctile([flat.(rowVar)], [5 95]);
rowInds = discretize([flat.(rowVar)], linspace(lims(1), lims(2), numRows+1));
plotBigStepKin(flat, rowInds);
set(gcf, 'Name', 'baseline')
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'baseline', 'baseline_bigStepKinControl.fig'))




