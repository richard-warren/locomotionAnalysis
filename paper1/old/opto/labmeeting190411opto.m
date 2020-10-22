%% -------------- 
% INITIALIZATIONS
% ---------------

% settings
% varsToAvg = {'mouse'};
varsToAvg = {'mouse'};
exp = 'cerInt'; % session, mtc, mtcHighPower, vermis, cerInt, alm (set to 'session' if you want to analyze only a single session)
session = '190331_005';
minSpeed = .35;
% sessions = '190327_004'; % set to 'all' to analyze all sessions


% load session metadata
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'optoNotes');
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session) & [sessionInfo.include],:); % remove empty rows, not included sessions

% initialize sessions
switch exp
    case 'mtc'
        sessions = {'190324_005', '190324_006', '190326_004', '190327_005', '190327_006', '190328_004', '190330_003', '190331_004', ...
                    '190401_002', '190401_003', '190402_002', '190403_002', '190403_003', '190404_001', '190406_002'}; % 190406_002 has dft coordinates
    case 'mtcHighPower'
        sessions = {'190401_002', '190401_003', '190402_002', '190403_002', '190403_003', '190404_001', '190406_002'}; % mtc opto
    case 'vermis'
        sessions = {'190325_005', '190326_005', '190327_004', '190329_006', '190402_004'};
    case 'cerInt'
        sessions = {'190331_005', '190402_003', '190404_002'};
    case 'alm'
        sessions = {'190407_001', '190407_002', '190408_001', '190409_001', '190409_002'};
    case 'session'
        sessions = {session};
        sessionInfo = sessionInfo(strcmp(sessionInfo.session, sessions),:);
        data = getExperimentData(sessionInfo, 'all');
end
sessionInfo = sessionInfo(ismember(sessionInfo.session, sessions), :); % only keep sessions to be analyzed
mice = unique(sessionInfo.mouse);
isSingleSession = isstr(sessions);


% set categorical vars
vars.paw = struct('name', 'paw', 'levels', 1:4, 'levelNames', {{'LH', 'LF', 'RF', 'RH'}});
vars.isOptoOn = struct('name', 'isOptoOn', 'levels', [0 1], 'levelNames', {{'no opto', 'opto'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'lead', 'lag'}});
vars.isContra = struct('name', 'isContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
vars.isFore = struct('name', 'isFore', 'levels', [0 1], 'levelNames', {{'hind', 'fore'}});
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.isModPawContra = struct('name', 'isModPawContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
vars.mouse = struct('name', 'mouse', 'levels', {mice}, 'levelNames', {mice});

% set conditionals
conditionals.lightOff = struct('name', 'isLightOn', 'condition', @(x) x==0);
conditionals.highSpeed = struct('name', 'trialVel', 'condition', @(x) x>.35);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals
% figConditionals = [conditionals.highSpeed];

isSided = strcmp(sessionInfo.side{1}, 'left') || strcmp(sessionInfo.side{1}, 'right');



%% load experiment data
load(fullfile(getenv('OBSDATADIR'), 'matlabData', [exp '_opto_data.mat']), 'data');
disp([exp ' opto data loaded!'])
% data = data(strcmp({data.mouse}, 'vgt6')); varsToAvg = {'session'}; % run this line to show bars for all sessions of a given mouse!

%% compute experiment data
loadOldData = true;
if exist('data', 'var') && loadOldData; data = getExperimentData(sessionInfo, 'all', data); else; data = getExperimentData(sessionInfo, 'all'); end
save(fullfile(getenv('OBSDATADIR'), 'matlabData', [exp '_opto_data.mat']), 'data'); disp('data saved')

%% compute experiment from scratch, in parallel
data = cell(1,length(mice));    
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data = cat(2,data{:});
save(fullfile(getenv('OBSDATADIR'), 'matlabData', [exp '_opto_data.mat']), 'data'); disp('all done!')

%% restrict analysis to single mouse
mouse = 'vgt6';
varsToAvg = {'session'};
data = data(strcmp({data.mouse}, mouse));

%% ----------
% PLOT THINGS
%  ----------

%% bar plots


% settings
rowVar = 4;
colVar = 4;


% initializations
figure('name', 'opto', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1800 900])
plotInd = 0;

% success (light, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isLightOn; vars.isOptoOn];
% conditions = [vars.isOptoOn];
dvMatrix = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'success rate', true)

% velocity
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isLightOn; vars.isOptoOn];
dvMatrix = getDvMatrix(data, 'trialVel', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'velocity (m/s)', true)

% body angle
if isSided
    plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
    conditions = [vars.isLightOn; vars.isOptoOn];
    dvMatrix = getDvMatrix(data, 'trialAngleContra', conditions, varsToAvg, figConditionals);
    barPlotRick(dvMatrix, {conditions.levelNames}, 'body angle (towards contra)', true)
end

% baseline step height
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
if isSided; conditions = [vars.isFore; vars.isContra; vars.isOptoOn]; else; conditions = [vars.isFore; vars.isOptoOn]; end
dvMatrix = getDvMatrix(data, 'baselineStepHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dvMatrix, {conditions.levelNames}, 'baseline step height (mm)', true)

% contra first rate (light, manip)
if isSided
    plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
    conditions = [vars.isOptoOn];
    dvMatrix = getDvMatrix(data, 'isContraFirst', conditions, varsToAvg, figConditionals);
    barPlotRick(dvMatrix, {conditions.levelNames}, 'contra paw first rate', true)
end

% penult step length (light, fore/hind, ipsi/contra, manip) - hgt, vel?
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd:plotInd);
conditions = [vars.isOptoOn; vars.isFore; vars.isLeading];
dvMatrix = getDvMatrix(data, 'penultStepLength', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'penultimate step length', true)

% paw error rate
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
if isSided; conditions = [vars.isFore; vars.isContra; vars.isOptoOn]; else; conditions = [vars.isFore; vars.isLeading; vars.isOptoOn]; end
dvMatrix = getDvMatrix(data, 'isPawSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'success rate', true)

% planting step distance (light, fore/hind, ipsi/contra, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
if isSided; conditions = [vars.isFore; vars.isContra; vars.isOptoOn]; else; conditions = [vars.isFore; vars.isLeading; vars.isOptoOn]; end
dvMatrix = getDvMatrix(data, 'stepOverStartingDistance', conditions, varsToAvg, [figConditionals; conditionals.isLagging])*-1000; % only take lagging paws
barPlotRick(dvMatrix, {conditions.levelNames}, 'planting foot distance (mm)', true)

% step over height (light, fore/hind, ipsi/contra, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
if isSided; conditions = [vars.isFore; vars.isContra; vars.isOptoOn]; else; conditions = [vars.isFore; vars.isLeading; vars.isOptoOn]; end
dvMatrix = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'step over anticipatory height', true)

% height shaping
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isFore; vars.isLeading; vars.isOptoOn];
dvMatrix = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, conditions, {'mouse'}, {'session'}, figConditionals); % perform regression for each session, then average slopes across sessions for each mouse
barPlotRick(dvMatrix, {conditions.levelNames}, 'height shaping', true)

% big step probability (light, modPawContra, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
if isSided; conditions = [vars.isModPawContra; vars.isOptoOn]; else; conditions = [vars.isOptoOn]; end
dvMatrix = getDvMatrix(data, 'isBigStep', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'big step probability', true)

% wisk contact position (light, manip)
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isOptoOn];
dvMatrix = getDvMatrix(data, 'wiskContactPosition', conditions, varsToAvg, figConditionals)*-1000;
barPlotRick(dvMatrix, {conditions.levelNames}, {'obstacle distance', 'to nose at contact (mm)'}, true)

% tail height
plotInd = plotInd+1; subplot(rowVar, colVar, plotInd);
conditions = [vars.isOptoOn];
dvMatrix = getDvMatrix(data, 'tailHgt', conditions, varsToAvg, figConditionals);
if strcmp(varsToAvg{1}, 'mouse')
    barPlotRick(dvMatrix, {conditions.levelNames}, 'tail height (m)', true, mice)
elseif strcmp(varsToAvg{1}, 'session')
    sessions = unique(squeeze(struct2cell(getNestedStructFields(data, 'session'))));
    barPlotRick(dvMatrix, {conditions.levelNames}, 'tail height (m)', true, sessions)
end


savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'opto', [exp '_opto_bars.fig']))


%% log plots

% settings
rows = 2;

% initializations
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', 'isOptoOn', 'obsHgt', ...
    'modPawPredictedDistanceToObs', 'isBigStep'});
conditions = [flat.isOptoOn]+1;
cols = ceil(length(mice)/rows);
figure('name', 'opto', 'color', 'white', 'menubar', 'none', 'position', [100 100 300*cols 200*rows])

for i = 1:length(mice)
    subplot(rows, cols,i)
    bins = strcmp({flat.mouse}, mice{i});
    logPlotRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).isBigStep], ...
        {'predicted distance to obstacle (m)', 'big step probability'}, conditions(bins))
    title(mice{i})
end
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'opto', 'opto_bigStepProbability_mice.fig'))

figure('name', 'opto', 'color', 'white', 'menubar', 'none', 'position', [2000 200 600 400])
logPlotRick([flat.modPawPredictedDistanceToObs], [flat.isBigStep], ...
    {'predicted distance to obstacle (m)', 'big step probability'}, conditions, vars.isOptoOn.levelNames)
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'opto', [exp '_opto_bigStepProbability.fig']))

%% speed vs. position / time plots

flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'conditionNum', 'isLightOn', 'isOptoOn', 'obsOnPositions', ...
    'velContinuousAtContact', 'velVsPosition', 'isWheelBreak', 'wiskContactPosition'});
flat = flat(~[flat.isWheelBreak]);
yLims = [.2 .8];

% speed vs position, control vs manip
figure('name', 'opto', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 500], 'inverthardcopy', 'off')

subplot(2,1,1)
plotDvPsth(flat, 'velVsPosition', [-.5 .2], 'isOptoOn')
line([0 0], yLims, 'color', [.5 .5 .5])
line(repmat(nanmean([flat.obsOnPositions]),1,2), yLims, 'color', [.5 .5 .5])
set(gca, 'YLim', yLims)
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')

subplot(2,1,2)
plotDvPsth(flat, 'velContinuousAtContact', [-.5 .5], 'isOptoOn')
line([0 0], yLims, 'color', [.5 .5 .5])
set(gca, 'YLim', yLims)
xlabel('time relative to whisker contact (s)')
ylabel('velocity (m/s)')    

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'opto', [exp '_opto_speed.fig']))


%% height shaping scatters

figure('name', 'opto', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1000 900])

% settings
if isSided; rowVar = vars.isContra; else; rowVar = vars.isFore; end
colVar = vars.isLeading;
scatVar = vars.isOptoOn;
xLims = [2 10];
yLims = [0 20];

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'isOptoOn'});
if ~isequal(rowVar, vars.isFore); flat = flat([flat.isFore]); end
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

savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'opto', [exp '_opto_heightShaping.fig']))

%% heat maps

% settings
if isSided; rowVar = vars.isModPawContra; else; rowVar = vars.isLightOn; end
colVar = vars.isOptoOn;
xLims = [-.03 .015];
yLims = [-.03 .03];

% predicted vs. actual mod paw distance
figure('name', 'opto', 'color', 'white', 'menubar', 'none', 'position', [2000 100 800 800])
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', 'sensoryCondition', 'isOptoOn', 'sessionNum', ...
    'modPawDistanceToObs', 'modPawPredictedDistanceToObs', 'isTrialSuccess', 'isModPawContra', 'conditionNum'});
rows = length(rowVar.levels);
cols = length(colVar.levels);

plotInd = 1;
for i = 1:rows
    for j = 1:cols
        subplot(rows, cols, plotInd)
        bins = cellfun(@(x) isequal(x, rowVar.levels(i)), {flat.(rowVar.name)}) & ...
               cellfun(@(x) isequal(x, colVar.levels(j)), {flat.(colVar.name)});
        if any(bins)
            heatmapRick([flat(bins).modPawPredictedDistanceToObs], [flat(bins).modPawDistanceToObs], ...
                {'predicted distance to obs', 'actual distance'}, xLims, yLims); hold on
            plot(xLims, xLims, 'color', [.6 .6 1], 'LineWidth', 2)
        end
        title(sprintf('%s, %s', rowVar.levelNames{i}, colVar.levelNames{j}))
        plotInd = plotInd + 1;
    end
end
savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'opto', [exp '_opto_heatMaps.fig']))


%% kinematics

% settings
% close all
rowVar = vars.isLeading;
colVar = vars.isFore;
plotVar = vars.isOptoOn;

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'isTrialSuccess', 'trial', 'isLightOn', 'isWheelBreak', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'isOptoOn', 'stepOverKinInterp', 'isBigStep', 'preObsKin'});
if isSided; figs = {'_ipsi', '_contra'}; else; figs = {''}; end % add conditionals here

% initializations
conditions = cellfun(@(x) find(ismember(plotVar.levels, x)), {flat.(plotVar.name)});
rows = length(rowVar.levels);
cols = length(colVar.levels);

% kinData = permute(cat(3, flat.stepOverKinInterp), [3,1,2]);
kinData = permute(cat(3, flat.preObsKin), [3,1,2]);
kinData = kinData(:,[1,3],:); % keep only x and z dimensions

for fig = 1:length(figs)
    
    figure('name', ['opto' figs{fig}], 'color', 'white', 'menubar', 'none', 'position', [2000 500*(fig-1)+100 1250 400])
    
    % select subset of figure trials
    if strcmp(figs{fig}, '_ipsi'); sideBins = logical(~[flat.isContra]);
    elseif strcmp(figs{fig}, '_contra'); sideBins = logical([flat.isContra]);
    else sideBins = true(1,size(kinData,1)); end
    
    plotInd = 1;
    for i = 1:rows
        for j = 1:cols
            subplot(rows, cols, plotInd)
            bins = cellfun(@(x) isequal(x, rowVar.levels(i)), {flat.(rowVar.name)}) & ...
                   cellfun(@(x) isequal(x, colVar.levels(j)), {flat.(colVar.name)}) & ...
                   sideBins;
            plotKinematics(kinData(bins,:,:), [flat(bins).obsHgt], conditions(bins), plotVar.levelNames)
            plotInd = plotInd+1;
            title(sprintf('%s, %s', rowVar.levelNames{i}, colVar.levelNames{j}))
        end
    end
    savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'opto', [exp '_opto_kinematics' figs{fig} '.fig']))
end

%% make example videos

% settings
% session = '190313_001';
lightTrialsToShow = 15;
overWriteVids = true;

for i = 1:length(data)
    for j = 1:length(data(i).sessions)
    
        condition = true(1,length([data(i).sessions(j).trials.isOptoOn]));
        isOptoOn = [data(i).sessions(j).trials.isOptoOn] & condition;
        optoTrials = find(isOptoOn);
        noOptoTrials = find(~isOptoOn);

        optoTrials = optoTrials(sort(randperm(length(optoTrials), lightTrialsToShow)));
        noOptoTrials = noOptoTrials(sort(randperm(length(noOptoTrials), lightTrialsToShow)));
        trials = sort([optoTrials, noOptoTrials]);

        file = fullfile(getenv('OBSDATADIR'), 'editedVid', 'opto', exp, ...
            [data(i).sessions(j).session '_' data(i).mouse '_' data(i).sessions(j).brainRegion '_' data(i).sessions(j).side '_' data(i).sessions(j).mW 'mW.avi']);
        if ~exist(file, 'file') || overWriteVids
            makeVidWisk(file, data(i).sessions(j).session, [-.45 .1], .15, trials, {'NO OPTO', 'OPTO'}, isOptoOn+1); % [-.05 .1] f
        end
    end
end





