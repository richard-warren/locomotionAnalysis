%% -------------- 
% INITIALIZATIONS
% ---------------

% settings
manipulation = 'muscimol';
brainRegion = 'mtc';

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
save(fullfile(getenv('OBSDATADIR'), 'matlabData', [brainRegion '_' manipulation '_data.mat']), 'data')

%% load experiment data
load(fullfile(getenv('OBSDATADIR'), 'matlabData', [brainRegion '_' manipulation '_data.mat']), 'data');
disp([manipulation ' data loaded!'])


%% ----------
% PLOT THINGS
%  ----------

%% bar plots


% settings
varsToAvg = {'mouse', 'session'};
rows = 4;
cols = 4;


% initializations
figure('name', [brainRegion '_' manipulation], 'color', 'white', 'menubar', 'none', 'position', [2000 50 1800 900])
plotInd = 0;

% success (light, manip)
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.isLightOn; vars.condition];
dvMatrix = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'success rate', true)

% velocity
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.isLightOn; vars.condition];
dvMatrix = getDvMatrix(data, 'trialVel', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'velocity (m/s)', true)

% body angle
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.isLightOn; vars.condition];
dvMatrix = getDvMatrix(data, 'trialAngleContra', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'body angle (towards contra)', true)

% baseline step height
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'baselineStepHgt', conditions, varsToAvg, figConditionals) * 1000;
barPlotRick(dvMatrix, {conditions.levelNames}, 'baseline step height (mm)', true)

% contra first rate (light, manip)
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.isLightOn; vars.condition];
dvMatrix = getDvMatrix(data, 'isContraFirst', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'contra paw first rate', true)

% penult step length (light, fore/hind, ipsi/contra, manip) - hgt, vel?
plotInd = plotInd+1; subplot(rows, cols, plotInd:plotInd);
conditions = [vars.condition; vars.isFore; vars.isLeading];
dvMatrix = getDvMatrix(data, 'penultStepLength', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'penultimate step length', true)

% paw error rate
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'isPawSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'success rate', true)

% planting step distance (light, fore/hind, ipsi/contra, manip)
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'stepOverStartingDistance', conditions, varsToAvg, [figConditionals; conditionals.isLagging])*-1000; % only take lagging paws
barPlotRick(dvMatrix, {conditions.levelNames}, 'planting foot distance (mm)', true)

% step over height (light, fore/hind, ipsi/contra, manip)
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.isFore; vars.isContra; vars.condition];
dvMatrix = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'step over anticipatory height', true)

% height shaping
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.isFore; vars.isLeading; vars.condition];
dvMatrix = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, conditions, {'mouse'}, {'session'}, figConditionals); % perform regression for each session, then average slopes across sessions for each mouse
barPlotRick(dvMatrix, {conditions.levelNames}, 'height shaping', true)

% big step probability (light, modPawContra, manip)
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.isLightOn; vars.isModPawContra; vars.condition];
dvMatrix = getDvMatrix(data, 'isBigStep', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'big step probability', true)

% wisk contact position (light, manip)
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.isLightOn; vars.condition];
dvMatrix = getDvMatrix(data, 'wiskContactPosition', conditions, varsToAvg, figConditionals)*-1000;
barPlotRick(dvMatrix, {conditions.levelNames}, {'obstacle distance', 'to nose at contact (mm)'}, true)

% tail height
plotInd = plotInd+1; subplot(rows, cols, plotInd);
conditions = [vars.isLightOn; vars.condition];
dvMatrix = getDvMatrix(data, 'tailHgt', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'tail height (m)', true)


% compare early and late success
if strcmp(manipulation, 'lesion')
    dvs = {'isTrialSuccess', 'trialVel'};
    for dv = dvs
        plotInd = plotInd+1; subplot(rows, cols, plotInd);
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

% big step prob by obs height (manip)

% big step prob by predicted distance to obs (manip)



%% big step kinematics

% !!! how to show this for different conditions?


%% sessions over time

% success, vel, body angle, baseline step height, 

%% speed vs. position / time plots

% speed vs position, control vs manip


% speed vs time centered around whisker contract


%% scatters

% obs hgt vs paw hgt (manip), plots for ipsi, contra, fore, hind


%% heat maps

% predicted vs actual planting distance, one map per paw


% delta mod step length vs predicted distance to obs (contra, manip)


%% kinematics

% fore/hind, contra, manip (leading only, light off only)


% fore/hind, leading, manip (contra only, light off only)













