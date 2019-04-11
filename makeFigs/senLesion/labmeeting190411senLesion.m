%% -------------- 
% INITIALIZATIONS
% ---------------

% settings
varsToAvg = {'mouse'};
maxEarlySession = 3;

% load session metadata
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'senLesionNotes');
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
conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
conditionals.isPre = struct('name', 'condition', 'condition', @(x) strcmp(x, 'pre'));
conditionals.isEarly = struct('name', 'conditionNum', 'condition', @(x) x<=maxEarlySession);
conditionals.isLate = struct('name', 'conditionNum', 'condition', @(x) x>=5 & x<=8);

% figConditionals = [conditionals.isEarly];
figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals


%% compute kinData for all sessions (only need to do once)
overwriteKindata = false;
sessions = unique(sessionInfo(logical([sessionInfo.include]),:).session);
parfor i = 1:length(sessions)
    if ~exist(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'kinData.mat'), 'file') || overwriteKindata
        getKinematicData5(sessions{i});
    end
end

%% load experiment data
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'senLesion_data.mat'), 'data');
disp('senLesion data loaded!')

%% compute new data and append to loaded data
loadOldData = true;
if exist('data', 'var') && loadOldData; data = getExperimentData(sessionInfo, 'all', data); else; data = getExperimentData(sessionInfo, 'all'); end
tic; save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'senLesion_data.mat'), 'data', '-v7.3'); disp('data saved'); toc

%% compute experiment from scratch, in parallel
data = cell(1,length(mice));
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data = cat(2,data{:});
save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'senLesion_data.mat'), 'data');

%% ----------
% PLOT THINGS
%  ----------

%% bar plots

close all; 
figure('name', 'senLesion', 'color', 'white', 'menubar', 'none', 'position', [[2206 154 334 718]])

% success
subplot(3, 1, 1);
conditions = [vars.conditionSub];
dvMatrix = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'success rate', true)

% height shaping
subplot(3, 1, 2);
conditions = [vars.conditionSub];
dvMatrix = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, conditions, {'mouse'}, {'session'}, [conditionals.isLeading; conditionals.isFore]); % perform regression for each session, then average slopes across sessions for each mouse
barPlotRick(dvMatrix, {conditions.levelNames}, 'height shaping', true)

% big step probability (light, modPawContra, manip)
subplot(3, 1, 3);
conditions = [vars.conditionSub];
dvMatrix = getDvMatrix(data, 'isBigStep', conditions, varsToAvg, figConditionals);
barPlotRick(dvMatrix, {conditions.levelNames}, 'big step probability', true)

print -clipboard -dbitmap


%% speed vs. position / time plots

close all
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'conditionNum', 'isLightOn', 'condition', 'obsOnPositions', ...
    'velContinuousAtContact', 'velVsPosition', 'isWheelBreak', 'wiskContactPosition'});
flat = flat(~[flat.isWheelBreak]);
% flat = flat(~([flat.conditionNum]>maxEarlySession & contains({flat.condition}, 'post'))); % only use first 3 lesion sessions
flat = flat(ismember({flat.condition}, {'pre', 'postContra', 'noWisk'}));
yLims = [.2 .8];

% speed vs position, control vs manip
figure('name', 'senLesion', 'Color', 'white', 'MenuBar', 'none', 'Position', [[2000 50 600 299]], 'inverthardcopy', 'off')

plotDvPsth(flat, 'velVsPosition', [-.5 .2], 'condition')
line([0 0], yLims, 'color', [.5 .5 .5])
line(repmat(nanmean([flat.obsOnPositions]),1,2), yLims, 'color', [.5 .5 .5])
set(gca, 'YLim', yLims)
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')

print -clipboard -dbitmap;


%% height shaping scatters

close all;
figure('name', 'senLesion', 'color', 'white', 'menubar', 'none', 'position', [2585 50 415 377])

% settings
xLims = [2 10];
yLims = [0 20];
conditionNames = {'noWisk', 'postContra', 'pre'};

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'condition'});
% flat = flat(~([flat.conditionNum]>maxEarlySession & contains({flat.condition}, 'post'))); % only use first 3 lesion sessions
% flat = flat(strcmp({flat.mouse})); % add conditionals here
flat = flat(ismember({flat.condition}, conditionNames) & [flat.isFore] & [flat.isLeading]);
obsHgts = [flat.obsHgt]*1000;
pawHgts = [flat.preObsHgt]*1000;

% initializations
conditions = cellfun(@(x) find(ismember(conditionNames, x)), {flat.condition});
scatterPlotRick(cat(1,obsHgts,pawHgts), {'obstacle height', 'paw height'}, conditions, conditionNames)

print -clipboard -dbitmap;



%% kinematics

% settings
rowVar = vars.isFore;
colVar = vars.isLeading;
plotVar = vars.condition;
conditionNames = {'noWisk', 'postContra', 'pre'};

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'isTrialSuccess', 'trial', 'isLightOn', 'isWheelBreak', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isContra', 'isLeading', 'condition', 'stepOverKinInterp', 'isBigStep', 'preObsKin'});
% flat = flat(~([flat.conditionNum]>maxEarlySession & contains({flat.condition}, 'post'))); % only use first 3 lesion sessions
% flat = flat(~[flat.isBigStep]); % add conditionals here
flat = flat(ismember({flat.condition}, conditionNames) & [flat.isLeading] & [flat.isFore]);

%% initializations
close all;
figure('name', 'senLesion', 'color', 'white', 'menubar', 'none', 'position', [2585 111 660 316])
conditions = cellfun(@(x) find(ismember(conditionNames, x)), {flat.condition});

% kinData = permute(cat(3, flat.stepOverKinInterp), [3,1,2]);
kinData = permute(cat(3, flat.preObsKin), [3,1,2]);
kinData = kinData(:,[1,3],:); % keep only x and z dimensions
        
plotKinematics(kinData, [flat.obsHgt], conditions, conditionNames)

print -clipboard -dbitmap;

%% render experiment vids

trialsPerSession = 15;
overWriteVids = false;
flat = getNestedStructFields(data, {'mouse', 'session', 'condition', 'trial', 'velAtWiskContact'});

for i = 1:length(mice)
    
    mouseConditions = unique({flat(strcmp({flat.mouse}, mice{i})).condition});
    for j = 1:length(mouseConditions)
        conditionSessions = unique({flat(strcmp({flat.mouse}, mice{i}) & ...
                                         strcmp({flat.condition}, mouseConditions{j})).session});
        session = conditionSessions{1};
        trials = sort(randperm(max([flat(strcmp({flat.session}, session)).trial]), trialsPerSession));
        
        file = fullfile(getenv('OBSDATADIR'), 'editedVid', 'senLesion', [mice{i} '_' session '_' mouseConditions{j}]);
        if ~exist(file, 'file') || overWriteVids
            makeVidWisk(file, session, [-.05 .1], .15, trials);
        end
    end    
end
disp('all done!')




