%% make fake nested struct data

mice = {'run1', 'run2', 'run3'};
sessions = {'180102_000', '180102_001', '180102_002'};
paws = {'LH', 'LF', 'RF', 'RH'};
trials = 100;

randLog = @(length) num2cell(logical(randi([0 1],1,length)));

data = struct('mouse', mice, 'sessions', ...
       struct('session', sessions, 'trials', ...
       struct('trialNum', num2cell(1:100), 'success', randLog(trials), ...
              'isLightOn', randLog(trials), 'paws', ...
       struct('pawNum', num2cell(1:4), 'height', num2cell(rand(1,4)), 'isLeading', randLog(4), ...
              'isIpsi', randLog(4), 'isFore', num2cell(logical([0 1 1 0]))))));
   
%% query things
dataOut = getNestedStructFields(data, {'session', 'mouse', 'isLightOn', 'pawNum', 'isLeading', 'isIpsi', 'height', 'success'});

% get height for all light on trials, leading, ipsi
dv = mean([dataOut([dataOut.isLightOn] & [dataOut.isLeading] & [dataOut.isIpsi]).height]);

% get success for all light on trials where LH is leading
dv = mean([dataOut([dataOut.pawNum]==1 & [dataOut.isLeading]).success]);

%% test out new stuff bro

% vars = {'mouse', 'session', 'paw', 'condition', 'side', 'brainRegion', 'isLightOn', 'isTrialSuccess', 'stepOverMaxHgt', 'obsHgt', ...
%         'isWheelBreak', 'velAtWiskContact', 'angleAtWiskContact', 'obsPosAtContact', 'trialVel', 'trialAngle', 'wiskContactPositions'};
vars = 'all';
data = getExperimentData(sessionInfo, vars);
dataFlat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'paw', 'trialAngle', 'trialAngleContra'});

%%

dv = 'penultStepLength';

paw = struct('name', 'paw', 'levels', 1:4, 'levelNames', {{'LH', 'LF', 'RF', 'RH'}});
pawType = struct('name', 'pawType', 'levels', {{'leadFore', 'lagFore', 'leadHind', 'lagHind'}}, 'levelNames', {{'leadFore', 'lagFore', 'leadHind', 'lagHind'}});
isContra = struct('name', 'isContra', 'levels', [0 1], 'levelNames', {{'ipsi', 'contra'}});
isWheelBreak = struct('name', 'isWheelBreak', 'levels', [0, 1], 'levelNames', {{'no break', 'break'}});
condition = struct('name', 'condition', 'levels', {{'saline', 'muscimol'}}, 'levelNames', {{'sal', 'mus'}});
isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'light', 'no light'}});



varsToAvg = {'mouse'};
vars = [condition; pawType];


dvMatrix = getDvMatrix(data, dv, vars, varsToAvg);
barPlotRick(dvMatrix, {vars.levelNames}, dv)

%%




