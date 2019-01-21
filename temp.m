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



dataOut = getStructFields(data, {'session', 'mouse', 'isLightOn', 'pawNum', 'isLeading', 'isIpsi', 'height', 'success'});

% get height for all light on trials, leading, ipsi
dv = mean([dataOut([dataOut.isLightOn] & [dataOut.isLeading] & [dataOut.isIpsi]).height]);

% get success for all light on trials where LH is leading
dv = mean([dataOut([dataOut.pawNum]==1 & [dataOut.isLeading]).success]);


