% settings
sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003', ...
            '180125_001', '180125_002', '180125_003'};


% initializations
data = getKinematicData(sessions);
dataNew = data([data.oneSwingOneStance]);





%% scatter

% settings
mouse = 'run8';

% initializations
inds = strcmp({dataNew.mouse}, mouse);

% predictors: wheel speed, previous stride length
lengths1 = cellfun(@(x) x(1,3), {dataNew(inds).controlSwingLengths});
vels2 = cellfun(@(x) x(2,3), {dataNew(inds).controlWheelVels});

% dependent variable: stride length
lengths2 = cellfun(@(x) x(2,3), {dataNew(inds).controlSwingLengths});

% collect into matrix
regData = cat(1, lengths1, vels2, lengths2)';