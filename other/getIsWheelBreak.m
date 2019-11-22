function isWheelBreak = getIsWheelBreak(session)

% for a session, returns binary vector indicating whether there has been a
% wheel break on a given trial

% initializations
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'obsOnTimes', 'obsOffTimes')
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mat'), 'breaks')

isWheelBreak = true(size(obsOnTimes));
for i = 1:length(obsOnTimes)
    isWheelBreak(i) = any(breaks.times>obsOnTimes(i) & breaks.times<obsOffTimes(i));
end