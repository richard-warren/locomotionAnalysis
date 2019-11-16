function frameTimes = getFrameTimes(ttlTimes, camTimes, camCounts, session)

% this is an attempt to make a getFrameTimes that works for both old and
% new sessions
%
% input    ttlTimes:      timestamps for times of camera exposure (vidTtl) recorded in Spike, which has temporal gaps in between trials
%          camTimes:      timestamps from camera metadata (s)
%          camCounts:     frame counter from camera metadata

% output   frameTimeStamps: timestamps for every frame collected in video


% settings
ttlGap = [.495 .505];  % the last frame is collected after a gap // look for a gap that is within these limits // this shojld match the setting in the Arduino script servoControl

% find ttlGaps in spike ttls camera metadata
ttlGaps = find(diff(ttlTimes)>ttlGap(1) & diff(ttlTimes)<ttlGap(2));
camGaps = find(diff(camTimes)>ttlGap(1) & diff(camTimes)<ttlGap(2));

% using the first ttlGap as a reference, shift camCounts such that the
% camCount corresponding to the first gap has same index as spike ttl
% corresponding to the first gap
camCountsShifted = camCounts - camCounts(camGaps(1)) + ttlGaps(1); 
frameTimes = ttlTimes(camCountsShifted);

% display number of missed frames
missedFrames = length(ttlTimes)-length(camTimes);
fprintf('  %s: %i missed frames, %i at the very beginning, %i elsewhere\n', ...
    session, missedFrames, camCountsShifted(1)-1, missedFrames - (camCountsShifted(1)-1));





