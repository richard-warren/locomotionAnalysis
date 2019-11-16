function frameTimes = getFrameTimes4(ttlTimes, camTimes, camCounts, session)

% finds the fimes at which camera frames wre acquired with respect to Spike clock
% does this by identifying gap in TTLs at the very end of the recording
%
% input    ttlTimes:      timestamps for times of camera exposure (vidTtl) recorded in Spike, which has temporal gaps in between trials
%          camTimes:      timestamps from camera metadata (s)
%          camCounts:     frame counter from camera metadata

% output   frameTimeStamps: timestamps for every frame collected in video


% TEMP (use for debugging)
% session = '191007_003';
% load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mat'), 'exposure')
% camMetadata = dlmread(fullfile(getenv('OBSDATADIR'), 'sessions', '191007_003', 'run.csv')); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
% ttlTimes = exposure.times;
% camTimes = timeStampDecoderFLIR(camMetadata(:,3));
% camCounts = camMetadata(:,2);


% settings
ttlGap = [.495 .505];  % the last frame is collected after a gap // look for a gap that is within these limits // this shojld match the setting in the Arduino script servoControl


% find second to last frame in camera metadata and video TTLs in spike
ttlGaps = find(diff(ttlTimes)>ttlGap(1) & diff(ttlTimes)<ttlGap(2));
camGaps = find(diff(camTimes)>ttlGap(1) & diff(camTimes)<ttlGap(2));


% only one TTL gap should be found in between the last two frames for both the camera metadata and video TTLs in spike
if length(ttlGaps)==1 && ttlGaps==length(ttlTimes)-1 && ...
   length(camGaps)==1 && camGaps==length(camTimes)-1
    
    camCountsShifted = camCounts - camCounts(end) + length(ttlTimes);  % shift camCounts such that last camCount has same index as last video TTL in spike
    frameTimes = ttlTimes(camCountsShifted);
else
    fprintf('  %s: WARNING! Problem finding TTL gaps. Should have found only one gap in between final two frames...\n', session);
end

% display number of missed frames
fprintf('  %s: %i missed frames, of which %i were at the very beginning\n', ...
    session, length(ttlTimes)-length(camTimes), camCountsShifted(1)-1);





