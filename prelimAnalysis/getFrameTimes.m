function frameTimes = getFrameTimes(ttlTimes, camTimes, camCounts, session)

% find the times of cameras frames with respect to spike clock. ttls
% trigger frames acquisition and are recorded in spike. given the times of
% these ttls, as well as the times of each frames wrt the cameras clock, as
% well as the frame counts, find times of each frame wrt to camera. does
% this by finding gaps in ttls, which serve as a reference point in both
% time systems, and aligning to this reference point
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

% find the interval between ttlGaps that is the closest in duration for camera and ttls
% this ensures we are aligning to the same ttl gap
if length(ttlGaps)==1 && length(camGaps)==1  % if same number of ttl gaps, just use the first one
    ttlGap = 1;
    camGap = 1;
else
    ttlGapDurations = diff(ttlTimes(ttlGaps));
    camGapDurations = diff(camTimes(camGaps));
    
    if length(ttlGaps)==length(camGaps)  % if same number of gaps, use the most similar gap (to avoid using a gap where a frame may be lost at the beginning or the end)
        [~, minInd] = min(abs(ttlGapDurations-camGapDurations));
        ttlGap = minInd;
        camGap = minInd;
    else  % otherwise try to find the best match (note that the following approach will fail for evenly spaced ttl gaps or if two gaps happen to be of similar duration) // would be smarter to loop across all entries, only keeping those with closely match durations
        fprintf('%s: WARNING! different number of ttl gaps detected in camera and spike!\n', session);
        [inds, diffs] = knnsearch(ttlGapDurations, camGapDurations);
        [~, minInd] = min(diffs);
        ttlGap = inds(minInd);
        camGap = minInd;
    end
end

% using the first ttlGap as a reference, shift camCounts such that the
% camCount corresponding to the first gap has same index as spike ttl
% corresponding to the first gap
camCountsShifted = camCounts - camCounts(camGaps(camGap)) + ttlGaps(ttlGap);

% if camera was started before spike, account for missing ttls at beginning
if any(camCountsShifted<1)
    beginningNans = nan(sum(camCountsShifted<1), 1);
    camCountsShifted = camCountsShifted(camCountsShifted>0);
    fprintf('%s: WARNING! %i camera frame(s) detected before spike recording turned on!\n', session, length(beginningNans));
end

% if camera was stopped after spike, account for extra ttls at endeginning
if any(camCountsShifted>length(ttlTimes))
    endNans = nan(sum(camCountsShifted>length(ttlTimes)), 1);
    camCountsShifted = camCountsShifted(camCountsShifted<=length(ttlTimes));
    fprintf('%s: WARNING! %i camera frame(s) detected after last spike ttl!\n', session, length(endNans));
end

frameTimes = ttlTimes(camCountsShifted);
if exist('beginningNans', 'var'); frameTimes = [beginningNans; frameTimes]; end
if exist('endNans', 'var'); frameTimes = [frameTimes; endNans]; end

% display number of missed frames
missedFrames = length(ttlTimes)-length(camTimes);
if exist('beginningNans', 'var'); missedFrames = missedFrames + length(beginningNans); end
if exist('endNans', 'var'); missedFrames = missedFrames + length(endNans); end
fprintf('%s: %i missed frame(s), %i at the very beginning, %i elsewhere\n', ...
    session, missedFrames, camCountsShifted(1)-1, missedFrames - (camCountsShifted(1)-1));





