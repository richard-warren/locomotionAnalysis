% test getFrameTimes2

while true

    % !!! breaks when misses frame right before reward...
    
    % settings
    frames = 250000;
    missedFrames = 100;
    rewards = 50;
    beginningDrop = 250;  % how many frames to drop at the beginning

    % simulate data
    ttlTimes = 100:(1/250):(frames/250+100);
    rewardInds = linspace(1, length(ttlTimes), rewards+2);
    rewardInds = round(rewardInds(2:end-1));
    deltas = zeros(size(ttlTimes));
    deltas(rewardInds) = .5;
    ttlTimes = ttlTimes + cumsum(deltas);

    frameTimesRaw = ttlTimes + 1000;  % simulate clock offsets
    frameCounts = 200:(200+length(frameTimesRaw));

    % drop frames
    bins = true(size(frameTimesRaw));
    missedInds = sort(randsample(length(frameTimesRaw)-beginningDrop+1, missedFrames)+beginningDrop+1);
    bins(missedInds) = false;
    frameTimesRaw = frameTimesRaw(bins);
    frameCounts = frameCounts(bins);

    % get ground truth frame times
    frameTimesTrue = ttlTimes';
    frameTimesTrue = frameTimesTrue(bins);  % remove missing frames

    % remove first second of frames
    frameTimesRaw = frameTimesRaw(beginningDrop+1:end);
    frameCounts = frameCounts(beginningDrop+1:end);
    frameTimesTrue = frameTimesTrue(beginningDrop+1:end);

%     frameTimesTrue(1:find(diff(frameTimesTrue)>.1,1,'first')) = nan; % replacing first trial with nans

    % run algorithm
%     frameTimes = getFrameTimes2(ttlTimes', frameTimesRaw', frameCounts', 'test');
%     frameTimes = getFrameTimes4(ttlTimes', frameTimesRaw', frameCounts', 'test');
    frameTimes = getFrameTimes(ttlTimes', frameTimesRaw', frameCounts', 'test');

    % show where they differ
    close all; figure('position', [2402.00 276.00 560.00 420.00]);
    plot(frameTimesTrue, 'lineWidth', 2); hold on; plot(frameTimes)
    frameTimesTrue(isnan(frameTimesTrue)) = 0; frameTimes(isnan(frameTimes)) = 0;  % in matlab nan~=nan, so turn nans to zeros here for sake of comparison
    differInds = find(frameTimesTrue~=frameTimes);
    scatter(differInds, frameTimesTrue(differInds))
    fprintf('differing frames: %i\n', length(differInds))
    
    if ~isequal(frameTimes, frameTimesTrue); break; end
end


%%