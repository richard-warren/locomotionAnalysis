function frameTimes = getFrameTimes3(ttlTimes, frameTimesRaw)

% given times of frames from camera metadata (frameTimesRaw) and ttls
% associated with these frames recorded in Spike (ttlTimes), assigns each
% frame in frameTimesRaw to the correct (hopefully) ttl in frameTimes //
% returns frameTimes, which are the times of each frame with respect to
% Spikes clock

% this is achieved by first stretching frameTimesRaw to compensate for
% differences in the rates of the clocks // then the lag that maximizes the
% xcorr between the two signals is determined // after these offset and
% drift corrections, the correct time for each frame is determined to be
% the ttlTime that is nearest the corrected frameTimesRaw, so long as it is
% within .001 seconds


% stretch frameTimesRaw to compensate for differences in clock speeds for spike and camera
medianDt = median(diff(ttlTimes)); % get frame rate

% get precise frame rate for camera and spike using average dt for all dt's
% wihtin .001 s of medianDt
ttlDiffs = diff(ttlTimes);
ttlDt = mean(ttlDiffs(abs(ttlDiffs-medianDt)<.001));
frameDiffs = diff(frameTimesRaw);
frameDt = mean(frameDiffs(abs(frameDiffs-medianDt)<.001));

frameTimesRaw = frameTimesRaw * (ttlDt/frameDt); % ratio between spike and camera Dts is used to stretch frameTimesRaw!


% find best lag by turning timestamps into vectors, then running cross-correlation

% turn into binary vectors
dt = .0001;
binEdges = min([ttlTimes; frameTimesRaw]):dt:max([ttlTimes; frameTimesRaw]);
ttlBins = histcounts(ttlTimes, binEdges);
frameBins = histcounts(frameTimesRaw, binEdges);

% convolve vectors with gaussian to make xcorr more sensible
sigma = median(diff(ttlTimes))*.5;
kernel = arrayfun(@(x) (1/(sigma*sqrt(2*pi))) * exp(-.5*(x/sigma)^2), -sigma*5:dt:sigma*5);


% determine lag that maximizes correlation
[r, lags] = xcorr(conv(ttlBins, kernel, 'same'), conv(frameBins, kernel, 'same'));
[~, maxLagInd] = max(r);
lag = lags(maxLagInd) * dt;


% plot camera and spike times after corrections
% figure;
% scatter(ttlTimes, zeros(1,length(ttlTimes)), 100); hold on
% scatter(frameTimesRaw+lag, zeros(1,length(frameTimesRaw)), 'filled');
% pimpFig; legend({'spike times', 'frame times'})


% assign each frame time to the nearest spike time

[ids, diffs] = knnsearch(ttlTimes, frameTimesRaw+lag);
frameTimes = ttlTimes(ids);
frameTimes(diffs>.001) = nan; % set to nan frames where there is no close ttl

if sum(~isnan(frameTimes)) ~= length(unique(frameTimes(~isnan(frameTimes))))
    disp('WARNING! same timestamp assigned to multiple frames!')
end










