function frameTimes = getFrameTimes(ttlTimes, frameTimesRaw, frameCounts)

% finds the fimes at which camera frames wre acquired with respect to Spike clock
% does this by finding breaks in camera exposures, which correspond to reward times, and matching camera metadata to exposure TLLs within every reward epoch
%
% input    ttlTimes:        timestamps for times of camera exposure (vidTtl) recorded in Spike, which has temporal gaps in between trials
%          frameTimesRaw:   timestamps from camera metadata (s)
%          frameCounts:     frame counter from camera metadata

% output   frameTimeStamps: timestamps for every frame collected in video


% initializations
ttlTrialStartInds = [1; find(diff(ttlTimes)>.1)+1; length(ttlTimes)];
frameTimeRawDiffs = abs(diff(frameTimesRaw));
frameCountDiffs = [abs(diff(frameCounts)); 1];
missedFrameInds = find(frameCountDiffs>1);
frameTrialStarts = [1; find(frameTimeRawDiffs>.250) + 1; length(frameTimesRaw)];
frameTimes = ttlTimes;


% ensure the correct number of trials were detected
if length(frameTrialStarts) ~= length(ttlTrialStartInds); disp('  same number of rewards not detected! problem!!!'); end


for i = 1:(length(ttlTrialStartInds)-1)

    % count vidTtls in trial and frames in trial (the latter determined from camera metadata)
    ttlsInTrial =   ttlTrialStartInds(i+1) - ttlTrialStartInds(i);
    framesInTrial = frameTrialStarts(i+1) - frameTrialStarts(i);

    % find indices of the missing frames in the trial
    trialMissedFrameInds = missedFrameInds(missedFrameInds>=frameTrialStarts(i) & missedFrameInds<frameTrialStarts(i+1));

    % count the number of detected missed frames
    missedInTrial = sum(frameCountDiffs(trialMissedFrameInds) - 1);

    % if all frames aren't accounted for, set all times to zero (these values will subsequently be changed to nan)
    if (ttlsInTrial - framesInTrial - missedInTrial) ~= 0
        fprintf('  cannot resolve timeStamps for trial: %i\n', i);
        frameTimes(ttlTrialStartInds(i):ttlTrialStartInds(i+1)-1) = 0;
        frameTimes = [frameTimes(1:ttlTrialStartInds(i)+framesInTrial-1); frameTimes(ttlTrialStartInds(i+1):end)]; % trim trial so there are only as many zeros as there are frames // there MUST be a more elegant way to do this

    % otherwise set all missing frame timestamps to nan (these values will subsequently be removed)
    else
        indexShift = 1;
        for j = 1:length(trialMissedFrameInds)
            numMissed = frameCountDiffs(trialMissedFrameInds(j))-1;
            startInd = ttlTrialStartInds(i) + (trialMissedFrameInds(j) - frameTrialStarts(i)) + indexShift;
            frameTimes(startInd : startInd+numMissed-1) = nan;
            indexShift = indexShift + numMissed;
        end
    end
end

% remove nan values
frameTimes = frameTimes(~isnan(frameTimes));

% set 0 values to nan to signify that timestamps are unavailable for frames
frameTimes(frameTimes==0) = nan;

% display number of missed frames
fprintf('  total missed frames: %i\n', length(missedFrameInds));











