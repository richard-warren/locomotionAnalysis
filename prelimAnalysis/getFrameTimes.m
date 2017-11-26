function frameTimes = getFrameTimes(ttlTimes, frameTimesRaw, frameCounts, session)

% finds the fimes at which camera frames wre acquired with respect to Spike clock
% does this by finding breaks in camera exposures, which correspond to reward times, and matching camera metadata to exposure TLLs within every reward epoch
%
% input    ttlTimes:        timestamps for times of camera exposure (vidTtl) recorded in Spike, which has temporal gaps in between trials
%          frameTimesRaw:   timestamps from camera metadata (s)
%          frameCounts:     frame counter from camera metadata

% output   frameTimeStamps: timestamps for every frame collected in video
%
% NOTE:    this code breaks if there are many adjacent missed frames - it will think this is a reward break in the camera when it is not really...
%          the break in ttls ocurring at reward times should be .5 seconds. gaps greater than .4 will be considered reward times


% initializations

ttlTrialStartInds = [1; find(diff(ttlTimes)>.4)+1; length(ttlTimes)];                 % inds of ttls at which trials start, plus extra ones at end and beginning
frameTrialStartInds = [1; find(diff(frameTimesRaw)>.45) + 1; length(frameTimesRaw)];   % inds of frames at which trials start, plus extra ones at end and beginning

frameCountDiffs = diff(frameCounts);
missedFrameInds = find(frameCountDiffs>1);
frameTimes = ttlTimes; % frameTimes will be changed as missed frames are removed...


% ensure the correct number of trials were detected
if length(frameTrialStartInds) ~= length(ttlTrialStartInds)
    keyboard
    fprintf('  %s: same number of rewards not detected! WTF!!! STOPPING ANALYSIS!!!', session)
%     figure; plot(diff(frameTimesRaw)) % plotting this may reveal that many adjacent frames were lost, causing the code to think based on the gap in frames that a reward was reached, when in fact it was not 
    return
    
end


for i = 1:(length(ttlTrialStartInds)-1)

    % count vidTtls in trial and frames in trial (the latter determined from camera metadata)
    ttlsInTrial =   ttlTrialStartInds(i+1) - ttlTrialStartInds(i);
    framesInTrial = frameTrialStartInds(i+1) - frameTrialStartInds(i);

    % find indices of the missing frames in the trial
    trialMissedFrameInds = missedFrameInds(missedFrameInds>=frameTrialStartInds(i) & missedFrameInds<frameTrialStartInds(i+1));

    % count the number of detected missed frames
    missedInTrial = sum(frameCountDiffs(trialMissedFrameInds) - 1);

    % if all frames aren't accounted for, set all times to zero (these values will subsequently be changed to nan)
    if (ttlsInTrial - framesInTrial - missedInTrial) ~= 0
        
        fprintf('  %s: cannot resolve timeStamps for trial: %i\n', session, i);
        frameTimes(ttlTrialStartInds(i):ttlTrialStartInds(i+1)-1) = 0;
        
        % trim trial so there are only as many zeros as there are frames // is there a better way to do this???
        frameTimes = [frameTimes(1:ttlTrialStartInds(i)+framesInTrial-1); frameTimes(ttlTrialStartInds(i+1):end)];

    % otherwise set all missing frame timestamps to nan (these values will subsequently be removed)
    else
        indexShift = 1;
        
        for j = 1:length(trialMissedFrameInds)
            
            numMissed = frameCountDiffs(trialMissedFrameInds(j))-1;
            startInd = ttlTrialStartInds(i) + (trialMissedFrameInds(j) - frameTrialStartInds(i)) + indexShift;
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
fprintf('  %s: total missed frames: %i\n', session, length(missedFrameInds));











