function obsPositionsFixed = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, noseX)


obsPositionsFixed = nan(size(obsPositions));

for i = 1:length(obsOnTimes)
    
    % get trial pixPositions and pixTimes
    trialFrameBins = (frameTimeStamps>=obsOnTimes(i)) & (frameTimeStamps<=obsOffTimes(i)) & ~isnan(obsPixPositions)';
    pixPositions = obsPixPositions(trialFrameBins);
    pixTimes = frameTimeStamps(trialFrameBins);
    
    % get obsPos at moment obs reaches nose
    if ~isempty(pixPositions) && length(unique(pixPositions))>1 % !!! sometimes pixPositions were all zero for reasons unknown
        
        noseTime = interp1(pixPositions, pixTimes, noseX);
        obsAtNosePos = interp1(obsTimes, obsPositions, noseTime);
        
        % get trial obsPos and subtract obsAtNosePos
        trialObsPosBins = (obsTimes>=obsOnTimes(i)) & (obsTimes<=obsOffTimes(i));
        obsPositionsFixed(trialObsPosBins) = obsPositions(trialObsPosBins) - obsAtNosePos;
    end
end