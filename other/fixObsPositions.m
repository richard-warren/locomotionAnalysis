function obsPositionsFixed = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, noseX)

% !!! what did this function do again??? i think it reformats obs positions
% so zero is where obs is under nose


obsPositionsFixed = nan(size(obsPositions));
epochTimes = [obsOnTimes; obsTimes(end)];

for i = 1:length(obsOnTimes)
    
    % get trial pixPositions and pixTimes
    trialFrameBins = (frameTimeStamps>=obsOnTimes(i)) & (frameTimeStamps<=obsOffTimes(i)) & ~isnan(obsPixPositions)';
    pixPositions = obsPixPositions(trialFrameBins);
    pixTimes = frameTimeStamps(trialFrameBins);
    
    % get obsPos at moment obs reaches nose
    if ~isempty(pixPositions) && length(unique(pixPositions))>1 % !!! sometimes pixPositions were all zero for reasons unknown
        
        % remove duplicate positional values
        [pixPositions, uniqueInds] = unique(pixPositions, 'stable');
        pixTimes = pixTimes(uniqueInds);
        noseTime = interp1(pixPositions, pixTimes, noseX);
        obsAtNosePos = interp1(obsTimes, obsPositions, noseTime);
        
        % get trial obsPos and subtract obsAtNosePos
%         trialObsPosBins = (obsTimes>=obsOnTimes(i)) & (obsTimes<=obsOffTimes(i));
%         obsPositionsFixed(trialObsPosBins) = obsPositions(trialObsPosBins) - obsAtNosePos;
        trialObsPosBins = (obsTimes>=epochTimes(i)) & (obsTimes<epochTimes(i+1));
        obsPositionsFixed(trialObsPosBins) = obsPositions(trialObsPosBins) - obsAtNosePos;
    end
end

% replace first obsPositions with first fixed obsPosition
firstRealInd = find(~isnan(obsPositionsFixed),1,'first');
obsPositionsFixed(1:firstRealInd-1) = obsPositionsFixed(firstRealInd);