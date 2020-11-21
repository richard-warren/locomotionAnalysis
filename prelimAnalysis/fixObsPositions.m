function obsPositionsFixed = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, noseX)

% obsPositionsFixed shifts obsPositions on a trial by trial basis such that
% 0 corresponds to the point at which the obstacle is beneath the nose of
% the mouse. in theory this could be done by simply adding a constant, but
% the rotary encoder measurements are subject (theoretically) to drift over
% time, so it is safer to do this trial by trial


obsPositionsFixed = nan(size(obsPositions));
epochTimes = [obsOnTimes; obsTimes(end)];
obsPixPositions = obsPixPositions(:)';  % enforce horizontal orientation
anyNotFixed = false;

for i = 1:length(obsOnTimes)
    
    % get trial pixPositions and pixTimes
    trialFrameBins = (frameTimeStamps>=obsOnTimes(i)) & (frameTimeStamps<=obsOffTimes(i)) & ~isnan(obsPixPositions)';
    pixPositions = obsPixPositions(trialFrameBins);
    pixTimes = frameTimeStamps(trialFrameBins);
    
    % get obsPos at moment obs reaches nose
    try 
        % remove duplicate positional values
        [pixPositions, uniqueInds] = unique(pixPositions, 'stable');
        pixTimes = pixTimes(uniqueInds);
        noseTime = interp1(pixPositions, pixTimes, noseX);  % time at which obstacle is underneath the nose
        obsAtNosePos = interp1(obsTimes, obsPositions, noseTime);
        
        % get trial obsPos and subtract obsAtNosePos
        trialObsPosBins = (obsTimes>=epochTimes(i)) & (obsTimes<epochTimes(i+1));
        obsPositionsFixed(trialObsPosBins) = obsPositions(trialObsPosBins) - obsAtNosePos;
    catch
        if ~anyNotFixed
            fprintf('WARNING! Could not fix obstacle positions for trial(s): %i', i)
            anyNotFixed = true;
        else
            fprintf(' %i', i);
        end
    end
end
if anyNotFixed; fprintf('\n'); end

% replace first obsPositions with first fixed obsPosition
firstRealInd = find(~isnan(obsPositionsFixed), 1, 'first');
obsPositionsFixed(1:firstRealInd-1) = obsPositionsFixed(firstRealInd);
