function obsPositionsFixed = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, noseX)



% temp
% session = '180123_001';
% load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
%     'obsPositions', 'obsTimes', 'obsPixPositions', 'frameTimeStamps', 'obsOnTimes', 'obsOffTimes', 'nosePos');
% noseX = nosePos(1);


obsPositionsFixed = nan(size(obsPositions));

for i = 1:length(obsOnTimes)
    
    % get trial pixPositions and pixTimes
    trialFrameBins = (frameTimeStamps>=obsOnTimes(i)) & (frameTimeStamps<=obsOffTimes(i)) & ~isnan(obsPixPositions)';
    pixPositions = obsPixPositions(trialFrameBins);
    pixTimes = frameTimeStamps(trialFrameBins);
    
    % get obsPos at moment obs reaches nose
    if ~isempty(pixPositions)
        try
        noseTime = interp1(pixPositions, pixTimes, noseX);
        obsAtNosePos = interp1(obsTimes, obsPositions, noseTime);
        catch; keyboard; end
    
        % get trial obsPos and subtract obsAtNosePos
        trialObsPosBins = (obsTimes>=obsOnTimes(i)) & (obsTimes<=obsOffTimes(i));
        obsPositionsFixed(trialObsPosBins) = obsPositions(trialObsPosBins) - obsAtNosePos;
    end
end


% close all; figure; plot(obsPositionsFixed); pimpFig





% obsPositionsFixed = obsPositions;
% 
% for i = 1:length(obsOnTimes)
%     
%     if i<length(obsOnTimes)
%         endTime = obsOnTimes(i+1);
%     else
%         endTime = max(obsTimes);
%     end
%     
%     trialInds = (obsTimes>obsOnTimes(i)) & (obsTimes<=endTime);
%     minPos = min(obsPositions(trialInds));
%     obsPositionsFixed(trialInds) = obsPositionsFixed(trialInds) - minPos;
%     
% end