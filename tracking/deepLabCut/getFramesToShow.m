function frameInds = getFramesToShow(session, onlyShowFramesNearObs)

% if onlyShowFramesNearObs is true, gets frame inds starting at reward preceding first obs trial to end, with periods surrounding reward editing out
% doesn't include frames starting at reward delivery and continuing until wheel has rotated postRewardDistance
% otherwise gets frames near obstacle

% settings
postRewardDistance = .5; % meters
obsPosRange = [-.05 .1]; % show frames where obs is obsPosRange meters in front of and behind the nose

% initializations
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
    'frameTimeStamps', 'rewardTimes', 'wheelPositions', 'wheelTimes', 'obsOnTimes', 'obsOffTimes', ...
    'obsPixPositions', 'obsTimes', 'obsPositions')
if ~exist('onlyShowFramesNearObs', 'var'); onlyShowFramesNearObs=false; end


% get frames when mouse is near obs
if onlyShowFramesNearObs
    
    % initialize obsPositions
    if exist('obsPixPositions', 'var')
        obsPositions = fixObsPositionsHacky(obsPositions);
    else
        obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));
    end
    
    frameInds = cell(1,length(obsOnTimes));
    for i = 1:length(obsOnTimes)
        % find trial indices
        startInd = find(obsTimes>obsOnTimes(i)  & obsPositions>=obsPosRange(1), 1, 'first');
        endInd   = find(obsTimes<obsOffTimes(i) & obsPositions<=obsPosRange(2), 1, 'last');
        frameInds{i} = find(frameTimeStamps>obsTimes(startInd) & frameTimeStamps<obsTimes(endInd));
    end
    frameInds = cat(1,frameInds{:});
    
    
    
% get all frames except when mouse is drinking reward    
else
    rewardTimes = rewardTimes(find(rewardTimes>obsOnTimes(1),1,'first'):end);
    toShowBins = false(size(frameTimeStamps));
    wheelPosInterp = interp1(wheelTimes, wheelPositions, frameTimeStamps);

    for i = 1:length(rewardTimes)-1
        if ~isnan(rewardTimes(i)) % this should never happen but i saw it happen once... wtf
            initialPos = wheelPosInterp(find(frameTimeStamps>rewardTimes(i),1,'first'));
            startInd = find(frameTimeStamps>rewardTimes(i) & wheelPosInterp>(initialPos+postRewardDistance), 1, 'first');
            endInd = find(frameTimeStamps<rewardTimes(i+1), 1, 'last');
            toShowBins(startInd:endInd) = true;
        end
    end

    frameInds = find(toShowBins);
end