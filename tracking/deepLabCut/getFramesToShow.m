function frameInds = getFramesToShow(session)

% gets frame inds starting at reward preceding first obs trial to end, with periods surrounding reward editing out
% doesn't include frames starting at reward delivery and continuing until wheel has rotated postRewardDistance

% settings
postRewardDistance = .5; % meters

% initializations
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
    'frameTimeStamps', 'rewardTimes', 'wheelPositions', 'wheelTimes', 'obsOnTimes')
rewardTimes = rewardTimes(find(rewardTimes<obsOnTimes(1),1,'last'):end); % start with reward preceding first obstacle trial
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