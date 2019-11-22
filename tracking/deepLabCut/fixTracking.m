function [locations, features, featurePairInds, isInterped, scores] = ...
    fixTracking(locationsTable, frameTimeStamps, pixelsPerM)


% todo: check vel thing // document // add options for each stage of analysis?

% settings
xDiffMax = .005;  % (m)
scoreThresh = .99;
pairNames = {'paw1', 'paw2', 'paw3', 'paw4'};  % two features both containing the same string in this array are considered the same feature in the top and bottom views
maxSpeed = 2;  % (m/s) tracked feature cannot move faster than this across adjacent frames
lookAheadFrames = 100;



% load data and convert from table to matrices
frameNum = height(locationsTable);
features = locationsTable(:,2:3:width(locationsTable)).Properties.VariableNames;
locations = nan(frameNum, 2, length(features));
scores = nan(frameNum, length(features));

for i = 1:length(features)
    locations(:,:,i) = table2array(locationsTable(:, (i-1)*3+2 : (i-1)*3+3));  % get x and y colums for feature i
    scores(:,i) = table2array(locationsTable(:, (i-1)*3+4));
end

% find inds of feature pairs
featurePairInds = nan(length(pairNames), 2); % featurePairs X 2 matrix containing pair of feature inds for bottom and top views
for i = 1:length(pairNames)
    featureBins = contains(features, pairNames{i});
    featurePairInds(i,1) = find(featureBins & contains(features, '_bot'));
    featurePairInds(i,2) = find(featureBins & contains(features, '_top'));
end

% check that vid has same number of frames as trackedFeaturesRaw has rows... fix if necessary
if length(frameTimeStamps) ~= frameNum
    fprintf('WARNING: %i frames in video and %i frames in trackedFeaturesRaw!\n', ...
        length(frameTimeStamps), height(locationsTable))
    
    if length(frameTimeStamps) > frameNum
        dif = length(frameTimeStamps)-frameNum;
        locations(end+dif,:,:) = nan;
    else
        locations = locations(1:length(frameTimeStamps),:,:);
    end
end




% 1) remove locations beneathe score thresh
lowScoreBins = permute(repmat(scores<scoreThresh,1,1,2), [1 3 2]);
locations(lowScoreBins) = nan;




% 2) remove inds that violate velocity constraint
maxSpeedPixels = maxSpeed * pixelsPerM;
for i = 1:length(features)

    speeds = sqrt(sum(diff(squeeze(locations(:,:,i)), 1, 1).^2,2)) ./ diff(frameTimeStamps);
    checkVelInds = find(diff(isnan(locations(:,1,i)))==1 | speeds>maxSpeedPixels) + 1;

    for j = checkVelInds'

        % find first frame where tracking doesn't violate max speed
        lastTrackedInd = j-1; while isnan(locations(lastTrackedInd,1,i)); lastTrackedInd = lastTrackedInd-1; end

        inds = lastTrackedInd+1 : min(lastTrackedInd+lookAheadFrames, size(locations,1));
        dp = sqrt(sum((locations(lastTrackedInd,:,i) - locations(inds,:,i)).^2,2));
        dt = frameTimeStamps(inds)-frameTimeStamps(lastTrackedInd);
        vels = dp ./ dt;

        nextTrackedInd = lastTrackedInd + find(vels<maxSpeedPixels,1,'first');

        % replace values violating speed constraint with nans
        locations(lastTrackedInd+1:nextTrackedInd-1,:,i) = nan;
    end
end




% 3) remove top locations where x is too far from bottom location
xDiffMaxPixels = xDiffMax * pixelsPerM;
for i = 1:size(featurePairInds)
    xDiffs = abs(diff(squeeze(locations(:,1,featurePairInds(i,:))), 1, 2));
    locations(xDiffs>xDiffMaxPixels, :, featurePairInds(i,2)) = nan;
end




% 4) fill in missing values AND median filter?
pawInds = find(contains(features, 'paw'));
isInterped = isnan(squeeze(locations(:,1,:)));

locations(:,1,pawInds) = medfilt2(squeeze(locations(:,1,pawInds)), [3, 1]);
locations(:,2,pawInds) = medfilt2(squeeze(locations(:,2,pawInds)), [3, 1]);

locations(:,1,pawInds) = fillmissing(squeeze(locations(:,1,pawInds)), 'pchip', 'endvalues', 'none');
locations(:,2,pawInds) = fillmissing(squeeze(locations(:,2,pawInds)), 'pchip', 'endvalues', 'none');
% locations(:,1,pawInds) = fillgaps(squeeze(locations(:,1,pawInds)), 1000, 50);
% locations(:,2,pawInds) = fillgaps(squeeze(locations(:,2,pawInds)), 1000, 50);




% 5) replace top x vals with bot x vals
% note: this is only for visualization purposes... in analysis i will use bot x vals only, but top x vals needed for making vids

% first find polyfit from bot x to top x values only using frames where neither view is interped)
% then replace top x values with predicted x values only for those frames that are interped in top view
for i = 1:size(featurePairInds,1)
    interpedBins = isInterped(:,featurePairInds(i,2));  % bins of frames where top is interped
    fitBins = ~any(isInterped(:,featurePairInds(i,:)),2);  % bins where neither top or bottom are interped, used to fit model mapping bottom to top x values
    
    x = locations(fitBins, 1, featurePairInds(i,1));
    y = locations(fitBins, 1, featurePairInds(i,2));
    fit = polyfit(x, y, 2);
    
    locations(interpedBins,1,featurePairInds(i,2)) = polyval(fit, locations(interpedBins,1,featurePairInds(i,1)));
end





