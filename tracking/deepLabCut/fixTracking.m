function [locations, features, featurePairInds, isInterped, scores] = fixTracking(locationsTable, frameTimeStamps)


% settings
xDiffMax = 20;
scoreThresh = .99;
featurePairNames = {'paw1', 'paw2', 'paw3', 'paw4'}; % two features both containing the same string in this array are considered the same feature in the top and bottom views...
maxSpeed = 5000; % pixels per second (5000 = 20 pix per frame at 250 fps)
lookAheadFrames = 100;
minimizeCorrections = false; % set to true to 'see' the tracking errors, which is useful for adding mistracked frames to training sets

% load data and convert from table to matrices
frameNum = height(locationsTable);
features = locationsTable(:,2:3:width(locationsTable)).Properties.VariableNames;

locations = nan(frameNum, 2, length(features));
scores = nan(frameNum, length(features));

for i = 1:length(features)
    locations(:,:,i) = table2array(locationsTable(:, (i-1)*3+2 : (i-1)*3+3));
    scores(:,i) = table2array(locationsTable(:, (i-1)*3+4));
end

% find inds of feature pairs
pawInds = find(~cellfun(@isempty, strfind(features, 'paw')));

featurePairInds = nan(length(featurePairNames),2); % !!! this code is hacky
for i = 1:length(featurePairNames)
    pair = find(~cellfun(@isempty, strfind(features, featurePairNames{i})));
    if strfind(features{pair(1)}, 'top'); pair = fliplr(pair); end % ensures the bot features come first
    if length(pair)==2; featurePairInds(i,:) = pair; end
end


% check that vid has same number of frames as trackedFeaturesRaw has rows... fix if necessary
if length(frameTimeStamps)~=height(locationsTable)
    fprintf('WARNING: %i frames in video and %i frames in trackedFeaturesRaw\n', ...
        length(frameTimeStamps), height(locationsTable))
    if length(frameTimeStamps) > height(locationsTable)
        dif = length(frameTimeStamps)-height(locationsTable);
        locations(end+dif,:,:) = nan;
    else
        locations = locations(1:length(frameTimeStamps),:,:);
    end
end



% remove locations beneathe score thresh
lowScoreBins = permute(repmat(scores<scoreThresh,1,1,2), [1 3 2]);
locations(lowScoreBins) = nan;


if ~minimizeCorrections
    % remove inds that violate velocity constraint
    for i = 1:length(features)

        speeds = sqrt(sum(diff(squeeze(locations(:,:,i)), 1, 1).^2,2)) ./ diff(frameTimeStamps);
        checkVelInds = find(diff(isnan(locations(:,1,i)))==1 | speeds>maxSpeed) + 1;

        for j = checkVelInds'

            % find first frame where tracking doesn't violate max speed
            lastTrackedInd = j-1; while isnan(locations(lastTrackedInd,1,i)); lastTrackedInd = lastTrackedInd-1; end

            inds = lastTrackedInd+1 : min(lastTrackedInd+lookAheadFrames, size(locations,1));
            dp = sqrt(sum((locations(lastTrackedInd,:,i) - locations(inds,:,i)).^2,2));
            dt = frameTimeStamps(inds)-frameTimeStamps(lastTrackedInd);
            vels = dp ./ dt;

            nextTrackedInd = lastTrackedInd + find(vels<maxSpeed,1,'first');

            % replace values violating speed constraint with nans
            locations(lastTrackedInd+1:nextTrackedInd-1,:,i) = nan;
        end
    end



    % remove top locations where x is too far from bottom location
    for i = 1:size(featurePairInds)
        xDiffs = abs(diff(squeeze(locations(:,1,featurePairInds(i,:))), 1, 2));
        locations(xDiffs>xDiffMax, :, featurePairInds(i,2)) = nan;
    end
end



% fill in missing values
isInterped = isnan(squeeze(locations(:,1,:)));

for i = pawInds % 1:length(features)
    locations(:,1,i) = fillmissing(locations(:,1,i), 'pchip', 'endvalues', 'none');
    locations(:,2,i) = fillmissing(locations(:,2,i), 'pchip', 'endvalues', 'none');
end



% replace top x vals with bot x vals
% note: this is only for visualization purposes... in analyses i will use bot x vals only, but top x vals needed for making vids

% first find polyfit from bot x to top x values only using frames where neither view is interped)
% then replace top x values with predicted x values only for those frames that are interped in top view
for i = 1:size(featurePairInds,1)
    interpedInds = isInterped(:,featurePairInds(i,2)); % inds of frames where top is interped
    notIInerpedInds = ~any(isInterped(:,featurePairInds(i,:)),2); % inds where neither top or bottom are interped    
    fit = polyfit(locations(notIInerpedInds,1,featurePairInds(i,1)), locations(notIInerpedInds,1,featurePairInds(i,2)), 2);
    locations(interpedInds,1,featurePairInds(i,2)) = polyval(fit,locations(interpedInds,1,featurePairInds(i,1)));
end





