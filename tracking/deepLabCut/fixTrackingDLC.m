function [locations, features, featurePairInds, isInterped] = fixTrackingDLC(locationsTable, frameTimeStamps)

% !!! currently assumes features for the top and second in the list of features, after features for the bottom... i should fix this yo

% settings
xDiffMax = 10;
scoreThresh = .95;
featurePairNames = {'paw1', 'paw2', 'paw3', 'paw4'}; % two features both containing the same string in this array are considered the same feature in the top and bottom views...
maxSpeed = 5000; % pixels per second (5000 = 20 pix per frame at 250 fps)
lookAheadFrames = 100;

% load data and convert from table to matrices
frameNum = height(locationsTable);
features = locationsTable(:,1:3:width(locationsTable)).Properties.VariableNames;

locations = nan(frameNum, 2, length(features));
scores = nan(frameNum, length(features));

for i = 1:length(features)
    locations(:,:,i) = table2array(locationsTable(:, (i-1)*3+1 : (i-1)*3+2));
    scores(:,i) = table2array(locationsTable(:, (i-1)*3+3));
end

% find inds of feature pairs
featurePairInds = nan(0,2);
for i = 1:length(featurePairNames)
    pair = find(~cellfun(@isempty, strfind(features, featurePairNames{i})));
    if length(pair)==2; featurePairInds(end+1,:) = pair; end
end

% check that vid has same number of frames as trackedFeaturesRaw has rows... fix if necessary
if length(frameTimeStamps)~=height(locationsTable)
    temp = table2array(locationsTable);
    duplicateInds = [];
    for i=1:size(temp,1)-1
        if all(temp(i,:) == temp(i+1,:)); duplicateInds(end+1)=i; end
    end
    
    fprintf('WARNING: %i frames in video and %i frames in trackedFeaturesRaw\n', ...
        length(frameTimeStamps), height(locationsTable))
    if (height(locationsTable)-length(duplicateInds))==length(frameTimeStamps)
        fprintf('  removed %i duplicate rows in trackedFeaturesRaw\n', length(duplicateInds))
    else
        fprintf('  COULD NOT RESOLVE DIFFERENCE IN FRAME NUMBERS! WTF!\n')
    end
    
    nonDuplicateBins = ~ismember(1:size(locations,1), duplicateInds);
    locations = locations(nonDuplicateBins, :, :);
    scores = scores(nonDuplicateBins, :);
    locationsTable = locationsTable(nonDuplicateBins,:);
end




% remove locations beneathe score thresh
lowScoreBins = permute(repmat(scores<scoreThresh,1,1,2), [1 3 2]);
locations(lowScoreBins) = nan;



% median filter
% locations = medfilt1(locations, 3, [], 'omitnan');


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



% fill in missing values
isInterped = isnan(squeeze(locations(:,1,:)));

for i = 1:length(features)
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





