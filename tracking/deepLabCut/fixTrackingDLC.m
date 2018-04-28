function [locations, features, featurePairInds] = fixTrackingDLC(session)

% !!! currently assumes features for the top and second in the list of features, after features for the bottom... i should fix this yo

% settings
xDiffMax = 10;
scoreThresh = .95;
featurePairNames = {'paw1', 'paw2', 'paw3', 'paw4', 'tailBase', 'tailMid'}; % two features both containing the same string in this array are considered the same feature in the top and bottom views...
maxSpeed = 5000; % pixels per second (5000 = 20 pix per frame at 250 fps)


% load data and convert from table to matrices
locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' session '\trackingRaw.csv']); % get tracking data
% load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'frameTimeStamps'); % get tracking data
frameTimeStamps = cumsum(ones(1,2500))*.004; % temp hack !!!
frameNum = height(locationsTable)-2;
features = unique(table2cell(locationsTable(1,2:end)), 'stable');
locations = nan(frameNum, 2, length(features));
scores = nan(frameNum, length(features));

for i = 1:length(features)
    locations(:,:,i) = cellfun(@str2num, locationsTable{3:end, 2+(i-1)*3 : 3+(i-1)*3});
    scores(:,i) = cellfun(@str2num, locationsTable{3:end, 4+(i-1)*3});
end

% find inds of feature pairs
featurePairInds = nan(0,2);
for i = 1:length(featurePairNames)
    pair = find(~cellfun(@isempty, strfind(features, featurePairNames{i})));
    if length(pair)==2; featurePairInds(end+1,:) = pair; end
end



% remove locations beneathe score thresh
lowScoreBins = permute(repmat(scores<scoreThresh,1,1,2), [1 3 2]);
locations(lowScoreBins) = nan;



% remove inds that jump too far

for i = 1:length(features)
    
    speeds = sqrt(sum(diff(squeeze(locations(:,:,i)), 1, 1).^2,2)) ./ diff(frameTimeStamps)';
    checkVelInds = find(diff(isnan(locations(:,1,i)))==1 | speeds>maxSpeed) + 1;
%     checkVelInds = find(speeds>maxSpeed) + 1;
    
    for j = checkVelInds'
        
        % find first frame where tracking doesn't violate max speed
        lastDetectedInd = j-1; while isnan(locations(lastDetectedInd,1,i)); lastDetectedInd = lastDetectedInd-1; end
        nextTrackedInd = j;
        nextTrackedSpeed = nan;
        while (nextTrackedSpeed>maxSpeed || isnan(nextTrackedSpeed)) && nextTrackedInd<size(locations,1)
            nextTrackedInd = nextTrackedInd + 1;
            nextTrackedSpeed = sqrt(sum((locations(nextTrackedInd,:,i) - locations(lastDetectedInd,:,i)).^2)) / ...
                (frameTimeStamps(nextTrackedInd)-frameTimeStamps(lastDetectedInd));
        end

        % replace values violating speed constraint with nans
        locations(j:nextTrackedInd-1,:,i) = nan;
    end
end



% remove top locations where x is too far from bottom location
for i = 1:size(featurePairInds)
    xDiffs = abs(diff(squeeze(locations(:,1,featurePairInds(i,:))), 1, 2));
    locations(xDiffs>xDiffMax,1,featurePairInds(i,2)) = nan;
end


% fill in missing values
for i = 1:length(features)
    locations(:,1,i) = fillmissing(locations(:,1,i), 'spline', 'endvalues', 'none');
    locations(:,2,i) = fillmissing(locations(:,2,i), 'spline', 'endvalues', 'none');
end























