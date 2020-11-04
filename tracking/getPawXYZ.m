function [pawXYZ, pawXYZ_pixels] = getPawXYZ(session)

% gets XYZ locations for paws by loading raw tracking data, applying
% fixTracking, and then saving everything in a nice table // pawXYZ is
% converted to meters, and adjusted such that z is relative to top of
% wheel, and y is relative to bodymidline // pawXYZ_pixels is in the
% original pixel coordinates

% todo: could rewrite this to accept as arguments all the things that are
% loaded up to increase the speed


% initializations
pawNames = {'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH'};

% load data and run fixTracking
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'wheelCenter', 'wheelRadius', 'pixelsPerM', 'nosePos')
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
scoreThresh = getScoreThresh(session, 'trackedFeaturesRaw_metadata.mat');  % scoreThresh depends on whether deeplabcut (old version) or deepposekit was used
[locations, locationsFeatures] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM, 'scoreThresh', scoreThresh);


[pawXYZ, pawXYZ_pixels] = deal(table(frameTimeStamps, 'VariableNames', {'t'}));

% stitch together xyz
botPawInds = contains(locationsFeatures, 'paw') & contains(locationsFeatures, '_bot');
topPawInds = contains(locationsFeatures, 'paw') & contains(locationsFeatures, '_top');
locationsPaws_pixels = nan(size(locations,1), 3, 4);  % frameNum X xyz X pawNum
locationsPaws_pixels(:,1:2,:) = locations(:,:,botPawInds);
locationsPaws_pixels(:,3,:) = locations(:,2,topPawInds);

locationsPaws = locationsPaws_pixels;
locationsPaws(:,2,:) = locationsPaws_pixels(:,2,:) - nosePos(2); % subtract midline from all y values
locationsPaws(:,3,:) = (wheelCenter(2)-wheelRadius) - locationsPaws(:,3,:); % flip z and set s.t. top of wheel is zero
locationsPaws = locationsPaws / pixelsPerM;  % convert to meters

% put results in table
for i = 1:length(pawNames)
    pawXYZ.(pawNames{i}) = locationsPaws(:,:,i);
    pawXYZ_pixels.(pawNames{i}) = locationsPaws_pixels(:,:,i);
end

