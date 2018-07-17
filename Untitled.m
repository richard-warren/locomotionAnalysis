session = '180714_000';

% settings
pythonPath = '\Users\rick\Anaconda3\envs\deepLabCut\python.exe';
confidenceThresh = .99;
proximityThresh = 15;

% run neural network classifier
% system([pythonPath ' pawContact\analyzeVideo.py ' getenv('OBSDATADIR') 'sessions ' session])
pawAnalyzed = readtable([getenv('OBSDATADIR') 'sessions\' session '\pawAnalyzed.csv']);
isAnyPawTouching = false(1,length(frameTimeStamps));
isAnyPawTouching(pawAnalyzed.framenum) = pawAnalyzed.forelimb>confidenceThresh | pawAnalyzed.hindlimb>confidenceThresh;

% get xz positions for paws
locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' session '\trackedFeaturesRaw.csv']); % get raw tracking data
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
    'frameTimeStamps', 'obsPixPositions', 'obsOnTimes', 'obsOffTimes')
[locations, features, ~, isInterped, scores] = fixTrackingDLC(locationsTable, frameTimeStamps);
xzInds = ismember(features, {});
pawXZ = nan(size(locations,1), 2, 4);
for i = 1:4
    pawXBin = contains(features, ['paw' num2str(i)]) & contains(features, '_bot');
    pawZBins = contains(features, ['paw' num2str(i)]) & contains(features, '_top');
    pawXZ(:,1,i) = locations(:,1,pawXBin);
    pawXZ(:,2,i) = locations(:,2,pawZBins);
end

%% get obs height in pixels for each trial
% (replace missing values with median of tracked values for each trial)
obsBin = ismember(features, 'obs_top');
obsHeights = nan(size(locations,1),1);
for i = 1:length(obsOnTimes)
    trialBins = frameTimeStamps>obsOnTimes(i) & frameTimeStamps<obsOffTimes(i);
    medianHgt = nanmedian(locations(trialBins,2,obsBin));
    obsHeights(trialBins) = locations(trialBins,2,obsBin);
    obsHeights(trialBins & isnan(locations(:,2,obsBin))) = medianHgt;
end
obsHeights = medfilt1(obsHeights,5);

%% get xz distance of paws to obs at all times

dx = squeeze(pawXZ(:,1,:)) - repmat(obsPixPositions',1,4);
dz = squeeze(pawXZ(:,2,:)) - repmat(obsHeights,1,4);
pawDistances = sqrt(dx.^2 + dz.^2);

%% determine which paws are touching obs
isTouching = repmat(isAnyPawTouching',1,4) & pawDistances<proximityThresh;






