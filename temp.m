session = '180808_001';

obsPos = -0.0087;
controlSteps = 2;
noObsSteps = 3;
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'],...
            'obsPositions', 'obsTimes', 'obsPixPositions', 'obsPixPositionsContinuous', 'frameTimeStamps', 'mToPixMapping', 'isLightOn', ...
            'obsOnTimes', 'obsOffTimes', 'nosePos', 'targetFs', 'wheelPositions', 'wheelTimes', 'targetFs', ...
            'wheelRadius', 'wheelCenter', 'obsHeightsVid', 'arePawsTouchingObs');
load([getenv('OBSDATADIR') 'sessions\' session '\run.mat'], 'breaks');
obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));
mToPixFactor = abs(mToPixMapping(1));
locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' session '\trackedFeaturesRaw.csv']); % get raw tracking data
[locations, features] = fixTrackingDLC(locationsTable, frameTimeStamps);
botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, wheelTimes, wheelCenter, wheelRadius, 250, mToPixMapping(1));

%%
contactPositions = ones(size(obsOnTimes))*obsPos;
        contactTimes = nan(size(obsOnTimes));
        
% get times when obs reaches obPos
for j = 1:length(obsOnTimes)
    indStart = find(obsPositions>=contactPositions(j) & obsTimes>obsOnTimes(j), 1, 'first');
    if ~isempty(indStart)
        contactTimes(j) = interp1(obsPositions(indStart-1:indStart), obsTimes(indStart-1:indStart), contactPositions(j));
    end
end

%%
[controlStepIdentities, modifiedStepIdentities, noObsStepIdentities] = ...
        getStepIdentities(stanceBins, locations(:,:,botPawInds), contactTimes, frameTimeStamps, ...
        obsOnTimes, obsOffTimes, obsPixPositions, obsPixPositionsContinuous, controlSteps, noObsSteps);


