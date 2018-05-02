% test speed of different analyses

session = '180122_001';
velPrePost = [.1 .1]; % compute trials velocity between these obstacle positions (relative to tip of mouse's nose)


load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'])
load([getenv('OBSDATADIR') 'sessions\' session '\wiskContactData.mat'], 'contactTimes')
obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));
mToPixFactor = median(mToPixMapping(:,1)); % get mapping from meters to pixels

%% fix tracking
[locations, features, featurePairInds, isInterped] = fixTrackingDLC(session);

%% vels
trialVels = getTrialVels(velPrePost, obsOnTimes, obsTimes, obsPositions);

%% stance bins

vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
wheelPoints = getWheelPoints(vidTop);
stanceBins = getStanceBins(frameTimeStamps, locations(:,:,8:11), wheelPositions, wheelTimes, wheelPoints, 250, mToPixFactor);

%% step segmentation

[controlStepIdentities, modifiedStepIdentities] = ...
    getStepIdentities(stanceBins, locations, contactTimes, frameTimeStamps, obsOnTimes, obsOffTimes, obsPixPositions);