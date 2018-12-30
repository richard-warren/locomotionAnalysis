% function getNeuralCodingData()

% prepares data struct used to determine relationship between behavior and
% other stimuli and neural activity // creates on row per session, and
% includes many kinematic variables and data about events in session
% (obstacle times, etc)

% TO DO: remove nan timestamps?

% temp
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'));
sessionInfo = sessionInfo(sessionInfo.numGoodUnits>0,:); % keep only sessions with usable units
i=1;
%%

% initializations
data = struct();

for i = 1:height(sessionInfo)
    
    % load session data
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessionInfo.session{i}, 'runAnalyzed.mat'));
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'neuralData.mat'), ...
        'spkTimes', 'unit_ids');
    load([getenv('OBSDATADIR') 'sessions\' session '\run.mat'], 'breaks');    
    locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', sessionInfo.session{i}, 'trackedFeaturesRaw.csv')); % get raw tracking data
    
    % session metadata
    data(i).session = sessionInfo.session{i};
    data(i).mouse = sessionInfo.mouse{i};
    
    % trial metadata
    data(i).is_light_on = isLightOn;
    data(i).obstacle_heights = obsHeightsVid;
    
    % neural data
    data(i).spike_times = spkTimes;
    data(i).unit_ids = unit_ids;
    
    
    % session events
    data(i).obstacle_on_times = obsOnTimes;
    data(i).whisker_contact_times = wiskContactTimes;
    data(i).wheel_break_times = breaks.times;
    data(i).reward_times = rewardTimes;
    
    
    % kinematics
    obsPixPositionsContinuous = getObsPixPositionsContinuous(...
        obsPosToWheelPosMappings, wheelTimes, wheelPositions, frameTimeStamps, ...
        obsPixPositions, obsPixPositionsUninterped, obsOnTimes, obsOffTimes);
    obsPositionsFixed = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));
    mToPixFactor = abs(mToPixMapping(1));
    [locations, features] = fixTrackingDLC(locationsTable, frameTimeStamps);
    botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
    topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
    stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, wheelTimes, wheelCenter, wheelRadius, 250, mToPixFactor);
    [controlStepIdentities, modifiedStepIdentities, noObsStepIdentities] = ...
        getStepIdentities(stanceBins, locations(:,:,botPawInds), wiskContactTimes, frameTimeStamps, ...
        obsOnTimes, obsOffTimes, obsPixPositions, obsPixPositionsContinuous, 2, 5); % last 2 arguments are number of control and no obstacle steps
    
    % put together xyz and convert to meters
    featuresTop = cellfun(@(x) x(1:end-4), features(contains(features, '_top')), 'UniformOutput', false);
    featuresBot = cellfun(@(x) x(1:end-4), features(contains(features, '_bot')), 'UniformOutput', false);
    featuresXyz = intersect(featuresTop, featuresBot);
    kinematics = nan(length(frameTimeStamps), 3, length(featuresXyz));d
    
    kinematics(:,1:2,:) = locations(:,:,contains(features, featuresXyz) & contains(features, '_bot')); % get xy values from bottom view
    kinematics(:,3,:) = locations(:,2,contains(features, featuresXyz) & contains(features, '_top')); % get z values from top view
    kinematics(:,2,:) = kinematics(:,2,:) - nosePos(2); % subtract midline from all y values
    kinematics(:,3,:) = (wheelCenter(2)-wheelRadius) - kinematics(:,3,:); % flip z and set s.t. top of wheel is zero
    kinematics = kinematics/ mToPixFactor; % convert to meters
    
    data(i).kinematics = kinematics;
    data(i).time_stamps = frameTimeStamps;
    data(i).is_paw_in_stance = stanceBins;
    data(i).body_angle = bodyAngles;
    
    
    % position and velocity
    data(i).wheel_velocity = getVelocity(wheelPositions, .02, targetFs); % second argument is seconds over which wheel speed is computed
    data(i).obstacle_position = obsPositionsFixed;
    
    
end
