function getNeuralCodingData(session)

% prepares data struct used to determine relationship between behavior and
% other stimuli and neural activity // creates on row per session, and
% includes many kinematic variables and data about events in session
% (obstacle times, etc)


% temp
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'));
sessionInfo = sessionInfo(strcmp(sessionInfo.session, session), :); % keep only sessions with usable units
% sessionInfo = sessionInfo(1:2,:); % temptemp



data = struct();

for i = 1:height(sessionInfo)
    
    fprintf('%s: preparing neural coding data...\n', sessionInfo.session{i});
    
    % load session data
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessionInfo.session{i}, 'runAnalyzed.mat'));
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessionInfo.session{i}, 'neuralData.mat'), ...
        'spkTimes', 'unit_ids');
    locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', sessionInfo.session{i}, 'trackedFeaturesRaw.csv')); % get raw tracking data
    validFrames = ~isnan(frameTimeStamps);
    
    % session metadata
    data(i).session = sessionInfo.session{i};
    data(i).mouse = sessionInfo.mouse{i};
    
    % trial metadata
    data(i).is_light_on = isLightOn;
    data(i).obstacle_heights = obsHeights';
    
    % neural data
    data(i).spike_times = spkTimes';
    data(i).unit_ids = unit_ids;
    
    
    % session events
    data(i).obstacle_on_times = obsOnTimes;
    data(i).obstacle_off_times = obsOffTimes;
    data(i).whisker_contact_times = wiskContactTimes';
    data(i).reward_times = rewardTimes;
    
    
    % kinematics
    [locations, features] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM);
%     botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
    topPawBins = contains(features, 'paw') & contains(features, '_top');
    stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawBins), wheelPositions, wheelTimes, wheelCenter, wheelRadius, 250, pixelsPerM);
%     [controlStepIdentities, modifiedStepIdentities, noObsStepIdentities] = ...
%         getStepIdentities(stanceBins, locations(:,:,botPawInds), wiskContactTimes, frameTimeStamps, ...
%         obsOnTimes, obsOffTimes, obsPixPositions, [], 2, 5); % last 2 arguments are number of control and no obstacle steps
    
    % put together xyz and convert to meters
    featuresTop = cellfun(@(x) x(1:end-4), features(contains(features, '_top')), 'UniformOutput', false);
    featuresBot = cellfun(@(x) x(1:end-4), features(contains(features, '_bot')), 'UniformOutput', false);
    featuresXyz = intersect(featuresTop, featuresBot);
    kinematics = nan(length(frameTimeStamps), 3, length(featuresXyz));
    
    for j = 1:length(featuresXyz)
        kinematics(:,1:2,j) = locations(:,:,contains(features, featuresXyz{j}) & contains(features, '_bot')); % get xy values from bottom view
        kinematics(:,3,j) = locations(:,2,contains(features, featuresXyz{j}) & contains(features, '_top')); % get z values from top view
    end
    kinematics(:,2,:) = kinematics(:,2,:) - nosePos(2); % subtract midline from all y values
    kinematics(:,1,:) = kinematics(:,1,:) - nosePos(1); % subtract nose position from all x values
    kinematics(:,3,:) = (wheelCenter(2)-wheelRadius) - kinematics(:,3,:); % flip z and set s.t. top of wheel is zero
    kinematics = kinematics/ pixelsPerM; % convert to meters
    
    data(i).kinematics = kinematics(validFrames,:,:);
    data(i).time_stamps = frameTimeStamps(validFrames);
    data(i).is_paw_in_stance = stanceBins(validFrames,:);
    data(i).body_angle = bodyAngles(validFrames);
    
    
    % osbtacle position, velocity, distance to reward
    temp = nan(size(obsPositionsFixed));
    for j = 1:length(obsOnTimes)
        trialBins = obsTimes>obsOnTimes(j) & obsTimes<obsOffTimes(j);
        temp(trialBins) = obsPositionsFixed(trialBins); % remove periods where obstacle is not one
    end
    obsPositionsFixed = interp1(obsTimes, temp, frameTimeStamps(validFrames)); % interpolate to frame times
    
    vel = getVelocity(wheelPositions, .02, targetFs);
    vel = interp1(wheelTimes, vel, frameTimeStamps(validFrames));
    
    distanceToReward = -wheelPositions;
    rewardEdges = [min(wheelTimes); rewardTimes; max(wheelTimes)];
    for j = 2:length(rewardEdges)
        epochBins = wheelTimes>=rewardEdges(j-1) & wheelTimes<rewardEdges(j);
        distanceToReward(epochBins) = distanceToReward(epochBins) - distanceToReward(find(epochBins,1,'last'));
    end
    distanceToReward(find(epochBins,1,'last'):end) = nan; % set to nan everything after last reward
    distanceToReward = interp1(wheelTimes, distanceToReward, frameTimeStamps(validFrames)); % interpolate to frame times
    
    data(i).wheel_velocity = vel; % second argument is seconds over which wheel speed is computed
    data(i).obstacle_position = obsPositionsFixed;
    data(i).distance_to_reward = distanceToReward;
    
    
    % paw contacts
    data(i).paw_touches = touchesPerPaw(validFrames,:);
    
end


% save
disp('saving results...')
body_parts = featuresXyz;
touch_class_names = touchClassNames(1:end-1); % remove no_touch class
fileName = fullfile(getenv('OBSDATADIR'), 'matlabData', 'ashokData', 'QZ_neuralData', [session, '.mat']);
save(fileName, 'data', 'body_parts', 'touch_class_names', '-v7.3', '-nocompression'); clear data;
disp('all done!')





