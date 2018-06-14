function data = getKinematicData(sessions, obsPos)

% note: only specify obPos if you would like to trigger analysis at specific obs position relative to mouse, as opposed to relative to time obs first contacts whiskers...

% settings
botPawInds = 1:4;
topPawInds = 8:11;
minVel = 0;
velPrePost = [-.1 .1]; % compute trials velocity between these obstacle positions (relative to tip of mouse's nose)
speedTime = .02; % compute velocity over this interval
interpSmps = 100; % strides are stretched to have same number of samples // interpSmps sets the number of samples per interpolated stride
swingMaxSmps = 50; % when averaging swing locations without interpolating don't take more than swingMaxSmps for each swing

% initializations
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx']);
controlSteps = 2; % !!! don't change this - currently requires this number is two
data = struct();
dataInd = 1;


% collect data for all trials
for i = 1:length(sessions)
    
    % report progress
    fprintf('%s: collecting data\n', sessions{i});
    
    
    % LOAD SESSION DATA (damn that's a lot of stuff homes)
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'],...
            'obsPositions', 'obsTimes', 'obsPixPositions', 'frameTimeStamps', 'mToPixMapping', 'isLightOn', ...
            'obsOnTimes', 'obsOffTimes', 'nosePos', 'targetFs', 'wheelPositions', 'wheelTimes', 'targetFs', ...
            'wheelRadius', 'wheelCenter');
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));
    mToPixFactor = median(mToPixMapping(:,1));
    locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' sessions{i} '\trackedFeaturesRaw2.csv']); % get raw tracking data
    [locations, features, featurePairInds, isInterped] = fixTrackingDLC(locationsTable, frameTimeStamps);
    trialVels = getTrialVels(velPrePost, obsOnTimes, obsTimes, obsPositions);
    vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runTop.mp4']);
    stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, wheelTimes, wheelPoints, 250, mToPixFactor);
    if exist('obsPos', 'var')
        contactPositions = ones(size(obsOnTimes))*obsPos;
        contactTimes = nan(size(obsOnTimes));
        
        % get times when obs reaches obPos
        for j = 1:length(obsOnTimes)
            indStart = find(obsPositions>=contactPositions(j) & obsTimes>obsOnTimes(j), 1, 'first');
            contactTimes(j) = interp1(obsPositions(indStart-1:indStart), obsTimes(indStart-1:indStart), contactPositions(j));
        end
    else
        load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\wiskContactData.mat'], 'contactTimes', 'contactPositions')
    end
    [controlStepIdentities, modifiedStepIdentities] = ...
        getStepIdentities(stanceBins, locations(:,:,botPawInds), contactTimes, frameTimeStamps, obsOnTimes, obsOffTimes, obsPixPositions);
    vel = getVelocity(wheelPositions, speedTime, targetFs);
    
    % put together xyz for paws only
    locationsPaws = nan(size(locations,1), 3, 4);
    locationsPaws(:,1:2,:) = locations(:,:,botPawInds);
    locationsPaws(:,3,:) = locations(:,2,topPawInds);
    locationsPaws(:,2,:) = locationsPaws(:,2,:) - (nosePos(2)+vidTop.Height); % subtract midline from all y values
    locationsPaws(:,3,:) = (wheelCenter(2)-wheelRadius) - locationsPaws(:,3,:); % flip z and set s.t. top of wheel is zero
    clear vidTop
    
    
    
    
    
    
    % COLLECT DATA FOR EACH TRIAL
    sessionVels = nan(1,length(obsOnTimes)); % vels at moment of contact
    trials = find(trialVels>minVel);
    
    for j = trials
        
        % GET TRIAL DATA
        % note: i pull out trial specific data because 'find' function works much quicker on the smaller data slices // however, this feels inelegant // is there a better way of doing this?)
        trialBins = frameTimeStamps>=obsOnTimes(j) & frameTimeStamps<=obsOffTimes(j) & ~isnan(obsPixPositions)';

        % get vel at moment of contact
        sessionVels(j) = interp1(wheelTimes, vel, contactTimes(j));
        
        % get time stamps relative to wisk contact
        trialTimeStamps = frameTimeStamps(trialBins)-contactTimes(j);
        [~, minInd] = min(abs(trialTimeStamps));
        trialTimeStampsInterp = trialTimeStamps - trialTimeStamps(minInd);

        % get trial data interpolated s.t. 0 is moment of wisk contact
        trialObsPixPositions = interp1(trialTimeStamps, obsPixPositions(trialBins), trialTimeStampsInterp);

        trialControlStepIds = nan(sum(trialBins), 4);
        trialModStepIds = nan(sum(trialBins), 4);
        trialLocations = nan(sum(trialBins), size(locationsPaws,2), size(locationsPaws,3));
        trialWheelVel = interp1(wheelTimes-contactTimes(j), vel, trialTimeStampsInterp);

        for k = 1:4
            trialControlStepIds(:,k) = interp1(trialTimeStamps, controlStepIdentities(trialBins,k), trialTimeStampsInterp, 'nearest');
            trialModStepIds(:,k) = interp1(trialTimeStamps, modifiedStepIdentities(trialBins,k), trialTimeStampsInterp, 'nearest');

            for m = 1:size(locationsPaws,2)
                trialLocations(:,m,k) = interp1(trialTimeStamps, locationsPaws(trialBins,m,k), trialTimeStampsInterp, 'linear', 'extrap');
            end
        end
        
        if any(all(isnan(trialControlStepIds),1) | all(isnan(trialModStepIds),1))
            fprintf('  missing steps in trial %i\n', j)
        else

            % correct x locations (transform them s.t. obs is always at position 0 and positions move forward as though there were no wheel)
            trialLocations(:,1,:) = trialLocations(:,1,:) - trialObsPixPositions;

            % convert to meters
            trialLocations = trialLocations / abs(mToPixFactor);

            % determine whether left and right forepaws are in swing at obsPos moment
            isLeftSwing = ~isnan(trialModStepIds(trialTimeStampsInterp==0,2));
            isRightSwing = ~isnan(trialModStepIds(trialTimeStampsInterp==0,3));
            oneSwingOneStance = xor(isLeftSwing, isRightSwing);

            % flip y values if the left fore is the swinging foot (thus making it the right paw)
            if oneSwingOneStance && isLeftSwing
                trialLocations = trialLocations(:,:,[4 3 2 1]);
                trialControlStepIds = trialControlStepIds(:,[4 3 2 1]);
                trialModStepIds = trialModStepIds(:,[4 3 2 1]);
                trialLocations(:,2,:) = -trialLocations(:,2,:);
                isFlipped = true;
            else
                isFlipped = false;
            end


            % get stance distance from obs
            stanceDistance = trialLocations(trialTimeStampsInterp==0,1,2); % left fore paw (2) is always the stance foot at this point after flipping y values above
            swingStartDistance = trialLocations(find(trialModStepIds(:,3)==1,1,'first'),1,3);


            % get mod and control step(s) length, duration, wheel velocity
            controlSwingLengths = nan(controlSteps,4);
            modifiedSwingLengths = nan(1,4);
            controlSwingDurations = nan(controlSteps,4);
            modifiedSwingDurations = nan(1,4);
            controlWheelVels = nan(controlSteps,4);
            modifiedWheelVels = nan(1,4);

            for k = 1:4

                % modified steps
                stepBins = trialModStepIds(:,k)==1;
                stepXLocations = trialLocations(stepBins,1,k);
                modifiedSwingLengths(k) = stepXLocations(end) - stepXLocations(1);
                stepTimes = trialTimeStampsInterp(stepBins);
                modifiedSwingDurations(1,k) = stepTimes(end) - stepTimes(1);
                modifiedWheelVels(k) = trialWheelVel(find(stepBins,1,'first'));

                % control steps
                for m = 1:controlSteps

                    stepBins = trialControlStepIds(:,k)==m;
                    stepXLocations = trialLocations(stepBins,1,k);
                    controlSwingLengths(m,k) = stepXLocations(end) - stepXLocations(1);
                    stepTimes = trialTimeStampsInterp(stepBins);
                    controlSwingDurations(m,k) = stepTimes(end) - stepTimes(1);
                    controlWheelVels(m,k) = trialWheelVel(find(stepBins,1,'first'));
                end
            end





            % GET CONTROL AND MOD PAW LOCATIONS (INTERP AND NON-INTERP)
            controlLocations = cell(1,4);
            modLocations = cell(1,4);
            controlLocationsInterp = cell(1,4);
            modLocationsInterp = cell(1,4);
            modStepNum = nan(1,4);
            pawObsPosIndInterp = nan(1,4);
            pawObsPosInd = nan(1,4);

            for k = 1:4

                % control
                stepNum = max(trialControlStepIds(:,k));
                pawControlLocations = nan(stepNum, 3, swingMaxSmps);
                pawControlLocationsInterp = nan(stepNum, 3, interpSmps);

                for m = 1:stepNum

                    % locations
                    startInd = find(trialControlStepIds(:,k)==m, 1, 'first');
                    stepIndsAll = startInd:min(startInd+swingMaxSmps-1, size(trialLocations,1)); % these inds continue past the end of swing !!! ideally they would stop at the start of the next swing
                    stepX = trialLocations(stepIndsAll,1,k);
                    stepY = trialLocations(stepIndsAll,2,k);
                    stepZ = trialLocations(stepIndsAll,3,k);
                    pawControlLocations(m,:,1:length(stepIndsAll)) = cat(1,stepX',stepY',stepZ');

                    % locations interp
                    stepBins = trialControlStepIds(:,k)==m;
                    xInterp = interp1(1:sum(stepBins), trialLocations(stepBins,1,k), linspace(1,sum(stepBins),interpSmps));
                    yInterp = interp1(1:sum(stepBins), trialLocations(stepBins,2,k), linspace(1,sum(stepBins),interpSmps));
                    zInterp = interp1(1:sum(stepBins), trialLocations(stepBins,3,k), linspace(1,sum(stepBins),interpSmps));
                    pawControlLocationsInterp(m,:,:) = cat(1,xInterp,yInterp,zInterp);
                end

                controlLocations{k} = pawControlLocations;
                controlLocationsInterp{k} = pawControlLocationsInterp;


                % modified
                modStepNum(k) = max(trialModStepIds(:,k));
                pawModifiedLocations = nan(modStepNum(k), 3, swingMaxSmps);
                pawModifiedLocationsInterp = nan(modStepNum(k), 3, interpSmps);

                for m = 1:modStepNum(k)

                    % locations
                    startInd = find(trialModStepIds(:,k)==m, 1, 'first');
                    stepIndsAll = startInd:min(startInd+swingMaxSmps-1, size(trialLocations,1));
                    stepX = trialLocations(stepIndsAll,1,k);
                    stepY = trialLocations(stepIndsAll,2,k);
                    stepZ = trialLocations(stepIndsAll,3,k);
                    pawModifiedLocations(m,:,1:length(stepIndsAll)) = cat(1,stepX',stepY', stepZ');

                    % locations interp
                    stepBins = trialModStepIds(:,k)==m;
                    xInterp = interp1(1:sum(stepBins), trialLocations(stepBins,1,k), linspace(1,sum(stepBins),interpSmps));
                    yInterp = interp1(1:sum(stepBins), trialLocations(stepBins,2,k), linspace(1,sum(stepBins),interpSmps));
                    zInterp = interp1(1:sum(stepBins), trialLocations(stepBins,3,k), linspace(1,sum(stepBins),interpSmps));
                    pawModifiedLocationsInterp(m,:,:) = cat(1,xInterp,yInterp,zInterp);

                    % get ind of obs hit in interpolated coordinates
                    if m==1
                        stepObsPosInd = find(trialTimeStampsInterp==0) - find(stepBins,1,'first') + 1;
                        pawObsPosIndInterp(k) = interp1(linspace(1,sum(stepBins),interpSmps), ...
                            1:interpSmps, stepObsPosInd, 'nearest');
                        pawObsPosInd(k) = find(trialTimeStampsInterp==0) - find(stepBins,1,'first') + 1;
                    end
                end

                modLocations{k} = pawModifiedLocations;
                modLocationsInterp{k} = pawModifiedLocationsInterp;

            end



            % STORE RESULTS
            sessionInfoBin = find(strcmp(sessionInfo.session, sessions{i}),1,'first');
            data(dataInd).mouse = sessionInfo.mouse{sessionInfoBin};
            data(dataInd).session = sessions{i};
            data(dataInd).trial = j;
            data(dataInd).isLightOn = isLightOn(j);

            data(dataInd).vel = sessionVels(j);  % mouse vel at moment of wisk contact
            data(dataInd).obsPos = contactPositions(j);       % position of obs relative to nose at moment of wisk contact
            data(dataInd).obsPosInd = find(trialTimeStampsInterp==0); % ind at which obs contacts wisks for trial
            data(dataInd).pawObsPosInd = pawObsPosInd;% ind at which obs contacts wisks for locations for each paw
            data(dataInd).pawObsPosIndInterp = pawObsPosIndInterp; % ind at which obs contacts wisks for interp locations for each paw
            data(dataInd).timeStamps = trialTimeStamps;
            data(dataInd).locations = trialLocations;
            data(dataInd).controlLocations = controlLocations;
            data(dataInd).modifiedLocations = modLocations;
            data(dataInd).controlLocationsInterp = controlLocationsInterp;
            data(dataInd).modifiedLocationsInterp = modLocationsInterp;
            data(dataInd).trialControlStepIdentities = trialControlStepIds;
            data(dataInd).modifiedStepIdentities = trialModStepIds;
            data(dataInd).modStepNum = modStepNum;
            data(dataInd).oneSwingOneStance = oneSwingOneStance;
            data(dataInd).stanceDistance = stanceDistance;
            data(dataInd).swingStartDistance = swingStartDistance;
            data(dataInd).isFlipped = isFlipped;

            data(dataInd).controlSwingLengths = controlSwingLengths;
            data(dataInd).modifiedSwingLengths = modifiedSwingLengths;
            data(dataInd).controlSwingDurations = controlSwingDurations;
            data(dataInd).modifiedSwingDurations = modifiedSwingDurations;
            data(dataInd).controlWheelVels = controlWheelVels;
            data(dataInd).modifiedWheelVels = modifiedWheelVels;

            dataInd = dataInd + 1;
        end
    end
end



% make model to predict would-be mod swing length using wheel vel and previous swing lengths as predictors
mice = unique({data.mouse});
models = cell(1,length(mice));

for i = 1:length(mice)
    
    mouseBins = strcmp({data.mouse}, mice{i});
    
    % make predictive model
    
    % predictors: wheel speed, previous stride length
    prevLengths = cellfun(@(x) x(1,3), {data(mouseBins).controlSwingLengths});
    vel = cellfun(@(x) x(2,3), {data(mouseBins).controlWheelVels});

    % dependent variable: stride length
    lengths = cellfun(@(x) x(2,3), {data(mouseBins).controlSwingLengths});
    
    % make linear model
    models{i} = fitlm(cat(1,prevLengths,vel)', lengths, 'Linear', 'RobustOpts', 'on');
    
    % generate control length predictions (this is used to validate method)
    predictedLengths = num2cell(predict(models{i}, cat(1,prevLengths,vel)'));
    [data(mouseBins).predictedControlLengths] = predictedLengths{:};
    
    % generate mod length predictions
    prevLengths = cellfun(@(x) x(2,3), {data(mouseBins).controlSwingLengths});
    vel = cellfun(@(x) x(1,3), {data(mouseBins).modifiedWheelVels});
    predictedLengths = num2cell(predict(models{i}, cat(1,prevLengths,vel)'));
    [data(mouseBins).predictedLengths] = predictedLengths{:};
        
end


fprintf('--- done collecting data ---\n');









