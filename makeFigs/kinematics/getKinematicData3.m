function data = getKinematicData3(sessions, obsPos)

% note: only specify obPos if you would like to trigger analysis at specific obs position relative to mouse, as opposed to relative to time obs first contacts whiskers...

% settings
minVel = 0; % only include trials where mouse is running at least this fast
velPrePost = [-.1 .1]; % compute trials velocity between these obstacle positions (relative to tip of mouse's nose)
speedTime = .02; % compute velocity over this interval
interpSmps = 100; % strides are stretched to have same number of samples // interpSmps sets the number of samples per interpolated stride
swingMaxSmps = 50; % when averaging swing locations without interpolating don't take more than swingMaxSmps for each swing
noObsSteps = 5;
controlSteps = 2; % needs to be at least 2

% initializations
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx']);
data = struct();
dataInd = 1;


% collect data for all trials
for i = 1:length(sessions)
    
    % report progress
    fprintf('%s: collecting data\n', sessions{i});
    
    
    % LOAD SESSION DATA (damn that's a lot of stuff)
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'],...
            'obsPositions', 'obsTimes', 'obsPixPositions', 'obsPixPositionsContinuous', 'frameTimeStamps', 'mToPixMapping', 'isLightOn', ...
            'obsOnTimes', 'obsOffTimes', 'nosePos', 'targetFs', 'wheelPositions', 'wheelTimes', 'targetFs', ...
            'wheelRadius', 'wheelCenter', 'obsHeightsVid');
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));
    mToPixFactor = abs(mToPixMapping(1));
    locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' sessions{i} '\trackedFeaturesRaw.csv']); % get raw tracking data
    [locations, features] = fixTrackingDLC(locationsTable, frameTimeStamps);
    botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
    topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
    trialVels = getTrialVels(velPrePost, obsOnTimes, obsTimes, obsPositions);
    stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, wheelTimes, wheelCenter, wheelRadius, 250, mToPixMapping(1));
    
    if exist('obsPos', 'var')
        contactPositions = ones(size(obsOnTimes))*obsPos;
        contactTimes = nan(size(obsOnTimes));
        
        % get times when obs reaches obPos
        for j = 1:length(obsOnTimes)
            indStart = find(obsPositions>=contactPositions(j) & obsTimes>obsOnTimes(j), 1, 'first');
            if ~isempty(indStart)
                contactTimes(j) = interp1(obsPositions(indStart-1:indStart), obsTimes(indStart-1:indStart), contactPositions(j));
            end
        end
    else
        load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\wiskContactData.mat'], 'contactTimes', 'contactPositions')
    end
    [controlStepIdentities, modifiedStepIdentities, noObsStepIdentities] = ...
        getStepIdentities(stanceBins, locations(:,:,botPawInds), contactTimes, frameTimeStamps, ...
        obsOnTimes, obsOffTimes, obsPixPositions, obsPixPositionsContinuous, controlSteps, noObsSteps);
    vel = getVelocity(wheelPositions, speedTime, targetFs);
    
    % put together xyz for paws only
    locationsPaws = nan(size(locations,1), 3, 4);
    locationsPaws(:,1:2,:) = locations(:,:,botPawInds);
    locationsPaws(:,3,:) = locations(:,2,topPawInds);
    locationsPaws(:,2,:) = locationsPaws(:,2,:) - nosePos(2); % subtract midline from all y values
    locationsPaws(:,3,:) = (wheelCenter(2)-wheelRadius) - locationsPaws(:,3,:); % flip z and set s.t. top of wheel is zero
    locationsPaws(:,1,:) = locationsPaws(:,1,:) - obsPixPositionsContinuous'; % unravel x coordinates s.t. 0 is position of obs and x coordinates move forward over time ('unheadfixing')
    locationsPaws = locationsPaws / mToPixFactor; % convert to meters
    
    
    
    
    
    
    % COLLECT DATA FOR EACH TRIAL
    sessionVels = nan(1,length(obsOnTimes)); % vels at moment of contact
    trials = find(trialVels>minVel);
    
    for j = trials
        try
            % GET TRIAL DATA
            % note: i pull out trial specific data because 'find' function works much quicker on the smaller data slices // however, this feels inelegant // is there a better way of doing this?)

            % find trial bins, beginning with obsOffSteps preceding
            % obstacle and ending with obs turning off
            if j==1; lastTrialEndTime=0; else; lastTrialEndTime = obsOffTimes(j-1); end

            startInd = find(frameTimeStamps>lastTrialEndTime & any(noObsStepIdentities==1,2), 1, 'first'); % first bin where there is any obsOff step preceding obstacle
            trialBins = (1:length(frameTimeStamps))'>=startInd & frameTimeStamps<=obsOffTimes(j);

            % get vel at moment of contact
            sessionVels(j) = interp1(wheelTimes, vel, contactTimes(j));

            % get time stamps relative to wisk contact
            trialTimeStamps = frameTimeStamps(trialBins)-contactTimes(j);
            [~, minInd] = min(abs(trialTimeStamps));
            trialTimeStampsInterp = trialTimeStamps - trialTimeStamps(minInd);

            % get trial data interpolated s.t. 0 is moment of wisk contact
            trialControlStepIds = nan(sum(trialBins), 4);
            trialModStepIds = nan(sum(trialBins), 4);
            trialNoObsStepIds = nan(sum(trialBins), 4);
            trialLocations = nan(sum(trialBins), size(locationsPaws,2), size(locationsPaws,3));
            trialWheelVel = interp1(wheelTimes-contactTimes(j), vel, trialTimeStampsInterp);

            for k = 1:4
                trialControlStepIds(:,k) = interp1(trialTimeStamps, controlStepIdentities(trialBins,k), trialTimeStampsInterp, 'nearest');
                trialModStepIds(:,k) = interp1(trialTimeStamps, modifiedStepIdentities(trialBins,k), trialTimeStampsInterp, 'nearest');
                trialNoObsStepIds(:,k) = interp1(trialTimeStamps, noObsStepIdentities(trialBins,k), trialTimeStampsInterp, 'nearest');

                for m = 1:size(locationsPaws,2)
                    trialLocations(:,m,k) = interp1(trialTimeStamps, locationsPaws(trialBins,m,k), trialTimeStampsInterp, 'linear', 'extrap');
                end
            end

            % check whether any control or modified or obsOff steps are missing
            missingControlStep = false;
            for k = 1:controlSteps
                if ~all(any(trialControlStepIds==k,1)); missingControlStep = true; end
            end

            missingObsOffStep = false; % !!! perhaps i shouldnt require that all obsOffSteps are present for a trial to be included...
            for k = 1:noObsSteps
                if ~all(any(trialNoObsStepIds==k,1)); missingObsOffStep = true; end
            end

            missingModStep = any(all(isnan(trialModStepIds),1));
            somethingWentWrong = missingModStep || missingControlStep || missingObsOffStep;
        catch
            somethingWentWrong = true;
        end


        
        
        
        % analyze trial if all steps are accounted for
        if somethingWentWrong
            fprintf('  missing steps in trial %i\n', j)
        else

            % determine whether left and right forepaws are in swing at obsPos moment
            isLeftSwingAtContact = ~isnan(trialModStepIds(trialTimeStampsInterp==0,2));
            isRightSwingAtContact = ~isnan(trialModStepIds(trialTimeStampsInterp==0,3));

            % determine paw that gets over obs first
            isStepping = ~isnan(trialModStepIds);
            lastModStepInds = table2array(rowfun(@(x)(find(x,1,'last')), table(isStepping')));
            [~, firstPawOver] = min(lastModStepInds .* [nan 1 1 nan]'); % mask out hind paws, inds 1 and 4
            
            % determine which is paw first modified step (forepaws only)
            %
            % if both one in swing and one in stance at contact, first mod step is one
            % in swing // if both in stance at moment of contact, first mod
            % step is first to enter swing // if both in swing, first mod
            % step is first paw to land on the other side
            if xor(isLeftSwingAtContact, isRightSwingAtContact)
                if isLeftSwingAtContact; firstModPaw = 2; else; firstModPaw = 3; end
            elseif isLeftSwingAtContact && isRightSwingAtContact
                firstModPaw = firstPawOver;
            elseif ~isLeftSwingAtContact && ~isRightSwingAtContact
                firstModStepInds = table2array(rowfun(@(x)(find(x,1,'first')), table(isStepping')));
                [~, firstModPaw] = min(firstModStepInds .* [nan 1 1 nan]'); % mask out hind paws, inds 1 and 4
            end


            % get stance distance from obs, but only for trials in which
            % both forepaws are not in swing at moment of contact
            if ~(isLeftSwingAtContact && isRightSwingAtContact)
                if firstModPaw==2; stancePaw=3; else; stancePaw=2; end
                stanceDistance = trialLocations(trialTimeStampsInterp==0,1,stancePaw);
            else
                stanceDistance = nan;
            end
            
            % get distance of firstModPaw to obs at the beginning of first mod step
            swingStartDistance = trialLocations(find(trialModStepIds(:,3)==1,1,'first'),1,firstModPaw);


            % get mod, control, and noObs step(s) length, duration, wheel velocity
            controlSwingLengths = nan(controlSteps,4);
            noObsSwingLengths = nan(noObsSteps,4);
            modifiedSwingLengths = nan(1,4);
            controlSwingDurations = nan(controlSteps,4);
            noObsSwingDurations = nan(noObsSteps,4);
            modifiedSwingDurations = nan(1,4);
            controlWheelVels = nan(controlSteps,4);
            noObsWheelVels = nan(noObsSteps,4);
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
                
                % no obs steps
                for m = 1:controlSteps
                    stepBins = trialNoObsStepIds(:,k)==m;
                    stepXLocations = trialLocations(stepBins,1,k);
                    noObsSwingLengths(m,k) = stepXLocations(end) - stepXLocations(1);
                    stepTimes = trialTimeStampsInterp(stepBins);
                    noObsSwingDurations(m,k) = stepTimes(end) - stepTimes(1);
                    noObsWheelVels(m,k) = trialWheelVel(find(stepBins,1,'first'));
                end
            end





            % GET CONTROL, NO OBS, AND MOD PAW LOCATIONS (INTERP AND NON-INTERP)
            controlLocations = cell(1,4);
            modLocations = cell(1,4);
            noObsLocations = cell(1,4);
            controlLocationsInterp = cell(1,4);
            modLocationsInterp = cell(1,4);
            noObsLocationsInterp = cell(1,4);
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


                % no obs
                stepNum = max(trialNoObsStepIds(:,k));
                pawNoObsLocations = nan(stepNum, 3, swingMaxSmps);
                pawNoObsLocationsInterp = nan(stepNum, 3, interpSmps);
                for m = 1:stepNum

                    % locations
                    startInd = find(trialNoObsStepIds(:,k)==m, 1, 'first');
                    stepIndsAll = startInd:min(startInd+swingMaxSmps-1, size(trialLocations,1)); % these inds continue past the end of swing !!! ideally they would stop at the start of the next swing
                    stepX = trialLocations(stepIndsAll,1,k);
                    stepY = trialLocations(stepIndsAll,2,k);
                    stepZ = trialLocations(stepIndsAll,3,k);
                    pawNoObsLocations(m,:,1:length(stepIndsAll)) = cat(1,stepX',stepY',stepZ');

                    % locations interp
                    stepBins = trialNoObsStepIds(:,k)==m;
                    xInterp = interp1(1:sum(stepBins), trialLocations(stepBins,1,k), linspace(1,sum(stepBins),interpSmps));
                    yInterp = interp1(1:sum(stepBins), trialLocations(stepBins,2,k), linspace(1,sum(stepBins),interpSmps));
                    zInterp = interp1(1:sum(stepBins), trialLocations(stepBins,3,k), linspace(1,sum(stepBins),interpSmps));
                    pawNoObsLocationsInterp(m,:,:) = cat(1,xInterp,yInterp,zInterp);
                end
                noObsLocations{k} = pawNoObsLocations;
                noObsLocationsInterp{k} = pawNoObsLocationsInterp;
                
                
                
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
            
            % trial metadata
            sessionInfoBin = find(strcmp(sessionInfo.session, sessions{i}),1,'first');
            data(dataInd).mouse = sessionInfo.mouse{sessionInfoBin};
            data(dataInd).session = sessions{i};
            data(dataInd).trial = j;
            data(dataInd).isLightOn = isLightOn(j);
            data(dataInd).obsHeightsVid = obsHeightsVid(j);

            % bunch of thangs
            data(dataInd).vel = sessionVels(j);  % mouse vel at moment of wisk contact
            data(dataInd).obsPos = contactPositions(j);       % position of obs relative to nose at moment of wisk contact
            data(dataInd).obsPosInd = find(trialTimeStampsInterp==0); % ind at which obs contacts wisks for trial
            data(dataInd).pawObsPosInd = pawObsPosInd;% ind at which obs contacts wisks for locations for each paw
            data(dataInd).pawObsPosIndInterp = pawObsPosIndInterp; % ind at which obs contacts wisks for interp locations for each paw
            data(dataInd).timeStamps = trialTimeStamps;
            data(dataInd).locations = trialLocations;
            data(dataInd).controlLocations = controlLocations;
            data(dataInd).modifiedLocations = modLocations;
            data(dataInd).noObsLocations = noObsLocations;
            data(dataInd).controlLocationsInterp = controlLocationsInterp;
            data(dataInd).modifiedLocationsInterp = modLocationsInterp;
            data(dataInd).noObsLocationsInterp = noObsLocationsInterp;
            data(dataInd).trialControlStepIdentities = trialControlStepIds;
            data(dataInd).modifiedStepIdentities = trialModStepIds;
            data(dataInd).modStepNum = modStepNum;
            data(dataInd).stanceDistance = stanceDistance;
            data(dataInd).swingStartDistance = swingStartDistance;
            
            % info about which paws did what
            data(dataInd).isLeftSwingAtContact = isLeftSwingAtContact;
            data(dataInd).isLeftRightAtContact = isRightSwingAtContact;
            data(dataInd).firstPawOver = firstPawOver;
            data(dataInd).firstModPaw = firstModPaw;

            data(dataInd).controlSwingLengths = controlSwingLengths;
            data(dataInd).modifiedSwingLengths = modifiedSwingLengths;
            data(dataInd).noObsSwingLengths = noObsSwingLengths;
            data(dataInd).controlSwingDurations = controlSwingDurations;
            data(dataInd).modifiedSwingDurations = modifiedSwingDurations;
            data(dataInd).noObsSwingDurations = noObsSwingDurations;
            data(dataInd).controlWheelVels = controlWheelVels;
            data(dataInd).modifiedWheelVels = modifiedWheelVels;
            data(dataInd).noObsWheelVels = noObsWheelVels;

            dataInd = dataInd + 1;
        end
    end
end



% make model to predict would-be mod swing length of first modified paw using wheel vel and previous swing lengths as predictors
mice = unique({data.mouse});
models = cell(1,length(mice));

for i = 1:length(mice)
    
    mouseBins = strcmp({data.mouse}, mice{i});
    mouseInds = find(mouseBins);
    
    % make predictive model
    % use second to last control length and last control wheel vel (at
    % first moment) to predict last control length
    % !!! is there a way to do this wiithout looping?
    prevLengths = nan(1, length(mouseInds));
    vel = nan(1, length(mouseInds));
    lengths = nan(1, length(mouseInds));
    for j = 1:length(mouseInds)
        ind = mouseInds(j);
        prevLengths(j) = data(ind).controlSwingLengths(end-1, data(ind).firstModPaw);
        vel(j) = data(ind).controlWheelVels(end, data(ind).firstModPaw);
        lengths(j) = data(ind).controlSwingLengths(end, data(ind).firstModPaw);
    end

    
    % make linear model
    validInds = ~isnan(prevLengths) & ~isnan(vel);
    models{i} = fitlm(cat(1,prevLengths(validInds),vel(validInds))', lengths, 'Linear', 'RobustOpts', 'on');
    
    % generate control length predictions (this is used to validate method)
    predictedLengths = num2cell(predict(models{i}, cat(1,prevLengths,vel)'));
    [data(mouseBins).predictedControlLengths] = predictedLengths{:};
    
    % generate mod length predictions
    prevLengths = nan(1, length(mouseInds));
    vel = nan(1, length(mouseInds));
    for j = 1:length(mouseInds)
        ind = mouseInds(j);
        prevLengths(j) = data(ind).controlSwingLengths(end, data(ind).firstModPaw);
        vel(j) = data(ind).modifiedWheelVels(1, data(ind).firstModPaw);
    end
    predictedLengths = num2cell(predict(models{i}, cat(1,prevLengths,vel)'));
    [data(mouseBins).predictedLengths] = predictedLengths{:};
        
end


fprintf('--- done collecting data ---\n');






