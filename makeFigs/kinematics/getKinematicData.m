function data = getKinematicData(sessions)

% settings
isObsPosStatic = false; % if true, assumes wisk contacts obstacle at median detected obsContactPosition PER SESSION
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
    
    
    % LOAD SESSION DATA
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'],...
            'obsPositions', 'obsTimes', 'obsPixPositions', 'frameTimeStamps', 'mToPixMapping', ...
            'obsOnTimes', 'obsOffTimes', 'nosePos', 'targetFs', 'wheelPositions', 'wheelTimes', 'targetFs');
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));
    mToPixMapping = median(mToPixMapping,1);
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\tracking\locationsBotCorrected.mat'], 'locations')
    locations = locations.locationsCorrected;
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\tracking\stepSegmentation.mat'], ...
        'stanceBins', 'controlStepIdentities', 'modifiedStepIdentities')
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\tracking\isExcluded.mat'], 'excludedTrials')
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\tracking\velocityInfo.mat'], 'trialVels', 'minVel')
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\wiskContactData.mat'], 'contactTimes', 'contactPositions')
    vel = getVelocity(wheelPositions, speedTime, targetFs);
    locations(:,2,:) = locations(:,2,:) - nosePos(2); % subtract midline from all y values
    
    
    
    
    
    % COLLECT DATA FOR EACH TRIAL
    sessionVels = nan(1,length(obsOnTimes));
    
    for j = 1:length(obsOnTimes)
        
        % only analyze trials that have kinematic data (no nans in locations), are not excluded, and meet velocity criterion
        trialBins = frameTimeStamps>=obsOnTimes(j) & frameTimeStamps<=obsOffTimes(j) & ~isnan(obsPixPositions)';
        
        if ~any(isnan(reshape(locations(trialBins,:,:), 1, []))) ... % !!! why do i need to check that there are no NaN values...
                && ~ismember(j, excludedTrials) ...
                && trialVels(j)>=minVel
            
            
            
            
            % GET TIME, OBS POSITION, AND SPEED AT MOMENT OF CONTACT
            if isObsPosStatic

                % set contact position as median of contact positiosn for session
                contactPositions(j) = nanmedian(contactPositions);

                % interpolate to find time at which obsPosition==contactPosition
                % note: this is done by finding the ind before and after threshold crossing, then interpolating between these two sample points!
                indStart = find(obsPositions>=contactPositions(j) & obsTimes>obsOnTimes(j), 1, 'first');
                contactTimes(j) = interp1(obsPositions(indStart-1:indStart), obsTimes(indStart-1:indStart), contactPositions(j));
            end
            sessionVels(j) = interp1(wheelTimes, vel, contactTimes(j));



            
            % GET TRIAL DATA
            % note: i pull out trial specific data because 'find' function works much quicker on the smaller data slices // however, this feels inelegant // is there a better way of doing this?)

            % get time stamps relative to wisk contact
            trialTimeStamps = frameTimeStamps(trialBins)-contactTimes(j);
            [~, minInd] = min(abs(trialTimeStamps));
            trialTimeStampsInterp = trialTimeStamps - trialTimeStamps(minInd);

            % get trial data interpolated s.t. 0 is moment of wisk contact
            trialObsPixPositions = interp1(trialTimeStamps, obsPixPositions(trialBins), trialTimeStampsInterp);

            trialControlStepIds = nan(sum(trialBins), 4);
            trialModStepIds = nan(sum(trialBins), 4);
            trialLocations = nan(sum(trialBins), size(locations,2), size(locations,3));
            trialWheelVel = interp1(wheelTimes-contactTimes(j), vel, trialTimeStampsInterp);
            
            for k = 1:4
                trialControlStepIds(:,k) = interp1(trialTimeStamps, controlStepIdentities(trialBins,k), trialTimeStampsInterp, 'nearest');
                trialModStepIds(:,k) = interp1(trialTimeStamps, modifiedStepIdentities(trialBins,k), trialTimeStampsInterp, 'nearest');
                
                for m = 1:size(locations,2)
                    trialLocations(:,m,k) = interp1(trialTimeStamps, locations(trialBins,m,k), trialTimeStampsInterp, 'linear', 'extrap');
                end
            end

            % correct x locations (transform them s.t. obs is always at position 0 and positions move forward as though there were no wheel)
            trialLocations(:,1,:) = trialLocations(:,1,:) - trialObsPixPositions;

            % convert to meters
            trialLocations = trialLocations / abs(mToPixMapping(1));
            
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
                pawControlLocations = nan(stepNum, 2, swingMaxSmps);
                pawControlLocationsInterp = nan(stepNum, 2, interpSmps);
                
                for m = 1:stepNum
                    
                    % locations
                    startInd = find(trialControlStepIds(:,k)==m, 1, 'first');
                    stepIndsAll = startInd:min(startInd+swingMaxSmps-1, size(trialLocations,1)); % these inds continue past the end of swing
                    stepX = trialLocations(stepIndsAll,1,k);
                    stepY = trialLocations(stepIndsAll,2,k);
                    pawControlLocations(m,:,1:length(stepIndsAll)) = cat(1,stepX',stepY');
                    
                    % locations interp
                    try
                    stepBins = trialControlStepIds(:,k)==m;
                    xInterp = interp1(1:sum(stepBins), trialLocations(stepBins,1,k), linspace(1,sum(stepBins),interpSmps));
                    yInterp = interp1(1:sum(stepBins), trialLocations(stepBins,2,k), linspace(1,sum(stepBins),interpSmps));
                    pawControlLocationsInterp(m,:,:) = cat(1,xInterp,yInterp);
                    catch; keyboard; end
                end
                
                controlLocations{k} = pawControlLocations;
                controlLocationsInterp{k} = pawControlLocationsInterp;
                
                
                % modified
                modStepNum(k) = max(trialModStepIds(:,k));
                pawModifiedLocations = nan(modStepNum(k), 2, swingMaxSmps);
                pawModifiedLocationsInterp = nan(modStepNum(k), 2, interpSmps);
                
                for m = 1:modStepNum(k)
                    
                    % locations
                    startInd = find(trialModStepIds(:,k)==m, 1, 'first');
                    stepIndsAll = startInd:min(startInd+swingMaxSmps-1, size(trialLocations,1));
                    stepX = trialLocations(stepIndsAll,1,k);
                    stepY = trialLocations(stepIndsAll,2,k);
                    pawModifiedLocations(m,:,1:length(stepIndsAll)) = cat(1,stepX',stepY');
                    
                    % locations interp
                    stepBins = trialModStepIds(:,k)==m;
                    xInterp = interp1(1:sum(stepBins), trialLocations(stepBins,1,k), linspace(1,sum(stepBins),interpSmps));
                    yInterp = interp1(1:sum(stepBins), trialLocations(stepBins,2,k), linspace(1,sum(stepBins),interpSmps));
                    pawModifiedLocationsInterp(m,:,:) = cat(1,xInterp,yInterp);
                    
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









