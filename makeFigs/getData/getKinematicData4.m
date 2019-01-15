function [data, stanceBins, touchClassNames] = getKinematicData4(sessions, sessionInfo, previousData, obsPos)

% note: only specify obPos if you would like to trigger analysis at specific obs position relative to mouse, as opposed to relative to time obs first contacts whiskers...
% note: stanceBins is only for the LAST analyzed session // only useful when calling this for a single session
% note: touchClassNames is for the last sessions, AND ASSUMES THEY ARE
% THE SAME FOR ALL SESSIONS ANALYZED

% settings
speedTime = .02; % compute velocity over this interval
interpSmps = 100; % strides are stretched to have same number of samples // interpSmps sets the number of samples per interpolated stride
swingMaxSmps = 50; % when averaging swing locations without interpolating don't take more than swingMaxSmps for each swing
noObsSteps = 3;
controlSteps = 2; % needs to be at least 2
contactPosLimits = [-.02 .02]; % whisker cant only contact obs this far in front of and behind nose
timeOperations = false;
numWorkers = 4;

% remove previously analyzed sessions from list of sessions
if ~isempty(previousData)
    previousSessions = unique({previousData.session});
    sessions = sessions(~ismember(sessions, previousSessions));
end


% collect data for all sessions
data = cell(1,length(sessions));
stanceBins = cell(1,length(sessions));
getDataForSessionHandle = @getDataForSession;
if any(strcmp(sessionInfo.Properties.VariableNames, 'notes'))
    sessionInfo = sessionInfo(:, ~strcmp(sessionInfo.Properties.VariableNames, 'notes')); % remove 'notes' field
end
metaDataFields = sessionInfo.Properties.VariableNames;
metaDataFields = cat(2, metaDataFields, {'sessionNum', 'conditionNum'});


parfor (i = 1:length(sessions), numWorkers)
% for i = 1:length(sessions)
    try
        % get metadata for sessions
        sessionInfoBin = strcmp(sessionInfo.session, sessions{i});
        sessionMetaData = table2struct(sessionInfo(sessionInfoBin,:));
        
        % get sesionNum and conditionNum for mouse
        mouseBins = strcmp(sessionInfo.mouse, sessionMetaData.mouse);
        if contains('condition', sessionInfo.Properties.VariableNames)
            conditionBins = strcmp(sessionInfo.condition, sessionMetaData.condition);
        else
            conditionBins = false(1,height(sessionInfo));
        end
        sessionMetaData.sessionNum = find(strcmp(sessionInfo.session(mouseBins), sessionMetaData.session)); % session num is 1,2,3... for sequential sessions for a given mouse // used to plot performace across days
        sessionMetaData.conditionNum = find(strcmp(sessionInfo.session(mouseBins & conditionBins), sessionMetaData.session)); % first session for condition pre and condition post is 1, second session for condition pre and condition post is 2, etc // used to restrict how many days post lesion to include in analysis
        
        % get session data
        if sessionMetaData.include % !!! restrict sessions if conditionNum too high here???
            fprintf('%s: collecting data...\n', sessions{i});
            [data{i}, stanceBins{i}] = feval(getDataForSessionHandle, sessions{i}, sessionMetaData);
        else
            fprintf('%s: skipped\n', sessions{i})
        end
    catch
        fprintf('%s: unable to analyze session!\n', sessions{i});
    end
end

% concatenate data across sessions
try
    data = cat(2,data{:});
catch
    disp('error concatenating data across sessions! entering debug mode...')
end


% make model to predict would-be mod swing length of first modified paw using wheel vel and previous swing lengths as predictors
mice = unique({data.mouse});
models = cell(1,length(mice));

for i = 1:length(mice)
    
    mouseBins = strcmp({data.mouse}, mice{i});
    mouseInds = find(mouseBins);
    
    % make predictive model
    % use second to last control length and last control wheel vel (at first moment) to predict last control length
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

% concatenate new and previous data, first removing fields that don't align
if ~isempty(previousData)
    mismatches = setdiff(fieldnames(data), fieldnames(previousData));
    for i = 1:length(mismatches)
        fprintf('WARNING: remove field "%s" due to mismatch in previous and newly computed data...', mismatches{i})
        if ismember(mismatches{i}, fieldnames(data))
            data = rmfield(data, mismatches{i});
        else
            previousData = rmfield(previousData, mismatches{i});
        end
    end
    data = cat(2, previousData, data);
end






% ---------
% FUNCTIONS
% ---------

function [sessionData, stanceBins] = getDataForSession(session, sessionMetaData)
    
    sessionData = struct();
    dataInd = 1;

    % LOAD SESSION DATA (damn that's a lot of stuff)
    if timeOperations; tic; end
    if ismember('side', metaDataFields) % if a neural manipulation occured, this stores information about the side of the brain
        if strcmp(sessionMetaData.side, 'left'); contraLimb = 3;
        elseif strcmp(sessionMetaData.side, 'right'); contraLimb = 2;
        else; contraLimb = nan; end
    end
    load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'],...
            'obsPositions', 'obsTimes', 'obsPixPositions', 'obsPixPositionsUninterped', 'frameTimeStamps', 'mToPixMapping', 'isLightOn', ...
            'obsOnTimes', 'obsOffTimes', 'nosePos', 'targetFs', 'wheelPositions', 'wheelTimes', 'targetFs', 'obsPosToWheelPosMappings', ...
            'wheelRadius', 'wheelCenter', 'obsHeightsVid', 'touchesPerPaw', 'touchClassNames', 'wiskContactFrames', 'wiskContactTimes', 'wiskContactPositions', ...
            'frameTimeStampsWisk', 'bodyAngles', 'obsPositionsFixed');
    load([getenv('OBSDATADIR') 'sessions\' session '\run.mat'], 'breaks');
    obsPixPositionsContinuous = getObsPixPositionsContinuous(...
        obsPosToWheelPosMappings, wheelTimes, wheelPositions, frameTimeStamps, ...
        obsPixPositions, obsPixPositionsUninterped, obsOnTimes, obsOffTimes);
    mToPixFactor = abs(mToPixMapping(1));
    locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' session '\trackedFeaturesRaw.csv']); % get raw tracking data
    [locations, features] = fixTrackingDLC(locationsTable, frameTimeStamps);
    botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
    topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
    stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, wheelTimes, wheelCenter, wheelRadius, 250, mToPixMapping(1));
    if timeOperations; fprintf('loading session data: %i seconds\n', round(toc)); end
    
        
    
    % get positions and times when obs reaches obPos or touches the whisker
    if timeOperations; tic; end
    if exist('obsPos', 'var')
        contactPositions = ones(size(obsOnTimes))*obsPos;
        contactTimes = nan(size(obsOnTimes));
    else
        contactPositions = wiskContactPositions;
        contactTimes = wiskContactTimes;
    end
    
    for j = 1:length(obsOnTimes)
        
        % if obsPos is defined by user, find times at which obs reaches this position
        if exist('obsPos', 'var')
            indStart = find(obsPositionsFixed>=contactPositions(j) & obsTimes>obsOnTimes(j) & obsTimes<obsOffTimes(j), 1, 'first');
            if ~isempty(indStart)
                contactTimes(j) = interp1(obsPositionsFixed(indStart-1:indStart), obsTimes(indStart-1:indStart), contactPositions(j));
                contactTimes(j) = frameTimeStamps(knnsearch(frameTimeStamps, contactTimes(j))); % force time to be in frameTimeStamps
            end
        end
    end
    
    % remove invalid contact positions
    if ~exist('obsPos', 'var')
        contactPositions(contactPositions<contactPosLimits(1) | contactPositions>contactPosLimits(2)) = nan;
    end
    if timeOperations; fprintf('getting wisk contact pos/times: %i seconds\n', round(toc)); end
    
    
    % get step identities and wheel velocity
    if timeOperations; tic; end
    [controlStepIdentities, modifiedStepIdentities, noObsStepIdentities] = ...
        getStepIdentities(stanceBins, locations(:,:,botPawInds), contactTimes, frameTimeStamps, ...
        obsOnTimes, obsOffTimes, obsPixPositions, obsPixPositionsContinuous, controlSteps, noObsSteps);
    vel = getVelocity(wheelPositions, speedTime, targetFs);
    if timeOperations; fprintf('getting step identities: %i seconds\n', round(toc)); end
    
    
    % put together xyz for paws only
    if timeOperations; tic; end
    locationsPaws = nan(size(locations,1), 3, 4);
    locationsPaws(:,1:2,:) = locations(:,:,botPawInds);
    locationsPaws(:,3,:) = locations(:,2,topPawInds);
    locationsPaws(:,2,:) = locationsPaws(:,2,:) - nosePos(2); % subtract midline from all y values
    locationsPaws(:,3,:) = (wheelCenter(2)-wheelRadius) - locationsPaws(:,3,:); % flip z and set s.t. top of wheel is zero
    
    
    % unravel x coordinates s.t. 0 is position of obs and x coordinates move forward over time ('unheadfixing')
    % this is done trial by trial, which each trial ending with the last ind of the last mod step over the obs
    prevInd = 1;
    for j = 1:length(obsOnTimes)
        finalInd = find(frameTimeStamps>obsOffTimes(j) & all(isnan(modifiedStepIdentities),2), 1, 'first');
        if j==length(obsOnTimes); finalInd = length(frameTimeStamps); end
        locationsPaws(prevInd:finalInd,1,:) = locationsPaws(prevInd:finalInd,1,:) - obsPixPositionsContinuous(j, prevInd:finalInd)';
        prevInd = finalInd+1;
    end
    locationsPaws = locationsPaws / mToPixFactor; % convert to meters
    if timeOperations; fprintf('getting xyz coords: %i seconds\n', round(toc)); end
    
    
    
    
    
    % COLLECT DATA FOR EACH TRIAL
    if timeOperations; tic; end
    
    for j = 1:length(obsOnTimes)
        try
            % GET TRIAL DATA
            % note: i pull out trial specific data because 'find' function works much quicker on the smaller data slices // however, this feels inelegant // is there a better way of doing this?)

            % get trial bins
            if j==1; lastTrialEndTime=0; else; lastTrialEndTime = obsOffTimes(j-1); end
            startInd = find(frameTimeStamps>lastTrialEndTime & any(noObsStepIdentities==1,2), 1, 'first'); % start at first obs off step after previous trial
            endInd = find(all(isnan(modifiedStepIdentities),2) & frameTimeStamps>obsOffTimes(j), 1, 'first'); % end when all mod steps in trial are over
            trialInds = startInd:endInd;


            % get time stamps relative to wisk contact
            trialTimeStamps = frameTimeStamps(trialInds)-contactTimes(j);
            
            % get trial data
            trialControlStepIds = controlStepIdentities(trialInds,:);
            trialModStepIds = modifiedStepIdentities(trialInds,:);
            trialNoObsStepIds = noObsStepIdentities(trialInds,:);
            trialLocations = locationsPaws(trialInds,:,:);
            trialWheelVel = vel(trialInds);
            trialTouchesPerPaw = touchesPerPaw(trialInds, :);
            
            % check whether any control or modified or obsOff steps are missing
            missingControlStep = false;
            for k = 1:controlSteps
                if ~all(any(trialControlStepIds==k,1)); missingControlStep = true; end
            end

            missingObsOffStep = false; % !!! perhaps i should not require that all obsOffSteps are present for a trial to be included...
            for k = 1:noObsSteps
                if ~all(any(trialNoObsStepIds==k,1)); missingObsOffStep = true; end
            end

            missingModStep = any(all(isnan(trialModStepIds),1));
            
            somethingWentWrong = missingModStep || missingControlStep || missingObsOffStep || isnan(contactPositions(j));
        catch
            somethingWentWrong = true;
        end


        % analyze trial if all steps are accounted for
        if somethingWentWrong
%             fprintf('  %s: missing steps in trial %i\n', session, j)
        else
            
            % find whether and where obstacle was toucheed
            isWheelBreak = any(breaks.times>obsOnTimes(j) & breaks.times<obsOffTimes(j));
            trialBinsTemp = frameTimeStamps>obsOnTimes(j) & frameTimeStamps<obsOffTimes(j);
            totalTouchFramesPerPaw = sum(touchesPerPaw(trialBinsTemp,:)>0,1); % get total number of obs touches in trial per paw
            
            % get trial velocity and body angle
            if isWheelBreak
                endTime = breaks.times(find(breaks.times>obsOnTimes(j),1,'first'));
            else
                endTime = obsOffTimes(j);
            end
            wheelInds = find(wheelTimes>=obsOnTimes(j) & wheelTimes<endTime);
            dp = wheelPositions(wheelInds(end)) - wheelPositions(wheelInds(1));
            dt = wheelTimes(wheelInds(end)) - wheelTimes(wheelInds(1));
            avgVel = dp/dt;
            avgAngle = nanmedian(bodyAngles(frameTimeStamps>obsOnTimes(j) & frameTimeStamps<endTime));
            
            
            % determine whether left and right forepaws are in swing at obsPos moment
            isLeftSwingAtContact = ~isnan(trialModStepIds(trialTimeStamps==0,2));
            isRightSwingAtContact = ~isnan(trialModStepIds(trialTimeStamps==0,3));

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
            % both only if at least one paw is in stance at moment of contact
            if ~(isLeftSwingAtContact && isRightSwingAtContact)
                if firstModPaw==2; stancePaw=3; else; stancePaw=2; end
                stanceDistance = trialLocations(trialTimeStamps==0,1,stancePaw);
            else
                stanceDistance = nan;
                stancePaw = nan;
            end
            
            % get distance of firstModPaw to obs at the beginning of first mod step
            swingStartDistance = trialLocations(find(trialModStepIds(:,firstModPaw)==1,1,'first'),1,firstModPaw);


            % get mod, control, and noObs step(s) length, duration, wheel velocity
            maxModSteps = max(trialModStepIds(:));
            
            controlSwingLengths = nan(controlSteps,4);
            noObsSwingLengths = nan(noObsSteps,4);
            modifiedSwingLengths = nan(maxModSteps,4);
            controlSwingDurations = nan(controlSteps,4);
            noObsSwingDurations = nan(noObsSteps,4);
            modifiedSwingDurations = nan(maxModSteps,4);
            controlWheelVels = nan(controlSteps,4);
            noObsWheelVels = nan(noObsSteps,4);
            modifiedWheelVels = nan(maxModSteps,4);
            

            for k = 1:4

                % modified steps
                for m = 1:maxModSteps
                    stepBins = trialModStepIds(:,k)==m;
                    if any(stepBins)
                        stepXLocations = trialLocations(stepBins,1,k);
                        modifiedSwingLengths(m,k) = stepXLocations(end) - stepXLocations(1);
                        stepTimes = trialTimeStamps(stepBins);
                        modifiedSwingDurations(m,k) = stepTimes(end) - stepTimes(1);
                        modifiedWheelVels(m,k) = trialWheelVel(find(stepBins,1,'first'));
                    end
                end

                % control steps
                for m = 1:controlSteps
                    stepBins = trialControlStepIds(:,k)==m;
                    stepXLocations = trialLocations(stepBins,1,k);
                    controlSwingLengths(m,k) = stepXLocations(end) - stepXLocations(1);
                    stepTimes = trialTimeStamps(stepBins);
                    controlSwingDurations(m,k) = stepTimes(end) - stepTimes(1);
                    controlWheelVels(m,k) = trialWheelVel(find(stepBins,1,'first'));
                end
                
                % no obs steps
                for m = 1:controlSteps
                    stepBins = trialNoObsStepIds(:,k)==m;
                    stepXLocations = trialLocations(stepBins,1,k);
                    noObsSwingLengths(m,k) = stepXLocations(end) - stepXLocations(1);
                    stepTimes = trialTimeStamps(stepBins);
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
            
            allLocations = {controlLocations, noObsLocations, modLocations};
            allLocationsInterp = {controlLocationsInterp, noObsLocationsInterp, modLocationsInterp};
            allIds = {trialControlStepIds, trialNoObsStepIds, trialModStepIds};

            for stepType = 1:3
                for k = 1:4

                    stepNum = max(allIds{stepType}(:,k));
                    pawLocations = nan(stepNum, 3, swingMaxSmps);
                    pawLocationsInterp = nan(stepNum, 3, interpSmps);

                    for m = 1:stepNum

                        % locations
                        stepBins = allIds{stepType}(:,k)==m;
                        startInd = find(stepBins, 1, 'first');
                        endInd = find(stepBins, 1, 'last');
                        if endInd-startInd>=swingMaxSmps; endInd = startInd+swingMaxSmps-1; end % make sure swing is not too long
        
                        
                        stepEndInd = endInd-startInd+1;
                        stepX = nan(1,swingMaxSmps); stepY = nan(1,swingMaxSmps); stepZ = nan(1,swingMaxSmps);
                        stepX(1:stepEndInd) = trialLocations(startInd:endInd,1,k); stepX(stepEndInd:end) = stepX(stepEndInd); % the latter statement ensures the kinematics don't bleed into the subsequent step (the step is 'frozen' at the moment swing ends)
                        stepY(1:stepEndInd) = trialLocations(startInd:endInd,2,k); stepY(stepEndInd:end) = stepY(stepEndInd);
                        stepZ(1:stepEndInd) = trialLocations(startInd:endInd,3,k); stepZ(stepEndInd:end) = stepZ(stepEndInd);
                        pawLocations(m,:,:) = cat(1,stepX,stepY,stepZ);

                        % locations interp
                        xInterp = interp1(1:sum(stepBins), trialLocations(stepBins,1,k), linspace(1,sum(stepBins),interpSmps));
                        yInterp = interp1(1:sum(stepBins), trialLocations(stepBins,2,k), linspace(1,sum(stepBins),interpSmps));
                        zInterp = interp1(1:sum(stepBins), trialLocations(stepBins,3,k), linspace(1,sum(stepBins),interpSmps));
                        pawLocationsInterp(m,:,:) = cat(1,xInterp,yInterp,zInterp);
                        
                        % things to do only for modified locations
                        if stepType==3
                            
                            modStepNum(k) = stepNum;
                            
                            % get ind of obs hit in interpolated coordinates
                            if m==1
                                stepObsPosInd = find(trialTimeStamps==0) - find(stepBins,1,'first') + 1;
                                pawObsPosIndInterp(k) = interp1(linspace(1,sum(stepBins), interpSmps), ...
                                    1:interpSmps, stepObsPosInd, 'nearest');
                                pawObsPosInd(k) = find(trialTimeStamps==0) - find(stepBins,1,'first') + 1;
                            end
                        end
                    end

                    allLocations{stepType}{k} = pawLocations;
                    allLocationsInterp{stepType}{k} = pawLocationsInterp;
                end
            end

            
            % STORE RESULTS
            
            % session metadata (determined from sessionInfo table)
            for k = metaDataFields
                sessionData(dataInd).(k{1}) = sessionMetaData.(k{1});
            end
            
            % trial metadata
            sessionData(dataInd).trial = j;
            sessionData(dataInd).isLightOn = isLightOn(j);
            sessionData(dataInd).obsHeightsVid = obsHeightsVid(j);
            if ismember('side', metaDataFields) % if a neural manipulation occured, this stores information about the side of the brain
                sessionData(dataInd).contraPawFirst = firstPawOver==contraLimb;
                sessionData(dataInd).ipsiPawFirst = firstPawOver~=contraLimb && ~isnan(contraLimb);
            end

            

            % bunch of thangs
            sessionData(dataInd).vel = interp1(wheelTimes, vel, contactTimes(j));  % mouse vel at moment of wisk contact
            sessionData(dataInd).angle = interp1(frameTimeStamps(~isnan(frameTimeStamps)), bodyAngles(~isnan(frameTimeStamps)), contactTimes(j));  % body angle at moment of wisk contact
            sessionData(dataInd).avgVel = avgVel;
            sessionData(dataInd).avgAngle = avgAngle;
            sessionData(dataInd).obsPos = contactPositions(j);       % position of obs relative to nose at moment of wisk contact
            sessionData(dataInd).obsPosInd = find(trialTimeStamps==0); % ind at which obs contacts wisks for trial
            sessionData(dataInd).pawObsPosInd = pawObsPosInd;% ind at which obs contacts wisks for locations for each paw
            sessionData(dataInd).pawObsPosIndInterp = pawObsPosIndInterp; % ind at which obs contacts wisks for interp locations for each paw
            sessionData(dataInd).timeStamps = trialTimeStamps;
            sessionData(dataInd).frameTimeStamps = frameTimeStamps(trialInds); % in original frameTimeStamps reference point
            sessionData(dataInd).locations = trialLocations;
            sessionData(dataInd).controlLocations = allLocations{1};
            sessionData(dataInd).noObsLocations = allLocations{2};
            sessionData(dataInd).modifiedLocations = allLocations{3};
            sessionData(dataInd).controlLocationsInterp = allLocationsInterp{1};
            sessionData(dataInd).noObsLocationsInterp = allLocationsInterp{2};
            sessionData(dataInd).modifiedLocationsInterp = allLocationsInterp{3};
            sessionData(dataInd).trialControlStepIdentities = trialControlStepIds;
            sessionData(dataInd).modifiedStepIdentities = trialModStepIds;
            sessionData(dataInd).modStepNum = modStepNum;
            sessionData(dataInd).stanceDistance = stanceDistance;
            sessionData(dataInd).stancePaw = stancePaw;
            sessionData(dataInd).swingStartDistance = swingStartDistance;
            
            % trial touch info
            sessionData(dataInd).totalTouchFramesPerPaw = totalTouchFramesPerPaw;
            sessionData(dataInd).trialTouchesPerPaw = trialTouchesPerPaw;
            sessionData(dataInd).isWheelBreak = isWheelBreak;
            
            
            % info about which paws did what
            sessionData(dataInd).isLeftSwingAtContact = isLeftSwingAtContact;
            sessionData(dataInd).isRightSwingAtContact = isRightSwingAtContact;
            sessionData(dataInd).firstPawOver = firstPawOver;
            sessionData(dataInd).firstModPaw = firstModPaw;

            % swing characteristics
            sessionData(dataInd).controlSwingLengths = controlSwingLengths;
            sessionData(dataInd).modifiedSwingLengths = modifiedSwingLengths;
            sessionData(dataInd).noObsSwingLengths = noObsSwingLengths;
            sessionData(dataInd).controlSwingDurations = controlSwingDurations;
            sessionData(dataInd).modifiedSwingDurations = modifiedSwingDurations;
            sessionData(dataInd).noObsSwingDurations = noObsSwingDurations;
            sessionData(dataInd).controlWheelVels = controlWheelVels;
            sessionData(dataInd).modifiedWheelVels = modifiedWheelVels;
            sessionData(dataInd).noObsWheelVels = noObsWheelVels;

            dataInd = dataInd + 1;
        end
    end
    if timeOperations; fprintf('getting data for all trials: %i seconds\n', round(toc)); end
    fprintf('%s: finished with %.2f of trials successfully analyzed\n', ...
        session, length(sessionData)/length(obsOnTimes));
end
disp('all done!');
end








