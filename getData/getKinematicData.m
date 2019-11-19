function [kinData, stanceBins, models] = getKinematicData(session, obsPos)


% TO DO: deal with unexpected occurences

% settings
speedTime = .05; % compute velocity over this interval
interpSmps = 100; % strides are stretched to have same number of samples // interpSmps sets the number of samples per interpolated stride
swingMaxSmps = 50; % when averaging swing locations without interpolating don't take more than swingMaxSmps for each swing
noObsSteps = 3;
controlSteps = 2; % needs to be at least 2
contactPosLimits = [-.02 .02]; % whisker cant only contact obs this far in front of and behind nose
timeOperations = false;


% load session data
fprintf('%s: getting kinematic data\n', session);
if timeOperations; tic; end
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'))
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
kinData(length(obsOnTimes)) = struct();

[locations, features] = fixTracking(locationsTable, frameTimeStamps);
botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, wheelTimes, wheelCenter, wheelRadius, 250, pixelsPerM);
if timeOperations; fprintf('loading session data and fixing tracking: %i seconds\n', round(toc)); end



% determine whisker contact positions and times (only if obsPos provided, which means whiskers were trimmed for the sessions // otherwise, they are inherited from runAnalyzed)
if timeOperations; tic; end
if exist('obsPos', 'var')
    wiskContactPositions = ones(size(obsOnTimes))*obsPos;
    wiskContactTimes = nan(size(obsOnTimes));
    
    for j = 1:length(obsOnTimes)
        indStart = find(obsPositionsFixed>=wiskContactPositions(j) & obsTimes>obsOnTimes(j) & obsTimes<obsOffTimes(j), 1, 'first');
        if ~isempty(indStart)
            wiskContactTimes(j) = interp1(obsPositionsFixed(indStart-1:indStart), obsTimes(indStart-1:indStart), wiskContactPositions(j));
            wiskContactTimes(j) = frameTimeStamps(knnsearch(frameTimeStamps, wiskContactTimes(j))); % force time to be in frameTimeStamps
        end
    end
end
wiskContactPositions(wiskContactPositions<contactPosLimits(1) | wiskContactPositions>contactPosLimits(2)) = nan; % remove invalid positions
if timeOperations; fprintf('getting wisk contact pos/times: %i seconds\n', round(toc)); end



% get step identities and wheel velocity
if timeOperations; tic; end
obsPixPositionsContinuous = getObsPixPositionsContinuous(...
    wheelToObsPixPosMappings, wheelTimes, wheelPositions, frameTimeStamps, ...
    obsPixPositions, obsPixPositionsUninterped, obsOnTimes, obsOffTimes);
[controlStepIdentities, modifiedStepIdentities, noObsStepIdentities] = getStepIdentities(...
    stanceBins, locations(:,:,botPawInds), wiskContactTimes, frameTimeStamps, ...
    obsOnTimes, obsOffTimes, obsPixPositions, obsPixPositionsContinuous, controlSteps, noObsSteps);
wheelVel = getVelocity(wheelPositions, speedTime, targetFs);
if timeOperations; fprintf('getting step identities: %i seconds\n', round(toc)); end


% put together xyz for paws
if timeOperations; tic; end
locationsPaws = nan(size(locations,1), 3, 4);
locationsPaws(:,1:2,:) = locations(:,:,botPawInds);
locationsPaws(:,3,:) = locations(:,2,topPawInds);
locationsPaws(:,2,:) = locationsPaws(:,2,:) - nosePos(2); % subtract midline from all y values
locationsPaws(:,3,:) = (wheelCenter(2)-wheelRadius) - locationsPaws(:,3,:); % flip z and set s.t. top of wheel is zero

% put together xyz for tail
botTailBins = contains(features, 'tail') & contains(features, '_bot');
topTailBins = contains(features, 'tail') & contains(features, '_top');
locationsTail = nan(size(locations,1), 3, 2); % frameNum X xyz X tailbase/mid
locationsTail(:,1:2,:) = locations(:,:,botTailBins);
locationsTail(:,3,:) = locations(:,2,topTailBins);
locationsTail(:,2,:) = locationsTail(:,2,:) - nosePos(2); % subtract midline from all y values
locationsTail(:,3,:) = (wheelCenter(2)-wheelRadius) - locationsTail(:,3,:); % flip z and set s.t. top of wheel is zero

% unravel x coordinates s.t. 0 is position of obs and x coordinates move forward over time ('unheadfixing')
% this is done trial by trial, with each trial ending at the last ind of the last mod step over the obs
prevInd = 1;
for j = 1:length(obsOnTimes)
    finalInd = find(frameTimeStamps>obsOffTimes(j) & all(isnan(modifiedStepIdentities),2), 1, 'first');
    if j==length(obsOnTimes); finalInd = length(frameTimeStamps); end
    
    locationsPaws(prevInd:finalInd,1,:) = locationsPaws(prevInd:finalInd,1,:) - obsPixPositionsContinuous(j, prevInd:finalInd)';
    locationsTail(prevInd:finalInd,1,:) = locationsTail(prevInd:finalInd,1,:) - obsPixPositionsContinuous(j, prevInd:finalInd)';
    prevInd = finalInd+1;
end
locationsPaws = locationsPaws / pixelsPerM; % convert to meters
locationsTail = locationsTail / pixelsPerM; % convert to meters
if timeOperations; fprintf('getting xyz coords: %i seconds\n', round(toc)); end







% collect data for each trial
if timeOperations; tic; end
isTrialAnalyzed = false(1,length(kinData)); % keeps track of whereh kinematic data could be analyzed for each trial
for j = 1:length(obsOnTimes)
        
    % first pull out trial specific data because 'find' function works much quicker on the smaller data slices

    % get trial inds
    if j==1; lastTrialEndTime=0; else; lastTrialEndTime = obsOffTimes(j-1); end
    startInd = find(frameTimeStamps>lastTrialEndTime & any(noObsStepIdentities==1,2), 1, 'first'); % start at first obs off step after previous trial
    endInd = find(all(isnan(modifiedStepIdentities),2) & frameTimeStamps>obsOffTimes(j), 1, 'first'); % end when all mod steps in trial are over
    trialInds = startInd:endInd;

    % get trial data
    trialControlStepIds = controlStepIdentities(trialInds,:);
    trialModStepIds = modifiedStepIdentities(trialInds,:);
    trialNoObsStepIds = noObsStepIdentities(trialInds,:);
    trialLocations = locationsPaws(trialInds,:,:);
    trialLocationsTail = locationsTail(trialInds,:,:);
    trialWheelVel = interp1(wheelTimes, wheelVel, frameTimeStamps(trialInds));

    % determine whether kinematic data can be analyzed for trial (all )
    [missingControlStep, missingObsOffStep] = deal(false);
    for k = 1:controlSteps; if ~all(any(trialControlStepIds==k,1)); missingControlStep = true; end; end
    for k = 1:noObsSteps; if ~all(any(trialNoObsStepIds==k,1)); missingObsOffStep = true; end; end
    missingModStep = any(all(isnan(trialModStepIds),1));
    wiskContactFrameFound = any(frameTimeStamps(trialInds)==wiskContactTimes(j)); % 
    isTrialAnalyzed(j) = ~missingModStep && ~missingControlStep && ~missingObsOffStep && ...
                         ~isnan(wiskContactPositions(j)) && wiskContactFrameFound;
    
    % ANALYZE KINEMATIC DATA
    if isTrialAnalyzed(j)
        % determine whether left and right forepaws are in swing at obsPos moment
        contactInd = find(frameTimeStamps(trialInds)==wiskContactTimes(j)); % ind within trial at which contact occurs
        isLeftSwingAtContact = ~isnan(trialModStepIds(contactInd,2));
        isRightSwingAtContact = ~isnan(trialModStepIds(contactInd,3));

        % determine sequence with which paws cross obstacle
        isStepping = ~isnan(trialModStepIds);
        lastModStepInds = table2array(rowfun(@(x)(find(x,1,'last')), table(isStepping'))); % final ind of final modified step for each paw
        [~, pawOverSequence] = sort(lastModStepInds);


        % determine which is paw first modified step (forepaws only)
        %
        % if both one in swing and one in stance at contact, first mod step is one
        % in swing // if both in stance at moment of contact, first mod
        % step is first to enter swing // if both in swing, first mod
        % step is first paw to land on the other side
        if xor(isLeftSwingAtContact, isRightSwingAtContact)
            if isLeftSwingAtContact; firstModPaw = 2; else; firstModPaw = 3; end
        elseif isLeftSwingAtContact && isRightSwingAtContact
            firstModPaw = pawOverSequence(1);
        elseif ~isLeftSwingAtContact && ~isRightSwingAtContact
            firstModStepInds = table2array(rowfun(@(x)(find(x,1,'first')), table(isStepping')));
            [~, firstModPaw] = min(firstModStepInds .* [nan 1 1 nan]'); % mask out hind paws, inds 1 and 4
        end


        % get mod, control, and noObs step(s) length, duration, wheel velocity
        maxModSteps = max(trialModStepIds(:));
        [controlSwingLengths, controlSwingDurations, controlWheelVels] = deal(nan(controlSteps,4));
        [noObsSwingLengths, noObsSwingDurations, noObsWheelVels] = deal(nan(noObsSteps,4));
        [modifiedSwingLengths, modifiedSwingDurations, modifiedWheelVels] = deal(nan(maxModSteps,4));

        for k = 1:4

            % modified steps
            for m = 1:maxModSteps
                stepBins = trialModStepIds(:,k)==m;
                if any(stepBins)
                    stepXLocations = trialLocations(stepBins,1,k);
                    modifiedSwingLengths(m,k) = stepXLocations(end) - stepXLocations(1);
                    stepTimes = frameTimeStamps(trialInds(stepBins));
                    modifiedSwingDurations(m,k) = stepTimes(end) - stepTimes(1);
                    modifiedWheelVels(m,k) = trialWheelVel(find(stepBins,1,'first'));
                end
            end

            % control steps
            for m = 1:controlSteps
                stepBins = trialControlStepIds(:,k)==m;
                stepXLocations = trialLocations(stepBins,1,k);
                controlSwingLengths(m,k) = stepXLocations(end) - stepXLocations(1);
                stepTimes = frameTimeStamps(trialInds(stepBins));
                controlSwingDurations(m,k) = stepTimes(end) - stepTimes(1);
                controlWheelVels(m,k) = trialWheelVel(find(stepBins,1,'first'));
            end

            % no obs steps
            for m = 1:controlSteps
                stepBins = trialNoObsStepIds(:,k)==m;
                stepXLocations = trialLocations(stepBins,1,k);
                noObsSwingLengths(m,k) = stepXLocations(end) - stepXLocations(1);
                stepTimes = frameTimeStamps(trialInds(stepBins));
                noObsSwingDurations(m,k) = stepTimes(end) - stepTimes(1);
                noObsWheelVels(m,k) = trialWheelVel(find(stepBins,1,'first'));
            end
        end



        % get noObs, control, and mod paw locations (interpolated and non-interpolated)
        [pawObsPosInd, pawObsPosIndInterp] = deal(nan(1,4));
        [allLocations, allLocationsInterp] = deal({cell(1,4), cell(1,4), cell(1,4)}); % three cell arrays containing data for control, noObs, and modLocations, respectively
        allIds = {trialControlStepIds, trialNoObsStepIds, trialModStepIds};

        for stepType = 1:3
            for k = 1:4

                stepNum = max(allIds{stepType}(:,k));
                pawLocations = nan(stepNum, 3, swingMaxSmps);
                pawLocationsInterp = nan(stepNum, 3, interpSmps);

                for m = 1:stepNum

                    try
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
                        if stepType==3 && m==1 % get ind of obs hit in interpolated coordinates only for first mod steps
                            pawObsPosInd(k) = contactInd - find(stepBins,1,'first') + 1;
                            pawObsPosIndInterp(k) = interp1(linspace(1,sum(stepBins), interpSmps), ...
                                1:interpSmps, pawObsPosInd(k), 'nearest');
                        end
                    catch
                        isTrialAnalyzed(j) = false;
                        fprintf('WARNING! Unable to analyze trial %i!\n', j);
                    end
                end

                allLocations{stepType}{k} = pawLocations;
                allLocationsInterp{stepType}{k} = pawLocationsInterp;
            end
        end


        % save trial kinematics
        kinData(j).trialInds = trialInds;
        kinData(j).locations = trialLocations;
        kinData(j).locationsTail = trialLocationsTail;
        kinData(j).wiskContactPositions = wiskContactPositions(j);
        kinData(j).wiskContactTimes = wiskContactTimes(j);
        kinData(j).contactInd = contactInd;
        
        kinData(j).controlLocations = allLocations{1};
        kinData(j).noObsLocations = allLocations{2};
        kinData(j).modifiedLocations = allLocations{3};
        kinData(j).controlLocationsInterp = allLocationsInterp{1};
        kinData(j).noObsLocationsInterp = allLocationsInterp{2};
        kinData(j).modifiedLocationsInterp = allLocationsInterp{3};
        kinData(j).controlStepIdentities = trialControlStepIds;
        kinData(j).noObsStepIdentities = trialNoObsStepIds;
        kinData(j).modifiedStepIdentities = trialModStepIds;
        kinData(j).stanceBins = stanceBins(trialInds,:);
        kinData(j).pawObsPosInd = pawObsPosInd;% ind in first mod paws at which obs contacts wisks for locations for each paw
        kinData(j).pawObsPosIndInterp = pawObsPosIndInterp; % ind in first mod paws at which obs contacts wisks for interp locations for each paw

        % info about which paws did what
        kinData(j).isLeftSwingAtContact = isLeftSwingAtContact;
        kinData(j).isRightSwingAtContact = isRightSwingAtContact;
        kinData(j).pawOverSequence = pawOverSequence;
        kinData(j).firstModPaw = firstModPaw;

        % save swing characteristics
        kinData(j).controlSwingLengths = controlSwingLengths;
        kinData(j).modifiedSwingLengths = modifiedSwingLengths;
        kinData(j).noObsSwingLengths = noObsSwingLengths;
        kinData(j).controlSwingDurations = controlSwingDurations;
        kinData(j).modifiedSwingDurations = modifiedSwingDurations;
        kinData(j).noObsSwingDurations = noObsSwingDurations;
        kinData(j).controlWheelVels = controlWheelVels;
        kinData(j).modifiedWheelVels = modifiedWheelVels;
        kinData(j).noObsWheelVels = noObsWheelVels;
    end
    kinData(j).isTrialAnalyzed = isTrialAnalyzed(j);
    if ~isTrialAnalyzed(j) % still useful to compute trialInds even when kinematics can't be analyzed...
        kinData(j).trialInds = find(frameTimeStamps>obsOnTimes(j) & frameTimeStamps<obsOffTimes(j));
    end
end
if timeOperations; fprintf('getting data for all trials: %i seconds\n', round(toc)); end
fprintf('%s: got kinematicData with %.2f of trials successfully analyzed\n', ...
    session, mean(isTrialAnalyzed));




try
    % MAKE SWING LENGTH PREDICTIVE MODELS
    % make model to predict would-be mod swing length of first modified paw using wheel vel and previous swing lengths as predictors

    % generate models
    models = cell(1,4);
    controlVels = cat(1,kinData.controlWheelVels);
    controlLengths = cat(1,kinData.controlSwingLengths);
    for paw = 1:4
        validInds = ~isnan(controlVels(:,paw)) & ~isnan(controlLengths(:,paw));
        models{paw} = fitlm(controlVels(validInds,paw), controlLengths(validInds,paw), 'Linear', 'RobustOpts', 'on');
    end

    % make predictions
    for j = find(isTrialAnalyzed)
        kinData(j).controlPredictedLengths = nan(size(kinData(j).controlWheelVels));
        kinData(j).modPredictedLengths = nan(size(kinData(j).modifiedWheelVels));
        for paw = 1:4
            kinData(j).controlPredictedLengths(:,paw) = predict(models{paw}, kinData(j).controlWheelVels(:,paw));
            kinData(j).modPredictedLengths(:,paw) = predict(models{paw}, kinData(j).modifiedWheelVels(:,paw));
        end
    end

    save(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'kinData.mat'), 'kinData', 'stanceBins', 'models', '-v7.3')


    % uncomment the following to show fits for step length for each paw
%     figure('name', [session ' paw length fits'], 'Position', [2000 300 700 600], 'color', 'white', 'menubar', 'none')
%     colors = hsv(4);
%     for i = 1:4
%         scatter(controlVels(:,i), controlLengths(:,i), 20, colors(i,:), 'filled', 'markerfacealpha', .4); hold on;
%         xLims = get(gca, 'XLim');
%         coefs = models{i}.Coefficients.Estimate;
%         plot(xLims, xLims*coefs(2)+coefs(1), 'Color', colors(i,:), 'LineWidth', 2);
%     end
%     xlabel('velocity (m/s)')
%     ylabel('step length')
%     for i = 1:4; lines(i) = plot([nan nan], 'color', colors(i,:), 'LineWidth', 2); end % create dummy lines
%     legend(lines, {'LH', 'LF', 'RF', 'RH'}, 'Location', 'best', 'box', 'off');

catch
    fprintf('%s: failed to make swing length model!\n', session)
    save(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'kinData.mat'), ...
        'kinData', 'stanceBins', 'models', '-v7.3')
end








