function [kinData, stanceBins, models] = getKinematicData(session)


% after analyzeSession computes low level information for each session,
% getKinematicData is used to extract more processed kinematic information,
% e.g. parsing each trial into interpolated kinematics for steps over the
% obstacle, steps preceding these steps, and steps when the obstacle is not
% engaged


% settings
speedTime = .05;         % compute velocity over this interval
interpSmps = 100;        % strides are stretched to have same number of samples // interpSmps sets the number of samples per interpolated stride
swingMaxSmps = 50;       % when averaging swing locations without interpolating don't take more than swingMaxSmps for each swing
noObsSteps = 3;          % how many steps per trial per paw to take before the obstacle becomes engaged
controlSteps = 2;        % needs to be at least 2 // how many steps per trial per paw to take before the first modified step
timeOperations = false;  % whether to report time it takes to compute differnt parts of this script
showStepLengthFits = false;  % whether to show linear first for step lengths for each paw


% load session data
fprintf('%s: getting kinematic data\n', session);
if timeOperations; tic; end
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'))
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
kinData(length(obsOnTimes)) = struct();
[locations, features] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM);
botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, wheelTimes, wheelCenter, wheelRadius, 250, pixelsPerM);
if timeOperations; fprintf('loading session data and fixing tracking: %i seconds\n', round(toc)); end

% get step identities and wheel velocity
if timeOperations; tic; end
obsPixPositionsContinuous = getObsPixPositionsContinuous(...
    wheelToObsPixPosMappings, wheelTimes, wheelPositions, frameTimeStamps, ...
    obsPixPositions, obsPixPositionsUninterped, obsOnTimes, obsOffTimes);
[controlStepIdentities, modifiedStepIdentities, noObsStepIdentities, trialStartInds] = getStepIdentities(...
    stanceBins, locations(:,:,botPawInds), wiskContactTimes, frameTimeStamps, ...
    obsOnTimes, obsOffTimes, obsPixPositionsContinuous, controlSteps, noObsSteps);
wheelVel = getVelocity(wheelPositions, speedTime, targetFs);
if timeOperations; fprintf('getting step identities: %i seconds\n', round(toc)); end

% put together xyz for paws
if timeOperations; tic; end
locationsPaws = nan(size(locations,1), 3, 4);  % frameNum X xyz X pawNum
locationsPaws(:,1:2,:) = locations(:,:,botPawInds);
locationsPaws(:,3,:) = locations(:,2,topPawInds);
locationsPaws(:,2,:) = locationsPaws(:,2,:) - nosePos(2); % subtract midline from all y values
locationsPaws(:,3,:) = (wheelCenter(2)-wheelRadius) - locationsPaws(:,3,:); % flip z and set s.t. top of wheel is zero

% put together xyz for tail
botTailBins = contains(features, 'tail') & contains(features, '_bot');
topTailBins = contains(features, 'tail') & contains(features, '_top');
locationsTail = nan(size(locations,1), 3, 2);  % frameNum X xyz X tailbase/mid
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
        
    % first get trial specific data because 'find' function works much quicker on the smaller data slices
    endInd = find(all(isnan(modifiedStepIdentities),2) & frameTimeStamps>obsOffTimes(j), 1, 'first'); % end when all mod steps in trial are over
    trialInds = trialStartInds(j):endInd;

    % get trial data
    trialTimes = frameTimeStamps(trialInds);
    trialControlStepIds = controlStepIdentities(trialInds,:);
    trialModStepIds = modifiedStepIdentities(trialInds,:);
    trialNoObsStepIds = noObsStepIdentities(trialInds,:);
    trialLocations = locationsPaws(trialInds,:,:);
    trialLocationsTail = locationsTail(trialInds,:,:);
    trialWheelVel = interp1(wheelTimes, wheelVel, frameTimeStamps(trialInds));

    % determine whether kinematic data can be analyzed for trial
    missingControlStep = ~all(max(trialControlStepIds, [], 1) >= controlSteps);
    missingNoObsStep = ~all(max(trialNoObsStepIds, [], 1) >= noObsSteps);
    missingModStep = any(all(isnan(trialModStepIds),1));
    wiskContactFrameFound = any(trialTimes==wiskContactTimes(j));
    isTrialAnalyzed(j) = ~missingModStep && ~missingControlStep && ~missingNoObsStep && ...
                         ~isnan(wiskContactPositions(j)) && wiskContactFrameFound;
    
    % ANALYZE KINEMATIC DATA
    if isTrialAnalyzed(j)
        try
            % determine whether left and right forepaws are in swing at obsPos moment
            contactInd = find(trialTimes==wiskContactTimes(j)); % ind within trial at which contact occurs
            isSwinging = ~isnan(trialModStepIds);
            isLeftSwingAtContact = isSwinging(contactInd,2);
            isRightSwingAtContact = isSwinging(contactInd,3);
            

            % determine sequence with which paws cross obstacle
            lastModStepInds = nan(4,1);
            for k = 1:4; lastModStepInds(k) = find(isSwinging(:,k),1,'last'); end
            [~, pawOverSequence] = sort(lastModStepInds);
            

            % determine 'first modified paw' (forepaws only)
            %
            % if one in swing and one in stance at contact, first mod step is 
            % one in swing // otherwise first mod paw is first paw to land on
            % the other side
            if xor(isLeftSwingAtContact, isRightSwingAtContact)
                if isLeftSwingAtContact; firstModPaw = 2; else; firstModPaw = 3; end
            else
                firstModPaw = pawOverSequence(1);
            end
            

            % get mod, control, and noObs step(s) length, duration, wheel velocity
            maxModSteps = max(trialModStepIds(:));
            allIds = {trialControlStepIds, trialNoObsStepIds, trialModStepIds};
            [allLengths, allDurations, allVels] = deal({nan(controlSteps,4), nan(noObsSteps,4), nan(maxModSteps,4)});

            for stepType = 1:3
                for k = 1:4
                    stepNum = max(allIds{stepType}(:,k));
                    for m = 1:stepNum
                        stepBins = allIds{stepType}(:,k)==m;
                        if any(stepBins)
                            allLengths{stepType}(m,k) = range(trialLocations(stepBins,1,k));
                            allDurations{stepType}(m,k) = range(trialTimes(stepBins));
                            allVels{stepType}(m,k) = trialWheelVel(find(stepBins,1,'first'));
                        end
                    end
                end
            end


            % get noObs, control, and mod paw locations (interpolated and non-interpolated)
            [pawObsPosInd, pawObsPosIndInterp] = deal(nan(1,4));
            [allLocations, allLocationsInterp] = deal({cell(1,4), cell(1,4), cell(1,4)});  % three cell arrays containing data for control, noObs, and modLocations, respectively

            for stepType = 1:3
                for k = 1:4

                    stepNum = max(allIds{stepType}(:,k));
                    pawLocations = nan(stepNum, 3, swingMaxSmps);
                    pawLocationsInterp = nan(stepNum, 3, interpSmps);

                    for m = 1:stepNum
                        % find inds for stepType, paw, and stepNum
                        stepBins = allIds{stepType}(:,k)==m;
                        startInd = find(stepBins, 1, 'first');
                        endInd = find(stepBins, 1, 'last');
                        if endInd-startInd>=swingMaxSmps; endInd = startInd+swingMaxSmps-1; end % make sure swing is not too long

                        % locations
                        stepEndInd = endInd-startInd+1;
                        [stepX, stepY, stepZ] = deal(nan(1,swingMaxSmps));
                        stepX(1:stepEndInd) = trialLocations(startInd:endInd,1,k);
                        stepY(1:stepEndInd) = trialLocations(startInd:endInd,2,k);
                        stepZ(1:stepEndInd) = trialLocations(startInd:endInd,3,k);
                        pawLocations(m,:,:) = fillmissing([stepX; stepY; stepZ], 'previous', 2);  % fillmissing replaces the nans at the end with the last non-nan entry, 'freezing' the kinematics at this value

                        % locations interpolated
                        xInterp = interp1(1:sum(stepBins), trialLocations(stepBins,1,k), linspace(1,sum(stepBins),interpSmps));
                        yInterp = interp1(1:sum(stepBins), trialLocations(stepBins,2,k), linspace(1,sum(stepBins),interpSmps));
                        zInterp = interp1(1:sum(stepBins), trialLocations(stepBins,3,k), linspace(1,sum(stepBins),interpSmps));
                        pawLocationsInterp(m,:,:) = [xInterp; yInterp; zInterp];

                        % things to do only for modified locations
                        if stepType==3 && m==1 % get ind of obs hit in interpolated coordinates only for first mod steps
                            pawObsPosInd(k) = contactInd - startInd + 1;
                            pawObsPosIndInterp(k) = round(interp1(linspace(1,sum(stepBins),interpSmps), ...
                                1:interpSmps, pawObsPosInd(k), 'nearest'));
                        end
                    end

                    allLocations{stepType}{k} = pawLocations;
                    allLocationsInterp{stepType}{k} = pawLocationsInterp;
                end
            end


            % save trial kinematics
            kinData(j).trialInds = trialInds;  % inds for frames included in trial
            kinData(j).locations = trialLocations;  % paw kinematics for trial (frame X xyz X pawNum)
            kinData(j).locationsTail = trialLocationsTail;  % tail kinematics for trial (frame X xyz X tailPosition)
            kinData(j).wiskContactPositions = wiskContactPositions(j);  % position relative to nose at which whiskers touch obstacle
            kinData(j).wiskContactTimes = wiskContactTimes(j);  % time at which whiskers touch obstacle
            kinData(j).contactInd = contactInd;  % frame index within trial at which whiskers contact obstacle

            % kinematics
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
            kinData(j).pawObsPosInd = pawObsPosInd;  % ind in first mod paws at which obs contacts wisks for locations for each paw
            kinData(j).pawObsPosIndInterp = pawObsPosIndInterp;  % ind in first mod paws at which obs contacts wisks for interp locations for each paw

            % which paws did what
            kinData(j).isLeftSwingAtContact = isLeftSwingAtContact;
            kinData(j).isRightSwingAtContact = isRightSwingAtContact;
            kinData(j).pawOverSequence = pawOverSequence;
            kinData(j).firstModPaw = firstModPaw;

            % swing characteristics
            kinData(j).controlSwingLengths = allLengths{1};
            kinData(j).noObsSwingLengths = allLengths{2};
            kinData(j).modifiedSwingLengths = allLengths{3};
            kinData(j).controlSwingDurations = allDurations{1};
            kinData(j).noObsSwingDurations = allDurations{2};
            kinData(j).modifiedSwingDurations = allDurations{3};
            kinData(j).controlWheelVels = allVels{1};
            kinData(j).noObsWheelVels = allVels{2};
            kinData(j).modifiedWheelVels = allVels{3};
        catch
            isTrialAnalyzed(j) = false;
            fprintf('WARNING! Unable to analyze trial %i!\n', j);
            kinData(j).trialInds = find(frameTimeStamps>obsOnTimes(j) & frameTimeStamps<obsOffTimes(j)); % still useful to compute trialInds even when kinematics can't be analyzed...
        end
    end
    
    kinData(j).isTrialAnalyzed = isTrialAnalyzed(j);
end
if timeOperations; fprintf('getting data for all trials: %i seconds\n', round(toc)); end
fprintf('%s: got kinematicData with %.2f of trials successfully analyzed\n', ...
    session, mean(isTrialAnalyzed));




try
    % MAKE SWING LENGTH PREDICTIVE MODELS
    % make model to predict would-be mod swing length from wheel vel at start of swing

    % generate models
    models = cell(1,4);
    controlVels = cat(1, kinData.controlWheelVels);
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

    % show model fits
    if showStepLengthFits
        figure('name', [session ' paw length fits'], 'Position', [2000 300 700 600], 'color', 'white')
        colors = hsv(4);
        for i = 1:4
            scatter(controlVels(:,i), controlLengths(:,i), 20, colors(i,:), 'filled', 'markerfacealpha', .4); hold on;
            xLims = get(gca, 'XLim');
            coefs = models{i}.Coefficients.Estimate;
            plot(xLims, xLims*coefs(2)+coefs(1), 'Color', colors(i,:), 'LineWidth', 2);
        end
        xlabel('velocity (m/s)')
        ylabel('step length')
        for i = 1:4; lines(i) = plot([nan nan], 'color', colors(i,:), 'LineWidth', 2); end % create dummy lines
        legend(lines, {'LH', 'LF', 'RF', 'RH'}, 'Location', 'best', 'box', 'off');
        pause(.001)
    end
catch
    fprintf('%s: failed to make swing length model!\n', session)
end

save(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'kinData.mat'), 'kinData', 'stanceBins', 'models', '-v7.3')



