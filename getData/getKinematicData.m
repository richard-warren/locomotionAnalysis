function [kinData, stanceBins, models] = getKinematicData(session, varargin)

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
s.showStepLengthFits = false;  % whether to show linear first for step lengths for each paw


% load session data
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
fprintf('%s: getting kinematic data\n', session);
if timeOperations; tic; end
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'))
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
kinData(length(obsOnTimes)) = struct();
scoreThresh = getScoreThresh(session, 'trackedFeaturesRaw_metadata.mat');  % scoreThresh depends on whether deeplabcut (old version) or deepposekit was used
[locations, features] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM, 'scoreThresh', scoreThresh);
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
locationsPawsPixels = locationsPaws;  % copy version in pixel coordinates before making subsequent transformations
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
if timeOperations; fprintf('getting xyz coords: %i seconds\n', round(toc)); end




% collect data for each trial
if timeOperations; tic; end
isTrialAnalyzed = true(1,length(kinData)); % keeps track of whereh kinematic data could be analyzed for each trial

for j = 1:length(obsOnTimes)
        
    % find trial indices
    if j==length(obsOnTimes); maxTime=frameTimeStamps(end); else; maxTime = obsOnTimes(j+1); end
    endInd = find(frameTimeStamps>obsOffTimes(j) & ...  % after obs is off
                  frameTimeStamps<maxTime & ...         % before next obs is on
                  ~any(controlStepIdentities==1 | noObsStepIdentities==1, 2), 1, 'last');  % before subsequent trial steps apear
    trialInds = trialStartInds(j):endInd;
    
    % unravel x coordinates s.t. 0 is position of obs and x coordinates move forward over time ('unheadfixing')
    % this is done trial by trial, with each trial ending at the last ind of the last mod step over the obs
    trialLocations = locationsPaws(trialInds,:,:);
    trialLocationsTail = locationsTail(trialInds,:,:);

    trialLocations(:,1,:) = trialLocations(:,1,:) - obsPixPositionsContinuous(j, trialInds)';
    trialLocationsTail(:,1,:) = trialLocationsTail(:,1,:) - obsPixPositionsContinuous(j, trialInds)';
    
    trialLocations = trialLocations / pixelsPerM; % convert to meters
    trialLocationsTail = trialLocationsTail / pixelsPerM; % convert to meters


    % first get trial specific data because 'find' function works much quicker on the smaller data slices
    trialTimes = frameTimeStamps(trialInds);
    trialControlStepIds = controlStepIdentities(trialInds,:);
    trialModStepIds = modifiedStepIdentities(trialInds,:);
    trialNoObsStepIds = noObsStepIdentities(trialInds,:);
    trialWheelVel = interp1(wheelTimes, wheelVel, frameTimeStamps(trialInds));
    trialObsPixPos = obsPixPositionsContinuous(j, trialInds);

    
    % get kinematic data for trial
    try
        % determine whether left and right forepaws are in swing at obsPos moment
        contactInd = knnsearch(trialTimes, wiskContactTimes(j)); % ind within trial at which contact occurs
        if isnan(wiskContactTimes(j)); contactInd = nan; end  % this is a hack to account for the fact that knnsearch returns 1 for nan values... this assignment will throw an error for the trial, which is the desired behavior
        isSwinging = ~isnan(trialModStepIds);
        isLeftSwingAtContact = isSwinging(contactInd,2);
        isRightSwingAtContact = isSwinging(contactInd,3);


        % determine sequence with which paws cross obstacle
        lastModStepInds = nan(4,1);
        for k = 1:4; lastModStepInds(k) = find(isSwinging(:,k),1,'last'); end
        [~, pawOverSequence] = sort(lastModStepInds);


        % determine 'first modified paw' (forepaws only)
        
        % if one swing, one stance: pick the one in swing
        if xor(isLeftSwingAtContact, isRightSwingAtContact)
            if isLeftSwingAtContact; firstModPaw = 2; else; firstModPaw = 3; end
        
        % if both in swing: pick the paw closest to the hurdle
        elseif (isLeftSwingAtContact && isRightSwingAtContact)
            if trialLocations(contactInd,1,2) > trialLocations(contactInd,1,3)
                firstModPaw = 2;
            else
                firstModPaw = 3;
            end
            
        % if both in stance: pick first paw to enter swing
        else
            if find(~isnan(trialModStepIds(:,2)), 1, 'first') < find(~isnan(trialModStepIds(:,3)), 1, 'first')
                firstModPaw = 2;
            else
                firstModPaw = 3;
            end
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
        kinData(j).locationsPix = locationsPawsPixels(trialInds,:,:);  % paw kinematics for trial in original pixel coordinates
        kinData(j).locationsTail = trialLocationsTail;  % tail kinematics for trial (frame X xyz X tailPosition)
        kinData(j).wiskContactPositions = wiskContactPositions(j);  % position relative to nose at which whiskers touch obstacle
        kinData(j).wiskContactTimes = wiskContactTimes(j);  % time at which whiskers touch obstacle
        kinData(j).contactInd = contactInd;  % frame index within trial at which whiskers contact obstacle
        kinData(j).obsPixPos = trialObsPixPos;  % positions of the obstacle in pixels coordinates, prior to shifting kinematics such that obstacle is 'stationary' // use this to 'undo' the unravelling of paw kinematics in other code as necessary

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
        if ~any(~isTrialAnalyzed)
            fprintf('WARNING! Unable to analyze trial(s): %i', j);
        else
            fprintf(' %i', j);
        end
        isTrialAnalyzed(j) = false;
        kinData(j).trialInds = find(frameTimeStamps>obsOnTimes(j) & frameTimeStamps<obsOffTimes(j)); % still useful to compute trialInds even when kinematics can't be analyzed...
    end
    
    kinData(j).isTrialAnalyzed = isTrialAnalyzed(j);
end
if any(~isTrialAnalyzed); fprintf('\n'); end
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
    controlVels = controlVels(2:2:end, :);  % take only second control step vels
    controlLengths = controlLengths(2:2:end, :);  % predict second control step length
    
    
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
    if s.showStepLengthFits
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

save(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'kinData.mat'), ...
    'kinData', 'stanceBins', 'models', '-v7.3')


