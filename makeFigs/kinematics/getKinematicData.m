function data = getKinematicData(sessions)

% settings
obsPos = -.0073;
speedTime = .02; % compute velocity over this interval
controlSteps = 2;
interpSmps = 100;
swingMin = .005; % % swings must be at least this long to be included in analysis (meters)

% initializations
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx']);
data = struct();
dataInd = 1;


% collect data for all trials
for i = 1:length(sessions)
    
    % report progress
    fprintf('%s: collecting data\n', sessions{i});
    
    % load session data
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'],...
            'obsPositions', 'obsTimes', 'obsPixPositions', 'frameTimeStamps', 'mToPixMapping', ...
            'obsOnTimes', 'obsOffTimes', 'nosePos', 'targetFs', 'wheelPositions', 'wheelTimes');
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));
    mToPixMapping = median(mToPixMapping,1);
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\tracking\locationsBotCorrected.mat'], 'locations')
    locations = locations.locationsCorrected;
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\tracking\stanceBins.mat'], 'stanceBins')
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\tracking\isExcluded.mat'], 'isExcluded')
    
    % get velocities for all trials in session
    sessionVels = getTrialSpeedsAtObsPos(obsPos, wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, speedTime, targetFs);
    
    % normalize y values
    locations(:,2,:) = locations(:,2,:) - nosePos(2); % subtract midline from all y values
    
    % get stance identities
    swingBins = ~stanceBins;
    swingIdentities = nan(size(swingBins));
    for j = 1:4
        
        % get start and end of swings
        swingStartInds = find([0; diff(swingBins(:,j))==1]');
        swingEndInds = find([diff(swingBins(:,j))==-1; 0]');
        
        % make sure find ind is a start and last ind is an end
        swingStartInds = swingStartInds(swingStartInds<swingEndInds(end));
        swingEndInds = swingEndInds(swingEndInds>swingStartInds(1));
        
        % remove swings that are too short
        validBins = locations(swingEndInds,1,j) - locations(swingStartInds,1,j) > swingMin;
        swingStartInds = swingStartInds(validBins);
        swingEndInds = swingEndInds(validBins);
        
        swingCount = 1;
        for k = 1:length(swingStartInds)
            swingIdentities(swingStartInds(k):swingEndInds(k),j) = swingCount;
            swingCount = swingCount + 1;
        end
    end
    
    
    % collect data for all trials within session
    for j = 1:length(obsOnTimes)-1
        
        % get trial bins, locations, and swingIdentities
        trialBins = frameTimeStamps>=obsOnTimes(j) & frameTimeStamps<=obsOffTimes(j) & ~isnan(obsPixPositions)';
        trialLocations = locations(trialBins,:,:);
        trialSwingIdentities = swingIdentities(trialBins,:);
        trialTimeStamps = frameTimeStamps(trialBins);
        trialObsPixPositions = obsPixPositions(trialBins);
        trialIsExcluded = isExcluded(trialBins);
        
        if any(~isnan(trialLocations(:))) && ~any(trialIsExcluded) % !!! this is a hack // should check that velocity criteria is met AND that the locations have in fact been analyzed for the session
        
            % get frame ind at which obs reaches obsPos
            obsPosTime = obsTimes(find(obsPositions>=obsPos & obsTimes>obsOnTimes(j), 1, 'first'));
            obsPosInd = knnsearch(trialTimeStamps, obsPosTime);
            
            % get trial swing identities and define control and modified steps
            controlStepIdentities = nan(size(trialSwingIdentities));
            modifiedStepIdentities = nan(size(trialSwingIdentities));

            for k = 1:4

                overObsInd = find(trialLocations(:,1,k)>trialObsPixPositions' & trialTimeStamps>obsOnTimes(j), 1, 'first');
                swingOverObsIdentity = trialSwingIdentities(overObsInd, k);
                firstModifiedIdentitiy = trialSwingIdentities(find(~isnan(trialSwingIdentities(:,k))' & 1:size(trialSwingIdentities,1)>=obsPosInd, 1, 'first'), k);

                modifiedBins = (trialSwingIdentities(:,k) >= firstModifiedIdentitiy) & (trialSwingIdentities(:,k) <= swingOverObsIdentity);
                controlBins = (trialSwingIdentities(:,k) >= (firstModifiedIdentitiy-controlSteps)) & (trialSwingIdentities(:,k) < firstModifiedIdentitiy);

                modifiedStepIdentities(:,k) = cumsum([0; diff(modifiedBins)==1]);
                modifiedStepIdentities(~modifiedBins,k) = nan;
                if ~any(~isnan(modifiedStepIdentities(:,k))); keyboard; end
                controlStepIdentities(:,k) = cumsum([0; diff(controlBins)==1]);
                controlStepIdentities(~controlBins,k) = nan;

            end
            
            % determine whether left and right forepaws are in swing at obsPos moment
            isLeftSwing = ~isnan(modifiedStepIdentities(obsPosInd,2));
            isRightSwing = ~isnan(modifiedStepIdentities(obsPosInd,3));
            oneSwingOneStance = xor(isLeftSwing, isRightSwing);
            
            % flip y values if the left fore is the swinging foot (thus making it the right paw)
            isFlipped = false;
            if oneSwingOneStance && isLeftSwing
                trialLocations = trialLocations(:,:,[4 3 2 1]);
                controlStepIdentities = controlStepIdentities(:,[4 3 2 1]);
                modifiedStepIdentities = modifiedStepIdentities(:,[4 3 2 1]);
                trialLocations(:,2,:) = -trialLocations(:,2,:);
                isFlipped = true;
            end
            
            % correct x locations (transform them s.t. obs is always at position 0 and positions move forward as though there were no wheel)
            trialLocations(:,1,:) = trialLocations(:,1,:) - trialObsPixPositions';           
            
            % convert to meters
            trialLocations = trialLocations / abs(mToPixMapping(1));
            
            % get stance distance from obs
            stanceDistance = trialLocations(obsPosInd,1,2); % left fore paw (2) is always the stance foot at this point after flipping y values above
            
            % get control step(s) length
            controlSwingLengths = nan(controlSteps,4);
            controlSwingDurations = nan(controlSteps,4);
            for k = 1:4
                for m = 1:controlSteps
                    stepBins = controlStepIdentities(:,k)==m;
                    stepXLocations = trialLocations(stepBins,1,k);
                    controlSwingLengths(m,k) = stepXLocations(end) - stepXLocations(1);
                    stepTimes = trialTimeStamps(stepBins);
                    controlSwingDurations(m,k) = stepTimes(end) - stepTimes(1);
                end
            end
            
            % get first modified step length for swing foot
            modifiedSwingLengths = nan(1,4);
            modifiedSwingDurations = nan(1,4);
            for k = 1:4
                stepBins = modifiedStepIdentities(:,k)==1;
                stepXLocations = trialLocations(stepBins,1,k);
                modifiedSwingLengths(k) = stepXLocations(end) - stepXLocations(1);
                stepTimes = trialTimeStamps(stepBins);
                modifiedSwingDurations(m,k) = stepTimes(end) - stepTimes(1);
            end
            
            % get interpolated control and modified step locations
            leftControlLocations = cell(1,4);
            leftModLocations = cell(1,4);
            modStepNum = nan(1,4);
            for k = 1:4
                
                % control
                stepNum = max(controlStepIdentities(:,k));
                pawControlLocations = nan(stepNum, 2, interpSmps);
                
                for m = 1:stepNum
                    stepInds = controlStepIdentities(:,k)==m;
                    stepX = trialLocations(stepInds,1,k);
                    stepY = trialLocations(stepInds,2,k);
                    xInterp = interp1(1:length(stepX), stepX, linspace(1,length(stepX),interpSmps));
                    yInterp = interp1(1:length(stepY), stepY, linspace(1,length(stepY),interpSmps));
                    pawControlLocations(m,:,:) = cat(1,xInterp,yInterp);
                end
                
                leftControlLocations{k} = pawControlLocations;
                
                % modified
                modStepNum(k) = max(modifiedStepIdentities(:,k));
                pawModifiedLocations = nan(modStepNum(k), 2, interpSmps);
                
                for m = 1:modStepNum(k)
                    stepInds = modifiedStepIdentities(:,k)==m;
                    stepX = trialLocations(stepInds,1,k);
                    stepY = trialLocations(stepInds,2,k);
                    xInterp = interp1(1:length(stepX), stepX, linspace(1,length(stepX),interpSmps));
                    yInterp = interp1(1:length(stepY), stepY, linspace(1,length(stepY),interpSmps));
                    pawModifiedLocations(m,:,:) = cat(1,xInterp,yInterp);
                end
                
                leftModLocations{k} = pawModifiedLocations;
            end



            % store results
            sessionInfoBin = find(strcmp(sessionInfo.session, sessions{i}));
            data(dataInd).mouse = sessionInfo.mouse{sessionInfoBin};
            data(dataInd).session = sessions{i};
            data(dataInd).vel = sessionVels(j);
            data(dataInd).obsPosInd = obsPosInd;
            data(dataInd).timeStamps = trialTimeStamps;
            data(dataInd).locations = trialLocations;
            data(dataInd).controlLocations = leftControlLocations;
            data(dataInd).modifiedLocations = leftModLocations;
            data(dataInd).controlStepIdentities = controlStepIdentities;
            data(dataInd).modifiedStepIdentities = modifiedStepIdentities;
            data(dataInd).modStepNum = modStepNum;
            data(dataInd).oneSwingOneStance = oneSwingOneStance;
            data(dataInd).stanceDistance = stanceDistance;
            data(dataInd).controlSwingLengths = controlSwingLengths;
            data(dataInd).modifiedSwingLengths = modifiedSwingLengths;
            data(dataInd).controlSwingDurations = controlSwingDurations;
            data(dataInd).modifiedSwingDurations = modifiedSwingDurations;
            data(dataInd).isFlipped = isFlipped;
            dataInd = dataInd + 1;
        end
    end
end

fprintf('--- done collecting data ---\n');