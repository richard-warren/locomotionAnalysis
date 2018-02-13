% function getObsTrajectories(sessions)

% temp
sessions = {'180122_001', '180122_002', '180122_003'};%, ...
%             '180123_001', '180123_002', '180123_003', ...
%             '180124_001', '180124_002', '180124_003'}; ...
%             '180125_001', '180125_002', '180125_003'};

% settings
obsPos = -.0073;
speedTime = .02; % compute velocity over this interval
controlSteps = 2;

% initializations
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx']);
data = struct();
dataInd = 1;

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
    
    % get velocities for all trials in session
    sessionVels = getTrialSpeedsAtObsPos(obsPos, wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, speedTime, targetFs);
    
    % normalize y values
    locations(:,2,:) = locations(:,2,:) - nosePos(2); % subtract midline from all y values
    
    % get stance identities
    swingBins = ~stanceBins;
    swingIdentities = nan(size(swingBins));
    for j = 1:4
        swingStartBins = [0; diff(swingBins(:,j))==1]';
        swingIdentities(:,j) = cumsum(swingStartBins);
    end
    swingIdentities(stanceBins | isnan(squeeze(locations(:,1,:)))) = nan;
    
    
    for j = 1:length(obsOnTimes)-1
        disp(j)
        
        % get trial bins, locations, and swingIdentities
        trialBins = frameTimeStamps>=obsOnTimes(j) & frameTimeStamps<=obsOffTimes(j) & ~isnan(obsPixPositions)';
        trialLocations = locations(trialBins,:,:);
        trialSwingIdentities = swingIdentities(trialBins,:);
        trialTimeStamps = frameTimeStamps(trialBins);
        trialObsPixPositions = obsPixPositions(trialBins);
        
        if any(~isnan(trialLocations(:))) % !!! this is a hack // should check that velocity criteria is met AND that the locations have in fact been analyzed for the session
        
            % get frame ind at which obs reaches obsPos
            obsPosTime = obsTimes(find(obsPositions>=obsPos & obsTimes>obsOnTimes(j), 1, 'first'));
            obsPosInd = knnsearch(trialTimeStamps, obsPosTime);

            % correct x locations (transform them s.t. obs is always at position 0 and positions move forward as though there were no wheel)
            trialLocations = trialLocations - trialObsPixPositions';           
            
            % convert to meters
            trialLocations = trialLocations / abs(mToPixMapping(1));
            keyboard
            
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
                controlStepIdentities(:,k) = cumsum([0; diff(controlBins)==1]);
                controlStepIdentities(~controlBins,k) = nan;

            end




            % store results
            sessionInfoBin = find(strcmp(sessionInfo.session, sessions{i}));
            data(dataInd).mouse = sessionInfo.mouse{sessionInfoBin};
            data(dataInd).session = sessions{i};
            data(dataInd).vel = sessionVels(j);
            data(dataInd).obsPosInd = obsPosInd;
            data(dataInd).locations = trialLocations;
            data(dataInd).controlStepIdentities = controlStepIdentities;
            data(dataInd).modifiedStepIdentities = modifiedStepIdentities;
            dataInd = dataInd + 1;
        end
    end
end

fprintf('--- done collecting data ---\n');


%% plot some thangs








