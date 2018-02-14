% function getObsTrajectories(sessions)

% temp
sessions = {'180122_001', '180122_002', '180122_003' ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003'};
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
            
            % get stance distance
            stanceDistance = trialLocations(obsPosInd,1,2); % left fore paw (2) is always the stance foot at this point after flipping y values above



            % store results
            sessionInfoBin = find(strcmp(sessionInfo.session, sessions{i}));
            data(dataInd).mouse = sessionInfo.mouse{sessionInfoBin};
            data(dataInd).session = sessions{i};
            data(dataInd).vel = sessionVels(j);
            data(dataInd).obsPosInd = obsPosInd;
            data(dataInd).timeStamps = trialTimeStamps;
            data(dataInd).locations = trialLocations;
            data(dataInd).controlStepIdentities = controlStepIdentities;
            data(dataInd).modifiedStepIdentities = modifiedStepIdentities;
            data(dataInd).oneSwingOneStance = oneSwingOneStance;
            data(dataInd).stanceDistance = stanceDistance;
            data(dataInd).isFlipped = isFlipped;
            dataInd = dataInd + 1;
        end
    end
end

fprintf('--- done collecting data ---\n');


%% plot some thangs

% settings
phaseBinNum = 3;
speedBinNum = 3;
tracesPerPlot = 20;
yLim = [-.1 .1];

% initializations
dataNew = data([data.oneSwingOneStance]);
phaseBinEdges = prctile([dataNew.stanceDistance], linspace(0,100,phaseBinNum+1));
phaseBins = discretize([dataNew.stanceDistance], phaseBinEdges);
speedBinEdges = prctile([dataNew.vel], linspace(0,100,speedBinNum+1));
speedBins = discretize([dataNew.vel], speedBinEdges);

close all; figure; pimpFig




for g = 1:speedBinNum
    for h = 1:phaseBinNum

        subaxis(speedBinNum, phaseBinNum ,sub2ind([speedBinNum phaseBinNum], g, h), ...
            'spacing', .01, 'padding', 0, 'margin', 0);
        
        binInds = find(speedBins==g & phaseBins==h);
        plotTraces = min(tracesPerPlot, length(binInds));
        dataInds = randperm(length(binInds), tracesPerPlot);
        dataInds = binInds(dataInds);

        for i = dataInds
            for j = 2:3

                realInds = ~isnan(dataNew(i).modifiedStepIdentities(:,j));
                steps = unique(dataNew(i).modifiedStepIdentities(realInds,j));
                colors = winter(length(steps));

                for k = steps'

                    % plot x and y trajectories
                    trialInds = dataNew(i).modifiedStepIdentities(:,j)==k;
                    x = dataNew(i).locations(trialInds,1,j);
                    y = dataNew(i).locations(trialInds,2,j);
                    plot(y, x, 'color', colors(k,:)); hold on

                    % scatter dots at start of each swing
                    scatter(y(end), x(end), 100, colors(k,:), 'filled'); hold on

                    % scatter position of swing foot at obsPos
                    if j==3
                        scatter(dataNew(i).locations(dataNew(i).obsPosInd,2,j), dataNew(i).locations(dataNew(i).obsPosInd,1,j), ...
                            100, [0 0 0], 'x'); hold on
                    end
                end
            end
        end

        % pimp figs
        set(gca, 'dataaspectratio', [1 1 1], 'ylim', yLim, 'box', 'off', 'tickdir', 'out');
            line(get(gca,'xlim'), [0 0], 'color', [0 0 0], 'linewidth', 3)'
            ax = gca; ax.YAxis.Visible = 'off';
        if h>1
            ax = gca;
            ax.XAxis.Visible = 'off';
            set(gca, 'xticklabel', [], 'yticklabel', [])
        end
    end
end





