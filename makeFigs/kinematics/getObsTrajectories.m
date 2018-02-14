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
    
    
    % collect data for all trials within session
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
            
            % get stance distance from obs
            stanceDistance = trialLocations(obsPosInd,1,2); % left fore paw (2) is always the stance foot at this point after flipping y values above
            
            % get control step(s) length
            controlSwingLengths = nan(controlSteps,4);
            for k = 1:4
                for m = 1:controlSteps
                    stepXLocations = trialLocations(controlStepIdentities(:,k)==m,1,k);
                    controlSwingLengths(m,k) = stepXLocations(end) - stepXLocations(1);
                end
            end
            
            % get first modified step length for swing foot
            modifiedSwingLengths = nan(1,4);
            for k = 1:4
                stepXLocations = trialLocations(modifiedStepIdentities(:,k)==1,1,k);
                modifiedSwingLengths(k) = stepXLocations(end) - stepXLocations(1);
            end



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
            data(dataInd).controlSwingLengths = controlSwingLengths;
            data(dataInd).modifiedSwingLengths = modifiedSwingLengths;
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

% initializations
dataNew = data([data.oneSwingOneStance]);

% get speed and phase bins
phaseBinEdges = prctile([dataNew.stanceDistance], linspace(0,100,phaseBinNum+1));
phaseBins = discretize([dataNew.stanceDistance], phaseBinEdges);
speedBinEdges = prctile([dataNew.vel], linspace(0,100,speedBinNum+1));
speedBins = discretize([dataNew.vel], speedBinEdges);

% create speed and phase labels
phaseLabels = cell(1,phaseBinNum);
for i = 1:phaseBinNum
    phaseLabels{i} = sprintf('%.3f', mean([dataNew(phaseBins==i).stanceDistance]));
end

speedLabels = cell(1,speedBinNum);
for i = 1:speedBinNum
    speedLabels{i} = sprintf('%.3f', mean([dataNew(speedBins==i).vel]));
end


%% sperm plots

% settings
xLims = [-.1 .1];
tracesPerPlot = 20;

close all; figure; pimpFig

for g = 1:speedBinNum
    for h = 1:phaseBinNum

        axInd = sub2ind([phaseBinNum speedBinNum], h, g);
        ax = subaxis(speedBinNum, phaseBinNum ,axInd, ...
            'spacing', .01, 'padding', .01, 'margin', .01);
        
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
                    plot(x, y, 'color', colors(k,:)); hold on

                    % scatter dots at start of each swing
                    scatter(x(end), y(end), 100, colors(k,:), 'filled'); hold on

                    % scatter position of swing foot at obsPos
                    if j==3
                        scatter(dataNew(i).locations(dataNew(i).obsPosInd,1,j), dataNew(i).locations(dataNew(i).obsPosInd,2,j), ...
                            100, [0 0 0], 'x'); hold on
                    end
                end
            end
        end

        % set appearance
        set(gca, 'dataaspectratio', [1 1 1], 'xlim', xLims, 'box', 'off', 'tickdir', 'out', 'ydir', 'reverse', ...
            'xtick', [], 'ytick', []);
        line([0 0], get(gca,'ylim'), 'color', [0 0 0], 'linewidth', 3)
        if g==speedBinNum; xlabel(['stance foot distance (m): ' phaseLabels{h}]); end
        if h==1; ylabel(['speed (m/s): ' speedLabels{g}]); end
    end
end



%% histograms

% settings
xLims = [.02 .12];
yLims = [0 .4];
binWidth = .005;

modifiedSwingLengths = {dataNew.modifiedSwingLengths}; modifiedSwingLengths = cat(1, modifiedSwingLengths{:});
controlSwingLengths = {dataNew.controlSwingLengths}; controlSwingLengths = cat(1, controlSwingLengths{:});

figure; pimpFig

for g = 1:speedBinNum
    for h = 1:phaseBinNum
        
        axInd = sub2ind([phaseBinNum speedBinNum], h, g);
        ax = subaxis(speedBinNum, phaseBinNum , axInd);%, ...
%             'spacing', .01, 'padding', 0, 'margin', 0);
        binInds = find(speedBins==g & phaseBins==h);
        
        h1 = histogram(modifiedSwingLengths(binInds,3)); hold on;
        h2 = histogram(controlSwingLengths(binInds,3));
        
        for i = {h1,h2}
            set(i{1}, 'normalization', 'probability', 'binwidth', binWidth);
        end
        
        % set apearance
        set(ax, 'box', 'off', 'xlim', xLims, 'ylim', yLims, 'tickdir', 'out')
        if g==speedBinNum; xlabel(['stance foot distance (m): ' phaseLabels{h}]); end
        if g<speedBinNum; set(ax, 'xticklabel', []); end
        if h==1; ylabel(['speed (m/s): ' speedLabels{g}]); end
        if h>1; set(ax, 'ytick', [], 'ylabel', []); ax.YAxis.Visible = 'off'; end
    end
end

legend('modified swing lengths', 'control swing lengths')

%% heat map

% settings
distanceLims = [-.06 -.02];
velLims = [.2 .7];
smps = 100;

% initializations
modifiedSwingLengths = {dataNew.modifiedSwingLengths}; modifiedSwingLengths = cat(1, modifiedSwingLengths{:});

controlSwingLengths = cellfun(@(x) mean(x,1), {dataNew.controlSwingLengths}, 'uniformoutput', 0);
controlSwingLengths = cat(1, controlSwingLengths{:});
deltaLength = abs(modifiedSwingLengths(:,3) - controlSwingLengths(:,3));

validInds = deltaLength<.08;

[xq, yq] = meshgrid(linspace(distanceLims(1),distanceLims(2),smps), linspace(velLims(1),velLims(2),smps));
heatMap = griddata([dataNew(validInds).stanceDistance], [dataNew(validInds).vel], deltaLength(validInds), xq, yq);

figure('color', [1 1 1]);
imagesc('xdata', xq(1,:), 'ydata', yq(:,1), 'cdata', heatMap);
set(gca, 'xlim', distanceLims, 'ylim', velLims)
xlabel('stance paw distance (m)');
ylabel('speed (m/s)');


%% scatter

% settings


% initializations
modifiedSwingLengths = {dataNew.modifiedSwingLengths}; modifiedSwingLengths = cat(1, modifiedSwingLengths{:});

controlSwingLengths = cellfun(@(x) mean(x,1), {dataNew.controlSwingLengths}, 'uniformoutput', 0);
controlSwingLengths = cat(1, controlSwingLengths{:});
deltaLength = abs(modifiedSwingLengths(:,3) - controlSwingLengths(:,3));

validInds = deltaLength<.08;




figure('color', [1 1 1]);
imagesc('xdata', xq(1,:), 'ydata', yq(:,1), 'cdata', heatMap);
set(gca, 'xlim', distanceLims, 'ylim', velLims)
xlabel('stance paw distance (m)');
ylabel('speed (m/s)');
















