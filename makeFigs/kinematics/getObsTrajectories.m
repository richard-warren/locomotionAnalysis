% function getObsTrajectories(sessions)

% temp
sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003'};
%             '180125_001', '180125_002', '180125_003'};

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
%     swingIdentities(stanceBins | isnan(squeeze(locations(:,1,:)))) = nan;
    
    
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
        
        bins = find(speedBins==g & phaseBins==h);
        plotTraces = min(tracesPerPlot, length(bins));
        dataInds = randperm(length(bins), tracesPerPlot);
        dataInds = bins(dataInds);

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
colors = winter(2);
controlColor = [.65 .65 .65];

% initializations
numModSteps = reshape([dataNew.modStepNum],4,length(dataNew))';
modifiedSwingLengths = {dataNew.modifiedSwingLengths}; modifiedSwingLengths = cat(1, modifiedSwingLengths{:});
controlSwingLengths = {dataNew.controlSwingLengths}; controlSwingLengths = cat(1, controlSwingLengths{:});

figure; pimpFig

for g = 1:speedBinNum
    for h = 1:phaseBinNum
        
        axInd = sub2ind([phaseBinNum speedBinNum], h, g);
        ax = subaxis(speedBinNum, phaseBinNum , axInd);
        bins = (speedBins==g & phaseBins==h)';
        oneStepBins = bins & numModSteps(:,3)==1;
        twoStepBins = bins & numModSteps(:,3)==2;
        
        if any(oneStepBins); h1 = histogram(modifiedSwingLengths(oneStepBins,3), 'facecolor', colors(2,:)); hold on; end
        h2 = histogram(modifiedSwingLengths(twoStepBins,3), 'facecolor', colors(1,:)); hold on;
        h3 = histogram(controlSwingLengths(bins,3), 'facecolor', controlColor); hold on;
        
        for i = {h1,h2,h3}
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

validBins = deltaLength<.08;

[xq, yq] = meshgrid(linspace(distanceLims(1),distanceLims(2),smps), linspace(velLims(1),velLims(2),smps));
heatMap = griddata([dataNew(validBins).stanceDistance], [dataNew(validBins).vel], deltaLength(validBins), xq, yq);

figure('color', [1 1 1]);
imagesc('xdata', xq(1,:), 'ydata', yq(:,1), 'cdata', heatMap);
set(gca, 'xlim', distanceLims, 'ylim', velLims)
xlabel('stance paw distance (m)');
ylabel('speed (m/s)');


%% average trajectories


% settings
controlColor = [.65 .65 .65];
xLims = [-.1 .05];
linWid = 4;
colors = winter(2);

% initializations
close all; figure; pimpFig;
numModSteps = reshape([dataNew.modStepNum],4,length(dataNew))';


for g = 1:speedBinNum
    for h = 1:phaseBinNum
        
        % get subplot bins
        axInd = sub2ind([phaseBinNum speedBinNum], h, g);
        ax = subaxis(speedBinNum, phaseBinNum , axInd, ...
            'spacing', .01, 'padding', .01, 'margin', .01);
        bins = (speedBins==g & phaseBins==h)';

        % get subplot bins for different conditions
        controlBins = bins;
        leftModBins = bins & numModSteps(:,2)==1;
        rightModOneStepBins = bins & (numModSteps(:,3)==1);
        rightModTwoStepBins = bins & (numModSteps(:,3)==2);
        oneTwoRatio = sum(rightModOneStepBins) / sum(rightModTwoStepBins); % ratio of trials in which swing foot takes one large step to those in which an additional step is taken
        oneTwoRatio = oneTwoRatio * 2 - 1; % scale from -1 to 1

        % get left and right control locations
        controlLocations = {dataNew(controlBins).controlLocations};
        leftControlLocations = cellfun(@(x) x{2}, controlLocations, 'uniformoutput', 0);
        leftControlLocations = cat(1,leftControlLocations{:});
        rightControlLocations = cellfun(@(x) x{3}, controlLocations, 'uniformoutput', 0);
        rightControlLocations = cat(1,rightControlLocations{:});

        % get left modified locations
        leftModLocations = {dataNew(leftModBins).modifiedLocations};
        leftModLocations = cellfun(@(x) x{2}, leftModLocations, 'uniformoutput', 0);
        leftModLocations = cat(1,leftModLocations{:});

        % get right modified (one step) locations
        rightModOneStepLocations = {dataNew(rightModOneStepBins).modifiedLocations};
        rightModOneStepLocations = cellfun(@(x) x{3}, rightModOneStepLocations, 'uniformoutput', 0);
        rightModOneStepLocations = cat(1,rightModOneStepLocations{:});

        % get right modified (two step) locations
        rightModTwoStepLocations = {dataNew(rightModTwoStepBins).modifiedLocations};
        rightModTwoStepLocations = cellfun(@(x) x{3}(1,:,:), rightModTwoStepLocations, 'uniformoutput', 0);
        rightModTwoStepLocations = cat(1,rightModTwoStepLocations{:});


        % plot control left
        x = squeeze(leftControlLocations(:,1,:));
        x = x - mean(x(:,1)) + mean(squeeze(leftModLocations(:,1,1)));
        y = squeeze(leftControlLocations(:,2,:));
        plot(mean(x,1), mean(y,1), 'color', controlColor, 'linewidth', linWid); hold on;

        % plot control right
        x = squeeze(rightControlLocations(:,1,:));
        x = x - mean(x(:,1)) + mean(squeeze(rightModTwoStepLocations(:,1,1))); % !!! make average of one and two step?
        y = squeeze(rightControlLocations(:,2,:));
        plot(mean(x,1), mean(y,1), 'color', controlColor, 'linewidth', linWid); hold on;
        
        % plot mod left
        x = squeeze(leftModLocations(:,1,:));
        y = squeeze(leftModLocations(:,2,:));
        plot(mean(x,1), mean(y,1), 'color', colors(1,:), 'linewidth', linWid); hold on;

        % plot mod right, one step
        if ~isempty(rightModOneStepLocations)
            x = squeeze(rightModOneStepLocations(:,1,:));
            y = squeeze(rightModOneStepLocations(:,2,:));
            plot(mean(x,1), mean(y,1), 'color', colors(2,:), 'linewidth', linWid + oneTwoRatio*linWid); hold on;
        end

        % plot mod right, two step
        x = squeeze(rightModTwoStepLocations(:,1,:));
        y = squeeze(rightModTwoStepLocations(:,2,:));
        plot(mean(x,1), mean(y,1), 'color', colors(1,:), 'linewidth', linWid + -oneTwoRatio*linWid); hold on;



        % set appearance
        set(gca, 'dataaspectratio', [1 1 1], 'xlim', xLims, 'box', 'off', 'tickdir', 'out', 'ydir', 'reverse', ...
            'xtick', [], 'ytick', []);
        line([0 0], get(gca,'ylim'), 'color', [0 0 0], 'linewidth', 3)
        if g==speedBinNum; xlabel(['stance foot distance (m): ' phaseLabels{h}]); end
        if h==1; ylabel(['speed (m/s): ' speedLabels{g}]); end
    end
end













