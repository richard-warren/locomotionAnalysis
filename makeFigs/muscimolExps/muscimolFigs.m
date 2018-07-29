%% load session info

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'muscimolNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);



%% get kinematic data
obsPos = -0.0087;
kinData = getKinematicData4(sessionInfo.session, obsPos);
data = kinData; save([getenv('OBSDATADIR') 'kinematicData\muscimolKinematicData.mat'], 'data');



%% get avoidance and speed data

speedAvoidanceData = getSpeedAndObsAvoidanceData(sessionInfo.session, false);
data = speedAvoidanceData; save([getenv('OBSDATADIR') 'kinematicData\muscimolSpeedAvoidanceData.mat'], 'data');

%% plot baseline kinematics

mice = unique(sessionInfo.mouse);

for j = 1:length(mice)
    
    mouseBins = strcmp(sessionInfo.mouse, mice{j});
    brainRegion = unique(sessionInfo.brainRegion(mouseBins));
    musBins = strcmp(sessionInfo.injectedSubstance, 'muscimol');
    vehBins = strcmp(sessionInfo.injectedSubstance, 'saline');
    sides = unique(sessionInfo.side(mouseBins & (musBins | vehBins)));
    
    for k = 1:length(sides)
        sideBins = strcmp(sessionInfo.side, sides{k});
        musSessions = sessionInfo.session(mouseBins & sideBins & musBins);
        vehSessions = sessionInfo.session(mouseBins & sideBins & vehBins);
        
        conditionBins = nan(1,length(kinData));
        conditionBins(ismember({kinData.session}, musSessions)) = 1;
        conditionBins(ismember({kinData.session}, vehSessions)) = 2;
        
        conditionLabels = {['muscimol - ' sides{k} ' ' brainRegion{1}], ['vehicle - ' sides{k} ' ' brainRegion{1}]};
        plotBaselineKinematics(kinData, conditionBins, conditionLabels, mice{j})
    end
end


%% plot dv averages

% settings
conditions = {'saline', 'muscimol'};
brainRegions = {'sen', 'mtc'};
dvs = {'success rate', ...
       'speed (m/s)', ...
       {'body angle towards contra', '(or right) side'}, ...
       'contra (or right) error rate', ...
       'ipsi (or left) error rate'};
dvYLims = [0 1; 0 .8; -15 15; 0 1; 0 1];
minTrial = 0;
% validTrials = ~[speedAvoidanceData.isLightOn] & [speedAvoidanceData.trialNum]>=minTrial;
validTrials = [speedAvoidanceData.trialNum]>=minTrial;


% initializations
figure('name', sprintf('min trial: %i', minTrial), 'Color', 'white', 'MenuBar', 'none', 'Position', [-1396 50 500 900])
dims = [length(dvs), length(brainRegions)]; % subplot dimensions

% loop over brain regions
for i = 1:length(brainRegions)
    
    brainRegionBins = strcmp(sessionInfo.brainRegion, brainRegions{i});
    mice = unique(sessionInfo.mouse(brainRegionBins));
    xJitters = linspace(-.1,.1,length(mice)); xJitters = xJitters-mean(xJitters); % jitters x position of scatter points
    colors = winter(length(mice));
    
    % containers for averages for each mouse for each condition across all sessions
    speeds = nan(length(mice), length(conditions)); % rows are mice, columns are conditions (saline, muscimol)
    successes = nan(length(mice), length(conditions));
    contraBodyAngles = nan(length(mice), length(conditions)); % angling towards contra side of body
    contraErrRates = nan(length(mice), length(conditions));
    ipsiErrRates = nan(length(mice), length(conditions));
    
    for j = 1:length(mice)
        for k = 1:length(conditions)
            
            conditionBins = brainRegionBins & ...
                            strcmp(sessionInfo.injectedSubstance, conditions{k}) & ...
                            strcmp(sessionInfo.mouse, mice{j});
            sessions = unique(sessionInfo.session(conditionBins));
            
            % containers for session averages for all session for a given mouse in a given condition
            sessionSuccesses = nan(1, length(sessions));
            sessionSpeeds = nan(1, length(sessions));
            sessionContraBodyAngles = nan(1, length(sessions));
            sessionContraErrRates = nan(1, length(sessions));
            sessionIpsiErrRates = nan(1, length(sessions));
            
            for m = 1:length(sessions)
                % get speed and success rate
                sessionBins = strcmp({speedAvoidanceData.session}, sessions{m}) & validTrials;
                sessionSuccesses(m) = mean([speedAvoidanceData(sessionBins).isObsAvoided]);
                sessionSpeeds(m) = mean([speedAvoidanceData(sessionBins).avgVel]);
                
                % get body angle
                load([getenv('OBSDATADIR') 'sessions\' sessions{m} '\runAnalyzed.mat'], ...
                'bodyAngles', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'arePawsTouchingObs');
                load([getenv('OBSDATADIR') 'sessions\' sessions{m} '\run.mat'], 'breaks');
                bodyAngles = getTrialBodyAngles(bodyAngles, obsOnTimes, obsOffTimes, frameTimeStamps, breaks);
                sessionContraBodyAngles(m) = nanmedian(bodyAngles);
                sideOfBrain = sessionInfo.side{strcmp(sessionInfo.session, sessions{m})};
                if strcmp(sideOfBrain, 'left'); sessionContraBodyAngles(m) = -sessionContraBodyAngles(m); end
                
                % get contra err rate
                leftErrors = 0;
                rightErrors = 0;
                for n = 1:length(obsOnTimes)
                    trialInds = frameTimeStamps>obsOnTimes(n) & frameTimeStamps<obsOffTimes(n);
                    leftErrors = leftErrors + any(any(arePawsTouchingObs(trialInds,[1 2])));
                    rightErrors = rightErrors + any(any(arePawsTouchingObs(trialInds,[3 4]),2)>0);
                end
                if strcmp(sideOfBrain, 'left')
                    sessionContraErrRates = rightErrors / length(obsOnTimes);
                    sessionIpsiErrRates = leftErrors / length(obsOnTimes);
                else
                    sessionContraErrRates = leftErrors / length(obsOnTimes);
                    sessionIpsiErrRates = rightErrors / length(obsOnTimes);
                end
            end
            
            successes(j,k) = nanmean(sessionSuccesses);
            speeds(j,k) = nanmean(sessionSpeeds);
            contraBodyAngles(j,k) = nanmean(sessionContraBodyAngles);
            contraErrRates(j,k) = nanmean(sessionContraErrRates);
            ipsiErrRates(j,k) = nanmean(sessionIpsiErrRates);
        end
        
        % plot mouse success
        subplot(dims(1), dims(2), i)
        line([1:length(conditions)] + xJitters(j), successes(j,:), 'color', [.5 .5 .5]); hold on
        scatter([1:length(conditions)] + xJitters(j), successes(j,:), 50, colors(j,:), 'filled');
        
        % plot mouse speed
        subplot(dims(1), dims(2), i+1*length(brainRegions))
        line([1:length(conditions)] + xJitters(j), speeds(j,:), 'color', [.5 .5 .5]); hold on
        scatter([1:length(conditions)] + xJitters(j), speeds(j,:), 50, colors(j,:), 'filled');
        
        % plot contra body angle
        subplot(dims(1), dims(2), i+2*length(brainRegions))
        line([1:length(conditions)] + xJitters(j), contraBodyAngles(j,:), 'color', [.5 .5 .5]); hold on
        scatter([1:length(conditions)] + xJitters(j), contraBodyAngles(j,:), 50, colors(j,:), 'filled');
        
        % plot contra err rates
        subplot(dims(1), dims(2), i+3*length(brainRegions))
        line([1:length(conditions)] + xJitters(j), contraErrRates(j,:), 'color', [.5 .5 .5]); hold on
        scatter([1:length(conditions)] + xJitters(j), contraErrRates(j,:), 50, colors(j,:), 'filled');
        
        % plot ipsi err rates
        subplot(dims(1), dims(2), i+4*length(brainRegions))
        line([1:length(conditions)] + xJitters(j), ipsiErrRates(j,:), 'color', [.5 .5 .5]); hold on
        scatter([1:length(conditions)] + xJitters(j), ipsiErrRates(j,:), 50, colors(j,:), 'filled');
    end
    
    % plot condition means
    for k = 1:length(conditions)
        % success
        subplot(dims(1), dims(2), i)
        avg = nanmean(successes(:,k));
        line([k-.1 k+.1], [avg avg], 'linewidth', 3, 'color', 'black')
        
        % speed
        subplot(dims(1), dims(2), i+1*length(brainRegions))
        avg = nanmean(speeds(:,k));
        line([k-.1 k+.1], [avg avg], 'linewidth', 3, 'color', 'black')
                
        % contra body angles
        subplot(dims(1), dims(2), i+2*length(brainRegions))
        avg = nanmean(contraBodyAngles(:,k));
        line([k-.1 k+.1], [avg avg], 'linewidth', 3, 'color', 'black')
        
        % contra err rates
        subplot(dims(1), dims(2), i+3*length(brainRegions))
        avg = nanmean(contraErrRates(:,k));
        line([k-.1 k+.1], [avg avg], 'linewidth', 3, 'color', 'black')
        
        % ipsi err rates
        subplot(dims(1), dims(2), i+4*length(brainRegions))
        avg = nanmean(ipsiErrRates(:,k));
        line([k-.1 k+.1], [avg avg], 'linewidth', 3, 'color', 'black')
    end
    
    % add mouse labels
    xLims = get(gca, 'xlim');
    xs = linspace(xLims(1)*1.2, xLims(2)*.8, length(mice));
    for j = 1:length(mice)
        text(xs(j), dvYLims(end,1)+(dvYLims(end,2)-dvYLims(end,1))*.2, mice{j}, 'Color', colors(j,:));
    end

end


% pimp figs
ind = 1;
for i = 1:length(dvs)
    for j = 1:length(brainRegions)
        subplot(dims(1), dims(2), ind);
        
        set(gca, 'xlim', [0.75 length(conditions)+0.25], 'xtick', 1:length(conditions), 'XTickLabel', conditions, ...
            'YLim', dvYLims(i,:));
        if i==1; title(brainRegions{j}); end
        if j==1; ylabel(dvs{i}); end
        ind = ind+1;
    end
end


saveas(gcf, [getenv('OBSDATADIR') 'figures\muscimolBarPlots.png']);
savefig([getenv('OBSDATADIR') 'figures\muscimolBarPlots.fig'])



%% plot dvs over time within sessions

% settings
trialBinSize = 40;
conditions = {'saline', 'muscimol'};
brainRegions = {'sen', 'mtc'};
dvs = {'success rate', ...
       'speed (m/s)', ...
       {'body angle towards contra', '(or right) side'}};
dvYLims = [0 1; 0 .8; -15 15];
xMax = 240;


% initializations
figure('Color', 'white', 'MenuBar', 'none', 'Position', [-1396 200 800 600])
dims = [length(dvs), length(brainRegions)]; % subplot dimensions
colors = winter(length(conditions));

% loop over brain regions
for i = 1:length(brainRegions)
    
    % get max trial num
    brainRegionBins = ismember({speedAvoidanceData.session}, sessions);
    maxTrial = max([speedAvoidanceData(brainRegionBins).trialNum]);
    maxTrial = min(maxTrial, xMax);
    maxTrial = floor(maxTrial/trialBinSize) * trialBinSize;
    trialBinEdges = 0:trialBinSize:maxTrial;
    xInds = trialBinEdges(2:end); % the x positions will be centered within the bins
    
    
    for j = 1:length(conditions)
        
        sessions = unique(sessionInfo.session(strcmp(sessionInfo.brainRegion, brainRegions{i}) ...
                                            & strcmp(sessionInfo.injectedSubstance, conditions{j})));
        
        % containers for bins averages across all sessions
        successes = nan(length(sessions), length(trialBinEdges)-1);
        speeds = nan(length(sessions), length(trialBinEdges)-1);
        contraBodyAngles = nan(length(sessions), length(trialBinEdges)-1);
        
        for k = 1:length(sessions)
            
            load([getenv('OBSDATADIR') 'sessions\' sessions{k} '\runAnalyzed.mat'], ...
                'bodyAngles', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps');
            load([getenv('OBSDATADIR') 'sessions\' sessions{k} '\run.mat'], 'breaks');
            bodyAngles = getTrialBodyAngles(bodyAngles, obsOnTimes, obsOffTimes, frameTimeStamps, breaks);
            sideOfBrain = sessionInfo.side{strcmp(sessionInfo.session, sessions{k})};
            if strcmp(sideOfBrain, 'left'); bodyAngles = -bodyAngles; end
            
            
            for m = 1:length(trialBinEdges)-1
            
                trialBinBins = strcmp({speedAvoidanceData.session}, sessions{k}) & ...
                               [speedAvoidanceData.trialNum]>trialBinEdges(m) & ...
                               [speedAvoidanceData.trialNum]<trialBinEdges(m+1);

                % get speed, success rate, and body angle
                successes(k,m) = mean([speedAvoidanceData(trialBinBins).isObsAvoided]);
                speeds(k,m) = mean([speedAvoidanceData(trialBinBins).avgVel]);
                bodyAngleBins = 1:length(bodyAngles) > trialBinEdges(m) & 1:length(bodyAngles) < trialBinEdges(m+1);
                contraBodyAngles(k,m) = nanmean(bodyAngles(bodyAngleBins));
            end
            
            
            % plot session data
            % plot session success
            subplot(dims(1), dims(2), i)
            line(xInds, successes(k,:), 'color', colors(j,:)); hold on

            % plot session speed
            subplot(dims(1), dims(2), i+1*length(brainRegions))
            line(xInds, speeds(k,:), 'color', colors(j,:)); hold on

            % plot session body angle
            subplot(dims(1), dims(2), i+2*length(brainRegions))
            line(xInds, contraBodyAngles(k,:), 'color', colors(j,:)); hold on
        end
        
        
        % plot condition means
        
        % mean successes
        subplot(dims(1), dims(2), i)
        line(xInds, nanmean(successes,1), 'color', colors(j,:), 'linewidth', 4); hold on

        % mean speeds
        subplot(dims(1), dims(2), i+1*length(brainRegions))
        line(xInds, nanmean(speeds,1), 'color', colors(j,:), 'linewidth', 4); hold on
        
        % mean session body angle
        subplot(dims(1), dims(2), i+2*length(brainRegions))
        line(xInds, nanmean(contraBodyAngles,1), 'color', colors(j,:), 'linewidth', 4); hold on
    end
    
    % add condition labels
    xs = linspace(trialBinEdges(2)*1.2, trialBinEdges(end)*.8, length(conditions));
    for j = 1:length(conditions)
        text(xs(j), dvYLims(end,1)+(dvYLims(end,2)-dvYLims(end,1))*.2, conditions{j}, 'Color', colors(j,:));
    end
end
    


% pimp figs
ind = 1;
for i = 1:length(dvs)
    for j = 1:length(brainRegions)
        subplot(dims(1), dims(2), ind);
        
        set(gca, 'xlim', [trialBinEdges(2) trialBinEdges(end)], 'xtick', trialBinEdges(2:end), 'YLim', dvYLims(i,:));
        if i==1; title(brainRegions{j}); end
        if j==1; ylabel(dvs{i}); end
        ind = ind+1;
    end
end


saveas(gcf, [getenv('OBSDATADIR') 'figures\muscimolSessionProgressPlots.png']);
savefig([getenv('OBSDATADIR') 'figures\muscimolSessionProgressPlotsPlots.fig'])


%% plot dvs across sessions

% settings
conditions = {'saline', 'muscimol'};
brainRegions = {'sen', 'mtc'};
dvs = {'success rate', ...
       'speed (m/s)', ...
       {'body angle towards contra', '(or right) side'}};
dvYLims = [.4 1; .1 .6; -15 15];
minTrial = 0;


% initializations
figure('name', sprintf('min trial: %i', minTrial), 'Color', 'white', 'MenuBar', 'none', 'Position', [-1396 200 800 600])
dims = [length(dvs), length(brainRegions)]; % subplot dimensions

% loop over brain regions
for i = 1:length(brainRegions)
    
    brainRegionBins = strcmp(sessionInfo.brainRegion, brainRegions{i});
    mice = unique(sessionInfo.mouse(brainRegionBins));
    colors = winter(length(mice));
    
    
    for j = 1:length(mice)
        conditionBins = brainRegionBins & ...
                        strcmp(sessionInfo.mouse, mice{j});
        sessions = unique(sessionInfo.session(conditionBins));

        % containers for session averages for all session for a given mouse in a given condition
        sessionSuccesses = nan(1, length(sessions));
        sessionSpeeds = nan(1, length(sessions));
        sessionContraBodyAngles = nan(1, length(sessions));
        musBins = false(1, length(sessions));

        for m = 1:length(sessions)
            % get speed and success rate
            sessionBins = strcmp({speedAvoidanceData.session}, sessions{m}) & ...
                          [speedAvoidanceData.trialNum]>=minTrial;
            sessionSuccesses(m) = mean([speedAvoidanceData(sessionBins).isObsAvoided]);
            sessionSpeeds(m) = mean([speedAvoidanceData(sessionBins).avgVel]);
            
            sideOfBrain = sessionInfo.side{strcmp(sessionInfo.session, sessions{m})};

            % get body angle
            load([getenv('OBSDATADIR') 'sessions\' sessions{m} '\runAnalyzed.mat'], ...
                'bodyAngles', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps');
            load([getenv('OBSDATADIR') 'sessions\' sessions{m} '\run.mat'], 'breaks');
            bodyAngles = getTrialBodyAngles(bodyAngles, obsOnTimes, obsOffTimes, frameTimeStamps, breaks);
            sessionContraBodyAngles(m) = nanmedian(bodyAngles);
            if strcmp(sideOfBrain, 'left'); sessionContraBodyAngles(m) = -sessionContraBodyAngles(m); end
            
            % get condition whether or not muscimol
            sessionBin = strcmp(sessionInfo.session, sessions{m});
            musBins(m) = strcmp(sessionInfo.injectedSubstance(sessionBin), 'muscimol');
        end
        
        
        % plot mouse data
        xInds = 1:length(sessions);
        
        % plot mouse success
        subplot(dims(1), dims(2), i)
        line(1:length(sessions), sessionSuccesses, 'color', [.5 .5 .5]); hold on
        scatter(xInds(musBins), sessionSuccesses(musBins), 50, colors(j,:), 'filled'); % muscimol
        scatter(xInds(~musBins), sessionSuccesses(~musBins), 50, colors(j,:)); % saline
        
        % plot mouse speed
        subplot(dims(1), dims(2), i+1*length(brainRegions))
        line(1:length(sessions), sessionSpeeds, 'color', [.5 .5 .5]); hold on
        scatter(xInds(musBins), sessionSpeeds(musBins), 50, colors(j,:), 'filled'); % muscimol
        scatter(xInds(~musBins), sessionSpeeds(~musBins), 50, colors(j,:)); % saline
        
        % plot contra body angle
        subplot(dims(1), dims(2), i+2*length(brainRegions))
        line(1:length(sessions), sessionContraBodyAngles, 'color', [.5 .5 .5]); hold on
        scatter(xInds(musBins), sessionContraBodyAngles(musBins), 50, colors(j,:), 'filled'); % muscimol
        scatter(xInds(~musBins), sessionContraBodyAngles(~musBins), 50, colors(j,:)); % saline;
    end
    
    % add mouse labels
    xLims = get(gca, 'xlim');
    xs = linspace(1, xLims(2)*.8, length(mice));
    for j = 1:length(mice)
        text(xs(j), dvYLims(end,1)+(dvYLims(end,2)-dvYLims(end,1))*.2, mice{j}, 'Color', colors(j,:));
    end
end
    


% pimp figs
ind = 1;
for i = 1:length(dvs)
    for j = 1:length(brainRegions)
        subplot(dims(1), dims(2), ind);
        
        xLims = get(gca, 'xlim');
        set(gca, 'xlim', [0.5 xLims(2)+.5], 'xtick', 1:xLims(2), 'YLim', dvYLims(i,:));
        if i==1; title(brainRegions{j}); end
        if j==1; ylabel(dvs{i}); end
        ind = ind+1;
    end
end

saveas(gcf, [getenv('OBSDATADIR') 'figures\muscimolWithinSessionProgressPlots.png']);
savefig([getenv('OBSDATADIR') 'figures\muscimolWithinSessionProgressPlotsPlots.fig'])

%% plot kinematics



% incorporate condition information into kinData struct
for i = 1:length(kinData)
    sessionInfoBin = strcmp(sessionInfo.session, kinData(i).session);
    brainRegion = sessionInfo.brainRegion{sessionInfoBin};
    injectedSubstance = sessionInfo.injectedSubstance{sessionInfoBin};
    if strcmp(sessionInfo.side{sessionInfoBin}, 'left')
        contraLimb = 3;
    elseif strcmp(sessionInfo.side{sessionInfoBin}, 'right')
        contraLimb = 2;
    else
        contraLimb = nan;
    end
    
    kinData(i).brainRegion = brainRegion;
    kinData(i).injectedSubstance = injectedSubstance;
    kinData(i).contraLimb = contraLimb;
end

% get trial bins

mice = unique(sessionInfo.mouse);

firstPawContraBins = [kinData.firstPawOver]==[kinData.contraLimb];
firstPawIpsiBins = [kinData.firstPawOver]~=[kinData.contraLimb];
musBins = strcmp({kinData.injectedSubstance}, 'muscimol');
vehBins = strcmp({kinData.injectedSubstance}, 'saline');
lightOnBins = [kinData.isLightOn];
mtcBins = strcmp({kinData.brainRegion}, 'mtc');
senBins = strcmp({kinData.brainRegion}, 'sen');
% validBins = ~[kinData.isWheelBreak] & [kinData.isObsAvoided] & ismember({kinData.mouse}, mice);
validBins = ~[kinData.isWheelBreak] & ~[kinData.isObsAvoided] & ismember({kinData.mouse}, mice);

% mtc, light on
bins = zeros(1,length(kinData));
figBins = mtcBins & validBins;
bins(firstPawIpsiBins & musBins & figBins) = 1; % ipsi
bins(firstPawContraBins & musBins & figBins) = 2; % contra
bins(vehBins & figBins) = 3; % vehicle
binLabels = {'ipsi', 'contra', 'vehicle ipsi'};
plotObsHeightTrajectories(kinData, bins, binLabels, 'motor cortex lesions, light on')

% mtc, light off
bins = zeros(1,length(kinData));
figBins = mtcBins & ~lightOnBins & validBins;
bins(firstPawIpsiBins & musBins & figBins) = 1; % ipsi
bins(firstPawContraBins & musBins & figBins) = 2; % contra
bins(vehBins & figBins) = 3; % vehicle
binLabels = {'ipsi', 'contra', 'vehicle'};
plotObsHeightTrajectories(kinData, bins, binLabels, 'motor cortex lesions, light off')

% % sen, light on
% bins = zeros(1,length(kinData));
% figBins = senBins & lightOnBins & validBins;
% bins(musBins & figBins) = 1; % mus
% bins(vehBins & figBins) = 2; % vehicle
% binLabels = {'muscimol', 'vehicle'};
% plotObsHeightTrajectories(kinData, bins, binLabels, 'sensory cortex lesions, light on')
% 
% % sen, light off
% bins = zeros(1,length(kinData));
% figBins = senBins & ~lightOnBins & validBins;
% bins(musBins & figBins) = 1; % mus
% bins(vehBins & figBins) = 2; % vehicle
% binLabels = {'muscimol', 'vehicle'};
% plotObsHeightTrajectories(kinData, bins, binLabels, 'sensory cortex lesions, light off')















