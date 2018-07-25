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


%% plot dvs averaged across sessions

% settings
conditions = {'saline', 'muscimol'};
brainRegions = {'sen', 'mtc'};
dvs = {'success rate', ...
       'speed (m/s)', ...
       {'contra (or right)', 'forepaw first rate'}, ...
       {'body angle towards contra', '(or right) side'}};
dvYLims = [0 1; 0 .8; 0 1; -15 15];
minTrial = 0;


% initializations
figure('name', sprintf('min trial: %i', minTrial), 'Color', 'white', 'MenuBar', 'none', 'Position', [-1396 200 500 800])
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
    contraFirstRates = nan(length(mice), length(conditions));
    contraBodyAngles = nan(length(mice), length(conditions)); % angling towards contra side of body
    
    for j = 1:length(mice)
        for k = 1:length(conditions)
            
            conditionBins = brainRegionBins & ...
                            strcmp(sessionInfo.injectedSubstance, conditions{k}) & ...
                            strcmp(sessionInfo.mouse, mice{j});
            sessions = unique(sessionInfo.session(conditionBins));
            
            % containers for session averages for all session for a given mouse in a given condition
            sessionSuccesses = nan(1, length(sessions));
            sessionSpeeds = nan(1, length(sessions));
            sessionContraFirstRates = nan(1, length(sessions));
            sessionContraBodyAngles = nan(1, length(sessions));
            
            for m = 1:length(sessions)
                % get speed and success rate
                sessionBins = strcmp({speedAvoidanceData.session}, sessions{m}) & ...
                              [speedAvoidanceData.trialNum]>=minTrial;
                sessionSuccesses(m) = mean([speedAvoidanceData(sessionBins).isObsAvoided]);
                sessionSpeeds(m) = mean([speedAvoidanceData(sessionBins).avgVel]);
                
                % get contra first paw over rate
                sideOfBrain = sessionInfo.side(strcmp(sessionInfo.session, sessions{m}));
                sessionKinBins = strcmp({kinData.session}, sessions{m}) & [kinData.trial]>minTrial; % inds for kinematic data struct
                sessionContraFirstRates(m) = mean([kinData(sessionKinBins).firstPawOver]-2); % assumes left forepaw is 2 and right forepaw is 3
                if strcmp(sideOfBrain, 'right'); sessionContraFirstRates(m) = 1-sessionContraFirstRates(m); end
                
                % get body angle
                load([getenv('OBSDATADIR') 'sessions\' sessions{m} '\runAnalyzed.mat'], 'bodyAngles');
                framesToAnalyze = getFramesToShow(sessions{m}, true);
                sessionContraBodyAngles(m) = nanmedian(bodyAngles(framesToAnalyze));
                if strcmp(sideOfBrain, 'left'); sessionContraBodyAngles(m) = -sessionContraBodyAngles(m); end
            end
            
            successes(j,k) = nanmean(sessionSuccesses);
            speeds(j,k) = nanmean(sessionSpeeds);
            contraFirstRates(j,k) = nanmean(sessionContraFirstRates);
            contraBodyAngles(j,k) = nanmean(sessionContraBodyAngles);
        end
        
        % plot mouse success
        subplot(dims(1), dims(2), i)
        line([1:length(conditions)] + xJitters(j), successes(j,:), 'color', [.5 .5 .5]); hold on
        scatter([1:length(conditions)] + xJitters(j), successes(j,:), 50, colors(j,:), 'filled');
        
        % plot mouse speed
        subplot(dims(1), dims(2), i+1*length(brainRegions))
        line([1:length(conditions)] + xJitters(j), speeds(j,:), 'color', [.5 .5 .5]); hold on
        scatter([1:length(conditions)] + xJitters(j), speeds(j,:), 50, colors(j,:), 'filled');
        
        % plot mouse fiirst paw over
        subplot(dims(1), dims(2), i+2*length(brainRegions))
        line([1:length(conditions)] + xJitters(j), contraFirstRates(j,:), 'color', [.5 .5 .5]); hold on
        scatter([1:length(conditions)] + xJitters(j), contraFirstRates(j,:), 50, colors(j,:), 'filled');
        
        % plot contra body angle
        subplot(dims(1), dims(2), i+3*length(brainRegions))
        line([1:length(conditions)] + xJitters(j), contraBodyAngles(j,:), 'color', [.5 .5 .5]); hold on
        scatter([1:length(conditions)] + xJitters(j), contraBodyAngles(j,:), 50, colors(j,:), 'filled');
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
        
        % first paw over
        subplot(dims(1), dims(2), i+2*length(brainRegions))
        avg = nanmean(contraFirstRates(:,k));
        line([k-.1 k+.1], [avg avg], 'linewidth', 3, 'color', 'black')
        
        % contra body angles
        subplot(dims(1), dims(2), i+3*length(brainRegions))
        avg = nanmean(contraBodyAngles(:,k));
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
        if i==1; title(brainRegions{j}); set(gca, 'ylim', [0 1]); end
        if j==1; ylabel(dvs{i}); end
        ind = ind+1;
    end
end


saveas(gcf, [getenv('OBSDATADIR') 'figures\muscimolBarPlots.png']);
savefig([getenv('OBSDATADIR') 'figures\muscimolBarPlots.fig'])



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
figure('name', sprintf('min trial: %i', minTrial), 'Color', 'white', 'MenuBar', 'none', 'Position', [-1396 200 800 800])
dims = [length(dvs), length(brainRegions)]; % subplot dimensions

% loop over brain regions
for i = 1:length(brainRegions)
    
    brainRegionBins = strcmp(sessionInfo.brainRegion, brainRegions{i});
    mice = unique(sessionInfo.mouse(brainRegionBins));
    colors = winter(length(mice));
    
    
    for j = 1:length(mice)
        disp(mice{j})
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
            load([getenv('OBSDATADIR') 'sessions\' sessions{m} '\runAnalyzed.mat'], 'bodyAngles');
            framesToAnalyze = getFramesToShow(sessions{m}, true);
            sessionContraBodyAngles(m) = nanmedian(bodyAngles(framesToAnalyze));
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
    xs = linspace(xLims(1)*1.2, xLims(2)*.8, length(mice));
    for j = 1:length(mice)
        text(xs(j), dvYLims(end,1)+(dvYLims(end,2)-dvYLims(end,1))*.2, mice{j}, 'Color', colors(j,:));
    end
end
    


%% pimp figs
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


saveas(gcf, [getenv('OBSDATADIR') 'figures\muscimolSessionProgressPlots.png']);
savefig([getenv('OBSDATADIR') 'figures\muscimolSessionProgressPlotsPlots.fig'])


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
validBins = ~[kinData.isWheelBreak] & [kinData.isObsAvoided] & ismember({kinData.mouse}, mice);
% validBins = ismember({kinData.mouse}, mice);
mtcBins = strcmp({kinData.brainRegion}, 'mtc');
senBins = strcmp({kinData.brainRegion}, 'sen');

% mtc, light on
bins = zeros(1,length(kinData));
figBins = mtcBins & lightOnBins & validBins;
bins(firstPawIpsiBins & musBins & figBins) = 1; % ipsi
bins(firstPawContraBins & musBins & figBins) = 2; % contra
bins(vehBins & figBins) = 3; % vehicle
binLabels = {'ipsi', 'contra', 'vehicle'};
plotObsHeightTrajectories(kinData, bins, binLabels, 'motor cortex lesions, light on')

% mtc, light off
bins = zeros(1,length(kinData));
figBins = mtcBins & ~lightOnBins & validBins;
bins(firstPawIpsiBins & musBins & figBins) = 1; % ipsi
bins(firstPawContraBins & musBins & figBins) = 2; % contra
bins(vehBins & figBins) = 3; % vehicle
binLabels = {'ipsi', 'contra', 'vehicle'};
plotObsHeightTrajectories(kinData, bins, binLabels, 'motor cortex lesions, light off')

% sen, light on
bins = zeros(1,length(kinData));
figBins = senBins & lightOnBins & validBins;
bins(musBins & figBins) = 1; % mus
bins(vehBins & figBins) = 2; % vehicle
binLabels = {'muscimol', 'vehicle'};
plotObsHeightTrajectories(kinData, bins, binLabels, 'sensory cortex lesions, light on')

% sen, light off
bins = zeros(1,length(kinData));
figBins = senBins & ~lightOnBins & validBins;
bins(musBins & figBins) = 1; % mus
bins(vehBins & figBins) = 2; % vehicle
binLabels = {'muscimol', 'vehicle'};
plotObsHeightTrajectories(kinData, bins, binLabels, 'sensory cortex lesions, light off')















