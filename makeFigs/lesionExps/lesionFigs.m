%% load session info

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'lesionNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);



%% get kinematic data
obsPos = -0.0087;
loadPreviousData = true;

if loadPreviousData
    load([getenv('OBSDATADIR') 'kinematicData\lesionKinematicData.mat'], 'data');
    kinData = getKinematicData4(sessionInfo.session, data, obsPos);
else
    kinData = getKinematicData4(sessionInfo.session, [], obsPos);
end
data = kinData; save([getenv('OBSDATADIR') 'kinematicData\lesionKinematicData.mat'], 'data');



%% get avoidance and speed data

speedAvoidanceData = getSpeedAndObsAvoidanceData(sessionInfo.session, false);
data = speedAvoidanceData; save([getenv('OBSDATADIR') 'kinematicData\lesionSpeedAvoidanceData.mat'], 'data');

%% plot baseline kinematics

mice = unique(sessionInfo.mouse);

for j = 1:length(mice)
    
    mouseBins = strcmp(sessionInfo.mouse, mice{j});
    brainRegion = unique(sessionInfo.brainRegion(mouseBins));
    preBins = strcmp(sessionInfo.preOrPost, 'pre');
    postBins = strcmp(sessionInfo.preOrPost, 'post');
    sides = unique(sessionInfo.side(mouseBins));
    
    for k = 1:length(sides)
        sideBins = strcmp(sessionInfo.side, sides{k});
        preSessions = sessionInfo.session(mouseBins & sideBins & preBins);
        postSessions = sessionInfo.session(mouseBins & sideBins & postBins);
        
        conditionBins = nan(1,length(kinData));
        conditionBins(ismember({kinData.session}, preSessions)) = 1;
        conditionBins(ismember({kinData.session}, postSessions)) = 2;
        
        conditionLabels = {['pre - ' sides{k} ' ' brainRegion{1}], ['post - ' sides{k} ' ' brainRegion{1}]};
        plotBaselineKinematics(kinData, conditionBins, conditionLabels, mice{j})
    end
end


%% plot dv averages

% settings
conditions = {'pre', 'post'};
brainRegions = {'sen'};
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
figure('name', sprintf('min trial: %i', minTrial), 'Color', 'white', 'MenuBar', 'none', 'Position', [0 50 250*length(brainRegions) 900])
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
                            strcmp(sessionInfo.preOrPost, conditions{k}) & ...
                            strcmp(sessionInfo.mouse, mice{j});
            sessions = unique(sessionInfo.session(conditionBins));
%             if strcmp(conditions{k}, 'post'); sessions = sessions(end); end
            
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


saveas(gcf, [getenv('OBSDATADIR') 'figures\lesionBarPlots.png']);
savefig([getenv('OBSDATADIR') 'figures\lesionBarPlots.fig'])




%% plot dvs across sessions

% settings
conditions = {'pre', 'post'};
brainRegions = {'sen'};
dvs = {'success rate', ...
       'speed (m/s)', ...
       {'body angle towards contra', '(or right) side'}};
dvYLims = [.4 1; .1 .8; -15 15];
minTrial = 0;
validTrials = ~[speedAvoidanceData.isLightOn] & [speedAvoidanceData.trialNum]>=minTrial;


% initializations
figure('name', sprintf('min trial: %i', minTrial), 'Color', 'white', 'MenuBar', 'none', 'Position', [0 200 400*length(brainRegions) 600])
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
        preBins = false(1, length(sessions));

        for m = 1:length(sessions)
            % get speed and success rate
            sessionBins = strcmp({speedAvoidanceData.session}, sessions{m}) & ...
                          validTrials;
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
            preBins(m) = strcmp(sessionInfo.preOrPost(sessionBin), 'pre');
        end
        
        
        % plot mouse data
        xInds = 1:length(sessions);
        
        % plot mouse success
        subplot(dims(1), dims(2), i)
        line(1:length(sessions), sessionSuccesses, 'color', [.5 .5 .5]); hold on
        scatter(xInds(preBins), sessionSuccesses(preBins), 50, colors(j,:), 'filled'); % muscimol
        scatter(xInds(~preBins), sessionSuccesses(~preBins), 50, colors(j,:)); % saline
        
        % plot mouse speed
        subplot(dims(1), dims(2), i+1*length(brainRegions))
        line(1:length(sessions), sessionSpeeds, 'color', [.5 .5 .5]); hold on
        scatter(xInds(preBins), sessionSpeeds(preBins), 50, colors(j,:), 'filled'); % muscimol
        scatter(xInds(~preBins), sessionSpeeds(~preBins), 50, colors(j,:)); % saline
        
        % plot contra body angle
        subplot(dims(1), dims(2), i+2*length(brainRegions))
        line(1:length(sessions), sessionContraBodyAngles, 'color', [.5 .5 .5]); hold on
        scatter(xInds(preBins), sessionContraBodyAngles(preBins), 50, colors(j,:), 'filled'); % muscimol
        scatter(xInds(~preBins), sessionContraBodyAngles(~preBins), 50, colors(j,:)); % saline;
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

saveas(gcf, [getenv('OBSDATADIR') 'figures\lesionWithinSessionProgressPlots.png']);
savefig([getenv('OBSDATADIR') 'figures\lesionWithinSessionProgressPlotsPlots.fig'])

%% plot kinematics



% incorporate condition information into kinData struct
for i = 1:length(kinData)
    sessionInfoBin = strcmp(sessionInfo.session, kinData(i).session);
    brainRegion = sessionInfo.brainRegion{sessionInfoBin};
    preOrPost = sessionInfo.preOrPost{sessionInfoBin};
    if strcmp(sessionInfo.side{sessionInfoBin}, 'left')
        contraLimb = 3;
    elseif strcmp(sessionInfo.side{sessionInfoBin}, 'right')
        contraLimb = 2;
    else
        contraLimb = nan;
    end
    
    kinData(i).brainRegion = brainRegion;
    kinData(i).preOrPost = preOrPost;
    kinData(i).contraLimb = contraLimb;
end

% get trial bins

mice = unique(sessionInfo.mouse);

firstPawContraBins = [kinData.firstPawOver]==[kinData.contraLimb];
firstPawIpsiBins = [kinData.firstPawOver]~=[kinData.contraLimb];
preBins = strcmp({kinData.preOrPost}, 'pre');
postBins = strcmp({kinData.preOrPost}, 'post');
lightOnBins = [kinData.isLightOn];
mtcBins = strcmp({kinData.brainRegion}, 'mtc');
senBins = strcmp({kinData.brainRegion}, 'sen');
validBins = ~[kinData.isWheelBreak] & [kinData.isObsAvoided] & ismember({kinData.mouse}, mice);

% % mtc, light on
% bins = zeros(1,length(kinData));
% figBins = mtcBins & validBins;
% bins(firstPawIpsiBins & preBins & figBins) = 1; % ipsi
% bins(firstPawContraBins & preBins & figBins) = 2; % contra
% bins(postBins & figBins) = 3; % vehicle
% binLabels = {'ipsi', 'contra', 'vehicle ipsi'};
% plotObsHeightTrajectories(kinData, bins, binLabels, 'motor cortex lesions, light on')
% 
% % mtc, light off
% bins = zeros(1,length(kinData));
% figBins = mtcBins & ~lightOnBins & validBins;
% bins(firstPawIpsiBins & preBins & figBins) = 1; % ipsi
% bins(firstPawContraBins & preBins & figBins) = 2; % contra
% bins(postBins & figBins) = 3; % vehicle
% binLabels = {'ipsi', 'contra', 'vehicle'};
% plotObsHeightTrajectories(kinData, bins, binLabels, 'motor cortex lesions, light off')

% sen, light on
bins = zeros(1,length(kinData));
figBins = senBins & lightOnBins & validBins;
bins(preBins & figBins) = 1;
bins(postBins & figBins) = 2;
binLabels = {'pre', 'post'};
plotObsHeightTrajectories(kinData, bins, binLabels, 'sensory cortex lesions, light on')

% sen, light off
bins = zeros(1,length(kinData));
figBins = senBins & ~lightOnBins & validBins;
bins(preBins & figBins) = 1;
bins(postBins & figBins) = 2;
binLabels = {'pre', 'post'};
plotObsHeightTrajectories(kinData, bins, binLabels, 'sensory cortex lesions, light off')









