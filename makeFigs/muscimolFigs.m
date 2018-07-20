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

%% plot vehicle vs mus baseline kinematics for all mice and all brain regions per mouse


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


%% plot avoidance


% settings
conditions = {'saline', 'muscimol'};
brainRegions = {'sen', 'mtc'};
dvs = {'success rate', 'speed (m/s)'};
speedYLims = [0 .8];



% initializations
close all; figure('Color', 'white', 'MenuBar', 'none', 'Position', [-1396 457 560 420])
dims = [length(dvs), length(brainRegions)]; % subplot dimensions

% loop over brain regions
for i = 1:length(brainRegions)
    
    
    brainRegionBins = strcmp(sessionInfo.brainRegion, brainRegions{i});
    mice = unique(sessionInfo.mouse(brainRegionBins));
    colors = winter(length(mice));
    
    speeds = nan(length(mice), length(conditions)); % rows are mice, columns are conditions (saline, muscimol)
    successes = nan(length(mice), length(conditions));
    
    for j = 1:length(mice)
        for k = 1:length(conditions)
            
            conditionBins = brainRegionBins & ...
                            strcmp(sessionInfo.injectedSubstance, conditions{k}) & ...
                            strcmp(sessionInfo.mouse, mice{j});
            sessions = unique(sessionInfo.session(conditionBins));
            
            sessionSuccesses = nan(1, length(sessions));
            sessionSpeeds = nan(1, length(sessions));
            for m = 1:length(sessions)
                sessionBins = strcmp({speedAvoidanceData.session}, sessions{m});
                sessionSuccesses(m) = mean([speedAvoidanceData(sessionBins).isObsAvoided]);
                sessionSpeeds(m) = mean([speedAvoidanceData(sessionBins).avgVel]);
            end
            
            successes(j,k) = mean(sessionSuccesses);
            speeds(j,k) = mean(sessionSpeeds);
        end
        
        % plot mouse success
        subplot(dims(1), dims(2), sub2ind(dims, i, 1))
        line(1:length(conditions), successes(j,:), 'color', [.5 .5 .5]); hold on
        scatter(1:length(conditions), successes(j,:), 50, colors(j,:), 'filled');
        
        % plot mouse speed
        subplot(dims(1), dims(2), sub2ind(dims, i, 2))
        line(1:length(conditions), speeds(j,:), 'color', [.5 .5 .5]); hold on
        scatter(1:length(conditions), speeds(j,:), 50, colors(j,:), 'filled');
    end
    
    % plot condition means
    for k = 1:length(conditions)
        % success
        subplot(dims(1), dims(2), sub2ind(dims, i, 1))
        avg = nanmean(successes(:,k));
        line([k-.1 k+.1], [avg avg], 'linewidth', 3, 'color', 'black')
        
        % speed
        subplot(dims(1), dims(2), sub2ind(dims, i, 2))
        avg = nanmean(speeds(:,k));
        line([k-.1 k+.1], [avg avg], 'linewidth', 3, 'color', 'black')
    end
end



% pimp figs
ind = 1;
for i = 1:length(dvs)
    for j = 1:length(brainRegions)
        subplot(dims(1), dims(2), sub2ind(dims,j,i));
        
        set(gca, 'xlim', [0.75 length(conditions)+0.25], 'xtick', 1:length(conditions), 'XTickLabel', conditions);
        if i==1; title(brainRegions{j}); set(gca, 'ylim', [0 1]); end
        if i==2; set(gca, 'ylim', speedYLims); end
        if j==1; ylabel(dvs{i}); end
        ind = ind+1;
    end
end




