function manipulationAcrossSessions(data, conditions, figTitle)

% settings
dvs = {'success rate', ...
       'speed (m/s)', ...
       {'body angle towards contra', '(or right) side'}};
dvYLims = [0 1; .1 .8; -15 15];
minTrial = 0;
validBins = [data.trialNum]>=minTrial & ~[data.isLightOn];
touchThresh = 5;

% initializations
brainRegions = unique({data.brainRegion});
isSuccess = cellfun(@sum, {data.totalTouchFramesPerPaw}) < touchThresh;
dims = [length(dvs), length(brainRegions)]; % subplot dimensions
figure('name', figTitle, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 200 800 600], 'InvertHardcopy', 'off')

% loop over brain regions
for i = 1:length(brainRegions)
    
    brainRegionBins = strcmp({data.brainRegion}, brainRegions{i});
    mice = unique({data(brainRegionBins).mouse});
    colors = hsv(length(mice));
    
    
    for j = 1:length(mice)
        conditionBins = brainRegionBins & ...
                        strcmp({data.mouse}, mice{j});
        sessions = unique({data(conditionBins).session});

        % containers for session averages for all session for a given mouse in a given condition
        sessionSuccesses = nan(1, length(sessions));
        sessionSpeeds = nan(1, length(sessions));
        sessionContraBodyAngles = nan(1, length(sessions));
        conditionOneBins = false(1, length(sessions));
        sessionNums = nan(1,length(sessions));

        for m = 1:length(sessions)
            
            % get speed and success rate
            sessionBins = strcmp({data.session}, sessions{m}) & validBins;
            sessionSuccesses(m) = nanmean(isSuccess(sessionBins));
            sessionSpeeds(m) = mean([data(sessionBins).avgVel]);
            sideOfBrain = unique({data(strcmp({data.session}, sessions{m})).side});
            sessionNums(m) = unique([data(sessionBins).sessionNum]);

            % get body angle
            load([getenv('OBSDATADIR') 'sessions\' sessions{m} '\runAnalyzed.mat'], ...
                'bodyAngles', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps');
            load([getenv('OBSDATADIR') 'sessions\' sessions{m} '\run.mat'], 'breaks');
            bodyAngles = getTrialBodyAngles(bodyAngles, obsOnTimes, obsOffTimes, frameTimeStamps, breaks);
            sessionContraBodyAngles(m) = nanmedian(bodyAngles);
            if strcmp(sideOfBrain, 'left'); sessionContraBodyAngles(m) = -sessionContraBodyAngles(m); end
            
            % get condition whether or not muscimol
            sessionBins = strcmp({data.session}, sessions{m});
            conditionOneBins(m) = strcmp(unique({data(sessionBins).condition}), conditions{1});
        end
        
        
        % plot mouse data
        allDvs = {sessionSuccesses, sessionSpeeds, sessionContraBodyAngles};
        for k = 1:length(allDvs)
            subplot(dims(1), dims(2), i+(k-1)*length(brainRegions))
            line(sessionNums, allDvs{k}, 'color', [.5 .5 .5]); hold on
            scatter(sessionNums(conditionOneBins), allDvs{k}(conditionOneBins), 50, colors(j,:), 'filled'); % muscimol
            scatter(sessionNums(~conditionOneBins), allDvs{k}(~conditionOneBins), 50, colors(j,:), 'LineWidth', 2); % saline
        end
        

    end
    
%     % add mouse labels
%     xLims = get(gca, 'xlim');
%     xs = linspace(1, xLims(2)*.8, length(mice));
%     for j = 1:length(mice)
%         text(xs(j), dvYLims(end,1)+(dvYLims(end,2)-dvYLims(end,1))*.2, mice{j}, 'Color', colors(j,:));
%     end
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

blackenFig

