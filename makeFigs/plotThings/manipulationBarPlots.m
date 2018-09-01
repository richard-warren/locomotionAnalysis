function manipulationBarPlots(data, conditions, figTitle)



% settings
touchThresh = 5; % if paw contacts obs for more than touchThresh frames, trial is considered touching
dvs = {'success rate', ...
       'speed (m/s)', ...
       {'body angle towards contra', '(or right) side'}, ...
       'contra (or right) error rate', ...
       'ipsi (or left) error rate'};
dvYLims = [0 1; 0 .8; -15 15; 0 1; 0 1];
minTrial = 0;
validBins = [data.trialNum]>=minTrial & ~[data.isLightOn];


% initializations
brainRegions = unique({data.brainRegion});
isSuccess = cellfun(@sum, {data.totalTouchFramesPerPaw}) < touchThresh;
dims = [length(dvs), length(brainRegions)]; % subplot grid
figure('name', figTitle, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 length(brainRegions)*250 900], 'inverthardcopy', 'off')



% loop over brain regions
for i = 1:length(brainRegions)
    
    brainRegionBins = strcmp({data.brainRegion}, brainRegions{i});
    mice = unique({data(brainRegionBins).mouse});
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
                            strcmp({data.condition}, conditions{k}) & ...
                            strcmp({data.mouse}, mice{j});
            sessions = unique({data(conditionBins).session});
            
            % containers for session averages for all session for a given mouse in a given condition
            sessionSuccesses = nan(1, length(sessions));
            sessionSpeeds = nan(1, length(sessions));
            sessionContraBodyAngles = nan(1, length(sessions));
            sessionContraErrRates = nan(1, length(sessions));
            sessionIpsiErrRates = nan(1, length(sessions));
            
            for m = 1:length(sessions)
                % get speed and success rate
                sessionBins = strcmp({data.session}, sessions{m}) & validBins;
                sessionSuccesses(m) = mean(isSuccess(sessionBins));
                sessionSpeeds(m) = mean([data(sessionBins).avgVel]);
                
                % get body angle
                load([getenv('OBSDATADIR') 'sessions\' sessions{m} '\runAnalyzed.mat'], ...
                    'bodyAngles', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'touchesPerPaw');
                load([getenv('OBSDATADIR') 'sessions\' sessions{m} '\run.mat'], 'breaks');
                bodyAngles = getTrialBodyAngles(bodyAngles, obsOnTimes, obsOffTimes, frameTimeStamps, breaks);
                sessionContraBodyAngles(m) = nanmedian(bodyAngles);
                sideOfBrain = unique({data(strcmp({data.session}, sessions{m})).side});
                if strcmp(sideOfBrain, 'left'); sessionContraBodyAngles(m) = -sessionContraBodyAngles(m); end
                
                % get contra err rate
                leftErrors = 0;
                rightErrors = 0;
                for n = 1:length(obsOnTimes)
                    trialInds = frameTimeStamps>obsOnTimes(n) & frameTimeStamps<obsOffTimes(n);
                    leftErrors = leftErrors + any(any(touchesPerPaw(trialInds,[1 2])));
                    rightErrors = rightErrors + any(any(touchesPerPaw(trialInds,[3 4]),2)>0);
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
        
        
        allDvs = {successes, speeds, contraBodyAngles, contraErrRates, ipsiErrRates};
        for k = 1:length(allDvs)
            subplot(dims(1), dims(2), i+(k-1)*length(brainRegions))
            line([1:length(conditions)] + xJitters(j), allDvs{k}(j,:), 'color', [.2 .2 .2]); hold on
            scatter([1:length(conditions)] + xJitters(j), allDvs{k}(j,:), 50, colors(j,:), 'filled');
        end
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
    
%     % add mouse labels
%     xLims = get(gca, 'xlim');
%     xs = linspace(xLims(1)*1.2, xLims(2)*.8, length(mice));
%     for j = 1:length(mice)
%         text(xs(j), dvYLims(end,1)+(dvYLims(end,2)-dvYLims(end,1))*.2, mice{j}, 'Color', colors(j,:));
%     end

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

blackenFig