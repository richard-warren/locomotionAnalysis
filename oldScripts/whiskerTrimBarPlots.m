function whiskerTrimBarPlots(data, figTitle)




% settings
poolSessions = false; % pool trials across sessions for given mouse if true
touchThresh = 5; % if paw contacts obs for more than touchThresh frames, trial is considered touching
dvs = {'success rate', ...
       'speed (m/s)', ...
       {'body angle towards contra', '(or right) side'}, ...
       'contra (or right) error rate', ...
       'ipsi (or left) error rate'};
dvYLims = [0 1; 0 .8; -15 15; 0 .8; 0 .8];
minTrial = 0;
validBins = [data.trialNum]>=minTrial;
columns = 3;
minConditionNum = 2; % only use a condition after the minConditionNum day of that condition


% initializations
isSuccess = cellfun(@(x) sum(any(x,2)), {data.trialTouchesPerPaw}) < touchThresh;
conditions = unique({data.condition});
rows = ceil(length(dvs)/columns);
figure('name', figTitle, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 300*columns 250*rows], 'inverthardcopy', 'off')



% loop over brain regions
mice = unique({data.mouse});
xJitters = linspace(-.1,.1,length(mice)); xJitters = xJitters-mean(xJitters); % jitters x position of scatter points
colors = hsv(length(mice));

% containers for averages for each mouse for each condition across all sessions
successes = nan(length(mice), length(conditions));
speeds = nan(length(mice), length(conditions)); % rows are mice, columns are conditions (saline, muscimol)
contraBodyAngles = nan(length(mice), length(conditions)); % angling towards contra side of body
contraErrRates = nan(length(mice), length(conditions));
ipsiErrRates = nan(length(mice), length(conditions));

for j = 1:length(mice)
    for k = 1:length(conditions)

        mouseBins = strcmp({data.mouse}, mice{j});
        conditionBins = strcmp({data.condition}, conditions{k}) & ...
                        mouseBins & ...
                        [data.conditionNum]>=minConditionNum & ...
                        validBins;
        sessions = unique({data(conditionBins).session});


        if ~isempty(sessions)
            
            % containers for averages for all sessions of a given mouse in a given condition
            if poolSessions; iters=1; else; iters=length(sessions); end
            [sessionSuccesses, sessionSpeeds, sessionContraBodyAngles, ...
                sessionContraErrRates, sessionIpsiErrRates] = deal(nan(1, iters));

            % get means
            for m = 1:iters

                if poolSessions; bins = conditionBins; else; bins = strcmp({data.session}, sessions{m}); end

                % get speed, success rate
                sessionSuccesses(m) = nanmean(isSuccess(bins));
                sessionSpeeds(m) = nanmean([data(bins).avgVel]);

                % get body angle
                sessionContraBodyAngles(m) = nanmean([data(bins).avgAngle]);
                untrimmedSide = unique({data(strcmp({data.session}, sessions{m})).side});
                if strcmp(untrimmedSide, 'left'); sessionContraBodyAngles(m) = -sessionContraBodyAngles(m); end

                % get contra err rate
                leftErrorRate = nanmean(cellfun(@(x) sum(any(x(:,[1 2]),2))>=touchThresh, {data(bins).trialTouchesPerPaw}));
                rightErrorRate = nanmean(cellfun(@(x) sum(any(x(:,[3 4]),2))>=touchThresh, {data(bins).trialTouchesPerPaw}));
                if strcmp(untrimmedSide, 'left')
                    [sessionContraErrRates, sessionIpsiErrRates] = deal(rightErrorRate, leftErrorRate);
                else
                    [sessionContraErrRates, sessionIpsiErrRates] = deal(leftErrorRate, rightErrorRate);
                end
            end

            successes(j,k) = nanmean(sessionSuccesses);
            speeds(j,k) = nanmean(sessionSpeeds);
            contraBodyAngles(j,k) = nanmean(sessionContraBodyAngles);
            contraErrRates(j,k) = nanmean(sessionContraErrRates);
            ipsiErrRates(j,k) = nanmean(sessionIpsiErrRates);
        end
    end
    

    % plot mouse means
    allDvs = {successes, speeds, contraBodyAngles, contraErrRates, ipsiErrRates};
    for k = 1:length(allDvs)
        subplot(rows, columns, k)
        line([1:length(conditions)] + xJitters(j), allDvs{k}(j,:), 'color', [.2 .2 .2]); hold on
        scatter([1:length(conditions)] + xJitters(j), allDvs{k}(j,:), 50, colors(j,:), 'filled');
    end
end


% plot condition means
for k = 1:length(conditions)
    for m = 1:length(allDvs)
        subplot(rows, columns, m)
        avg = nanmean(allDvs{m}(:,k));
        line([k-.1 k+.1], [avg avg], 'linewidth', 3, 'color', 'black')
        [~,p] = ttest(allDvs{m}(:,1), allDvs{m}(:,2));
        xlabel(sprintf('p=%.3f', p))
    end
end

% add mouse labels
xLims = get(gca, 'xlim');
xs = linspace(xLims(1), xLims(2)*.8, length(mice));
for j = 1:length(mice)
    text(xs(j), dvYLims(end,2), mice{j}, 'Color', colors(j,:));
end


% pimp figs
for i = 1:length(dvs)
    subplot(rows, columns, i);
    set(gca, 'xlim', [0.75 length(conditions)+0.25], 'xtick', 1:length(conditions), 'XTickLabel', conditions, ...
        'YLim', dvYLims(i,:));
    ylabel(dvs{i});
end

% blackenFig