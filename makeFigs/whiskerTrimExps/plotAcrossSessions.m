function whiskerTrimAcrossSessions(data, figTitle)

% settings
dvs = {'success rate', ...
       'speed (m/s)', ...
       {'body angle towards contra', '(or right) side'}, ...
       'contra (or right) error rate', ...
       'ipsi (or left) error rate'};
dvYLims = [0 1; .1 .8; -15 15; 0 .8; 0 .8];
minTrial = 0;
validBins = [data.trialNum]>=minTrial;
touchThresh = 5;
columns = 1;
mouseSymbols = {'o', '+', '*', '.', 'x', 's', 'd'};

% initializations
isSuccess = cellfun(@sum, {data.totalTouchFramesPerPaw}) < touchThresh;
rows = ceil(length(dvs) / columns);
figure('name', figTitle, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600*columns 200*rows], 'InvertHardcopy', 'off')
mice = unique({data.mouse});
conditions = unique({data.condition});
colors = hsv(length(conditions));
scatters = nan(1,length(mice));

[successes, speeds, contraBodyAngles, contraErrRate, ipsiErrRate, conditionNum, sessionNums] = ...
    deal(nan(length(mice), max([data.sessionNum])));

for j = 1:length(mice)
    conditionBins = strcmp({data.mouse}, mice{j});
    sessions = unique({data(conditionBins).session});

    % containers for session averages for all session for a given mouse in a given condition
    

    for m = 1:length(sessions)

        % get speed and success rate
        sessionBins = strcmp({data.session}, sessions{m}) & validBins;
        successes(j,m) = nanmean(isSuccess(sessionBins));
        speeds(j,m) = mean([data(sessionBins).avgVel]);
        untrimmedSide = unique({data(strcmp({data.session}, sessions{m})).side});
        sessionNums(j,m) = unique([data(sessionBins).sessionNum]);

        % get body angle
        contraBodyAngles(j,m) = nanmean([data(sessionBins).avgAngle]);
        if strcmp(untrimmedSide, 'left'); contraBodyAngles(m) = -contraBodyAngles(m); end
        
        % get contra and ipsi err rate
        leftErrorRate = nanmean(cellfun(@(x) sum(any(x(:,[1 2]),2))>=touchThresh, {data(sessionBins).trialTouchesPerPaw}));
        rightErrorRate = nanmean(cellfun(@(x) sum(any(x(:,[3 4]),2))>=touchThresh, {data(sessionBins).trialTouchesPerPaw}));
        if strcmp(untrimmedSide, 'left')
            [contraErrRate(j,m), ipsiErrRate(j,m)] = deal(rightErrorRate, leftErrorRate);
        else
            [contraErrRate(j,m), ipsiErrRate(j,m)] = deal(leftErrorRate, rightErrorRate);
        end

        % get condition whether or not muscimol
        sessionBins = strcmp({data.session}, sessions{m});
        conditionNum(j,m) = find(strcmp(unique({data(sessionBins).condition}), conditions));
    end


    % plot mouse data
    allDvs = {successes, speeds, contraBodyAngles, contraErrRate, ipsiErrRate};
    for k = 1:length(allDvs)
        subplot(rows, columns, k)
        line(sessionNums(j,:), allDvs{k}(j,:), 'color', [.5 .5 .5]); hold on
        scatters(j) = scatter(sessionNums(j,:), allDvs{k}(j,:), 50, colors(conditionNum(j,:),:), mouseSymbols{j});
    end
end


% add mouse labels
xLims = get(gca, 'xlim');
xs = linspace(1, xLims(2)*.4, length(mice));
for j = 1:length(conditions)
    text(xs(j), dvYLims(end,2), conditions{j}, 'Color', colors(j,:));
end
    


% pimp figs
for i = 1:length(dvs)
    subplot(rows, columns, i);
    xLims = get(gca, 'xlim');
    set(gca, 'xlim', [0.5 xLims(2)+.5], 'xtick', 1:xLims(2), 'YLim', dvYLims(i,:));
    ylabel(dvs{i});
end
xlabel('session number')
legend(scatters, mice, 'location', 'northwest', 'Color', 'none')

% blackenFig

