finalSessionsToShow = 4;

sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'sheet', 'sessions');
mice = {'sen2', 'sen3'};

mouseBins = false(1,height(sessionInfo));

for i = 1:length(mice)
    mouseBins(find(strcmp(sessionInfo.mouse, mice{i}), finalSessionsToShow, 'last')) = true;
end

data = getSpeedAndObsAvoidanceData(sessionInfo.session(mouseBins), sessionInfo, false);


%%
% prepare figure
figure('name', 'obsAvoidanceLearningSummary', 'menubar', 'none', 'units', 'pixels', ...
    'position', [100 300 900 600], 'color', [1 1 1], 'inverthardcopy', 'off');
xInds = 1:finalSessionsToShow;
mouseColors = hsv(length(mice));
mouseScatSize = 100;
subplotNames = {'success', 'speed'};

ventralTouchInds = find(contains(data(1).touchClassNames, {'fore_ventral', 'hind_ventral_low'}));
isSuccess = cellfun(@(x) sum(ismember(x, ventralTouchInds)), {data.trialTouches}) < 1;
avgVels = nan(length(mice), finalSessionsToShow, 2); % last dimension is whether light is on (1) or off (2)
successRates = nan(length(mice), finalSessionsToShow, 2);
for i = 1:length(mice)
    
    mouseSessions = unique({data(strcmp({data.mouse}, mice{i})).session});
    
    % get avoidance and avg vel for all sessions, separeted by light on and light off
    for j = 1:length(mouseSessions)
        dataBins = strcmp(mouseSessions{j}, {data.session});
        onBins = dataBins & [data.isLightOn];
        offBins = dataBins & ~[data.isLightOn];
        
        avgVels(i,j,1) = nanmean([data(onBins).avgVel]);
        avgVels(i,j,2) = nanmean([data(offBins).avgVel]);
        successRates(i,j,1) = mean(isSuccess(onBins));
        successRates(i,j,2) = mean(isSuccess(offBins));
    end
    
    % plot avoidance
    subplot(2,1,1)
    plot(xInds, successRates(i,:,1), 'Color', mouseColors(i,:)); hold on
    plot(xInds, successRates(i,:,2), 'Color', mouseColors(i,:));
    scatter(xInds, successRates(i,:,2), mouseScatSize, mouseColors(i,:), ...
        'filled', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', mouseColors(i,:), 'LineWidth', 1.5); % light off
    scatter(xInds, successRates(i,:,1), mouseScatSize, mouseColors(i,:), 'filled', 'MarkerFaceAlpha', .8) % light on
    
    % plot speed
    subplot(2,1,2)
    plot(xInds, avgVels(i,:,1), 'Color', mouseColors(i,:)); hold on
    plot(xInds, avgVels(i,:,2), 'Color', mouseColors(i,:));
    scatter(xInds, avgVels(i,:,2), mouseScatSize, 'filled', ...
        'MarkerFaceColor', 'white', 'MarkerEdgeColor', mouseColors(i,:), 'LineWidth', 1.5); % light off
    scatter(xInds, avgVels(i,:,1), mouseScatSize, mouseColors(i,:), 'filled', 'MarkerFaceAlpha', .8); hold on % light on
    
end

% pimp fig
for i = 1:2
    subplot(2,1,i)
    set(gca, 'box', 'off', 'xtick', xInds, 'xlim', [xInds(1)-.5 xInds(end)], 'ylim', [0 1])
    ylabel(subplotNames{i}, 'fontweight', 'bold')
end
xlabel('session #', 'fontweight', 'bold');

% add mouse labels
ys = fliplr(linspace(.2, .8, length(mice)));
for i = 1:length(mice)
    text(xInds(end), ys(i), mice{i}, 'Color', mouseColors(i,:));
end