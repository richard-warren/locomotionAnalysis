function data = plotLearning(mice)
% make plots to show how avoidance and speed progresses over time, before
% and after adding wheel break

% settings
touchThresh = 5;
subplotNames = {'success rate', 'velocity (m/s)'};
noBreakExpName = 'obsNoBrBar';
breakExpName = 'obsBrBar';
noBrSessions = 2; % uses the most recent noBrSessions 
brSessions = 5; % uses the first (oldest) brSessions
mouseScatSize = 50;

% initializations
xInds = 1:(noBrSessions + brSessions); % session inds for plots
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'sheet', 'sessions');
sessionBins = ismember(sessionInfo.mouse, mice) &...
              ismember(sessionInfo.experiment, {noBreakExpName, breakExpName}) &...
              sessionInfo.include;
sessionInfo = sessionInfo(sessionBins, :);
mouseColors = hsv(length(mice));



% determine which sessions to include (get last noBrSessions without wheel break and first brSessions with wheel break)
mouseSessions = cell(1,length(mice)); % each entry will be a cell array of session names

for i = 1:length(mice)
    % get inds for no break sessions
    indsNoBr = find(strcmp(sessionInfo.mouse, mice{i}) & ...
                    strcmp(sessionInfo.experiment, noBreakExpName), ...
                    noBrSessions, 'last');
    % get inds for break sessions
    indsBr = find(strcmp(sessionInfo.mouse, mice{i}) & ...
                    strcmp(sessionInfo.experiment, breakExpName), ...
                    brSessions, 'first');
    % get session names for this mouse
    mouseSessions{i} = sessionInfo.session([indsNoBr; indsBr]);
end



% get session data
allSessions = vertcat(mouseSessions{1,:});
data = getSpeedAndObsAvoidanceData(allSessions, sessionInfo, false);
isSuccess = cellfun(@sum, {data.totalTouchFramesPerPaw}) < touchThresh;


% prepare figure
figure('name', 'obsAvoidanceLearningSummary', 'menubar', 'none', 'units', 'pixels', ...
    'position', [100 300 900 600], 'color', [1 1 1], 'inverthardcopy', 'off');



avgVels = nan(length(mice), noBrSessions+brSessions, 2); % last dimension is whether light is on (1) or off (2)
successRates = nan(length(mice), noBrSessions+brSessions, 2);
for i = 1:length(mice)
    
    % get avoidance and avg vel for all sessions, separeted by light on and light off
    for j = 1:length(mouseSessions{i})
        dataBins = strcmp(mouseSessions{i}{j}, {data.session});
        onBins = dataBins & [data.isLightOn];
        offBins = dataBins & ~[data.isLightOn];
        
        avgVels(i,j,1) = nanmean([data(onBins).avgVel]);
        avgVels(i,j,2) = nanmean([data(offBins).avgVel]);
        successRates(i,j,1) = mean(isSuccess(onBins));
        successRates(i,j,2) = mean(isSuccess(offBins));;
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
    line([noBrSessions+.5 noBrSessions+.5], [0 1], 'lineWidth', 1.5, 'color', get(gca, 'xcolor'))
end
xlabel('session #', 'fontweight', 'bold');

% add mouse labels
ys = fliplr(linspace(.2, .8, length(mice)));
for i = 1:length(mice)
    text(xInds(end), ys(i), mice{i}, 'Color', mouseColors(i,:));
end

savefig([getenv('OBSDATADIR') 'learningProgress.fig'])
saveas(gcf, [getenv('OBSDATADIR') 'figures\learningProgress.png']);


