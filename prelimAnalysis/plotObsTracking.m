function plotObsTracking(session, obsTracking)

% analyzeSession generates obsTracking data struct that stores movement of
% wheel and obs for each trial // this function plots these on top of one
% another to visualize how well the movement of the obstalces tracks that
% of the wheel // if obsTracking is not provided, attempts to load
% obsTracking from runAnalyzed.mat, stored in session folder

% settings
plotRows = 9;

% initializations
if ~exist('obsTracking', 'var')
    load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'obsTracking');
end
figure('name', [session ': obstacle tracking of wheel velocity'], ...
    'color', 'white', 'position', [2002 41 1462 940])
plotCols = ceil(length(obsTracking) / plotRows);

% plot
for i = 1:plotRows
    subplot(plotRows, 1, i)
    rowInds = (i-1)*plotCols+1:min(i*plotCols+1, length(obsTracking));
    wheelVelRow = [obsTracking(rowInds).wheelVel];
    obsVelRow = [obsTracking(rowInds).obsVel];

    plot(wheelVelRow); hold on;
    plot(obsVelRow);

    trialStartInds = cumsum(cellfun(@length, {obsTracking(rowInds).times}));
    trialStartInds = [1 trialStartInds(1:end-1)];
    try  % the following breaks on session 180722_005
        scatter(trialStartInds, wheelVelRow(trialStartInds), 20, 'filled', 'black')
        set(gca, 'xlim', [0, length(wheelVelRow)], 'ylim', [0 1], 'box', 'off', 'xtick', [], 'TickDir', 'out')
        text(trialStartInds, zeros(1, length(trialStartInds)), sprintfc('%d', rowInds))
    end
end

legend({'wheel', 'obstacle'})
pause(.001)