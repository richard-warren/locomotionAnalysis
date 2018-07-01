function plotAllSessionSpeeds(experiments, posRange)

% settings
velThresh = .4;
yLims = [0 .8]; % m/s
scatterSize = 15;
numRows = 14;

% reasonable posRange is the following: [-.08 .08] + 0.3820

% get velocity data
data = getAllSessionTrialSpeeds(experiments, posRange);
data = data(~[data.lightOn]);
sessions = unique({data.session});


% initializations
figure('color', 'white', 'units', 'pixels', 'position', [100 100 1600 900],...
    'name', 'allSessionSpeeds', 'menubar', 'none');
numCols = ceil(length(sessions) / numRows);

% sort by mouse
mice = cell(1,length(sessions));
for i = 1:length(sessions)
    sessionBins = strcmp(sessions{i}, {data.session});
    mice(i) = unique({data(sessionBins).mouse});
end
[~, sortInds] = sort(mice);
sessions = sessions(sortInds);


% plot everything
for i = 1:length(sessions)
    
    subaxis(numCols, numRows, i, 'spacing', 0.05, 'margin', .05)
    sessionBins = strcmp(sessions{i}, {data.session});
    cmap = winter(sum(sessionBins));
    
    scatter(zeros(1,sum(sessionBins)), [data(sessionBins).vel], scatterSize, cmap, 'filled', 'jitter', 'on');
    line([-.5 .5], [velThresh velThresh], 'linewidth', 2, 'color', mean(cmap,1));
    
    % pimp fig
    firstTrialInd = find(sessionBins,1,'first');
    title(sprintf('%s (%s)', data(firstTrialInd).mouse, data(firstTrialInd).session(1:end-4)), 'fontweight', 'normal')
    set(gca, 'ylim', yLims, 'xcolor', 'none')
    
end

