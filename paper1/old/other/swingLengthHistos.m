%% swing length histograms

% settings
xGridLims = [.02 .12];
yLims = [0 .3];
binWidth = .005;
colors = winter(2);
controlColor = [.65 .65 .65];

% initializations
numModSteps = reshape([data.modStepNum],4,length(data))';
modifiedSwingLengths = {data.modifiedSwingLengths}; modifiedSwingLengths = cat(1, modifiedSwingLengths{:});
controlSwingLengths = {data.controlSwingLengths}; controlSwingLengths = cat(1, controlSwingLengths{:});

figure('color', 'white', 'menubar', 'none', 'position', [150 400 300*binNum 350]);

for h = 1:binNum

    ax = subaxis(1, binNum , h, 'spacing', .02);
    binBins = (bins==h)';
    oneStepBins = binBins & numModSteps(:,3)==1;
    twoStepBins = binBins & numModSteps(:,3)==2;
    oneTwoRatio = sum(oneStepBins) / (sum(oneStepBins) + sum(twoStepBins));

    % one step histo
    if any(oneStepBins)
        h1 = histogram(modifiedSwingLengths(oneStepBins,3), 'binwidth', binWidth); hold on
        counts = get(h1,'bincounts');
        set(h1, 'facecolor', colors(2,:), 'normalization', 'count', ...
            'bincounts', (counts/sum(counts)) * oneTwoRatio);
    end

    % two step histo
    h2 = histogram(modifiedSwingLengths(twoStepBins,3), 'binwidth', binWidth); hold on;        
    counts = get(h2,'bincounts');
    set(h2, 'facecolor', colors(1,:), 'normalization', 'count', ... 
        'bincounts', (counts/sum(counts)) * (1-oneTwoRatio));

    % control histo
    h3 = histogram(controlSwingLengths(binBins,3), 'binwidth', binWidth); hold on;
    counts = get(h3,'bincounts');
    set(h3, 'facecolor', controlColor, 'normalization', 'count', ...
        'bincounts', (counts/sum(counts)));

    % set apearance
    set(ax, 'box', 'off', 'xlim', xGridLims, 'ylim', yLims, 'tickdir', 'out')
    set(ax, 'ytick', [], 'ylabel', []); ax.YAxis.Visible = 'off';
end

legend('modified swing lengths (lengthened)', 'modified swing lengths (shortened)', 'control swing lengths')

saveas(gcf, [getenv('OBSDATADIR') 'figures\swingLengthHistograms.png']);


