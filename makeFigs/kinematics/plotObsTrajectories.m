

% settings
sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003', ...
            '180125_001', '180125_002', '180125_003'};
phaseBinNum = 3;
speedBinNum = 3;


% initializations
data = getKinematicData(sessions);
dataNew = data([data.oneSwingOneStance]);

% get speed and phase bins
phaseBinEdges = prctile([dataNew.swingStartDistance], linspace(0,100,phaseBinNum+1));
phaseBins = discretize([dataNew.swingStartDistance], phaseBinEdges);
speedBinEdges = prctile([dataNew.vel], linspace(0,100,speedBinNum+1));
speedBins = discretize([dataNew.vel], speedBinEdges);

% create speed and phase labels
phaseLabels = cell(1,phaseBinNum);
for i = 1:phaseBinNum; phaseLabels{i} = sprintf('%.3f', mean([dataNew(phaseBins==i).swingStartDistance])); end
speedLabels = cell(1,speedBinNum);
for i = 1:speedBinNum; speedLabels{i} = sprintf('%.3f', mean([dataNew(speedBins==i).vel])); end


%% sperm plots

% settings
xLims = [-.1 .1];
tracesPerPlot = 15;

figure; pimpFig

for g = 1:speedBinNum
    for h = 1:phaseBinNum

        axInd = sub2ind([phaseBinNum speedBinNum], h, g);
        ax = subaxis(speedBinNum, phaseBinNum ,axInd, ...
            'spacing', .01, 'padding', .01, 'margin', .01);
        
        bins = find(speedBins==g & phaseBins==h);
        plotTraces = min(tracesPerPlot, length(bins));
        dataInds = randperm(length(bins), tracesPerPlot);
        dataInds = bins(dataInds);

        for i = dataInds
            for j = 2:3

                realInds = ~isnan(dataNew(i).modifiedStepIdentities(:,j));
                steps = unique(dataNew(i).modifiedStepIdentities(realInds,j));
                colors = winter(length(steps));

                for k = steps'

                    % plot x and y trajectories
                    trialInds = dataNew(i).modifiedStepIdentities(:,j)==k;
                    x = dataNew(i).locations(trialInds,1,j);
                    y = dataNew(i).locations(trialInds,2,j);
                    plot(x, y, 'color', colors(k,:)); hold on

                    % scatter dots at start of each swing
                    scatter(x(end), y(end), 100, colors(k,:), 'filled'); hold on

                    % scatter position of swing foot at obsPos
                    if j==3
                        scatter(dataNew(i).locations(dataNew(i).obsPosInd,1,j), dataNew(i).locations(dataNew(i).obsPosInd,2,j), ...
                            100, [0 0 0], 'x'); hold on
                    end
                end
            end
        end

        % set appearance
        set(gca, 'dataaspectratio', [1 1 1], 'xlim', xLims, 'box', 'off', 'tickdir', 'out', 'ydir', 'reverse', ...
            'xtick', [], 'ytick', []);
        line([0 0], get(gca,'ylim'), 'color', [0 0 0], 'linewidth', 3)
        if g==speedBinNum; xlabel(['swing start distance (m): ' phaseLabels{h}]); end
        if h==1; ylabel(['speed (m/s): ' speedLabels{g}]); end
    end
end

saveas(gcf, [getenv('OBSDATADIR') 'figures\trialKinematics.png']);


%% histograms

% settings
xLims = [.02 .12];
yLims = [0 .4];
binWidth = .005;
colors = winter(2);
controlColor = [.65 .65 .65];

% initializations
numModSteps = reshape([dataNew.modStepNum],4,length(dataNew))';
modifiedSwingLengths = {dataNew.modifiedSwingLengths}; modifiedSwingLengths = cat(1, modifiedSwingLengths{:});
controlSwingLengths = {dataNew.controlSwingLengths}; controlSwingLengths = cat(1, controlSwingLengths{:});

figure; pimpFig

for g = 1:speedBinNum
    for h = 1:phaseBinNum
        
        axInd = sub2ind([phaseBinNum speedBinNum], h, g);
        ax = subaxis(speedBinNum, phaseBinNum , axInd);
        bins = (speedBins==g & phaseBins==h)';
        oneStepBins = bins & numModSteps(:,3)==1;
        twoStepBins = bins & numModSteps(:,3)==2;
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
        h3 = histogram(controlSwingLengths(bins,3), 'binwidth', binWidth); hold on;
        counts = get(h3,'bincounts');
        set(h3, 'facecolor', controlColor, 'normalization', 'count', ...
            'bincounts', (counts/sum(counts)));
        
        % set apearance
        set(ax, 'box', 'off', 'xlim', xLims, 'ylim', yLims, 'tickdir', 'out')
        if g==speedBinNum; xlabel(['stance foot distance (m): ' phaseLabels{h}]); end
        if g<speedBinNum; set(ax, 'xticklabel', []); end
        if h==1; ylabel(['speed (m/s): ' speedLabels{g}]); end
        if h>1; set(ax, 'ytick', [], 'ylabel', []); ax.YAxis.Visible = 'off'; end
    end
end

legend('modified swing lengths (lengthened)', 'modified swing lengths (shortened)', 'control swing lengths')

saveas(gcf, [getenv('OBSDATADIR') 'figures\swingLengthHistograms.png']);

%% scatter

% settings
distanceLims = [phaseBinEdges(1) phaseBinEdges(2)];
velLims = [speedBinEdges(1) speedBinEdges(2)];
smps = 9;

% initializations
modifiedSwingLengths = {dataNew.modifiedSwingLengths}; modifiedSwingLengths = cat(1, modifiedSwingLengths{:});

controlSwingLengths = cellfun(@(x) mean(x,1), {dataNew.controlSwingLengths}, 'uniformoutput', 0);
controlSwingLengths = cat(1, controlSwingLengths{:});
deltaLength = modifiedSwingLengths(:,3) - controlSwingLengths(:,3);

% validBins = deltaLength<.08;

[xq, yq] = meshgrid(linspace(distanceLims(1),distanceLims(2),smps), linspace(velLims(1),velLims(2),smps));
heatMap = griddata([dataNew.swingStartDistance], [dataNew.vel], ...
    deltaLength, xq, yq);

% close all; figure('color', [1 1 1]);
% imagesc('xdata', xq(1,:), 'ydata', yq(:,1), 'cdata', heatMap);
% set(gca, 'xlim', distanceLims, 'ylim', velLims)
% xlabel('stance paw distance (m)');
% ylabel('speed (m/s)');
% pimpFig

[~, sortInds] = sort(deltaLength);
colors = winter(length(dataNew));
close all; figure; scatter([dataNew(sortInds).swingStartDistance], [dataNew(sortInds).vel], ...
    100, colors, 'filled');
xlabel('swing start distance')
ylabel('speed')

xLims = get(gca, 'xlim'); yLims = get(gca, 'ylim');
for i = 2:phaseBinNum; line([phaseBinEdges(i) phaseBinEdges(i)],yLims); end
for i = 2:speedBinNum; line(xLims,[speedBinEdges(i) speedBinEdges(i)]); end
set(gca,'xlim',xLims,'ylim',yLims,'ydir','reverse')
pimpFig

saveas(gcf, [getenv('OBSDATADIR') 'figures\deltaSwingScatter.png']);


%% average trajectories


% settings
controlColor = [.65 .65 .65];
xLims = [-.1 .05];
linWid = 4;
colors = winter(2);

% initializations
figure; pimpFig;
numModSteps = reshape([dataNew.modStepNum],4,length(dataNew))';
obsPosIndInterps = reshape([dataNew.obsPosIndInterp],4,length(dataNew))';


for g = 1:speedBinNum
    for h = 1:phaseBinNum
        
        % get subplot bins
        axInd = sub2ind([phaseBinNum speedBinNum], h, g);
        ax = subaxis(speedBinNum, phaseBinNum , axInd, ...
            'spacing', .01, 'padding', .01, 'margin', .01);
        bins = (speedBins==g & phaseBins==h)';

        % get subplot bins for different conditions
        controlBins = bins;
        leftModBins = bins & numModSteps(:,2)==1;
        rightModOneStepBins = bins & (numModSteps(:,3)==1);
        rightModTwoStepBins = bins & (numModSteps(:,3)==2);
        oneTwoRatio = sum(rightModOneStepBins) / (sum(rightModOneStepBins) + sum(rightModTwoStepBins)); % ratio of trials in which swing foot takes one large step to those in which an additional step is taken
        oneTwoRatio = oneTwoRatio * 2 - 1; % scale from -1 to 1
        

        % get left and right control locations
        controlLocations = {dataNew(controlBins).controlLocations};
        leftControlLocations = cellfun(@(x) x{2}, controlLocations, 'uniformoutput', 0);
        leftControlLocations = cat(1,leftControlLocations{:});
        rightControlLocations = cellfun(@(x) x{3}, controlLocations, 'uniformoutput', 0);
        rightControlLocations = cat(1,rightControlLocations{:});

        % get left modified locations
        leftModLocations = {dataNew(leftModBins).modifiedLocations};
        leftModLocations = cellfun(@(x) x{2}, leftModLocations, 'uniformoutput', 0);
        leftModLocations = cat(1,leftModLocations{:});

        % get right modified (one step) locations
        rightModOneStepLocations = {dataNew(rightModOneStepBins).modifiedLocations};
        rightModOneStepLocations = cellfun(@(x) x{3}, rightModOneStepLocations, 'uniformoutput', 0);
        rightModOneStepLocations = cat(1,rightModOneStepLocations{:});

        % get right modified (two step) locations
        rightModTwoStepLocations = {dataNew(rightModTwoStepBins).modifiedLocations};
        rightModTwoStepLocations = cellfun(@(x) x{3}(1,:,:), rightModTwoStepLocations, 'uniformoutput', 0);
        rightModTwoStepLocations = cat(1,rightModTwoStepLocations{:});


        % plot control left
        x = squeeze(leftControlLocations(:,1,:));
        x = x - mean(x(:,1)) + mean(squeeze(leftModLocations(:,1,1)));
        y = squeeze(leftControlLocations(:,2,:));
        plot(mean(x,1), mean(y,1), 'color', controlColor, 'linewidth', linWid); hold on;

        % plot control right
        x = squeeze(rightControlLocations(:,1,:));
        x = x - mean(x(:,1)) + mean(squeeze(rightModTwoStepLocations(:,1,1))); % !!! make average of one and two step?
        y = squeeze(rightControlLocations(:,2,:));
        plot(mean(x,1), mean(y,1), 'color', controlColor, 'linewidth', linWid); hold on;
        
        % plot mod left
        x = squeeze(leftModLocations(:,1,:));
        y = squeeze(leftModLocations(:,2,:));
        plot(mean(x,1), mean(y,1), 'color', colors(1,:), 'linewidth', linWid); hold on;

        % plot mod right, one step
        if ~isempty(rightModOneStepLocations)
            x = squeeze(rightModOneStepLocations(:,1,:));
            y = squeeze(rightModOneStepLocations(:,2,:));
            plot(mean(x,1), mean(y,1), 'color', colors(2,:), 'linewidth', linWid + oneTwoRatio*linWid); hold on;
        end

        % plot mod right, two step
        x = squeeze(rightModTwoStepLocations(:,1,:));
        y = squeeze(rightModTwoStepLocations(:,2,:));
        plot(mean(x,1), mean(y,1), 'color', colors(1,:), 'linewidth', linWid + -oneTwoRatio*linWid); hold on;
        
        % mark avg position of obsPos
        if any(rightModOneStepBins)
            oneStepObsPos = round(mean(obsPosIndInterps(rightModOneStepBins,3)));
            scatter(mean(squeeze(rightModOneStepLocations(:,1,oneStepObsPos))), ...
                mean(squeeze(rightModOneStepLocations(:,2,oneStepObsPos))), 200, colors(2,:), 'filled'); hold on
        end
        twoStepObsPos = round(nanmean(obsPosIndInterps(rightModTwoStepBins,3)));
        scatter(mean(squeeze(rightModTwoStepLocations(:,1,twoStepObsPos))), ...
                mean(squeeze(rightModTwoStepLocations(:,2,twoStepObsPos))), 200, colors(1,:), 'filled'); hold on



        % set appearance
        set(gca, 'dataaspectratio', [1 1 1], 'xlim', xLims, 'box', 'off', 'tickdir', 'out', 'ydir', 'reverse', ...
            'xtick', [], 'ytick', []);
        line([0 0], get(gca,'ylim'), 'color', [0 0 0], 'linewidth', 3)
        if g==speedBinNum; xlabel(['stance foot distance (m): ' phaseLabels{h}]); end
        if h==1; ylabel(['speed (m/s): ' speedLabels{g}]); end
    end
end

saveas(gcf, [getenv('OBSDATADIR') 'figures\meanKinematics.png']);

