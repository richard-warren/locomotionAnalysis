

% to do: change colors on sperm // select sperm in proportion to strategy // merge interp and sperm? // fix green dots on interp plot

% settings
sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003', ...
            '180125_001', '180125_002', '180125_003'};


% initializations
data = getKinematicData(sessions);
tic; save([getenv('OBSDATADIR') 'kinematicData.mat'], 'data'); toc;
data = data([data.oneSwingOneStance]);

%%
binNum = 5;

% get prediceted distance to obs bins
predictedDistToObs = [data.swingStartDistance] + [data.predictedLengths];
binEdges = prctile(predictedDistToObs, linspace(0,100,binNum+1));
bins = discretize(predictedDistToObs, binEdges);

% create speed and phase labels
binLabels = cell(1,binNum);
for i = 1:binNum; binLabels{i} = sprintf('%.3f', mean(predictedDistToObs(bins==i))); end


%% sperm plots

% settings
yGridLims = [-.1 .1];
tracesPerPlot = 15;

figure('color', 'white', 'menubar', 'none', 'position', [400 200 250*binNum 800]);

for h = 1:binNum

    ax = subaxis(1, binNum , h, 'spacing', .02, 'padding', .04, 'margin', .02);

    binInds = find(bins==h);
    plotTraces = min(tracesPerPlot, length(binInds));
    binIndsInds = randperm(length(binInds), tracesPerPlot);
    dataInds = binInds(binIndsInds);

    for i = dataInds
        for j = 2:3

            realInds = ~isnan(data(i).modifiedStepIdentities(:,j));
            steps = unique(data(i).modifiedStepIdentities(realInds,j));
            colors = winter(length(steps));

            for k = steps'

                % plot x and y trajectories
                trialInds = data(i).modifiedStepIdentities(:,j)==k;
                x = data(i).locations(trialInds,1,j);
                y = data(i).locations(trialInds,2,j);
                plot(y, x, 'color', colors(k,:)); hold on

                % scatter dots at start of each swing
                scatter(y(end), x(end), 100, colors(k,:), 'filled'); hold on

                % scatter position of swing foot at obsPos
                if j==3
                    scatter(data(i).locations(data(i).obsPosInd,2,j), data(i).locations(data(i).obsPosInd,1,j), ...
                        100, [0 0 0], 'x'); hold on
                end
            end
        end
    end

    % set appearance
    set(gca, 'dataaspectratio', [1 1 1], 'ylim', yGridLims, 'box', 'off', 'tickdir', 'out', 'xtick', [], 'ytick', []);
    line(get(gca,'xlim'), [0 0], 'color', [0 0 0], 'linewidth', 3)
    xlabel(['predicted dist to obs (m): ' binLabels{h}]);
end

saveas(gcf, [getenv('OBSDATADIR') 'figures\trialKinematics.png']);


%% average interpolated trajectories


% settings
controlColor = [.65 .65 .65];
yGridLims = [-.09 .05];
linWid = 4;
colors = winter(2);

% initializations
figure('color', 'white', 'menubar', 'none', 'position', [400 200 250*binNum 800]);
numModSteps = reshape([data.modStepNum],4,length(data))';
obsPosIndInterps = reshape([data.pawObsPosIndInterp],4,length(data))';


for h = 1:binNum

    % get subplot bins
    ax = subaxis(1, binNum , h, 'spacing', .02, 'padding', .04, 'margin', .02);
    binBins = (bins==h)';

    % get subplot bins for different conditions
    controlBins = binBins;
    leftModBins = binBins & numModSteps(:,2)==1;
    rightModOneStepBins = binBins & (numModSteps(:,3)==1);
    rightModTwoStepBins = binBins & (numModSteps(:,3)==2);
    oneTwoRatio = sum(rightModOneStepBins) / (sum(rightModOneStepBins) + sum(rightModTwoStepBins)); % ratio of trials in which swing foot takes one large step to those in which an additional step is taken
    oneTwoRatio = oneTwoRatio * 2 - 1; % scale from -1 to 1


    % get left and right control locations
    controlLocations = {data(controlBins).controlLocationsInterp};
    leftControlLocations = cellfun(@(x) x{2}, controlLocations, 'uniformoutput', 0);
    leftControlLocations = cat(1,leftControlLocations{:});
    rightControlLocations = cellfun(@(x) x{3}, controlLocations, 'uniformoutput', 0);
    rightControlLocations = cat(1,rightControlLocations{:});

    % get left modified locations
    leftModLocations = {data(leftModBins).modifiedLocationsInterp};
    leftModLocations = cellfun(@(x) x{2}, leftModLocations, 'uniformoutput', 0);
    leftModLocations = cat(1,leftModLocations{:});

    % get right modified (one step) locations
    rightModOneStepLocations = {data(rightModOneStepBins).modifiedLocationsInterp};
    rightModOneStepLocations = cellfun(@(x) x{3}, rightModOneStepLocations, 'uniformoutput', 0);
    rightModOneStepLocations = cat(1,rightModOneStepLocations{:});

    % get right modified (two step) locations
    rightModTwoStepLocations = {data(rightModTwoStepBins).modifiedLocationsInterp};
    rightModTwoStepLocations = cellfun(@(x) x{3}(1,:,:), rightModTwoStepLocations, 'uniformoutput', 0);
    rightModTwoStepLocations = cat(1,rightModTwoStepLocations{:});


    % plot control left
    x = squeeze(leftControlLocations(:,1,:));
    x = x - mean(x(:,1)) + mean(squeeze(leftModLocations(:,1,1)));
    y = squeeze(leftControlLocations(:,2,:));
    plot(mean(y,1), mean(x,1), 'color', controlColor, 'linewidth', linWid); hold on;

    % plot control right
    x = squeeze(rightControlLocations(:,1,:));
    x = x - mean(x(:,1)) + mean(squeeze(rightModTwoStepLocations(:,1,1))); % !!! make average of one and two step?
    y = squeeze(rightControlLocations(:,2,:));
    plot(mean(y,1), mean(x,1), 'color', controlColor, 'linewidth', linWid); hold on;

    % plot mod left
    x = squeeze(leftModLocations(:,1,:));
    y = squeeze(leftModLocations(:,2,:));
    plot(mean(y,1), mean(x,1), 'color', colors(1,:), 'linewidth', linWid); hold on;

    % plot mod right, one step
    if ~isempty(rightModOneStepLocations)
        x = squeeze(rightModOneStepLocations(:,1,:));
        y = squeeze(rightModOneStepLocations(:,2,:));
        plot(mean(y,1), mean(x,1), 'color', colors(2,:), 'linewidth', linWid + oneTwoRatio*linWid); hold on;
    end

    % plot mod right, two step
    x = squeeze(rightModTwoStepLocations(:,1,:));
    y = squeeze(rightModTwoStepLocations(:,2,:));
    plot(mean(y,1), mean(x,1), 'color', colors(1,:), 'linewidth', linWid + -oneTwoRatio*linWid); hold on;

    % mark avg position of obsPos
    if any(rightModOneStepBins)
        oneStepObsPos = round(mean(obsPosIndInterps(rightModOneStepBins,3)));
        scatter(mean(squeeze(rightModOneStepLocations(:,2,oneStepObsPos))), ...
                mean(squeeze(rightModOneStepLocations(:,1,oneStepObsPos))), 200, colors(2,:), 'filled'); hold on
    end
    twoStepObsPos = round(nanmean(obsPosIndInterps(rightModTwoStepBins,3)));
    scatter(mean(squeeze(rightModTwoStepLocations(:,2,twoStepObsPos))), ...
            mean(squeeze(rightModTwoStepLocations(:,1,twoStepObsPos))), 200, colors(1,:), 'filled'); hold on



    % set appearance
    set(gca, 'dataaspectratio', [1 1 1], 'ylim', yGridLims, 'xlim', [-.02 .02], ...
        'box', 'off', 'tickdir', 'out', 'xtick', [], 'ytick', []);
    line(get(gca,'xlim'), [0 0], 'color', [0 0 0], 'linewidth', 3)
    xlabel(['predicted dist to obs (m): ' binLabels{h}]);
end

saveas(gcf, [getenv('OBSDATADIR') 'figures\meanKinematics.png']);

%% histograms

% settings
xGridLims = [.02 .12];
yGridLims = [0 .3];
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
    set(ax, 'box', 'off', 'xlim', xGridLims, 'ylim', yGridLims, 'tickdir', 'out')
    set(ax, 'ytick', [], 'ylabel', []); ax.YAxis.Visible = 'off';
end

legend('modified swing lengths (lengthened)', 'modified swing lengths (shortened)', 'control swing lengths')

saveas(gcf, [getenv('OBSDATADIR') 'figures\swingLengthHistograms.png']);





%% heat map

% settings
% (x is predicted position of paw relative to obs, y is swing length)

xLims = [-.03 .015];
yLims = [-.03 .04];
xGridLims = [-.03 .02];
dX = .001;
xWindowSize = .01;
yGridLims = [-.03 .08];
dY = .001;
yKernelSig = .01;
probColor = [0 .7 1];


% initializations
xWindowSmps = ceil(xWindowSize/dX) - (mod(xWindowSize/dX,2)==0); % round to nearest odd number
deltaLengths = cellfun(@(x) x(1,3), {data.modifiedSwingLengths}) - [data.predictedLengths];
modStepNum = cellfun(@(x) x(1,3), {data.modStepNum});
windowShift = floor(xWindowSmps/2);
xGrid = xGridLims(1):dX:xGridLims(2);
yGrid = yGridLims(1):dY:yGridLims(2);
kernel = arrayfun(@(x) (1/(yKernelSig*sqrt(2*pi))) * exp(-.5*(x/yKernelSig)^2), ...
    -yKernelSig*5:dY:yKernelSig*5);
kernel = kernel / sum(kernel);

heatMap = nan(length(yGrid), length(xGrid));
probs = nan(1, length(xGrid));

for i = 1:length(xGrid)
    
    % get data within bin
    xBinLims = [xGrid(max(1,i-windowShift)) xGrid(min(length(xGrid),i+windowShift))];
    dataInds = find(predictedDistToObs>=xBinLims(1) & predictedDistToObs<=xBinLims(2));
    deltaLengthsSub = deltaLengths(dataInds);
    
    binCounts = histogram(deltaLengthsSub, [yGrid yGrid(end)+dY] - .5*dY); % last argument changes bin centers to bin edges
    binCounts = binCounts.Values;
    
    histoConv = conv(binCounts, kernel, 'same');
    histoConv = histoConv / sum(histoConv);
    heatMap(:, i) = histoConv;
    
    probs(i) = sum(modStepNum(dataInds)==1) / length(dataInds);
    
end



close all;
figure('color', 'white', 'menubar', 'none', 'position', [1943 616 560 420])
colormap hot
imagesc(xGrid, yGrid, heatMap)
line(get(gca, 'xlim'), [0 0], 'color', 'white', 'linewidth', 3, 'linestyle', ':')
line([0 0], get(gca, 'ylim'), 'color', 'white', 'linewidth', 3, 'linestyle', ':')

set(gca, 'ydir', 'normal', 'box', 'off', 'xlim', xLims, 'ylim', yLims)
xlabel('predicted distance to obs (m)')
ylabel('\Deltax (m)')

yyaxis right
plot(xGrid, probs, 'color', probColor, 'linewidth', 5)
ylabel('probability of taking one big step')
set(gca, 'ycolor', probColor)

saveas(gcf, [getenv('OBSDATADIR') 'figures\deltaLengthHeatMap.png']);





