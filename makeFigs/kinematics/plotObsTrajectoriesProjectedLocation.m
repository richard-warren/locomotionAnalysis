

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
load([getenv('OBSDATADIR') 'kinematicData.mat'], 'data')
data = data([data.oneSwingOneStance]);
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
yLims = [-.1 .1];
xLims = [-.02 .02];
tracesPerPlot = 15;
colors = winter(2)*.85;
controlColor = [.65 .65 .65];
linWid = 4;
circSize = 150;
scaleBarSize = .01;

% initializations
modStepNum = cellfun(@(x) x(1,3), {data.modStepNum});
close all;
figure('color', 'white', 'menubar', 'none', 'position', [400 200 250*binNum 700]);

for h = 1:binNum

    ax = subaxis(1, binNum , h, 'spacing', .02, 'padding', .04, 'margin', .02);

    
    % get inds for one and two step trials within bin
    % (making sure the number of trials of each type reflects the proportion of each type across all trials)
    binBins = (bins==h);
    oneStepPortion = sum(modStepNum(binBins)==1) / sum(binBins);
    oneTwoStepTrialNums = [round(tracesPerPlot*oneStepPortion) round(tracesPerPlot*(1-oneStepPortion))];
    
    if oneTwoStepTrialNums(1)>0
        oneStepInds = find(binBins & modStepNum==1);
        oneStepInds = oneStepInds(randperm(length(oneStepInds), oneTwoStepTrialNums(1)));
    else
        oneStepInds = [];
    end
    
    twoStepInds = find(binBins & modStepNum>1);
    twoStepInds = twoStepInds(randperm(length(twoStepInds), oneTwoStepTrialNums(2)));
    allInds = [oneStepInds twoStepInds];
    colorInds = [ones(1,length(oneStepInds)) ones(1,length(twoStepInds))*2];
    
    
    % plot individual trials
    for i = [oneStepInds twoStepInds]
        for j = 2:3

            % plot x and y trajectories
            locationInds = data(i).modifiedStepIdentities(:,j)==1;
            x = data(i).locations(locationInds,1,j);
            y = data(i).locations(locationInds,2,j);
            if data(i).modStepNum(1,j)~=1; colorInd=1; else; colorInd=2; end
            plot(y, x, 'color', colors(colorInd,:)); hold on

            % scatter dots at start of each swing
            scatter(y(end), x(end), 100, colors(colorInd,:), 'filled'); hold on

            % scatter position of swing foot at obsPos
            if j==3
                scatter(data(i).locations(data(i).obsPosInd,2,j), data(i).locations(data(i).obsPosInd,1,j), ...
                    100, [0 0 0], 'x'); hold on
            end
        end
    end
    
    
    % get right control locations
    controlLocations = cellfun(@(x) x{3}(end,:,:), ...
        {data(binBins).controlLocationsInterp}, 'uniformoutput', 0); % only take last control step
    controlLocations = cat(1,controlLocations{:});
    
    % get x offset (mean starting x pos of modified steps)
    modLocations = cellfun(@(x) x{3}(1,:,:), {data(binBins).modifiedLocationsInterp}, 'uniformoutput', 0);
    modLocations = cat(1,modLocations{:});
    xOffset = mean(squeeze(modLocations(:,1,1)));
    
    x = squeeze(controlLocations(:,1,:));
    x = x - (mean(x(:,1)) - xOffset);
    y = squeeze(controlLocations(:,2,:));
    plot(mean(y,1), mean(x,1), 'color', controlColor, 'linewidth', linWid); hold on;
    
    
    

    % set appearance
    set(gca, 'dataaspectratio', [1 1 1], 'xlim', [-.02 .02], 'ylim', yLims);
    axis off
    line(get(gca,'xlim'), [0 0], 'color', [0 0 0], 'linewidth', 3)
    xlabel(['predicted dist to obs (m): ' binLabels{h}]);
end

% add scale bar
line([xLims(2)-scaleBarSize xLims(2)], repmat(yLims(1),1,2), 'linewidth', 3, 'color', 'black')
text(xLims(2)-.5*scaleBarSize, yLims(1)+.005, sprintf('%i mm', scaleBarSize*1000), 'horizontalalignment', 'center')

saveas(gcf, [getenv('OBSDATADIR') 'figures\trialKinematics.png']);




% AVERAGE INTERPOLATE TRAJECTORIES

% initializations
figure('color', 'white', 'menubar', 'none', 'position', [400 200 250*binNum 700]);
numModSteps = reshape([data.modStepNum],4,length(data))';
obsPosIndInterps = cellfun(@(x) x(1,3), {data.pawObsPosIndInterp});


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
    leftControlLocations = cellfun(@(x) x{2}(end,:,:), controlLocations, 'uniformoutput', 0);
    leftControlLocations = cat(1,leftControlLocations{:});
    rightControlLocations = cellfun(@(x) x{3}(end,:,:), controlLocations, 'uniformoutput', 0);
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
    
    % get left and right x offsets
    modLocations = cellfun(@(x) x{2}(1,:,:), {data(binBins).modifiedLocationsInterp}, 'uniformoutput', 0);
    modLocations = cat(1,modLocations{:});
    xOffsetLeft = mean(squeeze(modLocations(:,1,1)));
    modLocations = cellfun(@(x) x{3}(1,:,:), {data(binBins).modifiedLocationsInterp}, 'uniformoutput', 0);
    modLocations = cat(1,modLocations{:});
    xOffsetRight = mean(squeeze(modLocations(:,1,1)));


    % plot control left
    x = squeeze(leftControlLocations(:,1,:));
    x = x - (mean(x(:,1)) - xOffsetLeft);
    y = squeeze(leftControlLocations(:,2,:));
    plot(mean(y,1), mean(x,1), 'color', controlColor, 'linewidth', linWid); hold on;

    % plot control right
    x = squeeze(rightControlLocations(:,1,:));
    x = x - (mean(x(:,1)) - xOffsetRight);
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
        oneStepObsPos = round(mean(obsPosIndInterps(rightModOneStepBins)));
        scatter(mean(squeeze(rightModOneStepLocations(:,2,oneStepObsPos))), ...
                mean(squeeze(rightModOneStepLocations(:,1,oneStepObsPos))), circSize + oneTwoRatio*circSize, colors(2,:), 'filled'); hold on
    end
    twoStepObsPos = round(nanmean(obsPosIndInterps(rightModTwoStepBins)));
    scatter(mean(squeeze(rightModTwoStepLocations(:,2,twoStepObsPos))), ...
            mean(squeeze(rightModTwoStepLocations(:,1,twoStepObsPos))), circSize + -oneTwoRatio*circSize, colors(1,:), 'filled'); hold on



    % set appearance
    set(gca, 'dataaspectratio', [1 1 1], 'ylim', yLims, 'xlim', xLims);
    axis off
    line(get(gca,'xlim'), [0 0], 'color', [0 0 0], 'linewidth', 3)
    xlabel(['predicted dist to obs (m): ' binLabels{h}]);
end

% add scale bar
line([xLims(2)-scaleBarSize xLims(2)], repmat(yLims(1),1,2), 'linewidth', 3, 'color', 'black')
text(xLims(2)-.5*scaleBarSize, yLims(1)+.005, sprintf('%i mm', scaleBarSize*1000), 'horizontalalignment', 'center')

saveas(gcf, [getenv('OBSDATADIR') 'figures\meanKinematics.png']);

%% histograms

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





%% heat map

% settings
% (x is predicted position of paw relative to obs, y is swing length)

xLims = [-.03 .015];
yLims = [-.03 .04];
xGridLims = [-.03 .02];
dX = .001;
xWindowSize = .008;
yLims = [-.03 .08];
dY = .001;
yKernelSig = .008;
probColor = [0 .7 1];


% initializations
xWindowSmps = ceil(xWindowSize/dX) - (mod(xWindowSize/dX,2)==0); % round to nearest odd number
deltaLengths = cellfun(@(x) x(1,3), {data.modifiedSwingLengths}) - [data.predictedLengths];
modStepNum = cellfun(@(x) x(1,3), {data.modStepNum});
windowShift = floor(xWindowSmps/2);
xGrid = xGridLims(1):dX:xGridLims(2);
yGrid = yLims(1):dY:yLims(2);
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





