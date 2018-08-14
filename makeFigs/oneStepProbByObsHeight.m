% convolve instead of bin // figure out what is going on with longer
% trials...

% settings
xVar = [data.swingStartDistance] + [data.predictedLengths]; % predicted distance to ob
validBins = ~[data.isLightOn];

heightBinNum = 5;
yLims = [0 1];
xLims = [-.03 .03];
dx = .001;
xWindowSize = .01;
heightLims = [3.175 10];
colors = copper(heightBinNum);



% initializations
xWindowSmps = ceil(xWindowSize/dx) - (mod(xWindowSize/dx,2)==0); % samples within an x window
windowShift = floor(xWindowSmps/2); % amount to shift the window left and right of center point

modStepNum = cellfun(@(x,ind) x(data(ind).firstModPaw), ...
    {data.modStepNum}, num2cell(1:length(data)));
obsHeights = [data.obsHeightsVid];

% !!! select subset of all trials here???

heightBinEdges = linspace(heightLims(1), heightLims(2), heightBinNum+1);
heightBins = discretize(obsHeights, heightBinEdges);
binLabels = cell(1,heightBinNum); for i = 1:heightBinNum; binLabels{i} = sprintf('%.2f', mean(obsHeights(heightBins==i))); end



xGrid = xLims(1):dx:xLims(2);
probs = nan(heightBinNum, length(xGrid));

for i = 1:heightBinNum
    for j = 1:length(xGrid)

        xBinLims = [xGrid(max(1,j-windowShift)) xGrid(min(length(xGrid),j+windowShift))];
        dataBins = xVar>=xBinLims(1) & xVar<=xBinLims(2) & heightBins==i & validBins;
        probs(i,j) = sum(modStepNum(dataBins)==1) / sum(dataBins);

    end
end


%%
close all;

figure('menubar', 'none', 'color', 'white');

for i = 1:heightBinNum
    plot(xGrid, probs(i,:), 'linewidth', 3, 'Color', colors(i,:)); hold on
end

set(gca, 'box', 'off')
xlabel('predicted distance to obs (m)')
ylabel('big step probability')
legend(binLabels, 'Location', 'northwest')



% saveas(gcf, [getenv('OBSDATADIR') 'figures\deltaLengthHeatMap.png']);









