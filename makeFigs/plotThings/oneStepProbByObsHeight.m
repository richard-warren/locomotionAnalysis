function oneStepProbByObsHeight(data)

% settings
xVar = [data.swingStartDistance] + [data.predictedLengths]; % predicted distance to ob
validBins = ~[data.isLightOn] & ~[data.isWheelBreak];
% color = [237 9 127] / 255;
% colorFading = 0.5;

heightBinNum = 3;
yLims = [0 1];
xLims = [-.02 .02];
xRes = .0001;
heightLims = [3 11];
gausKernelSig = .004; % (m)


% initializations
% colors = repmat(color,heightBinNum,1) .* repmat(linspace(colorFading, 1, heightBinNum)',1,3);
colors = spring(heightBinNum);

modStepNum = cellfun(@(x,ind) x(data(ind).firstModPaw), ...
    {data.modStepNum}, num2cell(1:length(data)));
obsHeights = [data.obsHeightsVid];

kernel = arrayfun(@(x) (1/(gausKernelSig*sqrt(2*pi))) * exp(-.5*(x/gausKernelSig)^2), ...
    -gausKernelSig*5:xRes:gausKernelSig*5);
kernel = kernel / sum(kernel);
% xBinCenters = xLims(1) : xRes : xLims(2);
xBinCenters = xLims(1)-5*gausKernelSig : xRes : xLims(2)+5*gausKernelSig;
xBinEdges = [xBinCenters xBinCenters(end)+xRes] - .5*xRes;


heightBinEdges = linspace(heightLims(1), heightLims(2), heightBinNum+1);
heightBins = discretize(obsHeights, heightBinEdges);
binLabels = cell(1,heightBinNum); for i = 1:heightBinNum; binLabels{i} = sprintf('%.1f', mean(obsHeights(heightBins==i))); end



xGrid = xBinCenters(1):xRes:xBinCenters(end);
probs = nan(heightBinNum, length(xGrid));

for i = 1:heightBinNum
    
    conditionBins = heightBins==i & validBins;
    binCounts = histcounts(xVar(modStepNum==1 & conditionBins), xBinEdges) ./ ...
                histcounts(xVar(conditionBins), xBinEdges);
    binCounts = fillmissing(binCounts, 'linear');
    histoConv = conv(binCounts, kernel, 'same');
%     histoConv = histoConv / sum(conditionBins);
    
    probs(i,:) = histoConv;
end



figure('menubar', 'none', 'color', 'white', 'Position', [600 400 475 325], 'InvertHardcopy', 'off');

for i = 1:heightBinNum
    plot(xBinCenters, probs(i,:), 'linewidth', 3, 'Color', colors(i,:)); hold on
end

set(gca, 'box', 'off', 'XLim', xLims, 'YLim', yLims)
xlabel('predicted distance to obs (m)')
ylabel('big step probability')
legend(binLabels, 'Location', 'northwest')
blackenFig








