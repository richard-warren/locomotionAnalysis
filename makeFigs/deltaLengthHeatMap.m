function deltaLengthHeatMap(data, binVar)

% settings
% (x is binning variable, y is delta swing length)
xLims = [-.03 .015];
yLims = [-.03 .04];
xGridLims = [-.03 .02];
dX = .001;
xWindowSize = .008;
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
    dataInds = find(binVar>=xBinLims(1) & binVar<=xBinLims(2));
    deltaLengthsSub = deltaLengths(dataInds);
    
    binCounts = histogram(deltaLengthsSub, [yGrid yGrid(end)+dY] - .5*dY); % last argument changes bin centers to bin edges
    binCounts = binCounts.Values;
    
    histoConv = conv(binCounts, kernel, 'same');
    histoConv = histoConv / sum(histoConv);
    heatMap(:, i) = histoConv;
    
    probs(i) = sum(modStepNum(dataInds)==1) / length(dataInds);
    
end



close all;
figure('color', 'white', 'menubar', 'none', 'position', [1943 616 560 420], 'InvertHardcopy', 'off')
colormap hot
imagesc(xGrid*100, yGrid*100, heatMap)
line(get(gca, 'xlim'), [0 0], 'color', 'white', 'linewidth', 3)
line([0 0], get(gca, 'ylim'), 'color', 'white', 'linewidth', 3)

set(gca, 'ydir', 'normal', 'box', 'off', 'xlim', xLims*100, 'ylim', yLims*100, 'box', 'off', 'tickdir', 'out', ...
    'xtick', (xLims(1):.01:xLims(2))*100, 'ytick', (yLims(1):.01:yLims(2))*100)
xlabel('predicted distance to obs (cm)')
ylabel('\Deltax (cm)')

% yyaxis right
% plot(xGrid, probs, 'color', probColor, 'linewidth', 5)
% ylabel('probability of taking one big step')
% set(gca, 'ycolor', probColor, 'box', 'off')

saveas(gcf, [getenv('OBSDATADIR') 'figures\deltaLengthHeatMap.png']);




