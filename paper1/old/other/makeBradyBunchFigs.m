% function makeBradyBunchFigs

% settings
session = '180125_001';
cols = 6;
aspectRatio = 9/16; % height / width
yLims = [60 120]; % used to crop image
xLims = [50 360]; % used to crop image
contrastLims = [.2 1]; % pixels at these proportional values are mapped to 0 and 255
cmapRes = 1000;
gausKernelSig = .02; % (m)
velRes = .001; % (m)
velMin = .15;



% initializations
cmap = jet(cmapRes);
dims = [diff(yLims)+1, diff(xLims)+1]; % height, width
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'frameTimeStamps');
load([getenv('OBSDATADIR') 'sessions\' session '\wiskContactData.mat'], 'contactTimes');
load([getenv('OBSDATADIR') 'sessions\' session '\tracking\velocityInfo.mat'], 'trialVels');

validInds = trialVels>velMin;
trialVels = trialVels(validInds);
contactTimes = contactTimes(validInds);
trialVels(end) = nan; % !!! this is a hack because last trial always read zero vel for some reason... wtf


width = dims(2)*cols;
rows = floor( (width*aspectRatio) / dims(1));
trials = randperm(length(contactTimes)-1, rows*cols);
minVel = min(trialVels(trials));
maxVel = max(trialVels(trials));




% make histogram
close all;
figure('color', 'white', 'menubar', 'none', 'position', [100 100 dims(2)*2 dims(1)*4]*.75, 'InvertHardcopy', 'off')
kernel = arrayfun(@(x) (1/(gausKernelSig*sqrt(2*pi))) * exp(-.5*(x/gausKernelSig)^2), ...
    -gausKernelSig*5:velRes:gausKernelSig*5);
kernel = kernel / sum(kernel);

velBinEdges = velMin : velRes : max(trialVels)+3*gausKernelSig;
velBinCenters = velBinEdges(1:end-1) + .5*velRes;
binCounts = histogram(trialVels, velBinEdges);
binCounts = binCounts.Values;

histoConv = conv(binCounts, kernel, 'same');
histoConv = histoConv / sum(histoConv);
histoPlot = plot(velBinCenters, histoConv, 'linewidth', 5);

lineCmap = [uint8(jet(length(histoConv))*255) uint8(ones(length(histoConv),1))].'; pause(.01)
set(histoPlot.Edge, 'ColorBinding', 'interpolated', 'ColorData', lineCmap)
ax = gca;
set(ax, 'box', 'off', 'xlim', [velBinCenters(1) velBinCenters(end)], 'position', [0.05 .25 .9 .75], 'xcolor', [0 0 0]);
ax.YAxis.Visible = 'off';
xlabel('velocity (m/s)')
blackenFig;
savefig([getenv('OBSDATADIR') 'figures\bradyBunchVelHisto.fig'])
saveas(gcf, [getenv('OBSDATADIR') 'figures\bradyBunchVelHisto.png'])
print('-clipboard', '-dmeta')
%%



collage = uint8(nan(dims(1)*rows, dims(2)*cols));
collageColored = uint8(nan(dims(1)*rows, dims(2)*cols, 3));

for i = 1:length(trials)
    
    % get trial frame
    trialInd = knnsearch(frameTimeStamps, contactTimes(trials(i)));
    frame = rgb2gray(read(vid, trialInd));
    frame = frame(yLims(1):yLims(2), xLims(1):xLims(2));
    frame = imadjust(frame, contrastLims, [0 1]);
    
    % color image according to speed
    colorInd = round(interp1([velBinCenters(1) velBinCenters(end)], [1 cmapRes],  trialVels(trials(i))));
    color = cmap(colorInd,:);
    frameColored = cat(3, frame*color(1), frame*color(2), frame*color(3));
    
    % get collage inds
    [row, col] = ind2sub([rows cols], i);
    rowInd = (row-1)*dims(1)+1;
    colInd = (col-1)*dims(2)+1;
    
    collage(rowInd:rowInd+dims(1)-1, colInd:colInd+dims(2)-1) = frame;
    collageColored(rowInd:rowInd+dims(1)-1, colInd:colInd+dims(2)-1, :) = frameColored;
    
end

imwrite(collage, [getenv('OBSDATADIR') 'figures\bradyBunch.png'])
imwrite(collageColored, [getenv('OBSDATADIR') 'figures\bradyBunchColored.png'])










