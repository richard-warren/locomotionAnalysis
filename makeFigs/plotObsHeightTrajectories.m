function plotObsHeightTrajectories(data, bins, binNames, figTitle)

% eventually fcn will be given bins, vector of zeros for trials not to be
% included, and 1..n for trials to be included in rows 1...n in the plot,
% allowing independent plots to be created for dft conditions, such as
% speed, light on vs off, etc.


% settings
trialsPerPlot = 10;
maxSmpsThoughObs = 2;
heightBinNum = 3;
obsDiam = 3.175; % (mm)
xLims = [-.06 .04];
zLims = [0 .015];
lineWid = 2.5;

% initializations
avgColors = parula(heightBinNum);
trialColors = parula(trialsPerPlot);
allInds = find(bins~=0);
obsX = ([0 obsDiam]-.5*obsDiam) / 1000;

locationSmps = size(data(1).modifiedLocationsInterp{1},3);
locationsMod = nan(length(data), 3, locationSmps);
locationsControl = nan(length(data), 3, locationSmps);
throughObsBins = nan(length(data), locationSmps);

for i = allInds
    
    firstPawOverInd = data(i).firstPawOver;
    locationsMod(i,:,:) = squeeze(data(i).modifiedLocationsInterp{firstPawOverInd}(end,:,:));
    locationsControl(i,:,:) = squeeze(data(i).controlLocationsInterp{firstPawOverInd}(end,:,:));
    
    % check if trial passes through obstacle
    throughObsBins(i,:) = squeeze(locationsMod(i,3,:))<(data(i).obsHeightsVid/1000) & ... 
                          squeeze(locationsMod(i,1,:))>obsX(1) & ...
                          squeeze(locationsMod(i,1,:))<obsX(2);
    
    
    if sum(throughObsBins(i,:)) > maxSmpsThoughObs
        locationsMod(i,:,:) = nan;
        locationsControl(i,:,:) = nan;
    end
end
%
validTrials = (sum(throughObsBins,2)<=maxSmpsThoughObs)';
fprintf('%.2f of trials passed though obs\n', ...
    sum(sum(throughObsBins,2)>maxSmpsThoughObs) / size(throughObsBins,1));



figure('name', figTitle, 'menubar', 'none', 'color', 'white', 'position', [0 500 1500 450], 'InvertHardCopy', 'off');

for i = 1:length(binNames)
    
    % PLOT INDIVIDUAL TRIALS
    subaxis(length(binNames),2,(i-1)*2+1)
    
    xs = squeeze(locationsMod(bins==i & validTrials,1,:));
    zs = squeeze(locationsMod(bins==i & validTrials,3,:));
    
    trialInds = randperm(size(xs,1), min(trialsPerPlot, size(xs,1)));
    [~, sortInds] = sort([data(trialInds).obsHeightsVid]);
    trialInds = trialInds(sortInds);
    for j = 1:length(trialInds)
        plot(xs(trialInds(j),:), zs(trialInds(j),:), ...
            'color', trialColors(j,:), 'linewidth', 1); hold on
    end
    
    % plot control mean
    xOffset = nanmean(squeeze(locationsMod(bins==i,1,1)));
    xsCtl = squeeze(locationsControl(bins==i,1,:));
    zsCtl = squeeze(locationsControl(bins==i,3,:));
    plot(nanmean(xsCtl,1) - nanmean(xsCtl(:,1)) + xOffset, nanmean(zsCtl,1), ...
        'color', [0 0 0], 'linewidth', lineWid); hold on
    
    % pimp appearance
    daspect([1 1 1]);
    set(gca, 'xlim', xLims, 'zlim', zLims);
    axis off
    
    
    
    
    
    % PLOT HEIGHT BINNED MEANS
    subaxis(length(binNames),2,(i-1)*2+2)
    plot(nanmean(xsCtl,1) - nanmean(xsCtl(:,1)) + xOffset, nanmean(zsCtl,1), ...
        'color', [0 0 0], 'linewidth', lineWid); hold on
    title(binNames{i})

    % get bins
    heights = [data.obsHeightsVid];
    binEdges = linspace(min(heights), max(heights), heightBinNum+1);
    heightBins = discretize(heights, binEdges);
    binLabels = cell(1,heightBinNum);
    for j = 1:heightBinNum; binLabels{j} = sprintf('%.1f', mean(heights(heightBins==j))); end
    
    for j = 1:heightBinNum
        xs = squeeze(locationsMod(bins==i & heightBins==j & validTrials,1,:));
        zs = squeeze(locationsMod(bins==i & heightBins==j & validTrials,3,:));
        
        plot(nanmean(xs,1), nanmean(zs,1), ...
            'color', avgColors(j,:), 'linewidth', lineWid); hold on
        
        % add cylinder
        rad = obsDiam/1000/2;
        z = mean(heights(heightBins==j & bins==i)) / 1000;
        circ = rectangle('position', [0-rad, z-2*rad, 2*rad, 2*rad], ...
            'curvature', [1 1], 'facecolor', [avgColors(j,:) .8], 'edgecolor', 'none');
    end
    
    % pimp appearance
    daspect([1 1 1]);
    set(gca, 'xlim', xLims, 'ylim', zLims);
    axis off
end


% blackenFig;
set(gcf, 'menubar', 'figure')
saveas(gcf, [getenv('OBSDATADIR') 'figures\obsHeightTrajectories.png']);
savefig([getenv('OBSDATADIR') 'figures\obsHeightTrajectories.fig'])



