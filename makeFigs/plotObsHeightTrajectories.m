function plotObsHeightTrajectories(data, bins, binNames, figTitle)



% settings
validTrials = ~[data.isWheelBreak];
trialsPerPlot = 100;
maxSmpsThoughObs = 2;
heightBinNum = 3;
heightLims = [3.175 10];
obsDiam = 3.175; % (mm)
xLims = [-.06 .04];
zLims = [0 .015];
lineWid = 2.5;
colors = [.25 1 .25; 1 .25 .25]; % start and stop colors of gradient
% colors = [0 0 1; 0 1 0];



% initializations
avgColors = interp2(1:3, 1:2, colors, 1:3, linspace(1,2,heightBinNum)');
trialColors = nan(trialsPerPlot, 3); for i = 1:3; trialColors(:,i) = linspace(colors(1,i), colors(2,i), trialsPerPlot)'; end
allInds = find(bins~=0);
obsX = [-.5*obsDiam .5*obsDiam] / 1000;

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

validTrials = validTrials & (sum(throughObsBins,2)<=maxSmpsThoughObs)';
bins(~validTrials) = 0;
fprintf('%.2f of trials passed though obs\n', ...
    sum(sum(throughObsBins,2)>maxSmpsThoughObs) / size(throughObsBins,1));



figure('name', figTitle, 'menubar', 'none', 'color', 'white', 'position', [100 100 1500 150*length(binNames)], 'InvertHardCopy', 'off');

for i = 1:length(binNames)
    
    % PLOT INDIVIDUAL TRIALS
    subaxis(length(binNames),2,(i-1)*2+1)
    
    xs = squeeze(locationsMod(bins==i,1,:));
    zs = squeeze(locationsMod(bins==i,3,:));
    
    trialInds = randperm(size(xs,1), min(trialsPerPlot, size(xs,1)));
    for j = 1:length(trialInds)
        hgt = data(trialInds(j)).obsHeightsVid; hgt = min(hgt, heightLims(2)); hgt = max(hgt, heightLims(1));
        c = interp2(heightLims, 1:3, colors', hgt, 1:3);
        plot(xs(trialInds(j),:), zs(trialInds(j),:), ...
            'color', [c' 0.4], 'linewidth', 1); hold on
    end
    
    % plot control mean
    xOffset = nanmean(squeeze(locationsMod(bins==i,1,1)));
    xsCtl = squeeze(locationsControl(bins==i,1,:));
    zsCtl = squeeze(locationsControl(bins==i,3,:));
    plot(nanmean(xsCtl,1) - nanmean(xsCtl(:,1)) + xOffset, nanmean(zsCtl,1), ...
        'color', [0 0 0], 'linewidth', lineWid); hold on
    
    % add floor
    line(xLims, [0 0], 'linewidth', 1, 'color', [0 0 0]);
    
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
    binEdges = linspace(heightLims(1), heightLims(2), heightBinNum+1);
    heightBins = discretize(heights, binEdges);
    binLabels = cell(1,heightBinNum);
    for j = 1:heightBinNum; binLabels{j} = sprintf('%.1f', mean(heights(heightBins==j))); end
    
    for j = 1:heightBinNum
        xs = squeeze(locationsMod(bins==i & heightBins==j,1,:));
        zs = squeeze(locationsMod(bins==i & heightBins==j,3,:));
        
        plot(nanmean(xs,1), nanmean(zs,1), ...
            'color', avgColors(j,:), 'linewidth', lineWid); hold on
        
        % add cylinder
        rad = obsDiam/1000/2;
        z = mean(heights(heightBins==j & bins==i)) / 1000;
        circ = rectangle('position', [0-rad, z-2*rad, 2*rad, 2*rad], ...
            'curvature', [1 1], 'facecolor', [avgColors(j,:) 1], 'edgecolor', 'none');
    end
    
    % add floor
    line(xLims, [0 0], 'linewidth', 1, 'color', [0 0 0]);
    
    % pimp appearance
    daspect([1 1 1]);
    set(gca, 'xlim', xLims, 'ylim', zLims);
    axis off
end


blackenFig;



