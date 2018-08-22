function plotObsHeightTrajectories(data, bins, binNames, figTitle)


% settings
validTrials = ~[data.isWheelBreak];
trialsPerPlot = 100;
heightBinNum = 3;
heightLims = [3.175 10];
obsDiam = 3.175; % (mm)
xLims = [-.06 .04];
zLims = [0 .015];
lineWid = 2.5;
colors = [.25 1 .25; 1 .25 .25]; % start and stop colors of gradient




% initializations
avgColors = interp2(1:3, 1:2, colors, 1:3, linspace(1,2,heightBinNum)');

% exclude invalid trails, including those where max paw height is less than obs height
obsHgts = [data.obsHeightsVid];
pawHgts = cell2mat(cellfun(@(x,ind) max(squeeze(x{data(ind).firstPawOver}(end,3,:))), ...
          {data.modifiedLocations}, num2cell(1:length(data)), 'UniformOutput', false)) * 1000;
bins(obsHgts>pawHgts) = 0;
bins(~validTrials) = 0;
fprintf('%.3f of trials paw was not higher than obs..\n', ...
    sum(obsHgts>pawHgts)/length(data));
  

allInds = find(bins~=0);
locationSmps = size(data(1).modifiedLocationsInterp{1},3);
locationsMod = nan(length(data), 3, locationSmps);
locationsControl = nan(length(data), 3, locationSmps);
for i = allInds
    firstPawOverInd = data(i).firstPawOver;
    locationsMod(i,:,:) = squeeze(data(i).modifiedLocationsInterp{firstPawOverInd}(end,:,:));
    locationsControl(i,:,:) = squeeze(data(i).controlLocationsInterp{firstPawOverInd}(end,:,:));
end




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



