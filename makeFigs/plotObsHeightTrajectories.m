% function plotObsHeightTrajectories(data)

% eventually fcn will be given bins, vector of zeros for trials not to be
% included, and 1..n for trials to be included in rows 1...n in the plot,
% allowing independent plots to be created for dft conditions, such as
% speed, light on vs off, etc.

% temp
bins = zeros(1,length(data));
bins([data.isLightOn]) = 1;
bins(~[data.isLightOn]) = 2;
binNames = {'light on', 'light off'};
numModSteps = reshape([data.modStepNum],4,length(data))';
isOneStepTrial = numModSteps(:,3)==1;



% settings
trialsPerPlot = 10;
maxSmpsThoughObs = 2;
heightBinNum = 3;
azimuth = 0;
elevation = 0;
obsDiam = 3.175; % (mm)
xLims = [-.05 0];
% xLims = [-.05 .04];
zLims = [0 .015];
lineWid = 2.5;

% initializations
plotNum = max(bins);
avgColors = parula(heightBinNum);
trialColors = parula(trialsPerPlot);
allInds = find(bins~=0 & isOneStepTrial' & ~[data.isFlipped]);
obsX = ([0 obsDiam]-.5*obsDiam) / 1000;

locationSmps = size(data(1).modifiedLocationsInterp{1},3);
locationsMod = nan(length(data), 3, locationSmps);
locationsControl = nan(length(data), 3, locationSmps);
throughObsBins = nan(length(data), locationSmps);

for i = allInds
    
    % find fore paw that is first to get over osbtacle
    isStepping = ~isnan(data(i).modifiedStepIdentities(:,:));
    lastModStepInds = table2array(rowfun(@(x)(find(x,1,'last')), table(isStepping')));
    [~, firstPawOverInd] = min(lastModStepInds .* [nan 1 1 nan]'); % mask out hind paws, inds 1 and 4
    
    locationsMod(i,:,:) = squeeze(data(i).modifiedLocationsInterp{firstPawOverInd}(end,:,:));
    locationsControl(i,:,:) = squeeze(data(i).controlLocationsInterp{firstPawOverInd}(end,:,:));
    
    % flip around midline if left paw is selected
    if firstPawOverInd==2 % fix this... it is not working
        locationsMod(i,2,:) = -locationsMod(i,2,:);
        locationsControl(i,2,:) = -locationsControl(i,2,:);
    end
    
    % check if trial passes through obstacle
    throughObsBins(i,:) = squeeze(locationsMod(i,3,:))<(data(i).obsHeightsVid/1000) & ... 
                          squeeze(locationsMod(i,1,:))>obsX(1) & ...
                          squeeze(locationsMod(i,1,:))<obsX(2);
    
    
    if sum(throughObsBins(i,:)) > maxSmpsThoughObs
%         fprintf('trial %i passed though obstacle\n', i);
        locationsMod(i,:,:) = nan;
        locationsControl(i,:,:) = nan;
    end
end
%
validTrials = (sum(throughObsBins,2)<=maxSmpsThoughObs)';
fprintf('%.2f of trials passed though obs\n', ...
    sum(sum(throughObsBins,2)>maxSmpsThoughObs) / size(throughObsBins,1));



close all;
figure('menubar', 'none', 'color', 'white', 'position', [-1500 500 1500 450], 'InvertHardCopy', 'off');

for i = 1:length(binNames)
    
    % PLOT INDIVIDUAL TRIALS
    subaxis(length(binNames),2,(i-1)*2+1)
    
    xs = squeeze(locationsMod(bins==i & validTrials,1,:));
%     ys = squeeze(locationsMod(bins==i & validTrials,2,:));
    zs = squeeze(locationsMod(bins==i & validTrials,3,:));
    
    trialInds = randperm(size(xs,1), trialsPerPlot);
    [~, sortInds] = sort([data(trialInds).obsHeightsVid]);
    trialInds = trialInds(sortInds);
    for j = 1:length(trialInds)
%         plot3(xs(trialInds(j),:), ys(trialInds(j),:), zs(trialInds(j),:), ...
%             'color', trialColors(j,:), 'linewidth', 1); hold on
        plot(xs(trialInds(j),:), zs(trialInds(j),:), ...
            'color', trialColors(j,:), 'linewidth', 1); hold on
    end
    
    % plot control mean
    xOffset = nanmean(squeeze(locationsMod(bins==i,1,1)));
    xsCtl = squeeze(locationsControl(bins==i,1,:));
%     ysCtl = squeeze(locationsControl(bins==i,2,:));
    zsCtl = squeeze(locationsControl(bins==i,3,:));
%     plot3(nanmean(xsCtl,1) - nanmean(xsCtl(:,1)) + xOffset, nanmean(ysCtl,1), nanmean(zsCtl,1), ...
%         'color', [0 0 0], 'linewidth', lineWid); hold on
    plot(nanmean(xsCtl,1) - nanmean(xsCtl(:,1)) + xOffset, nanmean(zsCtl,1), ...
        'color', [0 0 0], 'linewidth', lineWid); hold on
    
    % pimp appearance
    daspect([1 1 1]);
%     set(gca, 'xlim', xLims, 'zlim', zLims, 'view', [azimuth elevation], 'YDir', 'reverse');
    set(gca, 'xlim', xLims, 'zlim', zLims);
    axis off
    
    
    
    
    
    % PLOT HEIGHT BINNED MEANS
    subaxis(length(binNames),2,(i-1)*2+2)
%     plot3(nanmean(xsCtl,1) - nanmean(xsCtl(:,1)) + xOffset, nanmean(ysCtl,1), nanmean(zsCtl,1), ...
%         'color', [0 0 0], 'linewidth', lineWid); hold on
    plot(nanmean(xsCtl,1) - nanmean(xsCtl(:,1)) + xOffset, nanmean(zsCtl,1), ...
        'color', [0 0 0], 'linewidth', lineWid); hold on

    % get bins
    heights = [data.obsHeightsVid];
    binEdges = linspace(min(heights), max(heights), heightBinNum+1);
    heightBins = discretize(heights, binEdges);
    binLabels = cell(1,heightBinNum);
    for j = 1:heightBinNum; binLabels{j} = sprintf('%.1f', mean(heights(heightBins==j))); end
    
    for j = 1:heightBinNum
        xs = squeeze(locationsMod(bins==i & heightBins==j & validTrials,1,:));
%         ys = squeeze(locationsMod(bins==i & heightBins==j & validTrials,2,:));
        zs = squeeze(locationsMod(bins==i & heightBins==j & validTrials,3,:));
        
%         plot3(nanmean(xs,1), nanmean(ys,1), nanmean(zs,1), ...
%             'color', avgColors(j,:), 'linewidth', lineWid); hold on
        plot(nanmean(xs,1), nanmean(zs,1), ...
            'color', avgColors(j,:), 'linewidth', lineWid); hold on
        
        % add cylinder
        rad = obsDiam/1000/2;
        z = mean(heights(heightBins==j & bins==i)) / 1000;
%         circ = viscircles([0 z-obsDiam/1000/2], obsDiam/1000/2, 'color', avgColors(j,:));
        circ = rectangle('position', [0-rad, z-2*rad, 2*rad, 2*rad], ...
            'curvature', [1 1], 'facecolor', [avgColors(j,:) .8], 'edgecolor', 'none');
%         xyz1 = [0 yLims(1) z];
%         xyz2 = [0 yLims(2) z];
%         Cylinder(xyz1, xyz2, obsDiam/1000/2, 40, avgColors(j,:) , 1, 0);
    end
    
    
    
    
    
    
    
    % pimp appearance
    daspect([1 1 1]);
%     set(gca, 'xlim', xLims, 'zlim', zLims, 'zlim', zLims, 'view', [azimuth elevation], 'YDir', 'reverse');
    set(gca, 'xlim', xLims, 'zlim', zLims);
    axis off
end


% %% plot control trajectory
% xOffset = nanmean(squeeze(locationsMod(:,1,1)));
% xs = squeeze(locationsControl(allBins,1,:));
% ys = squeeze(locationsControl(allBins,2,:));
% zs = squeeze(locationsControl(allBins,3,:));
% plot3(nanmean(xs,1) - nanmean(xs(:,1)) + xOffset, nanmean(ys,1), nanmean(zs,1), ...
%     'color', [0 0 0], 'linewidth', lineWid); hold on
% 

% 
% for i = 1:length(conditionNames)+1
%     
%     subaxis(length(conditionNames)+1,1,i)
%     
%     daspect([1 1 1]);
%     set(gca, 'xlim', xLims, 'ylim', yLims, 'zlim', zLims, ...
%         'view', [azimuth elevation], 'YDir', 'reverse');
%     axis off
% 
% 
%     patch('Vertices', vertices, 'Faces', faces);
% 
%     % add lines for sides of wheel
%     line([xLims; xLims]', [yLims; yLims], zeros(2,2), ...
%         'color', [0 0 0], 'linewidth', 2) % x1
% end

% legend(gca, conditionNames, 'location', 'northeast');
% pos = get(leg, 'position');
% set(leg, 'position', [xLims(2)-pos(3) zLims(2)-pos(4) pos(3) pos(4)])

% blackenFig;
set(gcf, 'menubar', 'figure')
saveas(gcf, [getenv('OBSDATADIR') 'figures\obsHeightTrajectories.png']);
savefig([getenv('OBSDATADIR') 'figures\obsHeightTrajectories.fig'])



