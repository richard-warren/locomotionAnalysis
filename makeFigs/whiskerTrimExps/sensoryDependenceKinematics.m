% function sensoryDependenceKinematics(data)

% temp
figTitle = 'temp';


% settings
validTrials = ~[data.isWheelBreak];
beforeObsOnly = true; % if beforeObsOnly is true, only take kinematics up until paw is near obs in x dimension, then interpolate these trajectories of common grid
heightBinNum = 3;
heightLims = [3.175 10];
obsDiam = 3.175; % (mm)
xLims = [-.06 .04];
zLims = [0 .015];
lineWid = 2.5;
colors = [.25 1 .25; 1 .25 .25]; % start and stop colors of gradient
% colors = [0 0 1; 0 1 0];
preObsLim = .005; % when before


% initialiations
if beforeObsOnly; xLims(2) = 0; end
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'whiskerTrimNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);
wiskSessions = unique(sessionInfo.session(strcmp(sessionInfo.preOrPost, 'pre')));
noWiskSessions = unique(sessionInfo.session(strcmp(sessionInfo.preOrPost, 'post')));

wiskLightBins = ismember({data.session}, wiskSessions) & [data.isLightOn] & validBins;
wiskBins = ismember({data.session}, wiskSessions) & ~[data.isLightOn] & validBins;
lightBins = ismember({data.session}, noWiskSessions) & [data.isLightOn] & validBins;
neitherBins = ismember({data.session}, noWiskSessions) & ~[data.isLightOn] & validBins;

conditionBins = {wiskLightBins, wiskBins, lightBins, neitherBins};
binNames = {'wisk+light', 'wisk', 'light', 'neither'};
bins = zeros(1,length(data));
for i = 1:4; bins(conditionBins{i}) = i; end


% exclude invalid trails, including those where max paw height is less than obs height
obsHgts = [data.obsHeightsVid];
pawHgts = cell2mat(cellfun(@(x,ind) max(squeeze(x{data(ind).firstPawOver}(end,3,:))), ...
          {data.modifiedLocations}, num2cell(1:length(data)), 'UniformOutput', false)) * 1000;
bins(obsHgts>pawHgts) = 0;
bins(~validTrials) = 0;
fprintf('%.3f of trials paw was not higher than obs..\n', ...
    sum(obsHgts>pawHgts)/length(data));


% initializations
avgColors = interp2(1:3, 1:2, colors, 1:3, linspace(1,2,heightBinNum)');
allInds = find(bins~=0);
obsX = [-.5*obsDiam .5*obsDiam] / 1000;
if beforeObsOnly; xLims(2) = 0; end

locationSmps = size(data(1).modifiedLocationsInterp{1},3);
locationsMod = nan(length(data), 3, locationSmps);
locationsControl = nan(length(data), 3, locationSmps);

for i = allInds
    
    firstPawOverInd = data(i).firstPawOver;
    
    if ~beforeObsOnly
        locationsMod(i,:,:) = squeeze(data(i).modifiedLocationsInterp{firstPawOverInd}(end,:,:));
        locationsControl(i,:,:) = squeeze(data(i).controlLocationsInterp{firstPawOverInd}(end,:,:));
    else
        % get locations up until paw is within preObsLim of obs
        mod = squeeze(data(i).modifiedLocations{firstPawOverInd}(end,:,:));
        control = squeeze(data(i).controlLocations{firstPawOverInd}(end,:,:));
        lastInd = find(mod(1,:)>-preObsLim,1,'first');
        
        % interpolate over common grid
        if lastInd>1
            modInterp = interp2(1:lastInd, [1:3]', mod(:,1:lastInd), linspace(1,lastInd,locationSmps), [1:3]');
            controlInterp = interp2(1:lastInd, [1:3]', control(:,1:lastInd), linspace(1,lastInd,locationSmps), [1:3]');
            locationsMod(i,:,:) = modInterp;
            locationsControl(i,:,:) = controlInterp;
        end
    end
    
    
end



% plot everyting
close all; figure('name', figTitle, 'menubar', 'none', 'color', 'white', 'position', [100 100 750 900], 'InvertHardCopy', 'off');

for i = 1:length(binNames)
    
    subaxis(length(binNames)+1,1,i)
    title(binNames{i})
    
    
    % plot control mean
    xOffset = nanmean(squeeze(locationsMod(bins==i,1,1)));
    xsCtl = squeeze(locationsControl(bins==i,1,:));
    zsCtl = squeeze(locationsControl(bins==i,3,:));
    plot(nanmean(xsCtl,1) - nanmean(xsCtl(:,1)) + xOffset, nanmean(zsCtl,1), ...
        'color', [0 0 0], 'linewidth', lineWid); hold on
    

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




