function sensoryDependenceKinematics(data)

% settings
validBins = ~[data.isWheelBreak];
beforeObsOnly = true; % if beforeObsOnly is true, only take kinematics up until paw is near obs in x dimension, then interpolate these trajectories of common grid
heightBinNum = 3;
heightLims = [3.175 10];
obsDiam = 3.175; % (mm)
xLims = [-.05 .04];
zLims = [0 .015];
lineWid = 2.5;
conditionColors = hsv(4)+.2; conditionColors(conditionColors>1)=1;
preObsLim = .005;
colorFading = 0.5;


% initialiations
obsRadius = obsDiam/1000/2;
if beforeObsOnly; xLims(2) = 0; end
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'whiskerTrimNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);

wiskSessions = unique(sessionInfo.session(strcmp(sessionInfo.preOrPost, 'pre')));
noWiskSessions = unique(sessionInfo.session(strcmp(sessionInfo.preOrPost, 'post')));
wiskLightBins = ismember({data.session}, wiskSessions) & [data.isLightOn];
wiskBins = ismember({data.session}, wiskSessions) & ~[data.isLightOn];
lightBins = ismember({data.session}, noWiskSessions) & [data.isLightOn];
neitherBins = ismember({data.session}, noWiskSessions) & ~[data.isLightOn];

conditionBins = {wiskLightBins, wiskBins, lightBins, neitherBins};
binNames = {'wisk+light', 'wisk', 'light', 'neither'};
bins = zeros(1,length(data));
for i = 1:4; bins(conditionBins{i}) = i; end


% exclude invalid trails, including those where max paw height is less than obs height
heights = [data.obsHeightsVid];
binEdges = linspace(heightLims(1), heightLims(2), heightBinNum+1);
heightBins = discretize(heights, binEdges);
pawHgts = cell2mat(cellfun(@(x,ind) max(squeeze(x{data(ind).firstPawOver}(end,3,:))), ...
          {data.modifiedLocations}, num2cell(1:length(data)), 'UniformOutput', false)) * 1000;
bins(heights>pawHgts) = 0;
bins(~validBins) = 0;
fprintf('%.3f of trials paw was not higher than obs..\n', ...
    sum(heights>pawHgts)/length(data));


% initializations
if beforeObsOnly; xLims(2) = 0; end
allInds = find(bins~=0);
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



% PLOT EACH CONDITION HEIGHT TRAJECTORIES
figure('name', 'sensoryDependence', 'menubar', 'none', 'color', 'white', 'position', [100 100 500 900], 'InvertHardCopy', 'off');

for i = 1:length(binNames)
    
    subaxis(length(binNames)+1,1,i)
    colors = interp2(1:3, 1:2, cat(1,conditionColors(i,:)*colorFading,conditionColors(i,:)), 1:3, linspace(1,2,heightBinNum)');
    
    
    % plot control mean
    xOffset = nanmean(squeeze(locationsMod(bins==i,1,1)));
    xsCtl = squeeze(locationsControl(bins==i,1,:));
    zsCtl = squeeze(locationsControl(bins==i,3,:));
    plot(nanmean(xsCtl,1) - nanmean(xsCtl(:,1)) + xOffset, nanmean(zsCtl,1), ...
        'color', [0 0 0], 'linewidth', lineWid); hold on
    title(binNames{i})

    
    
    for j = 1:heightBinNum
        xs = squeeze(locationsMod(bins==i & heightBins==j,1,:));
        zs = squeeze(locationsMod(bins==i & heightBins==j,3,:));
        
        plot(nanmean(xs,1), nanmean(zs,1), ...
            'color', colors(j,:), 'linewidth', lineWid); hold on
        
        % add cylinder
        z = mean(heights(heightBins==j & bins==i)) / 1000;
        circ = rectangle('position', [0-obsRadius, z-2*obsRadius, 2*obsRadius, 2*obsRadius], ...
            'curvature', [1 1], 'facecolor', [colors(j,:) 1], 'edgecolor', 'none');
    end
    
    % add floor
    line(xLims, [0 0], 'linewidth', 1, 'color', [0 0 0]);
    
    % pimp appearance
    daspect([1 1 1]);
    set(gca, 'xlim', xLims, 'ylim', zLims);
    axis off
end






% PLOT CONDITION COMPARISON
subaxis(length(binNames)+1,1,length(binNames)+1)

% plot control mean
xOffset = nanmean(squeeze(locationsMod(bins~=0,1,1)));
xsCtl = squeeze(locationsControl(bins~=0,1,:));
zsCtl = squeeze(locationsControl(bins~=0,3,:));
plot(nanmean(xsCtl,1) - nanmean(xsCtl(:,1)) + xOffset, nanmean(zsCtl,1), ...
    'color', [0 0 0], 'linewidth', lineWid); hold on

for i = 1:length(binNames)
    xs = squeeze(locationsMod(bins==i,1,:));
    zs = squeeze(locationsMod(bins==i,3,:));
    plot(nanmean(xs,1), nanmean(zs,1), ...
        'color', [conditionColors(i,:)], 'linewidth', lineWid); hold on
end

% add cylinder
z = mean(heights) / 1000;
circ = rectangle('position', [0-obsRadius, z-2*obsRadius, 2*obsRadius, 2*obsRadius], ...
    'curvature', [1 1], 'facecolor', [colors(j,:) 1], 'edgecolor', 'none');
    
% add floor
line(xLims, [0 0], 'linewidth', 1, 'color', [0 0 0]);

% pimp appearance
daspect([1 1 1]);
set(gca, 'xlim', xLims, 'ylim', zLims);
axis off

blackenFig;




