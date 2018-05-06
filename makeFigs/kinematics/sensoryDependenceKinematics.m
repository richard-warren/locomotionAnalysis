function sensoryDependenceKinematics(data, wiskSessions, noWiskSessions)


% validBins= cellfun(@(x) x(3), {data.modStepNum}) > 1;
validBins = true(1,length(data));
wiskLightBins = ismember({data.session}, wiskSessions) & [data.isLightOn] & validBins;
wiskBins = ismember({data.session}, wiskSessions) & ~[data.isLightOn] & validBins;
lightBins = ismember({data.session}, noWiskSessions) & [data.isLightOn] & validBins;
neitherBins = ismember({data.session}, noWiskSessions) & ~[data.isLightOn] & validBins;
allBins = wiskLightBins | wiskBins | lightBins | neitherBins;
allInds = find(allBins);
conditionBins = {wiskLightBins, wiskBins, lightBins, neitherBins};
conditionNames = {'wisk+light', 'wisk', 'light', 'neither'};
trials = 15;
maxSmpsThoughObs = 2;



% settings
azimuth = 0;
elevation = 0;
obsThickness = .003;
obsHeight = .009;
colors = hsv(4);
xLims = [-.05 .04];
yLims = [-0.0381 0.0381]; % 3 inches, which is width of wheel
zLims = [0 .015];
lineWid = 3;

% initializations
obsX = [0 obsThickness]-.5*obsThickness;
obsY = yLims;
obsZ = [0 obsHeight];

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
    if firstPawOverInd==2
        locationsMod(i,2,:) = -locationsMod(i,2,:);
        locationsControl(i,2,:) = -locationsControl(i,2,:);
    end
    
    % check if trial passes through obstacle
    throughObsBins(i,:) = squeeze(locationsMod(i,3,:)<obsHeight) & ... 
                          squeeze(locationsMod(i,1,:)>obsX(1)) & ...
                          squeeze(locationsMod(i,1,:)<obsX(2));
    if sum(throughObsBins(i,:)) > maxSmpsThoughObs
%         fprintf('trial %i passed though obstacle\n', i);
        locationsMod(i,:,:) = nan;
        locationsControl(i,:,:) = nan;
    end
end
fprintf('%.2f of trials passed though obs\n', ...
    sum(sum(throughObsBins,2)>maxSmpsThoughObs) / size(throughObsBins,1));



close all;
figure('menubar', 'none', 'color', 'white', 'position', [400 50 850 950], 'InvertHardCopy', 'off');

for i = 1:length(conditionNames)
    
    subaxis(length(conditionNames)+1,1,i)
    
    xs = squeeze(locationsMod(conditionBins{i},1,:));
    ys = squeeze(locationsMod(conditionBins{i},2,:));
    zs = squeeze(locationsMod(conditionBins{i},3,:));
    
    
    realInds = find(~isnan(xs(:,1)));
    trialInds = randperm(length(realInds), trials);
    for j = 1:length(trialInds)
        plot3(xs(realInds(trialInds(j)),:), ys(realInds(trialInds(j)),:), zs(realInds(trialInds(j)),:), ...
            'color', colors(i,:), 'linewidth', 1); hold on
    end
    
    % plot control mean
    xOffset = nanmean(squeeze(locationsMod(conditionBins{i},1,1)));
    xsCtl = squeeze(locationsControl(conditionBins{i},1,:));
    ysCtl = squeeze(locationsControl(conditionBins{i},2,:));
    zsCtl = squeeze(locationsControl(conditionBins{i},3,:));
    plot3(nanmean(xsCtl,1) - nanmean(xsCtl(:,1)) + xOffset, nanmean(ysCtl,1), nanmean(zsCtl,1), ...
        'color', [0 0 0], 'linewidth', lineWid); hold on
    
    % plot mean in final subplot
    subaxis(length(conditionNames)+1,1,length(conditionNames)+1)
    plot3(nanmean(xs,1), nanmean(ys,1), nanmean(zs,1), 'color', colors(i,:), 'linewidth', lineWid); hold on
end

% plot control trajectory
xOffset = nanmean(squeeze(locationsMod(:,1,1)));
xs = squeeze(locationsControl(allBins,1,:));
ys = squeeze(locationsControl(allBins,2,:));
zs = squeeze(locationsControl(allBins,3,:));
plot3(nanmean(xs,1) - nanmean(xs(:,1)) + xOffset, nanmean(ys,1), nanmean(zs,1), ...
    'color', [0 0 0], 'linewidth', lineWid); hold on

% pimp axes
vertices = [obsX(1) obsY(1) obsZ(1)
            obsX(1) obsY(2) obsZ(1)
            obsX(1) obsY(2) obsZ(2)
            obsX(1) obsY(1) obsZ(2)
            obsX(2) obsY(1) obsZ(1)
            obsX(2) obsY(2) obsZ(1)
            obsX(2) obsY(2) obsZ(2)
            obsX(2) obsY(1) obsZ(2)];
% specify which corners to connect for each of 6 faces (6 sides of obs)
faces = [1 2 3 4
         5 6 7 8
         1 2 6 5
         4 3 7 8
         2 6 7 3
         1 5 8 4];
for i = 1:length(conditionNames)+1
    
    subaxis(length(conditionNames)+1,1,i)
    
    daspect([1 1 1]);
    set(gca, 'xlim', xLims, 'ylim', yLims, 'zlim', zLims, ...
        'view', [azimuth elevation], 'YDir', 'reverse');
    axis off


    patch('Vertices', vertices, 'Faces', faces);

    % add lines for sides of wheel
    line([xLims; xLims]', [yLims; yLims], zeros(2,2), ...
        'color', [0 0 0], 'linewidth', 2) % x1
end

legend(gca, conditionNames, 'location', 'northeast');
% pos = get(leg, 'position');
% set(leg, 'position', [xLims(2)-pos(3) zLims(2)-pos(4) pos(3) pos(4)])

blackenFig;
saveas(gcf, [getenv('OBSDATADIR') 'figures\sensoryDependenceKinematics.png']);
savefig([getenv('OBSDATADIR') 'figures\sensoryDependenceKinematics.fig'])



