% function sensoryDependenceKinematics

% wiskSessions = {'171117_001', '171118_001', '180225_000', '180225_001', '180225_002', '180226_000', '180226_001', '180226_002'};
% noWiskSessions = {'171121_000', '171122_001', '180228_000', '180228_001', '180228_002', '180301_000', '180301_001', '180301_002'};

validBins= cellfun(@(x) x(3), {data.modStepNum}) > 1;
wiskLightBins = ismember({data.session}, wiskSessions) & [data.isLightOn] & validBins;
wiskBins = ismember({data.session}, wiskSessions) & ~[data.isLightOn] & validBins;
lightBins = ismember({data.session}, noWiskSessions) & [data.isLightOn] & validBins;
neitherBins = ismember({data.session}, noWiskSessions) & ~[data.isLightOn] & validBins;
allBins = wiskLightBins | wiskBins | lightBins | neitherBins;
allInds = find(allBins);
conditionBins = {wiskLightBins, wiskBins, lightBins, neitherBins};
conditionNames = {'wisk+light', 'wisk', 'light', 'neither'};



% settings
azimuth = 0;
elevation = 0;
obsThickness = .003;
obsHeight = .009;
colors = [hsv(4)];
controlColor = repmat(0, 1, 3);
xLims = [-.02 .06];
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

for i = allInds
    
    % find fore paw that is first to get over osbtacle
    isStepping = ~isnan(data(1).modifiedStepIdentities(:,:));
    lastModStepInds = table2array(rowfun(@(x)(find(x,1,'last')), table(isStepping')));
    [~, firstPawOverInd] = min(lastModStepInds .* [nan 1 1 nan]'); % mask out hind paws, inds 1 and 4
    
    locationsMod(i,:,:) = squeeze(data(i).modifiedLocationsInterp{firstPawOverInd}(end,:,:));
    locationsControl(i,:,:) = squeeze(data(i).controlLocationsInterp{firstPawOverInd}(end,:,:));
end


%%
close all;
figure('color', 'white', 'position', [500 500 900 300], 'InvertHardCopy', 'off');

for i = 1:length(conditionBins)
    
    xs = squeeze(locationsMod(conditionBins{i},1,:));
    ys = squeeze(locationsMod(conditionBins{i},2,:));
    zs = squeeze(locationsMod(conditionBins{i},3,:));
    
    plot3(mean(xs,1), mean(ys,1), mean(zs,1), 'color', colors(i,:), 'linewidth', lineWid); hold on
    
end

% plot control trajectory
xOffset = nanmean(squeeze(locationsMod(:,1,1)));
xs = squeeze(locationsControl(allBins,1,:));
ys = squeeze(locationsControl(allBins,2,:));
zs = squeeze(locationsControl(allBins,3,:));
plot3(mean(xs,1) - mean(xs(:,1)) + xOffset, mean(ys,1), mean(zs,1), 'color', [0 0 0], 'linewidth', lineWid); hold on


daspect([1 1 1]);
set(gca, 'xlim', xLims, 'ylim', yLims, 'zlim', zLims, ...
    'view', [azimuth elevation], 'YDir', 'reverse');
axis off

% add obstacle
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
patch('Vertices', vertices, 'Faces', faces);

% add lines for sides of wheel
line([xLims; xLims]', [yLims; yLims], zeros(2,2), ...
    'color', [0 0 0], 'linewidth', 2) % x1

legend(conditionNames)


pimpFig; blackenFig;



