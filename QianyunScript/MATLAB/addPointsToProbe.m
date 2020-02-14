function addPointsToProbe(session, LM, offset, opts)


if ~exist('offset', 'var'); offset = 0; end


% settings
s.showPCPoints = true;
s.showGCPoints = true;

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end


% get probe mapping file
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')
mapFile = ephysInfo.map{strcmp(session, ephysInfo.session)};

load(fullfile(getenv('OBSDATADIR'), 'ephys', 'channelMaps', 'kilosort', [mapFile '.mat']), ...
    'xcoords', 'ycoords', 'channelNum_OpenEphys', 'connected');

probeDepth = [];
for i = 1:64
    ind = find(channelNum_OpenEphys == i);
    probeDepth(ind, 1) = ind;
    probeDepth(ind, 2) = ycoords(i);
end


% get channel info of the recording session
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysChannelsInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysChannelsInfo.xlsx'), 'Sheet', 'EphysChannelsInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')


PCProbeDepth = ephysChannelsInfo.PCProbeDepth(strcmp(session, ephysChannelsInfo.session));
PCChannels = ephysChannelsInfo.PCChannels(strcmp(session, ephysChannelsInfo.session));
GCProbeDepth = ephysChannelsInfo.GCProbeDepth(strcmp(session, ephysChannelsInfo.session));
GCChannels = ephysChannelsInfo.GCChannels(strcmp(session, ephysChannelsInfo.session));

% get rid of nans 
PCProbeDepth = PCProbeDepth(~isnan(PCProbeDepth));
PCChannels = PCChannels(~isnan(PCChannels));
GCProbeDepth = GCProbeDepth(~isnan(GCProbeDepth));
GCChannels = GCChannels(~isnan(GCChannels));

% calculate the actual depth for Perkinje cell layer and channels that have
% good signals on.
PCDepth = zeros(size(PCProbeDepth, 1), 1);
for i = 1:length(PCProbeDepth)
    PCDepth(i, 1) = PCProbeDepth(i, 1) - probeDepth(PCChannels(i, 1), 2); 
end

GCDepth = zeros(size(GCProbeDepth, 1), 1);
for i = 1:length(GCProbeDepth)
    GCDepth(i, 1) = GCProbeDepth(i, 1) - probeDepth(GCChannels(i, 1), 2);
end


% getting parameters from the linear model
avg = LM.avg;
dirV = LM.dirVect;
if dirV(1, 2) > 0
    dirV = -dirV;
end

% determine the location of the brain surface cross point
MLCoord = avg(1, 1);
inds = find(LBSCoords(:, 1) == MLCood);
avgDVCoord = mean(LBSCoords(inds, 2));
tempDistance = (avg(1, 2) - avgDVCoord)/dirV(1, 2); % a hacky way of moving the avg point to where near the brain surface
tempBS_crossPoint = avg + dirV * tempDistance;

APCoord = tempBS_crossPoint(1, 2);
adjacentAPCoord = [LBSCoords(find(LBSCoord(:, 3) < APCoord, 'last'), 3), LBSCoords(find(LBSCoords(:, 3) > APCoord, 'first'), 3)];
MLCoord =  tempBS_crossPoint(1, 1);
DVCoords1 = LBSCoords(find((LBSCoords(:, 1) == MLCoord) & (LBSCoords(:, 3) == adjacentAPCoord(1))), 2);
DVCoords2 = LBSCoords(find((LBSCoords(:, 1) == MLCoord) & (LBSCoords(:, 3) == adjacentAPCoord(2))), 2);
DVCoords = cat(1, DVCoords1, DVCoords2);
DVCoord = mean(DVCoords);
clear DVCoords
distance = (avg(1, 2) - DVCoord)/dirV(1, 2);
BS_crossPoint = avg + dirV * distance;

APCoord = BS_crossPoint(1, 3);
weights = [abs(APCoord - adjacentAPCoord(1))/abs(diff(adjacentAPCoord)), abs(APCoord - adjacentAPCoord(2)/abs(diff(adjacentAPCoord))];
DVCoord = mean(DVCoords1)*weights(1) + mean(DVCoords2)*weights(2);   
distance = (avg(1, 2) - DVCoord)/dirV(1, 2);
BS_crossPoint = avg + dirV * distance;

% plot brain surface cross point
hold on
plot3(BS_crossPoint(:, 1), BS_crossPoint(:, 3), BS_crossPoint(:, 2), '.k', 'MarkerSize', 30)
points = [BS_crossPoint ; avg];
plot3(points(:,1),points(:,3),points(:,2),'-k','LineWidth',3) 



% plot PC layer cross points
if showPCPoints
    for i = 1:length(PCDepth)
        
        PC_crossPoint = (PCDepth(i) + offset)*dirV + BS_crossPoint;
        hold on
        plot3(PC_crossPoint(:, 1), PC_crossPoint(:, 3), PC_crossPoint(:, 2), '.r', 'MarkerSize', 30)        
    end
end


% plot good channel points
if showGCPoints
    for i = 1:length(GCDepth)
        
        goodChannel = (GCDepth(i) + offset)*dirV + BS_crossPoint;
        hold on
        plot3(goodChannel(:, 1), goodChannel(:, 3), goodChannel(:, 2), '.c', 'Markersize', 30);
        points = [goodChannel; avg];
        plot3(points(:,1),points(:,3),points(:,2),'-k','LineWidth',3)
        
    end
end



end