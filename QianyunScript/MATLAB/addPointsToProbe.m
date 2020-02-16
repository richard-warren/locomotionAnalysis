function addPointsToProbe(session, LM, BSCoord, offset, shankNum, opts)


if ~exist('offset', 'var'); offset = 0; end
if ~exist('shankNum', 'var'); shankNum = 1; end

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

PCShankNum =  ephysChannelsInfo.PCShankNum(strcmp(session, ephysChannelsInfo.session));
PCProbeDepth = ephysChannelsInfo.PCProbeDepth(strcmp(session, ephysChannelsInfo.session));
PCChannels = ephysChannelsInfo.PCChannels(strcmp(session, ephysChannelsInfo.session));

GCShankNum = ephysChannelsInfo.GCShankNum(strcmp(session, ephysChannelsInfo.session));
GCProbeDepth = ephysChannelsInfo.GCProbeDepth(strcmp(session, ephysChannelsInfo.session));
GCChannels = ephysChannelsInfo.GCChannels(strcmp(session, ephysChannelsInfo.session));

% get rid of nans 
PCShankNum = PCShankNum(~isnan(PCShankNum));
PCProbeDepth = PCProbeDepth(~isnan(PCProbeDepth));
PCChannels = PCChannels(~isnan(PCChannels));
GCShankNum = GCShankNum(~isnan(GCShankNum));
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
if dirV(1, 3) < 0
    dirV = -dirV; % make sure that the direction vector is pointing downward
end

hold on 
plot3(avg(:, 1), avg(:, 2), avg(:, 3), '.b', 'MarkerSize', 30)


% determine the location of the brain surface cross point
MLCoord = avg(1, 1);
MLCoord = BSCoords(find(LBSCoords(:, 1) > round(MLCoord), 1));
inds = find(BSCoords(:, 1) == MLCoord);
avgDVCoord = mean(BSCoords(inds, 3));
tempDistance = (avg(1, 3) - avgDVCoord)/dirV(1, 3); % a hacky way of moving the avg point to where near the brain surface
if tempDistance > 0; tempDistance = -tempDistance; end
tempBS_crossPoint = avg + dirV * tempDistance;
hold on
plot3(tempBS_crossPoint(:, 1), tempBS_crossPoint(:, 2), tempBS_crossPoint(:, 3), '.k', 'MarkerSize', 30)
points = [tempBS_crossPoint ; avg];
plot3(points(:,1),points(:,2),points(:,3),'-k','LineWidth',3) 

APCoord = tempBS_crossPoint(1, 2);
adjacentAPCoord = [BSCoords(find(LBSCoords(:, 2) < APCoord, 1, 'last'), 2), BSCoords(find(LBSCoords(:, 2) > APCoord, 1, 'first'), 2)];
MLCoord =  tempBS_crossPoint(1, 1);
clear tempBS_crossPoint inds avgDVCoord tempDistance
MLCoord = LBSCoords(find(LBSCoords(:, 1) > round(MLCoord), 1));
DVCoords1 = BSCoords(find((BSCoords(:, 1) == MLCoord) & (BSCoords(:, 2) == adjacentAPCoord(1))), 3);
DVCoords2 = BSCoords(find((BSCoords(:, 1) == MLCoord) & (BSCoords(:, 2) == adjacentAPCoord(2))), 3);
DVCoords = cat(1, DVCoords1, DVCoords2);
DVCoord = mean(DVCoords);
distance = (avg(1, 3) - DVCoord)/dirV(1, 3);
if distance > 0; distance = -distance; end
BS_crossPoint = avg + dirV * distance;
clear DVCoords MLCoord

APCoord = BS_crossPoint(1, 2);
weights = [abs(APCoord - adjacentAPCoord(1))/abs(diff(adjacentAPCoord)), abs(APCoord - adjacentAPCoord(2))/abs(diff(adjacentAPCoord)) ];
DVCoord = mean(DVCoords1)*weights(2) + mean(DVCoords2)*weights(1);   
distance = (avg(1, 3) - DVCoord)/dirV(1, 3);
if distance > 0; distance = -distance; end
BS_crossPoint = avg + dirV * distance;
clear weights APCoord DVCoord distance DVCoords1 DVCoords2 DVCoord

% plot brain surface cross point
plot3(BS_crossPoint(:, 1), BS_crossPoint(:, 2), BS_crossPoint(:, 3), '.k', 'MarkerSize', 30)
points = [BS_crossPoint ; avg];
plot3(points(:,1),points(:,2),points(:,3),'-k','LineWidth',3) 





% plot PC layer cross points
if showPCPoints
    disp('plotting PC cross points...');
    PC_inds = find(PCShankNum == shankNum);
    PCChannels = PCChannels(PC_inds, :);
    PCDepth = PCDepth(PC_inds, :);

    for i = 1:length(PCDepth)
        
        PC_crossPoint = (PCDepth(i) + offset)*dirV + BS_crossPoint;
        hold on
        plot3(PC_crossPoint(:, 1), PC_crossPoint(:, 2), PC_crossPoint(:, 3), '.r', 'MarkerSize', 30)        
    end
    
    disp('Done');    
end


% plot good channel points
if showGCPoints
    disp('plotting GC points...');
    
    inds = find(GCShankNum == shankNum);
    GCChannels = GCChannels(inds, :);
    GCDepth = GCDepth(inds, :);
    
    for i = 1:length(GCDepth)
        
        goodChannel = (GCDepth(i) + offset)*dirV + BS_crossPoint;
        hold on
        plot3(goodChannel(:, 1), goodChannel(:, 3), goodChannel(:, 2), '.c', 'Markersize', 30);
        points = [goodChannel; avg];
        plot3(points(:,1),points(:,3),points(:,2),'-k','LineWidth',3)
        
    end    
    disp('Done');
end



end