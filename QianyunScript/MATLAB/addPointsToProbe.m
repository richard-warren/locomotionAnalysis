function [BS_crossPoint, PC_crossPoint, goodChannelPoint, isStepModulated] = addPointsToProbe(session, LM, BSCoords, offset, shankNum, opts)


if ~exist('offset', 'var'); offset = 0; end
if ~exist('shankNum', 'var'); shankNum = 1; end

% settings
s.showPCPoints = true;
s.showGCPoints = true;
s.noDiIMode = false;

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end


% get probe mapping file
disp('getting probe mapping...');
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
disp('getting ephys channel info...');

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysChannelsInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysChannelsInfo.xlsx'), 'Sheet', 'EphysChannelsInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')

PCShankNum =  ephysChannelsInfo.PCShankNum(strcmp(session, ephysChannelsInfo.session));
PCProbeDepth = ephysChannelsInfo.PCProbeDepth(strcmp(session, ephysChannelsInfo.session));
PCChannels = ephysChannelsInfo.PCChannels(strcmp(session, ephysChannelsInfo.session));

GCShankNum = ephysChannelsInfo.GCShankNum(strcmp(session, ephysChannelsInfo.session));
GCProbeDepth = ephysChannelsInfo.GCProbeDepth(strcmp(session, ephysChannelsInfo.session));
GCChannels = ephysChannelsInfo.GCChannels(strcmp(session, ephysChannelsInfo.session));

isStepModulated = ephysChannelsInfo.stepModulated(strcmp(session, ephysChannelsInfo.session));


% get rid of nans 
PCShankNum = PCShankNum(~isnan(PCShankNum));
PCProbeDepth = PCProbeDepth(~isnan(PCProbeDepth));
PCChannels = PCChannels(~isnan(PCChannels));
GCShankNum = GCShankNum(~isnan(GCShankNum));
GCProbeDepth = GCProbeDepth(~isnan(GCProbeDepth));
GCChannels = GCChannels(~isnan(GCChannels));
isStepModulated = isStepModulated(~isnan(isStepModulated));

% calculate the actual depth for where I got Perkinje cell signals and for channels that have
% good signals.
PCDepth = zeros(size(PCProbeDepth, 1), 1);
for i = 1:length(PCProbeDepth)
    PCDepth(i, 1) = PCProbeDepth(i, 1) - probeDepth(PCChannels(i, 1), 2); 
end

GCDepth = zeros(size(GCProbeDepth, 1), 1);
for i = 1:length(GCProbeDepth)
    GCDepth(i, 1) = GCProbeDepth(i, 1) - probeDepth(GCChannels(i, 1), 2);
end


% plot brain surface cross point
disp('adding points...');
if ~s.noDiIMode
    
    % the code below deals with the reconstruction of dii-coated probe traces    
    % get avg point and direction vector from the linear model
    avg = LM.avg;
    dirV = LM.dirVect;
    if dirV(1, 3) < 0
        dirV = -dirV; % make sure that the direction vector is pointing downward
    end
    
%     % (optional) plotting the avg point
%     hold on
%     plot3(avg(:, 1), avg(:, 2), avg(:, 3), '.b', 'MarkerSize', 30)
    
    
    % coarsely determine the location of the brain surface cross point
    disp('adding brain surface cross point...');
    MLCoord = avg(1, 1);
    BS_MLCoords = unique(BSCoords(:, 1));
    MLCoord = BS_MLCoords(knnsearch(BS_MLCoords, MLCoord), 1);
    inds = find(BSCoords(:, 1) == MLCoord);
    avgDVCoord = mean(BSCoords(inds, 3));
    tempDistance = (avg(1, 3) - avgDVCoord)/dirV(1, 3); % a hacky way of moving the avg point to where near the brain surface
    if tempDistance > 0; tempDistance = -tempDistance; end
    tempBS_crossPoint = avg + dirV * tempDistance;
    
%     % (optional) plot the tempBS_crossPoint, just for quality checking
%     hold on
%     plot3(tempBS_crossPoint(:, 1), tempBS_crossPoint(:, 2), tempBS_crossPoint(:, 3), '.k', 'MarkerSize', 30)
%     points = [tempBS_crossPoint ; avg];
%     plot3(points(:,1),points(:,2),points(:,3),'-g','LineWidth',3)
%     
    
    % fine tune the DV coordinates for the brain surface cross point by
    % using the info from only the adjacent brain sections
    APCoord = tempBS_crossPoint(1, 2);
    if APCoord < 0;
        adjacentAPCoord = 0;
    else
        adjacentAPCoord = [BSCoords(find(BSCoords(:, 2) < APCoord, 1, 'last'), 2), BSCoords(find(BSCoords(:, 2) > APCoord, 1, 'first'), 2)];
    end
    
    MLCoord =  tempBS_crossPoint(1, 1);
    MLCoord = BSCoords(find(BSCoords(:, 1) > round(MLCoord), 1));
    
    % dealing with the situation if there is only one adjacent brain section (most anterior or posterior)
    if length(adjacentAPCoord) == 1 
        DVCoords = BSCoords(find((BSCoords(: , 1) == MLCoord) & (BSCoords(:, 2) == adjacentAPCoord)), 3);
    else
        DVCoords1 = BSCoords(find((BSCoords(: , 1) == MLCoord) & (BSCoords(:, 2) == adjacentAPCoord(1))), 3);
        DVCoords2 = BSCoords(find((BSCoords(:, 1) == MLCoord) & (BSCoords(:, 2) == adjacentAPCoord(2))), 3);
        DVCoords = cat(1, DVCoords1, DVCoords2);
    end
    
    DVCoord = mean(DVCoords);
    BS_distance = (avg(1, 3) - DVCoord)/dirV(1, 3);
    if BS_distance > 0; BS_distance = -BS_distance; end
    BS_crossPoint = avg + dirV * BS_distance;
    clear DVCoords MLCoord
    
    if length(adjacentAPCoord) ~= 1
        APCoord = BS_crossPoint(1, 2);
        adjacentAPCoord = [BSCoords(find(BSCoords(:, 2) < APCoord, 1, 'last'), 2), BSCoords(find(BSCoords(:, 2) > APCoord, 1, 'first'), 2)];
        weights = [abs(APCoord - adjacentAPCoord(1))/abs(diff(adjacentAPCoord)), abs(APCoord - adjacentAPCoord(2))/abs(diff(adjacentAPCoord)) ];
        DVCoord = mean(DVCoords1)*weights(2) + mean(DVCoords2)*weights(1);
        BS_distance = (avg(1, 3) - DVCoord)/dirV(1, 3);
        if BS_distance > 0; BS_distance = -BS_distance; end
        BS_crossPoint = avg + dirV * BS_distance;
        clear weights APCoord DVCoord distance DVCoords1 DVCoords2 DVCoord
    end
    
else
    
    % the code below applies for reconstructing no dii traces only
    tempBS_crossPoint = LM.BS_crossPoint;
    dirV = LM.dirV;
    if dirV(1, 3)<0; dirV = -dirV; end
    
    APCoord = tempBS_crossPoint(1, 2);
    adjacentAPCoord = [BSCoords(find(BSCoords(:, 2) < APCoord, 1, 'last'), 2), BSCoords(find(BSCoords(:, 2) > APCoord, 1, 'first'), 2)];
    MLCoord = tempBS_crossPoint(1, 1);
    MLCoord = BSCoords(find(BSCoords(:, 1) > round(MLCoord), 1));
    DVCoords1 = BSCoords(find((BSCoords(: , 1) == MLCoord) & (BSCoords(:, 2) == adjacentAPCoord(1))), 3);
    DVCoords2 = BSCoords(find((BSCoords(:, 1) == MLCoord) & (BSCoords(:, 2) == adjacentAPCoord(2))), 3);
    weights = [abs(APCoord - adjacentAPCoord(1))/abs(diff(adjacentAPCoord)), abs(APCoord - adjacentAPCoord(2))/abs(diff(adjacentAPCoord)) ];
    DVCoord = mean(DVCoords1)*weights(2) + mean(DVCoords2)*weights(1);
    BS_crossPoint = [tempBS_crossPoint(1, :), DVCoord];
    BS_distance = 0;   

end

% plot brain surface cross point
plot3(BS_crossPoint(:, 1), BS_crossPoint(:, 2), BS_crossPoint(:, 3), '.k', 'MarkerSize', 30)
if ~s.noDiIMode; points = [BS_crossPoint ; avg]; plot3(points(:,1), points(:,2), points(:,3),'-g','LineWidth',3); end



% plot PC layer cross points
if s.showPCPoints
    disp('plotting PC cross points...');
    PC_inds = find(PCShankNum == shankNum);
    PCChannels = PCChannels(PC_inds, :);
    PCDepth = PCDepth(PC_inds, :);
    if dirV(1, 3) < 0; dirV = -dirV; end
     
    PC_crossPoint = nan(length(PCDepth), 3);
    for i = 1:length(PCDepth)
        PC_crossPoint(i, :) = (PCDepth(i) + offset)*dirV + BS_crossPoint;
        hold on
        plot3(PC_crossPoint(i, 1), PC_crossPoint(i, 2), PC_crossPoint(i, 3), '.r', 'MarkerSize', 30)   
        points = [PC_crossPoint(i, :); BS_crossPoint];
        plot3(points(:,1),points(:,2),points(:,3),'-g','LineWidth',3)
    end

    disp('Done');    
end


% plot good channel points
if s.showGCPoints
    disp('plotting GC points...');
    
    inds = find(GCShankNum == shankNum);
    GCChannels = GCChannels(inds, :);
    GCDepth = GCDepth(inds, :);
    isStepModulated = isStepModulated(inds, :);
    
    goodChannelPoint = nan(length(GCDepth), 3);
    for i = 1:length(GCDepth)
        
        goodChannelPoint(i, :) = (GCDepth(i) + offset)*dirV + BS_crossPoint;
        hold on
        plot3(goodChannelPoint(i, 1), goodChannelPoint(i, 2), goodChannelPoint(i, 3), '.c', 'Markersize', 30);
        points = [goodChannelPoint(i, :); BS_crossPoint];
        plot3(points(:,1),points(:,2),points(:,3),'-g','LineWidth',3)
        
    end    
    disp('Done');
end



end