folder = uigetdir();
files = dir(fullfile(folder, '*.tif'));

includeDentate = false;
includeInt = true;
downsampleRate = 0.3;
thickness = 50;

if includeInt
    RICoords = importTiffStack(fullfile(folder, 'InterpositusRight.tif'), downsampleRate, thickness);
    LICoords = importTiffStack(fullfile(folder, 'InterpositusLeft.tif'), downsampleRate, thickness);
end

if includeDentate
    RDCoords = importTiffStack(fullfile(folder, 'DentateRight.tif'), downsampleRate, thickness);
    LDCoords = importTiffStack(fullfile(folder, 'DentateLeft.tif'), downsampleRate, thickness);

end

disp('All Done!');



%% import good channel locations

mouseID = 'cer12';
GC_columnInd = 6;
side_columnInd = 3;
stepModulated_columnInd = 7;

fileName = fullfile(getenv('OBSDATADIR'), 'ephys', 'ephysHistoData', 'ephysHistoData.mat');
temp = load(fileName);
ephysHistoData = temp.ephysHistoData;
clear temp

inds = find(strcmp(mouseID, ephysHistoData(:, 2)));
goodChannelRight = [];
goodChannelLeft = [];
isStepModulatedLeft = [];
isStepModulatedRight = [];

for i = 1:length(inds)
    if ~isempty(ephysHistoData{inds(i), GC_columnInd})
        if strcmp('right', ephysHistoData{inds(i), side_columnInd}) 
            goodChannelRight = cat(1, goodChannelRight, ephysHistoData{inds(i), GC_columnInd});
            isStepModulatedRight = cat(1, isStepModulatedRight, ephysHistoData{inds(i), stepModulated_columnInd});
        else
            goodChannelLeft = cat(1, goodChannelLeft, ephysHistoData{inds(i), GC_columnInd});
            isStepModulatedLeft = cat(1, isStepModulatedLeft, ephysHistoData{inds(i), stepModulated_columnInd});
        end            
    end
end


%% Plot!!!
markerSize = 25;
frontSize = 14;

% figure('Color', 'white', 'position', get(0,'ScreenSize'));
% 
% box off


hf = figure('Color', 'white', 'position', get(0,'ScreenSize'));
box off
hold on


% front view
h1 = subplot(2, 2, 1);
set(h1, 'position', [0.08, 0.5, 0.40, 0.5] );


xrange = [1400, max(LICoords(:, 1))*2+500];
yrange = [200, max(LICoords(:, 2))+300]; % just for visualizing purposes

% plot DCN and good units
if includeInt
    plotRegions3D(LICoords, 20, [1, 0.74, 0.35], xrange, yrange, 0.1);
    hold on
    plotRegions3D(RICoords, 20, [1, 0.74, 0.35], xrange, yrange, 0.1);
end

if includeDentate
    plotRegions3D(LDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange, 0.1);
    hold on
    plotRegions3D(RDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange, 0.1);
end

points_left = goodChannelLeft;
points_left(:, 2) = 400;
points_right = goodChannelRight;
points_right(:, 2) = 400;

plot3(points_left(:, 1), points_left(:, 2), points_left(:, 3), '.', 'Markersize', markerSize, 'Color', [0.4 0.69 0.90]);
plot3(points_right(:, 1), points_right, points_right(:, 3), '.', 'Markersize', markerSize, 'Color', [0.46 0.67 0.18]);

% adding numbers & color code units to indicate which is within and which is out of DCN
for i = 1:size(goodChannelLeft, 1)
    text(points_left(i, 1)+80, points_left(i, 2), points_left(i, 3), num2str(i), 'Color',  [0.46 0.67 0.18],...
        'FontWeight', 'bold', 'FontSize', frontSize);
    
    if ifInDCN(goodChannelLeft(i, :), LICoords)
        plot3(points_left(i, 1), points_left(i, 2), points_left(i, 3), '.', 'Markersize', 20, 'Color', [1 0.4 0.4]);
    else
        plot3(points_left(i, 1), points_left(i, 2), points_left(i, 3), '.', 'Markersize', 20, 'Color', [0.65 0.65 0.65]);
    end
    
%     if isStepModulatedLeft(i)
%         plot3(points_left(i, 1), points_left(i, 2), points_left(i, 3), '.', 'Markersize', 20, 'Color', [0.8 0.59 0.93]);
%     else
%         plot3(points_left(i, 1), points_left(i, 2), points_left(i, 3), '.', 'Markersize', 20, 'Color', [0.14 0.14 0.14]);
%     end
    
end

for i = 1:size(goodChannelRight, 1)
    text(points_right(i, 1)+80, points_right(i, 2), points_right(i, 3), num2str(i), 'Color',  [0.4 0.69 0.90],...
        'FontWeight', 'bold',  'FontSize', frontSize);
    
    if ifInDCN(goodChannelRight(i, :), RICoords)
        plot3(points_right(i, 1), points_right(i, 2), points_right(i, 3), '.', 'Markersize', 20, 'Color', [1 0.4 0.4]);
    else
        plot3(points_right(i, 1), points_right(i, 2), points_right(i, 3), '.', 'Markersize', 20, 'Color', [0.65 0.65 0.65]);
    end
    
%     if isStepModulatedRight(i)
%         plot3(points_right(i, 1), points_right(i, 2), points_right(i, 3), '.', 'Markersize', 20, 'Color',  [0.8 0.59 0.93]);
%     else
%         plot3(points_right(i, 1), points_right(i, 2), points_right(i, 3), '.', 'Markersize', 20, 'Color', [0.14 0.14 0.14]);
%     end
    
    
    
    
end


% determine the dorsal/ventral/lateral/medial edges for each brain section
% of DCN
if exist('LICoords', 'var')
    uniqueLICoords = unique(LICoords(:, 2));
    LIEdgeDorsal = nan(length(uniqueLICoords), 1);
    LIEdgeVentral = nan(length(uniqueLICoords), 1);
    LIEdgeLateral = nan(length(uniqueLICoords), 1);
    LIEdgeMedial = nan(length(uniqueLICoords), 1);
    
    for i = 1:length(uniqueLICoords)
        section = LICoords(find(LICoords(:, 2) == uniqueLICoords(i)), :);
        temp = sort(section(:, 3));
        LIEdgeDorsal(i, 1) = temp(1);
        LIEdgeVentral(i, 1) = temp(end);
        temp = sort(section(:, 1));
        LIEdgeLateral(i, 1) = temp(1);
        LIEdgeMedial(i, 1) = temp(end);        
    end    
    
end

if exist('RICoords', 'var')
    uniqueRICoords = unique(RICoords(:, 2));
    RIEdgeDorsal = nan(length(uniqueRICoords), 1);
    RIEdgeVentral = nan(length(uniqueRICoords), 1);
    RIEdgeLateral = nan(length(uniqueRICoords), 1);
    RIEdgeMedial = nan(length(uniqueRICoords), 1);
    
    for i = 1:length(uniqueRICoords)
        section = RICoords(find(RICoords(:, 2) == uniqueRICoords(i)), :);
        temp = sort(section(:, 3));
        RIEdgeDorsal(i, 1) = temp(1);
        RIEdgeVentral(i, 1) = temp(end);
        temp = sort(section(:, 1));
        RIEdgeLateral(i, 1) = temp(end);
        RIEdgeMedial(i, 1) = temp(1);        
    end  
end
view(0, 0)
daspect([1 1 1])
title(['Good Units Location ( ' mouseID ' front view )']);





% topdown view
h2 = subplot(2, 2, 2);
set(h2, 'position', [0.54, 0.51, 0.4, 0.5])


plotRegions3D(LICoords, 20, [1, 0.74, 0.35], xrange, yrange, 0.1);
hold on
plotRegions3D(RICoords, 20, [1, 0.74, 0.35], xrange, yrange, 0.1);
if includeDentate
    plotRegions3D(LDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange, 0.1);
    plotRegions3D(RDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange, 0.1);
end


points_left = goodChannelLeft;
points_left(:, 3) = 1500;
points_right = goodChannelRight;
points_right(:, 3) = 1500;

plot3(points_left(:, 1), points_left(:, 2), points_left(:, 3), '.', 'Markersize', markerSize, 'Color', [0.4 0.69 0.90]);
plot3(points_right(:, 1), points_right(:, 2), points_right(:, 3), '.', 'Markersize', markerSize, 'Color', [0.46 0.67 0.18]);

% adding numbers & color code units to indicate which is within and which is out of DCN
for i = 1:size(goodChannelLeft, 1)
    text(points_left(i, 1)+80, points_left(i, 2), points_left(i, 3), num2str(i), 'Color',  [0.46 0.67 0.18],...
        'FontWeight', 'bold', 'FontSize', frontSize);
    
    if ifInDCN(goodChannelLeft(i, :), LICoords)
        plot3(points_left(i, 1), points_left(i, 2), points_left(i, 3), '.', 'Markersize', 20, 'Color', [1 0.4 0.4]);
    else
        plot3(points_left(i, 1), points_left(i, 2), points_left(i, 3), '.', 'Markersize', 20, 'Color', [0.65 0.65 0.65]);
    end
    
%      if isStepModulatedLeft(i)
%         plot3(points_left(i, 1), points_left(i, 2), points_left(i, 3), '.', 'Markersize', 20, 'Color', [0.8 0.59 0.93]);
%     else
%         plot3(points_left(i, 1), points_left(i, 2), points_left(i, 3), '.', 'Markersize', 20, 'Color', [0.14 0.14 0.14]);
%     end
    
    
    
end

for i = 1:size(goodChannelRight, 1)
    text(points_right(i, 1)+80, points_right(i, 2), points_right(i, 3), num2str(i), 'Color',  [0.4 0.69 0.90],...
        'FontWeight', 'bold',  'FontSize', frontSize);
    if ifInDCN(goodChannelRight(i, :), RICoords)
        plot3(points_right(i, 1), points_right(i, 2), points_right(i, 3), '.', 'Markersize', 20, 'Color', [1 0.4 0.4]);
    else
        plot3(points_right(i, 1), points_right(i, 2), points_right(i, 3), '.', 'Markersize', 20, 'Color', [0.65 0.65 0.65]);
    end
    
    
%     if isStepModulatedRight(i)
%         plot3(points_right(i, 1), points_right(i, 2), points_right(i, 3), '.', 'Markersize', 20, 'Color',  [0.8 0.59 0.93]);
%     else
%         plot3(points_right(i, 1), points_right(i, 2), points_right(i, 3), '.', 'Markersize', 20, 'Color', [0.14 0.14 0.14]);
%     end
    
end


LI_lateralEdgePoint = [LIEdgeLateral, uniqueLICoords];
LI_medialEdgePoint = [LIEdgeMedial, uniqueLICoords];
plot(LI_lateralEdgePoint(:, 1), LI_lateralEdgePoint(:, 2), '-', 'Color', [0.96, 0.5, 0.3, 0.5], 'LineWidth', 2);
plot(LI_medialEdgePoint(:, 1), LI_medialEdgePoint(:, 2), '-', 'Color', [0.96, 0.5, 0.3, 0.5], 'LineWidth', 2);

RI_lateralPoint = [RIEdgeLateral, uniqueRICoords];
RI_medialEdgePoint = [RIEdgeMedial, uniqueRICoords];
plot(RI_lateralPoint(:, 1), RI_lateralPoint(:, 2), '-', 'Color', [0.96, 0.5, 0.3, 0.5], 'LineWidth', 2);
plot(RI_medialEdgePoint(:, 1), RI_medialEdgePoint(:, 2), '-', 'Color', [0.96, 0.5, 0.3, 0.5], 'LineWidth', 2);


view(2)
daspect([1 1 1])
title(['Good Units Location ( ' mouseID ' topdown view )']);





% side view, left side
h3 = subplot(2, 2, 3);
set(h3, 'position', [0.08, 0.11, 0.4, 0.4]);


plotRegions3D(LICoords, 20, [1, 0.74, 0.35], xrange, yrange, 0.1);
hold on

points_left = goodChannelLeft;
points_left(:, 1) = 2000;

plot3(points_left(:, 1), points_left(:, 2), points_left(:, 3), '.', 'Markersize', markerSize, 'Color', [0.4 0.69 0.90]);

for i = 1:size(goodChannelLeft, 1)
    text(points_left(i, 1)+80, points_left(i, 2), points_left(i, 3), num2str(i), 'Color',  [0.46 0.67 0.18],...
        'FontWeight', 'bold', 'FontSize', frontSize);
    
    if ifInDCN(goodChannelLeft(i, :), LICoords)
        plot3(points_left(i, 1), points_left(i, 2), points_left(i, 3), '.', 'Markersize', 20, 'Color', [1 0.4 0.4]);
    else
        plot3(points_left(i, 1), points_left(i, 2), points_left(i, 3), '.', 'Markersize', 20, 'Color', [0.65 0.65 0.65]);
    end
    
%     if isStepModulatedLeft(i)
%         plot3(points_left(i, 1), points_left(i, 2), points_left(i, 3), '.', 'Markersize', 20, 'Color', [0.8 0.59 0.93]);
%     else
%         plot3(points_left(i, 1), points_left(i, 2), points_left(i, 3), '.', 'Markersize', 20, 'Color', [0.14 0.14 0.14]);
%     end
%     
    
    
end



LI_dorsalEdgePoint = [LIEdgeMedial, uniqueLICoords, LIEdgeDorsal];
LI_ventralEdgePoint = [LIEdgeMedial, uniqueLICoords, LIEdgeVentral];
plot3(LI_dorsalEdgePoint(:, 1), LI_dorsalEdgePoint(:, 2), LI_dorsalEdgePoint(:, 3), '-', 'Color', [0.96, 0.5, 0.3, 0.5], 'LineWidth', 2);
plot3(LI_ventralEdgePoint(:, 1), LI_ventralEdgePoint(:, 2), LI_ventralEdgePoint(:, 3), '-', 'Color', [0.96, 0.5, 0.3, 0.5], 'LineWidth', 2);

view(270, 0)
daspect([1 1 1])
title(['Good Units Location ( ' mouseID ' leftside view )']);






% side view, right side
h4 = subplot(2, 2, 4);
set(h4, 'position', [0.53, 0.11, 0.4, 0.4]);


% plotRegions3D(LICoords, 20, [1, 0.74, 0.35], xrange, yrange, 0.1);

plotRegions3D(RICoords, 20, [1, 0.74, 0.35], xrange, yrange, 0.1);
hold on


points_right = goodChannelRight;
points_right(:, 1) = 7700;

plot3(points_right(:, 1), points_right(:, 2), points_right(:, 3), '.', 'Markersize', markerSize, 'Color', [0.4 0.69 0.90]);


for i = 1:size(goodChannelRight, 1)
    text(points_right(i, 1)+80, points_right(i, 2), points_right(i, 3), num2str(i), 'Color',  [0.4 0.69 0.90],...
        'FontWeight', 'bold',  'FontSize', frontSize);
    if ifInDCN(goodChannelRight(i, :), RICoords)
        plot3(points_right(i, 1), points_right(i, 2), points_right(i, 3), '.', 'Markersize', 20, 'Color', [1 0.4 0.4]);
    else
        plot3(points_right(i, 1), points_right(i, 2), points_right(i, 3), '.', 'Markersize', 20, 'Color', [0.65 0.65 0.65]);
    end
%     
%     if isStepModulatedRight(i)
%         plot3(points_right(i, 1), points_right(i, 2), points_right(i, 3), '.', 'Markersize', 20, 'Color',  [0.8 0.59 0.93]);
%     else
%         plot3(points_right(i, 1), points_right(i, 2), points_right(i, 3), '.', 'Markersize', 20, 'Color', [0.14 0.14 0.14]);
%     end
%     
    
end


RI_dorsalEdgePoint = [RIEdgeMedial, uniqueRICoords, RIEdgeDorsal];
RI_ventralEdgePoint = [RIEdgeMedial, uniqueRICoords, RIEdgeVentral];
plot3(RI_dorsalEdgePoint(:, 1), RI_dorsalEdgePoint(:, 2), RI_dorsalEdgePoint(:, 3), '-', 'Color', [0.96, 0.5, 0.3, 0.5], 'LineWidth', 2);
plot3(RI_ventralEdgePoint(:, 1), RI_ventralEdgePoint(:, 2), RI_ventralEdgePoint(:, 3), '-', 'Color', [0.96, 0.5, 0.3, 0.5], 'LineWidth', 2);


view(90, 0)
daspect([1 1 1])
title(['Good Units Location ( ' mouseID ' rightside view )']);







