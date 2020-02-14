
folder = uigetdir();
files = dir(fullfile(folder, '*.tif'));

% Right Side
display('Reformatting right side...');
RightDentate = reformatTiffStack(fullfile(folder, 'DentateRight.tif'));
RightInt = reformatTiffStack(fullfile(folder, 'InterpositusRight.tif'));

RightPCLayer = reformatTiffStack(fullfile(folder, 'PCLayerRight.tif'));
RightBrainSurface = reformatTiffStack(fullfile(folder, 'BrainSurfaceRight.tif'));

% Left Side
display('Reformatting left side...');
LeftInt = reformatTiffStack(fullfile(folder, 'InterpositusLeft.tif'));
LeftDentate = reformatTiffStack(fullfile(folder, 'DentateLeft.tif'));
% LeftFastigial = reformatTiffStack(fullfile(folder, 'FastigialLeft.tif'));
LeftPCLayer = reformatTiffStack(fullfile(folder, 'PCLayerLeftLess.tif'));
LeftBrainSurface = reformatTiffStack(fullfile(folder, 'BrainSurfaceLeft.tif'));

% probes
display('Reformatting probe traces...');
LeftProbe_1L = reformatTiffStack(fullfile(folder, 'ProbeLeft_1_LeftShank.tif'));
LeftProbe_1R = reformatTiffStack(fullfile(folder, 'ProbeLeft_1_RightShank.tif'));
LeftProbe_2L = reformatTiffStack(fullfile(folder, 'ProbeLeft_2_LeftShank.tif'));
LeftProbe_2R = reformatTiffStack(fullfile(folder, 'ProbeLeft_2_RightShank.tif'));

RightProbe_1L = reformatTiffStack(fullfile(folder, 'ProbeRight_1_LeftShank.tif'));
RightProbe_1R = reformatTiffStack(fullfile(folder, 'ProbeRight_1_RightShank.tif'));

display('All Done!');

%% Get xyzCoords 

thickness = 60;

RDCoords = getCoordinates(RightDentate, thickness);
RICoords = getCoordinates(RightInt, thickness);
RPCLCoords = getCoordinates(RightPCLayer, thickness);
RBSCoords = getCoordinates(RightBrainSurface, thickness);

RPCoords_1L = getCoordinates(RightProbe_1L, thickness);
RPCoords_1R = getCoordinates(RightProbe_1R, thickness);

LICoords = getCoordinates(LeftInt, thickness);
LDCoords = getCoordinates(LeftDentate, thickness);
% LFCoords = getCoordinates(LeftFastigial, thickness, thickness);
LPCLCoords = getCoordinates(LeftPCLayer, thickness);
LBSCoords = getCoordinates(LeftBrainSurface, thickness);

LPCoords_1L = getCoordinates(LeftProbe_1L, thickness);
LPCoords_1R = getCoordinates(LeftProbe_1R, thickness);
LPCoords_2L = getCoordinates(LeftProbe_2L, thickness);
LPCoords_2R = getCoordinates(LeftProbe_2R, thickness);


%% Plot brain regions + probe traces

figure('Color', 'white', 'position', get(0,'ScreenSize'));
box off
axis off

xrange = max(LICoords(:, 1))*2;
yrange = max(LDCoords(:, 3))+500; % +60 just for visualizing purposes

% Left side 
plotRegions3D(LDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange);
hold on;
plotRegions3D(LICoords, 20, [1, 0.74, 0.35], xrange, yrange);
% hold on 
% plotRegions3D(LFCoords, 20, [0.39, 0.83, 0.075], xrange, yrange);

hold on;
plotRegions3D(LPCLCoords, 10, [0.82, 0.56, 0.97], xrange, yrange, 0.1);
hold on;
plotRegions3D(LBSCoords, 10, [0.8 0.8 0.8], xrange, yrange, 0.1);

hold on;
plotRegions3D(LPCoords_1L, 10, [1.0, 0.43, 0.54], xrange, yrange);
hold on;
plotRegions3D(LPCoords_1R, 10, [1.0, 0.43, 0.54], xrange, yrange);
hold on;
plotRegions3D(LPCoords_2L, 10, [1.0, 0.43, 0.54], xrange, yrange);
hold on;
plotRegions3D(LPCoords_2R, 10, [1.0, 0.43, 0.54], xrange, yrange);

% Right side
plotRegions3D(RDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange);
hold on;
plotRegions3D(RICoords, 20, [1, 0.74, 0.35], xrange, yrange);
hold on;
plotRegions3D(RPCLCoords, 10, [0.82, 0.56, 0.97], xrange, yrange, 0.1);
hold on;
plotRegions3D(RBSCoords, 10, [0.8 0.8 0.8], xrange, yrange, 0.1);
hold on;
plotRegions3D(RPCoords_1L, 10, [1.0, 0.43, 0.54], xrange, yrange);
hold on;
plotRegions3D(RPCoords_1R, 10, [1.0, 0.43, 0.54], xrange, yrange);


title('3D Plot of Traced Features');
% legend('Dentate', 'Interpositus', 'Fastigial');
%% fit a line for the probe 

LM_Left_1L = LinearFit(LPCoords_1L);
LM_Left_1R = LinearFit(LPCoords_1R);
LM_Left_2L = LinearFit(LPCoords_2L);
LM_Left_2R = LinearFit(LPCoords_2R);
LM_Right_1L = LinearFit(RPCoords_1L);
LM_Right_1R = LinearFit(RPCoords_1R);
% legend('Dentate Nucleus', 'Purkinje Cell Layer', 'Brain Surface', ' ','Probe Traces', 'Fitted Probe Tracks', 'Location' , 'northeast');
% title('3D Plot of Traced Features with Fitted Probe Tracks'); 



%% For Probes with DiI

% probe Left 1L 
% BS crossPoint
avg = LM_Left_1L.avg;
dirV = LM_Left_1L.dirVect;
dirV = -dirV;
plot3(avg(:, 1), avg(:, 3), avg(:, 2), '.b', 'MarkerSize', 30)

d_BS = -1033; % Has to be manual determined!!
BS_crossPoint_left1L = (d_BS)*dirV + avg;
hold on
plot3(BS_crossPoint_left1L(:, 1), BS_crossPoint_left1L(:, 3), BS_crossPoint_left1L(:, 2), '.k', 'MarkerSize', 30)
% points = [BS_crossPoint_left1L ; avg];
% plot3(points(:,1),points(:,3),points(:,2),'-k','LineWidth',3) 


% PC cross point
offset = 70; % Has to be manual determined!!

% all the distance has to be manually calculated
d_PC1 = 545 + offset;
PC_crossPoint1_left1L = (d_PC1)*dirV + BS_crossPoint_left1L;
hold on
plot3(PC_crossPoint1_left1L(:, 1), PC_crossPoint1_left1L(:, 3), PC_crossPoint1_left1L(:, 2), '.r', 'MarkerSize', 30)

d_PC2 = 1604 + offset;
PC_crossPoint2_left1L = (d_PC2)*dirV + BS_crossPoint_left1L;
hold on
plot3(PC_crossPoint2_left1L(:, 1), PC_crossPoint2_left1L(:, 3), PC_crossPoint2_left1L(:, 2), '.r', 'MarkerSize', 30)



% Good Channels
% all the distance has to be manually calculated
d_GC1 = 1954.5 + offset;
goodChannel_Left1L_1 = (d_GC1)*dirV + BS_crossPoint_left1L;
hold on
plot3(goodChannel_Left1L_1(:, 1), goodChannel_Left1L_1(:, 3), goodChannel_Left1L_1(:, 2), '.c', 'Markersize', 30);

d_GC2 = 1854.5 + offset;
goodChannel_Left1L_2 = (d_GC2)*dirV + BS_crossPoint_left1L;
hold on
plot3(goodChannel_Left1L_2(:, 1), goodChannel_Left1L_2(:, 3), goodChannel_Left1L_2(:, 2), '.c', 'Markersize', 30);


d_GC3 = 1879.5 + offset;
goodChannel_Left1L_3 = (d_GC3)*dirV + BS_crossPoint_left1L;
hold on
plot3(goodChannel_Left1L_3(:, 1), goodChannel_Left1L_3(:, 3), goodChannel_Left1L_3(:, 2), '.c', 'Markersize', 30);


% Probe Left iR
% BS cross point
avg = LM_Left_1R.avg;
dirV = LM_Left_1R.dirVect;
dirV = -dirV;
plot3(avg(:, 1), avg(:, 3), avg(:, 2), '.b', 'MarkerSize', 30)

d_BS = -933.33; % Has to be manually determined!!
BS_crossPoint_left1R = (d_BS)*dirV + avg;
hold on
plot3(BS_crossPoint_left1R(:, 1), BS_crossPoint_left1R(:, 3), BS_crossPoint_left1R(:, 2), '.k', 'MarkerSize', 30)
% points = [BS_crossPoint_left1R ; avg];
% plot3(points(:,1),points(:,3),points(:,2),'-k','LineWidth',3) 


d_PC1 = 1360 + offset;
PC_crossPoint1_left1R = (d_PC1)*dirV + BS_crossPoint_left1R;
hold on
plot3(PC_crossPoint1_left1R(:, 1), PC_crossPoint1_left1R(:, 3), PC_crossPoint1_left1R(:, 2), '.r', 'MarkerSize', 30)

% d_PC2 = distanceConverter(1074.9 + offset, 0.9, 0.3);
% PC_crossPoint2_left1R = (d_PC2)*dirV + BS_crossPoint_left1R;
% hold on
% plot3(PC_crossPoint2_left1R(:, 1), PC_crossPoint2_left1R(:, 3), PC_crossPoint2_left1R(:, 2), '.r', 'MarkerSize', 30)
% 


% % Good Channels
% d_GC1 = distanceConverter(1954.5 + offset, 0.9, 0.3);
% goodChannel_Left1R_1 = (d_GC1)*dirV + BS_crossPoint_left1R;
% hold on
% plot3(goodChannel_Left1R_1(:, 1), goodChannel_Left1R_1(:, 3), goodChannel_Left1R_1(:, 2), '.c', 'Markersize', 30);
% 
% d_GC2 = distanceConverter(1854.5 + offset, 0.9, 0.3);
% goodChannel_Left1R_2 = (d_GC2)*dirV + BS_crossPoint_left1R;
% hold on
% plot3(goodChannel_Left1R_2(:, 1), goodChannel_Left1R_2(:, 3), goodChannel_Left1R_2(:, 2), '.c', 'Markersize', 30);
% 
% d_GC3 = distanceConverter(1879.5 + offset, 0.9, 0.3);
% goodChannel_Left1R_3 = (d_GC3)*dirV + BS_crossPoint_left1R;
% hold on
% plot3(goodChannel_Left1R_3(:, 1), goodChannel_Left1R_3(:, 3), goodChannel_Left1R_3(:, 2), '.c', 'Markersize', 30);
% 
% 
% points = [goodChannel_Left1R_1 ; avg];
% plot3(points(:,1),points(:,3),points(:,2),'-g','LineWidth',3) 


%% Probe Left 2

% probe left 2L
% BS cross point


avg = LM_Left_2L.avg;
dirV = LM_Left_2L.dirVect;
dirV = -dirV;
plot3(avg(:, 1), avg(:, 3), avg(:, 2), '.b', 'MarkerSize', 30)

d_BS = -1220; % Has to be manually determined!!
BS_crossPoint_left2L = (d_BS)*dirV + avg;
hold on
plot3(BS_crossPoint_left2L(:, 1), BS_crossPoint_left2L(:, 3), BS_crossPoint_left2L(:, 2), '.k', 'MarkerSize', 30)
points = [BS_crossPoint_left2L ; avg];
plot3(points(:,1),points(:,3),points(:,2),'-k','LineWidth',3) 


% PC cross point
offset = 350; % Has to be manual determined!!

% all the distance has to be manually calculated
d_PC1 = 1550 + offset;
PC_crossPoint1_left2L = (d_PC1)*dirV + BS_crossPoint_left2L;
hold on
plot3(PC_crossPoint1_left2L(:, 1), PC_crossPoint1_left2L(:, 3), PC_crossPoint1_left2L(:, 2), '.r', 'MarkerSize', 30)

% d_PC2 = distanceConverter(665 + offset, 0.9, 0.3);
% PC_crossPoint2_left2L = (d_PC2)*dirV + BS_crossPoint_left2L;
% hold on
% plot3(PC_crossPoint2_left2L(:, 1), PC_crossPoint2_left2L(:, 3), PC_crossPoint2_left2L(:, 2), '.r', 'MarkerSize', 30)

% Good Channels
% all the distance has to be manually calculated
d_GC1 = 1824.4 + offset;
goodChannel_Left2L_1 = (d_GC1)*dirV + BS_crossPoint_left2L;
hold on
plot3(goodChannel_Left2L_1(:, 1), goodChannel_Left2L_1(:, 3), goodChannel_Left2L_1(:, 2), '.c', 'Markersize', 30);

d_GC2 = 2174.4 + offset;
goodChannel_Left2L_2 = (d_GC2)*dirV + BS_crossPoint_left2L;
hold on
plot3(goodChannel_Left2L_2(:, 1), goodChannel_Left2L_2(:, 3), goodChannel_Left2L_2(:, 2), '.c', 'Markersize', 30);



% Probe Left iR
% BS cross point
avg = LM_Left_2R.avg;
dirV = LM_Left_2R.dirVect;
plot3(avg(:, 1), avg(:, 3), avg(:, 2), '.b', 'MarkerSize', 30)

d_BS = -1366.67; % Has to be manually determined!!
BS_crossPoint_left2R = (d_BS)*dirV + avg;
hold on
plot3(BS_crossPoint_left2R(:, 1), BS_crossPoint_left2R(:, 3), BS_crossPoint_left2R(:, 2), '.k', 'MarkerSize', 30)
points = [BS_crossPoint_left2R ; avg];
plot3(points(:,1),points(:,3),points(:,2),'-k','LineWidth',3) 


d_PC1 = 665 + offset;
PC_crossPoint1_left2R = (d_PC1)*dirV + BS_crossPoint_left2R;
hold on
plot3(PC_crossPoint1_left2R(:, 1), PC_crossPoint1_left2R(:, 3), PC_crossPoint1_left2R(:, 2), '.r', 'MarkerSize', 30)

d_PC2 = 1315 + offset;
PC_crossPoint2_left2R = (d_PC2)*dirV + BS_crossPoint_left2R;
hold on
plot3(PC_crossPoint2_left2R(:, 1), PC_crossPoint2_left2R(:, 3), PC_crossPoint2_left2R(:, 2), '.r', 'MarkerSize', 30)


% Good Channels
d_GC1 = 1874.4 + offset;
goodChannel_Left2R_1 = (d_GC1)*dirV + BS_crossPoint_left2R;
hold on
plot3(goodChannel_Left2R_1(:, 1), goodChannel_Left2R_1(:, 3), goodChannel_Left2R_1(:, 2), '.c', 'Markersize', 30);

% d_GC2 = distanceConverter(2149.4 + offset, 0.9, 0.3);
% goodChannel_Left2L_2 = (d_GC2)*dirV + BS_crossPoint_left2L;
% hold on
% plot3(goodChannel_Left2L_2(:, 1), goodChannel_Left2L_2(:, 3), goodChannel_Left2L_2(:, 2), '.c', 'Markersize', 30);

%% Add points for Probe Right 1

% Probe right 1L
% BS crossPoint
avg = LM_Right_1L.avg;
dirV = LM_Right_1L.dirVect;
dirV = -dirV;
plot3(avg(:, 1), avg(:, 3), avg(:, 2), '.b', 'MarkerSize', 30)

d_BS = -686.67;
BS_crossPoint_right1L = (d_BS)*dirV + avg;
hold on
plot3(BS_crossPoint_right1L(:, 1), BS_crossPoint_right1L(:, 3), BS_crossPoint_right1L(:, 2), '.k', 'MarkerSize', 30)
points = [BS_crossPoint_right1L ; avg];
plot3(points(:,1),points(:,3),points(:,2),'-k','LineWidth',3) 



% PC crossPoint
offset = 150; % Has to be manual determined!!

% all the distance has to be manually calculated
d_PC1 = 1154.3 + offset;
PC_crossPoint1_right1L = (d_PC1)*dirV + BS_crossPoint_right1L;
hold on
plot3(PC_crossPoint1_right1L(:, 1), PC_crossPoint1_right1L(:, 3), PC_crossPoint1_right1L(:, 2), '.r', 'MarkerSize', 30)

d_PC2 = 205 + offset;
PC_crossPoint2_right1L = (d_PC2)*dirV + BS_crossPoint_right1L;
hold on
plot3(PC_crossPoint2_right1L(:, 1), PC_crossPoint2_right1L(:, 3), PC_crossPoint2_right1L(:, 2), '.r', 'MarkerSize', 30)


% Good Channels
d_GC1 = 2016.8 + offset;
goodChannel_right1L_1 = (d_GC1)*dirV + BS_crossPoint_right1L;
hold on
plot3(goodChannel_right1L_1(:, 1), goodChannel_right1L_1(:, 3), goodChannel_right1L_1(:, 2), '.c', 'Markersize', 30);
points = [goodChannel_right1L_1 ; avg];
plot3(points(:,1),points(:,3),points(:,2),'-g','LineWidth',3) 

d_GC2 = 1566.8 + offset;
goodChannel_right1L_2 = (d_GC2)*dirV + BS_crossPoint_right1L;
hold on
plot3(goodChannel_right1L_2(:, 1), goodChannel_right1L_2(:, 3), goodChannel_right1L_2(:, 2), '.c', 'Markersize', 30);
points = [goodChannel_right1L_2 ; avg];
plot3(points(:,1),points(:,3),points(:,2),'-g','LineWidth',3) 



% Probe right 1R
% BS crossPoint 
avg = LM_Right_1R.avg;
dirV = LM_Right_1R.dirVect;
dirV = -dirV;
plot3(avg(:, 1), avg(:, 3), avg(:, 2), '.b', 'MarkerSize', 30)

d_BS = -800;
BS_crossPoint_right1R = (d_BS)*dirV + avg;
hold on
plot3(BS_crossPoint_right1R(:, 1), BS_crossPoint_right1R(:, 3), BS_crossPoint_right1R(:, 2), '.k', 'MarkerSize', 30)
points = [BS_crossPoint_right1R ; avg];
plot3(points(:,1),points(:,3),points(:,2),'-k','LineWidth',3) 


% all the distance has to be manually calculated
d_PC1 = 969.6 + offset;
PC_crossPoint1_right1R = (d_PC1)*dirV + BS_crossPoint_right1R;
hold on
plot3(PC_crossPoint1_right1R(:, 1), PC_crossPoint1_right1R(:, 3), PC_crossPoint1_right1R(:, 2), '.r', 'MarkerSize', 30)

% d_PC2 = distanceConverter(1247.1 + offset, 1, 0.3);
% PC_crossPoint2_right1R = (d_PC2)*dirV + BS_crossPoint_right1R;
% hold on
% plot3(PC_crossPoint2_right1R(:, 1), PC_crossPoint2_right1R(:, 3), PC_crossPoint2_right1R(:, 2), '.r', 'MarkerSize', 30)
% 




%% Backup codes

% Start Channel Point
d_startPoint = distanceConverter(140, 1, 0.3);
startPoint = (d_startPoint+d_PC)*dirV + avg;
hold on
plot3(startPoint(:, 1), startPoint(:, 3), startPoint(:, 2), '.m', 'MarkerSize', 20)


% End Channel Point
d_endPoint = distanceConverter(790, 1, 0.3);
endPoint = (d_PC + d_startPoint + d_endPoint) * dirV + avg;
hold on
points = [startPoint ; endPoint];
plot3(points(:,1),points(:,3),points(:,2),'-m','LineWidth',3) 
plot3(endPoint(:, 1), endPoint(:, 3), endPoint(:, 2), '.m', 'MarkerSize', 20)


% Good Channel Point
goodChannelNum = [63, 58, 61, 49, 47];
d_channels =  (goodChannelNum-32).*25;
d_channels = distanceConverter(d_channels, 1, 0.3);
for i = 1:length(goodChannelNum)
    channelPoint = (d_channels(i) + d_PC + d_startPoint)*dirV + avg; 
    hold on
    plot3(channelPoint(:, 1), channelPoint(:, 3), channelPoint(:, 2), '.c', 'Markersize', 30);
end








