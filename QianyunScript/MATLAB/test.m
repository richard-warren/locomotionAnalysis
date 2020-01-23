%% Test if we can directly apply linear regression without converting systems

histoCoords = [BS_crossPoint_left1L; BS_crossPoint_left2; BS_crossPoint_left3; BS_crossPoint_right1];
histoCoords = [histoCoords(:, 1) histoCoords(:, 3)];

manipuCoords = [-2.48, 6; -1.64, 6.16; -1.56 6.35; 2.17 6.19];

testmdl1 = fitlm([manipuCoords(:, 1) manipuCoords(:, 2)], histoCoords(:, 1))
testmdl2 = fitlm([manipuCoords(:, 1) manipuCoords(:, 2)], histoCoords(:, 2))

%% test - direct linear model for entry point, hold one, fit three


histoCoords = [BS_crossPoint_left1L; BS_crossPoint_left2; BS_crossPoint_left3; BS_crossPoint_right1];
histoCoords = [histoCoords(:, 1) histoCoords(:, 3)];

manipuCoords = [-2.48, 6; -1.64, 6.16; -1.56 6.35; 2.17 6.19];
fittedEntryPoint = nan(4, 3);

for i = 1:4

% Plot entry points - original data
figure;
plot(histoCoords(i, 1), histoCoords(i, 2), '.', 'MarkerSize', 30, 'Color', [0.49, 0.84, 0.97]);
hold on
plot(manipuCoords(i, 1), manipuCoords(i, 2), '.r', 'MarkerSize', 10);
xlabel('ML Axis');
ylabel('AP Axis');
title(['fitted point to probe traces', num2str(i)]);

tempManipu = manipuCoords;
tempHisto = histoCoords;
tempManipu(i, :) = [];
tempHisto(i, :) = [];
hold on
plot(tempManipu(:, 1), tempManipu(:, 2), '.m', 'MarkerSize', 10);
hold on
plot(tempHisto(:, 1), tempHisto(:, 2), '.b', 'MarkerSize', 20);


% build linear model using other probe traces
mdlML = fitlm([tempManipu(:, 1) tempManipu(:, 2)], tempHisto(:, 1))
fittedML = predict(mdlML, [manipuCoords(i, 1), manipuCoords(i, 2)])

mdlAP = fitlm([tempManipu(:, 1) tempManipu(:, 2)], tempHisto(:, 2))
fittedAP = predict(mdlAP, [manipuCoords(i, 1), manipuCoords(i, 2)])

fittedEntryPoint(i, :) = [fittedML, nan, fittedAP];

hold on
plot(fittedML, fittedAP, '.', 'MarkerSize', 30, 'Color', [0.84 0.55 0.73])
end

%% Test - generating fake entry points for no dii traces (topdown view)

histoCoords = [BS_crossPoint_left1L; BS_crossPoint_left2; BS_crossPoint_left3; BS_crossPoint_right1];
histoCoords = [histoCoords(:, 1) histoCoords(:, 3)];

manipuCoords = [-2.48, 6; -1.64, 6.16; -1.56 6.35; 2.17 6.19];
manipuCoords_noDiI = [1.5 6; 1.75 6; 2.11 6.11; 2.36 6.11];

plot(manipuCoords(:, 1), manipuCoords(:, 2), '.r', 'MarkerSize', 10);
hold on
plot(manipuCoords_noDiI(:, 1), manipuCoords_noDiI(:, 2), '.m', 'MarkerSize', 10);
hold on
plot(histoCoords(:, 1), histoCoords(:, 2), '.b', 'MarkerSize', 20);

mdlML_test = fitlm([manipuCoords(:, 1) manipuCoords(:, 2)], histoCoords(:, 1))
mdlAP_test = fitlm([manipuCoords(:, 1) manipuCoords(:, 2)], histoCoords(:, 2))


fittedML_nodii = predict(mdlML_test, [manipuCoords_noDiI(:, 1), manipuCoords_noDiI(:, 2)])
fittedAP_nodii = predict(mdlAP_test, [manipuCoords_noDiI(:, 1), manipuCoords_noDiI(:, 2)])

hold on
plot(fittedML_nodii, fittedAP_nodii, '.', 'MarkerSize', 20, 'Color', [1 0.74 0.35])

%% Test - generating fake probe traces for dii tracks. Hold one, fit three
figure;
xrange = size(LeftDentate, 2);
yrange = max(LDCoords(:, 3))+60; % +60 just for visualizing purposes
% 
PlotRegions3D(LDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange);
hold on;
PlotRegions3D(LICoords, 20, [1, 0.74, 0.35], xrange, yrange);

hold on;
PlotRegions3D(LPCLCoords, 10, [0.82, 0.56, 0.97], xrange, yrange);
hold on;
PlotRegions3D(LBSCoords, 10, [0.8 0.8 0.8], xrange, yrange);

hold on;
PlotRegions3D(LPCoords_1L, 10, [1.0, 0.43, 0.54], xrange, yrange);
hold on;
PlotRegions3D(LPCoords_1R, 10, [1.0, 0.43, 0.54], xrange, yrange);
hold on;
PlotRegions3D(LPCoords_2, 10, [1.0, 0.43, 0.54], xrange, yrange);
hold on;
PlotRegions3D(LPCoords_3, 10, [1.0, 0.43, 0.54], xrange, yrange);


PlotRegions3D(RDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange);
hold on;
PlotRegions3D(RICoords, 20, [1, 0.74, 0.35], xrange, yrange);
hold on;
PlotRegions3D(RPCLCoords, 10, [0.82, 0.56, 0.97], xrange, yrange);
hold on;
PlotRegions3D(RBSCoords, 10, [0.8 0.8 0.8], xrange, yrange);
hold on;
PlotRegions3D(RPCoords, 10, [1.0, 0.43, 0.54], xrange, yrange);



LM_Left_1L = LinearFit(LPCoords_1L);
LM_Left_1R = LinearFit(LPCoords_1R);
LM_Left_2 = LinearFit(LPCoords_2);
LM_Left_3 = LinearFit(LPCoords_3);
LM_Right_1 = LinearFit(RPCoords);



dirV = nan(4, 3);
dirV(1, :) = LM_Left_1L.dirVect;
dirV(2, :) = LM_Left_2.dirVect;
dirV(3, :) = LM_Left_3.dirVect;
dirV(4, :) = LM_Right_1.dirVect;



histoCoords = [BS_crossPoint_left1L; BS_crossPoint_left2; BS_crossPoint_left3; BS_crossPoint_right1];

hold on
plot3(histoCoords(:, 1), histoCoords(:, 3), histoCoords(:, 2), '.', 'MarkerSize', 30, 'Color', [0.6 0.33 0.1])
hold on
fittedEntryPoint(1:4, 2) = [158; 137; 136; 242.1994]; % get the depth for entry point manually
plot3(fittedEntryPoint(:, 1), fittedEntryPoint(:, 3), fittedEntryPoint(:, 2), '.', 'MarkerSize', 30, 'Color', [0.98 0.79 0.63])

manipuDepth = distanceConverter([2843, 2559, 2732, 2050] );

fittedEndPoint = nan(4, 3);

i = 1;
    
tempdirV = dirV;
tempdirV(i, :) = [];
meandirV = mean(tempdirV);
distance = manipuDepth(i);
fittedEndPoint(i, :) = fittedEntryPoint(i, :) + meandirV*distance;
points = [fittedEndPoint(i, :); fittedEntryPoint(i, :)];
hold on
plot3(points(:,1),points(:,3),points(:,2), '-', 'LineWidth', 2, 'Color', [1 0.42 0.21]);
PC_distance = distanceConverter(1990);
PC_crosspoint = fittedEntryPoint(i, :) + meandirV*PC_distance;
hold on
plot3(PC_crosspoint(:, 1), PC_crosspoint(:, 3), PC_crosspoint(:, 2), '.r', 'MarkerSize', 30);


% plot channel position
d_startPoint = distanceConverter(140, 1, 0.3);
startPoint = d_startPoint*meandirV + PC_crosspoint;
hold on
plot3(startPoint(:, 1), startPoint(:, 3), startPoint(:, 2), '.m', 'MarkerSize', 20)


% End Channel Point
d_endPoint = distanceConverter(790, 1, 0.3);
endPoint = (d_startPoint + d_endPoint) * meandirV + PC_crosspoint;
hold on
points = [startPoint ; endPoint];
plot3(points(:,1),points(:,3),points(:,2),'-m','LineWidth',3) 
plot3(endPoint(:, 1), endPoint(:, 3), endPoint(:, 2), '.m', 'MarkerSize', 20)


% Good Channel Point
goodChannelNum = [63, 58, 61, 49, 47];
d_channels =  (goodChannelNum-32).*25;
d_channels = distanceConverter(d_channels, 1, 0.3);
for i = 1:length(goodChannelNum)
    channelPoint = (d_channels(i) + d_startPoint)*meandirV + PC_crosspoint; 
    hold on
    plot3(channelPoint(:, 1), channelPoint(:, 3), channelPoint(:, 2), '.c', 'Markersize', 30);
end

%% Test - generating fake probe traces for no dii tracks (3d view)

figure;
xrange = size(LeftDentate, 2);
yrange = max(LDCoords(:, 3))+60; % +60 just for visualizing purposes
% 
PlotRegions3D(LDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange);
hold on;
PlotRegions3D(LICoords, 20, [1, 0.74, 0.35], xrange, yrange);

hold on;
PlotRegions3D(LPCLCoords, 10, [0.8471 0.7686 1.0000], xrange, yrange);
hold on;
PlotRegions3D(LBSCoords, 10, [0.8 0.8 0.8], xrange, yrange);

hold on;
PlotRegions3D(LPCoords_1L, 10, [1.0, 0.43, 0.54], xrange, yrange);
hold on;
PlotRegions3D(LPCoords_1R, 10, [1.0, 0.43, 0.54], xrange, yrange);
hold on;
PlotRegions3D(LPCoords_2, 10, [1.0, 0.43, 0.54], xrange, yrange);
hold on;
PlotRegions3D(LPCoords_3, 10, [1.0, 0.43, 0.54], xrange, yrange);


PlotRegions3D(RDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange);
hold on;
PlotRegions3D(RICoords, 20, [1, 0.74, 0.35], xrange, yrange);
hold on;
PlotRegions3D(RPCLCoords, 10, [0.8471 0.7686 1.0000], xrange, yrange);
hold on;
PlotRegions3D(RBSCoords, 10, [0.8 0.8 0.8], xrange, yrange);
hold on;
PlotRegions3D(RPCoords, 10, [1.0, 0.43, 0.54], xrange, yrange);



LM_Left_1L = LinearFit(LPCoords_1L);
LM_Left_1R = LinearFit(LPCoords_1R);
LM_Left_2 = LinearFit(LPCoords_2);
LM_Left_3 = LinearFit(LPCoords_3);
LM_Right_1 = LinearFit(RPCoords);



dirV = nan(4, 3);
dirV(1, :) = LM_Left_1L.dirVect;
dirV(2, :) = LM_Left_2.dirVect;
dirV(3, :) = LM_Left_3.dirVect;
dirV(4, :) = LM_Right_1.dirVect;
meandirV = mean(dirV);

histoCoords = [BS_crossPoint_left1L; BS_crossPoint_left2; BS_crossPoint_left3; BS_crossPoint_right1];

hold on
plot3(histoCoords(:, 1), histoCoords(:, 3), histoCoords(:, 2), '.', 'MarkerSize', 30, 'Color', [0.6 0.33 0.1])
hold on
fittedEntryPoint(5:8, 1) = fittedML_nodii;
fittedEntryPoint(5:8, 3) = fittedAP_nodii;
fittedEntryPoint(:, 2) = [158; 137; 136; 242.1994; 168; 183; 250; 275]; % get the depth for entry point manually
plot3(fittedEntryPoint(:, 1), fittedEntryPoint(:, 3), fittedEntryPoint(:, 2), '.', 'MarkerSize', 30, 'Color', [0.98 0.79 0.63])

manipuDepth_nodii = distanceConverter([2700, 2700, 2100, 2100]);

for i = 1:4
    
   fittedEndPoint(i+4, :) = fittedEntryPoint(i+4, :) + meandirV * manipuDepth_nodii(i); 
   points = [fittedEndPoint(i+4, :); fittedEntryPoint(i+4, :)];
   hold on
   plot3(points(:,1),points(:,3),points(:,2), '-', 'LineWidth', 2, 'Color', [0.478 1 0.827]);
   
end


%% Adjustin fake probe traces for Right 1R

% plot the end point
hold on
plot3(fittedEndPoint(6, 1), fittedEndPoint(6, 3), fittedEndPoint(6, 2), '.m', 'MarkerSize', 30);

% plot the PC layer cross points
PC_distance_right1R = distanceConverter(-675, 0.9, 0.3); % reference to end point
PC_distance_right1R2 = distanceConverter(-1305, 0.9, 0.3);
PC_crossPoint_right1R = fittedEndPoint(6, :) + meandirV * PC_distance_right1R;
PC_crossPoint_right1R2 = fittedEndPoint(6, :) + meandirV * PC_distance_right1R2;
plot3(PC_crossPoint_right1R(:, 1), PC_crossPoint_right1R(:, 3), PC_crossPoint_right1R(:, 2), '.r', 'MarkerSize', 30);
plot3(PC_crossPoint_right1R2(:, 1), PC_crossPoint_right1R2(:, 3), PC_crossPoint_right1R2(:, 2), '.r', 'MarkerSize', 30);

% plot good channel locations
goodChannel_right1R = [26 27 29 31 32];
distance_right1R = distanceConverter((32 - goodChannel_right1R) * (-25));

goodChannelPoint_right1R = nan(length(goodChannel_right1R), 3);
for i = 1:length(goodChannel_right1R)
    
   goodChannelPoint_right1R(i, :) = fittedEndPoint(6, :) + meandirV * distance_right1R(i);
    
end

hold on 
plot3(goodChannelPoint_right1R(:, 1), goodChannelPoint_right1R(:, 3), goodChannelPoint_right1R(:, 2), '.c', 'MarkerSize', 30);



%% Adjustin fake probe traces for Right 1L

hold on
plot3(fittedEndPoint(5, 1), fittedEndPoint(5, 3), fittedEndPoint(5, 2), '.m', 'MarkerSize', 30);


% plot the PC layer cross points
PC_distance_right1L = distanceConverter(-985, 0.9, 0.3); % reference to end point
PC_crossPoint_right1L = fittedEndPoint(5, :) + meandirV * PC_distance_right1R;

plot3(PC_crossPoint_right1L(:, 1), PC_crossPoint_right1L(:, 3), PC_crossPoint_right1L(:, 2), '.r', 'MarkerSize', 30);

% plot the good channel points
goodChannel_right1L = [43 55];
distance_right1L = distanceConverter((64 - goodChannel_right1L) * (-25));

goodChannelPoint_right1L = nan(length(goodChannel_right1L), 3);
for i = 1:length(goodChannel_right1L)
    
   goodChannelPoint_right1L(i, :) = fittedEndPoint(5, :) + meandirV * distance_right1L(i);
    
end

hold on 
plot3(goodChannelPoint_right1L(:, 1), goodChannelPoint_right1L(:, 3), goodChannelPoint_right1L(:, 2), '.c', 'MarkerSize', 30);

%% Adjustin fake probe traces for Right 2R


% plot the end point
hold on
plot3(fittedEndPoint(7, 1), fittedEndPoint(7, 3), fittedEndPoint(7, 2), '.m', 'MarkerSize', 30);

% plot the PC layer cross points
PC_distance_right2L = distanceConverter(-801, 0.9, 0.3); % reference to end point
PC_distance_right2L2 = distanceConverter(-1616, 0.9, 0.3);
PC_crossPoint_right2L = fittedEndPoint(7, :) + meandirV * PC_distance_right2L;
PC_crossPoint_right2L2 = fittedEndPoint(7, :) + meandirV * PC_distance_right2L2;
plot3(PC_crossPoint_right2L(:, 1), PC_crossPoint_right2L(:, 3), PC_crossPoint_right2L(:, 2), '.r', 'MarkerSize', 30);
plot3(PC_crossPoint_right2L2(:, 1), PC_crossPoint_right2L2(:, 3), PC_crossPoint_right2L2(:, 2), '.r', 'MarkerSize', 30);


% plot good channel locations
goodChannel_right2L = [56];
distance_right2L = distanceConverter((64 - goodChannel_right2L) * (-25));

goodChannelPoint_right2L = nan(length(goodChannel_right2L), 3);
for i = 1:length(goodChannel_right2L)
    
   goodChannelPoint_right2L(i, :) = fittedEndPoint(7, :) + meandirV * distance_right2L(i);
    
end

hold on 
plot3(goodChannelPoint_right2L(:, 1), goodChannelPoint_right2L(:, 3), goodChannelPoint_right2L(:, 2), '.c', 'MarkerSize', 30);


%% Adjustin fake probe traces for Right 2R

hold on
plot3(fittedEndPoint(8, 1), fittedEndPoint(8, 3), fittedEndPoint(8, 2), '.m', 'MarkerSize', 30);


% plot the PC layer cross points
PC_distance_right2R = distanceConverter(-925, 0.9, 0.3); % reference to end point
PC_crossPoint_right2R = fittedEndPoint(8, :) + meandirV * PC_distance_right2R;

plot3(PC_crossPoint_right2R(:, 1), PC_crossPoint_right2R(:, 3), PC_crossPoint_right2R(:, 2), '.r', 'MarkerSize', 30);

goodChannel_right2R = [9 10];
distance_right2R = distanceConverter((32 - goodChannel_right2R) * (-25));

goodChannelPoint_right2R = nan(length(goodChannel_right2R), 3);
for i = 1:length(goodChannel_right2R)
    
   goodChannelPoint_right2R(i, :) = fittedEndPoint(8, :) + meandirV * distance_right2R(i);
    
end

hold on 
plot3(goodChannelPoint_right2R(:, 1), goodChannelPoint_right2R(:, 3), goodChannelPoint_right2R(:, 2), '.c', 'MarkerSize', 30);


























