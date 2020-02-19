
folder = uigetdir();
files = dir(fullfile(folder, '*.tif'));

downsampleRate = 0.3;
thickness = 50;

% Right Side
disp('Reformatting right side...');
% RDCoords = importTiffStack(fullfile(folder, 'DentateRight.tif'), downsampleRate, thickness);
RICoords = importTiffStack(fullfile(folder, 'InterpositusRight.tif'), downsampleRate, thickness);

RPCLCoords = importTiffStack(fullfile(folder, 'PCLayerRight.tif'), downsampleRate, thickness);
RBSCoords = importTiffStack(fullfile(folder, 'BrainSurfaceRight.tif'), downsampleRate, thickness);

% Left Side
disp('Reformatting left side...');
LICoords = importTiffStack(fullfile(folder, 'InterpositusLeft.tif'), downsampleRate, thickness);
% LDCoords = importTiffStack(fullfile(folder, 'DentateLeft.tif'), downsampleRate, thickness);
% LeftFastigial = importTiffStack(fullfile(folder, 'FastigialLeft.tif'), downsampleRate, thickness);
LPCLCoords = importTiffStack(fullfile(folder, 'PCLayerLeft.tif'), downsampleRate, thickness);
LBSCoords = importTiffStack(fullfile(folder, 'BrainSurfaceLeft.tif'), downsampleRate, thickness);

% probes
disp('Reformatting probe traces...');
LeftProbe_1L = importTiffStack(fullfile(folder, 'ProbeLeft_1L.tif'), downsampleRate, thickness);
% LeftProbe_1M = importTiffStack(fullfile(folder, 'ProbeLeft_1M.tif'), downsampleRate, thickness);
LeftProbe_1R = importTiffStack(fullfile(folder, 'ProbeLeft_1R.tif'), downsampleRate, thickness);
LeftProbe_2L = importTiffStack(fullfile(folder, 'ProbeLeft_2L.tif'), downsampleRate, thickness);
% LeftProbe_2M = importTiffStack(fullfile(folder, 'ProbeLeft_2M.tif'), downsampleRate, thickness);
LeftProbe_2R = importTiffStack(fullfile(folder, 'ProbeLeft_2R.tif'), downsampleRate, thickness);

RightProbe_1L = importTiffStack(fullfile(folder, 'ProbeRight_1L.tif'), downsampleRate, thickness);
% RightProbe_1M = importTiffStack(fullfile(folder, 'ProbeRight_1M.tif'), downsampleRate, thickness);
RightProbe_1R = importTiffStack(fullfile(folder, 'ProbeRight_1R.tif'), downsampleRate, thickness);

disp('All Done!');



%% Plot brain regions + probe traces

figure('Color', 'white', 'position', get(0,'ScreenSize'));
box off
axis off

xrange = max(LICoords(:, 1))*2;
yrange = max(LICoords(:, 2))+500; % +60 just for visualizing purposes

% Left side 
% plotRegions3D(LDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange);
plotRegions3D(LICoords, 20, [1, 0.74, 0.35], xrange, yrange, 0.1);
% hold on 
% plotRegions3D(LFCoords, 20, [0.39, 0.83, 0.075], xrange, yrange);

hold on;
plotRegions3D(LPCLCoords, 10, [0.82, 0.56, 0.97], xrange, yrange, 0.1);
plotRegions3D(LBSCoords, 10, [0.8 0.8 0.8], xrange, yrange, 0.1);

plotRegions3D(LeftProbe_1L, 10, [1.0, 0.43, 0.54], xrange, yrange);
plotRegions3D(LeftProbe_1R, 10, [1.0, 0.43, 0.54], xrange, yrange);
plotRegions3D(LeftProbe_2L, 10, [1.0, 0.43, 0.54], xrange, yrange);
plotRegions3D(LeftProbe_2R, 10, [1.0, 0.43, 0.54], xrange, yrange);
% plotRegions3D(LeftProbe_2M, 10, [1.0, 0.43, 0.54], xrange, yrange);
% plotRegions3D(LeftProbe_1M, 10, [1.0, 0.43, 0.54], xrange, yrange);

% Right side
% plotRegions3D(RDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange);
hold on;
plotRegions3D(RICoords, 20, [1, 0.74, 0.35], xrange, yrange, 0.1);
plotRegions3D(RPCLCoords, 10, [0.82, 0.56, 0.97], xrange, yrange, 0.1);
plotRegions3D(RBSCoords, 10, [0.8 0.8 0.8], xrange, yrange, 0.1);

plotRegions3D(RightProbe_1L, 10, [1.0, 0.43, 0.54], xrange, yrange);
plotRegions3D(RightProbe_1R, 10, [1.0, 0.43, 0.54], xrange, yrange);
% plotRegions3D(RightProbe_1M, 10, [1.0, 0.43, 0.54], xrange, yrange);


title('3D Plot of Traced Features');
% legend('Dentate', 'Interpositus', 'Fastigial');
%% fit a line for the probe 

LM_Left_1L = LinearFit(LeftProbe_1L);
LM_Left_1R = LinearFit(LeftProbe_1R);
LM_Left_2L = LinearFit(LeftProbe_2L);
LM_Left_2R = LinearFit(LeftProbe_2R);
% LM_Left_1M = LinearFit(LeftProbe_1M);
% LM_Left_2M = LinearFit(LeftProbe_2M);


LM_Right_1L = LinearFit(RightProbe_1L);
LM_Right_1R = LinearFit(RightProbe_1R);
% LM_Right_1M = LinearFit(RightProbe_1M);
% legend('Dentate Nucleus', 'Purkinje Cell Layer', 'Brain Surface', ' ','Probe Traces', 'Fitted Probe Tracks', 'Location' , 'northeast');
% title('3D Plot of Traced Features with Fitted Probe Tracks'); 

%% For probes with DiI


[BS_crossPoint_left1L, BS_distance_left1L] = addPointsToProbe('200116_000', LM_Left_1L, LBSCoords, 70, 1);
[BS_crossPoint_left1R, BS_distance_left1R] = addPointsToProbe('200116_000', LM_Left_1R, LBSCoords, 70, 3);
% [BS_crossPoint_left1M, BS_distance_left1M] = addPointsToProbe('200201_000', LM_Left_1M, LBSCoords, 50, 2);

[BS_crossPoint_left2L, BS_distance_left2L] = addPointsToProbe('200117_000', LM_Left_2L, LBSCoords, 310, 1);
[BS_crossPoint_left2R, BS_distance_left2R] = addPointsToProbe('200117_000', LM_Left_2R, LBSCoords, 310, 3);
% [BS_crossPoint_left2M, BS_distance_left2M] = addPointsToProbe('200202_000', LM_Left_2M, LBSCoords, 50, 2);

[BS_crossPoint_right1L, BS_distance_right1L] = addPointsToProbe('200118_001', LM_Right_1L, RBSCoords, 150, 1);
[BS_crossPoint_right1R, BS_distance_right1R] = addPointsToProbe('200118_001', LM_Right_1R, RBSCoords, 150, 3);
% [BS_crossPoint_right1M, BS_distance_right1M] = addPointsToProbe('200131_000', LM_Right_1M, RBSCoords, -200, 2);

%% for no dii tracks, building linear models for locating the entry points

histoCoords = [BS_crossPoint_left1R; BS_crossPoint_left2R; BS_crossPoint_right1L];
histoCoords = [histoCoords(:, 1) histoCoords(:, 2)];
% BS_crossPoint_left1L; BS_crossPoint_left1M; BS_crossPoint_left1R; 


mouseID = 'cer11';

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')

diiMLCoords = ephysInfo.locationML( (strcmp(mouseID, ephysInfo.mouse)) & (ephysInfo.DiI == 1) );
diiAPCoords = ephysInfo.locationAP( (strcmp(mouseID, ephysInfo.mouse)) & (ephysInfo.DiI == 1) );

nodiiMLCoords = ephysInfo.locationML( (strcmp(mouseID, ephysInfo.mouse)) & (ephysInfo.DiI == 0) );
nodiiAPCoords = ephysInfo.locationAP( (strcmp(mouseID, ephysInfo.mouse)) & (ephysInfo.DiI == 0) );

% enter these info manually!!
manipuCoords = [diiMLCoords, diiAPCoords];
manipuCoords_noDiI = [nodiiMLCoords, nodiiAPCoords];


% plot to verify the histo info
figure;
plot(manipuCoords(:, 1), manipuCoords(:, 2), '.r', 'MarkerSize', 10);
hold on
plot(manipuCoords_noDiI(:, 1), manipuCoords_noDiI(:, 2), '.m', 'MarkerSize', 10);
plot(histoCoords(:, 1), histoCoords(:, 2), '.b', 'MarkerSize', 20);

mdlML_test = fitlm([manipuCoords(:, 1) manipuCoords(:, 2)], histoCoords(:, 1))
mdlAP_test = fitlm([manipuCoords(:, 1) manipuCoords(:, 2)], histoCoords(:, 2))

% fit the LM
fittedML_nodii = predict(mdlML_test, [manipuCoords_noDiI(:, 1), manipuCoords_noDiI(:, 2)])
fittedAP_nodii = predict(mdlAP_test, [manipuCoords_noDiI(:, 1), manipuCoords_noDiI(:, 2)])

% plot the fitted entry points for ni dii tracks
hold on
plot(fittedML_nodii, fittedAP_nodii, '.', 'MarkerSize', 20, 'Color', [1 0.74 0.35])


%% for ni dii tracks, reconstructing no-dii probe traces in 3d view

dirV = nan(6, 3);
dirV(1, :) = LM_Left_1L.dirVect;
dirV(2, :) = LM_Left_1R.dirVect;
dirV(3, :) = LM_Left_2L.dirVect;
dirV(4, :) = -LM_Left_2R.dirVect;
dirV(5, :) = LM_Right_1L.dirVect;
dirV(6, :) = LM_Right_1R.dirVect;
meandirV = mean(dirV);

nodiiLM.BS_crossPoint = [fittedML_nodii, fittedAP_nodii];
nodiiLM.dirV = meandirV;

nodiiLM.BS_crossPoint = [fittedML_nodii(1), fittedAP_nodii(1)];
[BS_crossPoint_right2L, ~] = addPointsToProbe('200113_000', nodiiLM, RBSCoords, 0, 1, {'showPCPoints', true, 'showGCPoints', true, 'noDiIMode', true});

nodiiLM.BS_crossPoint = [fittedML_nodii(2), fittedAP_nodii(2)];
[BS_crossPoint_right2R, ~] = addPointsToProbe('200113_000', nodiiLM, RBSCoords, 0, 2, {'showPCPoints', true, 'showGCPoints', true, 'noDiIMode', true});

nodiiLM.BS_crossPoint = [fittedML_nodii(3), fittedAP_nodii(3)];
[BS_crossPoint_right3L, ~] = addPointsToProbe('200114_000', nodiiLM, RBSCoords, 0, 1, {'showPCPoints', true, 'showGCPoints', true, 'noDiIMode', true});

nodiiLM.BS_crossPoint = [fittedML_nodii(4), fittedAP_nodii(4)];
[BS_crossPoint_right3R, ~] = addPointsToProbe('200114_000', nodiiLM, RBSCoords, 0, 2, {'showPCPoints', true, 'showGCPoints', true, 'noDiIMode', true});





























