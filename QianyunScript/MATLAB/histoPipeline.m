
folder = uigetdir();
files = dir(fullfile(folder, '*.tif'));

downsampleRate = 0.3;
thickness = 50;

includeDentate = false;
includeFastigial = false;


% Right Side 
disp('Reformatting right side...');
RICoords = importTiffStack(fullfile(folder, 'InterpositusRight.tif'), downsampleRate, thickness);
RPCLCoords = importTiffStack(fullfile(folder, 'PCLayerRight.tif'), downsampleRate, thickness);
RBSCoords = importTiffStack(fullfile(folder, 'BrainSurfaceRight.tif'), downsampleRate, thickness);

% Left Side
disp('Reformatting left side...');
LICoords = importTiffStack(fullfile(folder, 'InterpositusLeft.tif'), downsampleRate, thickness);
LPCLCoords = importTiffStack(fullfile(folder, 'PCLayerLeft.tif'), downsampleRate, thickness);
LBSCoords = importTiffStack(fullfile(folder, 'BrainSurfaceLeft.tif'), downsampleRate, thickness);

if includeDentate
    LDCoords = importTiffStack(fullfile(folder, 'DentateLeft.tif'), downsampleRate, thickness);
    RDCoords = importTiffStack(fullfile(folder, 'DentateRight.tif'), downsampleRate, thickness);
end

if includeFastigial
    LFCoords = importTiffStack(fullfile(folder, 'FastigialLeft.tif'), downsampleRate, thickness);
    RFCoords = importTiffStack(fullfile(folder, 'FastigialRight.tif'), downsampleRate, thickness);
end



% probes
disp('Reformatting probe traces...');
LeftProbe_1L = importTiffStack(fullfile(folder, 'ProbeLeft_1L.tif'), downsampleRate, thickness);
LeftProbe_1M = importTiffStack(fullfile(folder, 'ProbeLeft_1M.tif'), downsampleRate, thickness);
LeftProbe_1R = importTiffStack(fullfile(folder, 'ProbeLeft_1R.tif'), downsampleRate, thickness);
LeftProbe_2L = importTiffStack(fullfile(folder, 'ProbeLeft_2L.tif'), downsampleRate, thickness);
LeftProbe_2M = importTiffStack(fullfile(folder, 'ProbeLeft_2M.tif'), downsampleRate, thickness);
LeftProbe_2R = importTiffStack(fullfile(folder, 'ProbeLeft_2R.tif'), downsampleRate, thickness);

RightProbe_1L = importTiffStack(fullfile(folder, 'ProbeRight_1L.tif'), downsampleRate, thickness);
RightProbe_1M = importTiffStack(fullfile(folder, 'ProbeRight_1M.tif'), downsampleRate, thickness);
RightProbe_1R = importTiffStack(fullfile(folder, 'ProbeRight_1R.tif'), downsampleRate, thickness);

disp('All Done!');



%% Plot brain regions + probe traces

figure('Color', 'white', 'position', get(0,'ScreenSize'));



xrange = [1000, max(LICoords(:, 1))*2];
yrange = [0, max(LICoords(:, 2))+500]; % +60 just for visualizing purposes

% Left side 
plotRegions3D(LICoords, 20, [1, 0.74, 0.35], xrange, yrange, 0.1); 
hold on
plotRegions3D(LPCLCoords, 10, [0.82, 0.56, 0.97], xrange, yrange, 0.1);
plotRegions3D(LBSCoords, 10, [0.8 0.8 0.8], xrange, yrange, 0.1);

plotRegions3D(LeftProbe_1L, 10, [1.0, 0.43, 0.54], xrange, yrange);
plotRegions3D(LeftProbe_1R, 10, [1.0, 0.43, 0.54], xrange, yrange);
plotRegions3D(LeftProbe_2L, 10, [1.0, 0.43, 0.54], xrange, yrange);
plotRegions3D(LeftProbe_2R, 10, [1.0, 0.43, 0.54], xrange, yrange);
plotRegions3D(LeftProbe_2M, 10, [1.0, 0.43, 0.54], xrange, yrange);
plotRegions3D(LeftProbe_1M, 10, [1.0, 0.43, 0.54], xrange, yrange);

% Right side
plotRegions3D(RICoords, 20, [1, 0.74, 0.35], xrange, yrange, 0.1);
hold on
plotRegions3D(RPCLCoords, 10, [0.82, 0.56, 0.97], xrange, yrange, 0.1);
plotRegions3D(RBSCoords, 10, [0.8 0.8 0.8], xrange, yrange, 0.1);

plotRegions3D(RightProbe_1L, 10, [1.0, 0.43, 0.54], xrange, yrange);
plotRegions3D(RightProbe_1R, 10, [1.0, 0.43, 0.54], xrange, yrange);
plotRegions3D(RightProbe_1M, 10, [1.0, 0.43, 0.54], xrange, yrange);


if includeDentate
    plotRegions3D(LDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange, 0.1);
    plotRegions3D(RDCoords, 20, [0.3176, 0.8314, 0.9608], xrange, yrange, 0.1);
end

if includeFastigial
    plotRegions3D(LFCoords, 20, [0.39, 0.83, 0.075], xrange, yrange);
    plotRegions3D(RFCoords, 20, [0.39, 0.83, 0.075], xrange, yrange);    
end

view(3)


title('3D Plot of Traced Features');
% legend('Dentate', 'Interpositus', 'Fastigial');
%% fit a line for the probe 

LM_Left_1L = LinearFit(LeftProbe_1L);
LM_Left_1R = LinearFit(LeftProbe_1R);
LM_Left_2L = LinearFit(LeftProbe_2L);
LM_Left_2R = LinearFit(LeftProbe_2R);
LM_Left_1M = LinearFit(LeftProbe_1M);
LM_Left_2M = LinearFit(LeftProbe_2M);


LM_Right_1L = LinearFit(RightProbe_1L);
LM_Right_1R = LinearFit(RightProbe_1R);
LM_Right_1M = LinearFit(RightProbe_1M);
% legend('Dentate Nucleus', 'Purkinje Cell Layer', 'Brain Surface', ' ','Probe Traces', 'Fitted Probe Tracks', 'Location' , 'northeast');
% title('3D Plot of Traced Features with Fitted Probe Tracks'); 


LM = cell(1, 6);

LM{1} = LM_Left_1L;
LM{2} = LM_Left_1R;
LM{3} = LM_Left_2L;
LM{4} = LM_Left_2R;
LM{5} = LM_Right_1L;
LM{6} = LM_Right_1R;

% 3 shank version
LM = cell(1, 9);
LM{1} = LM_Right_1L;
LM{2} = LM_Right_1M;
LM{3} = LM_Right_1R;
LM{4} = LM_Left_1L;
LM{5} = LM_Left_1M;
LM{6} = LM_Left_1R;
LM{7} = LM_Left_2L;
LM{8} = LM_Left_2M;
LM{9} = LM_Left_2R;



%% For probes with DiI

[BS_crossPoint_left1L, PC_crossPoint_left1L, GC_left1L, ~] = addPointsToProbe('200201_000', LM_Left_1L, LBSCoords, 50, 1);
[BS_crossPoint_left1R, PC_crossPoint_left1R, GC_left1R, ~] = addPointsToProbe('200201_000', LM_Left_1R, LBSCoords, 50, 3);
[BS_crossPoint_left1M, PC_crossPoint_left1M, GC_left1M, ~] = addPointsToProbe('200201_000', LM_Left_1M, LBSCoords, 50, 2);

[BS_crossPoint_left2L, PC_crossPoint_left2L, GC_left2L, ~] = addPointsToProbe('200202_000', LM_Left_2L, LBSCoords, 0, 1);
[BS_crossPoint_left2R, PC_crossPoint_left2R, GC_left2R, ~] = addPointsToProbe('200202_000', LM_Left_2R, LBSCoords, 0, 3);
[BS_crossPoint_left2M, PC_crossPoint_left2M, GC_left2M, ~] = addPointsToProbe('200202_000', LM_Left_2M, LBSCoords, 0, 2);

[BS_crossPoint_right1L, PC_crossPoint_right1L, GC_right1L, ~] = addPointsToProbe('200131_000', LM_Right_1L, RBSCoords, -200, 1);
[BS_crossPoint_right1R, PC_crossPoint_right1R, GC_right1R, ~] = addPointsToProbe('200131_000', LM_Right_1R, RBSCoords, -200, 3);
[BS_crossPoint_right1M, PC_crossPoint_right1M, GC_right1M, ~] = addPointsToProbe('200131_000', LM_Right_1M, RBSCoords, -200, 2);


%% Save data to ephysHistoData

mouseID = 'cer12';

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysChannelsInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysChannelsInfo.xlsx'), 'Sheet', 'EphysChannelsInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')

shankNum = max(ephysChannelsInfo.PCShankNum(strcmp(mouseID, ephysChannelsInfo.mouse)));
sessions = ephysInfo.session((strcmp(mouseID, ephysInfo.mouse)) & (ephysInfo.DiI == 1) & (ephysInfo.include == 1));


colNum = 6;
ephysHistoData = cell(length(sessions)*shankNum, colNum);

for i = 1 : size(sessions, 1)
    sessionNum = sessions{i};
    offset = ephysInfo.offset((strcmp(sessionNum, ephysInfo.session)));
    if isnan(offset); offset = 0; end
    side = ephysInfo.side((strcmp(sessionNum, ephysInfo.session)));
    if strcmp(side, 'right')
        BSCoords = RBSCoords; 
    else
        BSCoords = LBSCoords; 
    end
    
    temp = (shankNum-1):-1:0;
    for j = 1:shankNum
        ind = i*shankNum - temp(j);
        ephysHistoData{ind, 1} = sessionNum;
        ephysHistoData{ind, 2} = mouseID;
        ephysHistoData{ind, 3} = side;
        [ephysHistoData{ind, 4}, ephysHistoData{ind, 5}, ephysHistoData{ind, 6}, ephysHistoData{ind, 7}] = ...
            addPointsToProbe(sessionNum, LM{ind}, BSCoords, offset, j);
    end
      
end

% check to update and save ephysHistoData
fileName = fullfile(getenv('OBSDATADIR'), 'ephys', 'ephysHistoData', 'ephysHistoData.mat');
if ~exist(fileName)
    save(fileName, 'ephysHistoData');
else
    temp = load(fileName);
    ephysHistoDataOld = temp.ephysHistoData;
    sessionsOld = ephysHistoDataOld(:, 1);
    sessions = ephysHistoData(:, 1);
    inds = zeros(length(sessionsOld), length(sessions));
    for i = 1:length(sessions)
        inds(:, i) = strcmp(sessions{i}, sessionsOld);
    end
    
    if sum(sum(inds)) == 0
        ephysHistoData = [ephysHistoDataOld; ephysHistoData];
        save(fileName, 'ephysHistoData');
    else
        updateIndsOld = [];
        updateIndsNew = [];
        [updateIndsOld, updateIndsNew] = find(inds == 1);
        for j = 1:length(updateIndsOld)
            ephysHistoDataOld(updateIndsOld(j), :) = ephysHistoData(updateIndsNew(j), :);
        end
        ephysHistoData(updateIndsNew, :) = [];
        ephysHistoData = [ephysHistoDataOld; ephysHistoData];
        save(fileName, 'ephysHistoData');
    end    
end

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

mdlML_test = fitlm(manipuCoords(:, 1), histoCoords(:, 1))
mdlAP_test = fitlm(manipuCoords(:, 2), histoCoords(:, 2))

% fit the LM
fittedML_nodii = predict(mdlML_test, manipuCoords_noDiI(:, 1));
fittedAP_nodii = predict(mdlAP_test, manipuCoords_noDiI(:, 2));
% manipuCoords_noDiI(:, 2)

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

mouseID = 'cer11';

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')

sessions = ephysInfo.session((strcmp(mouseID, ephysInfo.mouse)) & (ephysInfo.DiI == 0) & (ephysInfo.include == 1));

ephysHistoData = cell(length(sessions), 6);

for i = 1 : size(sessions, 1)
    sessionNum = sessions{i};
    offset = ephysInfo.offset((strcmp(sessionNum, ephysInfo.session)));
    if isnan(offset); offset = 0; end
    side = ephysInfo.side((strcmp(sessionNum, ephysInfo.session)));
    if strcmp(side, 'right')
        BSCoords = RBSCoords; 
    else
        BSCoords = LBSCoords; 
    end
    
    temp = [1, 0];
    for j = 1:2
        ind = i*2-temp(j);
        ephysHistoData{ind, 1} = sessionNum;
        ephysHistoData{ind, 2} = mouseID;
        ephysHistoData{ind, 3} = side;
        [ephysHistoData{ind, 4}, ephysHistoData{ind, 5}, ephysHistoData{ind, 6}] = ...
            addPointsToProbe(sessionNum, nodiiLM, BSCoords, offset, j,  {'showPCPoints', true, 'showGCPoints', true, 'noDiIMode', true});
    end
      
end

% check to update and save ephysHistoData
fileName = fullfile(getenv('OBSDATADIR'), 'ephys', 'ephysHistoData', 'ephysHistoData.mat');
if ~exist(fileName)
    save(fileName, 'ephysHistoData');
else
    temp = load(fileName);
    ephysHistoDataOld = temp.ephysHistoData;
    sessionsOld = ephysHistoDataOld(:, 1);
    sessions = ephysHistoData(:, 1);
    inds = zeros(length(sessionsOld), length(sessions));
    for i = 1:length(sessions)
        inds(:, i) = strcmp(sessions{i}, sessionsOld);
    end
 
    if sum(sum(inds)) == 0
        ephysHistoData = [ephysHistoDataOld; ephysHistoData];
        save(fileName, 'ephysHistoData');
    else
        [updateIndsOld, updateIndsNew] = find(inds == 1);
        for i = 1:length(updateIndsOld)
            ephysHistoDataOld(updateIndsOld(i), :) = ephysHistoData(updateIndsNew(i), :);
            ephysHistoData(updateIndsNew(i), :) = [];
        end
        ephysHistoData = [ephysHistoDataOld; ephysHistoData];
        save(fileName, 'ephysHistoData');
    end    
end







nodiiLM.BS_crossPoint = [fittedML_nodii(1), fittedAP_nodii(1)];
[BS_crossPoint_right2L, ] = addPointsToProbe('200113_000', nodiiLM, RBSCoords, 0, 1, {'showPCPoints', true, 'showGCPoints', true, 'noDiIMode', true});

nodiiLM.BS_crossPoint = [fittedML_nodii(2), fittedAP_nodii(2)];
[BS_crossPoint_right2R, ~] = addPointsToProbe('200113_000', nodiiLM, RBSCoords, 0, 2, {'showPCPoints', true, 'showGCPoints', true, 'noDiIMode', true});

nodiiLM.BS_crossPoint = [fittedML_nodii(3), fittedAP_nodii(3)];
[BS_crossPoint_right3L, ~] = addPointsToProbe('200114_000', nodiiLM, RBSCoords, 300, 1, {'showPCPoints', true, 'showGCPoints', true, 'noDiIMode', true});

nodiiLM.BS_crossPoint = [fittedML_nodii(4), fittedAP_nodii(4)];
[BS_crossPoint_right3R, ~] = addPointsToProbe('200114_000', nodiiLM, RBSCoords, 300, 2, {'showPCPoints', true, 'showGCPoints', true, 'noDiIMode', true});





























