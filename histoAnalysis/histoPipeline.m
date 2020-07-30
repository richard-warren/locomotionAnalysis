mouseID = 'cmu3';
folder = fullfile(getenv('OBSDATADIR'), 'histology', mouseID, 'TiffStack');
files = dir(fullfile(folder, '*.tif'));

downsampleRate = 0.3;
thickness = 50;

includeDentate = false;
includeInt = false;
includeFastigial = false;
includePCL_and_BS = true;

saveFiles = false;

% brain regions
if includeInt
    disp('Reformatting Int...');
    brain.RICoords = importTiffStack(fullfile(folder, 'InterpositusRight.tif'), downsampleRate, thickness);
    brain.LICoords = importTiffStack(fullfile(folder, 'InterpositusLeft.tif'), downsampleRate, thickness);
end

if includeDentate
    disp('Reformatting Dentate...');
    brain.LDCoords = importTiffStack(fullfile(folder, 'DentateLeft.tif'), downsampleRate, thickness);
    brain.RDCoords = importTiffStack(fullfile(folder, 'DentateRight.tif'), downsampleRate, thickness);
end

if includeFastigial
    disp('Reformatting Fastigial...');
    brain.LFCoords = importTiffStack(fullfile(folder, 'FastigialLeft.tif'), downsampleRate, thickness);
    brain.RFCoords = importTiffStack(fullfile(folder, 'FastigialRight.tif'), downsampleRate, thickness);
end

if includePCL_and_BS
    disp('Reformatting PC Layer and Brain Surface...');
    brain.RPCLCoords = importTiffStack(fullfile(folder, 'PCLayerRight.tif'), downsampleRate, thickness);
    brain.RBSCoords = importTiffStack(fullfile(folder, 'BrainSurfaceRight.tif'), downsampleRate, thickness);
    brain.LPCLCoords = importTiffStack(fullfile(folder, 'PCLayerLeft.tif'), downsampleRate, thickness);
    brain.LBSCoords = importTiffStack(fullfile(folder, 'BrainSurfaceLeft.tif'), downsampleRate, thickness);
end


% probes
disp('Reformatting probe traces...');

probeFiles = files(contains({files.name}, 'Probe'));
probe = struct();

for i = 1:length(probeFiles)
    probe(i).name = probeFiles(i).name;
    probe(i).coords = importTiffStack(fullfile(folder, probeFiles(i).name), downsampleRate, thickness);      
end


resultsFolder = fullfile(getenv('OBSDATADIR'), 'histology', mouseID, 'Reconstruction');
if saveFiles
    save(fullfile(resultsFolder, 'brain.mat' ));
end
disp('All Done!');



%% Plot brain regions + probe traces
mouseID = 'cmu3';
resultsFolder = fullfile('Z:\obstacleData\histology', mouseID, 'Reconstruction');
load(fullfile(resultsFolder, 'brain.mat'))

figure('Color', 'white', 'position', get(0,'ScreenSize'));

showDentate = false;
showInt = true;
showFastigial = true;
showPCL_and_BS = true;


xrange = [0, max(brain.LICoords(:, 1))*2.5];
yrange = [0, max(brain.LICoords(:, 2))+500]; % +60 just for visualizing purposes

PCL_color = [0.82, 0.56, 0.97];
BS_color = [0.8 0.8 0.8];
Int_color = [1, 0.74, 0.35];
Fastigial_color = [0.39, 0.83, 0.075];
Dentate_color = [0.3176, 0.8314, 0.9608];
Probe_color = [1.0, 0.43, 0.54];

transparency = 0.1;

if showPCL_and_BS
    plotRegions3D(brain.LBSCoords, 20, BS_color, xrange, yrange, transparency); hold on
    plotRegions3D(brain.RBSCoords, 20, BS_color, xrange, yrange, transparency);
    plotRegions3D(brain.LPCLCoords, 20, PCL_color, xrange, yrange, transparency);
    plotRegions3D(brain.RPCLCoords, 20, PCL_color, xrange, yrange, transparency);
end


if showInt
    plotRegions3D(brain.LICoords, 20, Int_color, xrange, yrange, transparency);
    plotRegions3D(brain.RICoords, 20, Int_color, xrange, yrange, transparency);
end


if showDentate
    plotRegions3D(brain.LDCoords, 20, Dentate_color, xrange, yrange, transparency);
    plotRegions3D(brain.RDCoords, 20, Dentate_color, xrange, yrange, transparency);
end

if showFastigial
    plotRegions3D(brain.LFCoords, 20, Fastigial_color, xrange, yrange, transparency);
    plotRegions3D(brain.RFCoords, 20, Fastigial_color, xrange, yrange, transparency);    
end


for i = 1:length(probe)
    
    plotRegions3D(probe(i).coords, 20, Probe_color, xrange, yrange, transparency);
    
end




view(3)


title('3D Plot of Traced Features');


%% fit a line for the probe 

for i = 1:length(probe)

    probe(i).LM = LinearFit(probe(i).coords);
    
    if contains(probe(i).name, 'ProbeLeft')
        probe(i).side = 'left';
    else
        probe(i).side = 'right';    
    end
    
    switch probe(i).name(end-4)
        case 'L'
            probe(i).shankNum = 1;
        case 'R'
            probe(i).shankNum = 2;
        case 'M'
            probe(i).shankNum = 3;
    end
    
end



%% For probes with DiI

mouseID = 'cmu3';

disp('getting sessions for this mouse...');
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')

mouse.sessions = {ephysInfo.session{strcmp(mouseID, ephysInfo.mouse) & (ephysInfo.DiI == 1)}};
mouse.DiINotes = {ephysInfo.DiINotes{strcmp(mouseID, ephysInfo.mouse) & (ephysInfo.DiI == 1)}};

% find session number associated with each probe traces
for i = 1:length(probe)
    probe(i).session = mouse.sessions{contains(mouse.DiINotes, probe(i).name(1:end-5))}; 
end

% plot BS_crossPoints, GC_crossPoints and good channels on the 3D map
% temporarily use 
for i = 1:length(probe)
    
    if strcmp(probe(i).side, 'left')
        [probe(i).BS_crossPoints, probe(i).PC_crossPoints, probe(i).GC_points] = ...
            addPointsToProbe(probe(i).session, probe(i).LM, brain.LBSCoords, probe(i).shankNum);
    elseif strcmp(probe(i).side, 'right')
        [probe(i).BS_crossPoints, probe(i).PC_crossPoints, probe(i).GC_points] = ...
            addPointsToProbe(probe(i).session, probe(i).LM, brain.RBSCoords, probe(i).shankNum);
    end

end


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

histoCoords = [BS_crossPoint_right1; BS_crossPoint_left1; BS_crossPoint_left2];
histoCoords = [histoCoords(:, 1) histoCoords(:, 2)];
% BS_crossPoint_left1L; BS_crossPoint_left1M; BS_crossPoint_left1R; 


mouseID = 'cer10';

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

dirV = nan(3, 3);

dirV(1, :) = LM_Right_1.dirVect;
dirV(2, :) = LM_Left_1.dirVect;
dirV(3, :) = LM_Left_2.dirVect;

% dirV(1, :) = LM_Left_1L.dirVect;
% dirV(2, :) = LM_Left_1R.dirVect;
% dirV(3, :) = LM_Left_2L.dirVect;
% dirV(4, :) = -LM_Left_2R.dirVect;
% dirV(5, :) = LM_Right_1L.dirVect;
% dirV(6, :) = LM_Right_1R.dirVect;

meandirV = mean(dirV);

nodiiLM.BS_crossPoint = [fittedML_nodii, fittedAP_nodii];
nodiiLM.dirV = LM_Right_1.dirVect;

mouseID = 'cer10';

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
[BS_crossPoint_right2, PC_crossPoint_right2, GC_right2, ~] = addPointsToProbe('200311_000', nodiiLM, RBSCoords, 0, 1, {'showPCPoints', true, 'showGCPoints', true, 'noDiIMode', true});

% nodiiLM.BS_crossPoint = [fittedML_nodii(2), fittedAP_nodii(2)];
% [BS_crossPoint_right2R, ~] = addPointsToProbe('200113_000', nodiiLM, RBSCoords, 0, 2, {'showPCPoints', true, 'showGCPoints', true, 'noDiIMode', true});
% 
% nodiiLM.BS_crossPoint = [fittedML_nodii(3), fittedAP_nodii(3)];
% [BS_crossPoint_right3L, ~] = addPointsToProbe('200114_000', nodiiLM, RBSCoords, 300, 1, {'showPCPoints', true, 'showGCPoints', true, 'noDiIMode', true});
% 
% nodiiLM.BS_crossPoint = [fittedML_nodii(4), fittedAP_nodii(4)];
% [BS_crossPoint_right3R, ~] = addPointsToProbe('200114_000', nodiiLM, RBSCoords, 300, 2, {'showPCPoints', true, 'showGCPoints', true, 'noDiIMode', true});
% 




























