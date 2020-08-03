%% import tiff stack mask files
mouseID = 'cer18';
folder = fullfile(getenv('OBSDATADIR'), 'histology', mouseID, 'TiffStack');
files = dir(fullfile(folder, '*.tif'));

downsampleRate = 0.3;
thickness = 50;

includeDentate = true;
includeInt = true;
includeFastigial = true;
includePCL_and_BS = true;

saveFiles = true;

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
    probe(i).name = probeFiles(i).name(1:end-4);
    probe(i).coords = importTiffStack(fullfile(folder, probeFiles(i).name), downsampleRate, thickness);   
    probe(i).isDiI = 1;
end


resultsFolder = fullfile(getenv('OBSDATADIR'), 'histology', mouseID, 'Reconstruction');
if saveFiles
    save(fullfile(resultsFolder, 'brain.mat' ));
end
disp('All Done!');



%% Plot brain regions + probe traces
mouseID = 'cer18';
resultsFolder = fullfile('Z:\obstacleData\histology', mouseID, 'Reconstruction');
load(fullfile(resultsFolder, 'brain.mat'))

figure('Color', 'white', 'position', get(0,'ScreenSize'));

showDentate = true;
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


%% fit a line for the probe traces

diI_inds = find([probe.isDiI]); 
for i = 1:length(diI_inds)

    probe(diI_inds(i)).LM = LinearFit(probe(diI_inds(i)).coords);
    
    if contains(probe(diI_inds(i)).name, 'ProbeLeft')
        probe(diI_inds(i)).side = 'left';
    else
        probe(diI_inds(i)).side = 'right';    
    end
    
    switch probe(diI_inds(i)).name(end)
        case 'L'
            probe(diI_inds(i)).shankNum = 1;
        case 'R'
            probe(diI_inds(i)).shankNum = 2;
        case 'M'
            probe(diI_inds(i)).shankNum = 3;
    end
    
end



%% AddPoints: For probes with DiI, calculate brain surface cross points and PC layer cross points
%  For probes with DiI, add good channel points 

mouseID = 'cer18';

disp('getting sessions for this mouse...');
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')

mouse.sessions = {ephysInfo.session{strcmp(mouseID, ephysInfo.mouse) & (ephysInfo.DiI == 1)}};
mouse.DiINotes = {ephysInfo.DiINotes{strcmp(mouseID, ephysInfo.mouse) & (ephysInfo.DiI == 1)}};

% find session number associated with each probe traces
for i = 1:length(diI_inds)
    probe(diI_inds(i)).session = mouse.sessions{contains(mouse.DiINotes, probe(diI_inds(i)).name(1:end-1))}; 
end

% plot BS_crossPoints, GC_crossPoints and good channels on the 3D map
% temporarily use 
for i = 1:length(diI_inds)
    
    if strcmp(probe(diI_inds(i)).side, 'left')
        [probe(diI_inds(i)).BS_crossPoints, probe(diI_inds(i)).PC_crossPoints, probe(diI_inds(i)).GC_points] = ...
            addPointsToProbe(probe(diI_inds(i)).session, probe(diI_inds(i)).LM, brain.LBSCoords, probe(diI_inds(i)).shankNum);
    elseif strcmp(probe(diI_inds(i)).side, 'right')
        [probe(diI_inds(i)).BS_crossPoints, probe(diI_inds(i)).PC_crossPoints, probe(diI_inds(i)).GC_points] = ...
            addPointsToProbe(probe(diI_inds(i)).session, probe(diI_inds(i)).LM, brain.RBSCoords, probe(diI_inds(i)).shankNum);
    end

end



%% for no dii tracks, building linear models for locating the entry points

% ------------------------- settings -------------------------------------
% HistoCoords DiI (get from probe structure)
temp = struct2table(probe(logical([probe.isDiI])));
histoPoints = table(temp.name, temp.BS_crossPoints(:, 1:2), 'VariableNames', {'HistoProbeNames', 'HistoCoords'}); 

% ManipulatorCoords (get from ephysInfo spreadsheet)
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')

ManipuProbeNames = ephysInfo.DiINotes((strcmp(mouseID, ephysInfo.mouse)));
ManipuCoords = [ephysInfo.locationML( (strcmp(mouseID, ephysInfo.mouse))), ephysInfo.locationAP( (strcmp(mouseID, ephysInfo.mouse)))];
manipuPoints = table(ManipuProbeNames, ManipuCoords, 'VariableNames', {'ManipuProbeNames', 'ManipuCoords'});

% Please review and modify these coordinates manually!


% --------------------------- processing ----------------------------------
% plot to verify the histo info
noDiI_inds = contains(manipuPoints.ManipuProbeNames, 'noDiI');
diI_inds = contains(manipuPoints.ManipuProbeNames, 'Probe');

figure;
plot(manipuPoints.ManipuCoords(diI_inds, 1), manipuPoints.ManipuCoords(diI_inds, 2), '.r', 'MarkerSize', 10);
hold on
plot(manipuPoints.ManipuCoords(noDiI_inds, 1), manipuPoints.ManipuCoords(noDiI_inds, 2), '.m', 'MarkerSize', 10);
plot(histoPoints.HistoCoords(:, 1), histoPoints.HistoCoords(:, 2), '.b', 'MarkerSize', 20);

mdlML_test = fitlm(manipuPoints.ManipuCoords(diI_inds, 1)*1000, histoPoints.HistoCoords(:, 1))
mdlAP_test = fitlm(manipuPoints.ManipuCoords(diI_inds, 2)*1000, histoPoints.HistoCoords(:, 2))

% fit the LM
fittedML_nodii = predict(mdlML_test, manipuPoints.ManipuCoords(noDiI_inds, 1)*1000);
fittedAP_nodii = predict(mdlAP_test, manipuPoints.ManipuCoords(noDiI_inds, 2)*1000);

% plot the fitted entry points for ni dii tracks
hold on
plot(fittedML_nodii, fittedAP_nodii, '.', 'MarkerSize', 20, 'Color', [1 0.74 0.35])
legend('manipulator coordinates (DiI)', 'manipulator coordinates (noDiI)', 'histo points (DiI)', 'predicted no DiI points');

%% for ni dii tracks, reconstructing no-dii probe traces in 3d view

dirV = nan(sum([probe.isDiI]), 3);
diI_inds = find([probe.isDiI]);
for i = 1:length(dirV)
    dirV(i, :) = probe(diI_inds(i)).LM.dirVect;
end
nodiiLM.dirV = mean(dirV);

noDiI_inds = find(contains(manipuPoints.ManipuProbeNames, 'noDiI'));
L = sum([probe.isDiI]);
for i = 1:length(noDiI_inds)
    tempName = char(manipuPoints.ManipuProbeNames(noDiI_inds(i)));
    probe(L+i).name = tempName;
    probe(L+i).session = ephysInfo.session{contains(ephysInfo.DiINotes, tempName(1:end-1)) & strcmp(ephysInfo.mouse, mouseID)};
    probe(L+i).side = ephysInfo.side{contains(ephysInfo.DiINotes, tempName(1:end-1))};
    probe(L+i).isDiI = 0;
    probe(L+i).coords = nan;
    probe(L+i).PC_crossPoints = [];
    probe(L+i).GC_points = [];
    
    tempLM = struct();
    tempLM.BS_crossPoints = [fittedML_nodii(noDiI_inds(i)), fittedAP_nodii(noDiI_inds(i))];
    tempLM.dirV = mean(dirV);
    probe(L+i).LM = tempLM;
    
    
    switch tempName(end)
        case 'L'
            probe(L+i).shankNum = 1;
        case 'R'
            probe(L+i).shankNum = 2;
        case 'M'
            probe(L+i).shankNum = 3;
    end
    
   
end

noDiI_inds = find([probe.isDiI] == 0);
for i = 1:length(noDiI_inds)
    
    if strcmp(probe(noDiI_inds(i)).side, 'right')
        [probe(noDiI_inds(i)).BS_crossPoints, probe(noDiI_inds(i)).PC_crossPoints, probe(noDiI_inds(i)).GC_points] = ...
            addPointsToProbe(probe(noDiI_inds(i)).session, probe(noDiI_inds(i)).LM, brain.RBSCoords, probe(noDiI_inds(i)).shankNum,...
            0, {'noDiIMode', true, 'showGCPoints', false});
    elseif strcmp(probe(noDiI_inds(i)).side, 'left')
         [probe(noDiI_inds(i)).BS_crossPoints, probe(noDiI_inds(i)).PC_crossPoints, probe(noDiI_inds(i)).GC_points] = ...
            addPointsToProbe(probe(noDiI_inds(i)).session, probe(noDiI_inds(i)).LM, brain.LBSCoords, probe(noDiI_inds(i)).shankNum,...
            0, {'noDiIMode', true, 'showGCPoints', false});       
    end
    
end

saveFiles = true;
resultsFolder = fullfile(getenv('OBSDATADIR'), 'histology', mouseID, 'Reconstruction');
if saveFiles
    save(fullfile(resultsFolder, 'probe.mat' ));
end
disp('All Done!');


%% Save data to ephysHistoData

mouseID = 'cer12';

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysChannelsInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysChannelsInfo.xlsx'), 'Sheet', 'EphysChannelsInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')

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























