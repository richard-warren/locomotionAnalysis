%% import tiff stack mask files
mouseID = 'cmu3';
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
    save(fullfile(resultsFolder, 'brain.mat' ), 'brain');
    save(fullfile(resultsFolder, 'probe.mat' ), 'probe');
end
disp('All Done!');



%% Plot brain regions + probe traces
mouseID = 'cer18';
resultsFolder = fullfile('Z:\obstacleData\histology', mouseID, 'Reconstruction');
load(fullfile(resultsFolder, 'brain.mat'))
load(fullfile(resultsFolder, 'probe.mat'))

figure('Color', 'white', 'position', get(0,'ScreenSize'));

showDentate = true;
showInt = true;
showFastigial = true;
showPCL_and_BS = true;


xrange = [0, max(brain.LICoords(:, 1))*2.5];
yrange = [0, max(brain.LICoords(:, 2))+500]; % +60 just for visualizing purposes

markerSize = 20;

PCL_color = [0.82, 0.56, 0.97];
BS_color = [0.8 0.8 0.8];
Int_color = [1, 0.74, 0.35];
Fastigial_color = [0.39, 0.83, 0.075];
Dentate_color = [0.3176, 0.8314, 0.9608];
Probe_color = [1.0, 0.43, 0.54];

transparency = 0.1;

if showPCL_and_BS
    plotRegions3D(brain.LBSCoords, markerSize, BS_color, xrange, yrange, transparency); hold on
    plotRegions3D(brain.RBSCoords, markerSize, BS_color, xrange, yrange, transparency);
    plotRegions3D(brain.LPCLCoords, markerSize, PCL_color, xrange, yrange, transparency);
    plotRegions3D(brain.RPCLCoords, markerSize, PCL_color, xrange, yrange, transparency);
end


if showInt
    plotRegions3D(brain.LICoords, markerSize, Int_color, xrange, yrange, transparency);
    plotRegions3D(brain.RICoords, markerSize, Int_color, xrange, yrange, transparency);
end


if showDentate
    plotRegions3D(brain.LDCoords, markerSize, Dentate_color, xrange, yrange, transparency);
    plotRegions3D(brain.RDCoords, markerSize, Dentate_color, xrange, yrange, transparency);
end

if showFastigial
    plotRegions3D(brain.LFCoords, markerSize, Fastigial_color, xrange, yrange, transparency);
    plotRegions3D(brain.RFCoords, markerSize, Fastigial_color, xrange, yrange, transparency);    
end

diI_inds = find([probe.isDiI]);
for i = 1:length(diI_inds)
    
    plotRegions3D(probe(diI_inds(i)).coords, markerSize, Probe_color, xrange, yrange, transparency);
    
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

% Please review and modify these coordinates manually! Then SAVE IT!
saveFiles = true;
resultsFolder = fullfile(getenv('OBSDATADIR'), 'histology', mouseID, 'Reconstruction');
if saveFiles
    save(fullfile(resultsFolder, 'histoPoints.mat' ), 'histoPoints');
    save(fullfile(resultsFolder, 'manipuPoints.mat' ), 'manipuPoints');
end
disp('File saved!');

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
    if dirV(i, 3) < 0; dirV(i, :) = -dirV(i, :); end
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
            {'noDiIMode', true, 'showGCPoints', true});
    elseif strcmp(probe(noDiI_inds(i)).side, 'left')
         [probe(noDiI_inds(i)).BS_crossPoints, probe(noDiI_inds(i)).PC_crossPoints, probe(noDiI_inds(i)).GC_points] = ...
            addPointsToProbe(probe(noDiI_inds(i)).session, probe(noDiI_inds(i)).LM, brain.LBSCoords, probe(noDiI_inds(i)).shankNum,...
            {'noDiIMode', true, 'showGCPoints', true});       
    end
    
end

saveFiles = true;
resultsFolder = fullfile(getenv('OBSDATADIR'), 'histology', mouseID, 'Reconstruction');
if saveFiles
    save(fullfile(resultsFolder, 'probe.mat' ));
end
disp('Probe.mat Saved!');


%% Save & Update data to ephysHistoTable

mouseID = 'cer18';

temp = struct2table(probe);

newEphysHistoTable = table();
newEphysHistoTable.session = temp.session;
newEphysHistoTable.mouseID(:) = {mouseID};
newEphysHistoTable.side = temp.side;
newEphysHistoTable.probeName = temp.name;
newEphysHistoTable.isDiI = temp.isDiI;
for i = 1:height(newEphysHistoTable)
    newEphysHistoTable.probeID(i) = ephysInfo.probeID(strcmp(ephysInfo.session, temp.session(i)));
end
newEphysHistoTable.shankNum = temp.shankNum;
newEphysHistoTable.BS_crossPoints = temp.BS_crossPoints;
newEphysHistoTable.PC_crossPoints = temp.PC_crossPoints;
newEphysHistoTable.GC_Points = temp.GC_points;


% check to update and save ephysHistoData
fileName = fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable.mat');
ephysHistoFolder = fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData');
if ~exist(fileName)
    ephysHistoTable = newEphysHistoTable;
    save(fileName, 'ephysHistoTable');
else
    old = load(fullfile(ephysHistoFolder, 'ephysHistoTable.mat'));
    old = old.ephysHistoTable;
    new = newEphysHistoTable;
    deleteInds = [];
    for i = 1:height(new)
        overlap = find(strcmp(old.session, new.session(i)) & strcmp(old.probeName, old.name(i)));
        if any(overlap) & length(overlap) == 1
            old(overlap) = new(i);
            deleteInds = [deleteInds, i];
        end
    end
    new(deleteInds, :) = [];
    
    if any(height(new))
        ephysHistoTable = [old; new];
    else
        ephysHistoTable = old;
    end  
    save(fileName, 'ephysHistoTable');
end

























