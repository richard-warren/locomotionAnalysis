function prepareHistoLabels(mouse, filename, scaling)

% inits
if ~exist('scaling', 'var'); scaling = .2; end
mmPerPixel = .002;  % this should be the same for all images
spreadsheet = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'Sheet', 'histology');
sectionThickness = spreadsheet.sectionThickness_um_(strcmp(spreadsheet.mouse, mouse)) / 1000;  % mm
loadArgs = {'scaling', scaling, 'dataType', 'logical'};
folder = fullfile(getenv('OBSDATADIR'), 'histology', mouse, 'TiffStack');

if isempty(sectionThickness)
    disp('WARNING! Could not determine sectionThickness from sessionInfo histology spreadsheet! Assuming 50 microns...')
    sectionThickness = .05;
end

% import neurotrace brain images
imgs = loadTiffStack(fullfile(getenv('OBSDATADIR'), 'histology', mouse, 'Neurotrace_final.tif'), 'scaling', scaling);

% import labels from .tiff files
nuclei = {'DentateLeft', 'DentateRight', 'InterpositusLeft', 'InterpositusRight', 'FastigialLeft', 'FastigialRight'};
for i = 1:length(nuclei)
    label = loadTiffStack(fullfile(folder, [nuclei{i} '.tif']), loadArgs{:});
    if i==1
        labels = uint8(label);
    else
        labels(label) = i;
    end
end

% load brain surface, purkinje cell layers, probe traces
[surface, pcLayers] = deal(false(size(labels)));
probes = uint8(zeros(size(labels)));

try; surface = surface | loadTiffStack(fullfile(folder, 'BrainSurfaceLeft.tif'), loadArgs{:}); catch; end
try; surface = surface | loadTiffStack(fullfile(folder, 'BrainSurfaceRight.tif'), loadArgs{:}); catch; end

try; pcLayers = pcLayers | loadTiffStack(fullfile(folder, 'PCLayerLeft.tif'), loadArgs{:}); catch; end
try; pcLayers = pcLayers | loadTiffStack(fullfile(folder, 'PCLayerRight.tif'), loadArgs{:}); catch; end

probeFiles = dir(fullfile(folder, 'Probe*'));
for i = 1:length(probeFiles)
    probes(loadTiffStack(fullfile(folder, probeFiles(i).name), loadArgs{:})) = i;
end

% define x, y, z grid
ap = (sectionThickness : sectionThickness : sectionThickness*size(labels,1));
ml = (mmPerPixel : mmPerPixel : mmPerPixel*size(labels,3)) / scaling;
dv = (mmPerPixel : mmPerPixel : mmPerPixel*size(labels,2)) / scaling;

save(filename, 'imgs', 'labels', 'surface', 'pcLayers', 'probes', ...
    'scaling', 'mouse', 'nuclei', 'ap', 'ml', 'dv')

% plot results
figure; montage(permute(imgs, [2 3 1])); colormap(gray)
figure; montage(permute(labels, [2 3 1])); colormap(lines)
