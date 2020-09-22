function prepareHistoLabels(mouse, scaling)

% inits
if ~exist('scaling', 'var'); scaling=.2; end
mmPerPixel = .002;  % this should be the same for all images
spreadsheet = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'Sheet', 'histology');
sectionThickness = spreadsheet.sectionThickness_um_(strcmp(spreadsheet.mouse, mouse)) / 1000;  % mm

% import neurotrace brain images
imgs = loadTiffStack(fullfile(getenv('OBSDATADIR'), 'histology', mouse, 'Neurotrace_final.tif'), 'scaling', scaling);

% import labels from .tiff files
nuclei = {'DentateLeft', 'DentateRight', 'InterpositusLeft', 'InterpositusRight', 'FastigialLeft', 'FastigialRight'};
for i = 1:length(nuclei)
    fileName = fullfile(getenv('OBSDATADIR'), 'histology', mouse, 'TiffStack', [nuclei{i} '.tif']);
    label = loadTiffStack(fileName, 'scaling', scaling, 'dataType', 'logical');
    if i==1
        labels = uint8(label);
    else
        labels(label) = i;
    end
end

% define x, y, z grid
ap = (0 : sectionThickness : sectionThickness*(size(labels,1)-1));
ml = (0 : mmPerPixel : mmPerPixel*(size(labels,3)-1)) / scaling;
dv = (0 : mmPerPixel : mmPerPixel*(size(labels,2)-1)) / scaling;

save(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']), ...
    'imgs', 'labels', 'scaling', 'mouse', 'nuclei', 'ap', 'ml', 'dv')

% plot results
figure; montage(permute(imgs, [2 3 1])); colormap(gray)
figure; montage(permute(labels, [2 3 1])); colormap(lines)
