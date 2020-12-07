function prepareHistoLabels(mouse, filename, scaling)

% inits
if ~exist('scaling', 'var'); scaling=.2; end
mmPerPixel = .002;  % this should be the same for all images
spreadsheet = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'Sheet', 'histology');
sectionThickness = spreadsheet.sectionThickness_um_(strcmp(spreadsheet.mouse, mouse)) / 1000;  % mm

if isempty(sectionThickness)
    disp('WARNING! Could not determine sectionThickness from sessionInfo histology spreadsheet! Assuming 50 microns...')
    sectionThickness = .05;
end

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

% define x, y, z grid (check this is correct)
ap = (sectionThickness : sectionThickness : sectionThickness*size(labels,1));
ml = (mmPerPixel : mmPerPixel : mmPerPixel*size(labels,3)) / scaling;
dv = (mmPerPixel : mmPerPixel : mmPerPixel*size(labels,2)) / scaling;

save(filename, 'imgs', 'labels', 'scaling', 'mouse', 'nuclei', 'ap', 'ml', 'dv')

% plot results
figure; montage(permute(imgs, [2 3 1])); colormap(gray)
figure; montage(permute(labels, [2 3 1])); colormap(lines)
