function prepareHistoLabels(mouse, scaling)

if ~exist('scaling', 'var'); scaling=.2; end

% import neurotrace brain images
brainImgs = loadTiffStack(fullfile(getenv('OBSDATADIR'), 'histology', mouse, 'Neurotrace_final.tif'), 'scaling', scaling);

% import labels from .tiff files
brainRegions = {'DentateLeft', 'DentateRight', 'InterpositusLeft', 'InterpositusRight', 'FastigialLeft', 'FastigialRight'};
for i = 1:length(brainRegions)
    fileName = fullfile(getenv('OBSDATADIR'), 'histology', mouse, 'TiffStack', [brainRegions{i} '.tif']);
    label = loadTiffStack(fileName, 'scaling', scaling, 'dataType', 'logical');
    if i==1
        labels = uint8(label);
    else
        labels(label) = i;
    end
end

save(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']), ...
    'brainImgs', 'labels', 'scaling', 'mouse', 'brainRegions')

figure; montage(permute(brainImgs, [2 3 1])); colormap(gray)
figure; montage(permute(labels, [2 3 1])); colormap(lines)