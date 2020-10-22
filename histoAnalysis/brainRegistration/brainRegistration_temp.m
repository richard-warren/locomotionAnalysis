% inits
mouse = 'cer12';  % cer12, cmu3, cer18 
scaling = .2;

% prepare histo labels (only need to do once per mouse)
prepareHistoLabels(mouse);


%% load neuron locations
load(fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable.mat'))
cellLocations = ephysHistoTable.GC_Points(strcmp(ephysHistoTable.mouseID, mouse));
cellLocations = cat(1, cellLocations{:}) / 1000;  % ml ap dv

% load mouse brain
data = load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']));

%% plot neurons on brain
figure('color', 'black'); hold on
imagesc(data.ml, data.dv, squeeze(data.imgs(15,:,:))); colormap gray
set(gca, 'ydir', 'reverse')
scatter(cellLocations(:,1), cellLocations(:,3), 20, 'yellow', 'filled')
axis equal
set(gca, 'visible', 'off')

%% load allen brain ccf

section = 470;

% display one coronal section
ccf = loadCCF();
figure('position', [47.00 435.00 1782.00 383.00]);
subplot(1,3,1); imagesc(squeeze(ccf.imgs(section,:,:))); colormap(gca, gray); axis equal 
subplot(1,3,2); imagesc(squeeze(ccf.labelsAll(section,:,:))); colormap(gca, lines); axis equal
subplot(1,3,3); imagesc(squeeze(ccf.labels(section,:,:)>0)); colormap(gca, lines); axis equal


%% test cropped warp
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 200;

labelsCcfCropped = ccf.labels(any(ccf.labels, [2 3]), any(ccf.labels, [1 3]), any(ccf.labels, [1 2]));
labelsCropped = data.labels(any(data.labels, [2 3]), any(data.labels, [1 3]), any(data.labels, [1 2]));
labelsCropped = imresize3(labelsCropped, size(labelsCcfCropped), 'nearest');

tform = imregtform(labelsCropped, labelsCcfCropped, 'affine', optimizer, metric, 'DisplayOptimization', true);
warped = imwarp(labelsCropped, tform, 'OutputView', imref3d(size(labelsCcfCropped)), 'interp', 'nearest');

close all; figure('color', 'white', 'position', [79.00 48.00 1794.00 928.00]); hold on
ax1 = subplot(1,2,1); hold on
plotLabels3D(labelsCcfCropped, 'surfArgs', {'FaceAlpha', .1});
plotLabels3D(labelsCropped);
ax2 = subplot(1,2,2); hold on
plotLabels3D(labelsCcfCropped, 'surfArgs', {'FaceAlpha', .1});
plotLabels3D(warped);
link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', link);

%% show nuclei traces



% mouse = 'cer12';  % cer12, cmu3, cer18 
% data = load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']));


close all; figure('color', 'white'); plotLabels2D(ccf.coarseLabels, 'dim', 'dv', 'smoothing', 3); set(gca, 'visible', 'off')


%% make ccf diagram

close all; figure('color', 'white', 'position', [200 166.00 867.00 851.00]);
ccf = loadCCF();
plotLabels3D(ccf.coarseLabels, 'method', 'contours', 'slices', {'ap'}, 'lineWidth', 2, 'contourAlpha', .1, 'colors', zeros(3,3));
plotLabels3D(ccf.labels, 'method', 'boundary', 'surfArgs', {'edgecolor', 'none', 'facealpha', 1}, 'colors', repelem(lines(3),2,1));
set(gca, 'visible', 'off');
view([0 0])










