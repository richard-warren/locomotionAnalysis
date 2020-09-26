%% inits
mouse = 'cer18';  % cer18
scaling = .2;

% prepare histo labels (only need to do once per mouse)
% prepareHistoLabels(mouse);


%% load neuron locations
load(fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable.mat'))
cellLocations = ephysHistoTable.GC_Points(strcmp(ephysHistoTable.mouseID, mouse));
cellLocations = cat(1, cellLocations{:}) / 1000;  % convert to mm

% load mouse brain
data = load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']));

% plot neurons on brain
figure; hold on
imagesc(data.ml, data.dv, squeeze(data.imgs(15,:,:))); colormap gray
set(gca, 'ydir', 'reverse')
scatter(cellLocations(:,1), cellLocations(:,3), 20, 'yellow', 'filled')
axis equal

%% load allen brain ccf

section = 470;

% display one coronal section
ccf = loadCCF();
figure('position', [47.00 435.00 1782.00 383.00]);
subplot(1,3,1); imagesc(squeeze(ccf.imgs(section,:,:))); colormap(gca, gray); axis equal 
subplot(1,3,2); imagesc(squeeze(ccf.labelsAll(section,:,:))); colormap(gca, lines); axis equal
subplot(1,3,3); imagesc(squeeze(ccf.labels(section,:,:)>0)); colormap(gca, lines); axis equal

%% test 3d imregister

% settings
mouse = 'cer18';

% load mouse brain and ccf
ccf = loadCCF();
data = load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']));

%% transformations: histo 1-> histoCropped 2-> histoResized 3-> ccfCropped 4-> ccf

% 1) histo to cropped histo coords
apOffset = find(any(data.labels, [2 3]), 1, 'first');
dvOffset = find(any(data.labels, [1 3]), 1, 'first');
mlOffset = find(any(data.labels, [1 2]), 1, 'first');
T1 = eye(4);
T1(end,1:3) = -[dvOffset apOffset mlOffset];  % dv ap ml

% 2) cropped histo to resized coords
apScale  = sum(any(ccf.labels, [2 3])) / sum(any(data.labels, [2 3]));
dvScale  = sum(any(ccf.labels, [1 3])) / sum(any(data.labels, [1 3]));
mlScale  = sum(any(ccf.labels, [1 2])) / sum(any(data.labels, [1 2]));
T2 = diag([dvScale apScale mlScale 1]);  % dv ap ml

% 3) cropped histo to cropped ccf
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 200;
labelsCcfCropped = ccf.labels(any(ccf.labels, [2 3]), any(ccf.labels, [1 3]), any(ccf.labels, [1 2]));
labelsCropped = data.labels(any(data.labels, [2 3]), any(data.labels, [1 3]), any(data.labels, [1 2]));
labelsCropped = imresize3(labelsCropped, size(labelsCcfCropped), 'nearest');
tform = imregtform(labelsCropped, labelsCcfCropped, 'affine', optimizer, metric, 'DisplayOptimization', true);
T3 = tform.T;

% 4) cropped ccf to normal ccf
apOffset = find(any(ccf.labels, [2 3]), 1, 'first');
dvOffset = find(any(ccf.labels, [1 3]), 1, 'first');
mlOffset = find(any(ccf.labels, [1 2]), 1, 'first');
T4 = eye(4);
T4(end,1:3) = [dvOffset apOffset mlOffset];  % dv ap ml

%% plot
close all; figure('color', 'white', 'position', [79.00 48.00 1794.00 928.00]); hold on

% full transform
T = T1 * T2 * T3 * T4;
tform = affine3d(T);
warped = imwarp(data.labels, tform, 'OutputView', imref3d(size(ccf.labels)), 'interp', 'nearest');
% ax1 = subplot(1,2,1); title('full transform'); hold on
plotLabels3D(ccf.labels, 'surfArgs', {'FaceAlpha', .1});
plotLabels3D(warped);
plotLabels3D(ccf.coarseLabels, 'downSampling', 8, 'colors', repmat([0 0 0],3,1), 'surfArgs', {'FaceAlpha', .02, 'EdgeAlpha', .2});

% % cropped only
% tform = affine3d(T3);
% warped = imwarp(labelsCropped, tform, 'OutputView', imref3d(size(labelsCcfCropped)), 'interp', 'nearest');
% ax2 = subplot(1,2,2); title('cropped only'); hold on
% plotLabels3D(labelsCcfCropped, 'surfArgs', {'FaceAlpha', .1});
% plotLabels3D(warped);

% link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
% setappdata(gcf, 'StoreTheLink', link);


%% test cropped warp only
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





