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

% load mouse brain and ccf
ccf = loadCCF();
data = load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']));

% todo: chop off ends of nuclei in ref

%% estimate initial matrix
apScale  = sum(any(ccf.labels, [2 3])) / sum(any(data.labels, [2 3]));
dvScale  = sum(any(ccf.labels, [1 3])) / sum(any(data.labels, [1 3]));
mlScale  = sum(any(ccf.labels, [1 2])) / sum(any(data.labels, [1 2]));

apOffset = find(any(ccf.labels, [2 3]), 1, 'first') - find(any(data.labels, [2 3]), 1, 'first')*apScale;
dvOffset = find(any(ccf.labels, [1 3]), 1, 'first') - find(any(data.labels, [1 3]), 1, 'first')*dvScale;
mlOffset = find(any(ccf.labels, [1 2]), 1, 'first') - find(any(data.labels, [1 2]), 1, 'first')*mlScale;

T = diag([dvScale, apScale, mlScale, 1]);
T(end, 1:3) = [dvOffset apOffset mlOffset];  % (dv, ap, ml)
tformInit = affine3d(T);

% find transformation
[optimizer, metric] = imregconfig('monomodal');
% optimizer.InitialRadius = 0.009;
% optimizer.Epsilon = 1.5e-4;
% optimizer.GrowthFactor = 1.01;
% optimizer.MaximumIterations = 20;
% optimizer.MaximumStepLength = 1e-2;

% create masks for left side of brain
% leftMask = false(size(data.labels)); leftMask(:, :, 1:round(size(data.labels,3)/2)) = true;
% leftMaskCcf = false(size(ccf.labels)); leftMaskCcf(:, :, 1:round(size(ccf.labels,3)/2)) = true;
% labelsAligned = cast(zeros(size(ccf.labels)), 'uint8');
% tforms = cell(1,2);

% preprocess images for alignment
% labels = imwarp(data.labels, tformInit, 'OutputView', imref3d(size(ccf.labels)), 'interp', 'nearest');  % apply initial transformation

% logical
% labels = uint8(labels>0);
% labelsRef = uint8(ccf.labels>0);

% cropped and size matched
labelsRef = ccf.labels(any(ccf.labels, [2 3]), any(ccf.labels, [1 3]), any(ccf.labels, [1 2]));
labels = data.labels(any(data.labels, [2 3]), any(data.labels, [1 3]), any(data.labels, [1 2]));
labels = imresize3(labels, size(labelsRef), 'nearest');

% smoothed double
% smoothing = 10;
% labels = imgaussfilt3(double(labels>0), smoothing);
% labelsRef = imgaussfilt3(double(labelsRef>0), smoothing);

% for isLeft = [true, false]
%     if isLeft; disp('registering left side...'); else; disp('registering right side...'); end
%     mask = leftMask==isLeft;
%     maskRef = leftMaskCcf==isLeft;
%     
% %     tform = imregtform(uint8(data.labels .* uint8(mask)), uint8(ccf.labels .* uint8(maskRef)), ...  % categorical
% %         'affine', optimizer, metric, 'InitialTransformation', tformInit);
% %     tform = imregtform(uint8(data.labels>0 .* mask*255), uint8(ccf.labels>0 .* maskRef), ...  % binary
% %         'affine', optimizer, metric, 'InitialTransformation', tformInit);     
%     tform = tformInit;
%     
%     warped = imwarp(data.labels .* uint8(mask), tform, 'OutputView', imref3d(size(ccf.labels)), 'interp', 'nearest');
%     labelsAligned(maskRef) = warped(maskRef);
%     tforms{isLeft+1} = tform;
% end

% optimize both sides at once
disp('registering both sides at once...')
tform = imregtform(labels, labelsRef, 'affine', optimizer, metric, 'DisplayOptimization', true);     
% warpedInt = imwarp(data.labels, tformInit, 'OutputView', imref3d(size(ccf.labels)), 'interp', 'nearest');
% warped = imwarp(warpedInt, tform, 'OutputView', imref3d(size(ccf.labels)), 'interp', 'nearest');

warped = imwarp(labels, tform, 'OutputView', imref3d(size(ccf.labels)), 'interp', 'nearest');


%% plot
close all; figure('color', 'white', 'position', [90.00 251.00 1763.00 607.00]); hold on

ax1 = subplot(1,2,1); title('original'); hold on
plotLabels3D(labelsRef, 'colors', repelem(lines(3),2,1)*.5, 'surfArgs', {'FaceAlpha', .1});
plotLabels3D(labels, 'colors', repelem(lines(3),2,1), 'surfArgs', {'FaceAlpha', .5});

ax2 = subplot(1,2,2); title('morphed'); hold on
plotLabels3D(labelsRef, 'colors', repelem(lines(3),2,1)*.5, 'surfArgs', {'FaceAlpha', .1});
plotLabels3D(warped, 'colors', repelem(lines(3),2,1), 'surfArgs', {'FaceAlpha', .5});

link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', link);

%% test plot3D function


% ccf = loadCCF();
close all; figure('color', 'white', 'position', [705.00 108.00 1176.00 837.00]); hold on
plotLabels3D(ccf.labels);







