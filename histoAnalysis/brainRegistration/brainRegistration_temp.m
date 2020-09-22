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

regionIds = [989, 91, 846];  % fastigial, interpositus, dentate
ccfSize = [528 320 456];     % don't change

% load ccf average brain
[brainImgsRef, labelsRef, ~, labelsAll] = loadCCF();

% display one coronal section
section = 470;
close all; figure('position', [47.00 435.00 1782.00 383.00]);
subplot(1,3,1); imagesc(squeeze(brainImgsRef(section,:,:))); colormap(gca, gray); axis equal 
subplot(1,3,2); imagesc(squeeze(labelsAll(section,:,:))); colormap(gca, lines); axis equal
subplot(1,3,3); imagesc(squeeze(labelsRef(section,:,:)>0)); colormap(gca, lines); axis equal

%% test 2d imregister

% settings
section = 16;
sectionRef = 450;
halfOnly = true;


% load mouse brain and ccf
[brainImgsRef, labelsRef] = loadCCF();
data = load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']));

label = squeeze(data.labels(section,:,:))>0;
labelRef = squeeze(labelsRef(sectionRef,:,:))>0;
img = squeeze(data.imgs(section,:,:));
imgRef = squeeze(brainImgsRef(sectionRef,:,:));

% resize images to help alignment algorithm...
label = imresize(label, size(labelRef));
img = imresize(img, size(labelRef));

if halfOnly
    inds = 1 : floor(size(label,2)/2);
%     inds = floor(size(label,2)/2) : size(label,2);
    label = label(:,inds);
    labelRef = labelRef(:,inds);
    img = img(:,inds);
    imgRef = imgRef(:,inds);
end


[optimizer, metric] = imregconfig('monomodal');
% optimizer.InitialRadius = 0.1;
% optimizer.Epsilon = 1.5e-3;
% optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 500;
% labelAligned = imregister(uint8(label), uint8(labelRef), 'affine', optimizer, metric);

tform = imregtform(uint8(label), uint8(labelRef), 'affine', optimizer, metric);
labelAligned = imwarp(label, tform, 'OutputView', imref2d(size(label)));
imgAligned = imwarp(img, tform, 'OutputView', imref2d(size(img)));

% plot
close all
figure('position', [110.00 230.00 1670.00 667.00])
subplot(2,4,1); imagesc(label); daspect([1 1 1]);
subplot(2,4,2); imagesc(labelRef); daspect([1 1 1])
subplot(2,4,3); imshowpair(labelRef, label, 'Scaling', 'joint'); daspect([1 1 1])
subplot(2,4,4); imshowpair(labelRef, labelAligned); daspect([1 1 1])

subplot(2,4,5); image(img); daspect([1 1 1]);
subplot(2,4,6); image(imgRef); daspect([1 1 1])
subplot(2,4,7); imshowpair(imgRef, img, 'Scaling', 'joint'); daspect([1 1 1])
subplot(2,4,8); imshowpair(imgRef, imgAligned, 'Scaling', 'joint'); daspect([1 1 1])


%% test 3d imregister

% settings

% load mouse brain and ccf
[brainImgsRef, labelsRef] = loadCCF();
data = load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']));
labels = data.labels;

% slice volumes to keep only sections containing nuclei
labelsRef = labelsRef(any(labelsRef, [2 3]), any(labelsRef, [1 3]), any(labelsRef, [1 2]));
labels = labels(any(labels, [2 3]), any(labels, [1 3]), any(labels, [1 2]));
% todo: incorporate user knowledge of final section in tissue... // store information about this slice, or don't slice at all (if that works) 

% resize images to match
labels = imresize3(labels, size(labelsRef), 'nearest');

% find transformation
[optimizer, metric] = imregconfig('monomodal');
% optimizer.InitialRadius = 0.1;
% optimizer.Epsilon = 1.5e-3;
% optimizer.GrowthFactor = 1.01;
% optimizer.MaximumIterations = 500;

% create mask for left side of brain
midLine = round(size(labels,3)/2);
leftMask = cast(zeros(size(labels)), 'uint8');
leftMask(:, :, 1:midLine) = 1;
labelsAligned = cast(zeros(size(labels)), 'uint8');

for mask = {leftMask, uint8(~leftMask)}
    tform = imregtform(uint8(labels.*mask{1}>0), uint8(labelsRef.*mask{1}>0), 'affine', optimizer, metric);
    warped = imwarp(labels.*mask{1}, tform, 'OutputView', imref3d(size(labels)), 'interp', 'nearest');
    labelsAligned(mask{1}>0) = warped(mask{1}>0);
end


%% plot
close all; figure('color', 'white', 'position', [90.00 251.00 1763.00 607.00]); hold on

ax1 = subplot(1,2,1); title('original'); hold on
plotLabels3D(labelsRef, 'colors', repelem(lines(3),2,1)*.5, 'surfArgs', {'FaceAlpha', .1});
plotLabels3D(labels, 'colors', repelem(lines(3),2,1), 'surfArgs', {'FaceAlpha', .5});

ax2 = subplot(1,2,2); title('morphed'); hold on
plotLabels3D(labelsRef, 'colors', repelem(lines(3),2,1)*.5, 'surfArgs', {'FaceAlpha', .1});
plotLabels3D(labelsAligned, 'colors', repelem(lines(3),2,1), 'surfArgs', {'FaceAlpha', .5});

link = linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
setappdata(gcf, 'StoreTheLink', link);

view(-45, 30)

%% test plot3D function

[~, labelsRef] = loadCCF();
data = load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']));

close all; figure('color', 'white', 'position', [800 228.00 987.00 750.00]); hold on
plotLabels3D(labelsRef, 'colors', repelem(lines(3),2,1));
plotLabels3D(data.labels, 'colors', repelem(lines(3),2,1));



