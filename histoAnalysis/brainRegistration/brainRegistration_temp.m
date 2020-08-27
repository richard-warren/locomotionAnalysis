%% settings
mouse = 'cer5';
scaling = .2;

%% import neuron locations
load(fullfile(getenv('OBSDATADIR'), 'histology', '0_ephysHistoData', 'ephysHistoTable_old.mat'))

cellLocations = ephysHistoTable.GC_Points(strcmp(ephysHistoTable.mouseID, mouse));
cellLocations = cat(1, cellLocations{:});
cellLocations = cellLocations / 2 * scaling;  % convert from microns to pixels

close all;
figure; hold on
image(squeeze(neurotrace(10,:,:))); colormap gray
set(gca, 'ydir', 'reverse')
scatter(cellLocations(:,1), cellLocations(:,3), 20, 'yellow', 'filled')

%% load allen brain ccf

regionIds = [989, 91, 846];  % fastigial, interpositus, dentate
ccfSize = [528 320 456];  % don't change

% load ccf average brain
fid = fopen(fullfile('histoAnalysis', 'brainRegistration', 'allenBrainCCF', 'atlasVolume.raw'), 'r', 'l');
brainImgsRef = fread(fid, prod(ccfSize), 'uint8');
fclose(fid);
brainImgsRef = reshape(brainImgsRef, ccfSize);

% load ccf annotations
fid = fopen(fullfile('histoAnalysis', 'brainRegistration', 'allenBrainCCF', 'annotation.raw'), 'r', 'l' );
labelsRef = fread(fid, prod(ccfSize), 'uint32' );
fclose(fid);
labelsRef = reshape(labelsRef, ccfSize);
 
% display one coronal section
section = 460;
close all; figure('position', [47.00 435.00 1782.00 383.00]);
subplot(1,3,1); imagesc(squeeze(brainImgsRef(section,:,:))); colormap(gca, gray); 
subplot(1,3,2); imagesc(squeeze(labelsRef(section,:,:))); colormap(gca, lines); 
subplot(1,3,3); imagesc(squeeze(ismember(labelsRef(section,:,:), regionIds))); colormap(gca, lines);

%% test imregister

mouse = 'cer5';

% load ccf
ccfSize = [528 320 456];  % don't change
fid = fopen(fullfile('histoAnalysis', 'brainRegistration', 'allenBrainCCF', 'atlasVolume.raw'), 'r', 'l');
brainImgsRef = fread(fid, prod(ccfSize), 'uint8');
fclose(fid);
brainImgsRef = reshape(brainImgsRef, ccfSize);

fid = fopen(fullfile('histoAnalysis', 'brainRegistration', 'allenBrainCCF', 'annotation.raw'), 'r', 'l' );
labelsRef = fread(fid, prod(ccfSize), 'uint32' );
fclose(fid);
labelsRef = reshape(labelsRef, ccfSize);
labelsRef = ismember(labelsRef, [989, 91, 846]);  % fastigial, interpositus, dentate

% load labels
load(fullfile(getenv('SSD'), 'paper2', 'histo', 'histoLabels', [mouse '_histoLabels.mat']), ...
    'brainImgs', 'labels', 'scaling', 'brainRegions')


%%

% settings
section = 16;
sectionRef = 450;
halfOnly = true;

label = squeeze(labels(section,:,:))>0;
labelRef = squeeze(labelsRef(sectionRef,:,:))>0;
img = squeeze(brainImgs(section,:,:));
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












