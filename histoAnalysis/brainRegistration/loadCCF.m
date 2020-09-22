function [imgs, labels, nuclei, labelsAll] = loadCCF()
% loads brain images and labels for allen brain common coordinate framework
% 'labels' is categorical mask for ONLY the 6 nuclei, whose identities are
% stored in 'nuclei' // 'labelsAll' has numerical code for all brain
% regions in the atlas // 3D dims are (AP, DV, ML)

% inits
nuclei = {'DentateLeft', 'DentateRight', 'InterpositusLeft', 'InterpositusRight', 'FastigialLeft', 'FastigialRight'};

% images
ccfSize = [528 320 456];  % don't change
fid = fopen(fullfile('histoAnalysis', 'brainRegistration', 'allenBrainCCF', 'atlasVolume.raw'), 'r', 'l');
imgs = fread(fid, prod(ccfSize), 'uint8');
fclose(fid);
imgs = reshape(imgs, ccfSize);

% labels
fid = fopen(fullfile('histoAnalysis', 'brainRegistration', 'allenBrainCCF', 'annotation.raw'), 'r', 'l' );
labelsAll = fread(fid, prod(ccfSize), 'uint32' );
fclose(fid);
labelsAll = reshape(labelsAll, ccfSize);

% nuclei labels
midLine = round(size(labelsAll,3)/2);
leftMask = false(size(labelsAll));
leftMask(:, :, 1:midLine) = true;

labels = cast(zeros(size(labelsAll)), 'uint8');
labels(labelsAll==846 & leftMask) = 1;   % dentate left
labels(labelsAll==846 & ~leftMask) = 2;  % dentate right
labels(labelsAll==91 & leftMask) = 3;    % interpositus left
labels(labelsAll==91 & ~leftMask) = 4;   % interpositus right
labels(labelsAll==989 & leftMask) = 5;   % fastigial left
labels(labelsAll==989 & ~leftMask) = 6;  % fastigial right

