function ccf = loadCCF()
% loads brain images and labels for allen brain common coordinate framework
% 'labels' is categorical mask for ONLY the 6 nuclei, whose identities are
% stored in 'nuclei' // 3D dims are (AP, DV, ML) // 'coarseLabels' has
% labels for cerebrum, cerebellum, and brainstem

% ccf documentation and downloads found here: http://help.brain-map.org/display/mousebrain/API
% ccf label identities found here: http://help.brain-map.org/download/attachments/2818171/MouseCCF.pdf

% todo: add downsamplin

% inits
ccf.nuclei = {'DentateLeft', 'DentateRight', 'InterpositusLeft', 'InterpositusRight', 'FastigialLeft', 'FastigialRight'};

% images
ccfSize = [528 320 456];  % don't change
fid = fopen(fullfile('histoAnalysis', 'brainRegistration', 'allenBrainCCF', 'atlasVolume.raw'), 'r', 'l');
ccf.imgs = fread(fid, prod(ccfSize), 'uint8');
fclose(fid);
ccf.imgs = reshape(ccf.imgs, ccfSize);


% labels
fid = fopen(fullfile('histoAnalysis', 'brainRegistration', 'allenBrainCCF', 'annotation.raw'), 'r', 'l' );
labelsAll = fread(fid, prod(ccfSize), 'uint32' );
fclose(fid);
labelsAll = reshape(labelsAll, ccfSize);
ccf.labelsAll = labelsAll;

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
ccf.labels = labels;

% coordinates
res = .025;
ccf.ap = res : res : size(labels,1)*res;
ccf.dv = res : res : size(labels,2)*res;
ccf.ml = res : res : size(labels,3)*res;

% get outlines of cerebellum, cerebrum, and brainstem
tree_json = fullfile('histoAnalysis', 'brainRegistration', 'allenBrainCCF', 'tree_structure', 'tree_structure.json');
% tree_json = 'D:\github\locomotionAnalysis\histoAnalysis\brainRegistration\allenBrainCCF\tree_structure\tree_structure.json';
ccf.coarseLabels = cast(zeros(size(labelsAll)), 'uint8');
ccf.coarseLabels(ismember(labelsAll, getAllDescenants(tree_json, 'Cerebrum'))) = 1;
ccf.coarseLabels(ismember(labelsAll, getAllDescenants(tree_json, 'Cerebellum'))) = 2;
ccf.coarseLabels(ismember(labelsAll, getAllDescenants(tree_json, 'Brain stem'))) = 3;



