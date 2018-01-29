function prepareTrainingDataForCnn(trainingData, targetSize)

% !!! need to document, but generally takes features matrix with columns as egs and rows as features and converts into image file in dft folders for each category
% this format can be used by matlabd imageDatastore function


% load features
featureDir = [getenv('OBSDATADIR') 'tracking\trainingData\' trainingData '\'];
load([featureDir 'labeledFeatures.mat'],...
    'features', 'labels', 'subFrameSize')

egNum = size(features, 2);
categories = unique(labels);


% make category folders (or overwrite if already exist)
for i = 1:length(categories)
    newDir = [featureDir num2str(categories(i))];
    if ~exist(newDir, 'dir')
        mkdir(newDir)
    end
end

% make training images
for i = 1:egNum
    
    disp(i/egNum)
    
    im = reshape(features(1:prod(subFrameSize),i), subFrameSize(1), subFrameSize(2));
    im = uint8(imresize(im, 'outputsize', targetSize));
    im = repmat(im, 1, 1, 3);
    imwrite(im, [featureDir num2str(labels(i)) '\img' num2str(i) '.tif'])
    
end

