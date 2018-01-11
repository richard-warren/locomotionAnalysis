function prepareTrainingData(trainingData, targetSize)

% !!! need to document


% load features
featureDir = [getenv('OBSDATADIR') 'svm\trainingData\' trainingData '\'];
load([featureDir 'labeledFeatures.mat'],...
    'features', 'labels', 'subFrameSize')

egNum = size(features, 2);
categories = unique(labels);


% make category folders
for i = 1:length(categories)
    newDir = [featureDir num2str(categories(i))];
    if ~exist(newDir, 'dir'); mkdir(newDir); end
end


% make training images
for i = 1:egNum
    
    im = reshape(features(:,i), subFrameSize(1), subFrameSize(2));
    im = uint8(imresize(im, 'outputsize', targetSize));
    imwrite(im, [featureDir num2str(labels(i)) '\img' num2str(i) '.png'])
    
end





% featuresReshaped = nan(targetSize(1), targetSize(2), 3, egs);
% 
% % !!! should try implementing this without a for loop, but with multidimensional interpolation...
% for i = 1:egs
%     
%     im = reshape(features(:,i), subFrameSize(1), subFrameSize(2));
%     im = imresize(im, targetSize);
%     featuresReshaped(:,:,:,i) = repmat(im,1,1,3);
%     
% end

