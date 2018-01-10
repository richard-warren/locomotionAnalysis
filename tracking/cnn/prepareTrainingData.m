function featuresReshaped = prepareTrainingData(features, subFrameSize, targetSize)

% !!! needs documentation, but generally restructures 2D feature matrix for grayscale images, where each col is training eg and each row is feature into 4D tensor for use in convolutional neural net
% dimensions of featuresReshaped are [imageHeight x imageWidth x 3 x numberOfExamples]
% the third dimension is redundant because there is no color dimension in the original gray images, but this dimension is added so data are properly structured for existing CNNs


% featuresReshaped = reshape(features, [subFrameSize(1), subFrameSize(1), size(features, 2)]);
% featuresReshaped = cat(4, featuresReshaped, featuresReshaped, featuresReshaped);
% featuresReshaped = permute(featuresReshaped, [1 2 4 3]);


egs = size(features, 2);
featuresReshaped = nan(targetSize(1), targetSize(2), 3, egs);

% !!! should try implementing this without a for loop, but with multidimensional interpolation...
for i = 1:egs
    
    im = reshape(features(:,i), subFrameSize(1), subFrameSize(2));
    im = imresize(im, targetSize);
    featuresReshaped(:,:,:,i) = repmat(im,1,1,3);
    
end