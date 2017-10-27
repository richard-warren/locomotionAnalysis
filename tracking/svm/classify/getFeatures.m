function [featureFrame, featureVector] = getFeatures(subFrame)
    % extracts features from subFrame
    % featureFrame is a visualization of the extracted features
    % featureVector are the features themselves
    
%     subFrame = medfilt2(subFrame, [3 3]);
%     featureFrame = imadjust(uint8(subFrame), [.05 .4], [0 1]);

%     subFrame = imadjust(subFrame, [.05 .5], [0 1]);
    
    featureFrame = subFrame;
    featureVector = subFrame(:)';
end