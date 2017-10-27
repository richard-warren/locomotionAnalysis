function featureFrame = getFeatures(frame)
    
    % extracts features from frame or sub-frame
    featureFrame = frame;
    
%     featureFrame = medfilt2(frame, [1 5], 'symmetric');
    featureFrame = imadjust(uint8(featureFrame), [.05 .4], [0 1]);
    
end