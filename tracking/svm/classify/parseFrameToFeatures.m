function [frameFeatures, xVals, yVals] = parseFrameToFeatures(frame, xRes, yRes, subHeight, subWidth, princComps)
    % given a single frame, parseFrame subsamples the frame into subframes
    % of size subHeightXsubWidth, and returns a matrix where is row is a
    % sub-frame WITH FEATURES EXTRACTED VIA getFeatures, and each column is
    % a feature of that sub-frame
    
    frameHeight = size(frame, 1);
    frameWidth = size(frame, 2);
    
    xVals = 1:xRes:frameWidth-subWidth;
    yVals = 1:yRes:frameHeight-subHeight;
    [~, temp] = getFeatures(ones(subHeight, subWidth, size(frame,3)));
    frameFeatures = nan(length(xVals)*length(yVals), length(temp));
    
    rowInd = 1;
    for i = yVals
        for j = xVals
            [~, frameFeatures(rowInd,:)] = getFeatures(frame(i:i+subHeight-1, j:j+subWidth-1, :));
            rowInd = rowInd+1;
        end
    end
end