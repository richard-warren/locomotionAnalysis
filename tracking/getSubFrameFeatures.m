function featureVector = getSubFrameFeatures(subFrame, location, includeLocation)
    
    % given a subFrame and location of subFrame, returns feature vector (column)
    % this could in principle be an arbitrary transformation, but for now just flattens subFrame and tacks location to the end
    % location should to be xy
    % only includes locations if includeLocations is true;
    
    featureVector = subFrame(:);
    
    if includeLocation
        featureVector = cat(1, featureVector, location(:));
    end
    
    
end