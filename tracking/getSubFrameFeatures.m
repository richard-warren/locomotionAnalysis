function [featureVector, imageVector] = getSubFrameFeatures(subFrame, location, featureSetting, neuralNetwork)
    
    % given a subFrame and location of subFrame, returns feature vector (column)
    % this could in principle be an arbitrary transformation, but for now just flattens subFrame and tacks location to the end
    % location should to be xy
    % only includes locations if includeLocations is true;
    
    imageVector = subFrame(:);
    
    switch featureSetting
        
        case 'imageOnly'
            featureVector = imageVector;
        
        case 'imageWithLocation'
            featureVector = cat(1, imageVector, location(:));
            
        case 'netFeaturesWithLocation'
            % !!!
            
    end
end