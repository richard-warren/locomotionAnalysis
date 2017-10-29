function pairwisePotentials = getPairwisePotentials(xy, xyLast, maxVelocity, velocityWeight)

% !!! need to document

% restructure matricies to allow direct subtration
xyLastReshaped = repelem(xyLast, size(xy,1), 1);
xyReshaped = repmat(xy, size(xyLast,1), 1);

% compute matrix of pairwise distances (distances is size xy by xyLast)
distances = sqrt(sum((xyReshaped - xyLastReshaped).^2, 2));
distances = reshape(distances, size(xy,1), size(xyLast,1));

% set invalid transitions to zero
distances(distances > maxVelocity) = 0;

% weight distances by velocityWeight
pairwisePotentials = distances * velocityWeight;