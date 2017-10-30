function pairwisePotentials = getPairwisePotentials(xy, xyLast, maxVelocity, velocityWeight, occludedCost)

% !!! need to document

% restructure matricies to allow direct subtration
xyLastReshaped = repelem(xyLast, size(xy,1), 1);
xyReshaped = repmat(xy, size(xyLast,1), 1);

% compute matrix of pairwise distances (distances is size xy by xyLast)
distances = sqrt(sum((xyReshaped - xyLastReshaped).^2, 2));
distances = reshape(distances, size(xy,1), size(xyLast,1));

% set normalize range and set invalid transitions to 0
pairwisePotentials = maxVelocity - distances;
pairwisePotentials = pairwisePotentials / maxVelocity;
pairwisePotentials(pairwisePotentials<0) = 0;

% weight distances by velocityWeight
pairwisePotentials = pairwisePotentials * velocityWeight;

% add single occluded state potentials
pairwisePotentials = [pairwisePotentials; ones(1,size(xyLast,1))*occludedCost];
pairwisePotentials = [pairwisePotentials, ones(size(xy,1)+1,1)*occludedCost];

% keyboard

