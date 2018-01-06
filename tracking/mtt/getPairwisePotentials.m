function pairwisePotentials = getPairwisePotentials(xy, xyLast, maxVelocity, velocityWeight, occludedCost, gridPts)

% given track locations in two adjacent frames, computes a matrix representing the likelihood of transitioning
% from every location in previous frame (xyLast) to every location in current frame (xy)
% this likelihood is based on the distance between track (ie the velocity necessary to move from one to
% another track


% initializations
numOccluded = size(gridPts,1);
numLast = size(xyLast,1);
numCurrent = size(xy,1);




% COST OF TRANSITIONING BETWEEN TRACKS

% restructure matrices to allow direct subtration
xyLastReshaped = repelem(xyLast, numCurrent, 1);
xyReshaped = repmat(xy, numLast, 1);

% compute matrix of pairwise distances (distances is size xy by xyLast)
distances = sqrt(sum((xyReshaped - xyLastReshaped).^2, 2));
distances = reshape(distances, size(xy,1), size(xyLast,1));

% normalize range and set invalid transitions to 0
transitionCosts = distances;
transitionCosts(transitionCosts > maxVelocity) = maxVelocity;
transitionCosts = (maxVelocity - transitionCosts) / maxVelocity;

% weight distances by velocityWeight
transitionCosts = transitionCosts * velocityWeight;




% COST OF BECOMING OCCLUDED

nearestOccludedInds = knnsearch(gridPts, xyLast);
occlusionCosts = zeros(numOccluded, numLast);
inds = sub2ind(size(occlusionCosts), nearestOccludedInds', 1:size(occlusionCosts,2));
occlusionCosts(inds) = occludedCost;




% COST OF RESURFACING FROM OCCLUSION
nearestResurfaceInds = knnsearch(gridPts, xy);
resurfaceCosts = zeros(numCurrent, numOccluded);
inds = sub2ind(size(resurfaceCosts), (1:numCurrent)', nearestResurfaceInds);
resurfaceCosts(inds) = occludedCost;


% resurfaceCosts = nan(numCurrent, numOccluded);
% 
% for i=1:numCurrent
%     distances = sqrt(sum((gridPts - repmat(xy(i,:),numOccluded,1)).^2, 2));
%     distances(distances > maxVelocity) = nan;
%     resurfaceCosts(i,:) = (max(distances) - distances) * ( occludedCost / max(distances));
%     resurfaceCosts(isnan(resurfaceCosts)) = 0;
% end



% COST OF REMAINING OCCLUDED
stayOccludedCosts = eye(numOccluded, numOccluded) * occludedCost;




% PUTTING IT ALL TOGETHER

pairwisePotentials = ([transitionCosts, resurfaceCosts; occlusionCosts, stayOccludedCosts]);






