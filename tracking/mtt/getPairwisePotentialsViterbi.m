function [pairwisePotentials, invalidTransitions] = getPairwisePotentialsViterbi(xy, xyLast, maxVelocity, velocityWeight, occludedScore, gridPts)

% given track locations in two adjacent frames, computes a matrix representing the likelihood of transitioning
% from every location in previous frame (xyLast) to every location in current frame (xy)
% this likelihood is based on the distance between track (ie the velocity necessary to move from one to
% another track


% initializations
numOccluded = size(gridPts,1);
numLast = size(xyLast,1);
numCurrent = size(xy,1);




% SCORE FOR TRANSITIONING BETWEEN TRACKS

% restructure matrices to allow direct subtration
xyLastReshaped = repelem(xyLast, numCurrent, 1);
xyReshaped = repmat(xy, numLast, 1);

% compute matrix of pairwise distances (distances is size xy by xyLast)
distances = sqrt(sum((xyReshaped - xyLastReshaped).^2, 2));
distances = reshape(distances, size(xy,1), size(xyLast,1));

% normalize range and set invalid transitions to 0
transitionScores = distances;
invalidTransitions = transitionScores > maxVelocity;
transitionScores(invalidTransitions) = maxVelocity;
transitionScores = (maxVelocity - transitionScores) / maxVelocity;

% weight distances by velocityWeight
transitionScores = transitionScores * velocityWeight;




% SCORE FOR BECOMING OCCLUDED

nearestOccludedInds = knnsearch(gridPts, xyLast);
occlusionScores = zeros(numOccluded, numLast);
inds = sub2ind(size(occlusionScores), nearestOccludedInds', 1:size(occlusionScores,2));
occlusionScores(inds) = occludedScore;




% SCORE FOR RESURFACING FROM OCCLUSION
nearestResurfaceInds = knnsearch(gridPts, xy);
resurfaceScores = zeros(numCurrent, numOccluded);
inds = sub2ind(size(resurfaceScores), (1:numCurrent)', nearestResurfaceInds);
resurfaceScores(inds) = occludedScore;


% SCORE FOR REMAINING OCCLUDED
stayOccludedScores = eye(numOccluded, numOccluded) * occludedScore;




% PUTTING IT ALL TOGETHER

pairwisePotentials = ([transitionScores, resurfaceScores; occlusionScores, stayOccludedScores]);






