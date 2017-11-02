function unaryPotentials = getUnaryPotentials(x, y, xyScores, frameWidth, frameHeight, anchorPointY, anchorPointX, numOccluded, unaryWeight)

% computes the prior likelihood that a paw exists in a given location
% simply computse distance of the paw to the anchorPoint (x,y)
% the likelihood is the inverse of the distance to the anchor point
% width and height are normalized from 0 to 1
%
% inputs        x:             vector of x locations of tracks
%               y:             vector of x locations of tracks
%               frameHeight:   height of frame (used to normalize y from 0 to 1)
%               frameWidth:    width of frame (used to normalize x from 0 to 1)
%               anchorPointX:  most likely x position of paw, expressed from 0 to 1
%               anchorPointY:  most likely y position of paw, expressed from 0 to 1
%               numOccluded:   number of occluded states (zeros will be added for the unary potentials of all occluded states)


% compute distant of points to anchorPoint
dx = (x / frameWidth) - anchorPointX;
dy = (y / frameHeight) - anchorPointY;
unaryPotentials = sqrt(2) - sqrt(dx.^2 + dy.^2); % sqrt(2) is the maximum possible distance, eg the distance from one corner to the opposite corner
unaryPotentials(abs(dx) > .75) = 0;

% multiply by location scores
% unaryPotentials = unaryPotentials .* xyScores;

% add zeros for occluded points
unaryPotentials = [unaryPotentials; zeros(numOccluded, 1)] * unaryWeight;
