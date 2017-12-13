function pairwisePotentials = getPairwisePotentials(x, y, prevX, prevY, dt, maxVel)

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
%               maxDistance:   scores are set to 0 if ddistances between anchor point and xy point is greater than max distance


% get velocities of all xy moving from prevXY
try; velocities = sqrt(sum(([x, y] - repmat([prevX, prevY], length(x), 1)).^2, 2)) / dt; catch; keyboard; end
velocities(velocities>maxVel) = maxVel;

% normalize (s.t. 0 velocities are pairwise potentials of 1 and velocities >= maxVel are pairwise potentials of 0
pairwisePotentials = (maxVel - velocities) / maxVel;
