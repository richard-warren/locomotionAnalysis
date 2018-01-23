function unaryPotentials = getUnaryPotentialsViterbi(x, y, frameWidth, frameHeight, anchorX, anchorY, maxDistanceX, maxDistanceY)

% takes the distance of x and y from anchor points as fraction of frame widht and height
% then normalizes st 0 corresponds to anything >= maxDistanceX and maxDistanceY and 1 is no distance
% finally averages xScore and yScore to get unary potentials
%
% inputs        x:             vector of x locations of tracks
%               y:             vector of x locations of tracks
%               frameHeight:   height of frame (used to normalize y from 0 to 1)
%               frameWidth:    width of frame (used to normalize x from 0 to 1)
%               anchorPointX:  most likely x position of paw, expressed from 0 to 1
%               anchorPointY:  most likely y position of paw, expressed from 0 to 1
%               maxDistance:   scores are set to 0 if ddistances between anchor point and xy point is greater than max distance


% compute distant of points to anchorPoint
dx = abs((x / frameWidth) - anchorX);
dy = abs((y / frameHeight) - anchorY);

xScore = dx;
yScore = dy;

xScore(xScore>maxDistanceX) = maxDistanceX;
yScore(yScore>maxDistanceY) = maxDistanceY;

xScore = (maxDistanceX - xScore) / maxDistanceX;
yScore = (maxDistanceY - yScore) / maxDistanceY;

unaryPotentials = (xScore + yScore) / 2;
invalidPositions = xScore==0 | yScore==0;
unaryPotentials(invalidPositions) = 0;




