function unaryPotentials = getUnaryPotentials(x, y, frameWidth, frameHeight, anchorPointX, anchorPointY)

% !!! need to document

dx = (x / frameWidth) - anchorPointX;
dy = (y / frameHeight) - anchorPointY;
unaryPotentials = 1 - sqrt(dx.^2 + dy.^2);




