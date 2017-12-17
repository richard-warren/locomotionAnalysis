function frameWithObs = addObsToFrame(frame, obsPos, obsThickness, yLims, color)

% !!! need to document
% generally adds a rectangle to frame to cover obstacle, either for masking or visualization purposes

% initializations
frameWithObs = frame;
obsPixLeft = floor(obsThickness/2);
obsPixRight = ceil(obsThickness/2);
yInds = yLims(1):yLims(2);
obsPixMinMax = round([obsPos - obsPixLeft, obsPos + obsPixRight - 1]);


if any(obsPixMinMax>0 & obsPixMinMax<=size(frame,2))
    
    obsPixMinMax(obsPixMinMax<1) = 1;
    obsPixMinMax(obsPixMinMax>size(frame,2)) = size(frame,2);
    frameWithObs(yInds, obsPixMinMax(1):obsPixMinMax(2)) = color;
    
end