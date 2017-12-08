function frameMasked = maskObs(frame, obsPixPosition)

% !!! need to document

% settings
obsMaskDepth = 0;
obsMaskWid = 18;

% initializations
obsPixLeft = floor(obsMaskWid/2);
obsPixRight = ceil(obsMaskWid/2);

obsPixMinMax = round([obsPixPosition - obsPixLeft, obsPixPosition + obsPixRight - 1]);

if any(obsPixMinMax>0 & obsPixMinMax<=size(frame,2))
    obsPixMinMax(obsPixMinMax<1) = 1;
    obsPixMinMax(obsPixMinMax>size(frame,2)) = size(frame,2);
    frame(:,obsPixMinMax(1):obsPixMinMax(2)) = frame(:,obsPixMinMax(1):obsPixMinMax(2)) .* obsMaskDepth;
end

frameMasked = frame;


