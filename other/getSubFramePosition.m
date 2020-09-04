function [yPos, xPos, scaling] = getSubFramePosition(frame, frameSub, scalings)

% !!! need to document, but basically takes frame and frameSub, where
% frameSub is located within frame, but needs to be translated and scaled
% to fit properly // this find the correct scaling and the x y positions
% corresponding to the position in frame where the top left corner of
% frameSub should live

% settings
scalings = .35 : .005 : .45; % all of the scales in this vector are attempted


% iterate through potential scalings
maxCorrs = nan(1, length(scalings));
maxPositions = nan(2, length(maxCorrs));


for i = 1:length(scalings)
    
   frameSubScaled = imresize(frameSub, scalings(i));
   
   % perform cross correlation
    ccorr = normxcorr2(frameSubScaled, frame);

    % get coordinates in frameTop whether top left corner of frameWisk should live
    [maxCorrs(i), imax] = max(abs(ccorr(:)));
    [yPos, xPos] = ind2sub(size(ccorr), imax(1)); % [y;x]
    
    maxPositions(1,i) = xPos - size(frameSubScaled,2);
    maxPositions(2,i) = yPos - size(frameSubScaled,1);
    
end

[~, maxInd] = max(maxCorrs);
xPos = maxPositions(1, maxInd);
yPos = maxPositions(2, maxInd);
scaling = scalings(maxInd);