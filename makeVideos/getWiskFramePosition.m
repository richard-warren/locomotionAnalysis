function [yWiskPos, xWiskPos] = getWiskFramePosition(frameTop, frameWisk, wiskScaling)

% !!! need to document


% rescale wisk view
frameWisk = imresize(frameWisk, wiskScaling);

% perform cross correlation
ccorr = normxcorr2(frameWisk, frameTop);

% get coordinates in frameTop whether top left corner of frameWisk should live
[~, imax] = max(abs(ccorr(:)));
[yWiskPos, xWiskPos] = ind2sub(size(ccorr),imax(1));
xWiskPos = xWiskPos-size(frameWisk,2);
yWiskPos = yWiskPos-size(frameWisk,1);