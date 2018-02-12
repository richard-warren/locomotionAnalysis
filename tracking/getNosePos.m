function [noseX, noseY, medianFrame] = getNosePos(vidBot)

% given bot vid of mouse, gets x and y positions of tip of nose
% does this by getting median of random selection of frames, then thresholding
% largest blob is taken to be mouse outline, and rightmost edge of blob is taken for noseX
% this blob is cropped to the rightmost noseSubFrameWidth pixels, then the y value of the centroid of this blob is taken as the noseY / midline



% settings
minFrame = 10000;  % when getting avg frame, only look for frames beyond this frame
threshFactor = 1.2; % 2 times the mean of the frame
frameNum = 500;    % number of random frames to use for min projection
noseSubFrameWidth = 10;


% get stack of random frames
frameInds = randperm(vidBot.NumberOfFrames-minFrame, frameNum) + minFrame;
frameInds = sort(frameInds);

frames = nan(vidBot.Height, vidBot.Width, frameNum);
for i = 1:length(frameInds)
    frames(:,:,i) = rgb2gray(read(vidBot, frameInds(i)));
end


% get median frame
medianFrame = uint8(median(frames, 3));
thresh = mean(medianFrame(:)) * threshFactor;


% get nose x position (rightmost edge of largest binary region)
threshed = medianFrame > thresh;
mouseOutlineInfo = regionprops(threshed, 'Area', 'BoundingBox', 'FilledImage');
[~, maxInd] = max([mouseOutlineInfo.Area]);
mouseOutlineInfo = mouseOutlineInfo(maxInd);
noseX = floor(mouseOutlineInfo.BoundingBox(1) + mouseOutlineInfo.BoundingBox(3));

% get nose y position / midline (y center of mass of largest binary region cropped to the rightmost noseSubFrameWidth pixels)
subFrame = mouseOutlineInfo.FilledImage(:, end-noseSubFrameWidth:end);
noseInfo = regionprops(subFrame, 'Centroid');
noseY = noseInfo.Centroid(2) + mouseOutlineInfo.BoundingBox(2);
