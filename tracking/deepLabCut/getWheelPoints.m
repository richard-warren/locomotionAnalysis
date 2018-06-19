function wheelPoints = getWheelPoints(vid)

% given a top view of the mouse running, returns three points along the top of the wheel
% does this by collecting a bunch of random frames from the vid, then taking the minimum projection across frames
% this preserves the wheel but hide the mouse's limb because they move around
% then treshold this image and get rid of connected regions that dont interest the bottom of the frame (as the wheel does)
% this leaves a blob containing only the wheel
% then for three xLocations find the topmost white pixel in the x colums and return wheelPoints



% settings
xLocations = [.1 .5 .9]; % x locations of circRoiPoints, expressed as fraction of vid width
frameNum = 100;
minFrame = 250*5*60; % five minutes in
threshFactor = .5; % threshFactor times average pixel vale



% initializations
wheelPoints = nan(3,2); % each row is x,y measured from top left of image
wheelPoints(:,1) = round(xLocations * vid.Width)';



% get stack of random frames
frameInds = randperm(vid.NumberOfFrames-minFrame, frameNum) + minFrame;
frameInds = sort(frameInds);

frames = nan(vid.Height, vid.Width, frameNum);
for i = 1:length(frameInds)
    frames(:,:,i) = rgb2gray(read(vid, frameInds(i)));
end
thresh = mean(frames(:))*threshFactor;

% get minimum project across all frames (this effectively removes the legs from above the wheel)
minProjection = uint8(min(frames, [], 3));



% threshold image and keep only connected regions that interest the bottom row of image
threshed = minProjection > thresh;
regionInfo = bwlabel(threshed);
botRowGroups = unique(regionInfo(end,:)); botRowGroups = botRowGroups(botRowGroups>0);
wheelOutline = ismember(regionInfo, botRowGroups);


% for each x value, find the top of the circle, which is the first white point in the wheelOutline (from the top to bot of column)
for i = 1:3
    x = wheelPoints(i,1);
    try
    wheelPoints(i,2) = find(wheelOutline(:,x), 1, 'first');
    catch; keyboard; end
end














