function wheelPoints = getWheelPoints(session)

% given a top view of the mouse running, returns three points along the top of the wheel
% does this by collecting a bunch of random frames from the vid, then taking the minimum projection across frames
% this preserves the wheel but hide the mouse's limb because they move around
% then treshold this image and get rid of connected regions that dont interest the bottom of the frame (as the wheel does)
% this leaves a blob containing only the wheel
% then for three xLocations find the topmost white pixel in the x colums and return wheelPoints


% settings
xLocations = [.1 .5 .9];  % x locations of circRoiPoints, expressed as fraction of vid width
frameNum = 100;  % number of frames to sample
minFrame = 250*5*60;  % start looking for frames this many seconds in (use to avoid periods before the mouse is put on the wheel)
threshFactor = .5; % threshFactor times average pixel vale


% initializations
try
    vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4'));
catch
    vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
end
wheelPoints = nan(3,2); % each row is x,y measured from top left of image
wheelPoints(:,1) = round(xLocations * vid.Width)';
numberOfFrames = floor(vid.Duration*vid.FrameRate);

% get stack of random frames
if minFrame>=numberOfFrames; minFrame = 0; end % a hack that allows me to analyze really short test sessions
frameInds = sort(randperm(numberOfFrames-minFrame, frameNum) + minFrame);
frames = nan(vid.Height, vid.Width, frameNum);
for i = 1:length(frameInds)
    frames(:,:,i) = rgb2gray(read(vid, frameInds(i)));
end
thresh = mean(frames(:))*threshFactor;

% get minimum project across all frames (this effectively removes the legs from above the wheel)
minProjection = uint8(min(frames, [], 3));


% threshold image and keep only connected regions that intersects the bottom row of image
threshed = minProjection > thresh;
threshed(end,:) = 0;  % this is a hack that allows the function to work with only the top view, or with top and bottom views concatenated
% [~, botRow] = min(diff(sum(threshed,2)));  % !!! switch to this line if switched completely away from unconcatenated vids // this finds the bottom row of the top view // the row after this has a sharp drop in brightness, because the wheel gets cut off abruptly
botRow = find(abs(diff(sum(threshed,2)))>60, 1, 'first');  % this finds the bottom row of the top view // the row after this has a sharp drop in brightness, because the wheel gets cut off abruptly
threshed = threshed(1:botRow, :);  % restrict to top camera view
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














