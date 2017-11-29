
% settings
rate = 8; % one of the following: 2 4 8 16 32

% initializations
anchorPtsBot = {[0 0], [0 1], [1 0], [1 1]};

vid = VideoReader(['C:\Users\rick\Google Drive\columbia\obstacleData\compressionTests\' num2str(rate) '\runBot.mp4']);
locationsBot.x = nan(vid.NumberOfFrames, 4);
locationsBot.y = nan(vid.NumberOfFrames, 4);
showLocations(vid, [], locationsBot, false, .02, anchorPtsBot, 1);