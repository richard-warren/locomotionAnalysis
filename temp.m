

% load files
vid = VideoReader('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\wiskTest5\runWisk.mp4');
showTracking = true;
load('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\wiskTest5\runAnalyzed.mat',...
     'obsOnTimes', 'obsOffTimes', 'frameTimeStampsWisk', 'wiskTouchSignal')

%% track wisk contacts
[isWiskTouching, contactPixels] = getWiskContacts(vid, showTracking, frameTimeStampsWisk, obsOnTimes, obsOffTimes);

