

% load files
vid = VideoReader('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\wiskTest5\runWisk.mp4');
showTracking = true;
load('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\wiskTest5\runAnalyzed.mat',...
     'obsOnTimes', 'obsOffTimes', 'frameTimeStampsWisk', 'wiskTouchSignal')

%% track wisk contacts
[isWiskTouching, contactPixels] = getWiskContacts(vid, showTracking, frameTimeStampsWisk, obsOnTimes, obsOffTimes);


%%
validInds = ~isnan(wiskTouchSignal);
normedTemp = zscore(wiskTouchSignal(validInds));
normed = nan(size(wiskTouchSignal));
normed(validInds) = normedTemp;

figure; yyaxis left; plot(wiskTouchSignal); yyaxis right; plot(normed); pimpFig
% close all; figure; plot(normed)


