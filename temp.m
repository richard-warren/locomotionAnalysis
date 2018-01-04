

vidOld = VideoReader('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\markerTest1\runBot.mp4');
vidNew = VideoReader('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\180102_002\runBot.mp4');

%%
frameOld = read(vidOld, 10000);
frameNew = read(vidNew, 10000);


%%

offset = uint8(mean(frameNew(:)) - mean(frameOld(:)));
scale = mean(frameNew(:)) / mean(frameOld(:));


%%
close all; figure;
subplot(2,1,1)
imshow(imadjust(frameOld, [.04 .3], [0 1]))
subplot(2,1,2)
imshow(frameNew)
pimpFig;