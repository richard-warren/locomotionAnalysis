

vid = VideoReader('C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\wiskTest2\runTop.mp4');
vidSub = VideoReader('C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\wiskTest2\runWisk.mp4');
load('C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\wiskTest2\runAnalyzed.mat')

%%
frameInd = 19247;
frame = rgb2gray(read(vid, frameInd));
subInd = find(frameTimeStampsWisk==frameTimeStamps(frameInd), 1, 'first');
frameSub = rgb2gray(read(vidSub, subInd));

%%
figure; subplot(2,1,1); image(frameTop); subplot(2,1,2); image(frameWisk)

getWiskFramePosition(frameTop, frameWisk);
