
implay('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\wiskTest\runBot.mp4', 50);
implay('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\wiskTest\runWisk.mp4', 50);

%%

spikeAnalysis('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\', 'wiskTest');
load('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\wiskTest\runAnalyzed.mat')

%%

makeVidWisk('wiskTest', [.25 .445], .1, .1);

%% test alignment

wiskScaling = 18/52;

% ind1 = 16091;
% ind2 = 16091;
ind1 = 15000;
ind2 = ind1;

vidTop = VideoReader('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\wiskTest\runTop.mp4');
vidWisk = VideoReader('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\wiskTest\runWisk.mp4');

top = rgb2gray(read(vidTop, ind1));
wisk = rgb2gray(read(vidWisk, ind2));
wisk = imresize(wisk, wiskScaling);

c = normxcorr2(wisk, top);
[max_c, imax] = max(abs(c(:)));
[ypeak, xpeak] = ind2sub(size(c),imax(1));
xpeak = xpeak-size(wisk,2);
ypeak = ypeak-size(wisk,1);
%%
topNew = [top, zeros(size(top,1), size(wisk,2))];
topInds = {ypeak:(ypeak+size(wisk,1)-1), xpeak:(xpeak+size(wisk,2)-1)};
% topNew(topInds{1}, topInds{2}) = (topNew(topInds{1}, topInds{2}) + wisk) / 2;
topNew(topInds{1}, topInds{2}) = wisk;

close all; figure; imshow(topNew); pimpFig



% close all; figure; subplot(2,1,1); imshow(top); subplot(2,1,2); imshow(wisk);




