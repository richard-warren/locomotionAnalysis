

session = '180122_003';

load(['C:\Users\rick\Google Drive\columbia\obstacleData\sessions\' session '\tracking\locationsBotCorrected.mat'], 'locations')
totalTrackingFrames = size(locations.locationsCorrected,1);

vid = VideoReader(['C:\Users\rick\Google Drive\columbia\obstacleData\sessions\' session '\runBot.mp4']);
totalFrames = vid.NumberOfFrames;

fprintf('frames: %i     tracked: %i\n', totalTrackingFrames, totalFrames)

