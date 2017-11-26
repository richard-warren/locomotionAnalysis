

% performs paw tracking

% settings
session = 'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171126_000\';
classifier = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\pawBot.mat';
showTracking = false;

% initializations
load(classifier, 'model', 'subHgt', 'subWid')
load([session 'runAnalyzed.mat'], 'obsPixPositions', 'frameTimeStamps', 'rewardTimes')
if ~exist([session '\tracking'], 'dir'); mkdir([session '\tracking']); end
vidBot = VideoReader([session '\runBot.mp4']);
startFrame = find(frameTimeStamps>rewardTimes(1));



%% perform tracking
tic
potentialLocationsBot = getPotentialLocationsBot(vidBot, model, subHgt, subWid, obsPixPositions, startFrame, true);
toc