

% performs paw tracking

% settings
session = 'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171126_000\';
classifier = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\pawBot.mat';
showPotentialLocations = true;

% initializations
load(classifier, 'model', 'subHgt', 'subWid')
load([session 'runAnalyzed.mat'], 'obsPixPositions', 'frameTimeStamps', 'rewardTimes')
if ~exist([session '\tracking'], 'dir'); mkdir([session '\tracking']); end
vidBot = VideoReader([session '\runBot.mp4']);
startFrame = find(frameTimeStamps>rewardTimes(1), 1, 'first');
anchorPtsBot = {[0 0], [0 1], [1 0], [1 1]};



%% perform tracking

potentialLocationsBot = getPotentialLocationsBot(vidBot, model, subHgt, subWid, obsPixPositions, startFrame, false);
save([session 'tracking\potentialLocationsBot.mat'], 'potentialLocationsBot');
%%
locationsBot = getLocationsBot(potentialLocationsBot, frameTimeStamps, vidBot.Width, vidBot.Height);
save([session 'tracking\locationsBot.mat'], 'locationsBot');
%%
showLocations(vidBot, potentialLocationsBot, fixTracking(locationsBot), showPotentialLocations, .02, anchorPtsBot, startFrame);