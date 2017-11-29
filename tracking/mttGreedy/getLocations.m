

% performs paw tracking

% settings
session = 'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171129_000\';
classifierBot = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\pawBot.mat';
classifierTop = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\pawTop.mat';
xMapping = 'C:\Users\rick\Desktop\github\locomotionAnalysis\xAlignment\xLinearMapping.mat';
showPotentialLocations = true;
fs = 250;

% initializations
load(xMapping, 'xLinearMapping');
load(classifierTop, 'model', 'subHgt', 'subWid');
modelTop = model;
load(classifierBot, 'model', 'subHgt', 'subWid');
modelBot = model;
load([session 'runAnalyzed.mat'], 'obsPixPositions', 'frameTimeStamps', 'rewardTimes')
if ~exist([session '\tracking'], 'dir'); mkdir([session '\tracking']); end
vidBot = VideoReader([session '\runBot.mp4']);
vidTop = VideoReader([session '\runTop.mp4']);
startFrame = find(frameTimeStamps>rewardTimes(1), 1, 'first');
anchorPtsBot = {[0 0], [0 1], [1 0], [1 1]};



% perform tracking

%% get potential locations for bottom
potentialLocationsBot = getPotentialLocationsBot(vidBot, model, subHgt, subWid, obsPixPositions, startFrame, false);
save([session 'tracking\potentialLocationsBot.mat'], 'potentialLocationsBot');

%% get locations for bottom
locationsBot = getLocationsBot(potentialLocationsBot, frameTimeStamps, vidBot.Width, vidBot.Height);
save([session 'tracking\locationsBot.mat'], 'locationsBot');
showLocations(vidBot, potentialLocationsBot, locationsBot, showPotentialLocations, .02, anchorPtsBot, startFrame*5);

%% get potential locations for top
potentialLocationsTop = getPotentialLocationsTop(vidTop, locationsBot, xLinearMapping, model, subHgt, subWid, startFrame, false);
save([session 'tracking\potentialLocationsTop.mat'], 'potentialLocationsTop');

%% get locations for top
locationsTop = getLocationsTop(potentialLocationsTop, locationsBot, frameTimeStamps, fs);
showLocations(vidTop, potentialLocationsTop, fixTracking(locationsTop), showPotentialLocations, .02, anchorPtsBot, startFrame*5);






