

% performs paw tracking

% settings
session = 'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171130_000\';
classifierBot = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\pawBotMarker.mat';
classifierTop = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\pawTopMarker.mat';
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
tic
close all
potentialLocationsBot = getPotentialLocationsBot(vidBot, modelBot, subHgt, subWid, obsPixPositions, startFrame, false);
save([session 'tracking\potentialLocationsBot.mat'], 'potentialLocationsBot');
toc

% track hind!!!
% close all
% load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\hind.mat', 'model', 'subHgt', 'subWid');
% potentialLocationsBot = getPotentialLocationsBot(vidBot, model, subHgt, subWid, obsPixPositions, startFrame, true);
% save([session 'tracking\potentialLocationsBot.mat'], 'potentialLocationsBot');

%% get locations for bottom
locationsBot = getLocationsBot(potentialLocationsBot, frameTimeStamps, vidBot.Width, vidBot.Height);
save([session 'tracking\locationsBot.mat'], 'locationsBot');
showLocations(vidBot, potentialLocationsBot, locationsBot, showPotentialLocations, .02, anchorPtsBot, 200000);

%% get potential locations for top
tic
potentialLocationsTop = getPotentialLocationsTop(vidTop, locationsBot, xLinearMapping, modelTop, subHgt, subWid, startFrame, true);
toc
save([session 'tracking\potentialLocationsTop.mat'], 'potentialLocationsTop');

%% get locations for top
locationsTop = getLocationsTop(potentialLocationsTop, locationsBot, frameTimeStamps, fs);
showLocations(vidTop, potentialLocationsTop, locationsTop, showPotentialLocations, .02, anchorPtsBot, 2155430);






