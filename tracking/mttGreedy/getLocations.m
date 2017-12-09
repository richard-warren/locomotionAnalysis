

% performs paw tracking

% settings
session = 'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\';
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
anchorPtsBot = {[0 0], [0 1], [1 0], [1 1]};



%% hand label paw locations

vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\runBot.mp4';
vid = VideoReader(vidFile);
frameInds = find(obsPixPositions>1 & obsPixPositions<vid.Width);
labelPawLocations('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\runBot.mp4', frameInds, 500);

%% create labeled set
makeLabeledSet('pawBot',...
               'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runBotHandLabeledLocations.mat', ...
               'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\runBot.mp4',...
               obsPixPositions)

viewTrainingSet('pawBot');

%% train svm

% train bot svm and reload
egPortion = .1; % only use this portion of the positive and negative examples
trainSVM('pawBot', egPortion);
load(classifierBot, 'model', 'subHgt', 'subWid');
modelBot = model;
close all; figure; imagesc(reshape(-model.Beta,51,51))

%% get potential locations for bottom

scoreThresh = .5;
showTracking = false;

fprintf('getting potential bottom locations... ')
close all
frameInds = find(~isnan(obsPixPositions));
tic; potentialLocationsBot = getPotentialLocationsBot(vidBot, modelBot, scoreThresh, subHgt, subWid,...
                                                      obsPixPositions, frameInds, showTracking);
save([session 'tracking\potentialLocationsBot.mat'], 'potentialLocationsBot');
fprintf('analysis time: %i minutes\n', toc/60)


%% get locations for bottom

locationsBot = getLocationsBot(potentialLocationsBot, frameTimeStamps, vidBot.Width, vidBot.Height, frameInds);
save([session 'tracking\locationsBot.mat'], 'locationsBot');
showLocations(vidBot, frameInds, potentialLocationsBot, (locationsBot), showPotentialLocations, .02, anchorPtsBot);


%% get potential locations for top
tic
potentialLocationsTop = getPotentialLocationsTop(vidTop, locationsBot, xLinearMapping, modelTop, subHgt, subWid, startFrame, false);
toc
save([session 'tracking\potentialLocationsTop.mat'], 'potentialLocationsTop');

% get locations for top
locationsTop = getLocationsTop(potentialLocationsTop, locationsBot, obsPixPositions, frameTimeStamps, fs);
showLocations(vidTop, potentialLocationsTop, locationsTop, showPotentialLocations, .02, anchorPtsBot, 200000);

%% make tracking vid
frameInds = find(~isnan(obsPixPositions));
makeTrackingVid(session, frameInds)




