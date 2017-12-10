

% performs paw tracking

% settings
session = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\171202_000\';
classifierBot = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\svm\classifiers\pawBot.mat';
classifierTop = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\svm\classifiers\pawTop.mat';
xMapping = 'C:\Users\LindseyBuckingham\Desktop\github\locomotionAnalysis\xAlignment\xLinearMapping.mat';
trainingData = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\svm\trainingData\pawBot\labeledFeatures.mat';
showPotentialLocations = true;
fs = 250;

% initializations
load(trainingData, 'features', 'labels', 'subFrameSize')
load(xMapping, 'xLinearMapping');
load(classifierTop, 'model');
modelTop = model;
load(classifierBot, 'model');
modelBot = model;
load([session 'runAnalyzed.mat'], 'obsPixPositions', 'frameTimeStamps', 'rewardTimes')
if ~exist([session '\tracking'], 'dir'); mkdir([session '\tracking']); end
vidBot = VideoReader([session '\runBot.mp4']);
vidTop = VideoReader([session '\runTop.mp4']);
anchorPtsBot = {[0 0], [0 1], [1 0], [1 1]};



%% hand label paw locations

vidFile = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\171202_000\runBot.mp4';
vid = VideoReader(vidFile);
frameInds = find(obsPixPositions>1 & obsPixPositions<vid.Width);
labelPawLocations('C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\171202_000\runBot.mp4', frameInds, 500);

%% create labeled set
makeLabeledSet('pawBot',...
               'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runBotHandLabeledLocations.mat', ...
               'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\171202_000\runBot.mp4',...
               obsPixPositions)

viewTrainingSet('pawBot');

%% create complex labeled set (paw vs not paw)

makeLabeledSetComplex('pawBotComplex',...
                      'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runBotHandLabeledLocations.mat', ...
                      'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\171202_000\runBot.mp4',...
                      obsPixPositions)

viewTrainingSet('pawBotComplex');

%% create complex labeled set (paw vs other)

paws = 1; % which paw to classify against other paws and negative examples

makeLabeledSetComplex('leftHindBot',...
                      'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runBotHandLabeledLocations.mat', ...
                      'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\171202_000\runBot.mp4',...
                      obsPixPositions, paws)

viewTrainingSet('leftHindBot');


%% train svm

% train bot svm and reload
trainSVM('pawBot');
load(classifierBot, 'model', 'subFrameSize');
modelBot = model;
close all; figure; imagesc(reshape(-model.Beta, subFrameSize(1), subFrameSize(2)))

%% get potential locations for bottom

scoreThresh = 0;
showTracking = true;

fprintf('getting potential bottom locations...\n')
close all
frameInds = find(~isnan(obsPixPositions));
tic; potentialLocationsBot = getPotentialLocationsBot(vidBot, modelBot, features, labels, scoreThresh, subFrameSize,...
                                                      obsPixPositions, frameInds, showTracking);
save([session 'tracking\potentialLocationsBot.mat'], 'potentialLocationsBot');
fprintf('analysis time: %i minutes\n', toc/60)


%% get locations for bottom

locationsBot = getLocationsBot(potentialLocationsBot, frameTimeStamps, vidBot.Width, vidBot.Height, frameInds);
save([session 'tracking\locationsBot.mat'], 'locationsBot');
showLocations(vidBot, frameInds, potentialLocationsBot, fixTracking(locationsBot), showPotentialLocations, .02, anchorPtsBot);


%% get potential locations for top
tic
potentialLocationsTop = getPotentialLocationsTop(vidTop, locationsBot, xLinearMapping, modelTop, subHgt, subWid, startFrame, false);
toc
save([session 'tracking\potentialLocationsTop.mat'], 'potentialLocationsTop');

%% get locations for top
locationsTop = getLocationsTop(potentialLocationsTop, locationsBot, obsPixPositions, frameTimeStamps, fs);
showLocations(vidTop, potentialLocationsTop, locationsTop, showPotentialLocations, .02, anchorPtsBot, 200000);

%% make tracking vid
frameInds = find(~isnan(obsPixPositions));
makeTrackingVid(session, frameInds)




