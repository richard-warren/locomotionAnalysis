

% performs paw tracking

% settings
session = 'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\';
xMapping = 'C:\Users\rick\Desktop\github\locomotionAnalysis\xAlignment\xLinearMapping.mat';

% initializations
load(xMapping, 'xLinearMapping');
load([session 'runAnalyzed.mat'], 'obsPixPositions', 'frameTimeStamps', 'rewardTimes')
anchorPtsBot = {[0 0], [0 1], [1 0], [1 1]};



%% hand label paw locations

vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\runBot.mp4';
vid = VideoReader(vidFile);
frameInds = find(obsPixPositions>1 & obsPixPositions<vid.Width);
labelPawLocations('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\runBot.mp4', frameInds, 500);

%% create labeled set, paw vs other

posEgs = 500;
negEgsPerEg = 5;
subFrameSize = [50 50];
includeLocation = false;
paws = 1:4;
class = 'pawBot';

makeLabeledSet(class,...
               'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runBotHandLabeledLocations.mat', ...
               'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\runBot.mp4',...
               subFrameSize, obsPixPositions, posEgs, negEgsPerEg, includeLocation, paws)

viewTrainingSet(class);


%% train svm

class = 'pawBot';

% train bot svm and reload
trainSVM(class);
load(['C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\' class], 'model', 'subFrameSize');

close all; figure; imagesc(reshape(-model.Beta, subFrameSize(1), subFrameSize(2)))

%% train second classifier

class = 'pawBot';

% train bot svm and reload
trainSecondClassifier(class); % this is saved as [class '2']


%% get potential locations for bottom


% settings
scoreThresh = 0;
showTracking = false;
model1 = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\pawBot';
model2 = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\pawBot2';

% initializations
load(model1, 'model', 'subFrameSize');
model1 = model; subFrameSize1 = subFrameSize;
load(model2, 'model', 'subFrameSize');
model2 = model; subFrameSize2 = subFrameSize;
vidBot = VideoReader([session '\runBot.mp4']);

fprintf('getting potential bottom locations...\n')
close all
frameInds = find(~isnan(obsPixPositions));
tic; potentialLocationsBot = getPotentialLocationsBot(vidBot, model1, model2, subFrameSize1, subFrameSize2,...
                                                      scoreThresh, obsPixPositions, frameInds, showTracking);
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




