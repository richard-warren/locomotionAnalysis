% PERFORM PAW TRACKING

% settings
session = '180122_001';
minVel = .4;
obsPrePost = [.2 .2];
velPositions = [-.08 .08] + 0.3820;
anchorPtsBot = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH // each entry is x,y pair measured from top left corner
colors = hsv(4); % red green blue purple

% initializations
trackingDir = [getenv('OBSDATADIR') 'sessions\' session '\tracking'];
if ~exist(trackingDir, 'dir'); mkdir(trackingDir); end % make tracking directory if it doesn't already exist

xMapping = [getenv('GITDIR') 'locomotionAnalysis\xAlignment\xLinearMapping.mat'];
load(xMapping, 'xLinearMapping');

load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsPixPositions', 'frameTimeStamps',...
    'wheelPositions', 'wheelTimes', 'obsPositions', 'obsTimes', 'obsOnTimes', 'obsOffTimes')
obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);

[frameInds, trialVels] = getTrialFrameInds(minVel, obsPrePost, velPositions, frameTimeStamps,...
    wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, obsOffTimes);

vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);



%% hand label paw bot locations

vidFile = [getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4'];
labelPawLocations(vidFile, frameInds, 100, anchorPtsBot, colors);


%% create bot labeled set, svm1

posEgs = 400;
negEgsPerEg = 10;
jitterNum = 0;
jitterPixels = 4;
subFrameSize1 = [45 45];
includeLocation = false;
paws = 1:4;
class = 'pawBot1';
maxOverlap = .5;
minBrightness = 2; % negative examples need to be minBrightness times the mean brightness of the current frame

makeLabeledSet(class,...
               [getenv('OBSDATADIR') 'sessions\' session '\tracking\runBotHandLabeledLocations.mat'], ...
               [getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4'],...
               subFrameSize1, obsPixPositions, posEgs, negEgsPerEg, includeLocation, paws,...
               jitterPixels, jitterNum, maxOverlap, minBrightness);

viewTrainingSet(class);


%% create bot labeled set, second classifier

posEgs = 400;
negEgsPerEg = 10;
jitterNum = 8;
jitterPixels = 2;
subFrameSize2 = [45 45];
includeLocation = 1;
paws = [1 4; 2 3]; % every row is a class // all paws in a row belong to that class (hind vs fore paws, for examlpe)
class = 'pawBot2';
maxOverlap = .5;
% targetSize = [227 227];
minBrightness = 2.5; % negative examples need to be minBrightness times the mean brightness of the current frame

makeLabeledSet(class,...
               [getenv('OBSDATADIR') 'sessions\' session '\tracking\runBotHandLabeledLocations.mat'], ...
               [getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4'],...
               subFrameSize2, obsPixPositions, posEgs, negEgsPerEg, includeLocation, paws,...
               jitterPixels, jitterNum, maxOverlap, minBrightness);

viewTrainingSet(class);

% prepareTrainingData(class, targetSize); % restructure features



%% train bot svm1

class = 'pawBot1';

% train bot svm and reload
trainSVM(class);
load(['C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\' class], 'model', 'subFrameSize');

close all; figure; imagesc(reshape(-model.Beta, subFrameSize(1), subFrameSize(2)))

%% train bot svm2

class = 'pawBot2';

% train bot svm2
trainSVM2(class); % this is saved as [class '2']

%% train cnn

% !!! use the script cnnTest to retrain alexnet on my data...



%% get bot potential locations


% settings
scoreThresh = 0;
showTracking = 1;
model1 = [getenv('OBSDATADIR') 'svm\classifiers\pawBot1'];
classNum = 2; % not included not paw class

% svm1
load(model1, 'model', 'subFrameSize');
model1 = model; subFrameSize1 = subFrameSize;

% svm2
load([getenv('OBSDATADIR') 'svm\classifiers\pawBot2'], 'model', 'subFrameSize');
model2 = model; subFrameSize2 = subFrameSize;

% cnn
% load([getenv('OBSDATADIR') 'svm\classifiers\pawBot2Cnn'], 'netTransfer');
% model2 = netTransfer; clear netTransfer;


tic
potentialLocationsBot = getPotentialLocationsBot(vidBot, model1, model2, classNum, ...
    subFrameSize1, subFrameSize2, scoreThresh, obsPixPositions, frameInds, showTracking);
                                                  
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsBotAll.mat'], 'potentialLocationsBot');
fprintf('potential locations bot analysis time: %i minutes\n', toc/60)


%% get bot locations

% settings
showPotentialLocations = true;

% initializations
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);

locationsBot = getLocationsBot(potentialLocationsBot, anchorPtsBot, frameTimeStamps, vidBot.Width, vidBot.Height, frameInds);
showLocations(vidBot, frameInds, potentialLocationsBot, fixTracking(locationsBot), showPotentialLocations, .02, anchorPtsBot, colors);
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBot.mat'], 'locationsBot');



%% hand label top locations

vidFile = [getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4'];
obsFrameInds = find(obsPixPositions>10 & obsPixPositions<vidTop.Width);
labelPawLocations(vidFile, obsFrameInds, 200);


%% create top labeled set, svm

posEgs = 400;
negEgsPerEg = 10;
jitterNum = 0;
jitterPixels = 2;
subFrameSize1 = [40 40];
includeLocation = false;
paws = 1:4;
class = 'pawTop1';
maxOverlap = .5;
minBrightness = 2.5; % negative examples need to be minBrightness times the mean brightness of the current frame

makeLabeledSet(class,...
               'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runTopHandLabeledLocations.mat', ...
               'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\runTop.mp4',...
               subFrameSize1, obsPixPositions, posEgs, negEgsPerEg, includeLocation, paws, threshIntensity, ...
               jitterPixels, jitterNum, maxOverlap, minBrightness);
viewTrainingSet(class);


%% create top labeled set, cnn

posEgs = 400;
negEgsPerEg = 10;
jitterNum = 8;
jitterPixels = 3;
subFrameSize2 = [40 40];
includeLocation = false;
paws = 1:4;
class = 'pawTop2';
maxOverlap = .5;
targetSize = [227 227];
minBrightness = 2.5; % negative examples need to be minBrightness times the mean brightness of the current frame

makeLabeledSet(class,...
               'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runTopHandLabeledLocations.mat', ...
               'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\runTop.mp4',...
               subFrameSize2, obsPixPositions, posEgs, negEgsPerEg, includeLocation, paws, threshIntensity, ...
               jitterPixels, jitterNum, maxOverlap, minBrightness);
viewTrainingSet(class);

% prepareTrainingData(class, targetSize); % restructure features

%% train top svm1

class = 'pawTop1';

% train bot svm and reload
trainSVM(class);
load([getenv('OBSDATADIR') 'svm\classifiers\' class], 'model', 'subFrameSize');

close all; figure; imagesc(reshape(-model.Beta, subFrameSize(1), subFrameSize(2)))


%% get potential locations for top (svm)

% settings
scoreThresh = 0;
showTracking = 0;
model1 = [getenv('OBSDATADIR') 'svm\classifiers\pawTop1'];
model2 = [getenv('OBSDATADIR') 'svm\classifiers\pawTop2Cnn'];
paws = 1:4;

% initializations
load(model1, 'model', 'subFrameSize');
model1 = model; clear model;
subFrameSize1 = subFrameSize; clear subFrameSize;
load(model2, 'netTransfer')
model2 = netTransfer; clear netTransfer;
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);

%
tic; potentialLocationsTop = getPotentialLocationsTop(vidTop, locationsBot, model1, model2, ...
    subFrameSize1, subFrameSize2, scoreThresh, frameInds, paws, showTracking);
fprintf('potential locations top analysis time: %i minutes\n', toc/60)

save([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsTop.mat'], 'potentialLocationsTop');



%% get potential locations for top (markers)

% settings
markerThresh = 150;
showTracking = false;

% initializations
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);

tic; potentialLocationsTop = getPotentialLocationsTopMarkers(vidTop, frameInds, markerThresh, showTracking);
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsTop.mat'], 'potentialLocationsTop');
fprintf('potential locations top analysis time: %i minutes\n', toc/60)


%% get locations for top (SVM)

% settings
showPotentialLocations = true;
paws = 1:4;
fs = 250;

% initializations
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);

% fix x alignment for bottom view
locationsBotFixed = fixTracking(locationsBot);
locationsBotFixed.x = locationsBotFixed.x*xLinearMapping(1) + xLinearMapping(2);


locationsTop = getLocationsTop(potentialLocationsTop, locationsBotFixed,...
    frameInds, obsPixPositions, frameTimeStamps, paws, fs);
showLocations(vidTop, frameInds, potentialLocationsTop, fixTracking(locationsTop),...
    showPotentialLocations, .08, anchorPtsBot, colors, locationsBotFixed);
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsTop.mat'], 'locationsTop');

%% get locations for top (markers)

% settings
hindOffset = -5;      % markers on hind legs are more anterior that foot pad that is tracked on the bot view, so bot view x values are shifted to the left by hindOffset
foreOffset = 5;

% initializations
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
locationsBotFixed = fixTracking(locationsBot);

locationsTop = getLocationsTopMarkers(potentialLocationsTop, locationsBotFixed, frameTimeStamps, xLinearMapping, frameInds);
showLocations(vidTop, frameInds(), potentialLocationsTop, (locationsTop), false, .02, anchorPtsBot, locationsBotFixed);
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsTop.mat'], 'locationsTop');


%% make tracking vid

makeTrackingVid(session, frameInds)




