

% performs paw tracking

% settings
session = 'markerTest1';

% initializations
xMapping = [getenv('GITDIR') 'locomotionAnalysis\xAlignment\xLinearMapping.mat'];
load(xMapping, 'xLinearMapping');
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsPixPositions', 'frameTimeStamps', 'rewardTimes')
frameInds = find(~isnan(obsPixPositions));
anchorPtsBot = {[0 0], [0 1], [1 0], [1 1]};



%% hand label paw bot locations

vidFile = [getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4'];
labelPawLocations(vidFile, frameInds, 50);


%% create bot labeled set, paw vs other

posEgs = 400;
negEgsPerEg = 5;
subFrameSize = [45 45];
includeLocation = false;
paws = 1:4;
class = 'pawBot';

makeLabeledSet(class,...
               [getenv('OBSDATADIR') 'sessions\' session '\tracking\runBotHandLabeledLocations.mat'], ...
               [getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4'],...
               subFrameSize, obsPixPositions, posEgs, negEgsPerEg, includeLocation, paws)

viewTrainingSet(class);


%% train bot svm1

class = 'pawBot';

% train bot svm and reload
trainSVM(class);
load(['C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\' class], 'model', 'subFrameSize');

close all; figure; imagesc(reshape(-model.Beta, subFrameSize(1), subFrameSize(2)))

%% train bot second classifier

class = 'pawBot';

% train bot svm and reload
trainSecondClassifier(class); % this is saved as [class '2']


%% get bot potential locations


% settings
scoreThresh = 0;
showTracking = false;
model1 = [getenv('OBSDATADIR') 'svm\classifiers\pawBot'];
model2 = [getenv('OBSDATADIR') 'svm\classifiers\pawBot2'];

% initializations
load(model1, 'model', 'subFrameSize');
model1 = model; subFrameSize1 = subFrameSize;
load(model2, 'model', 'subFrameSize');
model2 = model; subFrameSize2 = subFrameSize;
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);

tic; potentialLocationsBot = getPotentialLocationsBot(vidBot, model1, model2, subFrameSize1, subFrameSize2,...
                                                      scoreThresh, obsPixPositions, frameInds, showTracking);
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsBot.mat'], 'potentialLocationsBot');
fprintf('potential locations bot analysis time: %i minutes\n', toc/60)


%% get bot locations

% settings
showPotentialLocations = true;

% initializations
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);

locationsBot = getLocationsBot(potentialLocationsBot, frameTimeStamps, vidBot.Width, vidBot.Height, frameInds);
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBot.mat'], 'locationsBot');
showLocations(vidBot, frameInds, potentialLocationsBot, fixTracking(locationsBot), showPotentialLocations, .02, anchorPtsBot);


%% hand label top locations

vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\runTop.mp4';
vid = VideoReader(vidFile);
labelPawLocations(vidFile, frameInds, 200);


%% create top labeled set, paw vs other

posEgs = 400;
negEgsPerEg = 5;
subFrameSize = [40 40];
includeLocation = false;
paws = 1:4;
class = 'pawTop';

makeLabeledSet(class,...
               'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\tracking\runTopHandLabeledLocations.mat', ...
               'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171202_000\runTop.mp4',...
               subFrameSize, obsPixPositions, posEgs, negEgsPerEg, includeLocation, paws)

viewTrainingSet(class);


%% train top svm1

class = 'pawTop';

% train bot svm and reload
trainSVM(class);
load([getenv('OBSDATADIR') 'svm\classifiers\' class], 'model', 'subFrameSize');

close all; figure; imagesc(reshape(-model.Beta, subFrameSize(1), subFrameSize(2)))


%% get potential locations for top

% settings
scoreThresh = 0;
showTracking = true;
model = [getenv('OBSDATADIR') 'svm\classifiers\pawTop'];
paws = [2 4];

% initializations
load(model, 'model', 'subFrameSize');
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);

tic; potentialLocationsTop = getPotentialLocationsTop(vidTop, locationsBot, xLinearMapping, model, subFrameSize, scoreThresh,...
                                                      frameInds, paws, showTracking);
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsTop.mat'], 'potentialLocationsTop');
fprintf('potential locations top analysis time: %i minutes\n', toc/60)

%% get locations for top

% settings
showPotentialLocations = true;
paws = [2 4];
fs = 250;

locationsTop = getLocationsTop(potentialLocationsTop, locationsBot, xLinearMapping, frameInds, obsPixPositions, frameTimeStamps, paws, fs);
showLocations(vidTop, frameInds, potentialLocationsTop, (locationsTop), showPotentialLocations, .02, anchorPtsBot);
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsTop.mat'], 'locationsTop');

%% make tracking vid

makeTrackingVid(session, frameInds)




