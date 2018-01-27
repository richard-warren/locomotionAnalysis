% TRAIN CLASSIFIERS

% settings
session = '180122_001';
minVel = .4;
obsPrePost = [.2 .2];
velPositions = [-.08 .08] + 0.3820;
anchorPtsBot = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH // each entry is x,y pair measured from top left corner
colors = hsv(4); % red green blue purple

% initializations
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsPixPositions', 'frameTimeStamps',...
    'wheelPositions', 'wheelTimes', 'obsPositions', 'obsTimes', 'obsOnTimes', 'obsOffTimes', 'mToPixMapping', 'targetFs')
obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);

[frameInds, trialVels] = getTrialFrameInds(minVel, obsPrePost, velPositions, frameTimeStamps,...
    wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, obsOffTimes);


%% hand label paw bot locations

vidFile = [getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4'];
labelPawLocations(vidFile, frameInds, 100, anchorPtsBot, colors);


%% create bot labeled set, svm1

posEgs = 400;
negEgsPerEg = 10;
jitterNum = 0;
jitterPixels = 0;
subFrameSize1 = [45 45];
featureSetting = 'imageOnly';
paws = 1:4;
class = 'pawBot1';
maxOverlap = .5;
minBrightness = 2; % negative examples need to be minBrightness times the mean brightness of the current frame

makeLabeledSet(class,...
               [getenv('OBSDATADIR') 'sessions\' session '\tracking\runBotHandLabeledLocations.mat'], ...
               [getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4'],...
               subFrameSize1, obsPixPositions, posEgs, negEgsPerEg, featureSetting, paws,...
               jitterPixels, jitterNum, maxOverlap, minBrightness);

viewTrainingSet(class);


%% create bot labeled set 2, second classifier

posEgs = 400;
negEgsPerEg = 10;
jitterNum = 8;
jitterPixels = 3;
subFrameSize2 = [45 45];
featureSetting = 'imageOnly';  % !!!

paws = [1 4; 2 3]; % every row is a class // all paws in a row belong to that class (hind vs fore paws, for examlpe)
class = 'pawBot2'; % !!!
maxOverlap = .5;
minBrightness = 2.5; % negative examples need to be minBrightness times the mean brightness of the current frame

makeLabeledSet(class,...
               [getenv('OBSDATADIR') 'sessions\' session '\tracking\runBotHandLabeledLocations.mat'], ...
               [getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4'],...
               subFrameSize2, obsPixPositions, posEgs, negEgsPerEg, featureSetting, paws,...
               jitterPixels, jitterNum, maxOverlap, minBrightness);

viewTrainingSet(class);


%% train bot svm1

class = 'pawBot1';

% train bot svm and reload
trainSVM(class);
load([getenv('OBSDATADIR') 'svm\classifiers\' class], 'model', 'subFrameSize');

close all; figure; imagesc(reshape(-model.Beta, subFrameSize(1), subFrameSize(2)))

% %% train bot svm2
% 
% class = 'pawBot2';
% trainSVM2(class); % this is saved as [class '2']
% 
% 
% %% train nn
% 
% class = 'pawBot2';
% trainNn(class);

%% train alexnet cnn for top (transfer learning)

% settings
class = 'pawBot2';
network = 'alexNet';
% targetSize = [227 227]; prepareTrainingDataForCnn(class, targetSize); % uncomment to convert labeledFeatures.mat into folders containing images, a format appropriate for AlexNet
retrainCnn(class, network);


%% hand label top locations

vidFile = [getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4'];
obsFrameInds = find(obsPixPositions>10 & obsPixPositions<vidTop.Width);
labelPawLocations(vidFile, obsFrameInds, 100, anchorPtsBot, colors);

%% create top labeled set, svm

posEgs = 400;
negEgsPerEg = 10;
jitterNum = 0;
jitterPixels = 0;
subFrameSize1 = [35 35];
featureSetting = 'imageOnly';
paws = 1:4;
class = 'pawTop1';
maxOverlap = .5;
minBrightness = 2.5; % negative examples need to be minBrightness times the mean brightness of the current frame

makeLabeledSet(class,...
               [getenv('OBSDATADIR') 'sessions\' session '\tracking\runTopHandLabeledLocations.mat'], ...
               [getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4'],...
               subFrameSize1, obsPixPositions, posEgs, negEgsPerEg, featureSetting, paws, ...
               jitterPixels, jitterNum, maxOverlap, minBrightness);
viewTrainingSet(class);



%% create top labeled set, cnn

posEgs = 400;
negEgsPerEg = 10;
jitterNum = 8;
jitterPixels = 2;
subFrameSize2 = [35 35];
featureSetting = 'imageOnly';
paws = 1:4;
class = 'pawTop2';
maxOverlap = .5;
targetSize = [227 227];
minBrightness = 2.5; % negative examples need to be minBrightness times the mean brightness of the current frame

makeLabeledSet(class,...
               [getenv('OBSDATADIR') 'sessions\' session '\tracking\runTopHandLabeledLocations.mat'], ...
               [getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4'],...
               subFrameSize2, obsPixPositions, posEgs, negEgsPerEg, featureSetting, paws, ...
               jitterPixels, jitterNum, maxOverlap, minBrightness);
viewTrainingSet(class);



%% train top svm1

class = 'pawTop1';

% train bot svm and reload
trainSVM(class);
load([getenv('OBSDATADIR') 'svm\classifiers\' class], 'model', 'subFrameSize');

close all; figure; imagesc(reshape(-model.Beta, subFrameSize(1), subFrameSize(2)))


%% train alexnet cnn for top (transfer learning)

% settings
class = 'pawTop2';
network = 'alexNet';
% targetSize = [227 227]; prepareTrainingDataForCnn(class, targetSize); % uncomment to convert labeledFeatures.mat into folders containing images, a format appropriate for AlexNet
retrainCnn(class, network);




