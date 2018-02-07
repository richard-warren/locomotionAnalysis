
% TRAIN CLASSIFIERS


% settings
session = '180122_001';
minVel = .4;
obsPrePost = [0 0]; % how many m before and after obs on should be included in frameInds
velPositions = [-.08 .08] + 0.3820; % between what obs positions should trial velocity be computed
anchorPtsBot = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH // each entry is x,y pair measured from top left corner
colors = hsv(4); % red green blue purple

% initializations
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsPixPositions', 'frameTimeStamps',...
    'wheelPositions', 'wheelTimes', 'obsPositions', 'obsTimes', 'obsOnTimes', 'obsOffTimes', 'mToPixMapping', 'targetFs')
obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);

frameInds = getTrialFrameInds(minVel, obsPrePost, velPositions, frameTimeStamps,...
    wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, obsOffTimes);


%% hand label paw bot locations

labelPawLocations(session, 'Bot', frameInds, 100, anchorPtsBot, colors);


%% create bot labeled set, svm1

labeledDataSessions = {'171202_000', '180122_001'};
view = 'Bot'; % 'Bot' or 'Top'
posEgs = 500; % these positive egs are spread agross the classes in paws
negEgsPerEg = 10;
jitterNum = 0;
jitterPixels = 0;
flipBot = false;
subFrameSize1 = [45 45];
featureSetting = 'imageOnly';
pawGroups = [1 2 3 4];
class = 'pawBot1';
maxOverlap = .5;
minBrightness = 2; % negative examples need to be minBrightness times the mean brightness of the current frame

makeLabeledSet(class, labeledDataSessions, view,...
               subFrameSize1, posEgs, negEgsPerEg, featureSetting, pawGroups,...
               jitterPixels, jitterNum, maxOverlap, minBrightness, flipBot);

viewTrainingSet(class);


%% create bot labeled set 2, second classifier

labeledDataSessions = {'171202_000', '180122_001'};
view = 'Bot'; % 'Bot' or 'Top'
posEgs = 5000; % these positive egs are spread agross the classes in paws
negEgsPerEg = 10;
jitterNum = 0;
jitterPixels = 0;
flipBot = true;
subFrameSize2 = [45 45];
featureSetting = 'imageOnly';
pawGroups = [1;2;3;4];
class = 'pawBot2';
maxOverlap = .5;
minBrightness = 2.5; % negative examples need to be minBrightness times the mean brightness of the current frame

makeLabeledSet(class, labeledDataSessions, view,...
               subFrameSize2, posEgs, negEgsPerEg, featureSetting, pawGroups,...
               jitterPixels, jitterNum, maxOverlap, minBrightness, flipBot);

viewTrainingSet(class);

%% train bot svm1

class = 'pawBot1';

% train bot svm and reload
trainSVM(class);
load([getenv('OBSDATADIR') 'tracking\classifiers\' class], 'model', 'subFrameSize');

close all; figure; imagesc(reshape(-model.Beta, subFrameSize(1), subFrameSize(2)))

%% train bot svm2

% class = 'pawBot2';
% trainSVM2(class); % this is saved as [class '2']


%% train nn

% class = 'pawBot2';
% trainNn(class);

%% train alexnet cnn for top (transfer learning)

% settings
class = 'pawBot2';
network = 'alexNet';
% targetSize = [227 227]; prepareTrainingDataForCnn(class, targetSize); % uncomment to convert labeledFeatures.mat into folders containing images, a format appropriate for AlexNet
retrainCnn(class, network);


%% hand label top locations

session = '180123_003';

load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsPixPositions');
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
obsFrameInds = find(obsPixPositions>10 & obsPixPositions<vidTop.Width);
anchorPtsTop = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH // each entry is x,y pair measured from top left corner

labelPawLocations(session, 'Top', obsFrameInds, 50, anchorPtsTop, colors);

%% create top labeled set, svm

labeledDataSessions = {'180123_000', '180123_001', '180123_002', '180123_003'};
view = 'Top'; % 'Bot' or 'Top'
posEgs = 400; % these positive egs are spread agross the classes in paws
negEgsPerEg = 10;
jitterNum = 0;
jitterPixels = 0;
flipBot = false;
subFrameSize1 = [35 35];
featureSetting = 'imageOnly';
pawGroups = [1 2]; % top labeled set only distinguihes hind (1) and for (2) paws
class = 'pawTop1';
maxOverlap = .5;
minBrightness = 2; % negative examples need to be minBrightness times the mean brightness of the current frame

makeLabeledSet(class, labeledDataSessions, view,...
               subFrameSize1, posEgs, negEgsPerEg, featureSetting, pawGroups,...
               jitterPixels, jitterNum, maxOverlap, minBrightness, flipBot);

viewTrainingSet(class);


%% create top labeled set, cnn

labeledDataSessions = {'180123_000', '180123_001', '180123_002', '180123_003'};
view = 'Top'; % 'Bot' or 'Top'
posEgs = 400*8; % these positive egs are spread agross the classes in paws
negEgsPerEg = 10;
jitterNum = 8;
jitterPixels = 2;
flipBot = false;
subFrameSize1 = [35 35];
featureSetting = 'imageOnly';
pawGroups = [1 2]; % top labeled set only distinguihes hind (1) and for (2) paws
class = 'pawTop2';
maxOverlap = .5;
minBrightness = 2; % negative examples need to be minBrightness times the mean brightness of the current frame

makeLabeledSet(class, labeledDataSessions, view,...
               subFrameSize1, posEgs, negEgsPerEg, featureSetting, pawGroups,...
               jitterPixels, jitterNum, maxOverlap, minBrightness, flipBot);

viewTrainingSet(class);

%% train top svm1

class = 'pawTop1';

% train bot svm and reload
trainSVM(class);
load([getenv('OBSDATADIR') 'tracking\classifiers\' class], 'model', 'subFrameSize');

close all; figure; imagesc(reshape(-model.Beta, subFrameSize(1), subFrameSize(2)))


%% train alexnet cnn for top (transfer learning)

% settings
class = 'pawTop2';
network = 'alexNet';
% targetSize = [227 227]; prepareTrainingDataForCnn(class, targetSize); % uncomment to convert labeledFeatures.mat into folders containing images, a format appropriate for AlexNet
retrainCnn(class, network);




