%% initialize new training data structure

% settings
trainingSetFolder = 'topBotCat';
view = 'both';
sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003', ...
            '180125_001', '180125_002', '180125_003'};
frameNum = 1000;
obsPortion = .5; % portion of trials to include obstacle



trainingSetDir = [getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetFolder '\'];
if ~exist(trainingSetDir, 'dir'); mkdir(trainingSetDir); end

if ~exist([trainingSetDir 'trainingData.mat'], 'file')
    trainingData = createTrainingDataStruct(sessions, frameNum, obsPortion);
    save([trainingSetDir 'trainingData.mat'], 'trainingData', 'view')
else
    fprintf('%s already exists... did not create file\n', trainingSetFolder);
end

%% label things

trainingSetFolder = 'topBotCat';
trainingSetName = 'trainingData.mat';

trainingSetDir = [getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetFolder '\'];
features = {'paw1', 'paw2', 'paw3', 'paw4', 'gen', 'tailBase', 'tailMid', 'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH', 'tailBaseTop', 'tailMidTop'}; % with top view
labelFrames(trainingSetDir, trainingSetName, features);


%% add incorrect frames from tracked vid

session = '180123_000';
trainingSetFolder = 'topBotCat';
trainingSetName = 'trainingData.mat';

trainingDataPath = [getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetFolder '\' trainingSetName];
showTrackingDLC(session, .02, trainingDataPath)



%% prepare data for deepLabCut

trainingSetFolder = 'topBotCat';
scaling = 0.5;
trainingSetName = 'trainingData.mat';
features = {'paw1', 'paw2', 'paw3', 'paw4', 'gen', 'tailBase', 'tailMid', 'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH', 'tailBaseTop', 'tailMidTop'};

trainingDataPath = [getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetFolder '\' trainingSetName];
load(trainingDataPath, 'trainingData', 'view')
writeDir = [getenv('TRAININGEXAMPLESDIR') 'deepLabCut\' trainingSetFolder  'Scaling' num2str(scaling) '\']; % images will be written here
if ~exist(writeDir, 'dir'); mkdir(writeDir); end
prepareTrainingImages(writeDir, trainingData, view, features, scaling);


