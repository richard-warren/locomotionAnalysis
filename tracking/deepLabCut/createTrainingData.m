%% initialize new training data structure

% settings
trainingSetName = 'newTrainingData';
sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003', ...
            '180125_001', '180125_002', '180125_003'};
frameNum = 1000;
obsPortion = .5; % portion of trials to include obstacle


trainingSetDir = [getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetName '\'];

if ~exist(trainingSetDir, 'dir'); mkdir(trainingSetDir); end

if ~exist([trainingSetDir 'trainingData.mat'], 'file')
    trainingData = createTrainingDataStruct(sessions, frameNum, obsPortion);
    save([trainingSetDir 'trainingData.mat'], 'trainingData')
else
    fprintf('%s already exists... did not create file\n', trainingSetName);
end

%% label things

trainingSetName = 'newTrainingData';

trainingSetDir = [getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetName '\'];
features = {'pawTL', 'pawTR', 'pawBR', 'pawBL', 'gen', 'tailBase', 'tailMid', 'tailEnd'};
labelFrames(trainingSetDir, features);



%% prepare data for deepLabCut
trainingSetName = 'newTrainingData';
features = {'pawTL', 'pawTR', 'pawBR', 'pawBL','gen', 'tailBase', 'tailMid', 'tailEnd'}; % excluding genitals

trainingSetDir = [getenv('TRAININGEXAMPLESDIR') 'deepLabCut\' trainingSetName '\'];
if ~exist(trainingSetDir, 'dir'); mkdir(trainingSetDir); end
load([getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetName '\trainingData.mat'], 'trainingData')
prepareTrainingImages(trainingSetDir, trainingData, features);


