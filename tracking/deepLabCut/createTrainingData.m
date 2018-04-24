%% initialize new training data structure

% settings
trainingSetName = 'topBotCat';
view = 'both';
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
    save([trainingSetDir 'trainingData.mat'], 'trainingData', 'view')
else
    fprintf('%s already exists... did not create file\n', trainingSetName);
end

%% label things

trainingSetName = 'topBotCat';

trainingSetDir = [getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetName '\'];
% features = {'pawTL', 'pawTR', 'pawBR', 'pawBL', 'gen', 'tailBase', 'tailMid', 'tailEnd'};
features = {'paw1', 'paw2', 'paw3', 'paw4', 'gen', 'tailBase', 'tailMid', 'tailEnd', 'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH'}; % with top view
labelFrames(trainingSetDir, features);



%% prepare data for deepLabCut
trainingSetName = 'topBotCat';
% features = {'pawTL', 'pawTR', 'pawBR', 'pawBL','gen', 'tailBase', 'tailMid', 'tailEnd'}; % excluding genitals
features = {'paw1', 'paw2', 'paw3', 'paw4', 'gen', 'tailBase', 'tailMid', 'tailEnd', 'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH'}; % with top view

trainingSetDir = [getenv('TRAININGEXAMPLESDIR') 'deepLabCut\' trainingSetName '\'];
if ~exist(trainingSetDir, 'dir'); mkdir(trainingSetDir); end
load([getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetName '\trainingData.mat'], 'trainingData', 'view')
prepareTrainingImages(trainingSetDir, trainingData, view, features);


