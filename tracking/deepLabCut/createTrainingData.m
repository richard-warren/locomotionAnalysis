

% settings
trainingSetName = 'newTrainingData';
sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003', ...
            '180125_001', '180125_002', '180125_003'};
frameNum = 1000;
obsPortion = .5; % portion of trials to include obstacle
features = {'pawTL', 'pawTR', 'pawBR', 'pawBL', 'gen', 'tailBase', 'tailMid', 'tailEnd'};

% initializations
trainingSetDir = [getenv('TRAININGEXAMPLESDIR') 'deepLabCut\' trainingSetName '\'];
if ~exist(trainingSetDir, 'dir'); mkdir(trainingSetDir); end


%% create struct if doesn't already exist
if ~exist([trainingSetDir 'trainingData.mat'], 'file')
    trainingData = createTrainingDataStruct(sessions, frameNum, obsPortion);
    save([trainingSetDir 'trainingData.mat'], 'trainingData')
end

%% label things

trainingSetName = 'newTrainingData';
features = {'pawTL', 'pawTR', 'pawBR', 'pawBL', 'gen', 'tailBase', 'tailMid', 'tailEnd'};
labelFrames(trainingSetName, features);



%% create images