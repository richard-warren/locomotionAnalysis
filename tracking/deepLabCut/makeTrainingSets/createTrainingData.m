%% initialize new training data structure

% settings
trainingSetFolder = 'barObstacle';
view = 'both';
sessions = {'180605_000', '180605_001', '180605_002', ...
            '180609_000', '180609_001', '180609_002', '180609_003', '180609_004', '180609_005'};
frameNum = 200;
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

trainingSetFolder = 'barObstacle';
trainingSetName = 'trainingData.mat';

trainingSetDir = [getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetFolder '\'];
features = {'paw1LH_top', 'paw2LF_top', 'paw3RF_top', 'paw4RH_top', 'tailBase_top', 'tailMid_top', 'nose_top', 'obs_top', ...
            'paw1LH_bot', 'paw2LF_bot', 'paw3RF_bot', 'paw4RH_bot', 'tailBase_bot', 'tailMid_bot', 'nose_bot', 'obsHigh_bot', 'obsLow_bot'}; % with top view
labelFrames(trainingSetDir, trainingSetName, features);


%% add incorrect frames from tracked vid

session = '180623_006';
trainingSetFolder = 'barObstacle';
trainingSetName = 'trainingData.mat';

trainingDataPath = [getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetFolder '\' trainingSetName];
showTrackingDLC(session, .02, trainingDataPath)



%% prepare data for deepLabCut

trainingSetFolder = 'barObstacle';
scaling = 1;
trainingSetName = 'trainingData.mat';
% features = {'paw1', 'paw2', 'paw3', 'paw4', 'gen', 'tailBase', 'tailMid', 'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH', 'tailBaseTop', 'tailMidTop'};
features = {'paw1LH_top', 'paw2LF_top', 'paw3RF_top', 'paw4RH_top', 'tailBase_top', 'tailMid_top', 'nose_top', 'obs_top', ...
            'paw1LH_bot', 'paw2LF_bot', 'paw3RF_bot', 'paw4RH_bot', 'tailBase_bot', 'tailMid_bot', 'nose_bot', 'obsHigh_bot', 'obsLow_bot'}; % with top view

trainingDataPath = [getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetFolder '\' trainingSetName];
load(trainingDataPath, 'trainingData', 'view')
writeDir = [getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetFolder '\' trainingSetFolder 'Scaling' num2str(scaling) '\']; % images will be written here
if ~exist(writeDir, 'dir'); mkdir(writeDir); end
prepareTrainingImages(writeDir, trainingData, view, features, scaling);

% %% prepare data without obstacles in view for matthis paper
% 
% trainingSetFolder = 'barObstacle';
% xMax = 410; % a little larger than width of video
% targetSize = [554 554];
% trainingSetName = 'trainingData.mat';
% features = {'paw1LH_top', 'paw2LF_top', 'paw3RF_top', 'paw4RH_top', 'tailBase_top', 'tailMid_top', 'nose_top', ...
%             'paw1LH_bot', 'paw2LF_bot', 'paw3RF_bot', 'paw4RH_bot', 'tailBase_bot', 'tailMid_bot', 'nose_bot'}; % with top view
% trainingDataPath = [getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetFolder '\' trainingSetName];
% load(trainingDataPath, 'trainingData', 'view')
% 
% 
% % remove frames where obs is in view
% prevSes = '';
% obsLims = [.20 .45]; % !!! this is a hack, because obsPositions are unreliable and vary session to session
% validBins = false(1,length(trainingData));
% 
% for i = 1:length(trainingData)
%     
%     % if reached new session, load data
%     if ~strcmp(trainingData(i).session, prevSes)
%         sesData = load([getenv('OBSDATADIR') 'sessions\' trainingData(i).session '\runAnalyzed.mat'], ...
%             'obsTimes', 'obsPositions', 'frameTimeStamps');
%         prevSes = trainingData(i).session;
%     end
%     
%     % mark valid if obs is out of view
%     obsPos = sesData.obsPositions(find(sesData.obsTimes>sesData.frameTimeStamps(trainingData(i).frameNum),1,'first'));
%     if obsPos<obsLims(1) || obsPos>obsLims(2); validBins(i) = true; end
% 
% end
% sum(validBins)
% %
% writeDir = [getenv('OBSDATADIR') 'tracking\trainingData\deepLabCut\' trainingSetFolder '\' trainingSetFolder 'Scaling' num2str(scaling) 'NoObstacles\']; % images will be written here
% if ~exist(writeDir, 'dir'); mkdir(writeDir); end
% prepareTrainingImagesMathis(writeDir, trainingData(validBins), view, features, targetSize);








