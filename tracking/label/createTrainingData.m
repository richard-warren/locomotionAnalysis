%% initialize new training data structure

% settings
% view = 'run';  % 'run' or 'wisk'
% vidName = 'run_originalDimensions.mp4';
% trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_run.mat';

view = 'wisk';  % 'run' or 'wisk'
vidName = 'runWisk.mp4';
trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_wisk.mat';

frameNum = 40;  % total number of frames

% read sessions from spreadsheet
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'sheet', 'trainingSetSessions');
sessions = sessionInfo.session(sessionInfo.include==1);

if ~exist(trainingSet, 'file')
    trainingData = createTrainingDataStruct(sessions, view, vidName, frameNum, 'method', 'metadata');
    save(trainingSet, 'trainingData')
else
    fprintf('%s already exists... did not create file\n', trainingSet);
end

%% label dataset

% run
trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_run.mat';
skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv';  % skeletons follow the 'deepposekit' format
invert = false;

% wisk
% trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_wisk.mat';
% skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_wisk.csv';  % skeletons follow the 'deepposekit' format
% invert = true;


labelFrames(trainingSet, skeleton, 'vidScaling', 2, 'invertFrame', invert);

%% convert dataset to deepposekit format

% note: training data must be in sizes that work for deeppostkit // later
% on, when evaluating videos, the VideoReader can resize things as
% necessary on the fly

% settings
% trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_wisk.mat';
% skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_wisk.csv';  % skeletons follow the 'deepposekit' format
% imgDims = [352 384];  % images are resized to this dimension ([336 380] are the original dimensions) (must be valid DPK dimensions)

trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_run.mat';
skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv';  % skeletons follow the 'deepposekit' format
imgDims = [448 448];  % images are resized to this dimension (must be valid DPK dimensions)


% load data
load(trainingSet, 'trainingData');
trainingData = trainingData([trainingData.includeFrame]);  % restrict to frames to be included
skeletonTbl = readtable(skeleton);
features = skeletonTbl.name;
numEgs = length(trainingData);
numFeatures = length(features);

% create h5 file (for deepposekit) with datasets: annotated, annotations, images, skeleton

% note: matlab h5 files for some reason are read by python with the
% dimension order flipped, so i flip them here to make them work in python

annotated = true(numEgs, numFeatures);
annotations = nan(numEgs, numFeatures, 2);
for i = 1:numEgs; for j = 1:numFeatures; annotations(i,j,:) = trainingData(i).(features{j}); end; end


% resize images
for i = 1:length(trainingData)
    sz = size(trainingData(i).frame);  % original size of training image
    frame = zeros(imgDims);
    x = min(sz(2), imgDims(2));
    y = min(sz(1), imgDims(1));
    frame(1:y,1:x) = trainingData(i).frame(1:y,1:x);
    trainingData(i).frame = frame;
end
images = permute(cat(4, trainingData.frame), [3,2,1,4]);

skel = -ones(numFeatures, 2);  % first column in index (zero-based for python) of parent features; second column is zero-based index of 'swap' feature
for i = 1:numFeatures
    if ~isempty(skeletonTbl.parent{i}); skel(i,1) = find(strcmp(skeletonTbl.name, skeletonTbl.parent{i})) - 1; end  % get index for parent feature
    if ~isempty(skeletonTbl.swap{i}); skel(i,2) = find(strcmp(skeletonTbl.name, skeletonTbl.swap{i})) - 1; end  % get index for swap feature
end

% create h5 file
[filepath, name] = fileparts(trainingSet);
fileName = fullfile(filepath, [name '.h5']);
if exist(fileName, 'file'); delete(fileName); end

h5create(fileName, '/annotated', fliplr(size(annotated)), 'Datatype', 'uint8')  % matlab doesn't support boolean datatypes, what the fuck
h5create(fileName, '/annotations', fliplr(size(annotations)), 'Datatype', 'double')  % should really be <f8, whatever that means...
h5create(fileName, '/images', size(images), 'Datatype', 'uint8')
h5create(fileName, '/skeleton', fliplr(size(skel)), 'Datatype', 'int32')

h5write(fileName, '/annotated', permute(uint8(annotated), [2 1]))
h5write(fileName, '/annotations', permute(uint16(annotations), [3 2 1]))
h5write(fileName, '/images', images)
h5write(fileName, '/skeleton', permute(int32(skel), [2 1]))

h5disp(fileName) 
fprintf('created file: %s\n', fileName)


%% run deepposekit\training.py to train the network!


%% add incorrect frames from tracked vid

session = '999999_999';

% run
vid = 'run.mp4';  % run_originalDimensions;
trackedFeatures = 'trackedFeaturesRaw.csv';
trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_run.mat';
skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv';  % skeletons follow the 'deepposekit' format
invert = false;

% wisk
% vid = 'runWisk.mp4';  % run_originalDimensions;
% trackedFeatures = 'trackedFeatures_wisk.csv';
% trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_wisk.mat';
% skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_wisk.csv';  % skeletons follow the 'deepposekit' format
% invert = true;

addToTrainingSet(session, vid, trackedFeatures, trainingSet, 'skeleton', skeleton, 'invertFrame', invert, 'scoreThresh', 0);

%% analyze batch of videos

% read sessions from spreadsheet
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionsToLabel.xlsx'), 'sheet', 'runSessions');
groomingSessions = {'200113_000', '200116_000', '200117_000', '200114_000', '200131_000', '200202_000', '200130_000'};
otherSessions = sessionInfo.session(sessionInfo.include==1);
otherSessions = otherSessions(~ismember(otherSessions, groomingSessions));

% choose sessions to use
% sessions = otherSessions;
% sessions = groomingSessions;
sessions = {groomingSessions{:}, otherSessions{:}};

for i = 20:length(sessions)
    fprintf('%i/%i ', i, length(sessions))
    try
        dpkAnalysis(sessions{i}, 'verbose', false, ...
            'vid', 'run_originalDimensions.mp4', ...
            'model', 'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_run_StackedDenseNet.h5', ...
            'skeleton', 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv', ...
            'output', 'trackedFeatures_run.csv')
        dpkAnalysis(sessions{i}, 'verbose', false, ...
            'vid', 'run_originalDimensions.mp4', ...
            'model', 'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_run_DeepLabCut.h5', ...
            'skeleton', 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv', ...
            'output', 'trackedFeatures_runDLC.csv')
    catch
        fprintf('%s: problem with analysis!\n', sessions{i})
    end
end
disp('all done!')


% %% (temp) fix images in training set
% 
% % run
% % trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_run.mat';
% % vidName = 'run_originalDimensions.mp4';  % run_originalDimensions;
% 
% % wisk
% trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_wisk.mat';
% vidName = 'runWisk.mp4';  % run_originalDimensions;
% 
% 
% load(trainingSet)
% currentSession = '';
% for i = 1:length(trainingData)
%     
%     % load new session if necessary
%     if ~strcmp(trainingData(i).session, currentSession)
%         fprintf('%s: loading %s...\n', trainingData(i).session, vidName)
%         vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', trainingData(i).session, vidName));
%         currentSession = trainingData(i).session;
%     end
%     
%     % update frame
%     trainingData(i).frame = rgb2gray(read(vid, trainingData(i).frameNum));
% end
% 
% save(trainingSet, 'trainingData')
% disp('all done!')

%% (temp) find out which sessions have 'originalDimensions' versions

files = dir(fullfile(getenv('OBSDATADIR'), 'sessions'));
sessions = {files([files.isdir]).name};
origSessions = {};

fileExists = false(1,length(sessions));
fprintf('\n\n--------------------looking for originalDimensions--------------------\n')
for i = 1:length(sessions)
    dirSub = dir(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, '*.mp4'));
    bins = contains({dirSub.name}, 'originalDimensions');
    if any(bins)
        fprintf('%s: ', sessions{i})
        fprintf('%s ', dirSub(bins).name)
        fprintf('\n')
        origSessions{end+1} = sessions{i};
    end
end
disp('all done!')

%% (temp) test models on old frames with different dimensions

dpkAnalysis('190929_003', 'run', 'verbose', true, 'runModel', 'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_run_StackedDenseNet.h5', 'runOutput', 'trackedFeatures_run.csv', 'runVid', 'run.mp4')
% dpkAnalysis('190929_003', 'wisk', 'verbose', true, 'wiskModel', 'D:\github\locomotionAnalysis\tracking\deepposekit\models\model_wisk.h5', 'wiskOutput', 'trackedFeatures_wisk.csv')




