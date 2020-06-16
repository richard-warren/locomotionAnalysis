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
maxFiltering = 1;

% wisk
% trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_wisk.mat';
% skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_wisk.csv';  % skeletons follow the 'deepposekit' format
% invert = true;
% maxFiltering = 3;


labelFrames(trainingSet, skeleton, 'vidScaling', 2, 'invertFrame', invert, 'maxFiltering', maxFiltering);

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
fprintf('\ncreated file: %s\n', fileName)
fprintf('frames:       %i\n', size(images,4));


%% run deepposekit\training.py to train the network!


%% add incorrect frames from tracked vid

session = '181115_000';

% run
vid = 'run.mp4';  % run_originalDimensions;
trackedFeatures = 'trackedFeaturesRaw.csv';
trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_run.mat';
skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv';  % skeletons follow the 'deepposekit' format
invert = false;

% wisk
% vid = 'runWisk.mp4';  % run_originalDimensions;
% trackedFeatures = 'trackedFeaturesRaw_wisk.csv';
% trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_wisk.mat';
% skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_wisk.csv';  % skeletons follow the 'deepposekit' format
% invert = true;

addToTrainingSet(session, vid, trackedFeatures, trainingSet, 'skeleton', skeleton, 'invertFrame', invert, 'scoreThresh', .5, 'zoom', []);

