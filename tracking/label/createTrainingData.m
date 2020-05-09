%% initialize new training data structure

% settings
view = 'run';  % 'run' or 'wisk'
vidName = 'run_originalDimensions.mp4';
trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_run.mat';

% view = 'wisk';  % 'run' or 'wisk'
% vidName = 'runWisk.mp4';
% trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_wisk.mat';

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

% trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_wisk.mat';
% skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_wisk.csv';  % skeletons follow the 'deepposekit' format

trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_run.mat';
skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv';  % skeletons follow the 'deepposekit' format

labelFrames(trainingSet, skeleton, 'vidScaling', 2, 'invertFrame', false);

%% convert dataset to deepposekit format

% settings
% trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_wisk.mat';
% skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_wisk.csv';  % skeletons follow the 'deepposekit' format
% padDims = [];

trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_run.mat';
skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv';  % skeletons follow the 'deepposekit' format
padDims = [448 448];  % leave empty if you want to automatically determine the padding dimensions

padImgs = true;  % dpk only works with dimensions continuosly divisible by 2 // if true, pad images to the next highest number that is continuously divisible by 2


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

images = permute(cat(4, trainingData.frame), [4,1,2,3]);
if padImgs
    
    % make list of dimensions that are valid for use with deepposekit (only numbers that are repeatedly divisible by 2)
    validDims = 2.^(1:12)' * [1 3 5 7 11];
    validDims = unique(validDims(:));
    if isempty(padDims)
        hgt = validDims(find(validDims>=size(images,2),1,'first'));
        wid = validDims(find(validDims>=size(images,3),1,'first'));
    else
        hgt = padDims(1);
        wid = padDims(2);
    end
    images_temp = zeros(size(images,1), hgt, wid, size(images,4));
    images_temp(:,1:min(size(images,2),hgt),1:min(size(images,3),wid),:) = images(:,1:min(hgt,end),1:min(wid,end),:);
    images = images_temp;
    fprintf('padding images to size (h,w): %i, %i\n', hgt, wid)
    clear images_temp;
end

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
h5create(fileName, '/images', fliplr(size(images)), 'Datatype', 'uint8')
h5create(fileName, '/skeleton', fliplr(size(skel)), 'Datatype', 'int32')

h5write(fileName, '/annotated', permute(uint8(annotated), [2 1]))
h5write(fileName, '/annotations', permute(uint16(annotations), [3 2 1]))
h5write(fileName, '/images', permute(images, [4 3 2 1]))
h5write(fileName, '/skeleton', permute(int32(skel), [2 1]))

h5disp(fileName) 
fprintf('created file: %s\n', fileName)


%% run deepposekit\training.py to train the network!


%% add incorrect frames from tracked vid

session = '191113_002';
vid = 'runWisk.mp4';
trackedFeatures = 'trackedFeatures_wisk.csv';

trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_wisk.mat';
skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_wisk.csv';  % skeletons follow the 'deepposekit' format

% trainingSet = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_run.mat';
% skeleton = 'D:\github\locomotionAnalysis\tracking\label\training_sets\skeleton_run.csv';  % skeletons follow the 'deepposekit' format

addToTrainingSet(session, vid, trackedFeatures, trainingSet, 'skeleton', skeleton);

%% analyze batch of videos

% read sessions from spreadsheet
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'sheet', 'trainingSetSessions');
sessions = sessionInfo.session(sessionInfo.include==1);
sessions = sessions(1:2);  % !!! temp

%%
for i = 1:length(sessions)
    try
        dpkAnalysis(sessions{i}, 'run', 'verbose', false)
        dpkAnalysis(sessions{i}, 'wisk', 'verbose', false)
    catch
        pfrintf('%s: problem with analysis!\n', sessions{i})
    end
end
disp('all done!')


