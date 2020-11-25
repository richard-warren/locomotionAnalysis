%% load matlab training data struct

% fileName = 'Z:\loco\obstacleData\tracking\trainingData\deepLabCut\barObstacle\trainingData.mat';  % old dlc data
fileName = 'D:\github\locomotionAnalysis\tracking\label\training_sets\trainingset_run.mat';  % new dpk training set

load(fileName, 'trainingData');

%% copy all associated videos to folder

folder = 'D:\biderman_files\dpk';

% copy training data struct
save(fullfile(folder, 'trainingData.mat'), 'trainingData')

%% copy videos in dataset to folder
sessions = unique({trainingData.session});
for i = 1:length(sessions)
    fprintf('(%2i/%2i) copying session %s\n', i, length(sessions), sessions{i})
    source = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'run.mp4');
    if ~exist(source, 'file'); concatTopBotVids(sessions{i}); end
    destination = fullfile(folder, [sessions{i} '.mp4']);
    copyfile(source, destination);
end