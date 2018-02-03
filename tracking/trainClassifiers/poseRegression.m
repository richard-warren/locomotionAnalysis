

% settings
trainTestValPortions = [.6 .2 .2];
imgDir = 'C:\Users\rick\Desktop\trainingExamples\poseRegression\';


% initializations
imgs = imageDatastore([imgDir 'imgs'],...
    'IncludeSubfolders', true, 'FileExtensions', '.tif');
[trainImages, testImages, valImages] = splitEachLabel(imgs, ...
    trainTestValPortions(1), trainTestValPortions(2), trainTestValPortions(3), 'randomized');
load([imgDir 'pawLocations.mat'], 'locations')
%%

net = alexnet; % load alexNet

% get alexNet conv layers, and add new fully connected layers
layersTransfer = net.Layers(1:end-3);
numClasses = numel(categories(trainImages.Labels));
layers = [layersTransfer
          fullyConnectedLayer(numClasses, 'WeightLearnRateFactor', 20, 'BiasLearnRateFactor', 20)
          regressionLayer];

% set training parameters
miniBatchSize = 10;
numIterationsPerEpoch = floor(numel(trainImages.Labels)/miniBatchSize);
options = trainingOptions('sgdm',...
    'MiniBatchSize', miniBatchSize,...
    'MaxEpochs', 4,...
    'InitialLearnRate', 1e-4,...
    'Verbose', false,...
    'Plots', 'training-progress',...
    'ValidationData', valImages,...
    'ValidationFrequency', numIterationsPerEpoch);

% train!
convNetwork = trainNetwork(trainImages, locations, layers, options);
save([getenv('OBSDATADIR') 'tracking\classifiers\' class 'PoseRegressor.mat'], 'convNetwork', 'subFrameSize')

% classify
predictedLabels = classify(convNetwork, testImages);
fprintf('test accuracy: %f\n', mean(predictedLabels == testImages.Labels));


