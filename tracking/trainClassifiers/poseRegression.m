

% settings
trainTestValPortions = [.6 .2 .2];
imgDir = 'C:\Users\rick\Desktop\trainingExamples\poseRegression\';

% initializations
load([imgDir 'pawLocations.mat'], 'features', 'locations')

%%

net = alexnet; % load alexNet

% get alexNet conv layers, and add new fully connected layers
layersTransfer = net.Layers(1:16);
numOutputs = size(locations,2);
learningRateFactor = 10;
layers = [layersTransfer
          fullyConnectedLayer(4096, 'WeightLearnRateFactor', learningRateFactor, 'BiasLearnRateFactor', learningRateFactor)
          reluLayer
          dropoutLayer(.5)
          fullyConnectedLayer(4096, 'WeightLearnRateFactor', learningRateFactor, 'BiasLearnRateFactor', learningRateFactor)
          reluLayer
          dropoutLayer(.5)
          fullyConnectedLayer(numOutputs, 'WeightLearnRateFactor', learningRateFactor, 'BiasLearnRateFactor', learningRateFactor)
          regressionLayer];


% set training parameters
miniBatchSize = 32;
numIterationsPerEpoch = floor(size(locations,1)/miniBatchSize);
options = trainingOptions('sgdm',...
    'MiniBatchSize', miniBatchSize,...
    'MaxEpochs', 4,...
    'InitialLearnRate', 1e-4,...
    'Verbose', true,...
    'Plots', 'training-progress');

% train!
convNetwork = trainNetwork(features, layers, options);
save([getenv('OBSDATADIR') 'tracking\classifiers\botPoseRegressor.mat'], 'convNetwork')

% classify
predictedLabels = classify(convNetwork, testImages);
fprintf('test accuracy: %f\n', mean(predictedLabels == testImages.Labels));


