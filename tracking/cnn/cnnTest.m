


% settings
trainingData = 'pawBot';
targetSize = [227 227];
trainingPortion = .7;

% initializations
% prepareTrainingData(trainingData, targetSize); % restructure features
imgs = imageDatastore([getenv('OBSDATADIR') 'svm\trainingData\' trainingData],...
    'IncludeSubfolders', true, 'FileExtensions', '.tif', 'LabelSource', 'foldernames');
[trainingImages, testImages] = splitEachLabel(imgs, trainingPortion, 'randomized');

%% load alexNet
net = alexnet;

%% retrain alexNet on my data, omg


% get alexNet conv layers, and add new fully connected layers
layersTransfer = net.Layers(1:end-3);
numClasses = numel(categories(trainingImages.Labels));
layers = [layersTransfer
          fullyConnectedLayer(numClasses, 'WeightLearnRateFactor', 20, 'BiasLearnRateFactor', 20)
          softmaxLayer
          classificationLayer];

      
% set training parameters
miniBatchSize = 10;
numIterationsPerEpoch = floor(numel(trainingImages.Labels)/miniBatchSize);
options = trainingOptions('sgdm',...
    'MiniBatchSize', miniBatchSize,...
    'MaxEpochs', 4,...
    'InitialLearnRate', 1e-4,...
    'Verbose', false,...
    'Plots', 'training-progress',...
    'ValidationData', testImages,...
    'ValidationFrequency', numIterationsPerEpoch);

% train!
netTransfer = trainNetwork(trainingImages, layers, options);


%% classify

predictedLabels = classify(netTransfer, testImages);
fprintf('test accuracy: %f\n', mean(predictedLabels == testImages.Labels));


%% extract features

layer = 'fc7';
extractedFeatures = activations(net, featuresN, layer);







