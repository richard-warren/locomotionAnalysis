


% settings
trainingData = 'pawBot';
targetSize = [227 227];
trainingPortion = .7;

% initializations
% prepareTrainingData(trainingData); % restructure features
imgs = imageDatastore([getenv('OBSDATADIR') 'svm\trainingData\' trainingData],...
    'IncludeSubfolders', true, 'FileExtensions', '.png', 'LabelSource', 'foldernames');
featuresTrain = features(:,:,:,trainingBins);
featuresValidate = features(:,:,:,~trainingBins);

%% load alexNet
net = alexnet;

%% retrain alexNet on my data, omg


% get alexNet conv layers, and add new fully connected layers
layersTransfer = net.Layers(1:end-3);
numClasses = numel(length(unique(labels)));
layers = [layersTransfer
          fullyConnectedLayer(numClasses, 'WeightLearnRateFactor', 20, 'BiasLearnRateFactor', 20)
          softmaxLayer
          classificationLayer];

      
% set training parameters
miniBatchSize = 10;
numIterationsPerEpoch = floor(length(layers) / miniBatchSize);
options = trainingOptions('sgdm',...
    'MiniBatchSize', miniBatchSize,...
    'MaxEpochs', 4,...
    'InitialLearnRate', 1e-4,...
    'Verbose', false,...
    'Plots', 'training-progress',...
    'ValidationData', validationImages,...
    'ValidationFrequency', numIterationsPerEpoch);

% train!
netTransfer = trainNetwork(featuresN, labels, layers, options);


%% extract features

layer = 'fc7';
extractedFeatures = activations(net, featuresN, layer);







