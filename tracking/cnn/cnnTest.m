


% settings
class = 'pawTop2';
targetSize = [227 227];
trainTestValPortions = [.6 .2 .2];

% initializations
imgs = imageDatastore([getenv('OBSDATADIR') 'svm\trainingData\' class],...
    'IncludeSubfolders', true, 'FileExtensions', '.tif', 'LabelSource', 'foldernames');
[trainingImages, testImages, validationImages] = splitEachLabel(imgs, ...
    trainTestValPortions(1), trainTestValPortions(2), trainTestValPortions(3), 'randomized');

%% retrain alexnet on my data



% load alexNet
net = alexnet;


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
    'ValidationData', validationImages,...
    'ValidationFrequency', numIterationsPerEpoch);


% train!
netTransfer = trainNetwork(trainingImages, layers, options);
save([getenv('OBSDATADIR') 'svm\classifiers\' class 'Cnn.mat'], 'netTransfer')


% classify
predictedLabels = classify(netTransfer, testImages);
fprintf('test accuracy: %f\n', mean(predictedLabels == testImages.Labels));







