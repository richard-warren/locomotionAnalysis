

% settings
trainPortion = .9;
imgDir = 'C:\Users\rick\Desktop\trainingExamples\poseRegression\';

% initializations
originalImSize = [230 396];
load([imgDir 'pawLocations.mat'], 'features')
m = height(features);
allInds = randperm(m);
trainData = features(allInds(1:floor(m*trainPortion)),:);
valData = features(allInds(floor(m*trainPortion)+1:end),:);

%%

% settings
learningRateFactor = 5;
miniBatchSize = 32;

net = alexnet; % load alexNet

% get alexNet conv layers, and add new fully connected layers
layersTransfer = net.Layers(1:16);
numOutputs = size(trainData,2)-1;

layers = [layersTransfer
          fullyConnectedLayer(4096, 'WeightLearnRateFactor', learningRateFactor, 'BiasLearnRateFactor', learningRateFactor)
          reluLayer
          dropoutLayer(.5)
          fullyConnectedLayer(4096, 'WeightLearnRateFactor', learningRateFactor, 'BiasLearnRateFactor', learningRateFactor)
          reluLayer
          dropoutLayer(.5)
          fullyConnectedLayer(numOutputs, 'WeightLearnRateFactor', learningRateFactor, 'BiasLearnRateFactor', learningRateFactor)
          regressionLayer];
% layersTransfer = net.Layers(1:end-3);
% numOutputs = size(trainData,2)-1;
% learningRateFactor = 20;
% layers = [layersTransfer
%           fullyConnectedLayer(numOutputs, 'WeightLearnRateFactor', learningRateFactor, 'BiasLearnRateFactor', learningRateFactor)
%           regressionLayer];


% set training parameters
numIterationsPerEpoch = floor(size(trainData,1)/miniBatchSize);
options = trainingOptions('sgdm',...
    'MiniBatchSize', miniBatchSize,...
    'MaxEpochs', 20,...
    'InitialLearnRate', 1e-5,... % was originally 1e-4
    'Verbose', true,...
    'VerboseFrequency', 20,...
    'Plots', 'training-progress', ...
    'ValidationData', valData,...
    'ValidationFrequency', numIterationsPerEpoch);

% train!
convNetwork = trainNetwork(trainData, layers, options);
save([getenv('OBSDATADIR') 'tracking\classifiers\botPoseRegressor.mat'], 'convNetwork')

% classify
% predictedLabels = predict(convNetwork, testImages);
% fprintf('test accuracy: %f\n', mean(predictedLabels == testImages.Labels));

%% test on image

for imNum = randperm(height(features), 10)

    % get image and make predictions!
    img = imread(['C:\Users\rick\Desktop\trainingExamples\poseRegression\imgs\img' num2str(imNum) '.tif']);
    predictedLocations = predict(convNetwork, img);
%     imgResized = imresize(img, 'outputSize', originalImSize);

    % show results
    close all; figure('position', [1286 66 570 323]);
    imshow(img);

    hold on; scatter(predictedLocations([1 3 5 7])*size(img,2), predictedLocations([2 4 6 8])*size(img,1))
    hold on; scatter(table2array(features(imNum, [1 3 5 7]+1))*size(img,2), table2array(features(imNum, [2 4 6 8]+1))*size(img,1))
    pause(1)
end











