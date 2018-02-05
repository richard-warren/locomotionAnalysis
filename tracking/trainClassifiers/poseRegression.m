

% settings
trainPortion = .9;
imgDir = [getenv('TRAININGEXAMPLESDIR') 'poseRegression\lowRes\'];

% initializations
originalImSize = [230 396];
load([imgDir 'pawLocations.mat'], 'features')

% change file directory name if necessary
trainingDir = getenv('TRAININGEXAMPLESDIR');
sampleFileName = cell2mat(features.imgNames(1));

if ~strcmp(trainingDir, sampleFileName(1:length(trainingDir)))
    for i = 1:height(features)
        fileName = cell2mat(features.imgNames(i));
        fileName = [trainingDir fileName(40:end)]; % the 46 in this line is a super hacky temp
        features.imgNames(i) = {fileName};
    end
end

m = height(features);
allInds = randperm(m);
trainData = features(allInds(1:floor(m*trainPortion)),:);
valData = features(allInds(floor(m*trainPortion)+1:end),:);

%% retrain alexnet

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

%% train network from scratch

% settings
miniBatchSize = 64;
numOutputs = size(trainData,2)-1;

layers = [imageInputLayer([59 99 1])
    
    convolution2dLayer(3,32)
    reluLayer
    maxPooling2dLayer(2)
    
    convolution2dLayer(2,64)
    reluLayer
    maxPooling2dLayer(2)
    
    convolution2dLayer(2,128)
    reluLayer
    maxPooling2dLayer(2)
    
    fullyConnectedLayer(500)
    reluLayer
%     dropoutLayer(.5)
    
    fullyConnectedLayer(500)
    reluLayer
%     dropoutLayer(.5)
    
    fullyConnectedLayer(numOutputs)
    regressionLayer];

% set training parameters
numIterationsPerEpoch = floor(size(trainData,1)/miniBatchSize);
options = trainingOptions('sgdm',...
    'MiniBatchSize', miniBatchSize,...
    'MaxEpochs', 30,...
    'InitialLearnRate', 1e-2,... % was originally 1e-4
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.1, ...
    'LearnRateDropPeriod', 1, ... % epoches per drop
    'Verbose', true,...
    'VerboseFrequency', 50,...
    'Plots', 'training-progress', ...
    'ValidationData', valData,...
    'ValidationFrequency', numIterationsPerEpoch);

% train!
convNetwork = trainNetwork(trainData, layers, options);
% save([getenv('OBSDATADIR') 'tracking\classifiers\botPoseRegressorCustom.mat'], 'convNetwork')


%% test on image

% imNum = 1;
close all; figure();
img = imread([getenv('TRAININGEXAMPLESDIR') 'poseRegression\lowRes\imgs\img1.tif']);
preview = imshow(img);
hold on; scatterTruth = scatter(gca, [0 0 0 0], [0 0 0 0], 100, hsv(4), 'filled');
hold on; scatterPredict = scatter(gca, [0 0 0 0], [0 0 0 0], 200, hsv(4));
set(gcf, 'position', [646   173   984   750]);

% pimpFig;


for imNum = randperm(height(features), 10)
    disp(imNum)

    % get image and make predictions!
    img = imread([getenv('TRAININGEXAMPLESDIR') 'poseRegression\lowRes\imgs\img' num2str(imNum) '.tif']);
    predictedLocations = predict(convNetwork, img);
%     imgResized = imresize(img, 'outputSize', originalImSize);

    % show results
    
    set(preview, 'CData', img)
    set(scatterTruth, 'XData', table2array(features(imNum, [1 3 5 7]+1))*size(img,2), ...
        'YData', table2array(features(imNum, [2 4 6 8]+1))*size(img,1));
    set(scatterPredict, 'XData', predictedLocations([1 3 5 7])*size(img,2), ...
        'YData', predictedLocations([2 4 6 8])*size(img,1));
    
    pause(3)
end











