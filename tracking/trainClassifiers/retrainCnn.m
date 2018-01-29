function retrainCnn(class, network)


% settings
trainTestValPortions = [.6 .2 .2];



% initializations
imgs = imageDatastore([getenv('OBSDATADIR') 'tracking\trainingData\' class],...
    'IncludeSubfolders', true, 'FileExtensions', '.tif', 'LabelSource', 'foldernames');
[trainImages, testImages, valImages] = splitEachLabel(imgs, ...
    trainTestValPortions(1), trainTestValPortions(2), trainTestValPortions(3), 'randomized');
load([getenv('OBSDATADIR') 'tracking\trainingData\' class '\labeledFeatures.mat'], 'subFrameSize')


switch network
    
    case 'alexNet'
        
        net = alexnet; % load alexNet

        % get alexNet conv layers, and add new fully connected layers
        layersTransfer = net.Layers(1:end-3);
        numClasses = numel(categories(trainImages.Labels));
        layers = [layersTransfer
                  fullyConnectedLayer(numClasses, 'WeightLearnRateFactor', 20, 'BiasLearnRateFactor', 20)
                  softmaxLayer
                  classificationLayer];

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
        convNetwork = trainNetwork(trainImages, layers, options);
        save([getenv('OBSDATADIR') 'tracking\classifiers\' class 'AlexNet.mat'], 'convNetwork', 'subFrameSize')

        % classify
        predictedLabels = classify(convNetwork, testImages);
        fprintf('test accuracy: %f\n', mean(predictedLabels == testImages.Labels));
        
    
    case 'googleNet'
        
        net = googlenet; % load googlenet
        
        lgraph = layerGraph(net);
        lgraph = removeLayers(lgraph, {'loss3-classifier','prob','output'});

        numClasses = numel(categories(trainImages.Labels));
        newLayers = [
            fullyConnectedLayer(numClasses,'Name','fc','WeightLearnRateFactor',20,'BiasLearnRateFactor', 20)
            softmaxLayer('Name','softmax')
            classificationLayer('Name','classoutput')];
        lgraph = addLayers(lgraph,newLayers);
        
        lgraph = connectLayers(lgraph,'pool5-drop_7x7_s1','fc');
        
        options = trainingOptions('sgdm',...
            'MiniBatchSize',5,...
            'MaxEpochs',3,...
            'InitialLearnRate',1e-4,...
            'Plots', 'training-progress',...
            'Verbose', false, ...
            'ValidationData',valImages,...
            'ValidationFrequency',3);
        
        convNetwork = trainNetwork(trainImages,lgraph,options);
        save([getenv('OBSDATADIR') 'tracking\classifiers\' class 'GoogleNet.mat'], 'convNetwork', 'subFrameSize')
        
        predictedLabels = classify(convNetwork, testImages);
        fprintf('test accuracy: %f\n', mean(predictedLabels == testImages.Labels));


end




