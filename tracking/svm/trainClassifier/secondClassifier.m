function secondClassifier(className)

% !!! need to document


% user settings
dataDir = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\svm\';


% initializations
load([dataDir '\trainingData\' className '\labeledFeatures.mat'], 'features', 'labels', 'subFrameSize')


% train svm
model = fitcsvm(features', labels, 'KernelFunction', 'Gaussian', 'OptimizeHyperparameters', {'KernelScale', 'BoxConstraint'});
% modelCrossVal = crossval(model);
% fprintf('generalization loss: %f\n', kfoldLoss(modelCrossVal));


% save model
uisave ({'model', 'subFrameSize'}, [dataDir 'classifiers\' className '.mat']);



