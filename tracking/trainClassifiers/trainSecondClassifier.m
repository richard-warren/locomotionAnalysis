function trainSecondClassifier(className)

% !!! need to document


% user settings
dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\';


% initializations
load([dataDir '\trainingData\' className '\labeledFeatures.mat'], 'features', 'labels', 'subFrameSize')


% train second classifier
model = fitcsvm(features', labels, 'KernelFunction', 'Gaussian', 'KernelScale', sqrt(size(features,1)), 'Standardize', true);
% model = fitcsvm(features', labels, 'OptimizeHyperparameters', {'BoxConstraint', 'KernelScale'}, 'Standardize', true);
modelCrossVal = crossval(model);
fprintf('generalization loss: %f\n', kfoldLoss(modelCrossVal));


% save model
uisave ({'model', 'subFrameSize'}, [dataDir 'classifiers\' className '2.mat']);



