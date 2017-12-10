function trainSVM(className)

% !!! need to document


% user settings
dataDir = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\svm\';


% initializations
load([dataDir '\trainingData\' className '\labeledFeatures.mat'], 'features', 'labels', 'subFrameSize')


% train svm
model = fitcsvm(features', labels, 'KernelScale', 559.69);
modelCrossVal = crossval(model);
fprintf('generalization loss: %f\n', kfoldLoss(modelCrossVal));
% model.w = sum((modelRaw.Alpha .* modelRaw.Y(modelRaw.IsSupportVector)) .* modelRaw.SupportVectors, 1); % 267x1681


% save model
uisave ({'model', 'subFrameSize'}, [dataDir 'classifiers\' className '.mat']);



