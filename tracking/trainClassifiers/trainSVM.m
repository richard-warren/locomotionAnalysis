function trainSVM(className)

% !!! need to document


% user settings
dataDir = [getenv('OBSDATADIR') 'svm\'];


% initializations
load([dataDir '\trainingData\' className '\labeledFeatures.mat'], 'features', 'labels', 'subFrameSize')


% train svm
% model = fitcsvm(features', labels, 'KernelScale', 559.69);
model = fitcsvm(features', labels, 'KernelScale', 'Auto');
modelCrossVal = crossval(model);
fprintf('generalization loss: %f\n', kfoldLoss(modelCrossVal));


% save model
uisave ({'model', 'subFrameSize'}, [dataDir 'classifiers\' className '.mat']);



