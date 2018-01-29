function trainSVM(className)

% !!! need to document



% initializations
load([getenv('OBSDATADIR') 'tracking\trainingData\' className '\labeledFeatures.mat'], 'features', 'labels', 'subFrameSize')


% train svm
% model = fitcsvm(features', labels, 'KernelScale', 559.69);
model = fitcsvm(features', labels, 'KernelScale', 'Auto');
modelCrossVal = crossval(model);
fprintf('generalization loss: %f\n', kfoldLoss(modelCrossVal));


% save model
uisave ({'model', 'subFrameSize'}, [getenv('OBSDATADIR') 'tracking\classifiers\' className '.mat']);



