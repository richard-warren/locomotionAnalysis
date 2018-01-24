function trainSVM2(className)

% !!! need to document, but generally trains a second stage classifier, SVM, that classifies subframes identified by the first linear SVM convolution as potential paws
% i am currently seeing whether this classifier can distinguish WHICH paw it is, which would be pretty dope


% initializations
dataDir = [getenv('OBSDATADIR') 'svm\'];
load([dataDir '\trainingData\' className '\labeledFeatures.mat'], 'features', 'labels', 'subFrameSize')


% train second SVM classifier (old, single class)
% model = fitcsvm(features', labels, 'KernelFunction', 'Gaussian', 'KernelScale', sqrt(size(features,1)), 'Standardize', true);
% modelCrossVal = crossval(model);
% fprintf('generalization loss: %f\n', kfoldLoss(modelCrossVal));

% train second SVM classifier (new, multiclass)
template = templateSVM(...
    'KernelFunction', 'polynomial', ...
    'PolynomialOrder', 2, ... % maybe should be 3
    'Standardize', true);
%     'KernelScale', 'auto', ...
%     'BoxConstraint', 1, ...
    

model = fitcecoc(...
    features', ... % flip s.t. observations are in rows
    labels, ...
    'Learners', template, ...
    'Coding', 'onevsone', ...
    'ClassNames', (1:length(unique(labels)))', ...
    'OptimizeHyperParameters', {'BoxConstraint','KernelScale'});
modelCrossVal = crossval(model);
fprintf('generalization loss: %f\n', kfoldLoss(modelCrossVal));

% save model
uisave ({'model', 'subFrameSize'}, [dataDir 'classifiers\' className]);



