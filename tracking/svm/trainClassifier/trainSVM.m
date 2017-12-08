function trainSVM(className)
    
% user settings
dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\';
categories = {'positive', 'negative'};

% initializations
labels = [];
featuresAll = cell(1,length(categories));

% iterate through image categories
for i = 1:length(categories)

    % get list of files in category folder
    category = cell2mat(categories(i));
    files = dir([dataDir 'trainingImages\' className '\' category]);
    files = {files.name};
    files = files(3:end);
    
    % initialize category feature matrix
    load([dataDir 'trainingImages\' className '\' category '\' files{1}], 'img')
    featuresAll{i} = nan(length(files), numel(img));

    % get training examples
    for j = 1:length(files)
        
        % load image
        load([dataDir 'trainingImages\' className '\' category '\' files{j}], 'img')
                
        % extract and save features
        img = getFeatures(img);
        featuresAll{i}(j,:) = img(:);
    end
    
    % store category labels
    labels = vertcat(labels, ones(length(files),1)*(i));
end

% concatonate feature matrices
features = nan(0, numel(img));
for i = 1:length(categories)
    features = [features; featuresAll{i}];
end


% train classifer
fprintf('training classifier...\n');
tic
model = svmtrain(labels, features, '-t 0 -s svm_type 2');
fprintf('training time: %i minutes\n', toc\60);
model.w = model.sv_coef' * model.SVs;

% keyboard

% tic
% rng(1); % for reproducibility initialize random seed
% modelRaw = fitcsvm(features, labels, 'KernelFunction', 'linear', 'OptimizeHyperparameters',  {'BoxConstraint', 'KernelScale'})
% modelRaw = fitcsvm(features, labels, 'KernelFunction', 'polynomial', 'OptimizeHyperparameters',  {'BoxConstraint', 'KernelScale', 'PolynomialOrder'})
% modelRaw = fitcsvm(features, labels, 'KernelFunction', 'gaussian', 'OptimizeHyperparameters',  {'BoxConstraint', 'KernelScale'})
% toc
% modelCrossVal = crossval(modelRaw);
% model.classLoss = kfoldLoss(modelCrossVal);
% fprintf('generalization loss: %f', model.classLoss);


subHgt = size(img,1);
subWid = size(img,2);
uisave ({'model', 'subHgt', 'subWid'}, [dataDir 'classifiers\' className '.mat']);



