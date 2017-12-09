function trainSVM(className, egPortion)
    
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
    
    % randomize files list
    files = files(randperm(length(files)));
    
    % initialize category feature matrix
    load([dataDir 'trainingImages\' className '\' category '\' files{1}], 'img')
    egNum = round(length(files) * egPortion);
    featuresAll{i} = nan(egNum, numel(img));

    % get training examples
    for j = 1:egNum
        
        % load image
        load([dataDir 'trainingImages\' className '\' category '\' files{j}], 'img')
                
        % extract and save features
        img = getFeatures(img);
        featuresAll{i}(j,:) = img(:);
    end
    
    % store category labels
    labels = vertcat(labels, ones(egNum,1)*(i));
end

% concatonate feature matrices
features = nan(0, numel(img));
for i = 1:length(categories)
    features = [features; featuresAll{i}];
end


% train classifer
% fprintf('training classifier...\n');
% tic
% model = svmtrain(labels, features, '-t 0 -s svm_type 2');
% fprintf('training time: %i minutes\n', toc\60);
% model.w = model.sv_coef' * model.SVs;

% keyboard

% tic
% rng(1); % for reproducibility initialize random seed
model = fitcsvm(features, labels, 'KernelScale', 559.69);
modelCrossVal = crossval(model);
fprintf('generalization loss: %f\n', kfoldLoss(modelCrossVal));
% model.w = sum((modelRaw.Alpha .* modelRaw.Y(modelRaw.IsSupportVector)) .* modelRaw.SupportVectors, 1); % 267x1681

% modelRaw = fitcsvm(features, labels, 'KernelFunction', 'linear', 'OptimizeHyperparameters',  {'BoxConstraint', 'KernelScale'})
% modelRaw = fitcsvm(features, labels, 'KernelFunction', 'polynomial', 'OptimizeHyperparameters',  {'BoxConstraint', 'KernelScale', 'PolynomialOrder'})
% modelRaw = fitcsvm(features, labels, 'KernelFunction', 'gaussian', 'OptimizeHyperparameters',  {'BoxConstraint', 'KernelScale'})
% toc




subHgt = size(img,1);
subWid = size(img,2);
uisave ({'model', 'subHgt', 'subWid'}, [dataDir 'classifiers\' className '.mat']);



