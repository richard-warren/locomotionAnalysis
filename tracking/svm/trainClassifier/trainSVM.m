function trainSVM(className)
    
% user settings
dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\';
categories = {'negative', 'positive'};

% initializations
labels = [];


% iterate through image categories
for i = 1:length(categories)

    % get list of files in category folder
    category = cell2mat(categories(i));
    files = dir([dataDir 'trainingImages\' className '\' category]);
    files = {files.name};
    files = files(3:end);

    % get training examples
    for j = 1:length(files)

        load([dataDir 'trainingImages\' className '\' category '\' files{j}], 'imgTemp')

        % initialize storage variable on first pass
        if ~exist('features', 'var')
            features = nan(0, numel(imgTemp));
        end
        
        % extract and save features
        img = getFeatures(imgTemp);
        features(end+1,:) = img(:);

    end
    
    % store category labels
    labels = vertcat(labels, ones(length(files),1)*i);
end


% train classifer
model = svmtrain(labels, features, '-t 0');
model.w = model.sv_coef' * model.SVs;
subHgt = size(imgTemp,1);
subWid = size(imgTemp,2);
uisave ({'model', 'subHgt', 'subWid'}, [dataDir 'classifiers\' className '.mat']);



