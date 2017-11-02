function viewTrainingSet(view, type)

    % user settings
    dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\trainingImages\';
    
    % initializations
    fullDir = [dataDir view '\' type];
    files = dir(fullDir);
    files = {files.name};
    files = files(3:end);
    
    aspectRatio = 2/1;
    rows = ceil(sqrt((1/aspectRatio)*(length(files))));
    cols = ceil(length(files)/rows);
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'black', 'menubar', 'none')
    subaxis(rows, cols, 1, 'spacing', 0.02, 'margin', .02);
    
    for i = 1:length(files)
        load([fullDir '\' files{i}], 'img');
        subaxis(i)
        imshow(getFeatures(img));
        title(i, 'color', 'white')
    end
end

