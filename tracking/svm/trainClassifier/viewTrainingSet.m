function viewTrainingSet
    folder = uigetdir('trainClassifier\trainingImages');

    temp = dir(folder); files = {temp.name}; files = files(3:end);
    
    aspectRatio = 2/1;
    rows = ceil(sqrt((1/aspectRatio)*(length(files))));
    cols = ceil(length(files)/rows);
    
    figure('units', 'normalized', 'outerposition', [0 0 1 1])
    subaxis(rows, cols, 1, 'spacing', 0.01, 'margin', .01);
    
    for i = 1:length(files)
        load([folder '\' files{i}], 'imgTemp');
        subaxis(i)
        imshow(imgTemp);
        title(i)
    end
end

