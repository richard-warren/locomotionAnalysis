function viewTrainingSet(className)

% user settings
dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\trainingImages\';
figSize = [1000, 1800]; % h,w

% initializations
categories = {'positive', 'negative'};
colormap gray

for i = 1:length(categories)
    
    % get file names
    fullDir = [dataDir className '\' categories{i}];
    files = dir(fullDir);
    files = {files.name};
    files = files(3:end);
    
    % get img size
    load([fullDir '\' files{2}], 'img');
    imgSize = size(img);
    

    rows = floor(figSize(1) / imgSize(1));
    cols = floor(figSize(2) / imgSize(2));
    collage = nan(rows*imgSize(1), cols*imgSize(2));

    figure('units', 'pixels', 'position', [50 50 size(collage,2) size(collage,1)], 'color', 'black', 'menubar', 'none')
    
    for j = 1:rows
        for k = 1:cols
            
            % get random image
            fileInd = randi([1 length(files)]);
            load([fullDir '\' files{fileInd}], 'img');
            
            % incorporate image into collage
            rowInd =(j-1)*imgSize(1)+1;
            colInd =(k-1)*imgSize(2)+1; 
            collage(rowInd:rowInd+imgSize(1)-1, colInd:colInd+imgSize(2)-1) = img;
            
        end
    end
    
    % show collage
    colormap gray
    image(collage, 'CDataMapping', 'scaled');
    set(gca, 'units', 'normalized', 'position', [0 0 1 1], 'CLim', [0 255])
    
    
end

