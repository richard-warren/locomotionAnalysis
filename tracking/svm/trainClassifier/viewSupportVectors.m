function viewSupportVectors
    load('trainClassifier\classifiers\classifier_11-Feb-2017.mat');
    
    sptVecs = classifier.SupportVectors;
    
    aspectRatio = 2/1;
    rows = ceil(sqrt((1/aspectRatio)*(size(sptVecs,1))));
    cols = ceil(size(sptVecs,1)/rows);
    figure('units', 'normalized', 'outerposition', [0 .1 1 .9])
    subaxis(rows, cols, 1, 'spacing', 0.05, 'margin', .05);
    
    
    
    for i = 1:size(sptVecs,1)
        subaxis(rows, cols, i);
        imshow(reshape(sptVecs(i,:), sqrt(size(sptVecs,2)/3), sqrt(size(sptVecs,2)/3),3));
    end

end

