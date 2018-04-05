function viewTrainingSet(className)


% settings
% figSize = [1000, 1800]; % h,w
figSize = [300, 300]*2; % h,w

% initializations
load([getenv('OBSDATADIR') 'tracking\trainingData\' className '\labeledFeatures.mat'], 'images', 'labels', 'subFrameSize')


for i = unique(labels)
    
    % determine figure dimensions
    rows = floor(figSize(1) / subFrameSize(1));
    cols = floor(figSize(2) / subFrameSize(2));
    collage = nan(rows*subFrameSize(1), cols*subFrameSize(2));

    figure('units', 'pixels', 'position', [50 50 size(collage,2) size(collage,1)], 'color', 'black', 'menubar', 'none')
    
    inds = find(labels==i);
    inds = inds(randperm(length(inds)));
    currentInd = 1;
    
    for j = 1:rows
        for k = 1:cols
            if currentInd<=length(inds)
                
                % get random image
                img = images(1:prod(subFrameSize), inds(currentInd));
                img = reshape(img, subFrameSize(1), subFrameSize(2));

                % incorporate image into collage
                rowInd =(j-1) * subFrameSize(1) + 1;
                colInd =(k-1) * subFrameSize(2) + 1; 
                collage(rowInd:rowInd+subFrameSize(1)-1, colInd:colInd+subFrameSize(2)-1) = img;

                currentInd = currentInd + 1;
            end
        end
    end
    
    % show collage
    colormap gray;
    image(gca, collage, 'CDataMapping', 'scaled');
    set(gca, 'units', 'normalized', 'position', [0 0 1 1], 'CLim', [0 255])
end

