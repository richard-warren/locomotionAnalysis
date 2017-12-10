function potentialLocationsBot = getPotentialLocationsBot(vid, model, features, labels, scoreThresh, subFrameSize, obsPixPositions, frameInds, showTracking)

% !!! need to document


% settings
overlapThresh = .5; % used for non-maxima suppression // higher numbers = more tightly packed
objectNum = 4;      % number of paws

% initializations
sampleFrame = rgb2gray(read(vid,1));
totalFrames = vid.NumberOfFrames;
kernel = reshape(model.Beta, subFrameSize(1), subFrameSize(2));
bg = getBgImage(vid, 1000, false);
centPad = floor(subFrameSize / 2); % y,x


% prepare figure
if showTracking
    
    figure; imagesc(-kernel);

    figure('position', [680 144 698 834], 'menubar', 'none', 'color', 'black'); colormap gray

    rawAxis = subaxis(2,1,1, 'spacing', 0, 'margin', 0);
    rawIm = image(sampleFrame, 'parent', rawAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off', 'CLim', [0 255]);
    hold on; scatterPts = scatter(rawAxis, 0, 0, 50, [1 1 1], 'filled');
    hold on; scatterKnn = scatter(rawAxis, 0, 0, 150, [1 0 0], 'linewidth', 3);

    predictAxis = subaxis(2,1,2, 'spacing', 0.01, 'margin', .01);
    predictIm = image(sampleFrame, 'parent', predictAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off', 'CLim', [0 10]);
end


potentialLocationsBot = struct();

for i = frameInds
    
    disp(i/totalFrames)
    
    % get frame and subframes
    frame = rgb2gray(read(vid,i));
    frame = getFeatures(frame);
    frame = frame - bg;
    
    % mask obstacle
    frame = maskObs(frame, obsPixPositions(i));

    % filter with svm
    frameFiltered = -(conv2(double(frame)/model.KernelParameters.Scale, kernel, 'same') + model.Bias);
    
    frameFiltered(frameFiltered < scoreThresh) = scoreThresh;
    frameFiltered = frameFiltered - scoreThresh;
    [x, y, scores] = nonMaximumSupress(frameFiltered, subFrameSize, overlapThresh);
    
    % ensure only one location per blob
    if length(x)>objectNum

        % get blob labels for each point
        labelFrame = bwlabel(frameFiltered>0);
        labelInds = sub2ind(size(labelFrame), y, x);
        blobLabels = labelFrame(labelInds);

        % find blobs containing multiple points
        [counts, bins] = hist(blobLabels, 1:max(blobLabels(:)));
        blobsWithMultiples = bins(counts>1);

        if ~isempty(blobsWithMultiples)

            % keep only the most anterior point within each blob
            validInds = ~ismember(blobLabels, blobsWithMultiples);

            for j = 1:length(blobsWithMultiples)
                [~, anteriorInd] = max( x .* (blobLabels==blobsWithMultiples(j)));    
                validInds(anteriorInd) = 1;
            end

            x = x(validInds);
            y = y(validInds);
            scores = scores(validInds);
        end
    end
    
%     % knn search
%     
%     frameFeatures = nan(size(features,1), length(x));
%     
%     for j = 1:length(x)
%         
%         xy = [x(j) y(j)];
%         imgInds = {xy(2)-centPad(1):xy(2)+centPad(1)-1, xy(1)-centPad(2):xy(1)+centPad(2)-1};
%         
%         if ~any(imgInds{1}<1 | imgInds{1}>vid.Height) && ~any(imgInds{2}<1 | imgInds{2}>vid.Width) % !!! eventually pad images to deal with edges
%             img = frame(imgInds{1}, imgInds{2});
%             frameFeatures(:,j) = img(:);
%         end
%     end
%     
%     ids = knnsearch(features', frameFeatures', 'K', 1);
%     isPaw = labels(ids);
%     
    % store data
    potentialLocationsBot(i).x = x;
    potentialLocationsBot(i).y = y;
    potentialLocationsBot(i).scores = scores;
    
    
    if showTracking
        
        % update figure
        set(rawIm, 'CData', frame);
        set(predictIm, 'CData', frameFiltered)
        set(scatterPts, 'XData', x, 'YData', y);
%         set(scatterKnn, 'XData', x(isPaw==1), 'YData', y(isPaw==1));
        
        % pause to reflcet on the little things...
        pause(.001);
    end
end

close all




