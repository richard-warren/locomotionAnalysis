function potentialLocationsBot = getPotentialLocationsBot(vid, model, subHgt, subWid, obsPixPositions, startFrame, showTracking)

% !!! need to document


% settings
overlapThresh = .5; % used for non-maxima suppression // determines how tightly packed potential locations can be (higher numbers = more tightly packed) // ranged from 0 to 1
scoreThresh = .5;    % only pixels above scoreThresh are potential paw locations
objectNum = 4;      % number of paws
xMin = 35;
yMax = 220;
obsMaskDepth = .3;
obsMaskWid = 18;


% initializations
obsPixLeft = floor(obsMaskWid/2);
obsPixRight = ceil(obsMaskWid/2);
sampleFrame = rgb2gray(read(vid,1));
totalFrames = vid.NumberOfFrames;
kernel = reshape(model.w, subHgt, subWid);


% prepare figure
if showTracking
    
    figure; imagesc(kernel);

    figure('position', [680 144 698 834], 'menubar', 'none', 'color', 'black'); colormap gray

    rawAxis = subaxis(2,1,1, 'spacing', 0, 'margin', 0);
    rawIm = image(sampleFrame, 'parent', rawAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off', 'CLim', [0 255]);
    hold on; scatterPts = scatter(rawAxis, 0, 0, 200, 'filled');

    predictAxis = subaxis(2,1,2, 'spacing', 0.01, 'margin', .01);
    predictIm = image(sampleFrame, 'parent', predictAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off', 'CLim', [0 10]);
end


potentialLocationsBot = struct();

for i = startFrame:totalFrames
    
    disp(i/totalFrames)
    
    % get frame and subframes
    frame = rgb2gray(read(vid,i));
    frame = getFeatures(frame);
    
    % mask obstacle
    obsPixMinMax = [obsPixPositions(i) - obsPixLeft, obsPixPositions(i) + obsPixRight - 1];
    if any(obsPixMinMax>0 & obsPixMinMax<=vid.width)
        obsPixMinMax(obsPixMinMax<1) = 1;
        obsPixMinMax(obsPixMinMax>vid.Width) = vid.Width;
        frame(:,obsPixMinMax(1):obsPixMinMax(2)) = frame(:,obsPixMinMax(1):obsPixMinMax(2)) .* obsMaskDepth;
    end

    % filter with svm
    frameFiltered = (conv2(double(frame), kernel, 'same') - model.rho);
    frameFiltered(:, 1:xMin) = 0;
    frameFiltered(yMax:end, :) = 0;
    
    frameFiltered(frameFiltered < scoreThresh) = 0;
    [x, y, scores] = nonMaximumSupress(frameFiltered, [subHgt subWid], overlapThresh);
    
    % ensure only one location per blob
    if length(x)>objectNum

        % get blob labels for each point
        labelFrame = bwlabel(frameFiltered>0);
        labelInds = sub2ind(size(labelFrame), y, x);
        labels = labelFrame(labelInds);

        % find blobs containing multiple points
        [counts, bins] = hist(labels, 1:max(labels(:)));
        blobsWithMultiples = bins(counts>1);

        if ~isempty(blobsWithMultiples)

            % keep only the most anterior point within each blob
            validInds = ~ismember(labels, blobsWithMultiples);

            for j = 1:length(blobsWithMultiples)
                [~, anteriorInd] = max( x .* (labels==blobsWithMultiples(j)));    
                validInds(anteriorInd) = 1;
            end

            x = x(validInds);
            y = y(validInds);
            scores = scores(validInds);
        end
    end
    
    % store data
    potentialLocationsBot(i).x = x;
    potentialLocationsBot(i).y = y;
    potentialLocationsBot(i).scores = scores;
    
    
    if showTracking
        
        % update figure
        set(rawIm, 'CData', frame);
        set(predictIm, 'CData', frameFiltered)
        set(scatterPts, 'XData', x, 'YData', y);
        
        % pause to reflcet on the little things...
        pause(.001);
    end
end

close all




