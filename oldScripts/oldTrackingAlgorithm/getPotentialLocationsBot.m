function potentialLocationsBot = getPotentialLocationsBot(vid, model1, model2, classNum, subFrameSize1, subFrameSize2, ...
    scoreThresh, obsPixPositions, frameInds, trialIdentities, showTracking)

% !!! need to document


% settings
overlapThresh = .7; % used for non-maxima suppression // higher numbers = more tightly packed
% featureSetting = 'alexNetNoLocations';
% xNearness = 30;

% initializations
sampleFrame = rgb2gray(read(vid,1));
kernel = reshape(model1.Beta, subFrameSize1(1), subFrameSize1(2));
bg = getBgImage(vid, 1000, 120, 2*10e-4, false);
cmap = hsv(classNum);
% featureLength = length(getSubFrameFeatures(zeros(subFrameSize2), [0 0], featureSetting));


% prepare figure
if showTracking
    
%     figure(); imagesc(-kernel);

    figure('position', [680 144 698 834], 'menubar', 'none', 'color', 'black'); colormap gray

    rawAxis = subaxis(2,1,1, 'spacing', 0, 'margin', 0);
    rawIm = image(sampleFrame, 'parent', rawAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off', 'CLim', [0 255]);
    hold on; scatterAll = scatter(rawAxis, 0, 0, 50, [1 1 1], 'filled'); % shows results of svm1
    
    scatterPaws = cell(1, classNum); % shows results of second round of classification
    for i = 1:classNum
        hold on; scatterPaws{i} = scatter(rawAxis, 0, 0, 150, cmap(i,:), 'linewidth', 3);
    end

    predictAxis = subaxis(2,1,2, 'spacing', 0.01, 'margin', .01);
    predictIm = image(sampleFrame, 'parent', predictAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off', 'CLim', [0 10]);
end



w = waitbar(0, 'getting potentialLocationsBot...', 'position', [1500 50 270 56.2500]);
wInd = 0;
totalFrames = vid.NumberOfFrames;
potentialLocationsBot(totalFrames) = struct();

% create isAnalyzed field (records which frames were analyzed)
temp = false(1,totalFrames);
temp(frameInds) = true;
temp = mat2cell(temp, 1, ones(totalFrames,1));
[potentialLocationsBot.isAnalyzed] = temp{:};

% save trial identities
temp = nan(1,totalFrames);
temp(frameInds) = trialIdentities;
temp = mat2cell(temp, 1, ones(totalFrames,1));
[potentialLocationsBot.trialIdentities] = temp{:};

for i = frameInds
    
    % get frame and subframes
    frame = rgb2gray(read(vid,i));
    frame = frame - bg;
    
    % mask obstacle
    frame = maskObs(frame, obsPixPositions(i)); % !!! should replace this with addObsToFrame

    % filter with svm
    frameFiltered = -(conv2(double(frame)/model1.KernelParameters.Scale, kernel, 'same') + model1.Bias);
    
    % mask st only x locations close to x locations tracked in top view remain
%     if exist('locationsTop', 'var')
%         mask = zeros(size(frame));
%         for j = 1:length(locationsTop(i).x)
%             x = round(locationsTop(i).x(j));
%             startInd = max(1, x - xNearness);
%             endInd = min(vid.Width, x + xNearness);
%             mask(:, startInd:endInd) = 1;
%         end
%         frameFiltered = frameFiltered .* mask;
%     end
    
    
    frameFiltered(frameFiltered < scoreThresh) = scoreThresh;
    frameFiltered = frameFiltered - scoreThresh;
    [x, y, scores] = nonMaximumSupress(frameFiltered, subFrameSize1, overlapThresh);
    
%     % ensure only one location per blob
%     if length(x)>objectNum
% 
%         % get blob labels for each point
%         labelFrame = bwlabel(frameFiltered>0);
%         labelInds = sub2ind(size(labelFrame), y, x);
%         blobLabels = labelFrame(labelInds);
% 
%         % find blobs containing multiple points
%         [counts, bins] = hist(blobLabels, 1:max(blobLabels(:)));
%         blobsWithMultiples = bins(counts>1);
% 
%         if ~isempty(blobsWithMultiples)
% 
%             % keep only the most anterior point within each blob
%             validInds = ~ismember(blobLabels, blobsWithMultiples);
% 
%             for j = 1:length(blobsWithMultiples)
%                 [~, anteriorInd] = max( x .* (blobLabels==blobsWithMultiples(j)));    
%                 validInds(anteriorInd) = 1;
%             end
% 
%             x = x(validInds);
%             y = y(validInds);
%             scores = scores(validInds);
%         end
%     end
    
    
    
    % perform second round of classification (svm2, single class)
%     frameFeatures = nan(featureLength, length(x));
%     for j = 1:length(x)
%         img = getSubFrame(frame, [y(j) x(j)], subFrameSize2);
%         frameFeatures(:,j) = img(:);
%     end
%     
%     classes = predict(model2, frameFeatures');
%     isPaw = (classes==1);



    % perform second round of classification (svm2, multi class, with locations)
%     frameFeatures = nan(featureLength, length(x));
%     
%     for j = 1:length(x)
%         img = getSubFrame(frame, [y(j) x(j)], subFrameSize2);
%         frameFeatures(:,j) = getSubFrameFeatures(img, [x(j) y(j)], 'imageWithLocations');
%     end
%     
%     classes = predict(model2, frameFeatures');



    % perform second round of classification (svm trained on alexNet extracted features)
%     frameFeatures = nan(featureLength, length(x));
%     
%     for j = 1:length(x)
%         img = getSubFrame(frame, [y(j) x(j)], subFrameSize2);
%         frameFeatures(:,j) = getSubFrameFeatures(img, [x(j) y(j)], 'alexNetFeaturesWithLocation');
%     end
%     
%     classes = predict(model2, frameFeatures');
    
    
    % perform second round of classification (cnn)
    dims = model2.Layers(1).InputSize;
    frameFeatures = nan(dims(1), dims(2), 3, length(x));
    
    for j = 1:length(x)
        img = getSubFrame(frame, [y(j) x(j)], subFrameSize2);
        img = uint8(imresize(img, 'outputsize', model2.Layers(1).InputSize(1:2)));
        img = repmat(img, 1, 1, 3);
        frameFeatures(:,:,:,j) = img;
    end
    
    classes = uint8(classify(model2, frameFeatures, 'ExecutionEnvironment', 'gpu'));
    classProbs = predict(model2, frameFeatures);



    % perform second round of classification (nn)
%     frameFeatures = nan(featureLength, length(x));
%     
%     for j = 1:length(x)
%         img = getSubFrame(frame, [y(j) x(j)], subFrameSize2);
%         frameFeatures(:,j) = getSubFrameFeatures(img, [x(j) y(j)], 1);
%     end
%     
%     nnScores = model2(frameFeatures);
%     [~, classes] = max(nnScores, [], 1);


    
    
    % store data
    potentialLocationsBot(i).x = x;
    potentialLocationsBot(i).y = y;
    potentialLocationsBot(i).scores = scores;
    potentialLocationsBot(i).class = classes;
    potentialLocationsBot(i).classProbabilities = classProbs;
    
    
    if showTracking
        
        % update figure
        set(rawIm, 'CData', frame);
        set(predictIm, 'CData', frameFiltered)
        set(scatterAll, 'XData', x, 'YData', y);
        
        for j=1:classNum
            set(scatterPaws{j}, 'XData', x(classes==j), 'YData', y(classes==j));
        end
        
        % pause to reflcet on the little things...
        pause(.001);
    end
    
    wInd = wInd+1;
    waitbar(wInd / length(frameInds))
end

close(w)


