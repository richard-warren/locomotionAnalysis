function potentialLocationsBot = getPotentialLocationsBot(vid, model1, model2, subFrameSize1, subFrameSize2, ...
    scoreThresh, obsPixPositions, frameInds, showTracking, threshIntensity, locationsTop)

% !!! need to document


% settings
overlapThresh = .7; % used for non-maxima suppression // higher numbers = more tightly packed
% objectNum = 4;      % number of paws
xNearness = 30;

% initializations
sampleFrame = rgb2gray(read(vid,1));
totalFrames = vid.NumberOfFrames;
kernel = reshape(model1.Beta, subFrameSize1(1), subFrameSize1(2));
bg = getBgImage(vid, 1000, 120, 2*10e-4, false);


% prepare figure
if showTracking
    
%     figure(); imagesc(-kernel);

    figure('position', [680 144 698 834], 'menubar', 'none', 'color', 'black'); colormap gray

    rawAxis = subaxis(2,1,1, 'spacing', 0, 'margin', 0);
    rawIm = image(sampleFrame, 'parent', rawAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off', 'CLim', [0 255]);
    hold on; scatter1 = scatter(rawAxis, 0, 0, 50, [1 1 1], 'filled');
    hold on; scatter2 = scatter(rawAxis, 0, 0, 150, [1 0 0], 'linewidth', 3);
    hold on; scatter3 = scatter(rawAxis, 0, 0, 100, [0 1 0], 'linewidth', 3);

    predictAxis = subaxis(2,1,2, 'spacing', 0.01, 'margin', .01);
    predictIm = image(sampleFrame, 'parent', predictAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off', 'CLim', [0 10]);
end


potentialLocationsBot(max(frameInds)) = struct();

for i = frameInds
    
    disp(i/totalFrames)
    
    % get frame and subframes
    frame = rgb2gray(read(vid,i));
    frame = getFeatures(frame);
    frame = frame - bg;
    if exist('threshIntensity', 'var')
        frame(frame>threshIntensity) = threshIntensity; % a hack to limit influence of markers shining in bottom view
    end
    
    % mask obstacle
    frame = maskObs(frame, obsPixPositions(i)); % !!! should replace this with addObsToFrame

    % filter with svm
    frameFiltered = -(conv2(double(frame)/model1.KernelParameters.Scale, kernel, 'same') + model1.Bias);
    
    % mask st only x locations close to x locations tracked in top view remain
    if exist('locationsTop', 'var')
        mask = zeros(size(frame));
        for j = 1:length(locationsTop(i).x)
            x = round(locationsTop(i).x(j));
            startInd = max(1, x - xNearness);
            endInd = min(vid.Width, x + xNearness);
            mask(:, startInd:endInd) = 1;
        end
        frameFiltered = frameFiltered .* mask;
    end
    
    
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
    
    
    
    % perform second round of classification (svm)
%     frameFeatures = nan(prod(subFrameSize2), length(x));
%     for j = 1:length(x)
%         img = getSubFrame(frame, [y(j) x(j)], subFrameSize2);
%         frameFeatures(:,j) = img(:);
%     end
%     
%     classes = predict(model2, frameFeatures');
%     isPaw = (classes==1);
    
    
    
    % perform second round of classification (cnn)
%     dims = model2.Layers(1).InputSize;
%     frameFeatures = nan(dims(1), dims(2), 3, length(x));
    frameFeatures = nan(model2.inputs{1}.size, length(x));
    
    for j = 1:length(x)
        img = getSubFrame(frame, [y(j) x(j)], subFrameSize2);
%         img = uint8(imresize(img, 'outputsize', model2.Layers(1).InputSize(1:2)));
%         img = repmat(img, 1, 1, 3);
%         frameFeatures(:,:,:,j) = img;

        img(end, end-1:end) = [x(j) y(j)];
        frameFeatures(:,j) = img(:);
    end
    
%     classes = classify(model2, frameFeatures);
    [~, classes] = max(model2(frameFeatures), [], 1);
    isPaw  = (classes==3);
    isHind = (classes==1);
    isFore = (classes==2);
    
    
    % store data
    potentialLocationsBot(i).x = x(isPaw);
    potentialLocationsBot(i).y = y(isPaw);
    potentialLocationsBot(i).scores = scores(isPaw);
    
    
    if showTracking
        
        % update figure
        set(rawIm, 'CData', frame);
        set(predictIm, 'CData', frameFiltered)
        set(scatter1, 'XData', x, 'YData', y);
        set(scatter2, 'XData', x(isHind), 'YData', y(isHind));
        set(scatter3, 'XData', x(isFore), 'YData', y(isFore));
        
        % pause to reflcet on the little things...
        pause(.001);
    end
end




% % ---------
% % FUNCTIONS
% % ---------
% 
% % keypress controls
% function changeFrames(~,~)
%     
%     key = double(get(fig, 'currentcharacter'));
%     
%     if ~isempty(key) && isnumeric(key)
%         
%         if key==28                      % LEFT: move frame backward
%             pause(.001);
%             paused = true;
%             updateFrame(-1);
%         
%         elseif key==29                  % RIGHT: move frame forward
%             pause(.001);
%             paused = true;
%             updateFrame(1);
%         
%         elseif key==27                  % ESCAPE: close window
%             playing = false;
%             paused = false;
%         else                            % OTHERWISE: close window
%             paused = ~paused;
%         end
%     end
% end




