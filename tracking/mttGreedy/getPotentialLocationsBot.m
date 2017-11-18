function getPotentialLocationsBot(showTracking)

% settings
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\runBot.mp4';
classifier = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\pawBot.mat';
dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\';
objectNum = 4;

startFrame = 1;
overlapThresh = .8; % .5 for bottom // 
scoreThresh = 0;


% load classifier
load(classifier, 'model', 'subHgt', 'subWid')


% initializations
vid = VideoReader(vidFile);
sampleFrame = rgb2gray(read(vid,startFrame));
xMin = 20; % x and yMin are a temporary hack until i crop the videso properly
yMin = 15;
totalFrames = vid.NumberOfFrames;
cmap = winter(4);
kernel = reshape(model.w, subHgt, subWid);


% prepare figure
if showTracking

    figure('position', [680 144 698 834], 'menubar', 'none', 'color', 'black'); colormap gray

    rawAxis = subaxis(2,1,1, 'spacing', 0, 'margin', 0);
    rawIm = image(sampleFrame, 'parent', rawAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off');
    hold on; scatterPts = scatter(rawAxis, 0, 0, 200, 'filled');

    predictAxis = subaxis(2,1,2, 'spacing', 0.01, 'margin', .01);
    predictIm = image(sampleFrame, 'parent', predictAxis, 'CDataMapping', 'scaled');
    set(gca, 'visible', 'off');
end


potentialLocationsBot = struct();

for i = startFrame:totalFrames
    
    disp(i/totalFrames)
    
    % get frame and subframes
    frame = rgb2gray(read(vid,i));
    frame = getFeatures(frame);
    
    % filter with svm
    frameFiltered =  (conv2(double(frame), kernel, 'same') - model.rho);
    
    frameFiltered(frameFiltered < scoreThresh) = 0;
    frameFiltered(1:yMin,:) = 0;
    frameFiltered(:,1:xMin) = 0;
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
        pause(.05);
    end
end

save([dataDir 'potentialLocationsBot.mat'], 'potentialLocationsBot');
close all




