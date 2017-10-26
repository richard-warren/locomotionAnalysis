
% USER SETTINGS

% settings
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\171025_002\runBot.mp4';
classifier = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\classifiers\classifier_pawBot_26-Oct-2017';
xRes = 10; % sliding window resolution
yRes = 10; % sliding window resolution
startFrame = 6000;

% load sample frame and sample sub-frame
vid = VideoReader(vidFile);
sampleFrame = read(vid,startFrame-1);
frmHgt = size(sampleFrame,1);
frmWid = size(sampleFrame,2);

% load classifier
load(classifier)

% prepare figure
figure('units', 'normalized', 'position', [0 0 1 1]);
rawPlot = subaxis(2,1,1, 'spacing', 0.01, 'margin', .01);
predictPlot = subaxis(2,1,2);
%
for i = startFrame:vid.numberofframes
    % get frame and sub-frames
    frame = rgb2gray(read(vid,i));
    tic
    [frameFeatures xVals yVals] = parseFrameToFeatures(frame, xRes, yRes, subHgt, subWid);
    
    % get predicted labels and reshape
    predictions = predict(classifier, frameFeatures);
%     margins = margin(classifier, frameFeatures, predictions);
%     predictions = margins;
    predictions = reshape(predictions(:,1), length(xVals), length(yVals))';
    
    % interpolate to original spatial resolution
    [X, Y] = meshgrid(xVals+round(subWid/2), yVals+round(subHgt/2));
    [Xq, Yq] = meshgrid(1:frmWid, 1:frmHgt);
    
    predictions = interp2(X, Y, predictions, Xq, Yq)-1;
    predictions(isnan(predictions)) = 0;
    predictions = logical(round(predictions));
    toc
    
    % keep only the four largest blobs
    blobInfo = regionprops(predictions);
    
    % sort blobInfo by area
    if ~isempty(blobInfo)
        [~, inds] = sortrows(cat(1,blobInfo(:).Area),1);
        blobInfo = blobInfo(inds);
        if length(blobInfo)>4; blobInfo = blobInfo(1:4); end
    end

    % update figure
    imshow(getFeatures(frame), 'parent', rawPlot);
    imshow(predictions, [0 1], 'parent', predictPlot);
    
    % draw circles
    for j = 1:length(blobInfo)
        center = blobInfo(j).Centroid;
        rectangle('position', [center(1)-round(.5*subWid) center(2)-round(.5*subHgt) subWid subHgt], 'edgecolor', 'white', 'linewidth', 4, 'parent', rawPlot);
    end

    % pause to reflcet on the little things...
    pause(.001);
end





