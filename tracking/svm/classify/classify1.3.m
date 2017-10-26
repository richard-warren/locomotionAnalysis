
% USER SETTINGS
% sliding window resolution
xRes = 10;
yRes = 10;

% probability cloud thresholding settings
% thresh = -.1; % can go up to zero for more conservative detection
% erosion = 5; % pixel radius of erosion kernel shape
% erosionKernel = strel('diamond', erosion); % erosion kernel

% load sample frame and sample sub-frame
vid = VideoReader('video\topTest.mp4');
sampleFrame = read(vid,1);
frmHgt = size(sampleFrame,1);
frmWid = size(sampleFrame,2);

% load classifier
load('trainClassifier\classifiers\topPaws.mat')

% prepare figure
figure('units', 'normalized', 'position', [0 0 1 1]);
rawPlot = subaxis(2,1,1, 'spacing', 0.01, 'margin', .01);
predictPlot = subaxis(2,1,2);
%
for i = 1:vid.numberofframes
    % get frame and sub-frames
    frame = read(vid,i);
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





