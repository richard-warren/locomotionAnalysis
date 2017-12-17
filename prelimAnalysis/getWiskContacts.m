

vid = VideoReader('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\wiskTest4\runWisk.mp4');
load('C:\Users\rick\Google Drive\columbia\obstacleData\sessions\wiskTest4\runAnalyzed.mat');
bg = 255 - getBgImage(vid, 500, true);

%%

% settings
pixPos = [100 vid.Width];
bgThresh = 60;
wiskThresh = 30;
erosion = 3;
dilation = 5;

% initializations
frameInds = find(obsPixPositionsWisk>pixPos(1) & obsPixPositionsWisk<pixPos(2));
close all; figure; pimpFig

for i = frameInds
    
    frame = 255 - rgb2gray(read(vid, i));
    erodeKern = strel('disk', erosion);
    dilateKern = strel('disk', dilation);

    % raw
    subplot(3,3,1)
    imshow(frame);

    % bg mask
    subplot(3,3,2)
    bgMask = frame > bgThresh;
    bgMask = imerode(bgMask, erodeKern);
    bgMask = imdilate(bgMask, dilateKern);
    imshow(bgMask)

    % masked
    subplot(3,3,3)
    masked = frame .* uint8(~bgMask);
    imshow(masked)

    % obs mask
    subplot(3,3,4)
    labelFrame = bwlabel(bgMask>0);
    firstColInds = find(labelFrame(:,1));
    invalidLabels = unique(labelFrame(firstColInds,1));
    % !!! keep only largest blob that exceeds some minimum
    obsMask = ~ismember(labelFrame, invalidLabels) & labelFrame>0;
    imshow(obsMask);

    % wisk threshed
    subplot(3,3,5)
    maskedThreshed = masked > wiskThresh;
    imshow(maskedThreshed)

    % combined
    subplot(3,3,6)
    combined = uint8(obsMask);
    combined(maskedThreshed) = 2;
    imagesc(combined);
    
    pause(.1)
end


