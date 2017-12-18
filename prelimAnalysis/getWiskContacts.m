function [isWiskTouching, contactPixels] = getWiskContacts(vid, showTracking, frameTimeStampsWisk, obsOnTimes, obsOffTimes)

% !!! need to document

% settings
bgThresh = 60;
wiskThresh = 30;
erosion = 3;
dilation = 5;
minBlobArea = 2000;
borderThickness = 10;
borderThresh = .15;

% initializations
isWiskTouching = false(size(frameTimeStampsWisk));
contactPixels = cell(size(frameTimeStampsWisk));
erodeKern = strel('disk', erosion);
dilateKern = strel('disk', dilation);

% get inds of frames where obs is on
onInds = knnsearch(frameTimeStampsWisk, obsOnTimes);
offInds =  knnsearch(frameTimeStampsWisk, obsOffTimes);
isObsOn = zeros(size(frameTimeStampsWisk));
isObsOn(onInds) = 1;
isObsOn(offInds) = -1;
obsOnInds = find(logical(cumsum(isObsOn)))';

if showTracking
    
    figure('position', [1921 1 750 450], 'menubar', 'none', 'color', 'black');
    frame = rgb2gray(read(vid,1));
    bgMask = frame > bgThresh;
    
    subaxis(2,3,1, 'margin', 0, 'padding', 0, 'margin', 0)
    imRaw = imshow(frame);
    
    subaxis(2,3,2)
    imMasked = imshow(bgMask);
    
    subaxis(2,3,3)
    imObs = imshow(bgMask);
    
    subaxis(2,3,4)
    imWisk = imshow(bgMask);
    
    subaxis(2,3,5)
    imBorder = imshow(bgMask);
    
    subaxis(2,3,6)
    imTouching = imshow(bgMask);
end



for i = obsOnInds
    
%     disp(i/vid.NumberOfFrames)
    
    % raw
    frame = 255 - rgb2gray(read(vid, i));

    % masked by bg
    bgMask = frame > bgThresh;
    bgMask = imerode(bgMask, erodeKern);
    bgMask = imdilate(bgMask, dilateKern);
    masked = frame .* uint8(~bgMask);

    % obs mask
    tic
%     blobInfo = regionprops(bgMask, 'Area', 'BoundingBox', 'Image');
%     boundingBoxes = reshape([blobInfo.BoundingBox], 4, length(blobInfo));
%     validInds = ([blobInfo.Area] > minBlobArea) & (boundingBoxes(1,:)>1);
%     blobInfo = blobInfo(validInds);
%     obsMask = zeros(size(frame));
%     if ~isempty(blobInfo)
%         [~,maxInd] = max([blobInfo.Area]);
%         blobInfo = blobInfo(maxInd);
%         x = ceil(blobInfo.BoundingBox(1));
%         y = ceil(blobInfo.BoundingBox(2));
%         obsMask(y : y+blobInfo.BoundingBox(4)-1, x : x+blobInfo.BoundingBox(3)-1) = blobInfo.Image;
%     end

    keyboard    
    toc

    % wisk threshed
    wisks = masked > wiskThresh;
    
    % get mask for obstacle border
    borderMask = zeros(size(frame));
    
    if ~isempty(blobInfo)
        
        xP = x - borderThickness;
        yP = y - borderThickness;
        
        if xP>0 && yP>0
            obsMaskShifted = zeros(size(frame));
            obsMaskShifted(yP : yP+blobInfo.BoundingBox(4)-1, xP : xP+blobInfo.BoundingBox(3)-1) = blobInfo.Image;
            borderMask = xor(obsMask, obsMaskShifted) & ~obsMask;
        end
    end
    
    
    % get wisk touches, ie wisks masked by obstacle border
    wiskInBorder = wisks & borderMask;
    
    
    % store results
    if (sum(wiskInBorder(:)) / sum(borderMask(:))) > borderThresh
        isWiskTouching(i) = true;
    end
    contactPixels{i} = find(wiskInBorder);
    
    
    % update displays
    if showTracking
        
        if isWiskTouching(i)
            frame(wiskInBorder) = 255;
        end
        
        set(imRaw, 'CData', frame);
        set(imMasked, 'CData', masked);
        set(imObs, 'CData', obsMask);
        set(imWisk, 'CData', wisks);
        set(imBorder, 'CData', borderMask);
        set(imTouching, 'CData', wiskInBorder);
        
        if isempty(blobInfo)
            pause(.001)
        else
            pause(.05);
        end
    end
end


