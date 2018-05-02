function bgImage = getBgImage(vid, meanFrameNum, skipTime, diffThresh, showImage)
    % given a vid file, spits out a bg image that is the average of meanFrameNum frames
    % only uses frame where average pixel difference is > diffThresh to ensure it only takes frames where the mouse is moving

    
    % user settings
%     skipTime = 0; % time to skip at the beginning of the sesion, to account for periods where the mouse is being put on the wheel (s)
%     medFiltKernel = [10 5];
%     contrast = [.05 .5; 0 1]; % [low in, high in; low out, high out]
    maxTimeToSkip = 0; % skip this much time after a frame to avoid grabbing only successive frames
%     diffThresh = 2*10e-4;
    
    
    % initializations
    currentFrame = max(round(skipTime * vid.FrameRate), 2);
    frameLast = rgb2gray(read(vid, currentFrame-1));
    meanFrame = zeros(size(frameLast)); % re-computes mean every time it grabs a new frame to avoid having to store many frames in memory simultaneously
    meanFrameCount = 0;
    
%     w = waitbar(0, 'computing background frame...', 'position', [1500 50 270 56.2500]);
    
    
    
    
    
    % get mean background frame
    while meanFrameCount<meanFrameNum && currentFrame < vid.NumberOfFrames
        
        frame = rgb2gray(read(vid, currentFrame));
        frameDiff = abs(frame-frameLast);
        
        if mean(frameDiff(:))>diffThresh
            
%             frame = medfilt2(frame, medFiltKernel);
            
            meanFrameCount = meanFrameCount+1;
            meanFrame = ((meanFrame*(meanFrameCount-1)) + double(frame)) / meanFrameCount;
            
            % skip frames to prevent grabbing only successive frames
            if maxTimeToSkip>0
                timeToSkip = randi([1 maxTimeToSkip*1000])/1000;
                vid.CurrentTime = vid.CurrentTime + timeToSkip;
            end
            
%             waitbar(meanFrameCount/meanFrameNum)
        end

        frameLast = frame;
        currentFrame = currentFrame + 1;
    end
%     close(w)
%     meanFrame = getFeatures(meanFrame);
    bgImage = uint8(meanFrame);
    
    if showImage; figure('color', [0 0 0]); imshow(uint8(meanFrame)); set(gca, 'visible', 'off'); end
    if meanFrameCount<meanFrameNum; disp('COULD NOT FIND REQUESTED # OF FRAMES FOR BACKGROUND FRAME COMPUTATION!!!'); end        
end

    
    
    