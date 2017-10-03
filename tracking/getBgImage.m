function bgImage = getBgImage(vid, meanFrameNum, showImage)
    % given a vid file, spits out a bg image that is the average of meanFrameNum frames
    % only uses frame where average pixel difference is > diffThresh to ensure it only takes frames where the mouse is moving

    
    % user settings
    medFiltKernel = [10 10];
    contrast = [.05 .5; 0 1]; % [low in, high in; low out, high out]
    maxTimeToSkip = 0; % skip this much time after a frame to avoid grabbing only successive frames
    diffThresh = 2*10e-4;
    
    % get mean background frame
    frameLast = readFrame(vid);
    meanFrame = zeros(size(frameLast)); % re-computes mean every time it grabs a new frame to avoid having to store many frames in memory simultaneously
    meanFrameCount = 0;
    
    w = waitbar(0, 'computing background frame...');   
    
    while meanFrameCount<meanFrameNum && hasFrame(vid)
        frame = readFrame(vid);
        frameDiff = abs(frame-frameLast);
        
        if mean(frameDiff(:))>diffThresh
            meanFrameCount = meanFrameCount+1;
            meanFrame = ((meanFrame*(meanFrameCount-1)) + double(frame)) / meanFrameCount;
            
            % skip frames to prevent grabbing only successive frames
            if maxTimeToSkip>0
                timeToSkip = randi([1 maxTimeToSkip*1000])/1000;
                vid.CurrentTime = vid.CurrentTime + timeToSkip;
            end
            
            waitbar(meanFrameCount/meanFrameNum)
        end

        frameLast = frame;
    end
    close(w)
    bgImage = uint8(meanFrame);
    
    if showImage; figure('color', [0 0 0]); imagesc(uint8(meanFrame)); set(gca, 'visible', 'off'); end
    if meanFrameCount<meanFrameNum; disp('COULD NOT FIND REQUESTED # OF FRAMES FOR BACKGROUND FRAME COMPUTATION!!!'); end        
end

    
    
    