function markerTracking(session)

    % load file
    dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';
    vid = VideoReader([dataDir session '\runBot.mp4']);

    % settngs
    avgBgFrames = 1000;
    startTime = 120;
    
    % initializations
    thresh = 40;
    currentFrame = startTime * vid.FrameRate;
    frame = rgb2gray(read(vid, currentFrame));
%     circRoiPts = [15 190; 227 122; 404 162];
%     wheelMask = getWheelMask(circRoiPts);
    bg = getBgImage(vid, avgBgFrames, true);


    % prepare figure
    figure; pimpFig;
    rawAx = subaxis(2,1,1, 'spacing', 0.01, 'margin', .01);
    rawIm = imshow(frame);
    threshAx = subaxis(2,1,2, 'spacing', 0.01, 'margin', .01);
    threshIm = imshow(frame>thresh);
    sliderThresh = uicontrol('Style', 'slider', 'Position', [0 0 100 25],...
                             'Min', 0, 'Max', 255, 'Value', thresh, 'Callback', @threshUpdate);
    
    % main loop
    while true
        
        % acquire and process frame
        frame = rgb2gray(read(vid, currentFrame));% .* wheelMask;
        frame = frame - bg;
        frame = medfilt2(frame, [10 5]); % [h w]
        
        % update preview
        set(rawIm, 'CData', frame);
        set(threshIm, 'CData', frame>thresh);
        pause(.01)
        currentFrame = currentFrame + 1;
        
    end
    
    
    
    % FUNCTIONS
    
    function threshUpdate(~,~)
        thresh = get(sliderThresh, 'Value');
    end

    function wheelMask = getWheelMask(circRoiPts)
        
        [wheelRadius, wheelCenter] = fitCircle(circRoiPts);
        wheelMask = uint8(ones(vid.Height, vid.Width));
        
        for y=1:vid.Height
            for x=1:vid.Width
                if norm(wheelCenter-[x;y]) < wheelRadius
                    wheelMask(y,x) = 0;
                end
            end
        end
    end
end


