function trackObstacle(session)

    % load video
    dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';
    vid = VideoReader([dataDir session '\runBot.mp4']);
    
    % settings
    pixThresh = 128;
    obsThresh = .8;
    filtWidth = 10;
    
    % initialization
    numberOfFrames = vid.NumberOfFrames;
    currentFrame = 30000;
    obsPos = nan(1, numberOfFrames);
    
    
    % prepare figure
    close all;
    figure('units', 'normalized', 'position', [0.4714    0.1435    0.3422    0.7491]);
    
    subaxis(2,1,1, 'spacing', 0.01, 'margin', .01);
    frame = rgb2gray(read(vid, currentFrame));
    
    rawIm = imshow(frame);
    subaxis(2,1,2, 'spacing', 0.01, 'margin', .01);    
    threshIm = imshow(frame>pixThresh); hold on      
    
    obsLine = line([0 0], [1 size(frame,1)], 'Visible', 'off', 'lineWidth', 5, 'color', 'red');
    sliderThresh = uicontrol('Style', 'slider', 'Position', [0 0 100 25],...
                             'Min', 0, 'Max', 255, 'Value', pixThresh, 'Callback', @threshUpdate);


    

    for i = 1:numberOfFrames

        % acquire and process frame
        frame = rgb2gray(read(vid, currentFrame));
        frameThreshed = double(frame>pixThresh);
        
        % compute stats across rows
        colSum = sum(frameThreshed) / size(frame,1);
        colSum = smooth(colSum, filtWidth);
        obsLocation = median(find(colSum>obsThresh));
        
        
        % update preview
        set(rawIm, 'CData', frame);
        set(threshIm, 'CData', frameThreshed);
        
        currentFrame = currentFrame + 1;
        
        
        % update lines
        if ~isnan(obsLocation)
            set(obsLine, 'XData', [obsLocation obsLocation], 'Visible', 'on');
        else
            set(obsLine, 'Visible', 'off');
        end
        
        pause(.01)
        
        
        
    end
    
    
    function threshUpdate(~,~)
        keyboard
        pixThresh = get(sliderThresh, 'Value');
        fprintf('pixThresh set to: %.1f\n', pixThresh);
    end
end




