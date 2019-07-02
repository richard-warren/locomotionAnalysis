function getTrialsFromVideo(fileIn, fileOut, opts)

    % given a video with an LED that indicates trial start time, edits the
    % video to only include moments directly before and after start time
    
    % user settings
    s.timePre = .5;  % (s) show this much time before trial start
    s.timePost = .5;  % (s) show this much time after trial start
    s.fs = 250;  % original speed of video
    s.outputSpeed = .15;  % factor by which to slow down playback
    s.trialBuffer = 3; % disallow trials this many seconds within one another
    s.zScoreThresh = 6;  % how much does the ROI have to change for the trial to be detected
    s.roiPositions = []; % initial positions of ROI (Nx2 matrix of XY pairs describing positions of corners of polygon)
    s.previewAnalysis = true; % if true, shows the GUI, otherwise runs the analysis with specified settings
    
    % reassign settings contained in opts
    if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end

    % other variables
    currentFrame = 2;
    vid = VideoReader(fileIn);
    mask = ones(vid.Height,vid.Width);

    playing = true;
    diffHistory = [0 nan(1, vid.NumberOfFrames-1)];
    prevFrame = rgb2gray(read(vid,1));
    
    
    % set up figure
    previewFig = figure('Position', [100, 100, 400, 400], 'menubar', 'none', 'Visible', 'off');
    vidPreview = subaxis(1, 1, 1, 'spacing', 0.01, 'margin', .05);
    frame = read(vid,1);
    imPreview = imshow(rgb2gray(frame));
    
    if exist('prevSettings.mat') && isempty(s.roiPositions)
        load('prevSettings.mat', 'roiPositions');
        roi = impoly(vidPreview, roiPositions);
    elseif ~isempty(s.roiPositions)
            roi = impoly(vidPreview, s.roiPositions);
    else
        roi = impoly(vidPreview);
    end
    
    mask = createMask(roi);
    addNewPositionCallback(roi, @updateRois);
    
    % set up buttons
    uicontrol('Style', 'pushbutton', 'String', 'start analysis', 'Position', [0 0 100 20], 'Callback', @startAnalysis);
    uicontrol('Style', 'pushbutton', 'String', 'close', 'Position', [100 0 100 20], 'Callback', @closeAll);
    
    % preview analysis
    if s.previewAnalysis
        previewFig.Visible = 'on';
        while playing
            % get and show frame
            frame = rgb2gray(read(vid,currentFrame));
            set(imPreview, 'Cdata', frame);
            prevFrame = frame;
            if currentFrame<vid.NumberOfFrames; currentFrame=currentFrame+1; else; currentFrame = 1; end
            pause(.01);
        end
    else
        startAnalysis();
    end
    close(previewFig)
    
    
    % ---------
    % FUNCTIONS
    % ---------
    
    function closeAll(~,~)
        playing = false;
    end

    function updateRois(~)
        mask = createMask(roi);
    end

    function startAnalysis(~,~)
        
        % save roi position
        roiPositions = getPosition(roi);
        save('prevSettings.mat', 'roiPositions');
        
        playing = false;
        
        % get all frame differences
        f=1;
        prevFrame = rgb2gray(read(vid,f));
        w = waitbar(0, 'calculating frame differences...');
        while f<vid.NumberOfFrames
            f=f+1;
            frame = rgb2gray(read(vid,f));
            pixelDifs = double(frame).*double(mask) - double(prevFrame).*double(mask);
            frameDif =  sum(abs(pixelDifs(:))) / sum(mask(:));
            diffHistory(f) = frameDif;
            prevFrame = frame;
            waitbar(f/vid.NumberOfFrames)
        end
        close(w)
        
        % get trial start times
        trialStartInds = find(zscore(diffHistory)>s.zScoreThresh);
        trialStartInds = trialStartInds(logical([1 diff(trialStartInds/s.fs)>s.trialBuffer])); % remove trial start times that are too close together
        fprintf('number of trials found: %i\n', length(trialStartInds));
        
        % create video
        writerobj = VideoWriter(fileOut, 'MPEG-4');
        set(writerobj, 'FrameRate', round(s.fs*s.outputSpeed))
        open(writerobj);
        
        w = waitbar(0, 'writing video...');
        % write video to file
        for t = 1:length(trialStartInds)
            for f = round(trialStartInds(t)-s.timePre*s.fs):round(trialStartInds(t)+s.timePost*s.fs)
                frame = read(vid,f);
                frame = insertText(frame, [0 0], ['trial: ' num2str(t)]);
                writeVideo(writerobj, frame)
            end
            % add blank frame
            writeVideo(writerobj, zeros(size(frame)))
            waitbar(t/length(trialStartInds))
        end
        
        % close windows
        close(w)
        close(writerobj)
    end
end


