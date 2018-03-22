function getWiskContactTimes(sessions, thresh, showFrames)

% settings
rows = 5;
cols = 8;
contactPosLimits = [-.015 0]; % larger value 


for i = 1:length(sessions)
    
    fprintf('%s: getting wisk contact times, positions, and frames\n', sessions{i})
    
    vid = VideoReader([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runWisk.mp4']);
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'], ...
        'wiskTouchSignal', 'obsOnTimes', 'obsOffTimes', 'frameTimeStampsWisk', ...
        'obsPositions', 'obsTimes', 'obsPixPositions', 'nosePos', 'frameTimeStamps');
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));
    
    % convert wisk contacts to z scores
    realInds = ~isnan(wiskTouchSignal);
    normedReal = zscore(wiskTouchSignal(realInds));
    wiskTouchSignal = nan(size(wiskTouchSignal));
    wiskTouchSignal(realInds) = normedReal;
    
    % find first contact time and position for each trial and get sample frames
    contactTimes = nan(1, length(obsOnTimes));
    contactPositions = nan(1, length(obsOnTimes));
    preContactFrames = nan(vid.Height, vid.Width, length(obsOnTimes));
    contactFrames = nan(vid.Height, vid.Width, length(obsOnTimes));
    
        
    for j = 1:length(obsOnTimes)
        
        indStart = find(frameTimeStampsWisk>obsOnTimes(j) & frameTimeStampsWisk<obsOffTimes(j) & wiskTouchSignal>=thresh, 1, 'first');
        
        if ~isempty(indStart) && ~isnan(wiskTouchSignal(indStart-1))
                
            % get time of contact
            contactTimes(j) = interp1(wiskTouchSignal(indStart-1:indStart), frameTimeStampsWisk(indStart-1:indStart), thresh);

            % get position of contact
            if ~isnan(contactTimes(j))
                contactPositions(j) = interp1(obsTimes, obsPositions, contactTimes(j));

                % set to nan trials in which detected position is unreasonable
                if contactPositions(j)<contactPosLimits(1) || contactPositions(j)>contactPosLimits(2)
                    contactPositions(j) = nan;
                    contactTimes(j) = nan;
                end
            end

            % get contact frames
            if ~isnan(contactTimes(j))
                % find frame with time closest to contact time
                if abs(contactTimes(j)-frameTimeStampsWisk(indStart-1)) < abs(contactTimes(j)-frameTimeStampsWisk(indStart))
                    indStart = indStart-1;
                end
                preContactFrames(:,:,j) = rgb2gray(read(vid, indStart-1));
                contactFrames(:,:,j) = rgb2gray(read(vid, indStart));
            end        
        end
    end
    
    
    % save results
    save([getenv('OBSDATADIR') 'sessions\' sessions{i} '\wiskContactTimes.mat'], ...
        'contactPositions', 'contactTimes', 'preContactFrames', 'contactFrames', 'thresh');
    
    
    % show preContact contact frames
    if showFrames
        preContactPreview = nan(rows*vid.Height, cols*vid.Width);
        contactPreview = nan(rows*vid.Height, cols*vid.Width);
        imInds = randperm(size(contactFrames,3), rows*cols);
        imInd = 1;

        for j = 1:rows
            for k = 1:cols
                y = (j-1)*vid.Height + 1;
                x = (k-1)*vid.Width + 1;
                preContactPreview(y:y+vid.Height-1, x:x+vid.Width-1) = preContactFrames(:,:,imInds(imInd));
                contactPreview(y:y+vid.Height-1, x:x+vid.Width-1) = contactFrames(:,:,imInds(imInd));
                imInd = imInd+1;
            end
        end

        figure('name', [sessions{i} ' preContact']);
        imshow(uint8(preContactPreview)); pimpFig
        figure('name', [sessions{i} ' contact']);
        imshow(uint8(contactPreview)); pimpFig
    end
    fprintf('%s: getWiskContactTimes completed\n', sessions{i})
end



