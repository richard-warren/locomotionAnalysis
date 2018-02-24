function makeVidWisk(session, obsPosRange, playBackSpeed, trialProportion, trialLabels, trialInds)

% !!! need to document

% edits a video of mouse jumping over obstacles s.t. obstacle trials are
% kept and everything else is edited out. obsPosRange is in m and defines
% the start and end position of the obstacle along the track that each
% trial should include.


% settings
editedDir = [getenv('OBSDATADIR') 'editedVid\'];
maxTrialTime = 1.5; % trials exceeding maxTrialTime will be trimmed to this duration (s)
border = 4; % thickness (pixels) to draw around the wisk frame
scalings = .35 : .005 : .45; % the whisker vid is scaled by all of these values, and the scale that maximizes the correlation between the images is kept
obsBotThickness = 15;
wiskTouchThresh = -.8; % if wiskTouchSignal surpasses wiskTouchThresh, then wiskTouchPixels are drawn on wisk frame to show points of contact


% initializations
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
vidWisk = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runWisk.mp4']);
vidWeb = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\webCam.avi']);

frameDim = round([vidTop.Height + vidBot.Height, vidBot.Width + (vidWeb.Width * (vidBot.Height/vidWeb.Height))]);
webDim = [vidBot.Height, frameDim(2) - vidBot.Width];

vidSetting = 'MPEG-4';

fps = round(vidTop.FrameRate * playBackSpeed);
maxFps = 150; % fps > 150 can be accomplished using 'Motion JPEG AVI' as second argument to VideoWriter, but quality of video is worse

if fps>maxFps
    fprintf('WARNING: changing video mode to ''Motion JPEG AVI'' to acheive requested playback speed\n');
    vidSetting = 'Motion JPEG AVI';
end

vidWriter = VideoWriter(sprintf('%s%sspeed%.2f', editedDir, session, playBackSpeed), vidSetting);
set(vidWriter, 'FrameRate', fps)
if strcmp(vidSetting, 'MPEG-4'); set(vidWriter, 'Quality', 50); end
open(vidWriter)

load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsPositions', 'obsTimes', 'obsPixPositions',...
                                            'wheelPositions', 'wheelTimes',...
                                            'obsOnTimes', 'obsOffTimes',...
                                            'frameTimeStamps', 'frameTimeStampsWisk', 'webCamTimeStamps',...
                                            'touchSig', 'touchSigTimes',...
                                            'wiskTouchSignal', 'wiskTouchPixels', 'nosePos');


% convert wisk contacts to z scores
realInds = ~isnan(wiskTouchSignal);
normedReal = zscore(wiskTouchSignal(realInds));
wiskTouchSignal = nan(size(wiskTouchSignal));
wiskTouchSignal(realInds) = normedReal;

% obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes); % correct for drift in obstacle position readings
obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));

% get position where wisk frame should overlap with runTop frame
topInd = find(obsPixPositions>vidTop.Width-50 & obsPixPositions<vidTop.Width, 1, 'first'); % find a frame where the obstacle is at the right edge, eg within the whisker camera
frameTop = rgb2gray(read(vidTop, topInd));
wiskInd = find(frameTimeStampsWisk==frameTimeStamps(topInd), 1, 'first');
frameWisk = rgb2gray(read(vidWisk, wiskInd));
[yWiskPos, xWiskPos, wiskScaling] = getSubFramePosition(frameTop, frameWisk, scalings);

% edit video
w = waitbar(0, 'editing video...');

for i = 1 : round(1/trialProportion) : length(obsOnTimes)
    
    % find trial indices
    startInd = find(obsTimes>obsOnTimes(i)  & obsPositions>=obsPosRange(1), 1, 'first');
    endInd   = find(obsTimes<obsOffTimes(i) & obsPositions<=obsPosRange(2), 1, 'last');
    
    % get frame indices
    endTime = min(obsTimes(startInd)+maxTrialTime, obsTimes(endInd));
    frameInds = find(frameTimeStamps>obsTimes(startInd) & frameTimeStamps<endTime);
    
    % get webCame frame indices
    webFrameInds = find(webCamTimeStamps>obsTimes(startInd) & webCamTimeStamps<endTime);
    webFrames = read(vidWeb, [webFrameInds(1) webFrameInds(end)]);
    webFrames = squeeze(webFrames(:,:,1,:)); % collapse color dimension

    % interpolate webFrames to number of inds in frameInds
    webFramesInterp = interp1(webCamTimeStamps(webFrameInds), 1:length(webFrameInds), frameTimeStamps(frameInds), 'nearest', 'extrap');
    
    if isempty(frameInds) % if a block has NaN timestamps (which will happen when unresolved), startInd and endInd will be the same, and frameInds will be empty
        fprintf('skipping trial %i\n', i)
    else
        
        for j = 1:length(frameInds)
            
            frame = uint8(zeros(frameDim));
            
            % top
            frameTop = rgb2gray(read(vidTop, frameInds(j)));
            frame(1:vidTop.Height, 1:vidTop.Width) = frameTop;
            
            % bot
            frameBot = rgb2gray(read(vidBot, frameInds(j)));
            frameBot = addObsToFrame(frameBot, obsPixPositions(frameInds(j)), obsBotThickness, [1, size(frameBot,1)], 255);
            frame(vidTop.Height+1:end, 1:vidBot.Width) = frameBot;
            
            % wisk
            wiskFrameInd = find(frameTimeStampsWisk==frameTimeStamps(frameInds(j)), 1, 'first');
            
            if ~isempty(wiskFrameInd) % timeDif < .01 % only write wisk frame if it is temporally close to run frame
                
                % get wisk frame
                frameWisk = rgb2gray(read(vidWisk, wiskFrameInd));
                               
                % resize, adjust contrast, and draw border
                frameWisk = imresize(frameWisk, wiskScaling);
                frameWisk = imadjust(frameWisk, [.75 .95], [0 1]);
                frameWisk = 255 - frameWisk;
                frameWisk([1:border, end-border:end], :) = 255;
                frameWisk(:, [1:border, end-border:end]) = 255;
                
                % incorporate into frame
                frame(yWiskPos:yWiskPos+size(frameWisk,1)-1, xWiskPos:xWiskPos+size(frameWisk,2)-1) = frameWisk;
            end
            
            % webCam
            frameWeb = webFrames(:,:,webFramesInterp(j));
            frameWeb = imresize(frameWeb, webDim);
            frame(vidTop.Height+1:end, vidBot.Width+1:end) = frameWeb;          
            
            % add trial info text
            frame = insertText(frame, [size(frame,2) size(frame,1)], num2str(i),...
                               'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
            
            % add trial condition info
            if exist('trialLabels', 'var')
                if trialInds(i)==1
                    boxColor = 'yellow';
                    textColor = 'white';
                else
                    boxColor = 'blue';
                    textColor = 'black';
                end
                frame = insertText(frame, [size(frame,2), 0], trialLabels{trialInds(i)},...
                                   'BoxColor', boxColor, 'anchorpoint', 'RightTop', 'textcolor', textColor);
                frame = insertText(frame, [size(frame,2), size(frameTop,1)+size(frameBot,1)], trialLabels{trialInds(i)},...
                                   'BoxColor', boxColor, 'anchorpoint', 'RightTop', 'textcolor', textColor);
            end
            
            
            % change color of frame if touching
%             currentTouch = interp1(touchSigTimes, touchSig, frameTimeStamps(frameInds(j)));
%             if currentTouch
%                 frame(:,:,2) = frame(:,:,1)*.15;
%             end
            
            % make wisk view yellow if contacting obstacle
            if ~isempty(wiskFrameInd)
                if wiskTouchSignal(wiskFrameInd) > wiskTouchThresh
                    frame(yWiskPos:yWiskPos+size(frameWisk,1)-1, xWiskPos:xWiskPos+size(frameWisk,2)-1, 3) = ...
                        frame(yWiskPos:yWiskPos+size(frameWisk,1)-1, xWiskPos:xWiskPos+size(frameWisk,2)-1, 3) * .2;
                end
            end

            % write frame to video
            writeVideo(vidWriter, frame);
        end

        % add blank frame between trials
        writeVideo(vidWriter, zeros(size(frame)));
    end
    
    % update waitbar
    waitbar(i/length(obsOnTimes))
    
end


close(w)
close(vidWriter)


