function makeVidWisk(vidName, session, varargin)

% edits a video of mouse jumping over obstacles s.t. obstacle trials are
% kept and everything else is edited out. obsPosRange is in m and defines
% the start and end position of the obstacle along the track that each
% trial should include.
% namePrefix is added to the beginning of the file name... kind of weird i know, sorry bro


% settings
s.maxTrialTime = 10;  % (s) trials exceeding maxTrialTime will be trimmed to this duration
s.border = 4;  % (pidels) thickness of border surrounding whisker frame
s.contrastLims = [.1 .9];  % (0->1) contrast limits for video

s.includeWiskCam = true;  % whether to add whisker camera
s.showPawTouches = true;  % whether to color the frame when a paw contact occurs
s.showTrialInfo = false;  % whether to show trial metadata as text
s.showWiskTouches = true;  % whether to highlight when whisker contact occurs
s.showObsOn = false;  % whether to add text showing when obs turns on

s.obsPosRange = [-.05 .1];  % (m) for each trial, show when obs is within this range of the mouse's nose
s.playBackSpeed = .15;  % fraction of real time speed for playback
s.trialLabels = {};  % cell array of trial labels - can be useful for distinguishing different trial types // if providing trialLabels, you also need to provide 'trials', a list of the trial numbers you would like to show

s.trialNum = 10;  % number of trials (evenly spaced throughout session) to show
s.trials = [];  % array of specific trials to show // if provided, s.trialNum is ignored



% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4'));
if s.includeWiskCam; vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4')); end

load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'obsPositionsFixed', 'obsTimes', 'obsPixPositions', 'wheelPositions', 'wheelTimes', 'isLightOn', ...
    'obsOnTimes', 'frameTimeStamps', 'frameTimeStampsWisk');

if s.showPawTouches
    load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'touchesPerPaw')
    isTouching = any(touchesPerPaw,2);
end

if s.showWiskTouches; load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'wiskContactFrames'); end


% get position where wisk frame should overlap with runTop frame
if s.includeWiskCam
    [frame, yWiskPos, xWiskPos, wiskScaling] = ...
        getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, find(frameTimeStamps>obsOnTimes(1), 1, 'first'));  % use first frame where obstacle is on to ensure mouse is on the wheel when the whisker cam position is determined
else
    frame = read(vid, 1);
end
frameDim = [size(frame,1) size(frame,2) 3];


% determine video settings
vidSetting = 'MPEG-4';
fps = round(vid.FrameRate * s.playBackSpeed);
maxFps = 150; % fps > 150 can be accomplished using 'Motion JPEG AVI' as second argument to VideoWriter, but quality of video is worse

if fps>maxFps
    fprintf('WARNING: changing video mode to ''Motion JPEG AVI'' to acheive requested playback speed\n');
    vidSetting = 'Motion JPEG AVI';
end

vidWriter = VideoWriter(vidName);
set(vidWriter, 'FrameRate', fps)
if strcmp(vidSetting, 'MPEG-4'); set(vidWriter, 'Quality', 50); end
open(vidWriter)


% determine trials to include
textprogressbar([session ': editing video...']);

if isempty(s.trials)
    s.trials = floor(linspace(1, length(obsOnTimes), min(s.trialNum, length(obsOnTimes))));
end




% edit video
for i = s.trials
    
    % find trial inds
    obsAtNoseTime = obsTimes(find(obsPositionsFixed>=0 & obsTimes>obsOnTimes(i), 1, 'first'));
    obsAtNosePos = wheelPositions(find(wheelTimes>=obsAtNoseTime,1,'first'));
    inds = find((wheelPositions > obsAtNosePos+s.obsPosRange(1)) & (wheelPositions < obsAtNosePos+s.obsPosRange(2)));
    startInd = inds(1);
    endInd = inds(end);
    endTime = min(wheelTimes(startInd)+s.maxTrialTime, wheelTimes(endInd));
    trialInds = find(frameTimeStamps>wheelTimes(startInd) & frameTimeStamps<endTime);
    
    
    if isempty(trialInds) % if a block has NaN timestamps (which will happen when unresolved), startInd and endInd will be the same, and frameInds will be empty
        fprintf('skipping trial %i\n', i)
    else
        
        for j = trialInds'
            
            if s.includeWiskCam
                frame = getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, j, ...
                    'yWiskPos', yWiskPos, 'xWiskPos', xWiskPos, 'wiskScaling', wiskScaling, ...
                    'runContrast', s.contrastLims);
                frame = repmat(frame, 1, 1, 3);  % add color dimension
            else
                frame = read(vid, j);
            end
            
            % change color of frame if touching
            if s.showPawTouches
                if isTouching(j)
                	frame(:,:,3) = frame(:,:,3)*.2;
                end
            end
            
            % add trial info text
            if s.showTrialInfo
                if isLightOn(i); lightText = 'light on'; else; lightText = 'light off'; end
                wiskFrameInd = find(frameTimeStampsWisk==frameTimeStamps(j), 1, 'first');
                framesFromContact = nan;
                if ~isempty(wiskFrameInd); framesFromContact = wiskFrameInd - wiskContactFrames(i); end
                text = sprintf('%s, trial %i, touch frames %i, %s, %i', ...
                    session, i, sum(isTouching(trialInds)), lightText, framesFromContact);
                frame = insertText(frame, [size(frame,2) size(frame,1)], text,...
                                   'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
            end
            
%             % add trial condition info
%             if ~isempty(s.trialLabels)
%                 if trialInds(i)==1
%                     boxColor = 'yellow';
%                     textColor = 'white';
%                 else
%                     boxColor = 'blue';
%                     textColor = 'white';
%                 end
%                 frame = insertText(frame, [size(frame,2), 0], trialLabels{trialInds(i)},...
%                                    'BoxColor', boxColor, 'anchorpoint', 'RightTop', 'textcolor', textColor);
%                 frame = insertText(frame, [size(frame,2), size(frameTop,1)+size(frameBot,1)], trialLabels{trialInds(i)},...
%                                    'BoxColor', boxColor, 'anchorpoint', 'RightTop', 'textcolor', textColor);
%             end
            
            % show when obstacle turns on
            if s.showObsOn
                if frameTimeStamps(j)>=obsOnTimes(i)
                    frame = insertText(frame, [0 size(frame,1)], 'OBSTACLE ON', 'BoxOpacity', 1, ...
                                   'BoxColor', 'white', 'anchorpoint', 'LeftBottom', 'textcolor', 'black');
                end
            end

            % write frame to video
            writeVideo(vidWriter, frame);
        end

        % add blank frame between trials
        writeVideo(vidWriter, zeros(size(frame)));
    end
    
    % update progress bar
    textprogressbar((i/length(obsOnTimes)) * 100)
end


textprogressbar(' done')
close(vidWriter)


