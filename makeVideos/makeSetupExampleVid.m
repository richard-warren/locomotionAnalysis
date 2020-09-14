function makeSetupExampleVid(vidName, session, trials, varargin)
% make a wide vid with drawing of obs off screen to the right, intended to explain how obs appears and is controlled by the speed of the wheel
% should use trials in which there is variability in wheel speed, and ideally some in which mouse moves both forwards and backwards

% settings
s.prePostTime = [-.1 0];  % (s) time to add to beginning and end of a trial (before and after obs is engaged)
s.playBackSpeed = 0.15;
s.fps = 250;
s.obsRadius = 8;
s.contrastLims = [.2 1];  % pixels at these proportional values are mapped to 0 and 255
s.obsFadePixels = 20;     % when obs enter right side of vid frame, fade out drawing of obs over the course of this many pixels
s.obsColor = [1 1 1];
s.includeWiskCam = true;  % whether to add whisker camera
s.obsOffsetX = 0;         % (pixels) offset the x tracking of the obs in the top view to coimpensate for obs being tracking in the bottom view
s.obsOffsetZ = 0;         % (pixels) offset the z tracking of the obs to account for differences in the camera view...
s.text = '';              % text to include in bottom right corner of screen
s.textArgs = {};          % additional args to pass to the text() function


% initializations
fprintf('making setup example vid: %s...\n', vidName);
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'obsPixPositions', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'frameTimeStampsWisk');
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv'));
obsTopY = locationsTable.obs_top_1;  % used to determine height of obstacle
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4'));
if s.includeWiskCam; vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4')); end

% % get relationship between x pos of obs in top and bot views (this approach works with wisk cam, but fails with wisk cam because wisk cam ALSO has it's own x mapping...)
% obsBotX = mean([locationsTable.obsHigh_bot, locationsTable.obsLow_bot], 2);
% obsTopX = locationsTable.obs_top;
% bins = locationsTable.obs_top_2>=.99 & obsPixPositions'>vid.Width*.5 & obsPixPositions'<vid.Width;  % high confidence bins where obs is in right side of frame
% botToTopMapping = robustfit(obsBotX(bins), obsTopX(bins)); botToTopMapping = fliplr(botToTopMapping');
% obsPixPositionsTop = polyval(botToTopMapping, obsPixPositions);
% clear locationsTable

obsPixPositionsTop = obsPixPositions + s.obsOffsetX;

% get position where wisk frame should overlap with runTop frame
if s.includeWiskCam
    [frame, yWiskPos, xWiskPos, wiskScaling] = ...
        getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, find(frameTimeStamps>obsOnTimes(1), 1, 'first'));  % use first frame where obstacle is on to ensure mouse is on the wheel when the whisker cam position is determined
else
    frame = read(vid, 1);
end
frameDims = [size(frame,1), round(max(obsPixPositionsTop)*.85)];

% make video writer
vidWriter = VideoWriter(vidName, 'MPEG-4');
set(vidWriter, 'FrameRate', round(s.fps*s.playBackSpeed));
open(vidWriter);

% define lookup vector for obstacle opacity (it will fade out as it enters right side of screen)
transparencyGradient = ones(1, frameDims(2)*2); % make this unnecessarily long, trailing ones at the end
gradInds = round(vid.Width-.5*s.obsFadePixels) : round(vid.Width+.5*s.obsFadePixels);
transparencyGradient(1:gradInds(1)) = 0;
transparencyGradient(gradInds) = linspace(0,1,length(gradInds));


% iterate through trials
for i = 1:length(trials)
    trialBins = frameTimeStamps>=(obsOnTimes(trials(i))+s.prePostTime(1)) & ...
                frameTimeStamps<=(obsOffTimes(trials(i))+s.prePostTime(2));
    obsHeight = median(obsTopY(trialBins' & obsPixPositionsTop>vid.Width*.5 & obsPixPositionsTop<vid.Width));  % compute median when obs is near right side of frame
    obsHeight = obsHeight - s.obsOffsetZ;
    
    % iterate through frames within trial
    for j = find(trialBins)'
        
        % get vid frame and incorporate into wide frame
        if s.includeWiskCam
            vidFrame = getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, j, ...
                'yWiskPos', yWiskPos, 'xWiskPos', xWiskPos, 'wiskScaling', wiskScaling, ...
                'runContrast', s.contrastLims);
        else
            vidFrame = rgb2gray(read(vid, j));
            vidFrame = imadjust(vidFrame, s.contrastLims, [0 1]);
        end
        
        % insert frame into wider image
        frame = uint8(zeros(frameDims));
        frame(:,1:size(vidFrame,2)) = vidFrame;
        
        % add obs
        pixPos = round(obsPixPositionsTop(j));
        if pixPos>0 && pixPos<(frameDims(2)+2*s.obsRadius)
            opacity = transparencyGradient(pixPos);
            frame = insertShape(frame, 'FilledCircle', ...
                [pixPos obsHeight+s.obsRadius s.obsRadius], ...
                'opacity', opacity, 'color', s.obsColor*255);
        end
        
        % text
        if ~isempty(s.text)
            frame = insertText(frame, [0 size(frame,1)], s.text, ...  % [size(frame,2) size(frame,1)]
                'BoxColor', 'black', 'AnchorPoint', 'LeftBottom', 'TextColor', 'white', s.textArgs{:});
        end
        
        % write video
        writeVideo(vidWriter, frame);
    end
end
close(vidWriter)
disp('all done!')
