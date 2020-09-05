function makeVid(vidName, session, varargin)

% edits a video of mouse jumping over obstacles s.t. obstacle trials are
% kept and everything else is edited out // see abagillion option below

% todo:
% have text be a cell array so users can pass in trial specific text
% paw contact visualization
% whisker contact visualization


% settings
s.visible = 'on';           % ('on' or 'off') whether to show frames while writing video
s.maxTrialTime = 10;        % (s) trials exceeding maxTrialTime will be trimmed to this duration
s.border = 4;               % (pixels) thickness of border surrounding whisker frame
s.contrastLims = [.1 .9];   % (0->1) contrast limits for video
s.playBackSpeed = .15;      % fraction of real time speed for playback
s.dropFrames = 1;           % drop every s.dropFrames frames to keep file sizes down

s.includeWiskCam = true;    % whether to add whisker camera
s.text = '';                % text to add to bottom right corner

s.trialNum = 10;            % number of trials (evenly spaced throughout session) to show
s.trials = [];              % array of specific trials to show // if provided, s.trialNum is ignored
s.rewardTrials = false;     % whether trials are considered reward delivery // otherwise considered obstacles
s.obsPosRange = [-.05 .1];  % (m) for each trial, show when obs is within this range of the mouse's nose
s.rewardWindow = [2 .5];    % (s) if rewardTrials==true, seconds relative to previous reward for trial start, and relative to current reward to trial end
s.blankTime = .15;          % (s) how many seconds of black frames (or fadeout time) to put in between trials

s.showTracking = false;     % whether to overlay tracking
s.numCircles = 4;           % number of circles to show for each feature (there will be a trail of circles showing the tracking over time)
s.circSeparation = 2;       % how many frames separating circles in trail
s.circSz = 80;
s.featuresToShow = {'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH', 'tailBase', 'tailMid'};  % features to show (excluding _top and _bot suffix)
s.colors = [];

s.speedNearContact = [];    % playBackSpeed immediatly after whisker contact // use this to create a cool slow down effect as the obstacle jump occurs
s.contactWindow = [0 .2];   % (s) seconds relative to whisker contact to start and stop slow down effect


% initializations
fprintf('writing %s, trial: ', vidName)
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if isempty(s.speedNearContact); s.speedNearContact = s.playBackSpeed; end
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4'));
if s.includeWiskCam; vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4')); end
circIndOffsets = -(s.numCircles-1)*s.circSeparation : s.circSeparation : 0;

load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'obsPositionsFixed', 'obsTimes', 'wheelPositions', 'wheelTimes', 'rewardTimes', ...
    'obsOnTimes', 'frameTimeStamps', 'frameTimeStampsWisk', 'pixelsPerM', 'wiskContactTimes');


% get position where wisk frame should overlap with runTop frame
if s.includeWiskCam
    [frame, yWiskPos, xWiskPos, wiskScaling] = ...
        getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, find(frameTimeStamps>obsOnTimes(1), 1, 'first'));  % use first frame where obstacle is on to ensure mouse is on the wheel when the whisker cam position is determined
else
    frame = read(vid, 1);
end



% determine video settings
vidSetting = 'MPEG-4';
fps = round(vid.FrameRate * s.playBackSpeed / s.dropFrames);
maxFps = 150; % fps > 150 can be accomplished using 'Motion JPEG AVI' as second argument to VideoWriter, but quality of video is worse
blankFrames = round(fps * s.blankTime);  % how many black frames to put in between trials;

if fps>maxFps
    fprintf('WARNING: changing video mode to ''Motion JPEG AVI'' to acheive requested playback speed\n');
    vidSetting = 'Motion JPEG AVI';
end

vidWriter = VideoWriter(vidName);
set(vidWriter, 'FrameRate', fps)
if strcmp(vidSetting, 'MPEG-4'); set(vidWriter, 'Quality', 50); end
open(vidWriter)


% determine trials to include
if isempty(s.trials)
    if s.rewardTrials; maxNum = length(rewardTimes)-1; else; maxNum = length(obsOnTimes); end
    s.trials = floor(linspace(2, maxNum, min(s.trialNum, maxNum)));
end


% set up figure
fig = figure('name', session, 'color', [0 0 0], 'position', [800, 50, size(frame,2), size(frame,1)], 'menubar', 'none', 'visible', s.visible);
ax = axes('position', [0 0 1 1], 'CLim', [0 255]);
colormap gray
im = image(frame, 'CDataMapping', 'scaled'); hold on;
set(ax, 'visible', 'off')


% load tracking
if s.showTracking
    locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
    scoreThresh = getScoreThresh(session, 'trackedFeaturesRaw_metadata.mat');  % scoreThresh depends on whether deeplabcut (old version) or deepposekit was used
    [locations, features, ~, ~, scores] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM, 'scoreThresh', scoreThresh);
    
    % restrict features
    bins = contains(features, s.featuresToShow);
    locations = locations(:,:,bins); features = features(bins); scores = scores(bins);
    
    % define colors (s.t. same feature across views has same color)
    if ~isempty(s.colors); colorsTemp = s.colors; else; colorsTemp = hsv(length(s.featuresToShow)); end
    colors = nan(length(features),3);
    for i = 1:length(s.featuresToShow)
        bins =contains(features, s.featuresToShow{i});
        colors(bins,:) = repelem(colorsTemp(i,:),2,1);
    end
    
    colors = repelem(colors, s.numCircles, 1);
    colors = colors .* repmat([ones(1,s.numCircles-1)*.5 1], 1, length(features))';  % darken trailing circles
    circSizes = [linspace(s.circSz*.1, s.circSz*.5, s.numCircles-1) s.circSz];
    circSizes = repmat(circSizes,1,length(features));
    
    scat = scatter(nan(1,length(features)*s.numCircles), ...
        nan(1,length(features)*s.numCircles), ...
        circSizes, colors, 'filled', 'MarkerFaceAlpha', .8);
end








% edit video
for i = s.trials
    fprintf('%i ', i)
    
    % find trial inds
    if ~s.rewardTrials
        obsAtNoseTime = obsTimes(find(obsPositionsFixed>=0 & obsTimes>obsOnTimes(i), 1, 'first'));
        obsAtNosePos = wheelPositions(find(wheelTimes>=obsAtNoseTime,1,'first'));
        inds = find((wheelPositions > obsAtNosePos+s.obsPosRange(1)) & (wheelPositions < obsAtNosePos+s.obsPosRange(2)));
        startInd = inds(1);
        endInd = inds(end);
        endTime = min(wheelTimes(startInd)+s.maxTrialTime, wheelTimes(endInd));
        trialInds = find(frameTimeStamps>wheelTimes(startInd) & frameTimeStamps<endTime);
    else
        trialInds = find(frameTimeStamps>rewardTimes(i-1)+s.rewardWindow(1) & ...
                         frameTimeStamps<rewardTimes(i)+s.rewardWindow(2));
    end
    trialInds = trialInds(1) : s.dropFrames : trialInds(end);
    
    % get 'slow down' inds
    if s.playBackSpeed ~= s.speedNearContact
        contactInds = wiskContactTimes(wiskContactTimes>frameTimeStamps(trialInds(1)) & ...
                                       wiskContactTimes<frameTimeStamps(trialInds(end)));  % inds of whisker contacts within trial
        for j = contactInds
            tLims = j + s.contactWindow;                             % start and end time of slow down window
            times = tLims(1) : (s.speedNearContact/fps) : tLims(2);  % these are the 'desired' frame times given the specified frameRate
            trialSlowInds = knnsearch(frameTimeStamps, times');      % these are the inds closest to the desired frame times

            % splice in slow down inds
            trialInds = [trialInds(trialInds<trialSlowInds(1)) trialSlowInds' trialInds(trialInds>trialSlowInds(end))];
        end
    end
    
    
    if ~isempty(trialInds) % if a block has NaN timestamps (which will happen when unresolved), startInd and endInd will be the same, and frameInds will be empty
        for j = trialInds
            
            if s.includeWiskCam
                frame = getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, j, ...
                    'yWiskPos', yWiskPos, 'xWiskPos', xWiskPos, 'wiskScaling', wiskScaling, ...
                    'runContrast', s.contrastLims);
                frame = repmat(frame, 1, 1, 3);  % add color dimension
            else
                frame = read(vid, j);
            end
            
            % add text
            if ~isempty(s.text)
                frame = insertText(frame, [size(frame,2) size(frame,1)], s.text, ...
                                   'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
            end
            
            % update figure
            set(im, 'CData', frame);
            
            % tracking
            if s.showTracking
                x = squeeze(locations(j+circIndOffsets,1,:));
                y = squeeze(locations(j+circIndOffsets,2,:));
                set(scat, 'XData', x(:), 'YData', y(:));
            end
            
            % write frame to video
            frame = getframe(fig);
            writeVideo(vidWriter, frame);
        end

        % fade out last frame
        if blankFrames>0; for k = linspace(1,0,blankFrames); writeVideo(vidWriter, frame.cdata.*k); end; end
    end
end
fprintf('\nall done!\n')
close(vidWriter)
close(fig)


