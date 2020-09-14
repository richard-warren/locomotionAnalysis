function makeVidUnheadfixed(vidName, session, varargin)

% edits a video of mouse jumping over obstacles s.t. obstacle trials are
% kept and everything else is edited out // see abagillion option below


% settings
s.visible = 'on';           % ('on' or 'off') whether to show frames while writing video
s.border = 4;               % (pixels) thickness of border surrounding whisker frame
s.contrastLims = [.1 .9];   % (0->1) contrast limits for video
s.playBackSpeed = .15;      % fraction of real time speed for playback
s.dropFrames = 1;           % drop every s.dropFrames frames to keep file sizes down
s.includeWiskCam = true;    % whether to add whisker camera
s.text = '';                % text to add to bottom right corner
s.textArgs = {};            % additional args to pass to the text() function
s.insertTrialInfo = false;  % whether to insert trial info into bottom right corner
s.maxTrialTime = 2;         % (s) time trials that exceed this amount of time

s.trialNum = 10;            % number of trials (evenly spaced throughout session) to show
s.trials = [];              % array of specific trials to show // if provided, s.trialNum is ignored
s.xLims = [-.35 .15];       % (m) x limits relative to obstacle position
s.blankTime = .15;          % (s) how many seconds of black frames (or fadeout time) to put in between trials
s.noWheelBreaks = true;     % whether to exclude trials with wheel breaks (only applies when 'trials' NOT provided

s.obsColor = [1 .7 0];
s.obsDiam = .004;
s.topHgt = .045;            % (m) height of the top view in m (used to determinte where to draw the bottom obstacle)

s.showTracking = false;     % whether to overlay tracking
s.numCircles = 1;           % number of circles to show for each feature (there will be a trail of circles showing the tracking over time)
s.circSeparation = 2;       % how many frames separating circles in trail
s.circSz = 80;
s.lineWidth = 2;
s.featuresToShow = {'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH', 'tailBase', 'tailMid'};  % features to show (excluding _top and _bot suffix)
s.colors = [];



% initializations
fprintf('writing %s, trial: ', vidName)
if ~exist(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4'), 'file'); concatTopBotVids(session); end  % temp
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4'));
if s.includeWiskCam; vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4')); end
circIndOffsets = -(s.numCircles-1)*s.circSeparation : s.circSeparation : 0;

load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'wheelPositions', 'wheelTimes', 'wheelToObsPixPosMappings', ...
    'obsOnTimes', 'frameTimeStamps', 'frameTimeStampsWisk', 'pixelsPerM', ...
    'wheelCenter', 'wheelRadius', 'obsHeights', 'isWheelBreak');


% get position where wisk frame should overlap with runTop frame
if s.includeWiskCam
    [frame, yWiskPos, xWiskPos, wiskScaling] = ...
        getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, find(frameTimeStamps>obsOnTimes(1), 1, 'first'));  % use first frame where obstacle is on to ensure mouse is on the wheel when the whisker cam position is determined
else
    frame = read(vid, 1);
end
frameX = [0:size(frame,2)-1] / pixelsPerM;  % x grid for image in meters
frameY = [0:size(frame,1)-1] / pixelsPerM;  % x grid for image in meters

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
    if s.noWheelBreaks
        s.trials = sort(datasample(find(~isWheelBreak), min(s.trialNum, sum(~isWheelBreak)), 'replace', false))';  % excluding wheel break trials
    else
        s.trials = randsample(length(obsOnTimes), min(s.trialNum, length(obsOnTimes)))';  % excluding wheel break trials
    end
end


% set up figure
fig = figure('position', [0, 600, range(s.xLims) * pixelsPerM, size(frame,1)], ...
    'name', session, 'color', [0 0 0], 'menubar', 'none', 'visible', s.visible);
ax = axes('position', [0 0 1 1], 'CLim', [0 255]);
colormap gray
im = image(frame, 'CDataMapping', 'scaled', 'ydata', frameY); hold on;
set(ax, 'visible', 'off', 'XLim', s.xLims); hold on


% text
if ~isempty(s.text)
    text(s.xLims(1), frameY(end), s.text, 'color', 'white', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', s.textArgs{:});
end

if s.insertTrialInfo
    trialText = text(s.xLims(end), frameY(end), '', 'color', 'white', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', s.textArgs{:});
end


% load tracking
if s.showTracking
    locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
    scoreThresh = getScoreThresh(session, 'trackedFeaturesRaw_metadata.mat');  % scoreThresh depends on whether deeplabcut (old version) or deepposekit was used
    [locations, features] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM, 'scoreThresh', scoreThresh);
    
    % restrict features
    bins = contains(features, s.featuresToShow);
    locations = locations(:,:,bins); features = features(bins);
    
    % interpolate non-paw features
    nonPawInds = find(~contains(features, 'paw'));
    for i = nonPawInds; locations(:,:,i) = fillmissing(locations(:,:,i), 'spline'); end
    
    % define colors (s.t. same feature across views has same color)
    if ~isempty(s.colors); colorsTemp = s.colors; else; colorsTemp = hsv(length(s.featuresToShow)); end
    colors = nan(length(features),3);
    for i = 1:length(s.featuresToShow)
        bins =contains(features, s.featuresToShow{i});
        colors(bins,:) = repelem(colorsTemp(i,:),2,1);
    end

    circColors = repelem(colors, s.numCircles, 1);
    circColors = circColors .* repmat([ones(1,s.numCircles-1)*.5 1], 1, length(features))';  % darken trailing circles
    circSizes = [linspace(s.circSz*.1, s.circSz*.5, s.numCircles-1) s.circSz];
    circSizes = repmat(circSizes,1,length(features));
    
    scat = scatter(nan(1,length(features)*s.numCircles), ...
        nan(1,length(features)*s.numCircles), ...
        circSizes, circColors, 'filled', 'MarkerFaceAlpha', .8);
    
    lines = cell(1, length(features));
    for i = 1:length(features); lines{i} = plot(nan, nan, 'color', [colors(i,:) .8], 'LineWidth', s.lineWidth); end
end


% obstacle
obsTop = rectangle('position', [0-s.obsDiam/2, 0, s.obsDiam, s.obsDiam], 'curvature', [1 1], ...
    'facecolor', [s.obsColor .8], 'EdgeColor', 'none');
rectangle('position', [0-s.obsDiam/2, s.topHgt, s.obsDiam, frameY(end)-s.topHgt], ...
    'facecolor', [s.obsColor .8], 'EdgeColor', 'none');



% edit video
for i = s.trials
    fprintf('%i ', i)
    
    % load obstacle positions for trial
    obsPixPositions = polyval(wheelToObsPixPosMappings(i,:), wheelPositions);
    obsPixPositions = interp1(wheelTimes, obsPixPositions, frameTimeStamps);
    obsPositionsMeters = obsPixPositions / pixelsPerM;  % relative to left side of video frame
    
    % obstacle
    wheelTopZ = (wheelCenter(2)-wheelRadius)/pixelsPerM;
    obsZ = wheelTopZ - obsHeights(i)/1000;  % y axis is reversed
    try; set(obsTop, 'position', [0-s.obsDiam/2, obsZ, s.obsDiam, s.obsDiam]); catch; keyboard; end
    
    % find trial inds
    trialInds = find(obsPositionsMeters<(-s.xLims(1)+frameX(end)) & obsPositionsMeters>-s.xLims(2));
    trialInds = trialInds(1) : s.dropFrames : trialInds(end);
    trialInds = trialInds(frameTimeStamps(trialInds) - frameTimeStamps(trialInds(1)) < s.maxTrialTime);  % make sure doesn't exceed max trial time
    
    if s.showTracking
        trialLocations = locations/pixelsPerM;
        trialLocations(:,1,:) = trialLocations(:,1,:) - obsPositionsMeters;
    end
    
    % text
    if s.insertTrialInfo; set(trialText, 'String', sprintf('trial %i', i)); end
    
    for j = trialInds
        if s.includeWiskCam
            frame = getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, j, ...
                'yWiskPos', yWiskPos, 'xWiskPos', xWiskPos, 'wiskScaling', wiskScaling, ...
                'runContrast', s.contrastLims);
            frame = repmat(frame, 1, 1, 3);  % add color dimension
        else
            frame = read(vid, j);
        end

        % update figure
        set(im, 'XData', frameX - obsPositionsMeters(j), 'CData', frame);

        % tracking
        if s.showTracking
            % circles
            x = squeeze(trialLocations(j+circIndOffsets,1,:));
            y = squeeze(trialLocations(j+circIndOffsets,2,:));
            set(scat, 'XData', x(:), 'YData', y(:));
            
            % lines
            if j>trialInds(1)
                x = squeeze(trialLocations(trialInds(1):j,1,:));
                y = squeeze(trialLocations(trialInds(1):j,2,:));
                for k = 1:length(features); set(lines{k}, 'xdata', x(:,k), 'ydata', y(:,k)); end
            end
        end

        % write frame to video
        frame = getframe(fig);
        writeVideo(vidWriter, frame);
    end

    % fade out last frame
    if blankFrames>0; for k = linspace(1,0,blankFrames); writeVideo(vidWriter, frame.cdata.*k); end; end
    if s.showTracking; for k = 1:length(features); set(lines{k}, 'xdata', nan, 'ydata', nan); end; end  % clear tracking
end

fprintf('\nall done!\n')
close(vidWriter)
close(fig)


