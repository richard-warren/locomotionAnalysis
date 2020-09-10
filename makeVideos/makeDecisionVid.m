function makeDecisionVid(vidName, session, varargin)

% makes 'unheadfixed' video showing decision making strategies for a
% session // automatically selects exemplar trials for lengthened,
% shortened, and unmodified steps


% settings
s.visible = 'on';           % ('on' or 'off') whether to show frames while writing video
s.border = 4;               % (pixels) thickness of border surrounding whisker frame
s.contrastLims = [.1 .9];   % (0->1) contrast limits for video
s.playBackSpeed = .10;      % fraction of real time speed for playback
s.dropFrames = 1;           % drop every s.dropFrames frames to keep file sizes down
s.includeWiskCam = true;    % whether to add whisker camera
s.text = '';                % text to add to bottom right corner

s.xLims = [-.20 .10];       % (m) x limits relative to obstacle position
s.blankTime = .25;          % (s) how many seconds of black frames (or fadeout time) to put in between trials
s.contactPauseTime = 2;     % (s) time to pause at whisker contact

s.obsColor = [1 .7 0];
s.obsDiam = .004;           % (m)
s.topHgt = .045;            % (m) height of the top view in m (used to determinte where to draw the bottom obstacle)

s.colors = [hsv(2); .8 .8 .8];   % lengthened, shortened, unmodified



% initializations
fprintf('writing %s... ', vidName)
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4'));
if s.includeWiskCam; vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4')); end

load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'wheelPositions', 'wheelTimes', 'wheelToObsPixPosMappings', 'isLightOn', ...
    'obsOnTimes', 'frameTimeStamps', 'frameTimeStampsWisk', 'pixelsPerM', ...
    'wheelCenter', 'wheelRadius', 'obsHeights', 'wiskContactTimes');


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


% determine trials to include
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'kinData.mat'), 'kinData')
vars = getExperimentData(session, ...
    {'modPawPredictedDistanceToObs', 'modPawDeltaLength', 'modPawX', 'isBigStep', 'isTrialSuccess', 'trialVel'});
vars = vars(1).data(1).sessions(1).trials;

% only include trials where (1) left paw in stance
%                           (2) right in swing
%                           (3) horizontal position is in middle percentiles at contact
%                           (4) trial is successful (not too much paw contact)
%                           (5) light is off ?
%                           (6) vel is decent
bins = false(1, length(kinData));
bins([kinData.isTrialAnalyzed]) = [kinData.isRightSwingAtContact]==1 & [kinData.isLeftSwingAtContact]==0;
xPosLims = prctile([vars.modPawX], [10 80]);
bins = bins & [vars.modPawX]>xPosLims(1) & [vars.modPawX]<xPosLims(2);
bins = bins & [vars.isTrialSuccess] & [vars.trialVel]>.5;
% bins = bins & ~isLightOn';
if ~any(bins); disp('no trials met criteria!'); return; end

deltas = [vars.modPawDeltaLength];
deltas(~bins) = nan;

% among these trials, pick most shortened, most lengthened, and least changed steps
% todo: the following fails when everything is masked out...
[a, trialShortened]  = min(deltas);                              % most shortened step

deltasTemp = deltas; deltasTemp([vars.isBigStep]==0) = nan;
[b, trialLengthened] = max(deltasTemp);                          % most lengthened step (excluding little step trials)

deltasTemp = deltas; deltasTemp([vars.isBigStep]==1) = nan;
[c, trialNoChange]   = min(abs(deltasTemp));                     % step with least change (excluding big step trials)



% set up figure
close all
fig = figure('position', [0, 600, range(s.xLims) * pixelsPerM, size(frame,1)], ...
    'name', session, 'color', [0 0 0], 'menubar', 'none', 'visible', s.visible);
ax = axes('position', [0 0 1 1], 'CLim', [0 255]);
colormap gray
im = image(frame, 'CDataMapping', 'scaled', 'ydata', frameY); hold on;
set(ax, 'visible', 'off', 'XLim', s.xLims); hold on


% text
if ~isempty(s.text); text(s.xLims(2), frameY(end), s.text, 'color', 'white', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom'); end


% kinematics
modPlot = plot(nan(2,2), nan(2,2), 'LineWidth', 3);  % array of two lines (for top and bottom views)
ctlPlot = plot(nan(2,2), nan(2,2), 'LineWidth', 3, 'Color', s.colors(3,:), 'LineStyle', ':');


% obstacle
obsTop = rectangle('position', [0-s.obsDiam/2, 0, s.obsDiam, s.obsDiam], 'curvature', [1 1], ...
    'facecolor', [s.obsColor .8], 'EdgeColor', 'none');
rectangle('position', [0-s.obsDiam/2, s.topHgt, s.obsDiam, frameY(end)-s.topHgt], ...
    'facecolor', [s.obsColor .8], 'EdgeColor', 'none');


% get average control step
% (average the final control step across all trials)
controlKinAvg = nan(length(kinData), 3, 500);
for i = find([kinData.isTrialAnalyzed])
    bins = kinData(i).controlStepIdentities(:,3)==max(kinData(i).controlStepIdentities(:,3));  % bins of final control step
    trialKin = kinData(i).locationsPix(bins,:,3);   % right forepaw control kinematics
    
    % unravel kinematics
    obsPixPositions = polyval(wheelToObsPixPosMappings(i,:), wheelPositions);
    obsPixPositions = interp1(wheelTimes, obsPixPositions, frameTimeStamps);
    trialKin(:,1) = trialKin(:,1) - obsPixPositions(kinData(i).trialInds(bins));
    
    trialKin(:,1) = trialKin(:,1) - trialKin(1,1);  % make x start at 0
    controlKinAvg(i,:,:) = interp2(1:3, (1:size(trialKin,1))', trialKin, 1:3, linspace(1,size(trialKin,1),500)')';  % upsample
end
controlKinAvg = squeeze(nanmean(controlKinAvg, 1)) / pixelsPerM;
% figure; plot(controlKinAvg(1,:), controlKinAvg(2,:)); hold on; plot(controlKinAvg(1,:), controlKinAvg(3,:)); daspect([1 1 1])


% edit video
trials = [trialShortened, trialLengthened, trialNoChange];
names = {'shortened', 'lengthened', 'noChange'};


vidWriter = VideoWriter(vidName);
set(vidWriter, 'FrameRate', fps)
if strcmp(vidSetting, 'MPEG-4'); set(vidWriter, 'Quality', 50); end
open(vidWriter)

for i = find(~isnan([a b c]))  % only use trial types for which there were valid trials to try
    
    % load obstacle positions for trial
    obsPixPositions = polyval(wheelToObsPixPosMappings(trials(i),:), wheelPositions);
    obsPixPositions = interp1(wheelTimes, obsPixPositions, frameTimeStamps);
    obsPositionsMeters = obsPixPositions / pixelsPerM;  % relative to left side of video frame
    
    % get trial kinematics
    trialLocations = kinData(trials(i)).locationsPix;
    trialLocations = squeeze(trialLocations(:,:,3)) / pixelsPerM;  % right forepaw only
    trialLocations(:,1) = trialLocations(:,1) - obsPositionsMeters(kinData(trials(i)).trialInds);
    trialLocations(kinData(trials(i)).modifiedStepIdentities(:,3)~=1,:) = nan;  % mask out all but modified step
    
    % obstacle
    wheelTopZ = (wheelCenter(2)-wheelRadius)/pixelsPerM;
    obsZ = wheelTopZ - obsHeights(trials(i))/1000;  % y axis is reversed
    set(obsTop, 'position', [0-s.obsDiam/2, obsZ, s.obsDiam, s.obsDiam]);
    
    % find trial inds
    trialInds = find(obsPositionsMeters<(-s.xLims(1)+frameX(end)) & obsPositionsMeters>-s.xLims(2));
    trialInds = trialInds(1) : s.dropFrames : trialInds(end);
    contactInd = find(frameTimeStamps>=wiskContactTimes(trials(i)), 1, 'first');  % whisker contact ind
    contactInd = trialInds(knnsearch(trialInds', contactInd));  % makes sure contactInd is still found when s.dropFrames>1
    
    trialLocations = trialLocations(ismember(kinData(trials(i)).trialInds, trialInds), :);  % trial locations sliced the same as trialInds
    firstModInd = find(~isnan(trialLocations(:,1)), 1, 'first');  % index at which to start tracing kinematics of first mod paw
    lastModInd = find(~isnan(trialLocations(:,1)), 1, 'last');  % index at which to start tracing kinematics of first mod paw
    
    % reset kinematics
    for j = 1:2
        set(modPlot(j), 'xdata', nan, 'ydata', nan, 'Color', s.colors(i,:));
        set(ctlPlot(j), 'xdata', nan, 'ydata', nan);
    end
    
    for j = 1:length(trialInds)
        
        % image
        if s.includeWiskCam
            frame = getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, trialInds(j), ...
                'yWiskPos', yWiskPos, 'xWiskPos', xWiskPos, 'wiskScaling', wiskScaling, ...
                'runContrast', s.contrastLims);
            frame = repmat(frame, 1, 1, 3);  % add color dimension
        else
            frame = read(vid, trialInds(j));
        end
        set(im, 'XData', frameX - obsPositionsMeters(trialInds(j)), 'CData', frame);
        
        % tracking
        if j>firstModInd && j<=lastModInd
            xyz = trialLocations(firstModInd:j,:);
            xyz = interp2(1:3, (firstModInd:j)', xyz, 1:3, linspace(firstModInd,j,500)', 'spline');
            set(modPlot(1), 'xdata', xyz(:,1), 'ydata', xyz(:,3))
            set(modPlot(2), 'xdata', xyz(:,1), 'ydata', xyz(:,2))
        end

        % write frame to video
        writeVideo(vidWriter, getframe(fig));
        
        % pause at whisker contact
        if trialInds(j)==contactInd

            % stretch s.t. predicted length is correct
            controlKin = controlKinAvg;
            predictedLength = kinData(trials(i)).modPredictedLengths(1,3);
            controlKin(1,:) = controlKin(1,:) * (predictedLength/controlKin(1,end));
            
            % offset horizontally and vertically to match the position at beginning of swing
            controlKin(1,:) = controlKin(1,:)                   + trialLocations(firstModInd,1);
            
            % find ind at which x pos of controlKin matches x pos of paw at contact
            xInd = knnsearch(controlKin(1,:)', trialLocations(j,1));
            
            % shift y and z coords to match this position
            controlKin(2,:) = controlKin(2,:) - controlKin(2,xInd) + trialLocations(j,2);
            controlKin(3,:) = controlKin(3,:) - controlKin(3,xInd) + trialLocations(j,3);
            
            % draw
            drawInds = ceil(linspace(xInd, size(controlKin,2), round(fps*s.contactPauseTime)));
            for k = drawInds
                set(ctlPlot(1), 'xdata', controlKin(1,xInd:k), 'ydata', controlKin(3,xInd:k))
                set(ctlPlot(2), 'xdata', controlKin(1,xInd:k), 'ydata', controlKin(2,xInd:k))
                writeVideo(vidWriter, getframe(fig));
            end
        end
    end

    % fade out last frame
    if blankFrames>0; for k = linspace(1,0,blankFrames); writeVideo(vidWriter, getframe(fig).cdata*k); end; end
end

fprintf('\nall done!\n')
close(vidWriter)
close(fig)


