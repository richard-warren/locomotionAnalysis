function showTracking(session, varargin)

% to do: show confidence and gobal class in bottom left // think about how
% i should handle confidence thresholds...

% todo: whisker angle // arbitrary signal // arbitrary times? // whisker contacts // import skeleton?
% events: whisker contact // paw contact // licks // spike times

% settings
s.includeWiskCam = true;  % whether to add whisker camera
s.scoreThresh = .8;
s.vidDelay = .02;
s.showConfidence = false;
s.showPawTouch = true;
s.showPawTouchConfidence = true;
s.showStance = true;
s.circSize = 100;
s.vidScaling = 1.5;
s.colorMap = 'hsv';
s.faceColor = 'yellow';  % color for face tracking points
s.frameProps = {'isPaddingWhite', false, 'edgeFading', 5, 'border', 2};  % arguments to be passed to getFrameWithWisk (only when includeWiskCam=true)
connectedFeatures = {{'paw1LH_bot', 'paw1LH_top'}, ...
                     {'paw2LF_bot', 'paw2LF_top'}, ...
                     {'paw3RF_bot', 'paw3RF_top'}, ...
                     {'paw4RH_bot', 'paw4RH_top'}, ...
                     {'tailBase_bot', 'tailMid_bot'}, ...
                     {'tailBase_top', 'tailMid_top'}, ...
                     {'obsHigh_bot', 'obsLow_bot'}}; % features that are connected within a view (not across views)


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

% load video
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4'));
if s.includeWiskCam; vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4')); end

% get locations data and convert to 3d matrix
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'frameTimeStampsWisk', 'wheelPositions', 'wheelTimes', 'pixelsPerM', 'wiskContactTimes', 'rewardTimes', ...
    'wheelCenter', 'wheelRadius', 'touchesPerPaw', 'touchClassNames', 'touchConfidences', 'obsOnTimes', 'isLightOn', 'whiskerAngle');
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
[locations, features, ~, isInterped, scores] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM, 'scoreThresh', s.scoreThresh);
topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
if s.showStance
    stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, ...
        wheelTimes, wheelCenter, wheelRadius, 1/median(diff(frameTimeStamps)), pixelsPerM);
end


% get position where wisk frame should overlap with runTop frame
if s.includeWiskCam
    [frame, yWiskPos, xWiskPos, wiskScaling] = ...
        getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, find(frameTimeStamps>obsOnTimes(1), 1, 'first'), s.frameProps{:});  % use first frame where obstacle is on to ensure mouse is on the wheel when the whisker cam position is determined
else
    frame = read(vid, 1);
end
sz = [size(frame,1) size(frame,2) 3];


% set up figure
fig = figure('name', session, 'units', 'pixels', 'position', [600 100 sz(1)*s.vidScaling sz(2)*s.vidScaling],...
    'menubar', 'none', 'color', 'black', 'keypressfcn', @changeFrames);
colormap gray
imPreview = image(frame, 'CDataMapping', 'scaled'); hold on;
imAxis = gca;
set(imAxis, 'visible', 'off', 'units', 'pixels', 'units', 'normalized', 'position', [0 0 1 1]);


% draw circle at wheel location
viscircles(wheelCenter', wheelRadius, 'color', 'blue');


% add whisker features
if s.includeWiskCam
    wiskFeatures = {'jaw', 'tongue', 'nose'};
    
    % remove locations from run tracking that will be in wisk tracking
    locations(:,:,contains(features, wiskFeatures)) = nan;
    
    % load whisker tracking
    locationsWiskTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw_wisk.csv')); % get raw tracking data
    locationsWisk = nan(height(locationsWiskTable), 2, length(wiskFeatures));
    
    % find whisker pad and length
    pad = [median(locationsWiskTable.wisk_pad) median(locationsWiskTable.wisk_pad_1)];  % location of whisker pad
    pad = pad*wiskScaling + [xWiskPos yWiskPos];
    tipAvg = [median(locationsWiskTable.wisk_caudal) median(locationsWiskTable.wisk_caudal_1)]*wiskScaling + [xWiskPos yWiskPos];
    wiskLength = norm(tipAvg - pad);
    
    % fix tracking
    for i = 1:length(wiskFeatures)
        valid = locationsWiskTable.([wiskFeatures{i} '_2']) > s.scoreThresh;
        
        locationsWisk(valid,1,i) = locationsWiskTable.(wiskFeatures{i})(valid)*wiskScaling + xWiskPos;
        locationsWisk(valid,2,i) = locationsWiskTable.([wiskFeatures{i} '_1'])(valid)*wiskScaling + yWiskPos;
    end
    
    clear locationsWiskTable
    scatterWisk = scatter(imAxis, zeros(1,length(wiskFeatures)), zeros(1,length(wiskFeatures)), 20, s.faceColor, 'filled');
    
    % create whisker angle line
    wisk = plot([0 0], [0 0], 'color', s.faceColor, 'LineWidth', 2);
end


% set colors s.t. matching features in top and bot view have same color
cmap = eval(sprintf('%s(%i);', s.colorMap, length(features)));


% set up lines joining features within a view
connectedFeatureInds = cell(1,length(connectedFeatures));
linesConnected = cell(1,length(connectedFeatures));
for i = 1:length(connectedFeatures)
    connectedFeatureInds{i} = nan(1,2);
    for k = 1:length(connectedFeatures{i})
        connectedFeatureInds{i}(k) = find(ismember(features, connectedFeatures{i}(k)));
    end
    linesConnected{i} = line([0 0], [0 0], 'color', 'white');
end


% set up scatter points for tracked features
scatterLocations = scatter(imAxis, zeros(1,length(features)), zeros(1,length(features)),...
    s.circSize, cmap, 'linewidth', 3); hold on


% set up scatter points that will surround paw when it is touching obs
if s.showPawTouch
    obsTouchScatter = scatter(imAxis, [], [], s.circSize*3, [1 1 1], 'LineWidth', 2);
end


% set up stance scatter points
if s.showStance
    scatterStance = scatter(imAxis, ...
        zeros(1,length([botPawInds topPawInds])), zeros(1,length([botPawInds topPawInds])), ...
        s.circSize, cmap([botPawInds topPawInds],:), 'filled'); hold on
end


% set up text to show confidence
if s.showConfidence
    confidenceLabels = cell(1,length(features));
    for i = 1:length(features); confidenceLabels{i} = text(0,0,'', 'color', cmap(i,:)); end
end


% set up text to show touch tracking info
if s.showPawTouchConfidence
    touchConfidenceLabels = cell(1,length(features));
    for i = 1:4; touchConfidenceLabels{i} = text(0,0,'', 'color', [1 1 1], 'interpreter', 'none'); end
end

% set state variables
frameInd = 1;
playing = true;
paused = false;


% main loop
while playing
    while paused; pause(.001); end
    updateFrame(1);
end
close(fig)






% keypress controls
function changeFrames(~,~)
    
    key = double(get(fig, 'currentcharacter'));
    
    if ~isempty(key) && isnumeric(key)
        
        % left: move frame backward
        if key==28                      
            pause(.001);
            paused = true;
            updateFrame(-1);
        
        % right: move frame forward
        elseif key==29                  
            pause(.001);
            paused = true;
            updateFrame(1);
        
        % 'f': select frame
        elseif key==102                  
            pause(.001);
            paused = true;
            input = inputdlg('enter frame number');
            frameInd = str2num(input{1});
            updateFrame(0);
            
        % 't': go to specific trial
        elseif key==116
            pause(.001);
            paused = true;
            input = inputdlg('enter trial number');
            frameInd = find(frameTimeStamps>=obsOnTimes(str2num(input{1})),1,'first');
            updateFrame(0);
        
        % 'w': go to next water drop
        elseif key==119
            nextRewardTime = rewardTimes(find(rewardTimes>frameTimeStamps(frameInd), 1, 'first'));
            frameInd = find(frameTimeStamps>nextRewardTime, 1, 'first');
            updateFrame(0);
        
        % 'o': go to next whisker contact
        elseif key==111
            nextContactTime = wiskContactTimes(find(wiskContactTimes>frameTimeStamps(frameInd), 1, 'first'));
            frameInd = find(frameTimeStamps>nextContactTime, 1, 'first');
            updateFrame(0);
            
        % ESCAPE: close window
        elseif key==27                  
            playing = false;
            paused = false;
        
        % OTHERWISE: toggle pausing
        else                            
            paused = ~paused;
        end
    end
end



% update frame preview
function updateFrame(frameStep)
    
    frameInd = frameInd + frameStep;
    if frameInd < 1; frameInd = vid.NumberOfFrames;
    elseif frameInd > vid.NumberOfFrames; frameInd = 1; end
    
    % get frame and sub-frames
    if s.includeWiskCam
        frame = getFrameWithWisk(vid, vidWisk, frameTimeStamps, frameTimeStampsWisk, frameInd, ...
            'yWiskPos', yWiskPos, 'xWiskPos', xWiskPos, 'wiskScaling', wiskScaling, s.frameProps{:});
        frame = repmat(frame, 1, 1, 3);  % add color dimension
    else
        frame = read(vid, frameInd);
    end
    
	% add metadata
    trial = find(obsOnTimes>=frameTimeStamps(frameInd),1,'first')-1;
    if trial
        if isLightOn(trial); lightText = 'light on'; else; lightText = 'light off'; end
    else
        lightText = '';
    end
    frame = insertText(frame, [size(frame,2) size(frame,1)], ...
        sprintf('frame %i trial %i %s', ...
        frameInd, trial, lightText), ...
        'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
    
    % update figure
    set(imPreview, 'CData', frame);
    
    
    % lines connecting within view features
    for j = 1:length(connectedFeatures)
        set(linesConnected{j}, 'xdata', locations(frameInd,1,connectedFeatureInds{j}), ...
            'ydata', locations(frameInd,2,connectedFeatureInds{j}));
    end

    % upate scatter positions
    set(scatterLocations, 'XData', locations(frameInd,1,:), ...
        'YData', locations(frameInd,2,:), ...
        'SizeData', ones(1,length(features))*s.circSize - (ones(1,length(features)) ...
                    .* isInterped(frameInd,:)) * s.circSize * .9);
    
    % update scatter stance positions
    if s.showStance
        isStance = repmat(stanceBins(frameInd,:),1,2);
        set(scatterStance, ...
            'XData', squeeze(locations(frameInd,1,[botPawInds topPawInds])) .* isStance', ...
            'YData', squeeze(locations(frameInd,2,[botPawInds topPawInds])));
    end
    
    % update paw touch scatter
    if s.showPawTouch
        touchingBins = touchesPerPaw(frameInd,:)>0;
        x = locations(frameInd, 1, [topPawInds(touchingBins) botPawInds(touchingBins)]);
        y = locations(frameInd, 2, [topPawInds(touchingBins) botPawInds(touchingBins)]);
        set(obsTouchScatter, 'XData', x, 'YData', y);
    end

    % update scores text
    if s.showConfidence
        for j = 1:length(features)
            set(confidenceLabels{j}, 'position', [locations(frameInd,1,j)+10, locations(frameInd,2,j)], ...
                'string', sprintf('%.2f', scores(frameInd,j)));
        end
    end
    
    % update paw touch confidence text
    if s.showPawTouchConfidence
        for j = 1:4
            classInd = touchesPerPaw(frameInd,j);
            if classInd==0; classInd=find(strcmp(touchClassNames, 'no_touch')); end
            class = touchClassNames{classInd};
            confidence = touchConfidences(frameInd);
            set(touchConfidenceLabels{j}, ...
                'position', [locations(frameInd,1,topPawInds(j))+10, locations(frameInd,2,topPawInds(j))], ...
                'string', sprintf('%s (%.2f)', class, confidence));
            
        end
    end
    
    % update whisker view tracking
    if s.includeWiskCam
        set(scatterWisk, 'XData', locationsWisk(frameInd,1,:), 'YData', locationsWisk(frameInd,2,:));
        
        % whisker angle
        x = cosd(whiskerAngle(frameInd)) * wiskLength;
        y = -sind(whiskerAngle(frameInd)) * wiskLength;
        set(wisk, 'XData', [pad(1) pad(1)+x], 'YData', [pad(2) pad(2)+y])
    end

    % pause to reflcet on the little things...
    pause(s.vidDelay);
end



end