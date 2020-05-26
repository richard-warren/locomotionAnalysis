function showTracking(session, varargin)

% to do: show confidence and gobal class in bottom left // think about how
% i should handle confidence thresholds...

% additional features: paw touches // arbitrary signal // arbitrary times
% // wheel // post-processed tracking // whisker overlay // whisker
% contacts // confidences // go to water, whisker contact

% settings
s.scoreThresh = .8;
vidFs = 250;
vidDelay = .02;
showDlcScores = false;
showTouchData = true;
showStance = true;
circSize = 100;
vidSizeScaling = 1.5;
colorMap = 'hsv';
connectedFeatures = {{'paw1LH_bot', 'paw1LH_top'}, ...
                     {'paw2LF_bot', 'paw2LF_top'}, ...
                     {'paw3RF_bot', 'paw3RF_top'}, ...
                     {'paw4RH_bot', 'paw4RH_top'}, ...
                     {'tailBase_bot', 'tailMid_bot'}, ...
                     {'tailBase_top', 'tailMid_top'}, ...
                     {'obsHigh_bot', 'obsLow_bot'}}; % features that are connected within a view (not across views)



if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

% load video
vidName = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4');
if ~exist(vidName, 'file'); concatTopBotVids(session); end  % old sessions were recorded with separate top and bot views, which need to be concatenated
vid = VideoReader(vidName);

% get locations data and convert to 3d matrix
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'wheelPositions', 'wheelTimes', 'pixelsPerM', ...
    'wheelCenter', 'wheelRadius', 'touchesPerPaw', 'touchClassNames', 'touchConfidences', 'obsOnTimes', 'isLightOn');
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
[locations, features, ~, isInterped, scores] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM, 'scoreThresh', s.scoreThresh);
topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
if showStance
    stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, ...
        wheelTimes, wheelCenter, wheelRadius, vidFs, pixelsPerM);
end


% set up figure
figureName = session;
hgt = vid.Height;
fig = figure('name', figureName, 'units', 'pixels', 'position', [600 100 vid.Width*vidSizeScaling hgt*vidSizeScaling],...
    'menubar', 'none', 'color', 'black', 'keypressfcn', @changeFrames);
colormap gray
imPreview = image(zeros(hgt, vid.Width), 'CDataMapping', 'scaled'); hold on;
imAxis = gca;
set(imAxis, 'visible', 'off', 'units', 'pixels',...
    'position', [0 0 vid.Width*vidSizeScaling hgt*vidSizeScaling]);

% draw circle at wheel location
viscircles(wheelCenter', wheelRadius, 'color', 'blue');

% set colors s.t. matching features in top and bot view have same color
cmap = eval(sprintf('%s(%i);', colorMap, length(features)));


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
    circSize, cmap, 'linewidth', 3); hold on

% set up scatter points that will surround paw when it is touching obs
obsTouchScatter = scatter(imAxis, [], [], circSize*3, [1 1 1], 'LineWidth', 2);

% set up stance scatter points
if showStance
    scatterStance = scatter(imAxis, ...
        zeros(1,length([botPawInds topPawInds])), zeros(1,length([botPawInds topPawInds])), ...
        circSize, cmap([botPawInds topPawInds],:), 'filled'); hold on
end

% set up text to show dlc scores
scoreLabels = cell(1,length(features));
if showDlcScores
    for i = 1:length(features); scoreLabels{i} = text(0,0,'', 'color', cmap(i,:)); end
end

% set up text to show touch tracking info
touchScoreLabels = cell(1,length(features));
if showTouchData
    for i = 1:4; touchScoreLabels{i} = text(0,0,'', 'color', [1 1 1], 'interpreter', 'none'); end
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
    frame = rgb2gray(read(vid, frameInd));
    
	% add frame number
    trial = find(obsOnTimes>=frameTimeStamps(frameInd),1,'first')-1;
    if trial
        if isLightOn(trial); lightText = 'light on'; else; lightText = 'light off'; end
    else
        lightText = '';
    end
    frame = insertText(frame, [size(frame,2) size(frame,1)], ...
        sprintf('frame %i, trial %i, %s', ...
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
        'SizeData', ones(1,length(features))*circSize - (ones(1,length(features)) ...
                    .* isInterped(frameInd,:)) * circSize * .9);
    
    % update scatter stance positions
    if showStance
        isStance = repmat(stanceBins(frameInd,:),1,2);
        set(scatterStance, ...
            'XData', squeeze(locations(frameInd,1,[botPawInds topPawInds])) .* isStance', ...
            'YData', squeeze(locations(frameInd,2,[botPawInds topPawInds])));
    end
    
    % update paw touch scatter
    if exist('touchesPerPaw', 'var')
        touchingBins = touchesPerPaw(frameInd,:)>0;
        x = locations(frameInd, 1, [topPawInds(touchingBins) botPawInds(touchingBins)]);
        y = locations(frameInd, 2, [topPawInds(touchingBins) botPawInds(touchingBins)]);
        set(obsTouchScatter, 'XData', x, 'YData', y);
    end

    % update scores text
    if showDlcScores
        for j = 1:length(features)
            set(scoreLabels{j}, 'position', [locations(frameInd,1,j)+10, locations(frameInd,2,j)], ...
                'string', sprintf('%.2f', scores(frameInd,j)));
        end
    end
    
    if showTouchData
        for j = 1:4
            classInd = touchesPerPaw(frameInd,j);
            if classInd==0; classInd=find(strcmp(touchClassNames, 'no_touch')); end
            class = touchClassNames{classInd};
            confidence = touchConfidences(frameInd);
            set(touchScoreLabels{j}, ...
                'position', [locations(frameInd,1,topPawInds(j))+10, locations(frameInd,2,topPawInds(j))], ...
                'string', sprintf('%s (%.2f)', class, confidence));
            
        end
    end

    % pause to reflcet on the little things...
    pause(vidDelay);
end



end