function showTracking(session, trainingDataPath)

% to do: show confidence and gobal class in bottom left // think about how
% i should handle confidence thresholds...

% settings
onlyShowFramesNearObs = true;
vidFs = 250;
vidDelay = .02;
showDlcScores = false;
showTouchData = true;
showStance = true;
circSize = 100;
vidSizeScaling = 1.5;
colorMap = 'hsv';
scaling = 1.0; % network was trained on resolution of (saling)*(original resolution)
connectedFeatures = {{'paw1LH_bot', 'paw1LH_top'}, ...
                     {'paw2LF_bot', 'paw2LF_top'}, ...
                     {'paw3RF_bot', 'paw3RF_top'}, ...
                     {'paw4RH_bot', 'paw4RH_top'}, ...
                     {'tailBase_bot', 'tailMid_bot'}, ...
                     {'tailBase_top', 'tailMid_top'}, ...
                     {'obsHigh_bot', 'obsLow_bot'}}; % features that are connected within a view (not across views)



% initializations
frameInds = getFramesToShow(session, onlyShowFramesNearObs);
addingFrames = exist('trainingDataPath', 'var'); % keeps track of whether frames will be added to the training set
if addingFrames
    load(trainingDataPath, 'trainingData', 'view');
    framesAdded = 0;
end

% load video
vidName = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4');
if ~exist(vidName, 'file'); concatTopBotVids(session); end  % old sessions were recorded with separate top and bot views, which need to be concatenated
vid = VideoReader(vidName);

% get locations data and convert to 3d matrix
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'wheelPositions', 'wheelTimes', 'pixelsPerM', ...
    'wheelCenter', 'wheelRadius', 'touchesPerPaw', 'touchClassNames', 'touchConfidences', 'obsOnTimes', 'isLightOn');
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
[locations, features, ~, isInterped, scores] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM);
locations = locations / scaling; % bring back to original resolution
topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
if showStance
    stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, ...
        wheelTimes, wheelCenter, wheelRadius, vidFs, pixelsPerM);
end


% set up figure
if addingFrames; figureName = [session ', frames added: 0']; else; figureName = session; end
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
        
        % 'a': add frame to training set
        elseif key==97
            if addingFrames
                m = msgbox('adding frame to training set...'); pause(.5); close(m)
                ind = length(trainingData)+1;
                trainingData(ind).session = session;
                trainingData(ind).frameNum = frameInds(frameInd);
                trainingData(ind).includeFrame = false;
                for j = 1:length(features)
                    trainingData(ind).(features{j}) = squeeze(locations(frameInds(frameInd),:,j));
                end

                % resort the structure so like sessions stay together
                [~, sortInds] = sort({trainingData.session});
                trainingData = trainingData(sortInds);
                
                % update figure title
                framesAdded = framesAdded + 1;
                set(fig, 'name', sprintf('%s, frames added: %i', session, framesAdded))
            end
        
        % 's': save training set
        elseif key==115
            uisave({'trainingData', 'view'}, trainingDataPath)
        
        % 'f': select frame
        elseif key==102                  
            pause(.001);
            paused = true;
            input = inputdlg('enter frame number');
            frameInd = find(frameInds>=str2num(input{1}),1,'first');
            updateFrame(0);
            
        % 't': go to specific trial
        elseif key==116
            pause(.001);
            paused = true;
            input = inputdlg('enter trial number');
            frameInd = find(frameTimeStamps(frameInds)>=obsOnTimes(str2num(input{1})),1,'first');
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
    if frameInd < 1; frameInd = length(frameInds);
    elseif frameInd > length(frameInds); frameInd = 1; end
    
    
    % get frame and sub-frames
    frame = rgb2gray(read(vid, frameInds(frameInd)));
    
	% add frame number
    trial = find(obsOnTimes>=frameTimeStamps(frameInds(frameInd)),1,'first')-1;
    if isLightOn(trial); lightText = 'light on'; else; lightText = 'light off'; end
    frame = insertText(frame, [size(frame,2) size(frame,1)], ...
        sprintf('session %s, frame %i, trial %i, %s', ...
        session, frameInds(frameInd), trial, lightText), ...
        'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
    
    % update figure
    set(imPreview, 'CData', frame);
    
    
    % lines connecting within view features
    for j = 1:length(connectedFeatures)
        set(linesConnected{j}, 'xdata', locations(frameInds(frameInd),1,connectedFeatureInds{j}), ...
            'ydata', locations(frameInds(frameInd),2,connectedFeatureInds{j}));
    end

    % upate scatter positions
    set(scatterLocations, 'XData', locations(frameInds(frameInd),1,:), ...
        'YData', locations(frameInds(frameInd),2,:), ...
        'SizeData', ones(1,length(features))*circSize - (ones(1,length(features)) ...
                    .* isInterped(frameInds(frameInd),:)) * circSize * .9);
    
    % update scatter stance positions
    if showStance
        isStance = repmat(stanceBins(frameInds(frameInd),:),1,2);
        set(scatterStance, ...
            'XData', squeeze(locations(frameInds(frameInd),1,[botPawInds topPawInds])) .* isStance', ...
            'YData', squeeze(locations(frameInds(frameInd),2,[botPawInds topPawInds])));
    end
    
    % update paw touch scatter
    if exist('touchesPerPaw', 'var')
        touchingBins = touchesPerPaw(frameInds(frameInd),:)>0;
        x = locations(frameInds(frameInd), 1, [topPawInds(touchingBins) botPawInds(touchingBins)]);
        y = locations(frameInds(frameInd), 2, [topPawInds(touchingBins) botPawInds(touchingBins)]);
        set(obsTouchScatter, 'XData', x, 'YData', y);
    end

    % update scores text
    if showDlcScores
        for j = 1:length(features)
            set(scoreLabels{j}, 'position', [locations(frameInds(frameInd),1,j)+10, locations(frameInds(frameInd),2,j)], ...
                'string', sprintf('%.2f', scores(frameInds(frameInd),j)));
        end
    end
    
    if showTouchData
        for j = 1:4
            classInd = touchesPerPaw(frameInds(frameInd),j);
            if classInd==0; classInd=find(strcmp(touchClassNames, 'no_touch')); end
            class = touchClassNames{classInd};
            confidence = touchConfidences(frameInds(frameInd));
            set(touchScoreLabels{j}, ...
                'position', [locations(frameInds(frameInd),1,topPawInds(j))+10, locations(frameInds(frameInd),2,topPawInds(j))], ...
                'string', sprintf('%s (%.2f)', class, confidence));
            
        end
    end

    % pause to reflcet on the little things...
    pause(vidDelay);
end



end