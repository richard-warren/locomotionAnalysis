function showTrackingDLC(session, vidDelay, trainingDataPath)

% settings
onlyShowFramesNearObs = true;
vidFs = 250;
showScores = true;
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

% get videos
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);

% get locations data and convert to 3d matrix
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
    'frameTimeStamps', 'wheelPositions', 'wheelTimes', 'mToPixMapping', 'wheelCenter', 'wheelRadius', 'arePawsTouchingObs');
locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' session '\trackedFeaturesRaw.csv']); % get raw tracking data
[locations, features, ~, isInterped, scores] = fixTrackingDLC(locationsTable, frameTimeStamps);
locations = locations / scaling; % bring back to original resolution
topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
if showStance
    stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, ...
        wheelTimes, wheelCenter, wheelRadius, vidFs, mToPixMapping(1));
end


% set up figure
if addingFrames; figureName = [session ', frames added: 0']; else; figureName = session; end
hgt = (vidBot.Height+vidTop.Height);
fig = figure('name', figureName, 'units', 'pixels', 'position', [600 400 vidBot.Width*vidSizeScaling hgt*vidSizeScaling],...
    'menubar', 'none', 'color', 'black', 'keypressfcn', @changeFrames);
colormap gray
imPreview = image(zeros(hgt, vidBot.Width), 'CDataMapping', 'scaled'); hold on;
imAxis = gca;
set(imAxis, 'visible', 'off', 'units', 'pixels',...
    'position', [0 0 vidBot.Width*vidSizeScaling hgt*vidSizeScaling]);

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

% set up text to show scores
scoreLabels = cell(1,length(features));
if showScores
    for i = 1:length(features)
        scoreLabels{i} = text(0,0,'', 'color', cmap(i,:));
    end
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
            frameInd = find(frameInds==str2num(input{1}));
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
    frameBot = rgb2gray(read(vidBot, frameInds(frameInd)));
    frameTop = rgb2gray(read(vidTop, frameInds(frameInd)));
    frame = cat(1, frameTop, frameBot);
    
	% add frame number
    frame = insertText(frame, [size(frame,2) size(frame,1)], ...
        sprintf('session %s, frame %i', session, frameInds(frameInd)),...
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
    if exist('arePawsTouchingObs', 'var')
        touchingBins = arePawsTouchingObs(frameInds(frameInd),:);
        x = locations(frameInds(frameInd), 1, [topPawInds(touchingBins) botPawInds(touchingBins)]);
        y = locations(frameInds(frameInd), 2, [topPawInds(touchingBins) botPawInds(touchingBins)]);
        set(obsTouchScatter, 'XData', x, 'YData', y);
    end

    % update scores text
    if showScores
        for j = 1:length(features)
            set(scoreLabels{j}, 'position', [locations(frameInds(frameInd),1,j)+10, locations(frameInds(frameInd),2,j)], ...
                'string', sprintf('%.2f', scores(frameInds(frameInd),j)));
        end
    end

    % pause to reflcet on the little things...
    pause(vidDelay);
end



end