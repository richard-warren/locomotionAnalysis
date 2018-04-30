function showTrackingDLC(session, vidDelay, trainingDataPath)

% settings
circSize = 150;
vidSizeScaling = 1.25;
colorMap = 'hsv';
connectedFeatures = {{'gen', 'tailBase', 'tailMid'}, {'tailBaseTop', 'tailMidTop'}}; % features that are connected within a view (not across views)



% initializations
frameInds = getFramesToShow(session);
addingFrames = exist('trainingDataPath', 'var');
if addingFrames; load(trainingDataPath, 'trainingData', 'view'); end

% get videos
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);

% get locations data and convert to 3d matrix
[locations, features, featurePairInds, isInterped] = fixTrackingDLC(session);

% set up figure
hgt = (vidBot.Height+vidTop.Height);
fig = figure('units', 'pixels', 'position', [600 400 vidBot.Width*vidSizeScaling hgt*vidSizeScaling],...
    'menubar', 'none', 'color', 'black', 'keypressfcn', @changeFrames);
colormap gray
imPreview = image(zeros(hgt, vidBot.Width), 'CDataMapping', 'scaled'); hold on;
imAxis = gca;
set(imAxis, 'visible', 'off', 'units', 'pixels',...
    'position', [0 0 vidBot.Width*vidSizeScaling hgt*vidSizeScaling]);

% set colors s.t. matching features in top and bot view have same color
cmap = eval(sprintf('%s(%i);', colorMap, length(features)));
for i = 1:size(featurePairInds,1)
    cmap(featurePairInds(i,2),:) = cmap(featurePairInds(i,1),:);
end


% set up lines joining same featres in top and bot
lines = cell(size(featurePairInds,1),1);
for i = 1:length(lines)
    lines{i} = line([0 0], [0 0], 'color', cmap(featurePairInds(i),:));
end

% set up lines joing features within a view
connectedFeatureInds = cell(1,length(connectedFeatures));
for i = 1:length(connectedFeatures)
    connectedFeatureInds{i} = nan(1,length(connectedFeatures{i}));
    for k = 1:length(connectedFeatures{i})
        connectedFeatureInds{i}(k) = find(ismember(features, connectedFeatures{i}(k)));
    end
end

for i = 1:length(connectedFeatures)
    linesConnected{i} = line([0 0], [0 0], 'color', [1 1 1]);
end

% set up scatter points
scatterLocations = scatter(imAxis, zeros(1,length(features)), zeros(1,length(features)),...
    circSize, cmap, 'linewidth', 3); hold on

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
    if frameInds(frameInd) > vidBot.NumberOfFrames; frameInd = 1;
    elseif frameInd < 1; frameInd = length(frameInds); end
    
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
    
    % update vertical lines
    for j = 1:length(lines)
        set(lines{j}, ...
            'xdata', [locations(frameInds(frameInd), 1, featurePairInds(j,1)) locations(frameInds(frameInd), 1, featurePairInds(j,2))], ...
            'ydata', [locations(frameInds(frameInd), 2, featurePairInds(j,1)) locations(frameInds(frameInd), 2, featurePairInds(j,2))])
    end
    
    % lines connecting within view features
    for j = 1:length(connectedFeatures)
        set(linesConnected{j}, 'xdata', locations(frameInds(frameInd),1,connectedFeatureInds{j}), ...
            'ydata', locations(frameInds(frameInd),2,connectedFeatureInds{j}));
    end

    % upate scatter positions
    set(scatterLocations, 'XData', locations(frameInds(frameInd),1,:), ...
        'YData', locations(frameInds(frameInd),2,:), ...
        'SizeData', ones(1,length(features))*circSize - (ones(1,length(features)).*isInterped(frameInds(frameInd),:))*circSize*.9);

    % pause to reflcet on the little things...
    pause(vidDelay);
end



end