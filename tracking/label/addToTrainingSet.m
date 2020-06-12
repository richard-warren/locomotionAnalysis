function addToTrainingSet(session, vidName, trackedFeaturesName, trainingSet, varargin)

% show tracking results and add badly tracked frames to a training set

% settings
s.skeleton = [];
s.vidFs = 250;
s.vidDelay = .02;
s.showLabels = true;
s.circSize = 100;
s.vidScaling = 1.5;
s.lineColor = [.15 .15 1];
s.scoreThresh = .5;  % don't show locations when confidence is below scoreThresh
s.invertFrame = false;
s.textRotation = 45;
s.deltaThresh = 50;  % (pixels) // when looking for frames where paws in top are bottom views are misaligned, only include frames where the x difference is > deltaThresh



% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
load(trainingSet, 'trainingData');
framesAdded = 0;

% load data
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, vidName));
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, trackedFeaturesName)); % get raw tracking data
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'rewardTimes', 'frameTimeStamps', 'frameTimeStampsWisk', 'wiskContactTimes'); % get raw tracking data

% a hack to determine with frameTimes to use
if contains(vidName, 'isk'); times = frameTimeStampsWisk; else; times = frameTimeStamps; end

features = locationsTable.Properties.VariableNames(2:3:end);
scores = table2array(locationsTable(:,4:3:end));
locations = nan(height(locationsTable), 2, length(features));
for i = 1:length(features)
    locations(:,1,i) = table2array(locationsTable(:, (i-1)*3 + 2));
    locations(:,2,i) = table2array(locationsTable(:, (i-1)*3 + 3));
end
colors = hsv(length(features));

% find frames where x positions are misaligned
pawPairs = {{'paw1LH_top', 'paw1LH_bot'}, {'paw2LF_top', 'paw2LF_bot'}, {'paw3RF_top', 'paw3RF_bot'}, {'paw4RH_top', 'paw4RH_bot'}};
if any(ismember(features, cat(2,pawPairs{:})))  % if paws are actually in this video
    deltaX = nan(size(locations,1), 4);
    for i = 1:length(pawPairs)
        bins = ismember(features, pawPairs{i});
        deltaX(:,i) = diff(squeeze(locations(:,1,bins)),[],2);
    end
    isMisaligned = all(scores(:,contains(features, 'paw'))>s.scoreThresh,2) & ...  % all paws have high confidence
        any(deltaX>s.deltaThresh,2);  % at least one paw is misaligned
    isMisaligned = [false; diff(isMisaligned)==1];  % only keep first true is sequences of consecutive trues
    misalignedInds = find(isMisaligned);
else
    deltaX = zeros(size(locations,1), 4);
end


% set up figure
fig = figure('name', [session ', frames added: 0'], 'units', 'pixels', ...
    'position', [600 100 vid.Width*s.vidScaling vid.Height*s.vidScaling], ...
    'menubar', 'none', 'color', 'black', 'keypressfcn', @changeFrames);
colormap gray
[frame, frameRaw] = deal(zeros(vid.Height, vid.Width));
imPreview = image(frame, 'CDataMapping', 'scaled'); hold on;
imAxis = gca;
set(imAxis, 'visible', 'off', 'units', 'pixels',...
    'position', [0 0 vid.Width*s.vidScaling vid.Height*s.vidScaling]);


% set up scatter points for tracked features
scatterLocations = scatter(imAxis, zeros(1,length(features)), zeros(1,length(features)),...
    s.circSize, colors, 'linewidth', 3); hold on


% set up lines connecting features
if ~isempty(s.skeleton)
    skeletonTbl = readtable(s.skeleton);
    hasParentInds = find(~cellfun(@isempty, skeletonTbl.parent));  % bins of features with a parent feature, as determined by spreadsheet
    numEdges = length(hasParentInds);
    lineColors = hsv(numEdges);
    edges = nan(numEdges, 2);  % number of edges X 2 matrix storing connected features
    featurePairLines = cell(1,numEdges);

    lines = cell(1,numEdges);
    for i = 1:numEdges
        edges(i,1) = hasParentInds(i);
        edges(i,2) = find(strcmp(skeletonTbl.name, skeletonTbl.parent(hasParentInds(i))));  % index of parent feature
        featurePairLines{i} = line([0 0], [0 0], 'color', lineColors(i,:)); hold on
        uistack(featurePairLines{i}, 'bottom');
        uistack(featurePairLines{i}, 'up');
        
        lines{i} = line(0, 0, 'lineWidth', 2, 'color', lineColors(i,:));
    end
end


% set up text to show dlc scores
scoreLabels = cell(1,length(features));
if s.showLabels
    for i = 1:length(features); scoreLabels{i} = text(0,0,'', 'color', colors(i,:), 'Rotation', s.textRotation, 'interpreter', 'none'); end
end


% state variables
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
            trainingData(ind).frameNum = frameInd;
            trainingData(ind).frame = frameRaw;
            trainingData(ind).includeFrame = false;
            for j = 1:length(features)
                if scores(frameInd,j) > s.scoreThresh
                    trainingData(ind).(features{j}) = squeeze(locations(frameInd,:,j));
                else
                    trainingData(ind).(features{j}) = [nan nan];
                end
            end
            fprintf('added frame %i to trainingData\n', ind)

            % resort the structure so like sessions stay together
            [~, sortInds] = sort({trainingData.session});
            trainingData = trainingData(sortInds);

            % update figure title
            framesAdded = framesAdded + 1;
            set(fig, 'name', sprintf('%s, frames added: %i', session, framesAdded))
        
        % 't': go to specific time (in seconds)
        elseif key==116
            pause(.001);
            paused = true;
            input = inputdlg('enter time (seconds)');
            frameInd = find(times>=str2num(input{1}),1,'first');
            updateFrame(0);
            
        % 's': save training set
        elseif key==115
            uisave({'trainingData'}, trainingSet)
        
        % 'm': go to next misaligned frame
        elseif key==109
            paused = true;
            nextMisaligned = find(misalignedInds>frameInd,1,'first');
            if ~isempty(nextMisaligned)
                frameInd = misalignedInds(nextMisaligned);
            else
                disp('no more misaligned frames!')
            end
            updateFrame(0);
        
        % 'f': select frame
        elseif key==102                  
            pause(.001);
            paused = true;
            input = inputdlg('enter frame number');
            frameInd = str2num(input{1});
            updateFrame(0);
            
        % ESCAPE: close window
        elseif key==27                  
            playing = false;
            paused = false;
            
        % 'w': go to next water drop
        elseif key==119
            nextRewardTime = rewardTimes(find(rewardTimes>times(frameInd), 1, 'first'));
            frameInd = find(times>nextRewardTime, 1, 'first');
            updateFrame(0);
        
        % 'o': go to next whisker contact
        elseif key==111
            nextContactTime = wiskContactTimes(find(wiskContactTimes>times(frameInd), 1, 'first'));
            frameInd = find(times>nextContactTime, 1, 'first');
            updateFrame(0);
        
        % OTHERWISE: toggle pausing
        else                            
            paused = ~paused;
        end
    end
end



% update frame preview
function updateFrame(frameStep)
    
    frameInd = frameInd + frameStep;
    if frameInd < 1; frameInd = length(times);
    elseif frameInd > length(times); frameInd = 1; end
    
    % get frame
    frameRaw = rgb2gray(read(vid, frameInd));
    frame = frameRaw;
    if s.invertFrame; frame = 255-frame; end
    
	% add frame number
    frame = insertText(frame, [size(frame,2) size(frame,1)], ...
        sprintf('session %s, f=%i, t=%.2f', ...
        session, frameInd, times(frameInd)), ...
        'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
    
    % update figure
    set(imPreview, 'CData', frame);
    
    % determine which features to show
    showFeature = scores(frameInd,:) > s.scoreThresh;
    
    % upate scatter positions
    x = locations(frameInd,1,:);
    y = locations(frameInd,2,:);
    x(~showFeature) = nan;
    y(~showFeature) = nan;
    set(scatterLocations, 'XData', x, 'YData', y);
    
    % update scores text
    if s.showLabels
        for j = 1:length(features)
            set(scoreLabels{j}, 'position', [x(j)+10, y(j)], ...
                'string', sprintf('%s: %.2f', features{j}, scores(frameInd,j)));
        end
    end
    
    % update lines
    if exist('edges', 'var')
        for j = 1:numEdges
            set(lines{j}, 'XData', x(edges(j,:)), 'YData', y(edges(j,:)));
        end
    end
    
    % pause to reflcet on the little things...
    pause(s.vidDelay);
end



end