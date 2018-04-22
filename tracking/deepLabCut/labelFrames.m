function labelFrames(trainingSetDir, features)

% add way of excluding trials

% settings
vidScaling = 2;
textOffset = [5 0];
nanColor = [1 0 0];
labelledColor = [1 1 0];

% initializations
defaultPositions = repmat([10 10], length(features), 1) .* repmat(1:length(features),2,1)';
load([trainingSetDir 'trainingData.mat'], 'trainingData');
fields = fieldnames(trainingData);
structInd = 1;
selectedPoint = nan;



% create figure
currentSession = trainingData(structInd).session;
currentFrame = trainingData(structInd).frameNum;
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' currentSession '\runBot.mp4']);
fig = figure('name', sprintf('%s, frame %i', currentSession, currentFrame), ...
    'menubar', 'none', 'color', 'white', 'keypressfcn', @keypress, 'position', [200 200 vid.Width*vidScaling, vid.Height*vidScaling]);
colormap gray
im = image(rgb2gray(read(vid, currentFrame)), 'CDataMapping', 'scaled'); hold on
set(gca, 'position', [0 0 1 1], 'visible', 'off')
selectedCirc = scatter(0, 0, 200, [.5 .5 1], 'linewidth', 3, 'visible', 'off'); hold on



% add fields that have not yet been created, and creatable draggable objects
featurePoints = cell(1, length(features));
featureTexts = cell(1, length(features));


for i = 1:length(features)
    
    % initialize non-existent features
    if ~ismember(features{i}, fields)
        nanEntries = mat2cell(nan(length(trainingData),2), ones(1,length(trainingData)), 2);
        [trainingData.(features{i})] = nanEntries{:};
        fprintf('creating field: %s\n', features{i});
    end
    
    % create draggables for features
    featureTexts{i} = text(defaultPositions(i,1)+textOffset(1), defaultPositions(i,2)+textOffset(2), features{i});
    featurePoints{i} = impoint(gca, defaultPositions(i,:));
    addNewPositionCallback(featurePoints{i}, @(x) movePoint(i));
    setPositionConstraintFcn(featurePoints{i}, @constrainPosition)
    
    % add special constraint for end of tail
    if strcmp(features{i}, 'tailEnd')
        setPositionConstraintFcn(featurePoints{i}, @tailEndConstrainPosition)
    end
end
updateFrame(0); % set positions and colors of points based on whether the first frame has already been labeled


% ---------
% FUNCTIONS
% ---------


% respond to key press
function keypress(~,~)
        
    key = double(get(fig, 'currentcharacter'));

    if ~isempty(key) && isnumeric(key)
        switch key
            % LEFT: move frame backward
            case 28
                updateFrame(-1);

            % RIGHT: move frame forward
            case 29
                updateFrame(1);
        

            % ESCAPE: close window
            case 27
                close(fig)
                
            % d: delete selectetion
            case 100
                if ~isnan(selectedPoint)
                    trainingData(structInd).(features{selectedPoint}) = [nan nan];
                    updateFrame(0);
                end
                
            % n: move to next frame with NaN value
            case 110
                isntLabeled = false(length(trainingData), length(features));
                for j = 1:length(features)
                    x = cellfun(@(k) k(1), {trainingData.(features{j})});
                    nanBins = isnan(x);
                    isntLabeled(nanBins, j) = 1;
                end
                notFullyLabelledBins = any(isntLabeled,2);
                
                newStructInd = find(notFullyLabelledBins' & 1:length(trainingData)>structInd, 1, 'first'); % find first frame with nonlabelled part that is later than current frame
                if isempty(newStructInd); newStructInd = find(notFullyLabelledBins, 1, 'first'); end % if couldnt find any, search for nonlabelled part starting from beginning
                if ~isempty(newStructInd); structInd = newStructInd; end % otherwise, keep current frame
                    
                updateFrame(0);
    
                
            % i: include frame for analysis
            case 105
                trainingData(structInd).includeFrame = ~trainingData(structInd).includeFrame;
                updateFrame(0);
                
            % s: save progress
            case 115
                save([trainingSetDir 'trainingData.mat'], 'trainingData');
                m = msgbox('saving progress...'); pause(.5); close(m)
        end
    end
end




function updateFrame(frameStep)
    
    % move to next frame
    structInd = structInd + frameStep;
    if structInd==0; structInd=length(trainingData); end
    if structInd>length(trainingData); structInd=1; end
    currentFrame = trainingData(structInd).frameNum;
    
    % load new video if switching to new session
    if ~strcmp(currentSession, trainingData(structInd).session)
        currentSession = trainingData(structInd).session;
        m = msgbox(sprintf('loading session %s...', currentSession)); pause(.5); close(m);
        vid = VideoReader([getenv('OBSDATADIR') 'sessions\' currentSession '\runBot.mp4']);
    end
    
    % update frame
    frame = rgb2gray(read(vid, trainingData(structInd).frameNum));
    set(im, 'CData', frame);
    if trainingData(structInd).includeFrame; includedString = '(included)'; else; includedString = ''; end
    set(fig, 'name', sprintf('%s, frame %i %s', currentSession, currentFrame, includedString))
    
    % update colors and positions of points
    for j = 1:length(features)
        pos = trainingData(structInd).(features{j});
        
        if any(isnan(pos)) % not yet tracked
            set(featureTexts{j}, 'color', nanColor);
            setColor(featurePoints{j}, nanColor);
        else
            set(featureTexts{j}, 'position', pos, 'color', labelledColor);
            setColor(featurePoints{j}, labelledColor);
            setPosition(featurePoints{j}, pos);
        end
    end
    
    set(selectedCirc, 'visible', 'off')
    selectedPoint = nan;
end




function movePoint(featureInd)
    pos = getPosition(featurePoints{featureInd});
    trainingData(structInd).(features{featureInd}) = pos;
    setColor(featurePoints{featureInd}, labelledColor);
    set(featureTexts{featureInd}, 'position', pos + textOffset, 'color', labelledColor);
    
    set(selectedCirc, 'XData', pos(1), 'YData', pos(2), 'visible', 'on')
    selectedPoint = featureInd;
end




% make sure points don't go out of frame
function constrainedPos = constrainPosition(newPos)
    constrainedPos(1) = min(newPos(1), vid.Width);
    constrainedPos(1) = max(constrainedPos(1), 1);
    constrainedPos(2) = min(newPos(2), vid.Height);
    constrainedPos(2) = max(constrainedPos(2), 1); 
end


% make sure tail end stays on edge of train
function constrainedPos = tailEndConstrainPosition(newPos)
    
    if newPos(2)<1
        constrainedPos(1) = min(newPos(1), vid.Width);
        constrainedPos(1) = max(constrainedPos(1), 1);
        constrainedPos(2) = 1;
    elseif newPos(2)>vid.Height
        constrainedPos(1) = min(newPos(1), vid.Width);
        constrainedPos(1) = max(constrainedPos(1), 1);
        constrainedPos(2) = vid.Height;
    else
        constrainedPos(1) = 1;
        constrainedPos(2) = newPos(2);
    end
end

end

