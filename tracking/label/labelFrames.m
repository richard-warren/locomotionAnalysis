function labelFrames(trainingSet, skeleton, varargin)

% settings
s.vidScaling = 1.5;
s.textOffset = [5 0];
s.nanColor = [1 0 0];
s.lineColor = [0 0 1];
s.labeledColor = [1 1 0];
s.invertFrame = false;
s.textRotation = 45;
s.hideWhileMoving = true;  % whether to hide other tracked features while moving a single feature

% todo: hide other features while moving one feature // move unlabeled
% points to corner of screen?



% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
skeletonTbl = readtable(skeleton);
load(trainingSet, 'trainingData');
features = skeletonTbl.name;

hgt = max(cellfun(@(x) size(x,1), {trainingData.frame}));
wid = max(cellfun(@(x) size(x,2), {trainingData.frame}));
defaultPositions = [randi(wid, length(features), 1), randi(hgt, length(features), 1)];
fields = fieldnames(trainingData);  % columns already in trainingData
ind = 1;
selectedPoint = nan;



% create figure
fig = figure('name', sprintf('%s, frame %i', trainingData(ind).session, trainingData(ind).frameNum), ...
    'menubar', 'none', 'color', 'white', 'keypressfcn', @keypress, 'position', [200 0 wid*s.vidScaling, hgt*s.vidScaling]);
colormap gray
im = image(trainingData(ind).frame, 'CDataMapping', 'scaled'); hold on
set(gca, 'position', [0 0 1 1], 'visible', 'off')
selectedCirc = scatter(0, 0, 200, [.5 .5 1], 'linewidth', 3, 'visible', 'off'); hold on



% create arrays for text and impoints objects
points = cell(1, length(features));
texts = cell(1, length(features));

for i = 1:length(features)
    
    % initialize non-existent features
    if ~ismember(features{i}, fields)
        nanEntries = mat2cell(nan(length(trainingData),2), ones(1,length(trainingData)), 2);
        [trainingData.(features{i})] = nanEntries{:};
        fprintf('creating field: %s\n', features{i});
    end
    
    % create draggables for features
    texts{i} = text(defaultPositions(i,1)+s.textOffset(1), defaultPositions(i,2)+s.textOffset(2), features{i}, ...
        'interpreter', 'none', 'rotation', s.textRotation);
    points{i} = impoint(gca, defaultPositions(i,:));
    addNewPositionCallback(points{i}, @(x) movePoint(i));
    setPositionConstraintFcn(points{i}, @constrainPosition)
    
    % add special constraint for end of tail
    if strcmp(features{i}, 'tailEnd')
        setPositionConstraintFcn(points{i}, @tailEndConstrainPosition)
    end
end

% create lines connecting feature pairs
hasParentInds = find(~cellfun(@isempty, skeletonTbl.parent));  % bins of features with a parent feature, as determined by spreadsheet
numEdges = length(hasParentInds);
lineColors = hsv(numEdges);
edges = nan(numEdges, 2);  % number of edges X 2 matrix storing connected features
featurePairLines = cell(1,numEdges);

for i = 1:numEdges
    edges(i,1) = hasParentInds(i);
    edges(i,2) = find(strcmp(skeletonTbl.name, skeletonTbl.parent(hasParentInds(i))));  % index of parent feature
    featurePairLines{i} = line([0 0], [0 0], 'color', lineColors(i,:)); hold on
    uistack(featurePairLines{i}, 'bottom');
    uistack(featurePairLines{i}, 'up');
end

% move text to bottom of feature stack (so it doesn't prevent clicking on impoints)
for i = 1:length(features)
    uistack(texts{i}, 'bottom');
    uistack(texts{i}, 'up');
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
                
            % d: delete selection
            case 100
                if ~isnan(selectedPoint)
                    trainingData(ind).(features{selectedPoint}) = [nan nan];
                    updateFrame(0);
                end
                
            % h: toggle hiding of feature labels
            case 104
                for j = 1:length(texts); set(texts{j}, 'Color', 'none'); end
                
            % f: move to next frame that is not yet included
            case 102
                newStructInd = find(~[trainingData.includeFrame] & 1:length(trainingData)>ind, 1, 'first'); % find first frame with nonlabelled part that is later than current frame
                if isempty(newStructInd); newStructInd = find(~[trainingData.includeFrame], 1, 'first'); end % if couldnt find any, search for nonlabelled part starting from beginning
                if ~isempty(newStructInd); ind = newStructInd; end % otherwise, keep current frame
                updateFrame(0);
            
            % n: move to ind
            case 110
                input = inputdlg('enter new n...');
                ind = str2num(input{1});
                updateFrame(0);
                
            % ENTER: include frame for analysis
            case 13
                trainingData(ind).includeFrame = ~trainingData(ind).includeFrame;
                updateFrame(0);
                
            % s: save progress
            case 115
                save(trainingSet, 'trainingData');
                m = msgbox('saving progress...'); pause(.5); close(m)
        end
    end
end




function updateFrame(frameStep)
    
    % move to next frame
    ind = ind + frameStep;
    if ind==0; ind=length(trainingData); end
    if ind>length(trainingData); ind=1; end
    
    % update frame
    frame = trainingData(ind).frame;
    if s.invertFrame; frame = 255-frame; end
    set(im, 'CData', frame);
    set(gca, 'position', [0 0 1 1]);
    if trainingData(ind).includeFrame; includedStr = '(included)'; else; includedStr = ''; end
    set(fig, 'name', ...
        sprintf('%s, frame %i %s (n = %i/%i), %i unlabelled', ...
        trainingData(ind).session, trainingData(ind).frameNum, includedStr, ind, length(trainingData), sum([trainingData.includeFrame]==0)))
    
    % update colors and positions of points
    for j = 1:length(features)
        pos = trainingData(ind).(features{j});
        
        % update colors and positions
        if any(isnan(pos)) % not yet tracked
            set(texts{j}, 'color', s.nanColor);
            setColor(points{j}, s.nanColor);
        else
            set(texts{j}, 'position', pos, 'color', s.labeledColor);
            setColor(points{j}, s.labeledColor);
            setPosition(points{j}, pos);
        end
    end
    
    % update feature pair lines
    updateFeatureLines();
    
    set(selectedCirc, 'visible', 'off')
    selectedPoint = nan;
end




function movePoint(featureInd)
    
    pos = getPosition(points{featureInd});
    trainingData(ind).(features{featureInd}) = pos;
    setColor(points{featureInd}, s.labeledColor);
    set(texts{featureInd}, 'position', pos + s.textOffset, 'color', s.labeledColor);
    
    set(selectedCirc, 'XData', pos(1), 'YData', pos(2), 'visible', 'on')
    selectedPoint = featureInd;
    updateFeatureLines();
end




% make sure points don't go out of frame
function constrainedPos = constrainPosition(newPos)
    constrainedPos(1) = min(newPos(1), wid);
    constrainedPos(1) = max(constrainedPos(1), 1);
    constrainedPos(2) = min(newPos(2), hgt);
    constrainedPos(2) = max(constrainedPos(2), 1); 
end

% update feature pair lines
function updateFeatureLines()
    for j = 1:size(edges,1)
        inds = edges(j,:);
        pos1 = getPosition(points{inds(1)});
        pos2 = getPosition(points{inds(2)});
        set(featurePairLines{j}, 'xdata', [pos1(1) pos2(1)], 'ydata', [pos1(2) pos2(2)]);
    end
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
        constrainedPos(2) = hgt;
    else
        constrainedPos(1) = 1;
        constrainedPos(2) = newPos(2);
    end
end

end

