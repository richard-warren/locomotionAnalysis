function labelFrames(trainingSetName, features)

% store locations in struct // restrict point locations // ensure a point has been moved before saving it's value

% settings
vidScaling = 2;
textOffset = [5 0];
nanColor = [1 0 0];
labelledColor = [1 1 0];

% initializations
defaultPositions = repmat([10 10], length(features), 1) .* repmat(1:length(features),2,1)';
load([getenv('TRAININGEXAMPLESDIR') 'deepLabCut\' trainingSetName '\trainingData.mat'], 'trainingData');
fields = fieldnames(trainingData);
structInd = 1;



% create figure
currentSession = trainingData(structInd).session;
currentFrame = trainingData(structInd).frameNum;
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' currentSession '\runBot.mp4']);
fig = figure('name', sprintf('%s, frame %i', currentSession, currentFrame), ...
    'menubar', 'none', 'color', 'white', 'keypressfcn', @keypress, 'position', [200 200 vid.Width*vidScaling, vid.Height*vidScaling]);
colormap gray
im = image(rgb2gray(read(vid, currentFrame)), 'CDataMapping', 'scaled');
set(gca, 'position', [0 0 1 1], 'visible', 'off')



% add fields that have not yet been created, and creatable draggable objects
featurePoints = cell(1, length(features));
featureTexts = cell(1, length(features));

for i = 1:length(features)
    
    % initialize non-existent features
    if ~ismember(features{i}, fields)
        nanEntries = num2cell(nan(length(trainingData),1));
        [trainingData.(features{i})] = nanEntries{:};
        fprintf('creating field %s\n', features{i});
    end
    
    % create draggables for features
    featureTexts{i} = text(defaultPositions(i,1)+textOffset(1), defaultPositions(i,2)+textOffset(2), features{i});
    featurePoints{i} = impoint(gca, defaultPositions(i,:));
    addNewPositionCallback(featurePoints{i}, @(x) movePoint(i));
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
                
            case 115 % 's'
                save([getenv('TRAININGEXAMPLESDIR') 'deepLabCut\' trainingSetName '\trainingData.mat'], 'trainingData');
                m = msgbox('saving progress...'); pause(.5); close(m)
        end
    end
end


function updateFrame(frameStep)
    
    % move to next frame
    structInd = structInd + frameStep;
    currentFrame = trainingData(structInd).frameNum;
    
    % load new video if switching to new session
    if ~strcmp(currentSession, trainingData(structInd).session)
        currentSession = trainingData(structInd).session;
        vid = VideoReader([getenv('OBSDATADIR') 'sessions\' currentSession '\runBot.mp4']);
    end
    
    % update frame
    frame = rgb2gray(read(vid, trainingData(structInd).frameNum));
    set(im, 'CData', frame);
    set(fig, 'name', sprintf('%s, frame %i', currentSession, currentFrame))
    
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
end


function movePoint(featureInd)
    pos = getPosition(featurePoints{featureInd});
    trainingData(structInd).(features{featureInd}) = pos;
    setColor(featurePoints{featureInd}, labelledColor);
    set(featureTexts{featureInd}, 'position', pos + textOffset, 'color', labelledColor);
end

end

