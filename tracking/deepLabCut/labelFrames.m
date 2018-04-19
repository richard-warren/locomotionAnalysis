function labelFrames(trainingSetName, features)

% settings
vidScaling = 2;

% initializations
load([getenv('TRAININGEXAMPLESDIR') 'deepLabCut\' trainingSetName '\trainingData.mat'], 'trainingData');
fields = fieldnames(trainingData);
structInd = 1;



% create figure
close all
currentSession = trainingData(structInd).session;
currentFrame = trainingData(structInd).frameNum;
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' currentSession '\runBot.mp4']);
fig = figure('name', sprintf('%s, frame %i', currentSession, currentFrame), ...
    'menubar', 'none', 'color', 'white', 'keypressfcn', @keypress, 'position', [200 200 vid.Width*vidScaling, vid.Height*vidScaling]);
colormap gray
im = image(rgb2gray(read(vid, currentFrame)), 'CDataMapping', 'scaled');
set(gca, 'position', [0 0 1 1])



% add fields that have not yet been created, and creatable draggable objects
featurePoints = cell(1, length(features));

for i = 1:length(features)
    
    % initialize non-existent features
    if ~ismember(features{i}, fields)
        nanEntries = num2cell(nan(length(trainingData),1));
        [trainingData.(features{i})] = nanEntries{:};
        fprintf('creating field %s\n', features{i});
    end
    
    % create draggables for features
    featurePoints{i} = impoint(gca, [10 10]*i);
end


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
            case 27% save current progress
                close(fig)
                
            case 115 % 's'
                save([getenv('TRAININGEXAMPLESDIR') 'deepLabCut\' trainingSetName '\trainingData.mat'], 'trainingData');
                disp('data saved')
        end
    end
end


function updateFrame(frameStep)
    
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
end


end

