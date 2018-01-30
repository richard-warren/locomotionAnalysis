function correctTracking(vid, locations, frameInds, vidDelay, anchorPts, lineLocations)


% settings
circSize = 200;
vidSizeScaling = 1.5;
colors = hsv(4);


% initializations
playing = true;
paused = false;

currentFrameInd = frameInds(1);
sampleFrame = rgb2gray(read(vid,currentFrameInd));



% prepare figure
close all;
fig = figure('units', 'pixels', 'position', [600 400 vid.Width*vidSizeScaling vid.Height*vidSizeScaling],...
    'menubar', 'none', 'color', 'black', 'keypressfcn', @changeFrames);
colormap gray
preview = image(sampleFrame, 'CDataMapping', 'scaled'); hold on;
rawAxis = gca;
set(rawAxis, 'visible', 'off', 'units', 'pixels',...
    'position', [0 0 vid.Width*vidSizeScaling vid.Height*vidSizeScaling]);
circSizes = circSize * ones(1,length(anchorPts));


% prepare lines showing x locations of bottom tracked paws
lines = cell(4,1);
if exist('lineLocations', 'var')
    for i = 1:4
        lines{i} = line([0 0], [vid.Height vid.Height-50], 'color', cmap(i,:));
    end
end


% prepare impoints, draggable markers used to show / adjust tracking
for i = 1:4
    impoints{i} = impoint(gca, [10 10]*i);
    setColor(impoints{i}, colors(i,:));
end
keyboard


% load data
fields = fieldnames(locations);
dim2 = fields{2};



% main loop
while playing
    while paused; pause(.001); end
    updateFrame(1);
end
close(fig)


% ---------
% FUNCTIONS
% ---------

% keypress controls
function changeFrames(~,~)
    
    key = double(get(fig, 'currentcharacter'));
    
    if ~isempty(key) && isnumeric(key)
        
        if key==28                      % LEFT: move frame backward
            pause(.001);
            paused = true;
            updateFrame(-1);
        
        elseif key==29                  % RIGHT: move frame forward
            pause(.001);
            paused = true;
            updateFrame(1);
        
        elseif key==102                  % 'f': select frame
            pause(.001);
            paused = true;
            input = inputdlg('enter frame number');
            currentFrameInd = find(frameInds == str2double(input{1}));
            updateFrame(1);
        
        elseif key==27                  % ESCAPE: close window
            % !!! save results here!!!
            playing = false;
            paused = false;
        
        else                            % OTHERWISE: close window
            paused = ~paused;
        end
    end
end



% update frame preview
function updateFrame(frameStep)
    
    currentFrameInd = currentFrameInd + frameStep;
    if currentFrameInd > length(frameInds); currentFrameInd = 1;
    elseif currentFrameInd < 1; currentFrameInd = length(frameInds); end
    
    % get frame and sub-frames
    frame = rgb2gray(read(vid, frameInds(currentFrameInd)));
    
    
    % add vertical lines
    if exist('lineLocations', 'var')
        inds = lineLocations.x(frameInds(currentFrameInd),:);
        for j = 1:4
            set(lines{j}, 'XData', [inds(j) inds(j)])
        end
    end
    
    
    % add frame number
    frame = insertText(frame, [size(frame,2) size(frame,1)], num2str(frameInds(currentFrameInd)),...
                               'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
    
    % update figure
    set(preview, 'CData', frame);
    
    
    % !!! update positions of pointers
    for j = 1:4
        setPosition(impoints{j}, locations.x(frameInds(currentFrameInd), j), locations.y(frameInds(currentFrameInd), j));
    end

    
    
    % pause to reflcet on the little things...
    pause(vidDelay);
end


end













