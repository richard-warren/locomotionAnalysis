function showLocations(vid, frameInds, potentialLocations, locations,...
    showPotentialLocations, vidDelay, anchorPts, cmap, lineLocations)
    
% settings
circSize = 200;
vidSizeScaling = 1.5;
% lineMaskWid = 15;


% initializations
currentFrame = 1;
sampleFrame = rgb2gray(read(vid,currentFrame));
fields = fieldnames(locations);
dim2 = fields{2};


% prepare figure
close all;
fig = figure('units', 'pixels', 'position', [600 400 vid.Width*vidSizeScaling vid.Height*vidSizeScaling],...
    'menubar', 'none', 'color', 'black', 'keypressfcn', @changeFrames);
colormap gray


rawIm = image(sampleFrame, 'CDataMapping', 'scaled'); hold on;
rawAxis = gca;
set(rawAxis, 'visible', 'off', 'units', 'pixels',...
    'position', [0 0 vid.Width*vidSizeScaling vid.Height*vidSizeScaling]);
circSizes = circSize * ones(1,length(anchorPts));
% circSizes = linspace(100,500,4);

lines = cell(4,1);
if exist('lineLocations', 'var')
    for i = 1:4
        lines{i} = line([0 0], [vid.Height vid.Height-50], 'color', cmap(i,:));
    end
end
    

if showPotentialLocations
    scatterPotentialLocations = scatter(rawAxis, 0, 0, 50, 'white', 'filled', 'linewidth', 2);
end
scatterLocations = scatter(rawAxis, zeros(1,length(anchorPts)), zeros(1,length(anchorPts)),...
    circSizes, cmap, 'linewidth', 5); hold on
scatter(rawAxis, [anchorPts{1}(1) anchorPts{2}(1) anchorPts{3}(1) anchorPts{4}(1)] .* (vid.Width-1) + 1,...
                 [anchorPts{1}(2) anchorPts{2}(2) anchorPts{3}(2) anchorPts{4}(2)] .* (vid.Height-1) + 1,...
                 circSizes, cmap, 'filled', 'linewidth', 3);     % show anchor points

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
            currentFrame = find(frameInds == str2double(input{1}));
            updateFrame(1);
        
        elseif key==27                  % ESCAPE: close window
%             disp(key)
            playing = false;
            paused = false;
        
        else                            % OTHERWISE: close window
            paused = ~paused;
        end
    end
end



% update frame preview
function updateFrame(frameStep)
    
    currentFrame = currentFrame + frameStep;
    if currentFrame > length(frameInds); currentFrame = 1;
    elseif currentFrame < 1; currentFrame = length(frameInds); end
    
    % get frame and sub-frames
    frame = rgb2gray(read(vid,frameInds(currentFrame)));
    frame = imadjust(uint8(frame), [.05 1], [0 1]);
%     frame = getFeatures(frame);
    
    
    % add vertical lines
    if exist('lineLocations', 'var')
        inds = lineLocations.x(frameInds(currentFrame),:);
        for j = 1:4
            set(lines{j}, 'XData', [inds(j) inds(j)])
        end
%         inds = round(inds(~isnan(inds)));
%         frame(:, inds) = 255;
    end
    
    
    % add frame number
    frame = insertText(frame, [size(frame,2) size(frame,1)], num2str(frameInds(currentFrame)),...
                               'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
    
    % update figure
    set(rawIm, 'CData', frame);

    set(scatterLocations, 'XData', locations.x(frameInds(currentFrame),:), 'YData', locations.(dim2)((frameInds(currentFrame)),:), 'visible', 'on');

    if showPotentialLocations
        set(scatterPotentialLocations, 'XData', potentialLocations(frameInds(currentFrame)).x, 'YData', potentialLocations(frameInds(currentFrame)).(dim2));
    end
    
    % pause to reflcet on the little things...
    pause(vidDelay);
end


end