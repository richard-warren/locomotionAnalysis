function showLocations(vid, locations, labels, showPotentialLocations, vidDelay, anchorPts, startFrame, lineLocations)
    
% settings
circSize = 200;
vidSizeScaling = 1.5;


% initializations
if ~exist('startFrame', 'var'); startFrame = 1; end
currentFrame = startFrame;
sampleFrame = rgb2gray(read(vid,currentFrame));
cmap = winter(length(anchorPts));

% prepare figure
close all;
fig = figure('units', 'pixels', 'position', [600 400 vid.Width*vidSizeScaling vid.Height*vidSizeScaling],...
    'menubar', 'none', 'color', 'black', 'keypressfcn', @changeFrames);
colormap gray


rawIm = image(sampleFrame, 'CDataMapping', 'scaled'); hold on;
rawAxis = gca;
set(rawAxis, 'visible', 'off', 'units', 'pixels', 'position', [0 0 vid.Width*vidSizeScaling vid.Height*vidSizeScaling]);
circSizes = circSize * ones(1,length(anchorPts)); % linspace(50,500,4)

if showPotentialLocations
    scatterPotentialLocations = scatter(rawAxis, 0, 0, 50, 'white', 'filled', 'linewidth', 2);
end
scatterLocations =    scatter(rawAxis, zeros(1,length(anchorPts)), zeros(1,length(anchorPts)), circSizes, cmap, 'filled', 'linewidth', 3); hold on
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
        
        elseif key==27                  % ESCAPE: close window
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
    if currentFrame > vid.NumberOfFrames; currentFrame = 1;
    elseif currentFrame < 1; currentFrame = vid.NumberOfFrames; end
    
    % get frame and sub-frames
    frame = rgb2gray(read(vid,currentFrame));
    frame = imadjust(uint8(frame), [.05 1], [0 1]);
    frame = getFeatures(frame);
    
    
    % add vertical lines
    if exist('lineLocations', 'var')
        inds = lineLocations.x(currentFrame,:);
        inds = round(inds(~isnan(inds)));
        frame(:, inds) = 255;
    end
    
    
    % add frame number
    frame = insertText(frame, [size(frame,2) size(frame,1)], num2str(currentFrame),...
                               'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
    
    % update figure
    set(rawIm, 'CData', frame);
    
    xs = [locations(currentFrame).x; 0]; % add zero for occluded state
    ys = [locations(currentFrame).y; 0]; % add zero for occluded state

    inds = labels(currentFrame,:);
    inds(isnan(inds)) = length(xs);
    
    set(scatterLocations, 'XData', xs(inds), 'YData', ys(inds), 'visible', 'on');
    if showPotentialLocations
        set(scatterPotentialLocations, 'XData', locations(currentFrame).x, 'YData', locations(currentFrame).y);
    end
    
    % pause to reflcet on the little things...
    pause(vidDelay);
end


end