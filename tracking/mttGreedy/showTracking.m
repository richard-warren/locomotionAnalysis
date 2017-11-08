function showTracking(vid, locations, labels, vidDelay, paws)
    

currentFrame = 1;
sampleFrame = rgb2gray(read(vid,currentFrame));
totalFrames = vid.NumberOfFrames;
cmap = winter(length(paws));

% prepare figure
close all;
fig = figure('position', [567 383 698 400], 'color', 'black', 'keypressfcn', @changeFrames);
colormap gray


rawIm = image(sampleFrame, 'CDataMapping', 'scaled');
rawAxis = gca;
set(rawAxis, 'visible', 'off')
hold on;
scatterPtsAll = scatter(rawAxis, 0, 0, 50, 'red', 'filled', 'linewidth', 2);
scatterPts =    scatter(rawAxis, zeros(1,length(paws)), zeros(1,length(paws)), linspace(50,500,4), cmap, 'linewidth', 3); hold on
scatter(rawAxis, [1 vid.Width 1 vid.Width], [1 1 vid.Height vid.Height], linspace(50,500,4), cmap, 'linewidth', 3); hold on

playing = true;
paused = false;


% main loop
while playing
    while paused; pause(.001); end
    updateFrame(1);
end

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
            close(fig)
        else                            % OTHERWISE: close window
            paused = ~paused;
        end
    end
end

% update frame preview
function updateFrame(frameStep)
    
    currentFrame = currentFrame + frameStep;
    
    % get frame and sub-frames
    frame = rgb2gray(read(vid,currentFrame));
    frame = getFeatures(frame);
    
    % update figure
    set(rawIm, 'CData', frame);
    
    xs = [locations(currentFrame).x;
    ys = [locations(currentFrame).y;

    inds = labels(currentFrame,paws);
    
    set(scatterPts, 'XData', xs(inds), 'YData', ys(inds), 'visible', 'on');
    set(scatterPtsAll, 'XData', locations(currentFrame).x, 'YData', locations(currentFrame).y);
    
    % pause to reflcet on the little things...
    pause(vidDelay);
    if currentFrame==totalFrames; currentFrame = 0; end
end


end