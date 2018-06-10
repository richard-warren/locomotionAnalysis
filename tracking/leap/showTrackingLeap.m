function showTrackingLeap(session, vidDelay)

% settings
circSize = 150;
vidSizeScaling = 1.25;
colorMap = 'hsv';
trainingDims = [404 396]; % height width
showFeature = [1 1 1 1 0 0 0 1 1 1 1 1 0 0 1 0 0 0];



% get videos
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);

% get locations data and convert to 3d matrix
locations = h5read([getenv('OBSDATADIR') 'sessions\' session '\trackingDataRaw.preds.h5'], '/positions_pred');
confidences = h5read([getenv('OBSDATADIR') 'sessions\' session '\trackingDataRaw.preds.h5'], '/conf_pred');

% initializations
frameInds = getFramesToShow(session);
scaling = [vidTop.Width/trainingDims(2), ...
           (vidBot.Height+vidTop.Height)/trainingDims(1)];
cmap = eval(sprintf('%s(%i);', colorMap, size(locations,1)));

% upscale tracking to original resolution
locations = double(locations) .* repmat(scaling, size(locations,1), 1, size(locations,3));

% set up figure
hgt = (vidBot.Height+vidTop.Height);
fig = figure('units', 'pixels', 'position', [600 400 vidBot.Width*vidSizeScaling hgt*vidSizeScaling],...
    'menubar', 'none', 'color', 'black', 'keypressfcn', @changeFrames);
colormap gray
imPreview = image(zeros(hgt, vidBot.Width), 'CDataMapping', 'scaled'); hold on;
imAxis = gca;
set(imAxis, 'visible', 'off', 'units', 'pixels',...
    'position', [0 0 vidBot.Width*vidSizeScaling hgt*vidSizeScaling]);



% set up scatter points for tracked features
scatterLocations = scatter(imAxis, zeros(1,size(locations,1)), zeros(1,size(locations,1)),...
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
    if frameInd < 1; frameInd = length(frameInds);
    elseif frameInd > length(frameInds); frameInd = 1; end
    
    
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
    
    % upate scatter positions
%     validBins = confidences(:,frameInds(frameInd))>10 & showFeature';
    validBins = showFeature';
    set(scatterLocations, 'XData', locations(:,1,frameInds(frameInd)) .* double(validBins), ...
                          'YData', locations(:,2,frameInds(frameInd)));

    % pause to reflcet on the little things...
    pause(vidDelay);
end



end