function obsPixPositions = trackObstacles(sessionPath, obsOnTimes, obsOffTimes, frameTimeStamps, showTracking)

% !!! need to document


% user settings
pixThresh = 80;
obsMinThickness = 10;                                           



% initializations
vid = VideoReader([sessionPath '\runBot.mp4']);
totalFrames = vid.NumberOfFrames;
if mod(obsMinThickness,2)==1; obsMinThickness = obsMinThickness-1; end % ensure obsMinThickness is even, which ensures medFiltSize is odd // this way the filtered version doesn't shift by one pixel relative to the unfiltered version
medFiltSize = obsMinThickness*2+1;
if ~exist('showTracking', 'var')
    showTracking = false;
end

% prepare figure if showTracking enabled
if showTracking
    
    fig = figure('position', [1923 35 vid.Width vid.Height*2], 'menubar', 'none', 'color', 'black');
    frame = rgb2gray(read(vid,1));
    
    axRaw = subplot(2,1,1, 'units', 'pixels');
    imRaw = imshow(frame);
    
    axSum = subplot(2,1,2, 'units', 'pixels');
    imSum = imshow(frame);
    
    set(axRaw, 'position', [1 size(frame,1)+1 size(frame,2) size(frame,1)]);
    set(axSum, 'position', [1 1 size(frame,2) size(frame,1)]);
end


% iternate through all obstacle on epochs
obsPixPositions = nan(1, totalFrames);


for i = 1:length(obsOnTimes)

    % get frame indices for current obstacle epoch
    frameInds = find(frameTimeStamps>=obsOnTimes(i) & frameTimeStamps<=obsOffTimes(i));
    
    % iterate through all frames within epoch
    for j = 1:length(frameInds)
        
        % get frame
        frame = rgb2gray(read(vid, frameInds(j)));
        
        % sum across columns and normalize
        colSums = sum(frame,1) / size(frame,1);
        
        % threshold and median filter to remove thin threshold crossings with few adjacent columns
        colThreshed = colSums > pixThresh;
        colThreshed = medfilt1(double(colThreshed), medFiltSize);
        
        obsPixPositions(frameInds(j)) = mean(find(colThreshed));
        
        % update figure if showTracking enabled
        if showTracking
            threshFrame = frame .* uint8(repmat(colThreshed, size(frame,1), 1));
            colFrame = repmat(colSums, size(frame,1), 1);
            
            set(imRaw, 'CData', frame .* uint8(~threshFrame));
            set(imSum, 'CData', colFrame);
            
            pause(.001)
            
        end
        
    end
end

if showTracking; close(fig); end



