function obsPixPositions = trackObstacles(session, obsOnTimes, obsOffTimes, frameTimeStamps)

% !!! need to document





% user settings
dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';
pixThresh = 80;
obsMinThickness = 10;                                           



% initializations
vid = VideoReader([dataDir session '\runBot.mp4']);
totalFrames = vid.NumberOfFrames;
if mod(obsMinThickness,2)==1; obsMinThickness = obsMinThickness-1; end % ensure obsMinThickness is even, which ensures medFiltSize is odd // this way the filtered version doesn't shift by one pixel relative to the unfiltered version
medFiltSize = obsMinThickness*2+1;



% iternate through all obstacle on epochs
obsPixPositions = nan(1, totalFrames);


for i = 1:length(obsOnTimes)

    frameInds = find(frameTimeStamps>=obsOnTimes(i) & frameTimeStamps<=obsOffTimes(i));
    disp(i/length(obsOnTimes));
    
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
        
    end
end



