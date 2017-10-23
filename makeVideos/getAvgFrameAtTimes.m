function avgFrame = getAvgFrameAtTimes(session, view, times)

% returns the average of frames occuring at 'times'

% settings
fps = 250;
dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';

% initializations
minFrameDiff = (1/fps) * 2;
vid = VideoReader([dataDir session '\run' view '.mp4']);
load([dataDir session '\runAnalyzed.mat'], 'frameTimeStamps');


allFrames = nan(vid.Height, vid.Width, length(times));

for i = 1:length(times)
    
    [timeDif, frameInd] = min(abs(frameTimeStamps - times(i)));
    
    if timeDif < minFrameDiff
        frame = read(vid, frameInd);
        allFrames(:,:,i) = squeeze(frame(:,:,1));
    else
        fprintf('no frame found within %.2f seconds of %.2f\n', minFrameDiff, times(i))
    end
    
end

avgFrame = uint8(nanmean(allFrames, 3));