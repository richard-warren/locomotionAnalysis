% this make a wide vid with top view plus drawing of obs off screen, intended to explain how obs appears and is controlled by the speed of the wheel
% should use trials in which there is variability in wheel speed, and ideally some in which mouse moves both forwards and backwards

% settings
session = '180225_000';
trials = [3]; % 3 6
prePostTime = [-0.5 0]; % (s) time to add to beginning and end of a trial (before and after obs is engaged)
playbackSpeed = 0.1;
fps = 250;
obsYLims = [71 107];
obsThickness = 12;

% initializations
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
    'obsPixPositions', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps');
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
frameDims = [vid.Height round(max(obsPixPositions))];
vidWriter = VideoWriter([getenv('OBSDATADIR') 'editedVid\' sprintf('%strials%s', session, num2str(trials))], 'MPEG-4');
set(vidWriter, 'FrameRate', round(fps*playbackSpeed));
open(vidWriter);


% iterate through trials
for i = 1:length(trials)
    
    trialBins = frameTimeStamps>=(obsOnTimes(trials(i))+prePostTime(1)) & ...
                frameTimeStamps<=(obsOffTimes(trials(i))+prePostTime(2));
    trialInds = find(trialBins);
    
    % iterate through frames within trial
    for j = 1:length(trialInds)
        
        % get vid frame and incorporate into wide frame
        vidFrame = rgb2gray(read(vid, trialInds(j)));
        frame = uint8(zeros(frameDims));
        frame(:,1:vid.Width) = vidFrame;
        
        % add obs to frame
        obsPos = round(obsPixPositions(trialInds(j)));
        if obsPos>vid.Width
            frame = addObsToFrame(frame, obsPos, obsThickness, obsYLims, 255);
        end
        
        % write to video
        writeVideo(vidWriter, frame);
    end
end



close(vidWriter)
disp('all done')