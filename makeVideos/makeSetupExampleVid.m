function makeSetupExampleVid
% this make a wide vid with top view plus drawing of obs off screen, intended to explain how obs appears and is controlled by the speed of the wheel
% should use trials in which there is variability in wheel speed, and ideally some in which mouse moves both forwards and backwards

% settings
session = '180225_000';
trials = [34 50 29 60 110]; % back and forth trials: 34 3 6 57 87
prePostTime = [-.25 0]; % (s) time to add to beginning and end of a trial (before and after obs is engaged)
playbackSpeed = 0.15;
fps = 250;
obsYLims = [71 107];
obsThickness = 12;
constrastLims = [.2 1]; % pixels at these proportional values are mapped to 0 and 255
obsFadePixels = 50; % when obs enter right side of vid frame, fade out drawing of obs over the course of this many pixels

% initializations
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
    'obsPixPositions', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps');
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
frameDims = [vid.Height round(max(obsPixPositions))];
vidWriter = VideoWriter([getenv('OBSDATADIR') 'editedVid\' sprintf('setupExample%strials%s', session, num2str(trials))], 'MPEG-4');
set(vidWriter, 'FrameRate', round(fps*playbackSpeed));
open(vidWriter);


% iterate through trials
for i = 1%:length(trials)
    
    trialBins = frameTimeStamps>=(obsOnTimes(trials(i))+prePostTime(1)) & ...
                frameTimeStamps<=(obsOffTimes(trials(i))+prePostTime(2));
    trialInds = find(trialBins);
    
    % iterate through frames within trial
    for j = 1:length(trialInds)
        
        % get vid frame and incorporate into wide frame
        vidFrame = rgb2gray(read(vid, trialInds(j)));
        vidFrame = imadjust(vidFrame, constrastLims, [0 1]);
        frame = uint8(zeros(frameDims));
        frame(:,1:vid.Width) = vidFrame;
        
        % add obs to frame
        obsPos = round(obsPixPositions(trialInds(j)));
        if obsPos>=(vid.Width-obsFadePixels)
            brightness = ((obsPos)-(vid.Width-.5*obsFadePixels)) / obsFadePixels; % fade occurs on the left and right sides of frame edge
            brightness = min(brightness,1); brightness = max(brightness,0);
            frameObs = addObsToFrame(frame, obsPos, obsThickness, obsYLims, 255);
            frame = (1-brightness)*frame + (brightness)*frameObs;
        end
        
        % write to video
        writeVideo(vidWriter, frame);
%         if j==49; break; end
    end
end


close(vidWriter)
disp('all done')