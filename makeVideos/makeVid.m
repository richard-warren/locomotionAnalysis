function makeVid(session)

% user settings
dataDir = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\';
% obsPosRange = [.31 .45]; % each trial will capture frames in which the obstacle was between these positions on the track // meters
obsPosRange = [.31 .445];
maxTrialTime = 1; % trials exceeding maxTrialTime will be trimmed to this duration (s)
playbackSpeed = .1;



% initializations
vidTop = VideoReader([dataDir session '\runTop.mp4']);
vidBot = VideoReader([dataDir session '\runBot.mp4']);
vidWriter = VideoWriter([dataDir session '\edited.mp4'], 'MPEG-4');
set(vidWriter, 'FrameRate', round(vidTop.FrameRate * playbackSpeed))
open(vidWriter)

load([dataDir session '\run.mat'], 'ObsLight', 'touch');
load([dataDir session '\runAnalyzed.mat'], 'obsPositions', 'obsTimes', 'wheelPositions', 'wheelTimes', 'targetFs');
load([dataDir session '\frameTimeStamps.mat'], 'timeStamps')
maxFrames = vidTop.FrameRate * maxTrialTime;

obsOnTimes = ObsLight.times(logical(ObsLight.level));
obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes); % correct for drift in obstacle position readings



% edit video
w = waitbar(0, 'editing video...');

for i = 1:length(obsOnTimes)
    
    % find trial indices
    startInd = find(obsTimes>obsOnTimes(i) & obsPositions>obsPosRange(1), 1, 'first');
    endInd   = find(obsTimes>obsOnTimes(i) & obsPositions>obsPosRange(2), 1, 'first');
    
    % get frame indices
    frameInds = find(timeStamps>obsTimes(startInd) & timeStamps<obsTimes(endInd));
    if diff(frameInds) > maxFrames
        frameInds(2) = frameInds(1) + maxFrames;
    end
    
    if isempty(frameInds) % if a block has NaN timestamps (which will happen when unresolved), startInd and endInd will be the same, and frameInds will be empty
        fprintf('skipping trial %i due to unresolved timeStamps\n', i)
    else
        
        for f = frameInds'

            % put together top and bot frames
            frameTop = rgb2gray(read(vidTop, f));
            frameBot = rgb2gray(read(vidBot, f));
            frame = imadjust([frameTop; frameBot]);

            % add trial info text
            frame = insertText(frame, [0 0], ['trial: ' num2str(i)]);

            % write frame to video
            writeVideo(vidWriter, frame);
        end

        % add blank frame between trials
        writeVideo(vidWriter, zeros(size(frame)))
    end
    
    % update waitbar
    waitbar(i/length(obsOnTimes))
    
end



close(w)
close(vidWriter)


