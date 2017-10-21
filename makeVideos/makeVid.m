function makeVid(session)

% edits a video of mouse jumping over obstacles s.t. obstacle trials are
% kept and everything else is edited out. obsPosRange is in m and defines
% the start and end position of the obstacle along the track that each
% trial should include.


% user settings
dataDir = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\';
% obsPosRange = [.1 .5]; %[.31 .445]; % (m)
obsPosRange = [.25 .445]; % (m)
maxTrialTime = 2; % trials exceeding maxTrialTime will be trimmed to this duration (s)
playBackSpeed = .1;



% initializations
webCamExists = exist([dataDir session '\webCam.csv'], 'file');
vidTop = VideoReader([dataDir session '\runTop.mp4']);
vidBot = VideoReader([dataDir session '\runBot.mp4']);
if webCamExists; vidWeb = VideoReader([dataDir session '\webCam.avi']); end

fps = round(vidTop.FrameRate * playBackSpeed);
maxFps = 150; % fps > 150 can be accomplished using 'Motion JPEG AVI' as second argument to VideoWriter, but quality of video is worse
if fps>maxFps
    playBackSpeed = maxFps/vidTop.FrameRate;
    fprintf('WARNING: requested playback speed not possible. speed set to %.2f', playBackSpeed);
    fps = maxFps;
end

vidWriter = VideoWriter(sprintf('%s%s\\%sspeed%.2f', dataDir, session, session, playBackSpeed), 'MPEG-4');
set(vidWriter, 'FrameRate', fps, 'Quality', 50)
open(vidWriter)

load([dataDir session '\run.mat'], 'touch');
load([dataDir session '\runAnalyzed.mat'], 'obsPositions', 'obsTimes',...
                                           'wheelPositions', 'wheelTimes',...
                                           'obsOnTimes', 'obsOffTimes',...
                                           'frameTimeStamps', 'webCamTimeStamps');

obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes); % correct for drift in obstacle position readings



% edit video
w = waitbar(0, 'editing video...');

for i = 1:length(obsOnTimes)
    
    % find trial indices
    startInd = find(obsTimes>obsOnTimes(i)  & obsPositions>=obsPosRange(1), 1, 'first');
    endInd   = find(obsTimes<obsOffTimes(i) & obsPositions<=obsPosRange(2), 1, 'last');
    
    % get frame indices
    endTime = min(obsTimes(startInd)+maxTrialTime, obsTimes(endInd));
    frameInds = find(frameTimeStamps>obsTimes(startInd) & frameTimeStamps<endTime);
    
    if webCamExists
        % get webCame frame indices
        webFrameInds = find(webCamTimeStamps>obsTimes(startInd) & webCamTimeStamps<endTime);
        webFrames = read(vidWeb, [webFrameInds(1) webFrameInds(end)]);
        webFrames = squeeze(webFrames(:,:,1,:)); % collapse color dimension
        
        % interpolate webFrames to number of inds in frameInds
        webFramesInterp = interp1(webCamTimeStamps(webFrameInds), 1:length(webFrameInds), frameTimeStamps(frameInds), 'nearest', 'extrap');
    end
    
    if isempty(frameInds) % if a block has NaN timestamps (which will happen when unresolved), startInd and endInd will be the same, and frameInds will be empty
        fprintf('skipping trial %i\n', i)
    else
        
        for j = 1:length(frameInds)
            
            % put together top and bot frames
            frameTop = rgb2gray(read(vidTop, frameInds(j)));
            frameBot = rgb2gray(read(vidBot, frameInds(j)));
            frame = imadjust([frameTop; frameBot]);
            
            % add webCam view
            if webCamExists
                frameWeb = webFrames(:,:,webFramesInterp(j));
                frameWeb = imresize(frameWeb, (size(frame,2)/size(frameWeb,2)));
                frame = [frame; frameWeb];
            end

            % add trial info text
            frame = insertText(frame, [0 0], num2str(i), 'BoxColor', 'yellow');

            % write frame to video
            writeVideo(vidWriter, frame);
        end

        % add blank frame between trials
        writeVideo(vidWriter, zeros(size(frame)));
    end
    
    % update waitbar
    waitbar(i/length(obsOnTimes))
    
end


close(w)
close(vidWriter)


