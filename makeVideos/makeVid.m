function makeVid(session)

% user settings
dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';
obsPosRange = [.31 .445];
maxTrialTime = 1; % trials exceeding maxTrialTime will be trimmed to this duration (s)
playbackSpeed = .1;



% initializations
webCamExists = exist([dataDir session '\webCam.csv'], 'file');
vidTop = VideoReader([dataDir session '\runTop.mp4']);
vidBot = VideoReader([dataDir session '\runBot.mp4']);
if webCamExists; vidWeb = VideoReader([dataDir session '\webCam.avi']); end

vidWriter = VideoWriter([dataDir session '\edited.mp4'], 'MPEG-4');
set(vidWriter, 'FrameRate', round(vidTop.FrameRate * playbackSpeed))
open(vidWriter)

load([dataDir session '\run.mat'], 'touch');
load([dataDir session '\runAnalyzed.mat'], 'obsPositions', 'obsTimes',...
                                           'wheelPositions', 'wheelTimes',...
                                           'obsOnTimes', 'obsOffTimes');
load([dataDir session '\frameTimeStamps.mat'], 'timeStamps')
load([dataDir session '\webCamTimeStamps.mat'], 'webCamTimeStamps')
maxFrames = vidTop.FrameRate * maxTrialTime;

obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes); % correct for drift in obstacle position readings



% edit video
w = waitbar(0, 'editing video...');

for i = 1:length(obsOnTimes)
    
    % find trial indices
    startInd = find(obsTimes>obsOnTimes(i)  & obsPositions>=obsPosRange(1), 1, 'first');
    endInd   = find(obsTimes<obsOffTimes(i) & obsPositions<=obsPosRange(2), 1, 'last');
    
    % get frame indices
    endTime = min(obsTimes(startInd)+maxTrialTime, obsTimes(endInd));
    frameInds = find(timeStamps>obsTimes(startInd) & timeStamps<endTime);
    
    if webCamExists
        % get webCame frame indices
        webFrameInds = find(webCamTimeStamps>obsTimes(startInd) & webCamTimeStamps<endTime);
        webFrames = read(vidWeb, [webFrameInds(1) webFrameInds(end)]);
        webFrames = squeeze(webFrames(:,:,1,:)); % collapse color dimension

        % increase framerate using interpolation    
        webFrames = double(webFrames);
        webFramesInterp = nan(size(webFrames,1), size(webFrames,2), length(frameInds));

        for j = 1:size(webFrames,1)
            for k = 1:size(webFrames,2)

                webFramesInterp(j,k,:) = interp1(webCamTimeStamps(webFrameInds),...
                                                 squeeze(webFrames(j,k,:)),...
                                                 timeStamps(frameInds),...
                                                 'linear', 'extrap');        
            end
        end
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
                frameWeb = webFramesInterp(:,:,j);
                frameWeb = imresize(frameWeb, (size(frame,2)/size(frameWeb,2)));
                frame = [frame; frameWeb];
            end

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


