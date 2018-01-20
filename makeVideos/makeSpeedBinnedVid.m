function makeSpeedBinnedVid(trialNum, speedLimits, allTrialsData, pauseAtWisk)

% note that allTrialsData is genreated using getAllSessionTrialSpeeds function
% edits together a montage of trials across all sessions in experiments in which mouse was running within speedLimits around position of obstacle centered at wheel


% temp
% experiments = {'obsBr', 'obsWiskMrk', 'obsCompareOn', 'obsCompareOnOff'};
% trialNum = 10;
% speedLimits = [.2 .25];

% settings
obsPosRange = [.25 .445]; % edit video between these obsPoses
playBackSpeed = .1;
% obsPos = 0.3820; % obsPos at which obs is at center of wheel
% posRange = .06; % computer vel between obsPos plus or minus posRange
% allTrialsData = getAllSessionTrialSpeeds(experiments, obsPos, posRange);
originalFps = 250;
maxTrialTime = 1.5; % trials exceeding maxTrialTime will be trimmed to this duration (s)
obsBotThickness = 15;
wiskPos = .33; % this is rough estimate
pauseTime = 1.5;

% initializations
editedDir = [getenv('OBSDATADIR') 'editedVid\'];
vidWriter = VideoWriter(sprintf('%sspeedMontage%.2f-%.2f', editedDir, speedLimits(1), speedLimits(2)), 'MPEG-4');
set(vidWriter, 'FrameRate', round(originalFps * playBackSpeed), 'Quality', 50)
open(vidWriter)
pauseFrames = pauseTime * playBackSpeed*originalFps;


% only keep trials within posRange
validBins = [allTrialsData.vel]>=speedLimits(1) & [allTrialsData.vel]<=speedLimits(2) &...
            [allTrialsData.videoRecorded]; % only include trials where a video was actually recorded...
allValidTrials = allTrialsData(validBins);

% get random inds
inds = randperm(length(allValidTrials));
inds = inds(1 : min(trialNum, length(allValidTrials)));
inds = sort(inds);



% edit video
w = waitbar(0, 'editing video...');
lastSession = '';

for i = inds
    
    currentSession = allValidTrials(i).session;
    
    % if we are switching to a new session/vid, reload data
    if ~strcmp(currentSession, lastSession)
        
        vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' currentSession '\runTop.mp4']);
        vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' currentSession '\runBot.mp4']);
        
        load([getenv('OBSDATADIR') 'sessions\' currentSession '\runAnalyzed.mat'],...
            'obsPositions', 'obsTimes',...
            'wheelPositions', 'wheelTimes',...
            'obsOnTimes', 'obsOffTimes',...
            'obsPixPositions', 'frameTimeStamps');
        obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes); % correct for drift in obstacle position readings
        
        % some videos are cropped differently, so this allows me to reshape down the line to a uniform size
        if ~exist('frameSize', 'var'); frameSize = [vidTop.Height+vidBot.Height vidTop.Width]; end
    end
    
    lastSession = currentSession;
    trial = allValidTrials(i).trialNum;
    
    % find trial indices
    startInd = find(obsTimes>obsOnTimes(trial)  & obsPositions>=obsPosRange(1), 1, 'first');
    endInd   = find(obsTimes<obsOffTimes(trial) & obsPositions<=obsPosRange(2), 1, 'last');
    
    % get frame indices
    endTime = min(obsTimes(startInd)+maxTrialTime, obsTimes(endInd));
    frameInds = find(frameTimeStamps>obsTimes(startInd) & frameTimeStamps<endTime);
    
    % get wisk contact ind
    positionsAtFrameInds = interp1(obsTimes, obsPositions, frameTimeStamps(frameInds));
    [minDiff, minInd] = min(abs(positionsAtFrameInds - wiskPos));
    if minDiff<.005
        wiskTouchFrameInd = minInd;
    else
        wiskTouchFrameInd = [];
        disp('couldnt find wisk contact frame lololol')
    end
    
    
    if isempty(frameInds) % if a block has NaN timestamps (which will happen when unresolved), startInd and endInd will be the same, and frameInds will be empty
        fprintf('skipping trial %i\n', i)
    else
        
        for j = 1:length(frameInds)
            
            % top
            frameTop = rgb2gray(read(vidTop, frameInds(j)));
            
            % bot
            frameBot = rgb2gray(read(vidBot, frameInds(j)));
            frameBot = addObsToFrame(frameBot, obsPixPositions(frameInds(j)), obsBotThickness, [1, size(frameBot,1)], 255);
            
            % combine frames
            frame = cat(1, frameTop, frameBot);
            if any(size(frame) ~= frameSize)
            	frame = imresize(frame, frameSize);
                disp('resized frame! lol')
            end
            
            % add trial info text
            frame = insertText(frame, [size(frame,2) size(frame,1)], sprintf('%s: %s', allValidTrials(i).mouse, allValidTrials(i).session),...
                               'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
            
            % write frame to video
            writeVideo(vidWriter, frame);
            
            % add pause
            if pauseAtWisk && j==wiskTouchFrameInd
                for k = 1:pauseFrames; writeVideo(vidWriter, frame); end
            end
        end

        % add blank frame between trials
        writeVideo(vidWriter, zeros(size(frame)));
    end
    
    waitbar(find(i==inds) / length(inds))
end

close(vidWriter)
close(w)









