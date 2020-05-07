function trainingData = createTrainingDataStruct(sessions, view, vidName, frameNum, varargin)

% settings
s.method = 'metadata';  % method for frame sampling // 'uniform', 'random' or 'metadata' // 'metadata' uses information from the spike recording to make sure ceretain events of interest are represeted

% if s.method is 'metadata', then the following fractions of frames are included from the following epochs
s.ledOnFraction = .1;  %  fraction of trials where reward LED is on
s.obsVisibleOn = .5;  %  fraction of trials where obstable is visible and engaged (as opposed to visible but returning to home position)
s.obsVisibleOff = .1;  %  fraction of trials where obstacle is visible while returning to home position
s.nearReward = .1;  %  fraction of trials to include near the reward times

s.rewardProximity = [0 2];  % (s) 'nearReward' frames are selected randomly from this temporal interval surrounding reward delivery
s.obsVisibleRun = [.29 .39];  % (m) esimate of obstacle positions along tracking in which it is visible from the run camera (run 'getAvgFrameAtObsLocation' to get values for this)
s.obsVisibleWisk = [.28 .31];  % (m) esimate of obstacle positions along tracking in which it is visible from the whisker camera (run 'getAvgFrameAtObsLocation' to get values for this)


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
sesFrames = ones(1, length(sessions)) * round(frameNum / length(sessions));  % frames to get for each session
sesFrames(end) = sesFrames(end) + (frameNum-sum(sesFrames));
trainingData(frameNum,1) = struct();
ind = 1;  % index for trainingData


% get frames for each session
for i = 1:length(sessions)
    
    % load sesssion data
    fprintf('%s: loading session data...\n', sessions{i})
    
    % load session data
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
        'rewardTimes', 'frameTimeStamps', 'frameTimeStampsWisk', 'obsOnTimes', 'obsOffTimes', 'obsPositions', 'obsTimes');
    vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, vidName));
    totalFrames = floor(vid.Duration*vid.FrameRate);

    if strcmp(view, 'wisk')
        frameTimes = frameTimeStampsWisk;
        lims = s.obsVisibleWisk;
    
    elseif strcmp(view, 'run')
        frameTimes = frameTimeStamps;
        lims = s.obsVisibleRun;
    else
        disp('WARNING: ''view'' must be either ''wisk'' or ''run''')
    end
    
    
    % get frames
    if strcmp(s.method, 'uniform')
        fprintf('%s: finding uniformly spaced frames...\n', sessions{i})
        frames = round(linspace(1, totalFrames, sesFrames(i)+2));
        frames = frames(2:end-1)';
    
    elseif strcmp(s.method, 'random')
        fprintf('%s: finding random frames...\n', sessions{i})
        frames = randsample(totalFrames, sesFrames(i))';
    
    elseif strcmp(s.method, 'metadata')
        fprintf('%s: finding frames based on session metadata...\n', sessions{i})
        
        % determine number of frames for each frame type
        frames = [];
        ledFrameNum = round(sesFrames(i) * s.ledOnFraction);
        obsOnFrameNum = round(sesFrames(i) * s.obsVisibleOn);
        obsOffFrameNum = round(sesFrames(i) * s.obsVisibleOff);
        rewardFrameNum = round(sesFrames(i) * s.nearReward);
        randFrameNum = sesFrames(i) - ledFrameNum - obsOnFrameNum - obsOffFrameNum - rewardFrameNum;
        
        % get led on frames
        times = rewardTimes(randsample(length(rewardTimes), ledFrameNum));
        frames(end+1:end+ledFrameNum) = knnsearch(frameTimes, times)+1;  % add +1 to make sure we get frames during and not right before LED
        
        % create logical vector coding whether obstacle is on for all frames
        isObsOn = zeros(size(frameTimes));
        isObsOn(knnsearch(frameTimes, obsOnTimes)) = 1;
        isObsOn(knnsearch(frameTimes, obsOffTimes)) = -1;
        isObsOn = cumsum(isObsOn);
        obsPosInterp = interp1(obsTimes, obsPositions, frameTimes);  % interpolate onto same time grid as frames
        
        % get obs visible while engaged frames
        obsOnFrames = find(obsPosInterp>lims(1) & obsPosInterp<lims(2) & isObsOn);
        frames(end+1:end+obsOnFrameNum) = obsOnFrames(randsample(length(obsOnFrames), obsOnFrameNum));
        
        % get obs visible while returning home frames
        obsOffFrames = find(obsPosInterp>lims(1) & obsPosInterp<lims(2) & ~isObsOn);
        frames(end+1:end+obsOffFrameNum) = obsOffFrames(randsample(length(obsOffFrames), obsOffFrameNum));
        
        % get frames near reward delivery
        isNearReward = zeros(size(frameTimes));
        isNearReward(knnsearch(frameTimes, rewardTimes+s.rewardProximity(1))) = 1;
        isNearReward(knnsearch(frameTimes, rewardTimes+s.rewardProximity(2))) = -1;
        nearRewardFrames = find(cumsum(isNearReward));
        frames(end+1:end+rewardFrameNum) = nearRewardFrames(randsample(length(nearRewardFrames), rewardFrameNum));
        
        % get random frames
        frames(end+1:end+randFrameNum) = randsample(totalFrames, randFrameNum);
        
    end
    frames = sort(frames);
    
    for j = 1:length(frames)
        
        frame = read(vid, frames(j));
        trainingData(ind).session = sessions{i};
        trainingData(ind).frameNum = frames(j);
        trainingData(ind).frame = frame;
        trainingData(ind).includeFrame = false;
        ind = ind + 1;
    end
    
end

fprintf('%s: all done!\n', sessions{i})


