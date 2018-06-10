function trainingData = createTrainingDataStruct(sessions, frameNum, obsPortion)

% settings
obsLims = [.32 .39]; % !!! this is a hack, because obsPositions are unreliable and vary session to session // perhaps if I normalize by end position for each trial like i used to? // also, this will break if length of track changes
minFrame = 250*5*60; % 4 minutes at 250 fps

% initializations
sessionFrameNums = ones(1, length(sessions)) * round(frameNum / length(sessions));
sessionFrameNums(end) = sessionFrameNums(end) + (frameNum-sum(sessionFrameNums));
trainingData(frameNum,1) = struct();

for i = 1:length(sessions)
    
    % load vid and sesion data
    fprintf('loading video: %s\n', sessions{i})
    vid = VideoReader([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runBot.mp4']);
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'], ...
        'frameTimeStamps', 'obsTimes', 'obsPositions', 'obsOnTimes', 'obsOffTimes')
    obsPositionsInterped = interp1(obsTimes, obsPositions, frameTimeStamps);
    
    % get number of obs and obs frames for this session
    obsNum = round(sessionFrameNums(i)*obsPortion); % number of frames to get with obstacle
    noObsNum = sessionFrameNums(i) - obsNum;
    
    % find frames with obstacles
    obsBins = false(1, length(frameTimeStamps));
    for j = 1:length(obsOnTimes)
        trialObsInFrameBins = frameTimeStamps>obsOnTimes(j) & frameTimeStamps<obsOffTimes(j) & ...
                              obsPositionsInterped>obsLims(1) & obsPositionsInterped<obsLims(2);
        obsBins(trialObsInFrameBins) = true;
    end
    obsInds = find(obsBins);
    noObsInds = find(~obsBins & 1:length(frameTimeStamps)>minFrame);
    
    % find frames without obstacles
    obsFrames = obsInds(randperm(length(obsInds), obsNum)); % inds of randomly selected obstacle frames
    noObsFrames = noObsInds(randperm(length(noObsInds), noObsNum));
    
    frames = sort([obsFrames noObsFrames]); % all obstacle and non-obstacle frames
    frames = mat2cell(frames, 1, ones(sessionFrameNums(i),1));
    
    
    % get struct inds for session
    if i==1; startInd = 1; else; startInd = sum(sessionFrameNums(1:i-1))+1; end
    endInd = startInd + sessionFrameNums(i)-1;
    
    % insert session name and frame numbers into struct
    sessionName = cell(1, sessionFrameNums(i)); sessionName(:) = {sessions{i}};
    [trainingData(startInd:endInd).session] = sessionName{:};
    [trainingData(startInd:endInd).frameNum] = frames{:};
    
end

falseVector = num2cell(false(1, frameNum));
[trainingData(:).includeFrame] = falseVector{:};



