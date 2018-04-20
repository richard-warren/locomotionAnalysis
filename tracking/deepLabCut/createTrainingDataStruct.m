function trainingData = createTrainingDataStruct(sessions, frameNum, obsPortion)




% initializations
sessionFrameNums = ones(1, length(sessions)) * round(frameNum / length(sessions));
sessionFrameNums(end) = sessionFrameNums(end) + (frameNum-sum(sessionFrameNums));
trainingData(frameNum,1) = struct();


for i = 1:length(sessions)
    
    % load vid and sesion data
    fprintf('loading video: %s\n', sessions{i})
    vid = VideoReader([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runBot.mp4']);
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'], 'obsPixPositions')
    minFrame = find(~isnan(obsPixPositions),1,'first');
    
    
    % select random frame numbers
    obsNum = round(sessionFrameNums(i)*obsPortion); % number of frames to get with obstacle
    noObsNum = sessionFrameNums(i) - obsNum;
    
    obsBins = obsPixPositions>0 & obsPixPositions<=vid.Width;
    obsInds = find(obsBins & 1:length(obsPixPositions)>minFrame); % inds of all obstacle frames
    noObsInds = find(~obsBins & 1:length(obsPixPositions)>minFrame);
    
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