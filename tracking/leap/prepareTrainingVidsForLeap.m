function prepareTrainingVidsForLeap(sessions, vidDir, framesPerVid, outputDims)

% given a list of sessions, makes a series of h5 files, each containing a subset of frames from the video in each session
% each frame has top and bot views concatenated, and is resized to outputDims
% these h5 files can be used to create training sets in LEAP

% settings
deflation = 1;
obsLims = [.32 .39]; % !!! this is a hack, because obsPositions are unreliable and vary session to session // perhaps if I normalize by end position for each trial like i used to? // also, this will break if length of track changes

% load video
for i = 1:length(sessions)
    
    % load session data
    tic
    session = sessions{i};
    fprintf('%s: converting to h5... ', session)
    vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
    vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
    load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
        'frameTimeStamps', 'obsOnTimes', 'obsOffTimes', 'obsPositions', 'obsTimes')
    obsPositionsInterped = interp1(obsTimes, obsPositions, frameTimeStamps);

    % create h5 file
    h5create([vidDir session '.h5'], '/box', [outputDims(1) outputDims(2) 1 framesPerVid], 'Chunksize', [outputDims(1) outputDims(2) 1 2], ...
        'Datatype', 'uint8', 'deflate', deflation);
    
    % find frames with obs in it
    isObsInFrame = false(1, length(frameTimeStamps));
    for j = 1:length(obsOnTimes)
        trialObsInFrameBins = frameTimeStamps>obsOnTimes(j) & frameTimeStamps<obsOffTimes(j) & ...
                              obsPositionsInterped>obsLims(1) & obsPositionsInterped<obsLims(2);
        isObsInFrame(trialObsInFrameBins) = true;
    end
    obsInFrameInds = find(isObsInFrame);
    
    
    % convert to h5
    startInd = randi(length(obsInFrameInds)-framesPerVid);
    obsInFrameInds = obsInFrameInds(startInd:startInd+framesPerVid-1);
    
    for j = 1:framesPerVid
        frameTop = read(vidTop, obsInFrameInds(j));
        frameBot = read(vidBot, obsInFrameInds(j));
        frame = imresize(cat(1,frameTop,frameBot), outputDims);
        h5write([vidDir session '.h5'], '/box', frame(:,:,1), [1 1 1 j], [outputDims(1) outputDims(2) 1 1]);
%         disp(j/framesPerVid)
    end
    fprintf('finished in %.1f minutes\n', toc/60)
end
disp('all done!')











