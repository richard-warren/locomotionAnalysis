%% make wide view setup example view (shows obs coming into view from far off on the right side of the screen)

session = '180913_003';
trials = [7 10 30 50];
makeSetupExampleVid(['setupExample' session], session, trials);



%% make fast and slow matrix vids

session = '180715_004';
trials = 20:30; % 5 14 43
slowSpeed = .05;
makeMatrixVid2(sprintf('matrixVidFast%s', session), session, trials, false, 1)
makeMatrixVid2(sprintf('matrixVidSlow%s', session), session, trials, false, slowSpeed)


%% make slowed down example for tracking

session = '180715_004';
trials = 1:10; % 5 14 43
slowSpeed = .15;
makeMatrixVid2(sprintf('trackingEg%s', session), session, trials, true, slowSpeed)




%% make video showing how we can monitor wisk and paw contacts for talks

session = '180912_004';
obsPosRange = [-.05 .1];
playbackSpeed = .1;
trials = [1 6 11 20 31 35 41 46 51 56];


makeVidWisk([getenv('OBSDATADIR') 'editedVid\talkExamle' session], session, obsPosRange, playbackSpeed, trials);



%% make vid with obstacles

session = '180125_000';
trialNum = 5;


load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'rewardTimes', 'frameTimeStamps', 'obsOnTimes', 'obsOffTimes')

frameInds = [];
trials = randperm(length(obsOnTimes), trialNum);

for i = trials
    
    trialInds = find(frameTimeStamps>obsOnTimes(i) & frameTimeStamps<obsOffTimes(i));
    frameInds = cat(1, frameInds, trialInds);
end


makeTrackingVidDLC(session, frameInds);



%% make vids with no obstacle, only locomotion

session = '180125_000';
trialNum = 3;
obsOffBuffer = 1; % (s)
rewardBuffer = 4;


load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'rewardTimes', 'frameTimeStamps', 'obsPixPositions', 'obsOnTimes', 'obsOffTimes')
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);


frameInds = [];
trials = randperm(length(obsOnTimes)-2, trialNum);

for i = trials
    
    startInd = find(frameTimeStamps>obsOffTimes(i)+obsOffBuffer, 1, 'first'); % right after obs turns off, plus a buffer
    endInd = find(frameTimeStamps>obsOffTimes(i)+obsOffBuffer & ...
        frameTimeStamps<obsOffTimes(i+1) & obsPixPositions'<(vid.Width+20), 1, 'first'); % right before obs get into frame on subsequent trial
   
    % if a reward is in the trial start after the reward time
    rewardInTrialInd = find(rewardTimes>frameTimeStamps(startInd) & rewardTimes<frameTimeStamps(endInd));
    if ~(isempty(rewardInTrialInd))
        startInd = find(frameTimeStamps > rewardTimes(rewardInTrialInd)+rewardBuffer);
    end
    
    frameInds = cat(2, frameInds, startInd:endInd);
end


makeTrackingVidDLC(session, frameInds');