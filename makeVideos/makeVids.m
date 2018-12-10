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



%% make vids with no obstacle, only locomotion (for KineMouse Wheel demonstration)

session = '180125_000';
trialNum = 100;
maxFrames = 10000;
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


makeTrackingVidDLC(session, frameInds'); % make vid with markers


%% simply concatenate top and bot only including frameInds(1:maxFrames)
filename = ['noObsVid' session];
vidTop = VideoReader([getenv('OBSDATADIR') '\sessions\' session '\runTop.mp4']);
vidBot = VideoReader([getenv('OBSDATADIR') '\sessions\' session '\runBot.mp4']);
vidWrite = VideoWriter([getenv('OBSDATADIR') '\editedVid\' filename 'trackingSample.mp4'], 'MPEG-4');
set(vidWrite, 'Quality', 50, 'FrameRate', 100);
open(vidWrite);

for i = 1:maxFrames
    disp(i/maxFrames)
    frame = cat(1, read(vidTop, frameInds(i)), read(vidBot, frameInds(i)));
    writeVideo(vidWrite, frame);
end

close(vidWrite)


%% make matched trials vids for muscimol and lesion experiments


% settings
manipulation = 'lesion';
brainRegion = 'sen';
maxLesionSession = 3;
params = {'vel', 'avgVel', 'avgAngle', 'firstModPaw', 'swingStartDistance', 'obsPos'}; % params to match across trials
tolerances = [.05, .01, 2, .5, .01, .005]; % trials are only matched if they are within this range on a given parameter
maxTrials = 15;




% initializations
load([getenv('OBSDATADIR') 'matlabData\' manipulation 'KinematicData.mat'], 'data');
kinData = data; clear data;
if strcmp(manipulation, 'lesion'); kinData = kinData([kinData.conditionNum]<=maxLesionSession); end
disp([manipulation ' kinematic data loaded!'])

if strcmp(manipulation, 'muscimol'); conditions = {'saline', 'muscimol'};
elseif strcmp(manipulation, 'lesion'); conditions = {'pre', 'post'}; end
controlBins = strcmp({kinData.condition}, conditions{1});
manipBins = strcmp({kinData.condition}, conditions{2});
brainRegionBins = strcmp({kinData.brainRegion}, brainRegion);
mice = unique({kinData(brainRegionBins).mouse});

% get trial bins for each parameter
binIds = nan(length(kinData), length(params));
for i = 1:length(params)
    binEdges = min([kinData.(params{i})]) : tolerances(i) : max([kinData.(params{i})]+tolerances(i));
    [~,~,binIds(:,i)] = histcounts([kinData.(params{i})], binEdges);
end


%% get trial pairs
for i = 3:length(mice)
    
    mouseBins = strcmp({kinData.mouse}, mice{i});
    uniqueTrialTypes = unique(binIds(mouseBins,:), 'rows');
    
    sessions = cell(2,0);
    trials = nan(2,0);
    trialText = cell(2,0);
    
    for j = 1:size(uniqueTrialTypes,1)
        manipInds = find(ismember(binIds, uniqueTrialTypes(j,:), 'rows')' & mouseBins & manipBins);
        controlInds = find(ismember(binIds, uniqueTrialTypes(j,:), 'rows')' & mouseBins & controlBins);
        
        if ~isempty(manipInds) && ~isempty(controlInds)
            sessions(:,end+1) = {kinData(controlInds(1)).session; kinData(manipInds(1)).session};
            trials(:,end+1) = [kinData(controlInds(1)).trial; kinData(manipInds(1)).trial];
            
            control = kinData(controlInds(1));
            manip = kinData(manipInds(1));
            trialText(:,end+1) = ...
                {sprintf('%s, %s, %s, %s, %s, light:%i, trial:%i', control.mouse, manip.brainRegion, control.condition, control.side, control.session, control.isLightOn, manip.trial); ...
                sprintf('%s, %s, %s, %s, %s, light:%i, trial:%i', manip.mouse, manip.brainRegion, manip.condition, manip.side, manip.session, manip.isLightOn, manip.trial)};
            
            if size(sessions,2)>=maxTrials; break; end
        end
    end
    
    fprintf('making video for mouse %s...\n', mice{i})
    makeTrialPairsVid(fullfile(getenv('OBSDATADIR'), 'editedVid', 'matchedTrials', [mice{i} '_' manipulation '_' brainRegion '.mp4']), ...
        sessions, trials, trialText);
end



%% make video with neural firing for reward times

session = '181001_002';
unit_id = 6;
trialDuration = 5;


% load reward times
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'rewardTimes')
timeEpochs = cat(2, rewardTimes, rewardTimes+trialDuration);

% make video
makeUnitVid(session, unit_id, timeEpochs)








