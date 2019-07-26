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

session = '190308_000';
trials = 55:65; % 5 14 43
slowSpeed = .15;
makeMatrixVid2(fullfile(getenv('OBSDATADIR'), 'editedVid', session), session, trials, false, slowSpeed)



%% make video showing how we can monitor wisk and paw contacts for talks

session = '190327_003';
obsPosRange = [-.05 .1];
playbackSpeed = .1;
trials = [1 6 11 20 31 35 41 46 51 56];


makeVidWisk(fullfile(getenv('OBSDATADIR'), 'editedVid', session), session, obsPosRange, playbackSpeed, trials);

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



%% make videos with neural firing for obstacles and rewards

% settings
timePrePost = [-1 1.5];

ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
ephysInfo = ephysInfo(strcmp(ephysInfo.spikesSorted, 'yes') & ephysInfo.numGoodUnits>0, :);
close all

for i = 1:height(ephysInfo)
    
    % load session data
    s1 = load(fullfile(getenv('OBSDATADIR'), 'sessions', ephysInfo.session{i}, 'runAnalyzed.mat'));
    s2 = load(fullfile(getenv('OBSDATADIR'), 'sessions', ephysInfo.session{i}, 'neuralData.mat'));
    
    for j = s2.unit_ids'
        try
            fprintf('\n%s unit %i\n', ephysInfo.session{i}, j)
            
            % obstacle response vid
            fileName = fullfile(getenv('OBSDATADIR'), 'editedVid', 'vidsWithNeurons', 'obstacleResponses', ...
                [ephysInfo.session{i} 'unit' num2str(j) '.avi']);
            timeEpochs = cat(2, s1.obsOnTimes, s1.obsOffTimes);
            makeUnitVid(ephysInfo.session{i}, j, fileName, timeEpochs)

            % reward response vid
            fileName = fullfile(getenv('OBSDATADIR'), 'editedVid', 'vidsWithNeurons', 'rewardResponses', ...
                [ephysInfo.session{i} 'unit' num2str(j) '.avi']);
            timeEpochs = s1.rewardTimes + timePrePost;
            makeUnitVid(ephysInfo.session{i}, j, fileName, timeEpochs)
        catch
            fprintf('PROBLEM WITH %s unit %i\n', ephysInfo.session{i}, j)
        end
    end
end





%% make video with neural firing for reward times

session = '181004_003';
unit_id = 74;
timePrePost = [-1 1.5];


% load reward times
fileName = fullfile(getenv('OBSDATADIR'), 'editedVid', 'vidsWithNeurons', 'rewardResponses', [session 'unit' num2str(unit_id) '.avi']);
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'rewardTimes')
timeEpochs = rewardTimes + timePrePost;

% make video
makeUnitVid(session, unit_id, fileName, timeEpochs)


%% make video of all wheel break trials for a session

% sessions = selectSessions;
sessions = {'190405_000'};
for i = 1:length(sessions)
    wheelBreakTrials = find(getIsWheelBreak(sessions{i}))';
    makeVidWisk(fullfile(getenv('OBSDATADIR'), 'editedVid', 'wheelBreakOnlyVids', [sessions{i} 'WheelBreaks']), ...
        sessions{i}, [-.05 .1], .2, wheelBreakTrials);
end


%% make video comparing success and failure trials

% settings
sessions = {'180801_001', '180805_001', '180809_001', '180812_001', '190227_000'};
maxTrialsToShow = 15;

for i = 1:length(sessions)
    
    % identify success and failure trials
    sessionData = getExperimentData(sessions{i}, {'isTrialSuccess'});
    sessionData = getNestedStructFields(sessionData, {'mouse', 'session', 'trial', 'isTrialSuccess'});
    successTrials = find([sessionData.isTrialSuccess]);
    failTrials = find(~[sessionData.isTrialSuccess]);
    
    % limit to maxTrialsToShow per trial type
    successTrials = successTrials(sort(randperm(length(successTrials), min(maxTrialsToShow, length(successTrials)))));
    failTrials = failTrials(sort(randperm(length(failTrials), min(maxTrialsToShow, length(failTrials)))));
    
    % make vids
    makeVidWisk(fullfile(getenv('OBSDATADIR'), 'editedVid', 'failureTrialComparison', [sessions{i} '_successTrials']), ...
        sessions{i}, [-.05 .1], .2, successTrials);
    makeVidWisk(fullfile(getenv('OBSDATADIR'), 'editedVid', 'failureTrialComparison', [sessions{i} '_failTrials']), ...
        sessions{i}, [-.05 .1], .2, failTrials);    
end

%% make videos with whisker contact to manually estimate reaction times

% settings
experiment = 'baseline';
trialsPerSession = 10;
modPawLims = [5, 20]; % only include trials where contact occurs within modPawLims samples of swing start
minVel = .4;

% initializations
varsToGet = {'mouse', 'session', 'trial', 'preModPawKin', 'modPawKin', 'isBigStep', 'isLightOn', 'modPawContactInd', 'velAtWiskContact'};
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', [experiment 'Notes']);
sessionInfo = sessionInfo(~cellfun(@isempty, sessionInfo.session),:); % remove empty rows, not included sessions, and those without correct brain region
mice = unique(sessionInfo.mouse(~cellfun(@isempty, sessionInfo.session)));
data = getExperimentData(sessionInfo, varsToGet);
flat = getNestedStructFields(data, varsToGet);
flat = flat([flat.modPawContactInd]>=modPawLims(1) & [flat.modPawContactInd]<=modPawLims(2) & ...
            ~[flat.isLightOn] & ...
            [flat.velAtWiskContact]>minVel); % add conditionals here
flat = struct2table(flat);
spreadsheet = array2table(zeros(0,3), 'VariableNames', {'mouse', 'session', 'trial'});
%%
for i = 1:length(mice)
    session = sessionInfo.session{find(strcmp(sessionInfo.mouse, mice{i}),1,'first')};
    inds = find(strcmp(flat.session, session));
    inds = inds(sort(randperm(length(inds), min(trialsPerSession, length(inds)))));
    if ~isempty(inds)
        spreadsheet = cat(1, spreadsheet, flat(inds, 1:3));
        makeVidWisk(fullfile(getenv('OBSDATADIR'), 'editedVid', 'reactionTimeVids', session), ...
            session, [-.05 .1], .15, flat.trial(inds)');
    end
end

file = fullfile(getenv('OBSDATADIR'), 'editedVid', 'reactionTimeVids', [experiment '_reactionTimes.csv']);
if ~exist(file, 'file'); writetable(spreadsheet, file); else; disp('file already exists!'); end
disp('all done')


%% make vid, and tracking spreadsheet, with no obstacles (for hyun soo park)

session = '190222_002';
trialsToInclude = 20;

file = fullfile(getenv('OBSDATADIR'), 'other', 'hyunSooParkData', session);
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'obsPixPositions', 'frameTimeStamps', 'obsOnTimes', 'obsOffTimes')
trackingData = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv'));

% take only trialsToInclude fastest trials
[~, sortInds] = sort(obsOffTimes - obsOnTimes);
% sortInds = flipud(sortInds); % uncomment this line to take the slowest trials!
trials = sort(sortInds(1:trialsToInclude));

vidTop = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
vidBot = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runBot.mp4'));
vidWriter = VideoWriter(file);
set(vidWriter, 'FrameRate', 50)
open(vidWriter)

inds = cell(1,trialsToInclude);

for i = 1:trialsToInclude
    
    % get trial inds
    inds{i} = find(frameTimeStamps>obsOnTimes(trials(i)) & ...
                frameTimeStamps<obsOffTimes(trials(i)) & ...
                obsPixPositions'>(vidTop.Width+20));
    
    % write trial inds to video
    for j = inds{i}'
        frame = cat(1, read(vidTop,j), read(vidBot,j));
        writeVideo(vidWriter, frame);
    end
end
close(vidWriter)

frameInds = cat(1,inds{:});
noObsColBins = ~contains(trackingData.Properties.VariableNames, 'obs') & ...
               ~strcmp(trackingData.Properties.VariableNames, 'Var1');

% change column headings to make easier to understand
for i = find(noObsColBins)
    varEnd = trackingData.Properties.VariableNames{i}(end-1:end);
    if strcmp(varEnd, '_1')
        trackingData.Properties.VariableNames{i} = [trackingData.Properties.VariableNames{i}(1:end-1) 'y'];
    elseif strcmp(varEnd, '_2')
        trackingData.Properties.VariableNames{i} = [trackingData.Properties.VariableNames{i}(1:end-1) 'confidence'];
    else
        trackingData.Properties.VariableNames{i} = [trackingData.Properties.VariableNames{i} '_x'];
    end
end
writetable(trackingData(frameInds, noObsColBins), ...
    fullfile(getenv('OBSDATADIR'), 'other', 'hyunSooParkData', 'deepLabCutOutput.csv'));
disp('all done!')






