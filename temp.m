%% measure times of trials and other little thangs

fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

flat = flattenData(data, {'wiskContactTimes', 'obsOnTimes', 'obsOffTimes', 'lightOnTimes', 'isLightOn'});
flat = flat(~[flat.isLightOn]);
medTimeUntilContact = nanmedian([flat.wiskContactTimes] - [flat.obsOnTimes]);
medObsOnTime = nanmedian([flat.obsOffTimes] - [flat.obsOnTimes]);
medTrackingTime = nanmedian([flat.lightOnTimes] - [flat.obsOnTimes]);


%% measure time between rewards
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'baselineNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);

interRewardIntervals = nan(1, height(sessionInfo));
for i = 1:height(sessionInfo)
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessionInfo.session{i}, 'runAnalyzed.mat'), 'rewardTimes')
    interRewardIntervals(i) = nanmedian(diff(rewardTimes));
end

figure; histogram(interRewardIntervals)


%% check obstacle tracking...
% use as criterion number of samples where vel deviation is above some
% threshold

% settings
session = '190401_004';
velTime = .01;
obsOnBuffer = .2;
velTolerance = .02;


% initializations
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'obsOnTimes', 'obsOffTimes', 'wheelPositions', 'wheelTimes', 'obsPositions', 'obsTimes', 'targetFs')
wheelVel = getVelocity(wheelPositions, velTime, targetFs);
obsVel = getVelocity(obsPositions, velTime, targetFs);

% get trial data
data = struct();
rowInd = 1;
for i = 1:length(obsOnTimes)
    
    wheelBins = wheelTimes>(obsOnTimes(i)+obsOnBuffer) & wheelTimes<obsOffTimes(i);
    obsBins = obsTimes>obsOnTimes(i) & obsTimes<obsOffTimes(i);
    wheelVelTrial = wheelVel(wheelBins);
    obsVelTrial = obsVel(obsBins);
    times = wheelTimes(wheelBins);
    
    if ~isequal(times, obsTimes(obsBins))
        obsVelTrial = interp1(obsTimes(obsBins), obsVelTrial, times);
    end
    
    data(rowInd).wheelVel = wheelVelTrial;
    data(rowInd).obsVel = obsVelTrial;
    data(rowInd).times = times;
    data(rowInd).percentBadTracking = nanmean(abs(wheelVelTrial-obsVelTrial) > velTolerance);
    rowInd = rowInd + 1;
end

close all; figure; hold on;
plot([data.wheelVel]);
plot([data.obsVel])

%% compute percentage bad tracking for all past sessions!

sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'whiskerTrimNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);
% !!! start with 54 for lesionNotes
for i = 1:height(sessionInfo); spikeAnalysis2(sessionInfo.session{i}, 'obsTracking'); end
disp('all done!')


%% paw contact rate as function of obs height


fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

flat = flattenData(data, {'obsHgt', 'numTouchFrames', 'paw', 'isPawSuccess', 'mouse'});
anyTouches = num2cell([flat.numTouchFrames]>0);
[flat.anyTouches] = deal(anyTouches{:});

%%

close all; figure;
logPlotRick([flat.obsHgt]*1000, [flat.anyTouches], ...
    {'colors', 'hsv', 'conditions', [flat.paw], 'xlabel', 'obstacle height', 'ylabel', 'contact rate', 'plotMice', false, ...
     'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'conditionNames', {'LH', 'LF', 'RF', 'RH'}, ...
     'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'computeVariance', false, 'mouseNames', {flat.mouse}})













