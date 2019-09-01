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

% settings
session = '190814_002';  % .8 acceleration
% session = '190820_001';  % .6 acceleration
% session = '190822_002';  % .6 acceleration
% session = '190823_002';  % .6 acceleration
% session = '190821_000';  % .4 acceleration
% session = '190827_002';  % .8 acceleration // 1.8 current

load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'obsTracking')

figure;
thresholds = linspace(0,1,100);
percentBadTracking = nan(1,length(thresholds));
for i = 1:length(thresholds)
    percentBadTracking(i) = nanmean([obsTracking.percentBadTracking]>thresholds(i));
end

figure; plot(thresholds, percentBadTracking);
xlabel('percent threshold')
ylabel('percent trials with good tracking')

fprintf('%s: %.2f of trials with poor obstalce tracking...\n', session, nanmean([obsTracking.percentBadTracking]>.15))


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

close all; figure;
logPlotRick([flat.obsHgt]*1000, [flat.anyTouches], ...
    {'colors', 'hsv', 'conditions', [flat.paw], 'xlabel', 'obstacle height', 'ylabel', 'contact rate', 'plotMice', false, ...
     'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'conditionNames', {'LH', 'LF', 'RF', 'RH'}, ...
     'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'computeVariance', false, 'mouseNames', {flat.mouse}})
 
%% light spot diam vs fiber size and distance

distances = [2 4 6 8 10 12 16 20];
diams = [0.27 0.4 0.64 1.1 1.2 1.63 2.18 2.7];

distances = [4 6 8 10 12 16 20];
diamsOnSkull = [1.54 2.23 2.92 3.57 4.5 5.82 7.38];

figure; plot(diamsOnSkull, distances)
xlabel('diam')
ylabel('distance')

fit = polyfit(diamsOnSkull, distances, 1);

fitX = 0:.01:10;
fitY = polyval(fit, fitX);

hold on; plot(fitX, fitY);

polyval(fit, .4)


%% fraction of trials on which light terminates prematurely

flat = flattenData(data, {'obsOffTimes', 'earlyOptoTermination', 'session', 'optoPower'});
% flat = flat([flat.optoPower]>.9)












