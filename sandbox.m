%% demonsrate face tracking stuff

load(fullfile(getenv('OBSDATADIR'), 'sessions', '999999_999', 'runAnalyzed.mat'), ...
    'wheelTimes', 'wheelPositions', 'lickTimes', 'whiskerAngles', 'frameTimeStampsWisk', 'targetFs', 'rewardTimes')
wheelVel = getVelocity(wheelPositions, .05, targetFs);
whiskerAngle = fillmissing(whiskerAngle, 'pchip');
% whiskerAngle(whiskerAngle<-120) = nan;
%%
close all; figure('color', 'white', 'position', [286.00 456.00 1437.00 360.00]); hold on
yLims = [-2 2];

plot(wheelTimes, zscore(wheelVel))
plot(frameTimeStampsWisk, zscore(whiskerAngle));
plot([lickTimes, lickTimes], yLims, 'color', [.6 .6 .6])
plot([rewardTimes, rewardTimes], yLims, 'color', 'blue', 'LineWidth', 2)


% legend('velocity', 'whisker angle')
set(gca, 'xlim', [288.7246  304.5106], 'ylim', yLims, 'visible', 'off')


%% play around with whisker angle computation

% locationsWisk = readtable('Z:\loco\obstacleData\sessions\999999_999\trackedFeaturesRaw_wisk.csv');
% load('Z:\loco\obstacleData\sessions\999999_999\runAnalyzed.mat', 'lickTimes', 'frameTimeStampsWisk')

% settings
confidenceThresh = .5;
medianFiltering = 5;
smoothing = 5;
angleLims = [-120 -60];

pad = [median(locationsWisk.wisk_pad) median(locationsWisk.wisk_pad_1)];  % location of whisker pad
valid = locationsWisk.wisk_caudal_2>confidenceThresh | locationsWisk.wisk_rostral_2>confidenceThresh;

c = [medfilt1(locationsWisk.wisk_caudal,medianFiltering), medfilt1(locationsWisk.wisk_caudal_1,medianFiltering)];
r = [medfilt1(locationsWisk.wisk_rostral,medianFiltering), medfilt1(locationsWisk.wisk_rostral_1,medianFiltering)];

avg = mean(cat(3,c,r), 3);
avg = avg - pad;
avg(:,2) = -avg(:,2);
angles = rad2deg(atan2(avg(:,2), avg(:,1)));
angles = smooth(angles, smoothing);
% interpolation?

angles(angles<angleLims(1) | angles>angleLims(2)) = nan;


% close all; figure('position', [186.00 129.00 1601.00 849.00]); hold on
% plot(frameTimeStampsWisk, r(:,1))
% plot(frameTimeStampsWisk, c(:,1))
% plot(frameTimeStampsWisk, a(:,1), 'LineWidth', 2)
% plot([lickTimes, lickTimes], get(gca, 'YLim'), 'color', [.6 .6 .6])

% close all; figure('position', [186.00 129.00 1601.00 849.00]); hold on
% t = randi(size(a,1)); close all; figure; plot([0 a(t,1)], [0 a(t,2)]); set(gca, 'YDir', 'normal'); daspect([1 1 1])
% disp(angles(t))

close all; figure('position', [186.00 129.00 1601.00 849.00]); hold on
angles(~valid) = nan;
plot(frameTimeStampsWisk, angles)
plot([lickTimes, lickTimes], get(gca, 'YLim'), 'color', [.6 .6 .6])
set(gca, 'xlim', [445.1172  570.2924])


%% try peak finding algorithms on lick signal

locationsWisk = readtable('Z:\loco\obstacleData\sessions\999999_999\trackedFeaturesRaw_wisk.csv');

conf = .5;
valid = locationsWisk.tongue_2 > conf;

raw = locationsWisk.tongue_1;
sig = smooth(raw, 5);

sig(~valid) = nan;
raw(~valid) = nan;

lims = prctile(raw, [50 99]);
valid = sig<lims(2);

sig(~valid) = nan;
raw(~valid) = nan;


close all; figure('position', [186.00 129.00 1601.00 849.00]);

plot(raw, 'color', 'blue'); hold on;
plot(sig, 'LineWidth', 2); hold on;


set(gca, 'xlim', [178704 179657])

figure;
findpeaks(sig, 'MinPeakDistance', 10, 'MinPeakHeight', lims(1))

% check that: confidence high // peak above AND BELOW thresh


%% test hildebrand plots

sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'baselineNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);  % remove not included sessions and empty lines
sessionInfo = sessionInfo(1:6,:);  % only analyze a handful of sessions for testing purposes

plotHildebrands(sessionInfo, 'colors', stepColors, 'stepPercentiles', [5 95], 'plotMice', false)


%% check obs tracking confidence

session = '191120_000';
data = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv'));



%% compare whisker trim videos in old and new experiments

% settings
vids = 5;

info = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'sensoryDependenceNotes');
oldSessions = info.session(strcmp(info.whiskers, 'none'));

info = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'senLesionNotes');
newSessions = info.session(strcmp(info.condition, 'noWisk') & [info.include]==1);


fullfile(getenv('OBSDATADIR'), 'editedVid', 'whiskerTrimComparison', seesions{i})
makeVidWisk(vidName, session, obsPosRange, playBackSpeed, trialProportion, trialLabels, trialInds)




%% summarize step lengths vs speed

fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')
flat = flattenData(data, {'mouse', 'session', 'trial', 'paw', 'controlStepLength'});

figure('color', 'white'); hold on
for paws = [1,4; 2 3]'
    histogram([flat(ismember([flat.paw], paws)).controlStepLength] * 1000, 'Normalization', 'probability')
end
set(gca, 'YTick', [], 'xlim', [0 160])
xlabel('step length (mm)')

%% single session analysis

session = '180715_004';
% load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'kinData.mat'), 'kinData');




%% check out which trials are being excluded based on paw height

fprintf('loading data... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', ['mtc_muscimol_data.mat']), 'data'); fprintf('data loaded!\n')
flat = flattenData(data, {'mouse', 'session', 'trial', 'paw', 'isValidZ', 'stepOverMaxHgt', 'obsHgt', 'modPawKin'});
sessions = unique({flat.session});

%%
fprintf('INVALID-Z TRIALS:\n')
fprintf('-----------------\n')
for i = 1:length(sessions)
    fprintf('%s:  ', sessions{i})
    bins = strcmp({flat.session}, sessions{i}) & ~[flat.isValidZ] & [flat.paw]==2;
    fprintf('%i ', unique([flat(bins).trial]))
    fprintf('\n')
end



%% estimate tissue shrinkage landmarks
mmWidth = 11;
im = imread('Y:\obstacleData\papers\hurdles_paper1\figures\histology\rostalHippocampus_coronal.PNG');

mmPerPix = mmWidth / size(im,2);
figure; imshow(im)
ventricalScale = [626-266] * mmPerPix;
baseScale = [569-331] * mmPerPix;

%% test load times for large .mat files, FML!!!

tic; fprintf('\n\n-----SPEED TESTS-----\n')

% small data
tic; load('D:\matlabData\mtc_muscimol_data.mat'); fprintf('loading small data, solid state: %.2f\n', toc)
tic; load('Y:\obstacleData\matlabData\mtc_muscimol_data.mat'); fprintf('loading small data, network: %.2f\n', toc)
tic; m = matfile('D:\matlabData\mtc_muscimol_data.mat'); fprintf('loading small data, solid state, mafile(): %.2f\n', toc)

% large data
tic; load('D:\matlabData\senLesion_data.mat'); fprintf('loading large data, solid state: %.2f\n', toc)
tic; load('Y:\obstacleData\matlabData\senLesion_data.mat'); fprintf('loading large data, network: %.2f\n', toc)
tic; m = matfile('Y:\obstacleData\matlabData\senLesion_data.mat'); fprintf('loading large data, network, matfile(): %.2f\n', toc)

%% saving
tic; save('D:\matlabData\temp.mat', 'data', '-v7.3', '-nocompression'); fprintf('saving large data, solid state, no compression: %.2f\n', toc)
tic; load('D:\matlabData\temp.mat'); fprintf('loading large data, solid state, no compression: %.2f\n', toc)
tic; load('D:\matlabData\senLesion_data.mat'); fprintf('loading large data, solid state, with compression: %.2f\n', toc)


%% recompute frameTimes using new script

sessions = selectSessions;

%%

for i = 1:length(sessions)

    

    
    file = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'run.csv');
    if exist(file, 'file')    
        
        s1 = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
        'frameTimeStamps', 'frameTimeStampsWisk', 'exposure');
        s2 = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'run.mat'), ...
            'exposure');

        camMetadata = dlmread(file); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
        frameCounts = camMetadata(:,2);
        timeStampsFlir = timeStampDecoderFLIR(camMetadata(:,3));
        
        % recompute frameTimeStamps
        frameTimeStamps = getFrameTimes(s2.exposure.times, timeStampsFlir, frameCounts, sessions{i});
        bins = ~isnan(s1.frameTimeStamps);  % only compare these bins
        sum(frameTimeStamps(bins)~=s1.frameTimeStamps(bins));
        
        close all; figure('position', [2402.00 276.00 560.00 420.00]);
        plot(s1.frameTimeStamps, 'lineWidth', 2); hold on; plot(frameTimeStamps)
        differInds = find(s1.frameTimeStamps~=frameTimeStamps & bins);
        scatter(differInds, s1.frameTimeStamps(differInds))
        fprintf('differing frames: %i\n', sum(frameTimeStamps(bins)~=s1.frameTimeStamps(bins)))
        keyboard
        
    end
    
    
end