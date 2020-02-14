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