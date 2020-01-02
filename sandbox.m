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