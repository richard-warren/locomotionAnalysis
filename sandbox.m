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