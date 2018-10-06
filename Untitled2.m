sessions = {'180923_000', '180924_000', '180925_000', '180926_000', '180927_000', '180928_000', '180929_000', '180930_000', '181001_000', '181002_000', '181003_000', '181004_000', '180923_001', '180924_001', '180925_001', '180926_001', '180927_001', '180928_001', '180929_001', '180930_001', '181001_001', '181002_001', '181003_001', '181004_001'};
% sessions = {'180924_001', '180925_001', '180926_001', '180927_001', '180928_001', '180929_001', '180930_001', '181001_001', '181002_001', '181003_001', '181004_001'};

%%
for i = 1:length(sessions)
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
        'obsOnTimes', 'obsOffTimes',  'wiskContactFrames', 'frameTimeStampsWisk');
    
    wiskContactFramesTemp = wiskContactFrames;
    wiskContactFramesTemp(wiskContactFramesTemp<1) = 1; % this version is used to index into frameTimeStampsWisk, which can only be done if -1 values are replaced
    validBins = wiskContactFrames>0 & ...
                frameTimeStampsWisk(wiskContactFramesTemp)>obsOnTimes & ...
                frameTimeStampsWisk(wiskContactFramesTemp)<obsOffTimes;
    
%     fprintf('%s: %.2f of trials are invalid\n', sessions{i}, sum(~validBins)/length(validBins))
    isWiskFrameMissing = max(frameTimeStampsWisk) < max(obsOffTimes);
    fprintf('%i\n', isWiskFrameMissing)
    
    
end


%%

load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
        'obsOnTimes', 'obsOffTimes',  'wiskContactFrames', 'wiskContactFramesConfidences', 'frameTimeStamps', 'frameTimeStampsWisk', 'rewardTimes', 'isLightOn');