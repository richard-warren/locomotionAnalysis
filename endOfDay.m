% pick sessions

sessions = selectSessions;



%% run DeepLabCut analysis
dlcPath = 'C:\Users\rick\Desktop\github\DeepLabCut';
for i = 1:length(sessions)
    try
        fprintf('%s: conducting DeepLabCut analysis...\n', sessions{i})
        tic; [~,~] = system(['cd ' dlcPath ' && batchDLC.bat ' sessions{i}]);
        fprintf('%s: DeepLabCut analysis finished in %.1f hours\n', sessions{i}, toc/60/60)
    catch
        fprintf('%s: problem with DeepLabCut analysis!\n', sessions{i})
    end
end



%% analyze spike data
disp('starting to analyze sessions...')
for i = 1:length(sessions)
%     try
        spikeAnalysis2(sessions{i});
%     catch
%         fprintf('%s: problem with spike analysis!\n', sessions{i})
%     end
end
disp('all done!')



%% make video with trials labelled by condition
vidTrialProportion = 0.1;

for i = 1:length(sessions)
    try
        load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'], 'isLightOn');
        makeVidWisk('', sessions{i}, [-.05 .1], .15, vidTrialProportion, {'OFF', 'ON'}, isLightOn+1);
    catch
        fprintf('%s: problem editing video!\n', sessions{i})
    end
end


%% plot baseline

% mice = {'run6', 'run7', 'run8', 'sen1', 'sen2', 'sen3', 'sen4', 'sen5', 'sen6', 'mtc1', 'mtc2', 'mtc3', 'mtc4', 'mtc5', 'mtc6', 'den2', 'den4', 'den5'};
mice = {'sen1', 'sen2', 'sen3', 'sen4', 'sen5', 'sen6', 'mtc1', 'mtc2', 'mtc3', 'mtc4', 'mtc5', 'mtc6', 'den2', 'den4', 'den5'};
baselineSummary(mice);

%% plot learning progress
senMice = {'sen2', 'sen3', 'sen4', 'sen5', 'sen6'};
mtcMice = {'mtc1', 'mtc2', 'mtc3', 'mtc4', 'mtc5', 'mtc6'};
senData = makeSpeedAndAvoidanceFigs(senMice);
mtcData = makeSpeedAndAvoidanceFigs(mtcMice);


