% pick sessions

sessions = selectSessions;


%% analyze spike data

disp('starting to analyze sessions...')
for i = 1:length(sessions)    
    spikeAnalysis2(sessions{i}, {'arePawsTouchingObs'});
end
disp('all done!')




%% make video with trials labelled by condition

vidTrialProportion = 0.5;

for i = 1:length(sessions)
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'], 'isLightOn');
    makeVidWisk('', sessions{i}, [-.05 .1], .15, vidTrialProportion, {'OFF', 'ON'}, isLightOn+1);
end


%% plot baseline

% mice = {'run6', 'run7', 'run8', 'sen1', 'sen2', 'sen3', 'sen4', 'sen5', 'sen6', 'mtc1', 'mtc2', 'mtc3', 'mtc4', 'mtc5', 'mtc6', 'den2', 'den4', 'den5'};
mice = {'sen1', 'sen2', 'sen3', 'sen4', 'sen5', 'sen6', 'mtc1', 'mtc2', 'mtc3', 'mtc4', 'mtc5', 'mtc6', 'den2', 'den4', 'den5'};
baselineSummary(mice);

%% plot learning progress
% mice = {'sen2', 'sen3', 'sen4', 'sen5', 'sen6'};
mice = { 'mtc4', 'mtc5', 'mtc6'};
makeSpeedAndAvoidanceFigs(mice);