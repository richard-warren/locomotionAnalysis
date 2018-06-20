% spike preliminary analysis and video making!

% settings
sessions = selectSessions;



%%
for i = 1:length(sessions)
    spikeAnalysis2(sessions{i});
end
disp('all done!')



%% make video with trials labelled by condition
vidTrialProportion = 1;

for i = 1:length(sessions)
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'], 'isLightOn');
    makeVidWisk('', sessions{i}, [-.05 .1], .15, vidTrialProportion, {'OFF', 'ON'}, isLightOn+1);
end





