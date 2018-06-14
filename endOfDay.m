% spike preliminary analysis and video making!

% settings
sessions = selectSessions;



%%
for i = 1:length(sessions)
    spikeAnalysis2(sessions{i}, {'obsPixPositions'});
end
disp('all done!')



%% make video with trials labelled by condition
vidTrialProportion = .5;

for i = 1:length(sessions)
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'], 'isLightOn');
    makeVidWisk('', sessions{i}, [-.1 .1], .1, vidTrialProportion, {'OFF', 'ON'}, isLightOn+1);
end





