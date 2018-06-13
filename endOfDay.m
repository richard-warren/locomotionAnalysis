% spike preliminary analysis and video making!

% settings
sessions = selectSessions;



%%
for i = 1:length(sessions)
    spikeAnalysis2(sessions{i});
end
disp('all done!')



%% make video with trials labelled by condition
vidTrialProportion = .15;

for i = 1:length(sessionDirs)
    
    % load session data
    nameInd = find(sessionDirs{i}=='\',1,'last');
    session = sessionDirs{i}(nameInd+1:end);
    load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'isLightOn');
    makeVidWisk('', session, [-.1 .1], .1, vidTrialProportion, {'OFF', 'ON'}, isLightOn+1);
end





