% spike preliminary analysis and video making!

% settings
sessionDirs = uigetdir2([getenv('OBSDATADIR') 'sessions\'], 'select folders to analyze');



%%
parfor i = 1:length(sessionDirs)
    % spike analysis
    nameInd = find(sessionDirs{i}=='\',1,'last');
    spikeAnalysis2(sessionDirs{i}(1:nameInd), sessionDirs{i}(nameInd+1:end));
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





