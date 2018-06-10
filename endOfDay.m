% spike analysis, make videos, make figures

% settings
sessionDirs = uigetdir2([getenv('OBSDATADIR') 'sessions\'], 'select folders to analyze');
% miceToPlot = {'run6', 'run7', 'run8'};
vidTrialProportion = .15;
%%
parfor i = 1:length(sessionDirs)
    % spike analysis
    nameInd = find(sessionDirs{i}=='\',1,'last');
    spikeAnalysis2(sessionDirs{i}(1:nameInd), sessionDirs{i}(nameInd+1:end));
end
disp('all done!')
%% delete(gcp); % delete parallel pool

% generate plots
% baselineAnalysis('run6')
% baselineAnalysis('run7')
% baselineAnalysis('run8')
% obsAvoidanceLearning('run6', {'obsBr'})
% obsAvoidanceLearning('run7', {'obsBr'})
% obsAvoidanceLearning('run8', {'obsBr'})
% obsAvoidanceLearningSummary(miceToPlot)
% speedOverTimeSummary(miceToPlot)
% plotAllSessionSpeeds({'obsBr', 'obsCompareOnOff', 'obsCompareOff', 'obsWiskMrk'}, [-.08 .08] + 0.3820)
% plotAllSessionSpeeds({'obsWisk'}, [-.08 .08] + 0.3820)


% make video with trials labelled by condition
for i = 1:length(sessionDirs)
    
    % load session data
    nameInd = find(sessionDirs{i}=='\',1,'last');
    session = sessionDirs{i}(nameInd+1:end);
    load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'isLightOn');
    makeVidWisk('', session, [-.1 .1], .1, vidTrialProportion, {'OFF', 'ON'}, isLightOn+1);
end





