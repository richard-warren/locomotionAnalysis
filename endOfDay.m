% spike analysis, make videos, make figures

% settings
sessionDirs = uigetdir2([getenv('OBSDATADIR') 'sessions\'], 'select folders to analyze');
% miceToPlot = {'run6', 'run7', 'run8'};
% vidTrialProportion = .15;
%%
for i = 1:length(sessionDirs)
    
    % spike analysis
    nameInd = find(sessionDirs{i}=='\',1,'last');
    spikeAnalysis(sessionDirs{i}(1:nameInd), sessionDirs{i}(nameInd+1:end), {'isLightOn'});
    
end

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
for j = 1:length(sessionDirs)
    
    % load session data
    nameInd = find(sessionDirs{j}=='\',1,'last');
    session = sessionDirs{j}(nameInd+1:end);
    load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'],...
         'obsOnTimes', 'obsOffTimes',...
	     'obsLightOnTimes', 'obsLightOffTimes');

    % find obstacle light on trial inds
    trialConditions = zeros(1,length(obsOnTimes));

    for i =1:length(obsOnTimes)

        if min(abs(obsOnTimes(i) - obsLightOnTimes)) < .5
            trialConditions(i) = 2; % light is on
        else
            trialConditions(i) = 1; % light is off
        end
    end

    makeVidWisk('', session, [-.1 .1], .1, vidTrialProportion, {'OFF', 'ON'}, trialConditions);
end





