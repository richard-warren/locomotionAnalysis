% spike analysis, make videos, make figures

% settings
sessionDirs = uigetdir2([getenv('OBSDATADIR') 'sessions\'], 'select folders to analyze');
trialProportion = .15;

parfor i = 1:length(sessionDirs)
    
    nameInd = find(sessionDirs{i}=='\',1,'last');
    
    % spike analysis
    spikeAnalysis(sessionDirs{i}(1:nameInd), sessionDirs{i}(nameInd+1:end));
    
    % make video
%     makeVid(sessionDirs{i}(nameInd+1:end), [.25 .445], .1, trialProportion);

end

% delete(gcp); % delete parallel pool

% generate plots
% baselineAnalysis('run6')
% baselineAnalysis('run7')
% baselineAnalysis('run8')
% obsAvoidanceLight('run6', {'obsNoBr', 'obsBr'})
% obsAvoidanceLight('run7', {'obsNoBr', 'obsBr'})
% obsAvoidanceLight('run8', {'obsNoBr', 'obsBr'})
obsAvoidanceLearningSummary({'run6', 'run7', 'run8'})



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

    makeVidWisk(session, [.25 .445], .1, trialProportion, {'OFF', 'ON'}, trialConditions);
end





