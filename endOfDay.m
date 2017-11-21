% spike analysis, make videos, make figures

% settings
sessionDirs = uigetdir2('C:\Users\rick\Google Drive\columbia\obstacleData\sessions', 'select folders to analyze');
dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\sessions\';
trialProportion = .2;

parfor i=1:length(sessionDirs)
    
    nameInd = find(sessionDirs{i}=='\',1,'last');
    
    % spike analysis
    spikeAnalysis(sessionDirs{i}(1:nameInd), sessionDirs{i}(nameInd+1:end), {'obsPixPositions'});
    
    % make video
%     makeVid(sessionDirs{i}(nameInd+1:end), [.25 .445], .1, trialProportion);

end

delete(gcp); % delete parallel pool

% generate plots
obsAvoidanceLight2('run3', {'obsTestLight2'})
obsAvoidanceLight2('run5', {'obsTestLight2'})



%% make video with trials labelled by condition
for j = 1:length(sessionDirs)
    
    % load session data
    nameInd = find(sessionDirs{j}=='\',1,'last');
    session = sessionDirs{j}(nameInd+1:end);
    load([dataDir session '\runAnalyzed.mat'],...
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

    makeVid(session, [.25 .445], .1, trialProportion, {'OFF', 'ON'}, trialConditions);
%     makeVid(session, [0 .445], .5, trialProportion, {'OFF', 'ON'}, trialConditions);
end





