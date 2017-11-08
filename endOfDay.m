% spike analysis, make videos, make figures

% settings
sessionDirs = uigetdir2('C:\Users\rick\Google Drive\columbia\obstacleData\sessions', 'select folders to analyze');
trialProportion = .1;

for i=1:length(sessionDirs)
    
    nameInd = find(sessionDirs{i}=='\',1,'last');
    
    % spike analysis
    spikeAnalysis(sessionDirs{i}(1:nameInd), sessionDirs{i}(nameInd+1:end));
    
    % make video
    makeVid(sessionDirs{i}(nameInd+1:end), [.25 .445], .1, trialProportion);

end

delete(gcp); % delete parallel pool

% generate plots
% obsAvoidance2('run3', 'obsHgtTest')
% obsAvoidance2('run4', 'obsHgtTest')
obsAvoidance2('run5', 'obsHgtTest')



%% make video with trials labelled by condition

% settings
session = '171103_000';
dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';
conditions = {'light off', 'light on'};
trialProportion = .1;

% load session data
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





