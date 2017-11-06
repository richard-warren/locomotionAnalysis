
% preliminary spike analysis
spikeAnalysis('C:\Users\Rick\Google Drive\columbia\obstacleData\sessions')


% generate plots
obsAvoidance2('run3', 'obsHgtTest')
obsAvoidance2('run4', 'obsHgtTest')
obsAvoidance2('run5', 'obsHgtTest')


% make simple video
trialProportion = .1;
makeVid('171105_000', [.25 .445], .1, trialProportion);
makeVid('171105_001', [.25 .445], .1, trialProportion);
makeVid('171105_002', [.25 .445], .1, trialProportion);



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





