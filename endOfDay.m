
% preliminary spike analysis
spikeAnalysis('C:\Users\Rick\Google Drive\columbia\obstacleData\sessions')


% generate plots
obsAvoidanceLight2('run3', 'obsTestLight')
obsAvoidanceLight2('run4', 'obsTestLight')
obsAvoidance2('run5', 'obsTestBr')



%% make simple video

makeVid('171102_002', [.25 .445], .1);



%% make video with trials labelled by condition

% settings
session = '171102_000';
dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';
conditions = {'light off', 'light on'};

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


makeVid(session, [.25 .445], .1, {'OFF', 'ON'}, trialConditions);





