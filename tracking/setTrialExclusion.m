function setTrialExclusion(session, trials)

% given a session, creates isExcluded.mat variable in tracking folder for session with one boolean entry per frame
% recording whether that frame belongs to a trial that should be excluded
% note: if this file already exists it is overwritten! omg

% temp
% session = '180122_001';
% trials = [107];
% excludeTrials = [1];

load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBot.mat'], 'locations')
isExcluded = false(length(locations.isAnalyzed), 1);

for i = 1:length(trials)
    trialInds = (locations.trialIdentities==trials(i));
    isExcluded(trialInds) = true;
end

save([getenv('OBSDATADIR') 'sessions\' session '\tracking\isExcluded.mat'], 'isExcluded')