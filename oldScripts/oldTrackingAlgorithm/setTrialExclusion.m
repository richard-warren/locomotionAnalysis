function setTrialExclusion(session, excludedTrials)

% given a session, creates isExcluded.mat variable in tracking folder for session with one boolean entry per frame
% recording whether that frame belongs to a trial that should be excluded
% note: if this file already exists it is overwritten! omg


load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBot.mat'], 'locations')
isExcluded = ismember(locations.trialIdentities, excludedTrials);
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\isExcluded.mat'], 'isExcluded', 'excludedTrials')