function isExcludedNew = setTrialExclusion(isExcluded, trialIdentities, trials, excludeTrial)

% locations struct ahs isExcluded vector that keeps track of whether a frame belongs to a trial that should be excluded
% given isExcluded vector (on index per frame in video), a vector of trial numbers, and excludeTrial vector (same length as trials),
% this function returns new isExcluded vector in which entires corresponding to trials are set to values in ecludeTrial

isExcludedNew = isExcluded;

for i = 1:length(trials)
    
    trialInds = (trialIdentities==trials(i));
    isExcludedNew(trialInds) = excludeTrial(i);
    
end