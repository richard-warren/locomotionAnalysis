% function getTrialIndsFromStepIdentities(controlStepIdentities, modifiedStepIdentities)

% when correcting tracking on the top view, it will be useful to restrict correction to frames including no more than control and modified steps, as determined by getStepIdentities
% this function takes stepIdentities and returns inds of frames containing for each trial all frames between first control start and last mod step


