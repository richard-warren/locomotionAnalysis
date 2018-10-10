
% make speed montage videos


% get velocity data for all sessions
experiments = {'obsBr', 'obsWiskMrk', 'obsCompareOn', 'obsCompareOnOff'};
obsPos = 0.3820; % obsPos at which obs is at center of wheel
posRange = .06; % computer vel between obsPos plus or minus posRange
allTrialsDataRaw = getAllSessionTrialSpeeds(experiments, obsPos, posRange);
allTrialsData = allTrialsDataRaw(~[allTrialsDataRaw.lightOn]);

%% make videos in dft speed bins with no pauses


% settings
trialNum = 10;
binLims = [.15 .65];
binSeparation = .1;

bins = (binLims(1) : binSeparation : binLims(2)-binSeparation)';
bins = cat(2, bins, bins+binSeparation);

for i = 1:size(bins,1)
    makeSpeedBinnedVid(trialNum, bins(i,:), allTrialsData, false);
end

%% make version with pauses
trialNum = 50;
speedLims = [.4 .7];

makeSpeedBinnedVid(trialNum, speedLims, allTrialsData, true);
