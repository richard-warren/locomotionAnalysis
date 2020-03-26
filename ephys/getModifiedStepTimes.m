function times = getModifiedStepTimes(session, paw, opts)

% load runAnalyzed.mat
sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
load(fullfile(sessionFolder, 'runAnalyzed.mat'), 'frameTimeStamps');


% load kinData if it exists // otherwise compute kinData
if exist(fullfile(sessionFolder, 'kinData.mat'))
    load(fullfile(sessionFolder, 'kinData.mat'), 'kinData', 'stanceBins')
else
    [kinData, stanceBins] = getKinematicData(session);
end

% get the start & end times for modified steps
% times = cell(1,1);
% [times{:}] = deal(nan(length(kinData), 1));

times = nan(length(kinData), 2);

for i = find([kinData.isTrialAnalyzed])
    stepOverBins = kinData(i).modifiedStepIdentities(:,paw) == max(kinData(i).modifiedStepIdentities(:,paw));
    times(i, 1) = frameTimeStamps(kinData(i).trialInds([find(stepOverBins,1,'first')]));
    times(i, 2) = frameTimeStamps(kinData(i).trialInds([find(stepOverBins,1,'last')]));
    
end

end