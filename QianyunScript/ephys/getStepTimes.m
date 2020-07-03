function times = getStepTimes(session, paw, stepType, kinData)

% load runAnalyzed.mat
sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
load(fullfile(sessionFolder, 'runAnalyzed.mat'), 'frameTimeStamps');


% load kinData if it exists // otherwise compute kinData
if ~exist('kinData', 'var')
    if exist(fullfile(sessionFolder, 'kinData.mat'))
        load(fullfile(sessionFolder, 'kinData.mat'), 'kinData', 'stanceBins')
    else
        [kinData, stanceBins] = getKinematicData(session);
    end
end
% get the start & end times for modified steps
% times = cell(1,1);
% [times{:}] = deal(nan(length(kinData), 1));

switch stepType
    case 'stepOver'
        times = cell(length(kinData), 2);
        for i = find([kinData.isTrialAnalyzed])
            stepOverBins = kinData(i).modifiedStepIdentities(:,paw) == max(kinData(i).modifiedStepIdentities(:,paw));
            times(i, 1) = {frameTimeStamps(kinData(i).trialInds([find(stepOverBins,1,'first')]))};
            times(i, 2) = {frameTimeStamps(kinData(i).trialInds([find(stepOverBins,1,'last')]))};
        end
        
    case 'control'
        times = cell(length(kinData), 2);
        startTimes = [];
        endTimes = [];
        for i = find([kinData.isTrialAnalyzed])
            maxStepNum = max(kinData(i).controlStepIdentities(:,paw));
            for j = 1:maxStepNum
                stepBins = kinData(i).controlStepIdentities(:, paw) == j;
                startTimes(j) = frameTimeStamps(kinData(i).trialInds([find(stepBins,1,'first')]));
                endTimes(j) = frameTimeStamps(kinData(i).trialInds([find(stepBins,1,'last')]));
            end
            times(i, 1) = {startTimes};
            times(i, 2) = {endTimes};
        end                
        
    case 'noObs'
        times = cell(length(kinData), 2);
        startTimes = [];
        endTimes = [];
        for i = find([kinData.isTrialAnalyzed])
            maxStepNum = max(kinData(i).noObsStepIdentities(:,paw));
            for j = 1:maxStepNum
                stepBins = kinData(i).noObsStepIdentities(:, paw) == j;
                startTimes(j) = frameTimeStamps(kinData(i).trialInds([find(stepBins,1,'first')]));
                endTimes(j) = frameTimeStamps(kinData(i).trialInds([find(stepBins,1,'last')]));
            end
            times(i, 1) = {startTimes};
            times(i, 2) = {endTimes};
        end
end

times = cell2struct(times, {'startTimes', 'endTimes'}, 2);

end