function frameInds = getTrialIndsFromStepIdentities(controlStepIdentities, modifiedStepIdentities, ...
    frameTimeStamps, obsOnTimes, obsOffTimes, paws, buffer)

% when correcting tracking on the top view, it will be useful to restrict correction to frames including no more than control and modified steps, as determined by getStepIdentities
% this function takes stepIdentities and returns inds of frames containing for each trial all frames between first control start and last mod step
% buffer frames are added to the beginning and the end for every trial


frameBins = false(size(frameTimeStamps));
isSwing = ~isnan(controlStepIdentities) | ~isnan(modifiedStepIdentities);
isSwing = any(isSwing(:, paws),2); % true where there is a control or modified step in any of the paws
swingStarts = [0; diff(isSwing)==1];
swingEnds = [diff(isSwing)==-1; 0];

for i = 1:length(obsOnTimes)
    
    trialBins = frameTimeStamps>obsOnTimes(i) & frameTimeStamps<obsOffTimes(i);
    firstSwingInd = find(trialBins & swingStarts, 1, 'first') - buffer;
    lastSwingInd = find(trialBins & swingEnds, 1, 'last') + buffer;
    frameBins(firstSwingInd:lastSwingInd) = 1;
        
end

frameInds = find(frameBins);