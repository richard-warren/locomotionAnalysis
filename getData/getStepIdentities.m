function [controlStepIdentities, modifiedStepIdentities, noObsStepIdentities] = ...
    getStepIdentities(stanceBins, locationsBotPaws, contactTimes, frameTimeStamps, ...
    obsOnTimes, obsOffTimes, obsPixPositions, obsPixPositionsContinuous, controlSteps, noObsSteps)

% returns three identities for control, modified, and noObs steps //
% modified steps are those occuring during or after whisker contact up
% until the step over the obs // control and the 'controlSteps' steps
% before modified steps // noObs steps are the 'noObsSteps' steps occuring
% before the obstacle turns on // the output matrices are of size
% (numberOfFrames X 4), with all NaNs except for inds whether the paw is
% swinging for the given step type, in which case 1s indicated the first
% swing, 2s for the second swing, etc...


% settings
plotExample = true;


% give each swing a number across all trials, in ascending order
swingBins = ~stanceBins;
swingDeltas = double([zeros(1,4); diff(swingBins,1,1)]);  % matrix where 1 indcates swing start and -1 indicates stance start
swingDeltas = swingDeltas .* [zeros(1,4); cumsum(diff(swingBins,1,1)==1, 1)];  % matrix with n at swing start and -n at start of next stance, where n indicates the swing number
allSwingIdentities = cumsum(swingDeltas);  % matrix where inds corresponding to swing n have value n, and zeros elsewhere
allSwingIdentities(allSwingIdentities==0) = nan;


% get step identities
controlStepIdentities = nan(size(allSwingIdentities));
modifiedStepIdentities = nan(size(allSwingIdentities));
noObsStepIdentities = nan(size(allSwingIdentities));
keyboard

for i = 1:length(obsOnTimes)
    for j = 1:4
        try
            % .01 find id for swing that crosses obs
            overObsInd = find(locationsBotPaws(:,1,j)<obsPixPositionsContinuous(i,:)', 1, 'last') + 1; % finding last frame that is less than obsPos, rather than first frame that is greater than obsPos, prevents me from getting steps that go under the obs (because these have to return behind obsPos and then get over it again)
            swingOverObsIdentity = allSwingIdentities(overObsInd, j);  % step number 'swingOverObsIdentity' is the first to get over the obstacle
            
            % .06 find id of first swing during or after obs contact with wisk
            firstModifiedInd = find(~isnan(allSwingIdentities(:,j)) & frameTimeStamps>=contactTimes(i), 1, 'first');
            firstModifiedIdentity = allSwingIdentities(firstModifiedInd, j);
            
            % .06 find id of first swing during or after obs turns on
            firstObsOnInd = find(~isnan(allSwingIdentities(:,j)) & frameTimeStamps>=obsOnTimes(i), 1, 'first');
            firstObsOnIdentitiy = allSwingIdentities(firstObsOnInd, j);

        
            % .018
            modifiedBins = (allSwingIdentities(:,j) >= firstModifiedIdentity) & ...
                           (allSwingIdentities(:,j) <= swingOverObsIdentity);
            controlBins = (allSwingIdentities(:,j) >= (firstModifiedIdentity-controlSteps)) & ...
                          (allSwingIdentities(:,j) < firstModifiedIdentity);
            obsOffBins = (allSwingIdentities(:,j) >= (firstObsOnIdentitiy-noObsSteps)) & ...
                         (allSwingIdentities(:,j) < firstObsOnIdentitiy);
            
            % .002
            modifiedStepIdentities(modifiedBins, j) = allSwingIdentities(modifiedBins,j)-firstModifiedIdentity+1;
            controlStepIdentities(controlBins, j) = allSwingIdentities(controlBins,j)-(firstModifiedIdentity-controlSteps)+1;
            noObsStepIdentities(obsOffBins, j) = allSwingIdentities(obsOffBins,j)-(firstObsOnIdentitiy-noObsSteps)+1;
        catch
            fprintf('  problem with trial %i\n', i)
        end
    end
end


% plot control and modified step segmentation
if plotExample
    
    % get trial
    trial = randperm(length(obsOnTimes)-1, 1)+1;
    trialBins = frameTimeStamps>obsOffTimes(trial-1) & frameTimeStamps<obsOffTimes(trial);
    paws = [1 2 3 4];
    xLocations = squeeze(locationsBotPaws(:,1,:)) - obsPixPositionsContinuous(trial,:)';
    colors = hsv(4);

    figure;

    % plot x positions
    for i = 1:length(paws)

        % plot all x positions
        plot(frameTimeStamps(trialBins), xLocations(trialBins, paws(i)), ...
            'linewidth', 2, 'color', colors(paws(i),:)); hold on

        % highlight control swings
        for j = 1:max(controlStepIdentities(trialBins,paws(i)))
            controlBins = trialBins & controlStepIdentities(:,paws(i))==j;
            plot(frameTimeStamps(controlBins), xLocations(controlBins, paws(i)), ...
                'linewidth', 5, 'color', [0 0 0]); hold on
        end

        % highlight modified swings
        for j = 1:max(modifiedStepIdentities(trialBins,paws(i)))
            modBins = trialBins & modifiedStepIdentities(:,paws(i))==j;
            plot(frameTimeStamps(modBins), xLocations(modBins, paws(i)), ...
                'linewidth', 5, 'color', colors(paws(i),:)); hold on
        end
        
        % highlight no obs swings
        for j = 1:max(noObsStepIdentities(trialBins,paws(i)))
            noObsBins = trialBins & noObsStepIdentities(:,paws(i))==j;
            plot(frameTimeStamps(noObsBins), xLocations(noObsBins, paws(i)), ...
                'linewidth', 5, 'color', [.5 .5 .5]); hold on
        end
    end
    
    % add lines for obs position and obsContactTime
    line(get(gca,'xlim'), [0 0])
    line([contactTimes(trial) contactTimes(trial)], get(gca,'ylim'))
    line([obsOnTimes(trial) obsOnTimes(trial)], get(gca,'ylim'))

    pimpFig
end








