function [controlStepIdentities, modifiedStepIdentities, noObsStepIdentities, trialStartInds] = ...
    getStepIdentities(stanceBins, locationsBotPaws, contactTimes, frameTimeStamps, ...
    obsOnTimes, obsOffTimes, obsPixPositionsContinuous, controlSteps, noObsSteps)

% returns three identities for control, modified, and noObs steps //
% modified steps are those occuring during or after whisker contact up
% until the step over the obs // control and the 'controlSteps' steps
% before modified steps // noObs steps are the 'noObsSteps' steps occuring
% before the obstacle turns on // the output matrices are of size
% (numberOfFrames X 4), with all NaNs except for inds whether the paw is
% swinging for the given step type, in which case 1s indicated the first
% swing, 2s for the second swing, etc... // trialStartsInds records the
% index of the first step (noObsStep) for each trial


% settings
plotExample = false;


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
trialStartInds = knnsearch(frameTimeStamps, [0; obsOffTimes(1:end-1)]);  % default trial start ind is previous obsOffTime, but will be adjusted to the ind of the first noObs step

for i = 1:length(obsOnTimes)
    
    % extract trial data to make indexing operations faster
    if i==1; minTime=0; else; minTime=obsOffTimes(i-1); end
    if i==length(obsOnTimes); maxTime=frameTimeStamps(end); else; maxTime=obsOnTimes(i+1); end
    inds = find(frameTimeStamps>minTime & frameTimeStamps<maxTime);
    locationsTrial = locationsBotPaws(inds,:,:);
    obsPosTrial = obsPixPositionsContinuous(:,inds);
    swingIdsTrial = allSwingIdentities(inds,:);
    timesTrial = frameTimeStamps(inds);
    
    try
        for j = 1:4
            % find id for swing that crosses obs
            overObsInd = find(locationsTrial(:,1,j)<obsPosTrial(i,:)', 1, 'last') + 1; % finding last frame that is less than obsPos, rather than first frame that is greater than obsPos, prevents me from getting steps that go under the obs (because these have to return behind obsPos and then get over it again)
            swingOverObsIdentity = swingIdsTrial(overObsInd, j);  % step number 'swingOverObsIdentity' is the first to get over the obstacle
            
            % find id of first swing during or after obs contact with wisk
            firstModifiedInd = find(~isnan(swingIdsTrial(:,j)) & timesTrial>=contactTimes(i), 1, 'first');
            firstModifiedIdentity = swingIdsTrial(firstModifiedInd, j);
            
            % find id of first swing during or after obs turns on
            firstObsOnInd = find(~isnan(swingIdsTrial(:,j)) & timesTrial>=obsOnTimes(i), 1, 'first');
            firstObsOnIdentitiy = swingIdsTrial(firstObsOnInd, j);

            % find bins wihtin swingIdsTrial
            modifiedBins = (swingIdsTrial(:,j) >= firstModifiedIdentity) & ...
                           (swingIdsTrial(:,j) <= swingOverObsIdentity);
            controlBins = (swingIdsTrial(:,j) >= (firstModifiedIdentity-controlSteps)) & ...
                          (swingIdsTrial(:,j) < firstModifiedIdentity);
            obsOffBins = (swingIdsTrial(:,j) >= (firstObsOnIdentitiy-noObsSteps)) & ...
                         (swingIdsTrial(:,j) < firstObsOnIdentitiy);
            
            % store results
            modifiedStepIdentities(inds(modifiedBins), j) = swingIdsTrial(modifiedBins,j)-firstModifiedIdentity+1;
            controlStepIdentities(inds(controlBins), j) = swingIdsTrial(controlBins,j)-(firstModifiedIdentity-controlSteps)+1;
            noObsStepIdentities(inds(obsOffBins), j) = swingIdsTrial(obsOffBins,j)-(firstObsOnIdentitiy-noObsSteps)+1;
        end
        trialStartInds(i) = find(any(noObsStepIdentities,2) & frameTimeStamps>minTime, 1, 'first');
    catch
        fprintf('  problem with trial %i\n', i)
    end
end




% plot control and modified step segmentation
if plotExample
    
    % get trial
    trial = randperm(length(obsOnTimes)-1, 1)+1;
    trialBins = (1:length(frameTimeStamps)>trialStartInds(trial))' & frameTimeStamps<obsOffTimes(trial);  % plus one to be safe
    paws = [1 2 3 4];
    xLocations = squeeze(locationsBotPaws(:,1,:)) - obsPixPositionsContinuous(trial,:)';
    colors = hsv(4);

    figure('color', 'white', 'position', [1944.00 65.00 1335.00 838.00]);

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
    xLims = get(gca, 'xlim');
    line(xLims, [0 0])
    line([contactTimes(trial) contactTimes(trial)], get(gca,'ylim'))
    text(contactTimes(trial), min(get(gca, 'ylim')), 'WHISKER CONTACT', 'VerticalAlignment', 'bottom')
    line([obsOnTimes(trial) obsOnTimes(trial)], get(gca,'ylim'))
    text(obsOnTimes(trial), min(get(gca, 'ylim')), 'OBSTACLE ON', 'VerticalAlignment', 'bottom')

    set(gca, 'box', 'off', 'TickDir', 'out', 'xlim', xLims)
    ylabel('time (s)')
    ylabel('position relative to obstacle (m)')
    pause(.001)
end


