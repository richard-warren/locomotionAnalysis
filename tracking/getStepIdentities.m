function [controlStepIdentities, modifiedStepIdentities] = ...
    getStepIdentities(stanceBins, locationsBotPaws, contactTimes, frameTimeStamps, obsOnTimes, obsOffTimes, obsPixPositions)

% given the bins in which paws are in stance and time of wisk contact with obs, determines frames corresponding to control and modified swings
% modified swings are those occuring during or after obs contact with wisk, including but not after swing that actually gets over obs
% control swings are those before mod swings, for now two swings before mod swings for each paw


% settings
controlSteps = 2;
plotExample = false;

% give each swing a number across all trials, in ascending order
swingBins = ~stanceBins;
swingDeltas = [zeros(1,4); double(diff(swingBins,1,1)==1) - double(diff(swingBins,1,1)==-1)]; % matrix where 1 indcates swing start and -1 indicates stannce stance (not swing end, but stance stance, which is the next ind over)
swingDeltas = swingDeltas .* [zeros(1,4); cumsum(diff(swingBins,1,1)==1, 1)]; % matrix with n at swing start and -n at start of next stance, where n indicates the swing number
allSwingIdentities = cumsum(swingDeltas); % matrix where inds corresponding to swing n have value n, and zeros elsewhere
allSwingIdentities(allSwingIdentities==0) = nan;



% get trial swing identities and define control and modified steps

controlStepIdentities = nan(size(allSwingIdentities));
modifiedStepIdentities = nan(size(allSwingIdentities));


for i = 1:length(obsOnTimes)
    for j = 1:4

        % find id of swing that crosses obs
        overObsInd = find(frameTimeStamps>obsOnTimes(i) & ... 
                          frameTimeStamps<obsOffTimes(i) & ...
                          locationsBotPaws(:,1,j)>=obsPixPositions', 1, 'first');
        swingOverObsIdentity = allSwingIdentities(overObsInd, j);

        % find id of first swing during or after obs contact with wisk
        firstModifiedInd = find(~isnan(allSwingIdentities(:,j)) & frameTimeStamps>=contactTimes(i), 1, 'first');
        firstModifiedIdentitiy = allSwingIdentities(firstModifiedInd, j);

        try
            modifiedBins = (allSwingIdentities(:,j) >= firstModifiedIdentitiy) & ...
                           (allSwingIdentities(:,j) <= swingOverObsIdentity);
            controlBins = (allSwingIdentities(:,j) >= (firstModifiedIdentitiy-controlSteps)) & ...
                          (allSwingIdentities(:,j) < firstModifiedIdentitiy);

            steps = cumsum([0; diff(modifiedBins)==1]);
            modifiedStepIdentities(modifiedBins, j) = steps(modifiedBins);
            steps = cumsum([0; diff(controlBins)==1]);
            controlStepIdentities(controlBins, j) = steps(controlBins);
        catch
            fprintf('  problem with trial %i\n', i)
        end
        
    end
end


% plot control and modified step segmentation
if plotExample

    % get trial
    trial  = randperm(length(obsOnTimes), 1);
    trialBins = frameTimeStamps>obsOnTimes(trial) & frameTimeStamps<obsOffTimes(trial);
    paws = [1 2 3 4];
    xLocations = squeeze(locationsBotPaws(:,1,:)) - obsPixPositions';
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

    end
    % add lines for obs position and obsContactTime
    line(get(gca,'xlim'), [0 0])
    line([contactTimes(trial) contactTimes(trial)], get(gca,'ylim'))

    pimpFig
end










