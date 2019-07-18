function plotManyPSTHs(session, opts)

% for a given session, plots a series of PSTHs for each cell in the session


% settings
s.folder = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'PSTHs');  % folder in which the PSTHs will be saved
s.rows = 4;
s.cols = 5;
s.stepPercentiles = [40 60]; % only include steps with durations in between these percentile limits
s.pawNames = {'LH', 'LF', 'RF', 'RH'};
s.pawColors = hsv(4);
s.mainColor = [.8 .4 1];

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end


% initializations
folder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
load(fullfile(folder, 'runAnalyzed.mat'), ...
        'obsOnTimes', 'obsOffTimes',  'wiskContactFrames', 'frameTimeStamps', ...
        'frameTimeStampsWisk', 'rewardTimes', 'isLightOn', 'touches', 'touchesPerPaw', 'touchClassNames');

% load kinData if it exists // otherwise compute kinData    
if exist(fullfile(folder, 'kinData.mat'))
    load(fullfile(folder, 'kinData.mat'), 'kinData', 'stanceBins')
else
    [kinData, stanceBins] = getKinematicData5(session);
end

load(fullfile(folder, 'neuralData.mat'), 'unit_ids');

for cellNum = 1:length(unit_ids)
    
    fprintf('%s: plotting cell %i/%i\n', session, cellNum, length(unit_ids))
    figure('name', sprintf('%s - unit %i', session, unit_ids(cellNum)), ...
        'color', 'white', 'MenuBar', 'none', 'units', 'pixels', 'position', [2000 20 1800 1000]); hold on
    plotInd = 0;
    
    % reward delivery
    plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);
    plotPSTH2(session, cellNum, rewardTimes);
    xlabel('reward delivery')
    
    % reward delivery -> reward delivery
    plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);
    times = nan(length(rewardTimes)-1, 2);
    times(:,1) = rewardTimes(1:end-1);
    times(:,2) = rewardTimes(2:end);
    plotPSTH2(session, cellNum, times);
    xlabel('reward delivery -> reward delivery')
    
    % obs on
    plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);
    times = {obsOnTimes(isLightOn), obsOnTimes(~isLightOn)};
    plotPSTH2(session, cellNum, times, ...
        {'conditionNames', {'light on', 'light off'}});
    xlabel('obstacle turns on')
    
    % obs on -> obs off
    plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);
    times = {cat(2, obsOnTimes(isLightOn), obsOffTimes(isLightOn)), ...
             cat(2, obsOnTimes(~isLightOn), obsOffTimes(~isLightOn))};
    plotPSTH2(session, cellNum, times, ...
        {'conditionNames', {'light on', 'light off'}});
    xlabel('obstacle on -> obstacle off')
    
    % wisk contact
    plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);
    times = {frameTimeStampsWisk(wiskContactFrames(wiskContactFrames>0 & isLightOn)), ...
             frameTimeStampsWisk(wiskContactFrames(wiskContactFrames>0 & ~isLightOn))};
    plotPSTH2(session, cellNum, times, ...
        {'conditionNames', {'light on', 'light off'}});
    xlabel('whisker contact')
    
    % obsOn -> wisk contact
    plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);
    allTimes = nan(0,2);
    allTimesIsLightOn = false(0,0);
    wiskContactTimes = frameTimeStampsWisk(wiskContactFrames(wiskContactFrames>0));
	for i = 1:length(obsOnTimes)
        contactInd = find(wiskContactTimes>obsOnTimes(i) & wiskContactTimes<obsOffTimes(i), 1, 'first');
        if ~isempty(contactInd)
            allTimes(end+1,:) = [obsOnTimes(i) wiskContactTimes(contactInd)];
            allTimesIsLightOn(end+1) = isLightOn(i);
        end
    end
    times = {allTimes(allTimesIsLightOn,:), allTimes(~allTimesIsLightOn,:)};
    plotPSTH2(session, cellNum, times, ...
        {'conditionNames', {'light on', 'light off'}});
    xlabel('obstacle on -> whisker contact')
    
    
    % step tuning for all paws
    plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);
    times = cell(1,4);
    for paw = 1:4
        
        % get swing start and end times
        swingStartInds = find(diff(~stanceBins(:,paw))==1);
        swingStartTimes = frameTimeStamps(swingStartInds(1:end-1));
        swingEndTimes = frameTimeStamps(swingStartInds(2:end)-1);
        sessionEvents = cat(2, swingStartTimes, swingEndTimes);
        sessionEvents = sessionEvents(~isnan(sum(sessionEvents,2)),:); % remove nan entries

        % only take steps in middle of duration distribution
        durations = diff(sessionEvents,1,2);
        durationLimits = prctile(durations, s.stepPercentiles);
        sessionEvents = sessionEvents(durations>durationLimits(1) & durations<durationLimits(2), :);

        times{paw} = sessionEvents;
    end
    plotPSTH2(session, cellNum, times, ...
        {'conditionNames', s.pawNames, 'errorFcn', false});
    xlabel('swing start -> stance end')
    
    
    % step over vs. control steps
    for paw = 1:4
        times = cell(1,2);
        [times{:}] = deal(nan(length(kinData), 2));
        
        % get step start and stop times
        for i = 1:length(kinData)
            stepOverBins = kinData(i).modifiedStepIdentities(:,paw) == max(kinData(i).modifiedStepIdentities(:,paw));
            ctlBins = kinData(i).controlStepIdentities(:,paw) == max(kinData(i).controlStepIdentities(:,paw));
            times{1}(i,:) = frameTimeStamps(kinData(i).trialInds([find(stepOverBins,1,'first'), find(stepOverBins,1,'last')]));
            times{2}(i,:) = frameTimeStamps(kinData(i).trialInds([find(ctlBins,1,'first'), find(ctlBins,1,'last')]));
        end

        % plot it!
        plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);
        plotPSTH2(session, cellNum, times, ...
            {'conditionNames', {'step over', 'control'}, 'colors', [s.pawColors(paw,:); s.pawColors(paw,:)*.2]});
        xlabel(sprintf('%s: swing start -> swing end', s.pawNames{paw}))
    end
    
    
    % one vs two step for first modified paw, left and rightforepaws
    for paw = 2:3
        
        % get all step times
        numModSteps = cellfun(@(x) max(x(:,paw)), {kinData.modifiedStepIdentities});
        stepTimes = nan(length(kinData), 2);
        for i = 1:length(kinData)
            firstModBins = kinData(i).modifiedStepIdentities(:,paw) == 1;
            stepTimes(i,:) = frameTimeStamps(kinData(i).trialInds([find(firstModBins,1,'first'), find(firstModBins,1,'last')]));
        end
        times = cell(1,2);
        times{1} = stepTimes([kinData.firstModPaw]==paw & numModSteps==1, :); % one step times
        times{2} = stepTimes([kinData.firstModPaw]==paw & numModSteps==2, :); % two step times
        
        % plot it!
        plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);
        plotPSTH2(session, cellNum, times, ...
            {'conditionNames', {'one step', 'two step'}, 'colors', [s.pawColors(paw,:); s.pawColors(paw,:)*.2]});
        xlabel(sprintf('%s: swing start -> swing end', s.pawNames{paw}))
    end
    
    
    % ventral paw contacts
    minTouchInds = 1; % paw must be touching for at least this many inds for trial to be included
    contactsToInclude = {'fore_ventral', 'hind_ventral_low'};
%     contactsToInclude = {'hind_ventral_low'};
    classInds = find(ismember(touchClassNames, contactsToInclude));
    
    for paw = 1:4
        plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);
        times = [];
        
        % get first contact for paw on each trial
        for i = 1:length(obsOnTimes)
            trialBins = frameTimeStamps>obsOnTimes(i) & frameTimeStamps<obsOffTimes(i);
            pawTouchBins = trialBins & touchesPerPaw(:,paw) & ismember(touches, classInds)';
            pawTouchInd = find(pawTouchBins, 1, 'first');
            if sum(pawTouchBins)>=minTouchInds; times(end+1) = frameTimeStamps(pawTouchInd); end
        end
        
        plotPSTH2(session, cellNum, times, {'colors', s.pawColors(paw,:)});
        xlabel(sprintf('%s: ventral touch', s.pawNames{paw}))
    end
    
    
    % save
    pause(.01)
    fileName = fullfile(folder, [session 'unit' num2str(unit_ids(cellNum))]);
%     savefig(fileName)
    saveas(gcf, [fileName '.png'])
end
disp('all done!')

