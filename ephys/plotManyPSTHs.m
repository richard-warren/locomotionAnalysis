function plotManyPSTHs(session, cells)


% settings
rows = 4;
cols = 5;
stepPercentiles = [40 60]; % only include steps with durations in between these percentile limits
pawNames = {'LH', 'LF', 'RF', 'RH'};
pawColors = hsv(4);
mainColor = [.8 .4 1];



% initializations
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
        'obsOnTimes', 'obsOffTimes',  'wiskContactFrames', 'frameTimeStamps', ...
        'frameTimeStampsWisk', 'rewardTimes', 'isLightOn', 'touches', 'touchesPerPaw', 'touchClassNames');
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'Sheet', 'sessions');
[data, stanceBins] = getKinematicData4({session}, sessionInfo, []);
stanceBins = stanceBins{1}; % stanceBins will always have only a single entry when getKinematicData is called with a single session
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), 'unit_ids');
if ~exist('cells', 'var'); cells = 1:length(unit_ids); end

for cellNum = cells
    
    fprintf('%s: plotting cell %i/%i\n', session, cellNum, length(cells))
    figure('name', sprintf('%s - unit %i', session, unit_ids(cellNum)), ...
        'color', 'white', 'MenuBar', 'none', 'units', 'pixels', 'position', [2000 20 1800 1000]); hold on
    plotInd = 0;
    
    % reward delivery
    plotInd = plotInd + 1; subplot(rows, cols, plotInd);
    plotPSTH2(session, cellNum, rewardTimes, {}, mainColor);
    xlabel('reward delivery')
    
    % reward delivery -> reward delivery
    plotInd = plotInd + 1; subplot(rows, cols, plotInd);
    times = nan(length(rewardTimes)-1, 2);
    times(:,1) = rewardTimes(1:end-1);
    times(:,2) = rewardTimes(2:end);
    plotPSTH2(session, cellNum, times, {}, mainColor);
    xlabel('reward delivery -> reward delivery')
    
    % obs on
    plotInd = plotInd + 1; subplot(rows, cols, plotInd);
    times = {obsOnTimes(isLightOn), obsOnTimes(~isLightOn)};
    plotPSTH2(session, cellNum, times, {'light on', 'light off'}, [mainColor; mainColor*.25]);
    xlabel('obstacle turns on')
    
    % obs on -> obs off
    plotInd = plotInd + 1; subplot(rows, cols, plotInd);
    times = {cat(2, obsOnTimes(isLightOn), obsOffTimes(isLightOn)), ...
             cat(2, obsOnTimes(~isLightOn), obsOffTimes(~isLightOn))};
    plotPSTH2(session, cellNum, times, {'light on', 'light off'}, [mainColor; mainColor*.25]);
    xlabel('obstacle on -> obstacle off')
    
    % wisk contact
    plotInd = plotInd + 1; subplot(rows, cols, plotInd);
    times = {frameTimeStampsWisk(wiskContactFrames(wiskContactFrames>0 & isLightOn)), ...
             frameTimeStampsWisk(wiskContactFrames(wiskContactFrames>0 & ~isLightOn))};
    plotPSTH2(session, cellNum, times, {'light on', 'light off'}, [mainColor; mainColor*.25]);
    xlabel('whisker contact')
    
    % obsOn -> wisk contact
    plotInd = plotInd + 1; subplot(rows, cols, plotInd);
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
    plotPSTH2(session, cellNum, times, {'light on', 'light off'}, [mainColor; mainColor*.25]);
    xlabel('obstacle on -> whisker contact')
    
    
    % step tuning for all paws
    plotInd = plotInd + 1; subplot(rows, cols, plotInd);
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
        durationLimits = prctile(durations, stepPercentiles);
        sessionEvents = sessionEvents(durations>durationLimits(1) & durations<durationLimits(2), :);

        times{paw} = sessionEvents;
    end
    plotPSTH2(session, cellNum, times, pawNames, pawColors, false);
    xlabel('swing start -> stance end')
    
    
    % step over vs. control steps
    for paw = 1:4
        times = cell(1,2);
        
        % get mod step times
        stepTimes = cellfun(@(x,y) x(y(:,paw)==max(y(:,paw))), ...
            {data.frameTimeStamps}, {data.modifiedStepIdentities}, ...
            'UniformOutput', false);
        stepTimes = cellfun(@(x) [x(1) x(end)], stepTimes, 'UniformOutput', false);
        times{1} = cat(1, stepTimes{:});
        
        % get control step times
        stepTimes = cellfun(@(x,y) x(y(:,paw)==max(y(:,paw))), ...
            {data.frameTimeStamps}, {data.trialControlStepIdentities}, ...
            'UniformOutput', false);
        stepTimes = cellfun(@(x) [x(1) x(end)], stepTimes, 'UniformOutput', false);
        times{2} = cat(1, stepTimes{:});
        
        % plot it!
        plotInd = plotInd + 1; subplot(rows, cols, plotInd);
        plotPSTH2(session, cellNum, times, {'step over', 'control'}, [pawColors(paw,:); .6 .6 .6]);
        xlabel(sprintf('%s: swing start -> swing end', pawNames{paw}))
    end
    
    
    % one vs two step for first modified paw, left and rightforepaws
    for paw = 2:3
        
        % get all step times
        numModSteps = cellfun(@(x) x(1,paw), {data.modStepNum});
        stepTimes = cellfun(@(x,y) x(y(:,paw)==1), ...
            {data.frameTimeStamps}, {data.modifiedStepIdentities}, ...
            'UniformOutput', false);
        stepTimes = cellfun(@(x) [x(1) x(end)], stepTimes, 'UniformOutput', false);
        stepTimes = cat(1, stepTimes{:});
        
        times = cell(1,2);
        times{1} = stepTimes([data.firstModPaw]==paw & numModSteps==1, :); % one step times
        times{2} = stepTimes([data.firstModPaw]==paw & numModSteps==2, :); % two step times
        
        % plot it!
        plotInd = plotInd + 1; subplot(rows, cols, plotInd);
        plotPSTH2(session, cellNum, times, {'one step', 'two step'}, [pawColors(paw,:); .6 .6 .6]);
        xlabel(sprintf('%s: swing start -> swing end', pawNames{paw}))
    end
    
    % first vs second paw over
    times = cell(1,2);
    firstAndSecondPawsOver = {[data.firstPawOver], ~([data.firstPawOver]-2)+2};
    for i = 1:2
        stepTimes = cellfun(@(times,ids,paw) times(ids(:,paw)==max(ids(:,paw))), ...
                {data.frameTimeStamps}, {data.modifiedStepIdentities}, num2cell(firstAndSecondPawsOver{i}), ...
                'UniformOutput', false);
        stepTimes = cellfun(@(x) [x(1) x(end)], stepTimes, 'UniformOutput', false);
        stepTimes = cat(1, stepTimes{:});
        times{i} = stepTimes;
    end
    
    plotInd = plotInd + 1; subplot(rows, cols, plotInd);
    plotPSTH2(session, cellNum, times, {'first paw over', 'second paw over'}, [mainColor; .6 .6 .6]);
    xlabel(sprintf('swing start -> swing end'))
    
    % ventral paw contacts
    minTouchInds = 1; % paw must be touching for at least this many inds for trial to be included
    contactsToInclude = {'fore_ventral', 'hind_ventral_low'};
%     contactsToInclude = {'hind_ventral_low'};
    classInds = find(ismember(touchClassNames, contactsToInclude));
    
    for paw = 1:4
        plotInd = plotInd + 1; subplot(rows, cols, plotInd);
        times = [];
        
        % get first contact for paw on each trial
        for i = 1:length(obsOnTimes)
            trialBins = frameTimeStamps>obsOnTimes(i) & frameTimeStamps<obsOffTimes(i);
            pawTouchBins = trialBins & touchesPerPaw(:,paw) & ismember(touches, classInds)';
            pawTouchInd = find(pawTouchBins, 1, 'first');
            if sum(pawTouchBins)>=minTouchInds; times(end+1) = frameTimeStamps(pawTouchInd); end
        end
        
        plotPSTH2(session, cellNum, times, {}, pawColors(paw,:));
        xlabel(sprintf('%s: ventral touch', pawNames{paw}))
    end
    
    
    % save
    pause(.01)
    fileName = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'PSTHs', [session 'unit' num2str(unit_ids(cellNum))]);
%     savefig(fileName)
    saveas(gcf, [fileName '.png'])
end
disp('all done!')

