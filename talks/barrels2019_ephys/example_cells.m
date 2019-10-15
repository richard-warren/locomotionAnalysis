%% global settings

clear all
tcm190910_config;
figSettings = {};




%% whisker contact

close all; figure('position', [766.00 169.00 543.00 638.00], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false);

sessions = {'181020_001', '181020_001', '181103_000'};
cellNums = [37, 65, 140];


for i = 1:length(sessions)
    
    subplot(length(sessions), 1, i)
    
    sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i});
    load(fullfile(sessionFolder, 'neuralData.mat'), 'unit_ids');
    load(fullfile(sessionFolder, 'runAnalyzed.mat'), ...
            'obsOnTimes', 'obsOffTimes',  'wiskContactFrames', 'frameTimeStamps', 'wiskContactTimes', ...
            'frameTimeStampsWisk', 'rewardTimes', 'isLightOn', 'touches', 'touchesPerPaw', 'touchClassNames');

    % wisk contact
    plotPSTH2(sessions{i}, find(unit_ids==cellNums(i)), {wiskContactTimes}, ...
        {'xLims', [-.2 .4], 'yLimsNormalized', [-2 6], 'plotLegend', false, 'colors', [245 245 66]/255});
    
    if i==length(sessions)
        xlabel('time from whisker contact (s)')
        ylabel('spike rate')
    end
    
    yticks = get(gca, 'YTick');
    set(gca, 'XColor', axisColor, 'YColor', axisColor, 'Color', 'black', 'FontSize', fontSize, 'FontName', font, 'YTick', [yticks(1) yticks(end)])
    print -clipboard -dmeta
end


%% whisker ramps

close all; figure('position', [2500.00 253.00 543.00 538.00], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false);
sessions = {'181019_002', '181103_000'};
cellNums = [67, 140];


for i = 1:length(sessions)
    
    subplot(length(sessions), 1, i)
    
    sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i});
    load(fullfile(sessionFolder, 'neuralData.mat'), 'unit_ids');
    load(fullfile(sessionFolder, 'runAnalyzed.mat'), ...
            'obsOnTimes', 'obsOffTimes',  'wiskContactFrames', 'frameTimeStamps', 'wiskContactTimes', ...
            'frameTimeStampsWisk', 'rewardTimes', 'isLightOn', 'touches', 'touchesPerPaw', 'touchClassNames');

    % wisk contact
    plotPSTH2(sessions{i}, find(unit_ids==cellNums(i)), {wiskContactTimes}, ...
        {'xLims', [-1 .4], 'yLimsNormalized', [-2 4.2], 'plotLegend', false, 'colors', [245 245 66]/255});
    
    if i==length(sessions)
        xlabel('time from whisker contact (s)')
        ylabel('spike rate')
    end
    
    yticks = get(gca, 'YTick');
    set(gca, 'XColor', axisColor, 'YColor', axisColor, 'Color', 'black', 'FontSize', fontSize, 'FontName', font, 'YTick', [yticks(1) yticks(end)])
    print -clipboard -dmeta
end

%% step tuning


close all; figure('position', [711.00 294.00 633.00 438.00], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false);
sessions = {'180917_002', '180917_002', '180920_002', '180922_001', '181001_002', '181001_002', '181002_002', '181002_002', '181030_000'};
cellNums = [64, 82, 70, 73, 6, 93, 2, 23, 132];
stepPercentiles = [40 60];


for i = 1:length(sessions)
    
    subplot(3, 3, i)
    
    sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i});
    load(fullfile(sessionFolder, 'neuralData.mat'), 'unit_ids');
    load(fullfile(sessionFolder, 'runAnalyzed.mat'), ...
            'obsOnTimes', 'obsOffTimes',  'wiskContactFrames', 'frameTimeStamps', 'wiskContactTimes', ...
            'frameTimeStampsWisk', 'rewardTimes', 'isLightOn', 'touches', 'touchesPerPaw', 'touchClassNames');
    load(fullfile(sessionFolder, 'kinData.mat'), 'kinData', 'stanceBins')    
    

    % step tuning for all paws
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
    
    plotPSTH2(sessions{i}, find(unit_ids==cellNums(i)), times, ...
        {'xLims', [-1 .4], 'yLimsNormalized', [-1.5 1.5], 'errorFcn', false, 'plotLegend', false, 'colors', pawColors});
    
    if i==7
        xlabel('fraction of step cycle')
        ylabel('spike rate')
    end
    
    yticks = get(gca, 'YTick');
    set(gca, 'XColor', axisColor, 'YColor', axisColor, 'Color', 'black', ...
        'FontSize', fontSize, 'FontName', font, 'XTick', [], 'YTick', [yticks(1) yticks(end)])
    print -clipboard -dmeta
end





