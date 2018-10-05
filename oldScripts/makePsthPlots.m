%% select sessions

sessions = {'180917_002', '180920_002', '180922_001', '181001_002'};
sessions = {'180917_002'};


%% make PSTHs

eventName = 'obsOnToObsOff'; % rewardTimes, wiskContacts, obsOnTimes, obsLightOnTimes, obsOnToObsOff

eventTimes = cell(1,length(sessions));
for i = 1:length(sessions)
    
    if strcmp(eventName, 'wiskContacts')
        temp = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
            'wiskContactFrames', 'frameTimeStamps');
        sessionEvents = temp.frameTimeStamps(temp.wiskContactFrames(temp.wiskContactFrames>0));
    elseif strcmp(eventName, 'obsOnToObsOff')
        temp = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
            'obsOnTimes', 'obsOffTimes');
        sessionEvents = cat(2, temp.obsOnTimes, temp.obsOffTimes);
    else
        sessionEvents = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), eventName);
        sessionEvents = sessionEvents.(eventName);
    end
    
    eventTimes{i} = sessionEvents;
end

plotPSTH(sessions, eventTimes, eventName);


%% plot inter-reward epochs

eventTimes = cell(1,length(sessions));

for i = 1:length(sessions)        
    temp = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), 'rewardTimes');
    sessionEvents = nan(length(temp.rewardTimes)-1, 2);
    sessionEvents(:,1) = temp.rewardTimes(1:end-1);
    sessionEvents(:,2) = temp.rewardTimes(2:end);
    eventTimes{i} = sessionEvents;
end

plotPSTH(sessions, eventTimes, 'rewardEpochs');

%% plot time between obsOn and wisk contact

eventTimes = cell(1,length(sessions));

for i = 1:length(sessions)
    
    temp = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
            'wiskContactFrames', 'frameTimeStamps', 'obsOnTimes', 'isLightOn');
    sessionEvents = nan(0,2);
    
	for j = 1:length(temp.obsOnTimes)
        if temp.wiskContactFrames(j)>0
            wiskTime = temp.frameTimeStamps(temp.wiskContactFrames(j));
            if wiskTime>temp.obsOnTimes(j) && ~isLightOn(j)
                sessionEvents(end+1,:) = [temp.obsOnTimes(j) wiskTime];
            end
        end
    end
    eventTimes{i} = sessionEvents;
end

plotPSTH(sessions, eventTimes, 'obsOnToWiskContact');



%% plot step tuning

sessions = {'180917_002', '180920_002', '180922_001', '181001_002'};
sessions = {'181002_002'};
eventTimes = cell(1,length(sessions));
paw = 4;
percentileLimits = [40 60]; % only include swings within these percentile limits

for i = 1:length(sessions)
    
    fprintf('%s: loading data...\n', sessions{i})
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'],...
            'frameTimeStamps', 'wheelPositions', 'wheelTimes', 'wheelCenter', 'wheelRadius', 'mToPixMapping');
    locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' sessions{i} '\trackedFeaturesRaw.csv']); % get raw tracking data
    [locations, features] = fixTrackingDLC(locationsTable, frameTimeStamps);
    topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
    stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, wheelTimes, wheelCenter, wheelRadius, 250, mToPixMapping(1));
    
    % get swing start and end times
    swingStartInds = find(diff(~stanceBins(:,paw))==1);
    swingStartTimes = frameTimeStamps(swingStartInds(1:end-1));
    swingEndTimes = frameTimeStamps(swingStartInds(2:end)-1);
    sessionEvents = cat(2, swingStartTimes, swingEndTimes);
    sessionEvents = sessionEvents(~isnan(sum(sessionEvents,2)),:); % remove nan entries
    
    % remove steps that take to long
    durations = diff(sessionEvents,1,2);
    durationLimits = prctile(durations, percentileLimits);
    sessionEvents = sessionEvents(durations>durationLimits(1) & durations<durationLimits(2), :);
    
    eventTimes{i} = sessionEvents;
end


plotPSTH(sessions, eventTimes, sprintf('paw%iSwingStartToStanceEnd', paw));








