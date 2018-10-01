




%% make PSTHs

eventName = 'rewardTimes'; % rewardTimes, wiskContacts, obsOnTimes, obsLightOnTimes

sessions = {'180917_002', '180920_002', '180922_001'};
% sessions = {'180922_001'};
eventTimes = cell(1,length(sessions));
for i = 1:length(sessions)
    
    if strcmp(eventName, 'wiskContacts')
        load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
            'wiskContactFrames', 'frameTimeStamps');
        sessionEvents = frameTimeStamps(wiskContactFrames(wiskContactFrames>0));
        clear wiskContactFrames frameTimeStamps
    else
        sessionEvents = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), eventName);
        sessionEvents = sessionEvents.(eventName);
    end
    
    eventTimes{i} = sessionEvents;
end

plotPSTH(sessions, eventTimes, eventName);


%% plot inter-reward epoches

sessions = {'180917_002', '180920_002', '180922_001'};
eventTimes = cell(1,length(sessions));

for i = 1:length(sessions)    
    
    temp = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), 'rewardTimes');
    sessionEvents = nan(length(temp.rewardTimes)-1, 2);
    sessionEvents(:,1) = temp.rewardTimes(1:end-1);
    sessionEvents(:,2) = temp.rewardTimes(2:end);
    eventTimes{i} = sessionEvents;
end

plotPSTH(sessions, eventTimes, 'rewardEpochs');

