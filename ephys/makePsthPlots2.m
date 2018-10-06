%% select sessions

session = '181002_002';
cellNum = 2;
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx'], 'Sheet', 'sessions');

%% make PSTHs

eventName = 'obsOnToObsOff'; % rewardTimes, wiskContacts, obsOnTimes, obsLightOnTimes, obsOnToObsOff

    
if strcmp(eventName, 'wiskContacts')
    temp = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
        'wiskContactFrames', 'frameTimeStamps');
    sessionEvents = temp.frameTimeStampsWisk(temp.wiskContactFrames(temp.wiskContactFrames>0));
elseif strcmp(eventName, 'obsOnToObsOff')
    temp = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
        'obsOnTimes', 'obsOffTimes');
    sessionEvents = cat(2, temp.obsOnTimes, temp.obsOffTimes);
else
    sessionEvents = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), eventName);
    sessionEvents = sessionEvents.(eventName);
end

cellAxis = plotPSTH2(session, cellNum, {sessionEvents}, eventName);


%% test multi-condition PSTHs

temp = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'obsOnTimes', 'obsOffTimes');
sessionEvents = cell(1,2);
sessionEvents{1} = temp.obsOnTimes;
sessionEvents{2} = temp.obsOffTimes;


cellAxis = plotPSTH2(session, cellNum, sessionEvents, eventName);


%% plot inter-reward epochs

temp = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'rewardTimes');
sessionEvents = nan(length(temp.rewardTimes)-1, 2);
sessionEvents(:,1) = temp.rewardTimes(1:end-1);
sessionEvents(:,2) = temp.rewardTimes(2:end);

cellAxis = plotPSTH2(session, cellNum, sessionEvents);

%% plot time between obsOn and wisk contact

eventTimes = cell(1,length(sessions));

for i = 1:length(sessions)
    
    temp = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
            'wiskContactFrames', 'frameTimeStamps', 'obsOnTimes');
    sessionEvents = nan(0,2);
    
	for j = 1:length(temp.obsOnTimes)
        if temp.wiskContactFrames(j)>0
            wiskTime = temp.frameTimeStamps(temp.wiskContactFrames(j));
            if wiskTime>temp.obsOnTimes(j)
                sessionEvents(end+1,:) = [temp.obsOnTimes(j) wiskTime];
            end
        end
    end
    eventTimes{i} = sessionEvents;
end

plotPSTH(sessions, eventTimes, 'obsOnToWiskContact');



%% plot step tuning

paws = 1:4;
percentileLimits = [40 60]; % only include swings within these percentile limits

    
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'],...
        'frameTimeStamps', 'wheelPositions', 'wheelTimes', 'wheelCenter', 'wheelRadius', 'mToPixMapping');
locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' session '\trackedFeaturesRaw.csv']); % get raw tracking data
[locations, features] = fixTrackingDLC(locationsTable, frameTimeStamps);
topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, wheelTimes, wheelCenter, wheelRadius, 250, mToPixMapping(1));


eventTimes = cell(1,length(paws));
for paw = paws
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

    eventTimes{paw} = sessionEvents;
end

cellAxis = plotPSTH2(session, cellNum, eventTimes, 'stepTuning');


%% plot control vs mod step tuning

data = getKinematicData4({session}, sessionInfo, []);
%%
paw = 3;
buffer = [0 0];
% load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'],...
%             'obsPositions', 'obsTimes', 'obsPixPositions', 'obsPixPositionsContinuous', 'frameTimeStamps', 'mToPixMapping', 'isLightOn', ...
%             'obsOnTimes', 'obsOffTimes', 'nosePos', 'targetFs', 'wheelPositions', 'wheelTimes', 'targetFs', ...
%             'wheelRadius', 'wheelCenter', 'obsHeightsVid', 'touchesPerPaw', 'wiskContactFrames', 'frameTimeStampsWisk');
% load([getenv('OBSDATADIR') 'sessions\' session '\run.mat'], 'breaks');
% locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' session '\trackedFeaturesRaw.csv']); % get raw tracking data
% [locations, features] = fixTrackingDLC(locationsTable, frameTimeStamps);
% botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
% topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
% stanceBins = getStanceBins(frameTimeStamps, locations(:,:,topPawInds), wheelPositions, wheelTimes, wheelCenter, wheelRadius, 250, mToPixMapping(1));
% contactTimes = nan(size(obsOnTimes));
% [controlStepIdentities, modifiedStepIdentities, noObsStepIdentities] = ...
%         getStepIdentities(stanceBins, locations(:,:,botPawInds), contactTimes, frameTimeStamps, ...
%         obsOnTimes, obsOffTimes, obsPixPositions, obsPixPositionsContinuous, controlSteps, noObsSteps);



events = cell(1,2);
events{1} = nan(length(data),2);
events{2} = nan(length(data),2);
numModSteps = nan(1,length(data));

for i = 1:length(data)
    
%     paw = data(i).firstPawOver;
%     if paw==2; paw=3; else; paw=2; end % change to second paw over
    paw = data(i).firstModPaw;
    numModSteps(i) = data(i).modStepNum(paw);
    
    % get control start and stop time for paw
    bins = data(i).trialControlStepIdentities(:,paw) == max(data(i).trialControlStepIdentities(:,paw));
    controlTimes = data(i).frameTimeStamps(bins);
    bins = data(i).modifiedStepIdentities(:,paw) == min(data(i).modifiedStepIdentities(:,paw));
    modTimes = data(i).frameTimeStamps(bins);
    
    events{1}(i,:) = [controlTimes(1)-buffer(1) controlTimes(end)+buffer(2)];
    events{2}(i,:) = [modTimes(1)-buffer(1) modTimes(end)+buffer(2)];
end

% validBins = [data.firstPawOver]~=paw;
validBins = numModSteps==1;
for i = 1:length(events); events{i}=events{i}(validBins,:); end
cellAxis = plotPSTH2(session, cellNum, events, 'stepTuning');

