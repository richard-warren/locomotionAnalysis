function obsAvoidanceLearningSummary(mice)

% shows obstacle avoidance for all mice over time, with and without wheel break...
% assumes at least noBrSessions have been collected... otherwise will behave incorrectly
%
% input         mice:      name of mice to analyze


% user settings
minTouchTime = .05; % only touches count that are >= minTouchTime
lightVsNoLightSessions = 3; % for the plot comparing light vs no light speed, only take the most recent lightVsNoLightSessions for each mouse
conditionYAxes = {'(light)', '(no light)'};
experimentNames = {'obsNoBr', 'obsBr'};
frameEdges = [.336 .415]; % (m) % [(pos at which obs enters frame) (pos at which obs clears hind paws (this is to avoid including tail contacts))]
noBrSessions = 3; % uses the most recent noBrSessions 
brSessions = 7; % uses the first (oldest) brSessions
mouseScatSize = 25;
meanScatSize = 100;

% initializations
xInds = 1:(noBrSessions + brSessions);
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx']);

sessionBins = ismember(sessionInfo.mouse, mice) &...
              ismember(sessionInfo.experiment, experimentNames) &...
              sessionInfo.include;
sessions = sessionInfo(sessionBins, :);

data = struct(); % stores trial data for all sessions

sessionNum = length(mice);
cmap = winter(sessionNum*2);
condColors = [cmap(end,:); cmap(1,:)];
cmap = {cmap(sessionNum+1:end,:), cmap(1:sessionNum,:)};



% collect data
for i = 1:size(sessions,1)

    % load session data
    load([getenv('OBSDATADIR') 'sessions\' sessions.session{i} '\runAnalyzed.mat'],...
            'obsPositions', 'obsTimes',...
            'obsPixPositions', 'frameTimeStamps',...
            'obsOnTimes', 'obsOffTimes',...
            'obsLightOnTimes', 'obsLightOffTimes',...
            '');
    load([getenv('OBSDATADIR') 'sessions\' sessions.session{i} '\run.mat'], 'breaks', 'touch');
    
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes); % !!! very important to determine whether this is in fact a good idea
    
    
    % remove brief touches
    validLengthInds = (touchOffTimes - touchOnTimes) >= minTouchTime;
    touchOnTimes = touchOnTimes(validLengthInds);
    touchOffTimes = touchOffTimes(validLengthInds);
    
    
    % get touch positions and ensure all touches fall within frame
    touchPositions = interp1(obsTimes, obsPositions, touchOnTimes, 'linear');
    validPosInds = touchPositions>frameEdges(1) & touchPositions<frameEdges(2);
    touchOnTimes = touchOnTimes(validPosInds);
    touchOffTimes = touchOffTimes(validPosInds);
           
    
    
    isAvoided = nan(length(obsOnTimes), 1);
    isLightOn = false(length(obsOnTimes), 1);
    
    % iterate over all trials
    for j = 1:length(obsOnTimes)
        
        % find whether and where obstacle was toucheed
%         isAvoided(j) = ~any(touchOnTimes>obsOnTimes(j) & touchOnTimes<obsOffTimes(j));
        isAvoided(j) = ~any(touchOnTimes>obsOnTimes(j) & touchOnTimes<obsOffTimes(j));
        
        % find whether light was on
        isLightOn(j) = min(abs(obsOnTimes(j) - obsLightOnTimes)) < 1; % did the light turn on near whether the obstacle turned on
                
    end
    
    data(i).mouse = sessions.mouse{i};
    data(i).lightOnAvoidance = sum(isAvoided(isLightOn)) / sum(isLightOn);
    data(i).lightOffAvoidance = sum(isAvoided(~isLightOn)) / sum(~isLightOn);
    data(i).isWheelBreak = strcmp(sessions.experiment{i}, 'obsBr');
    
end


% determine which sessions to include (get last noBrSessions without wheel break and first brSessions with wheel break)
temp = num2cell(false(length(data)));
[data(1:length(data)).includeSessions] = temp{:};

for i = 1:length(mice)
    
    % get no break session nums
    indsNoBr = find(strcmp(mice{i}, {data.mouse}) & ~[data.isWheelBreak], noBrSessions, 'last');
    indsBr = find(strcmp(mice{i}, {data.mouse}) & [data.isWheelBreak], brSessions, 'first');
    temp = num2cell(ones(1,length([indsNoBr indsBr])));
    [data([indsNoBr indsBr]).includeSessions] = temp{:};
    
end

% ---------------
% plot everything
% ---------------

% prepare figure
figure('name', 'obsAvoidanceLearningSummary', 'menubar', 'none', 'units', 'pixels', 'position', [500 200 600 650], 'color', [1 1 1]);
fields = {'lightOnAvoidance', 'lightOffAvoidance'};

% plot light on and light off avoidance for each mouse

allAvoidanceData = nan(length(mice), (noBrSessions+brSessions), 2); % (mice, session, isLightOn)

for i = 1:2
    
    subplot(3,1,i)
    
    for j = 1:length(mice)
        
        bins = strcmp(mice{j}, {data.mouse}) & [data.includeSessions];
        
        plot(xInds(1:sum(bins)), [data(bins).(fields{i})], 'color', cmap{i}(j,:)); hold on
        scatter(xInds(1:sum(bins)), [data(bins).(fields{i})], mouseScatSize, cmap{i}(j,:), 'filled')
        
        allAvoidanceData(j, 1:sum(bins), i) = [data(bins).(fields{i})];
        
    end
    
    % pimp fig
    set(gca, 'box', 'off', 'xtick', xInds, 'xlim', [xInds(1)-.5 xInds(end)], 'ylim', [0 1])
    if i==2; xlabel('session', 'fontweight', 'bold'); end
    ylabel({'fraction avoided', conditionYAxes{i}}, 'fontweight', 'bold')
    line([noBrSessions+.5 noBrSessions+.5], [0 1], 'lineWidth', 2, 'color', get(gca, 'xcolor'))
    
end

% plot means
for i = 1:2
    
    subplot(3,1,i)
    meanAvoidance = nanmean(squeeze(allAvoidanceData(:,:,i)),1);
    plot(xInds, meanAvoidance, 'lineWidth', 3, 'color', get(gca, 'xcolor'))
    scatter(xInds, meanAvoidance, meanScatSize, get(gca, 'xcolor'), 'filled')
    
end

% plot light vs no light
lightMeans = nan(length(mice), 1);
lightOffMeans = nan(length(mice), 1);

for i = 1:length(mice)
    inds = find(strcmp(mice{i}, {data.mouse}), lightVsNoLightSessions, 'last');
    lightMeans(i) = nanmean([data(inds).lightOnAvoidance]);
    lightOffMeans(i) = nanmean([data(inds).lightOffAvoidance]);
end


subplot(3,1,3)
centers = [-.5 .5];
spread = .25;
border = .5;
jitterPos = linspace(-spread/2, spread/2, length(mice));
jitterPos = jitterPos(randperm(length(jitterPos)));
allMeans = {lightMeans, lightOffMeans};

for i = 1:length(mice)
    line(centers+jitterPos(i), [lightMeans(i) lightOffMeans(i)], 'linewidth', 1, 'color', [.5 .5 .5]); hold on
end

for i = 1:2
    scatter(centers(i) + jitterPos, allMeans{i}, 75, cmap{i}, 'filled'); hold on
    condMean = nanmean(allMeans{i});
    line([centers(i)-spread centers(i)+spread], [condMean condMean], 'linewidth', 3, 'color', get(gca, 'xcolor'));
end

posTemp = get(gca, 'position');
set(gca, 'xlim', [centers(1)-border centers(2)+border], 'ylim', [0 1], 'box', 'off',...
    'position', [posTemp(1) posTemp(2) .25 posTemp(4)], 'xtick', centers, 'xticklabel', {'light', 'no light'}, 'xcolor', [0 0 0])
ylabel({'fraction avoided', conditionYAxes{i}}, 'fontweight', 'bold')

% save fig
savefig([getenv('OBSDATADIR') 'figures\obsAvoidanceLearningSummary.fig'])










