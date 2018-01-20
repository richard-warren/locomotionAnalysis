function speedOverTimeSummary(mice)

% shows average speed as mice jump over obstacle for each session, showing (hopefully) that mice get faster over time
% one plot for light, one for no light
% also shows mean light vs no light for most recent sessions
% starts with the first session with wheel break
% note: for each session MEDIAN speed is computed (to minimize impact of very slow trials) - then averaging is preformed across sessions
%
% input         mice:      name of mice to analyze

% !!! right now if a session is not included it assumes the two adjacent sessions are on two adjacent days, when in fact they are separated by a day... should fix this...

% user settings
sessionsToPlot = 7;
lightVsNoLightSessions = 3; % for the plot comparing light vs no light speed, only take the most recent lightVsNoLightSessions for each mouse
obsPrePost = [.6 .25]; % plot this much before and after the obstacle reaches the mouse
posRes = .001; % resolution of x axis, in meters
yLims = [.1 .6]; % m/s
obsPos = .382; % m, position at which obstacle is in the middle of the frame // use getFrameTimes function to determine this value
conditionYAxes = {'(light)', '(no light)', ''};

% initializations
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx']);

sessionBins = ismember(sessionInfo.mouse, mice) &...
              ismember(sessionInfo.experiment, {'obsBr'}) &...
              sessionInfo.include;
sessions = sessionInfo(sessionBins, :);

posInterp = -obsPrePost(1) : posRes : obsPrePost(2); % velocities will be interpolated across this grid of positional values

data = struct(); % stores trial data for all sessions



% collect data
for i = 1:size(sessions,1)

    [sessionVels, isLightOn, obsOnPositions] = getSessionSpeedInfo(sessions.session{i}, posInterp, obsPos);
    
    data(i).mouse = sessions.mouse{i};
    data(i).lightOnVel = nanmedian(sessionVels(isLightOn,:),1);
    data(i).lightOffVel = nanmedian(sessionVels(~isLightOn,:),1);
    data(i).avgObsOnPos = nanmean(obsOnPositions);
    
end

% determine which sessions to include (get the first sessionsToPlot sessions with the wheel break)
for i = 1:length(mice)
    
    inds = find(strcmp(mice{i}, {data.mouse}), sessionsToPlot, 'first');
    temp = num2cell(1:min(length(inds), sessionsToPlot));
    [data(inds).sessionNum] = temp{:};
    
end

sessionNum = max([data.sessionNum]);
cmap = winter(sessionNum*2);
condColors = [cmap(end,:); cmap(1,:)];
cmap = {cmap(sessionNum+1:end,:), cmap(1:sessionNum,:)};

% plot everything

% prepare figure
figure('name', 'speedOverTimeSummary', 'menubar', 'none', 'units', 'pixels', 'position', [500 200 550 800], 'color', [1 1 1]);
fields = {'lightOnVel', 'lightOffVel'};

% plot light on and light off avoidance for each mouse

for i = 1:3
    
    subplot(3,1,i)
    
    % plot individual conditions over time
    if i<3
        for j = 1:max([data.sessionNum])

            bins = [data.sessionNum] == j;
            sessionVels = reshape([data(bins).(fields{i})], length(posInterp), sum(bins))';

            plot(posInterp, nanmean(sessionVels,1), 'color', cmap{i}(j,:), 'linewidth', 3); hold on

        end
        
    % plot condition comparison    
    else
        
        inds = nan(lightVsNoLightSessions * length(mice), 1);

        for j = 1:length(mice)
            indInds = (j-1)*lightVsNoLightSessions+1 : j*lightVsNoLightSessions;
            inds(indInds) = find(strcmp(mice{j}, {data.mouse}), lightVsNoLightSessions, 'last');
        end
        
        for j = 1:2
            allVels = reshape([data(inds).(fields{j})], length(posInterp), length(inds))';
            plot(posInterp, nanmean(allVels,1), 'color', condColors(j,:), 'linewidth', 3); hold on
        end
    end
    
    % pimp fig
    set(gca, 'box', 'off', 'xlim', [posInterp(1) posInterp(end)], 'ylim', yLims)
    if i==3; xlabel('position (m)', 'fontweight', 'bold'); end
    ylabel({'speed (m/s)', conditionYAxes{i}}, 'fontweight', 'bold')
    obsOnPos = -mean([data.avgObsOnPos]);   
    line([obsOnPos obsOnPos], yLims, 'color', get(gca, 'xcolor'), 'linewidth', 2)
    line([0 0], yLims, 'color', get(gca, 'xcolor'), 'linewidth', 2)
    
end



% save fig
savefig([getenv('OBSDATADIR') 'figures\speedOverTimeSummary.fig'])








