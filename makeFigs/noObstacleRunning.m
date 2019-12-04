%% todo: load previous sessions from spreadsheet rather than computing from scratch

%% find no obstacle sessions
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'Sheet', 'sessions');
sessions = dir(fullfile(getenv('OBSDATADIR'), 'sessions'));
sessions = {sessions([sessions.isdir]).name};  % restrcit to folders
sessions = sessions(cellfun(@(x) length(regexp(x, '\d\d\d\d\d\d_\d\d\d', 'match'))==1, sessions));  % restrict to folders with the proper formatting

sessionData = struct('mouse', 'session');
dataInd = 1;
for i = 1:length(sessions)
    try
        folder = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i});
        load(fullfile(folder, 'run.mat'), 'obsOn')
        if obsOn.length==0
            fprintf('adding session %i: %s...\n', dataInd, sessions{i})
            sessionData(dataInd).mouse = sessionInfo.mouse{strcmp(sessionInfo.session, sessions{i})};
            sessionData(dataInd).session = sessions{i};
            dataInd = dataInd+1;
        end
    catch
        fprintf('WARNING! Problem with session %i: %s...\n', dataInd, sessions{i})
    end
        
end

writetable(struct2table(sessionData), fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'noObstacleSessions.csv'));
disp('all done!')

%% compute velocity (with loops)

% settings
xGridPoints = 200;
velTime = .05;  % (s)
sessions = {sessionData.session};
xLims = [-5.4, 0];  % would be better to compute this emperically


mice = unique({sessionData.mouse});
vels = nan(length(mice), xGridPoints);

for i = 1:length(mice)
    fprintf('computing vels for %s\n', mice{i})
    mouseSessions = {sessionData(strcmp({sessionData.mouse}, mice{1})).session};
    mouseVels = cell(1, length(mouseSessions));  % one row per trial across all sessions
    
    for j = 1:length(mouseSessions)
        load(fullfile(getenv('OBSDATADIR'), 'sessions', mouseSessions{j}, 'runAnalyzed.mat'), ...
            'wheelPositions', 'wheelTimes', 'rewardTimes')
        mouseVels{j} = nan(length(rewardTimes)-1, xGridPoints);
        sessionVel = getVelocity(wheelPositions, velTime, 1/median(diff(wheelTimes)));
        
        for k = 1:length(rewardTimes)-1
            bins = wheelTimes>rewardTimes(k) & wheelTimes<rewardTimes(k+1);
            vel = sessionVel(bins);
            pos = wheelPositions(bins);
            pos = pos(end) - pos;  % express as distance to reward
            
            % remove duplicate positional values (would be better to average all values in a particular bin)
            [pos, uniqueInds] = unique(pos, 'stable');
            vel = vel(uniqueInds);

            mouseVels{j}(end+1,:) = interp1(pos, vel, linspace(pos(1), pos(end), xGridPoints));
        end
    end
    vels(i,:) = nanmean(cat(1, mouseVels{:}), 1);
end
disp('all done!')

%% compute velocity struct (which can be used as input to plotDvPsth)

% settings
xGridPoints = 200;
velTime = .05;  % (s)
sessions = {sessionData.session};
xLims = [-5.2, 0];  % would be better to compute this emperically

data = struct();
dataInd = 1;
mice = unique({sessionData.mouse});
vels = nan(length(mice), xGridPoints);
x = linspace(xLims(1), xLims(2), xGridPoints);

for i = 1:length(sessions)
    
    fprintf('computing vels for %s\n', sessions{i})
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
        'wheelPositions', 'wheelTimes', 'rewardTimes')
    sessionVel = getVelocity(wheelPositions, velTime, 1/median(diff(wheelTimes)));
    mouse = sessionInfo.mouse{strcmp(sessionInfo.session, sessions{i})};

    for k = 1:length(rewardTimes)-1
        bins = wheelTimes>rewardTimes(k) & wheelTimes<rewardTimes(k+1);
        vel = sessionVel(bins);
        pos = wheelPositions(bins);
        pos = pos - pos(end);  % express as distance to reward

        % remove duplicate positional values (would be better to average all values in a particular bin)
        [pos, uniqueInds] = unique(pos, 'stable');
        vel = vel(uniqueInds);
        vel = interp1(pos, vel, x);
        
        if range(pos) > (.95*range(xLims))  % only if the reward interval looks right
            data(dataInd).mouse = mouse;
            data(dataInd).session = sessions{i};
            data(dataInd).velX = x;
            data(dataInd).vel = vel;
            dataInd = dataInd + 1;
        end
    end
end
disp('all done!')


%% plot that ish

close all
close all; figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 300], 'inverthardcopy', 'off'); hold on
velData = plotDvPsth(data, 'vel', [], ...
    {'plotMouseAvgs', true, 'showLegend', false, 'conditionColors', [0 0 0], 'errorFcn', @(x) nanstd(x), 'xlim', [-.5 .2]});
    















