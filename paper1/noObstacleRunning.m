%% find no obstacle sessions (you can skip this by loading computed sessions below)
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'Sheet', 'sessions');
mouseInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'Sheet', 'mice');
sessions = dir(fullfile(getenv('OBSDATADIR'), 'sessions'));
sessions = {sessions([sessions.isdir]).name};  % restrict to folders
sessions = sessions(cellfun(@(x) length(regexp(x, '\d\d\d\d\d\d_\d\d\d', 'match'))==1, sessions));  % restrict to folders with the proper formatting

%%
sessionData = struct('mouse', 'session');
dataInd = 1;
for i = 1:length(sessions)
    try
        % find genotype
        mouse = sessionInfo.mouse{strcmp(sessionInfo.session, sessions{i})};
        mouseInfoInd = find(strcmp(mouseInfo.mouse, mouse));
        if isempty(mouseInfoInd)
            genotype = 'C57/B6';  % assume wild-type if not listed in spreadsheet
        else
            genotype = mouseInfo.line_Genotype{mouseInfoInd};
        end
        
        if strcmp(genotype, 'C57/B6')
            folder = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i});
            load(fullfile(folder, 'run.mat'), 'obsOn')
            if obsOn.length==0
                fprintf('adding session %i: %s...\n', dataInd, sessions{i})
                sessionData(dataInd).mouse = mouse;
                sessionData(dataInd).session = sessions{i};
                dataInd = dataInd+1;
            end
        else
            fprintf('skipping %s because genoype is %s\n', mouse, genotype);
        end
        
    catch
        fprintf('WARNING! Problem with session %i: %s...\n', dataInd, sessions{i})
    end
end

writetable(struct2table(sessionData), fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'noObstacleSessions.csv'));
disp('all done!')

%% load obstacle sessions

sessionData = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'noObstacleSessions.csv'), 'Delimiter', ',');
sessionData = table2struct(sessionData);

%% compute velocity struct (which can be used as input to plotDvPsth)
% (you can instead use the next block to load previously computed data)

% settings
xGridPoints = 200;
velTime = .01;  % (s)
xLims = [-5.5, 1];  % meters before and after reward delivery
smoothing = .5;  % (m)

data = struct();
dataInd = 1;
mice = unique({sessionData.mouse});
vels = nan(length(mice), xGridPoints);
x = linspace(xLims(1), xLims(2), xGridPoints);
smoothSmps = round(smoothing / mean(diff(x)));

for i = 1:length(sessionData)
    fprintf('computing vels for %s\n', sessionData(i).session)
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessionData(i).session, 'runAnalyzed.mat'), ...
        'wheelPositions', 'wheelTimes', 'rewardTimes')
    sessionVel = getVelocity(wheelPositions, velTime, 1/median(diff(wheelTimes)));

    for k = 1:length(rewardTimes)-1
        rewardPos = wheelPositions(find(wheelTimes>=rewardTimes(k),1,'first'));
        pos = wheelPositions - rewardPos;  % wheel positions shifted relative to position of current reward
        bins = pos>=xLims(1) & pos<=xLims(2);
        pos = pos(bins);
        vel = sessionVel(bins);

        % remove duplicate positional values (would be better to average all values in a particular bin)
        [pos, uniqueInds] = unique(pos, 'stable');
        vel = vel(uniqueInds);
        vel = interp1(pos, vel, x);
        vel = smooth(vel, smoothSmps);
        
        if range(pos) > (.9*range(xLims))  % only if the reward interval looks right
            data(dataInd).mouse = sessionData(i).mouse;
            data(dataInd).session = sessionData(i).session;
            data(dataInd).velX = x;
            data(dataInd).vel = vel;
            dataInd = dataInd + 1;
        end
    end
end
disp('all done!')

fprintf('saving...'); save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'noObstacleRunning_data.mat'), 'data'); disp('data saved!')

%% load previously computed data
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'noObstacleRunning_data.mat'))

%% plot that ish

yLims = [0 .85];
xLims = [-5.4 0];

global_config;  % load global settings
figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 300], 'inverthardcopy', 'off'); hold on
velData = plotDvPsth(data, 'vel', [], ...
    {'plotMouseAvgs', true, 'showLegend', false, 'conditionColors', [0 0 0], ...
    'errorFcn', @(x) nanstd(x), 'mouseAlpha', .3});
set(gca, 'xlim', xLims, 'xtick', -5:0, 'ylim', yLims, 'ytick', 0:.2:.8)
xlabel('distance to water reward (m)')
ylabel('velocity (m/s)')
line([0 0], [yLims(2)-.08*range(yLims) yLims(2)], 'linewidth', 3, 'color', waterColor)
text(0, .85, 'water', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'noObsRunning');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

% compute running speed statistics
vels = velData(:, floor(length(data(1).velX)/2));  % sample velocity in middle of x axis across all mice
fprintf('velocity: %.2f +- %.2f SEM\n', mean(vels), std(vels)/sqrt(length(vels)))













