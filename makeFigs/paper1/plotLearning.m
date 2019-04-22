

% settings
yLims = [.2 .6];



% load experiment data
disp('loading...'); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'learning_data.mat'), 'data'); disp('learning data loaded!')




flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'condition', 'conditionNum', 'sessionNum', ...
    'isLightOn', 'obsOnPositions', 'velVsPosition', 'isWheelBreak', 'isTrialSuccess'});
mice = {data.mouse};
colors = lines(length(mice));


%% speed vs position
figure('name', 'learning', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 700 600], 'inverthardcopy', 'off')
subplot(2,1,1)
plotDvPsth(flat(strcmp({flat.condition}, 'wheelBreaks')), 'velVsPosition', [-.5 .2], 'conditionNum')
% line(repmat(nanmean([flat.obsOnPositions]),1,2), yLims, 'color', [.5 .5 .5])
line([0 0], yLims, 'color', [.5 .5 .5])
set(gca, 'YLim', yLims);
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')


%% success
subplot(2,1,2); cla

% compute success per session
success = nan(length(data), max([flat.sessionNum]));
noBrSessions = max([flat(strcmp({flat.condition}, 'noWheelBreaks')).conditionNum]);
brSessions = max([flat(strcmp({flat.condition}, 'wheelBreaks')).conditionNum]);
xTicks = -noBrSessions:brSessions-1;

for i = 1:length(data) % iterate through mice
    for j = 1:length(data(i).sessions)
        bins = strcmp({flat.mouse}, data(i).mouse) & strcmp({flat.session}, data(i).sessions(j).session);
        success(i,j) = nanmean([flat(bins).isTrialSuccess]);
    end
    
    plot(xTicks, success(i,:)); hold on
    scatter(xTicks, success(i,:), 20, colors(i,:), 'filled');
end

% plot means
plot(xTicks, nanmean(success,1), 'linewidth', 3, 'color', 'black')
scatter(xTicks, nanmean(success,1), 50, 'black', 'filled')
for i = 1:size(success,2)
    line([xTicks(i) xTicks(i)], nanmean(success(:,i)) + [-1 1]*nanstd(success(:,i)), ...
        'color', 'black', 'linewidth', 2)
end

set(gca, 'box', 'off')
xlabel('session #')
ylabel('success rate')





