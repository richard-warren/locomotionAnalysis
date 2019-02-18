function plotSpeedVsPosition(sessionInfo, figTitle, conditions)



% settings
plotMouseAvgs = false;
yLims = [.3 .6];
xLims = [-.5 .2]; % (m)
posRes = .001; % (m) // resoultion of positional grid that velocities are computed over
errorFcn = @(x) nanstd(x)/sqrt(size(x,1)); % function for error bars
% errorFcn = @(x) nanstd(x);

% initializations
if ~exist('conditions', 'var'); conditions = unique(sessionInfo.condition); end
mice = unique(sessionInfo.mouse);
posInterp = xLims(1) : posRes : xLims(2); % velocities will be interpolated across this grid of positional values
lightConditions = {'light off', 'light on'};
colors = hsv(length(conditions));

% collect data for each mouse in each condition for both light on and light off trials
mouseAvgs = nan(length(conditions), length(mice), 2, length(posInterp)); % condition X mouse X light off/on X position
for i = 1:length(conditions)
    for j = 1:length(mice)
        for k = 1:2
            bins = strcmp({sessionInfo.condition}, conditions{i}) & ...
                   strcmp({sessionInfo.mouse}, mice{j}) & ...
                   [sessionInfo.isLightOn]+1 == k & ...
                   ~[sessionInfo.isWheelBreak];
            
            mouseAvgs(i,j,k,:) = mean(cat(1,sessionInfo(bins).trialVelInterp));
        end
    end
end


% plot for light on/off
figure('name', figTitle, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 500], 'inverthardcopy', 'off')

for i = 1:2
    subplot(2,1,i)
    
    for j = 1:length(conditions)
        shadedErrorBar(posInterp, squeeze(mouseAvgs(j,:,i,:)), {@nanmean, errorFcn}, ...
            'lineprops', {'linewidth', 3, 'color', colors(j,:)}, 'patchSaturation', .1); hold on;
        
        if plotMouseAvgs
            for k = 1:length(mice)
                plot(posInterp, squeeze(mouseAvgs(j,k,i,:)), ...
                    'LineWidth', 1, 'Color', [colors(j,:) .4]); hold on
            end
        end
        
        % vertical line at avg wisk contact position for condition
        contactPos = nanmean([sessionInfo([sessionInfo.isLightOn]+1==j & ~[sessionInfo.isWheelBreak]).wiskContactPositions]);
        line([contactPos contactPos], yLims, 'color', [colors(j,:) .5])
    end
    
    % vertical line at avg obs on position
    obsOnPos = nanmean([sessionInfo(~[sessionInfo.isWheelBreak]).obsOnPositions]);
    line([obsOnPos obsOnPos], yLims, 'color', mean(colors,1))
    
    % pimp fig
    set(gca, 'YLim', yLims, 'XLim', xLims, 'Box', 'off')
    title(lightConditions{i})
end

xlabel('osbtacle distance to nose (m)')
ylabel('velocity (m/s)')
for i = 1:length(conditions); lines(i) = plot([nan nan], 'color', colors(i,:), 'LineWidth', 2); end % create dummy lines
legend(lines, conditions, 'Location', 'northeast', 'Box', 'off')










