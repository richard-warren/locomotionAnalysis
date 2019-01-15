function plotSpeedAtWiskContact(sessionInfo, figTitle, conditions)


% to do: make work with condition number restriction // change data structs
% to make collecting data faster

% settings
yLims = [.2 .6];
timePrePost = [-.5 .5];
velocityDelta = .01;  % compute velocity over this interval
errorFcn = @(x) nanstd(x)/sqrt(size(x,1)); % function for error bars
plotMouseAvgs = false;
% errorFcn = @(x) nanstd(x);


% initializations
load(fullfile(getenv('OBSDATADIR'), 'sessions', sessionInfo.session{1}, 'runAnalyzed.mat'), 'wheelTimes')
dt = median(diff(wheelTimes)); % temporal resolution for wheel velocity
if ~exist('conditions', 'var'); conditions = unique(sessionInfo.condition(logical(sessionInfo.include))); end
colors = hsv(length(conditions));
lightConditions = {'light off', 'light on'};
mice = unique(sessionInfo.mouse);
times = timePrePost(1):dt:timePrePost(2);



% make struct where each row is a trial and contains vel surrounding
% moment of whisker contact
data = cell(0);
dataInd = 1;

for i = 1:height(sessionInfo)
    
    if sessionInfo.include(i)
        load(fullfile(getenv('OBSDATADIR'), 'sessions', sessionInfo.session{i}, 'runAnalyzed.mat'), ...
            'wiskContactTimes', 'wheelTimes', 'wheelPositions', 'isLightOn')
        wheelVel = getVelocity(wheelPositions, velocityDelta, 1/dt);

        for j = 1:length(wiskContactTimes)
            if ~isnan(wiskContactTimes(j))
                startInd = find(wheelTimes>=wiskContactTimes(j)+timePrePost(1), 1, 'first');
                contactVel = wheelVel(startInd:startInd+length(times)-1);

                if length(contactVel)~=length(times); keyboard; end

                data(dataInd).session = sessionInfo.session{i};
                data(dataInd).condition = sessionInfo.condition{i};
                data(dataInd).mouse = sessionInfo.mouse{i};
                data(dataInd).isLightOn = isLightOn(j);
                data(dataInd).contactVel = contactVel;
                dataInd = dataInd+1;
            end
        end
    end
end



% collect data for each mouse in each condition for both light on and light off trials
mouseAvgs = nan(length(conditions), length(mice), 2, length(times)); % condition X mouse X light off/on X position
for i = 1:length(conditions)
    for j = 1:length(mice)
        for k = 1:2
            bins = strcmp({data.condition}, conditions{i}) & ...
                   strcmp({data.mouse}, mice{j}) & ...
                   [data.isLightOn]+1 == k;
            
            mouseAvgs(i,j,k,:) = mean(cat(1,data(bins).contactVel));
        end
    end
end


% plot for light on/off
figure('name', figTitle, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 500], 'inverthardcopy', 'off')

for i = 1:2
    subplot(2,1,i)
    
    for j = 1:length(conditions)
        shadedErrorBar(times, squeeze(mouseAvgs(j,:,i,:)), {@nanmean, errorFcn}, ...
            'lineprops', {'linewidth', 3, 'color', colors(j,:)}, 'patchSaturation', .1); hold on;
        
        if plotMouseAvgs
            for k = 1:length(mice)
                plot(times, squeeze(mouseAvgs(j,k,i,:)), ...
                    'LineWidth', 1, 'Color', [colors(j,:) .4]); hold on
            end
        end
    end
    line([0 0], yLims, 'color', mean(colors,1))
    
    % pimp fig
    set(gca, 'YLim', yLims, 'XLim', timePrePost, 'Box', 'off')
    title(lightConditions{i})
end

xlabel('time from whisker contact (s)')
ylabel('velocity (m/s)')
for i = 1:length(conditions); lines(i) = plot([nan nan], 'color', colors(i,:), 'LineWidth', 2); end % create dummy lines
legend(lines, conditions, 'Location', 'northeast', 'Box', 'off')











