%% load experiment data
fprintf('loading...'); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')


%% ----------
% PLOT THINGS
%  ----------

%% speed vs. position plot

% settings
obsOnColor = [0 0 0];
mouseColors = lines(length(data));

% initializations
yLims = [0 .55];
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions' ...
    'velContinuousAtContact', 'velVsPosition', 'isWheelBreak', 'wiskContactPosition', 'isBigStep'});
flat = flat(~[flat.isWheelBreak]);
flat = flat([flat.isLightOn]);
mice = unique({flat.mouse});

allVels = cat(1,flat.velVsPosition);
xGrid = allVels(1,:);
allVels = allVels(1:2:end,:); % every other row conditions the x vector


%% collect data
vels = nan(length(mice), size(flat(1).velContinuousAtContact,2));
for i = 1:length(mice)
    
    sessions = unique({flat(strcmp({flat.mouse}, mice{i})).session});
    sesVels = nan(length(mice), size(flat(1).velContinuousAtContact,2);
    for j = 1:length(sessions)
        sesVels(j,:) = flat(strcmp({flat.session}, sessions{j})).velVsPosition;
        
    end
    
end


figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 300], 'inverthardcopy', 'off')

% plotDvPsth(flat, 'velVsPosition', [-.5 .2], 'isLightOn')
% line(repmat(nanmean([flat.obsOnPositions]),1,2), yLims, 'color', [.5 .5 .5])
% line([0 0], yLims, 'color', [.5 .5 .5])

% add shaded area where obstacle is on
x = [nanmean(trialObsOnPosits(:,i),1) nanmean(trialObsOffPosits(:,i),1)];
    rectangle('Position', [x(1) yLims(1) diff(x) diff(yLims)], ...
        'FaceColor', [obsOnColor obsOnAlpha], 'EdgeColor', 'none');

set(gca, 'YLim', yLims, 'YTicks', 0:.2:1);
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')


% savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'baseline', 'baseline_speed.fig'))

