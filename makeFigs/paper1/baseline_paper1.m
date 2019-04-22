%% load experiment data
fprintf('loading...'); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')


%% ----------
% PLOT THINGS
%  ----------

%% speed vs. position / time plots

yLims = [0 .55];
flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions' ...
    'velContinuousAtContact', 'velVsPosition', 'isWheelBreak', 'wiskContactPosition', 'isBigStep'});
flat = flat(~[flat.isWheelBreak]);
flat = flat([flat.isLightOn]);

% speed vs position, 
figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 300], 'inverthardcopy', 'off')



% plotDvPsth(flat, 'velVsPosition', [-.5 .2], 'isLightOn')
% line(repmat(nanmean([flat.obsOnPositions]),1,2), yLims, 'color', [.5 .5 .5])
% line([0 0], yLims, 'color', [.5 .5 .5])
rectangle
set(gca, 'YLim', yLims, 'YTicks', 0:.2:1);
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')


% savefig(fullfile(getenv('OBSDATADIR'), 'figures', 'baseline', 'baseline_speed.fig'))

