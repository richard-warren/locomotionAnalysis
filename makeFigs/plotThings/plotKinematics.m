function plotKinematics(trajectories, obsHgts, conditions, opts)

% plots kinematic trajectories for different conditions on top of one
% another, including obstacle position // trajectories is [number of trials
% X 2 (xz) X time] matrix with kinematics // obsHgts is height of obstacle
% for each trial // conditions are condition number for each trial, and
% conditionNames is names of conditions for legend // opts is cell array containing name value pairs
% of settings to adjust // note that all settings are in structure 's'


% settings
obsRadius = 3.175 / 2 / 1000; % (m)
s.colors = 'hsv';
s.conditionNames = {}; % user can specify this by passing in via 'opts'
s.obsAlpha = .8;
s.lineAlpha = 1;

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end

% initializations
if ischar(s.colors); s.colors = eval([s.colors '(max(conditions))']); end % set colorspace if color is specified as a string

for i = 1:max(conditions)
    
    % plot median trajectory for condition
    kinMedian = squeeze(nanmean(trajectories(conditions==i,:,:), 1));
    plot(kinMedian(1,:), kinMedian(2,:), ...
        'LineWidth', 2, 'Color', [s.colors(i,:) s.lineAlpha]); hold on
    
    % draw obstacle
    z = nanmean(obsHgts(conditions==i));
    rectangle('position', [0-obsRadius, z-2*obsRadius, 2*obsRadius, 2*obsRadius], ...
        'curvature', [1 1], 'facecolor', [s.colors(i,:) s.obsAlpha], 'edgecolor', 'none');
end


% pimp fig
line([-.1 .1], [0 0], 'color', 'black') % add line at top of wheel
set(gca, 'DataAspectRatio', [1 1 1], 'YLim', max(get(gca, 'YLim'), 0), ...
    'XColor', 'none', 'YColor', 'none')
if ~isempty(s.conditionNames); legend(s.conditionNames, 'box', 'off', 'Location', 'best'); end




