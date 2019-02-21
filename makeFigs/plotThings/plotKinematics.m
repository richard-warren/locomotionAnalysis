function plotKinematics(trajectories, obsHgts, conditions, conditionNames)

% plots kinematic trajectories for different conditions on top of one
% another, including obstacle position // trajectories is [number of trials
% X 2 (xz) X time] matrix with kinematics // obsHgts is height of obstacle
% for each trial // conditions are condition number for each trial, and
% conditionNames is names of conditions for legend

% temp
% flat = getNestedStructFields(data, {'mouse', 'session', 'conditionNum', 'trial', 'isLightOn', 'obsHgt', ...
%     'stepOverKin', 'isLeading', 'isFore'});
% flat = flat([flat.isFore] & ~[flat.isLeading]);
% trajectories = cellfun(@(x) x([1,3],:), {flat.stepOverKin}, 'UniformOutput', false);
% trajectories = cat(3, trajectories{:});
% trajectories = permute(trajectories, [3 1 2]);
% obsHgts = [flat.obsHgt];
% conditions = discretize(obsHgts, 3);
% conditionNames = {'low', 'med', 'high'};
% close all; figure('Position', [1943 618 687 252], 'menubar', 'none', 'color', 'white')

% settings
obsRadius = 3.175 / 2 / 1000; % (m)
colors = hsv(max(conditions));



for i = 1:max(conditions)
    
    % plot avg trajectory for condition
    kinMean = squeeze(nanmean(trajectories(conditions==i,:,:), 1));
    plot(kinMean(1,:), kinMean(2,:), ...
        'LineWidth', 2, 'Color', colors(i,:)); hold on
    
    % draw obstacle
    z = nanmean(obsHgts(conditions==i));
    rectangle('position', [0-obsRadius, z-2*obsRadius, 2*obsRadius, 2*obsRadius], ...
        'curvature', [1 1], 'facecolor', [colors(i,:) 1], 'edgecolor', 'none');    
end


% pimp fig
line(get(gca, 'XLim'), [0 0], 'color', 'black') % add line at top of wheel
set(gca, 'DataAspectRatio', [1 1 1], 'YLim', max(get(gca, 'YLim'), 0), ...
    'XColor', 'none', 'YColor', 'none')
if exist('conditionNames', 'var'); legend(conditionNames, 'box', 'off', 'Location', 'northeast'); end




