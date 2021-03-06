function plotKinematics(trajectories, obsHgts, conditions, varargin)

% plots kinematic trajectories for different conditions on top of one
% another, including obstacle position // trajectories is [number of trials
% X 2 (xz) X time] matrix with kinematics // obsHgts is height of obstacle
% for each trial // conditions are condition number for each trial // opts is cell array containing name value pairs
% of settings to adjust // note that all settings are in structure 's'

% TO DO: check if my error bars are correct... im computing variability in
% z without over the avg x grid... is this correct?


% settings
obsRadius = 3.175 / 2 / 1000; % (m)
s.colors = 'hsv';  % colors for kinematic traces
s.obsColors = [];
s.conditionNames = {}; % user can specify this by passing in via 'opts'
s.obsAlpha = .8; % transparency of obstacles
s.lineAlpha = 1; % transparency of plot lines
s.trialAlpha = .2; % transparency of individual trial lines
s.lineWidth = 2;
s.mouseNames = {}; % if this exists (is set in opts), kinematics are first averaged within, then across mice // cell array of mouse name corresponding to each trial
s.errorFcn = []; % if this is set to an anonymous function, e.g. "@(x) nanstd(x)/sqrt(size(x,1))", then error bars are plotted over traces
s.errorAlpha = .2; % transparency of error bars
s.isBotView = false; % if true, then plots are assumed to plot kinematics from the bottom view (xy)
s.trialsToOverlay = []; % if not empty, plot individual trials for each condition
s.yLimZero = true; % makes min of y axis zero so traces don't dip beneathe 'floor' on plots
s.lineColor = 'black';
s.plotObs = true;  % whether to plot the obstacle
s.scatterEnds = true;  % whether to add scatter points to the end of the individual trials
s.evenlySampleTrials = true;  % when showing individual trials, whether to evenly sample across conditions // otherwise samples the same number of trials per condition regardless of condition probability


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if ischar(s.colors); s.colors = eval([s.colors '(max(conditions))']); end % set colorspace if color is specified as a string
if ~isempty(s.mouseNames); mice = unique(s.mouseNames); end

% add line at top of wheel
if ~s.isBotView; line([-.2 .2], [0 0], 'color', s.lineColor); end

% draw obstacles for each condition
if isempty(s.obsColors); s.obsColors = s.colors; end
if s.plotObs
    for i = 1:max(conditions)
        if ~s.isBotView
            z = nanmean(obsHgts(conditions==i));
            rectangle('position', [0-obsRadius, z-2*obsRadius, 2*obsRadius, 2*obsRadius], ...
                'curvature', [1 1], 'facecolor', [s.obsColors(i,:) s.obsAlpha], 'edgecolor', [.4 .4 .4]); hold on
        else
            rectangle('position', [0-obsRadius, -.1, 2*obsRadius, .2], ...
                'facecolor', [s.obsColors(i,:) s.obsAlpha], 'edgecolor', 'none'); hold on
        end
    end
end

% plot kinematics for each condition
for i = unique(conditions)'
    
    % first avg within each mouse ONLY if s.mouseNames is provided by user
    if ~isempty(s.mouseNames)
        conditionTrajectories = nan(length(mice), 2, size(trajectories,3));
        for j = 1:length(mice)
            conditionTrajectories(j,:,:) = nanmean(trajectories(conditions==i & strcmp(s.mouseNames, mice{j}),:,:),1);
        end
    else
        conditionTrajectories = trajectories(conditions==i,:,:);
    end
    
    % plot mean trajectory for condition
    kinMean = squeeze(nanmean(conditionTrajectories, 1));
    if isempty(s.errorFcn)
        plot(kinMean(1,:), kinMean(2,:), ...
            'LineWidth', s.lineWidth, 'Color', [s.colors(i,:) s.lineAlpha]); hold on
    else
        shadedErrorBar(kinMean(1,:), squeeze(conditionTrajectories(:,2,:)), {@nanmean, s.errorFcn}, ...
            'lineprops', {'linewidth', s.lineWidth, 'color', [s.colors(i,:) s.lineAlpha]}, 'patchSaturation', s.errorAlpha); hold on;
    end
    
    % plot individual trials
    if ~isempty(s.trialsToOverlay)
        if s.evenlySampleTrials
            conditionTrialsToShow = min(s.trialsToOverlay, size(conditionTrajectories,1)) * mean(conditions==i);
        else
            conditionTrialsToShow = min(s.trialsToOverlay, size(conditionTrajectories,1)) / length(unique(conditions));
        end
        inds = randperm(size(conditionTrajectories,1), floor(conditionTrialsToShow));
        
        for j = inds
            plot(squeeze(conditionTrajectories(j,1,:)), squeeze(conditionTrajectories(j,2,:)), ...
                'LineWidth', 1, 'Color', [s.colors(i,:) s.trialAlpha]); hold on
        end
        if s.scatterEnds
            scatter(conditionTrajectories(inds,1,end), conditionTrajectories(inds,2,end), ...
                20, [s.colors(i,:)], 'filled', 'MarkerFaceAlpha', (1+s.trialAlpha)/2)
        end
    end
end


% pimp fig
set(gca, 'DataAspectRatio', [1 1 1], 'XColor', 'none', 'YColor', 'none')
% set(gca, 'DataAspectRatio', [1 1 1])
if s.yLimZero && ~s.isBotView; set(gca, 'YLim', max(get(gca, 'YLim'),0)); end
if ~isempty(s.conditionNames); legend(s.conditionNames, 'box', 'off', 'Location', 'best'); end

pause(.001)
