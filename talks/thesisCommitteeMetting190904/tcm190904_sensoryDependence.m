%% INITIALIZATIONS

% load global settings
clear all; tcm190910_config;

% load data
fprintf('loading data... '); 
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'sensoryDependence_data.mat'), 'data');
data.data = data.data(~ismember({data.data.mouse}, {'sen1'}));  % exclude sen1
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'isWheelBreak', 'wiskContactPosition', 'obsHgt', 'isPawSuccess', ...
    'isLeading', 'isFore', 'stepOverKinInterp', 'paw', 'controlStepKinInterp', 'preObsHgt', 'sensoryCondition', 'preObsKin'});
disp('data loaded!')

% conditions and conditionals
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore paw', 'hind paw'}});
vars.sensoryCondition = struct('name', 'sensoryCondition', 'levels', {{'WL', 'W', 'L', '-'}}, 'levelNames', {{'W+V', 'W', 'V', '-'}});

conditionals.lightOn = struct('name', 'isLightOn', 'condition', @(x) x==1);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
conditionals.none = struct('name', '', 'condition', @(x) x); % no conditional
conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);


%% SUCCESS

figure('units', 'inches', 'position', [6.98 5.07 3.29 figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
dv = getDvMatrix(data, 'isTrialSuccess', vars.sensoryCondition, {'mouse'}, conditionals.none);
barFancy(dv, 'ylabel', 'success rate', 'colors', sensoryColors, 'YLim', [0 1], barProperties{:}, ...
    'levelNames', {{'', '', '', ''}}, 'YTick', [0 .5 1])
set(gca, 'Position', [.2 .11 .78 .82])
print -clipboard -dmeta

%% SPEED

figure('units', 'inches', 'position', [6.98 5.07 3.29 figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
dv = getDvMatrix(data, 'velAtWiskContact', vars.sensoryCondition, {'mouse'}, conditionals.none);
barFancy(dv, 'ylabel', 'velocity (m/s)', 'colors', sensoryColors, barProperties{:}, ...
    'levelNames', {vars.sensoryCondition.levelNames}, 'YTick', [0 .3 .6])
set(gca, 'Position', [.2 .11 .78 .82])
print -clipboard -dmeta

%% SPEED VS POSITION

% initializations
yLims = [.3 .7];
figure('units', 'inches', 'position', [21.18 3.21 6.28 figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false); hold on
flat_sub = flat(~[flat.isWheelBreak]);

% add obstacle rectangle and lines
x = [nanmean([flat_sub.obsOnPositions]) nanmean([flat_sub.obsOffPositions])];
line([0 0], yLims, 'linewidth', 2, 'color', axisColor)
line([x(1) x(1)], yLims, 'linewidth', 2, 'color', axisColor)
text(x(1), yLims(2), '\itobstacle on', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSize, 'FontName', font)
text(0, yLims(2), '\itobstacle reached', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSize, 'FontName', font)

% speed vs position
plotDvPsth(flat_sub, 'velVsPosition', 'sensoryCondition', ...
    {'showLegend', false, 'conditionColors', sensoryColors([4 3 2 1],:), 'xlim', [-.5 .2], 'showErrorBars', false, ... 
     'plotConditions', vars.sensoryCondition.levels([4 3 2 1]), 'errorAlpha', .1, 'lineWidth', 5})
set(gca, 'YLim', yLims, 'YTick', linspace(yLims(1),yLims(2),3));
set(gca, 'YLim', yLims);
xlabel('distance (m)')
ylabel('velocity (m/s)')
set(gca, 'XColor', axisColor, 'YColor', axisColor, 'Color', 'black', 'FontSize', fontSize, 'FontName', font)
print -clipboard -dmeta

%% KINEMATICS

% settings
figure('units', 'inches', 'position', [10.67 3.16 5.7 7.5], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false); hold on
obsHgtBins = 4; % discretize obstacle heights into this many bins
fading = .5; % within a plot, each paw's color fades from fading*color -> color
xLims = [-.055 0];
yLims = [0 .016];

% initializations
obsHgtDiscretized = num2cell(discretize([flat.obsHgt], linspace(3.4, 10, obsHgtBins+1)/1000));
[flat.obsHgtDiscretized] = obsHgtDiscretized{:};
flat_sub = flat(~isnan([flat.obsHgtDiscretized]) & ...
            ~[flat.isWheelBreak] & ...
            [flat.isLeading] & ...
            [flat.isFore]);
[~, conditions] = ismember({flat_sub.sensoryCondition}, vars.sensoryCondition.levels);
kinData = permute(cat(3, flat_sub.preObsKin), [3,1,2]);
kinDataCtl = permute(cat(3, flat_sub.controlStepKinInterp), [3,1,2]);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1); % change the x starting x position of ctl steps to match steps over


colorsTemp = [sensoryColors(1:end-1,:); .8 .8 .8];
buf = .01;
for i = 1:4
    subplot(4,1,i); hold on
    bins = conditions==i;
    plotColor = repmat(colorsTemp(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    plotKinematics(kinData(bins,[1,3],:), [flat_sub(bins).obsHgt], [flat_sub(bins).obsHgtDiscretized], ...
        {'colors', plotColor, 'obsAlpha', 1, 'mouseNames', {flat_sub(bins).mouse}}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    set(gca, 'XLim', xLims, 'YLim', yLims, 'color', 'black', ...
        'position', [0 1-i*.24 1 .22], 'units', 'normalized')
end

print -clipboard -dmeta

%% SHAPING MOVING AVGS
flat_sub = flat(~isnan([flat.obsHgt]) & ...
                [flat.isLeading] & ...
                [flat.isFore]);
[~, conditions] = ismember({flat_sub.sensoryCondition}, vars.sensoryCondition.levels);
figure('units', 'inches', 'position', [4.46 1.36 4.25 figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
plot([0 100], [0 100], 'Color', [1 1 1]*.6, 'LineWidth', 3) % add unity line
logPlotRick([flat_sub.obsHgt]*1000, [flat_sub.preObsHgt]*1000, ...
    {'colors', sensoryColors, 'conditions', conditions, 'xlabel', 'obstacle height', 'ylabel', 'paw height (mm)',  ...
    'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat_sub.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'lineWidth', 3, 'errorAlpha', .4})
set(gca, 'color', 'black', 'xcolor', axisColor, 'YColor', axisColor, 'FontSize', fontSize, 'FontName', font)
print -clipboard -dmeta

%% SHAPING CORRELATIONS

figure('units', 'inches', 'position', [7.23 5.15 4.25 figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
corrs = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, vars.sensoryCondition, {}, {'mouse'}, ...
    [conditionals.isLeading; conditionals.isFore], 'corr');
barFancy(corrs, 'ylabel', 'paw obstacle correlation', 'colors', sensoryColors, barProperties{:}, ...
    'levelNames', {vars.sensoryCondition.levelNames}, 'YTick', -.2:.2:.6)
print -clipboard -dmeta









