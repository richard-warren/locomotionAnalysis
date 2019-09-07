%% GLOBAL SETTINGS

% first load data // each dataset has a corresponding flattened data set //
% incidual plots may take subsets of the flattened dataset

axisColor = [.95 .95 .95];
colors = [.6 .83 .54; .2 .8 1; 1 .87 .08; .5 .5 .5];
mouseColors = 'jet';

%% INITIALIZATIONS

% load data
fprintf('loading data... '); 
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'sensoryDependence_data.mat'), 'data');
data.data = data.data(~ismember({data.data.mouse}, {'sen1'}));  % exclude sen1
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'isWheelBreak', 'wiskContactPosition', 'obsHgt', 'isPawSuccess', ...
    'isLeading', 'isFore', 'stepOverKinInterp', 'paw', 'controlStepKinInterp', 'preObsHgt', 'sensoryCondition', 'preObsKin'});
disp('data loaded!')

% conditions and conditionals
conditionNames = {{'fore', 'hind'}, {'leading', 'trailing'}};
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

figure('units', 'inches', 'position', [6.98 5.07 3.29 2.90], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)

dv = getDvMatrix(data, 'isTrialSuccess', vars.sensoryCondition, {'mouse'}, conditionals.none);
barFancy(dv, 'axisColor', axisColor, 'ylabel', 'success rate', ...
    'colors', colors, 'scatterAlpha', .8, 'barAlpha', .4, 'labelSizePerFactor', .1, 'lineThickness', 2, 'YLim', [0 1])
set(gca, 'Position', [.2 .11 .78 .82])
print -clipboard -dmeta


%% SPEED VS POSITION

% initializations
yLims = [.3 .7];
close all
figure('units', 'inches', 'position', [21.18 3.21 6.28 3.20], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false); hold on

% add obstacle rectangle and lines
x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
line([0 0], yLims, 'linewidth', 2, 'color', axisColor)
line([x(1) x(1)], yLims, 'linewidth', 2, 'color', axisColor)
text(x(1), yLims(2), '\itobstacle on', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
text(0, yLims(2), '\itobstacle reached', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

% speed vs position
plotDvPsth(flat, 'velVsPosition', 'sensoryCondition', ...
    {'showLegend', false, 'conditionColors', colors([4 3 2 1],:), 'xlim', [-.5 .2], 'showErrorBars', false, ... 
     'plotConditions', vars.sensoryCondition.levels([4 3 2 1]), 'errorAlpha', .1, 'lineWidth', 3})
set(gca, 'YLim', yLims, 'YTick', linspace(yLims(1),yLims(2),3));
set(gca, 'YLim', yLims);
xlabel('distance (m)')
ylabel('velocity (m/s)')
set(gca, 'XColor', axisColor, 'YColor', axisColor, 'Color', 'black')
print -clipboard -dmeta

%% KINEMATICS

% settings
figure('units', 'inches', 'position', [10.67 3.16 4.51 6.06], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false); hold on
obsHgtBins = 4; % discretize obstacle heights into this many bins
fading = .5; % within a plot, each paw's color fades from fading*color -> color
xLims = [-.055 0];
yLims = [0 .016];

% initializations
obsHgtDiscretized = num2cell(discretize([flat.obsHgt], linspace(3.4, 10, obsHgtBins+1)/1000));
[flat.obsHgtDiscretized] = obsHgtDiscretized{:};
flat_sub = flat_sub(~isnan([flat_sub.obsHgtDiscretized]) & ...
            ~[flat_sub.isWheelBreak] & ...
            [flat_sub.isLeading] & ...
            [flat_sub.isFore]);
[~, conditions] = ismember({flat_sub.sensoryCondition}, vars.sensoryCondition.levels);
kinData = permute(cat(3, flat_sub.preObsKin), [3,1,2]);
kinDataCtl = permute(cat(3, flat_sub.controlStepKinInterp), [3,1,2]);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1); % change the x starting x position of ctl steps to match steps over


colorsTemp = [colors(1:end-1,:); .8 .8 .8];
buf = .01;
for i = 1:4
    subplot(4,1,i); hold on
    bins = conditions==i;
    plotColor = repmat(colorsTemp(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    plotKinematics(kinData(bins,[1,3],:), [flat_sub(bins).obsHgt], flat_subflat(bins).obsHgtDiscretized], ...
        {'colors', plotColor, 'obsAlpha', 1, 'lineAlpha', .8, 'mouseNames', {flat_sub(bins).mouse}}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    set(gca, 'XLim', xLims, 'YLim', yLims, 'color', 'black', ...
        'position', [0 1-i*.24 1 .22], 'units', 'normalized')
end

print -clipboard -dmeta

%% SHAPING MOVING AVGS
[~, conditions] = ismember({flat.sensoryCondition}, vars.sensoryCondition.levels);
figure('units', 'inches', 'position', [4.46 1.36 3.8 2.8], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
plot([0 xLims(2)], [0 xLims(2)], 'Color', [1 1 1]*.6, 'LineWidth', 3) % add unity line
logPlotRick([flat.obsHgt]*1000, [flat.preObsHgt]*1000, ...
    {'colors', colors, 'conditions', conditions, 'xlabel', 'obstacle height', 'ylabel', 'paw height (mm)',  ...
    'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'lineWidth', 2})
set(gca, 'color', 'black', 'xcolor', axisColor, 'YColor', axisColor)
print -clipboard -dmeta

%% SHAPING CORRELATIONS
figure('units', 'inches', 'position', [7.23 5.15 4.00 3.00], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
corrs = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, vars.sensoryCondition, {}, {'mouse'}, ...
    [conditionals.isLeading; conditionals.isFore], 'corr');
barFancy(corrs, 'axisColor', axisColor, 'ylabel', 'paw obstacle correlation', ...
    'colors', colors, 'scatterAlpha', .8, 'barAlpha', .4, 'labelSizePerFactor', .15, 'lineThickness', 2, 'barWidth', .75, ...
    'ylabelPosX', -.75)
print -clipboard -dmeta









