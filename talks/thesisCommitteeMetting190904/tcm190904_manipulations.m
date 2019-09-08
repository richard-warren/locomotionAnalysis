%% INITIALIZATIONS

% load global settings
clear all; tcm190910_config;


% settings
dataset = 'opto_mtc2mm';
maxEarlySessions = 3;  % only include this many days post lesion in lesions figs

% motor cortex muscimol
if strcmp(dataset, 'mtc_muscimol')
    colors = manipColors;
    vars.condition = struct('name', 'condition', 'levels', {{'saline', 'muscimol'}}, 'levelNames', {{'saline', 'muscimol'}});
    figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals
    
% motor cortex lesion
elseif strcmp(dataset, 'mtc_lesion')
    colors = manipColors;
    vars.condition = struct('name', 'condition', 'levels', {{'pre', 'post'}}, 'levelNames', {{'pre', 'post'}});
	figConditionals = struct('name', 'conditionNum', 'condition', @(x) x<=maxEarlySessions);

% barrel cortex lesion
elseif strcmp(dataset, 'senLesion')
    colors = manipColors;
    vars.condition = struct('name', 'condition', ...
        'levels', {{'pre', 'postBi', 'noWisk'}}, ...
        'levelNames', {{'pre', 'post', sprintf('no whiskers')}});
    figConditionals = struct('name', 'conditionNum', 'condition', @(x) x<=maxEarlySessions);

% opto (settings the same for all brain regions)
elseif contains(dataset, 'opto')
    colors = optoColors;
    vars.condition = struct('name', 'powerCondition', 'levels', [1 2 3 4], 'levelNames', {{'0', '1.4', '4.5', '15mW'}});
    figConditionals = struct('name', '', 'condition', @(x) x); % no conditionals
end


% load data
fprintf('loading data... '); 
load(fullfile(getenv('OBSDATADIR'), 'matlabData', [dataset '_data.mat']), 'data');
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'isWheelBreak', 'obsHgt', 'isPawSuccess', 'optoOnPositions', ...
    'isLeading', 'isFore', 'stepOverKinInterp', 'paw', 'controlStepKinInterp', 'preObsHgt', 'stepType', ...
    'isBigStep', 'preObsKin', 'conditionNum', 'isOptoOn', vars.condition.name});
if contains(dataset, {'lesion', 'Lesion'}); flat = flat([flat.conditionNum]<=maxEarlySessions); disp('restricting flat to early sessions...'); end
disp([dataset ' data loaded!'])

% define variables and conditionals to be used in bar plots
conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
figPadding = (length(vars.condition.levels)-2) * .25;

%% SUCCESS

figure('units', 'inches', 'position', [6.98 5.07 2.44+figPadding figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
dv = getDvMatrix(data, 'isTrialSuccess', vars.condition, {'mouse'}, figConditionals);
barFancy(dv, 'ylabel', 'success rate', 'levelNames', {vars.condition.levelNames}, 'colors', colors, 'YLim', [0 1], barProperties{:})
set(gca, 'position', [.2 .11 .78 .82])
print -clipboard -dmeta

%% SPEED VS POSITION

close all
yLims = [.2 .7];
figure('units', 'inches', 'position', [21.18 3.21 6.28 figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false); hold on

% add obstacle rectangle and lines
x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
line([0 0], yLims, 'linewidth', 2, 'color', axisColor)
line([x(1) x(1)], yLims, 'linewidth', 2, 'color', axisColor)
text(x(1), yLims(2), '\itobstacle on', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSize, 'FontName', font)
text(0, yLims(2), '\itobstacle reached', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSize, 'FontName', font)

if ~contains(dataset, 'opto')  % plots for non-opto experiments
    plotDvPsth(flat, 'velVsPosition', vars.condition.name, ...
        {'showLegend', false, 'conditionColors', colors, 'xlim', [-.5 .2], 'showErrorBars', true, ... 
         'plotConditions', vars.condition.levels, 'errorAlpha', .3, 'lineWidth', 3})

else  % plots for opto experiments
    
    % add light on, obs on, and obs at nose markers
    optoOnPos = nanmedian([flat.optoOnPositions]);
    x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
    line([optoOnPos optoOnPos], yLims, 'linewidth', 2, 'color', mean(colors(2:end,:),1))

    plotDvPsth(flat, 'velVsPosition', vars.condition.name, ...
        {'showLegend', false, 'conditionColors', colors, 'xlim', [-.9 .2], 'showErrorBars', false, ... 
         'plotConditions', vars.condition.levels, 'errorAlpha', .3, 'lineWidth', 3})
end

set(gca, 'YLim', yLims, 'YTick', linspace(yLims(1),yLims(2),3));
set(gca, 'YLim', yLims);
xlabel('distance (m)')
ylabel('velocity (m/s)')
set(gca, 'XColor', axisColor, 'YColor', axisColor, 'Color', 'black', 'FontSize', fontSize, 'FontName', font)
print -clipboard -dmeta

%% PAW HEIGHT

figure('units', 'inches', 'position', [6.98 5.07 2.44+figPadding figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
dv = getDvMatrix(data, 'preObsHgt', vars.condition, {'mouse'}, [figConditionals; conditionals.isFore])*1000;
barFancy(dv, 'ylabel', 'paw height (mm)', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:})
set(gca, 'position', [.2 .11 .78 .82])
print -clipboard -dmeta

%% PAW SHAPING

figure('units', 'inches', 'position', [6.98 5.07 2.44+figPadding figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, vars.condition, {'mouse'}, {'session'}, ...
    [figConditionals; conditionals.isFore; conditionals.isLeading], 'corr');
barFancy(dv, 'levelNames', {vars.condition.levelNames}, 'ylabel', 'paw-obstacle correlation', 'colors', colors, barProperties{:})
set(gca, 'position', [.2 .11 .78 .82])
print -clipboard -dmeta


%% SHAPING MOVING AVGS

close all;
flat_sub = flat(~isnan([flat.obsHgt]) & ...
                [flat.isLeading] & ...
                [flat.isFore]);
if ~contains(dataset, 'opto')
    [~, conditions] = ismember({flat_sub.(vars.condition.name)}, vars.condition.levels);  % convert condition name into condition number
else
    conditions = [flat_sub.isOptoOn]+1;
end
colorsTemp = manipColors;
figure('units', 'inches', 'position', [20.95 1.88 3.80 2.79], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
plot([0 100], [0 100], 'Color', axisColor*.5, 'LineWidth', 3) % add unity line
logPlotRick([flat_sub.obsHgt]*1000, [flat_sub.preObsHgt]*1000, ...
    {'colors', colorsTemp, 'conditions', conditions, 'xlabel', 'obstacle height', 'ylabel', 'paw height (mm)',  ...
    'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat_sub.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'lineWidth', 2})
set(gca, 'color', 'black', 'xcolor', axisColor, 'YColor', axisColor, 'FontSize', fontSize, 'FontName', font)
print -clipboard -dmeta

%% KINEMATICS

% settings
close all
obsHgtBins = 4; % discretize obstacle heights into this many bins
fading = .5; % within a plot, each paw's color fades from fading*color -> color
xLims = [-.058 0];
yLims = [0 .016];

% initializations
obsHgtDiscretized = num2cell(discretize([flat.obsHgt], linspace(3.4, 10, obsHgtBins+1)/1000));
[flat.obsHgtDiscretized] = obsHgtDiscretized{:};
flat_sub = flat(~isnan([flat.obsHgtDiscretized]) & ...
            ~[flat.isWheelBreak] & ...
            [flat.isLeading] & ...
            [flat.isFore]);
if ~contains(dataset, 'opto')
    [~, conditions] = ismember({flat_sub.(vars.condition.name)}, vars.condition.levels);  % convert condition name into condition number
else
    conditions = [flat_sub.isOptoOn]+1;
end
colorsTemp = manipColors;
figure('units', 'inches', 'position', [10.67 3.16 4.51 1.44*max(conditions)], ...
    'color', 'black', 'menubar', 'none', 'inverthardcopy', false); hold on
kinData = permute(cat(3, flat_sub.preObsKin), [3,1,2]);
kinDataCtl = permute(cat(3, flat_sub.controlStepKinInterp), [3,1,2]);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1); % change the x starting x position of ctl steps to match steps over

rows = max(conditions);
set(gca, 'visible', 'off')
for i = 1:rows
    ax = axes('position', [0 1-i*(1/rows) 1 1/rows ], 'color', 'black');
    bins = conditions==i;
    plotColor = repmat(colorsTemp(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    plotKinematics(kinData(bins,[1,3],:), [flat_sub(bins).obsHgt], [flat_sub(bins).obsHgtDiscretized], ...
        {'colors', plotColor, 'obsAlpha', 1, 'mouseNames', {flat_sub(bins).mouse}}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    set(gca, 'XLim', xLims, 'YLim', yLims, 'color', 'black')
end
print -clipboard -dmeta








