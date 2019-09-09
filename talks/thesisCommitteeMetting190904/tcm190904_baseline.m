%% INITIALIZATIONS

% load global settings
clear all; tcm190910_config;

% load data
fprintf('loading data... '); 
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data');
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'isWheelBreak', 'wiskContactPosition', 'obsHgt', 'isPawSuccess', ...
    'isLeading', 'isFore', 'stepOverKinInterp', 'paw', 'controlStepKinInterp', 'preObsHgt', 'stepType', ...
    'isBigStep'});
disp('data loaded!')

% conditions and conditionals
conditionNames = {{'fore', 'hind'}, {'leading', 'trailing'}};
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore paw', 'hind paw'}});

conditionals.lightOn = struct('name', 'isLightOn', 'condition', @(x) x==1);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
conditionals.none = struct('name', '', 'condition', @(x) x); % no conditionals

disp('all done!')


%% SPEED VS POSITION

% initializations
figure('units', 'inches', 'position', [3.83 3.61 7.5 figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false); hold on
yLims = [0 .8];
flat_sub = flat(~[flat.isWheelBreak]);

% add obstacle rectangle and lines
x = [nanmean([flat_sub.obsOnPositions]) nanmean([flat_sub.obsOffPositions])];
line([0 0], yLims, 'linewidth', 2, 'color', axisColor)
line([x(1) x(1)], yLims, 'linewidth', 2, 'color', axisColor)
text(x(1), yLims(2), '\itobstacle on', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSize, 'FontName', font)
text(0, yLims(2), '\itobstacle reached', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fontSize, 'FontName', font)

% plot
plotDvPsth(flat_sub, 'velVsPosition', [], ...
    {'plotMouseAvgs', true, 'showLegend', false, 'conditionColors', [1 1 1], 'mouseColors', mouseColors, ...
    'errorFcn', @(x) nanstd(x), 'xlim', [-.5 .2], 'mouseAlpha', .5, 'errorAlpha', .4})

% pimp fig
set(gca, 'YLim', [0 .8]);
xlabel('distance (m)')
ylabel('velocity (m/s)')
set(gca, 'XColor', axisColor, 'YColor', axisColor, 'Color', 'black', 'FontSize', fontSize, 'FontName', font)
print -clipboard -dmeta

%% BASELINE SUCCESS

figure('units', 'inches', 'position', [3.83 3.61 3.92 figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)

logPlotRick([flat.obsHgt]*1000, [flat.isPawSuccess], ...
    {'colors', pawColors, 'conditions', [flat.stepType], 'xlabel', 'obstacle height (mm)', 'ylabel', 'success rate', 'plotMice', false, ...
    'xlim', [3.4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'computeVariance', false, 'ylim', [0 1], 'plotMouseErrors', false, 'lineWidth', 3})
set(gca, 'xlim', [3.4 10], 'color', 'black', 'xcolor', axisColor, 'YColor', axisColor, 'FontSize', fontSize, 'FontName', font)

print -clipboard -dmeta

%% BASELINE KINEMATICS

% settings
obsHgtBins = 4; % discretize obstacle heights into this many bins

close all
figure('units', 'inches', 'position', [4.27 6.41 13.33 2.8], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
colNames = {'hind', 'fore'};
rowNames = {'trailing', 'leading'};
conditionSequence = [4 2 3 1]; % mapping between plot index and condition sequence, which is specified above
xLims = [-.05 .05];
yLims = [0 .016];
fading = .3; % within a plot, each paw's color fades from fading*color -> color
obsHgtDiscretized = num2cell(discretize([flat.obsHgt], linspace(3.4, 10, obsHgtBins+1)/1000));
[flat.obsHgtDiscretized] = obsHgtDiscretized{:};
flat_sub = flat(~isnan([flat.obsHgtDiscretized]) & ...
                ~[flat.isWheelBreak]);


if obsHgtBins==1; lineWidth = 5; else lineWidth = 3; end

% get kin data
kinData = permute(cat(3, flat_sub.stepOverKinInterp), [3,1,2]);
kinData = cat(1, kinData, permute(cat(3, flat_sub.controlStepKinInterp), [3,1,2])); % append with control steps temporarily

% flip y values s.t. leading is always right and lagging always left
bins = [flat_sub.paw]==1 & [flat_sub.isLeading]; kinData(bins,2,:) = -kinData(bins,2,:);
bins = [flat_sub.paw]==2 & [flat_sub.isLeading]; kinData(bins,2,:) = -kinData(bins,2,:);
bins = [flat_sub.paw]==3 & ~[flat_sub.isLeading]; kinData(bins,2,:) = -kinData(bins,2,:);
bins = [flat_sub.paw]==4 & ~[flat_sub.isLeading]; kinData(bins,2,:) = -kinData(bins,2,:);

% split real and control steps
kinDataCtl = kinData(length(flat_sub)+1:end,:,:);
kinData = kinData(1:length(flat_sub),:,:);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1);  % change the x starting x position of ctl steps to match steps over

plotInd = 1;
buf = .03;
wid = (1-3*buf)/2;
hgt = (1-buf)/2;
positions = [buf hgt; wid+2*buf hgt; buf 0; wid+2*buf 0];
for i = conditionSequence
    subplot(2,2,plotInd)
    plotColor = repmat(pawColors(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    bins = [flat_sub.stepType]==i;

    % plot control step
    plotKinematics(kinDataCtl(bins,[1,3],:), [flat_sub(bins).obsHgt], ones(1,sum(bins)), ...
        {'colors', axisColor, 'obsAlpha', 0, 'lineAlpha', .8, 'lineWidth', lineWidth*.5}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    
    % plot step over obstacle
    plotKinematics(kinData(bins,[1,3],:), [flat_sub(bins).obsHgt], [flat_sub(bins).obsHgtDiscretized], ...
        {'colors', plotColor, 'obsAlpha', 1, 'mouseNames', {flat_sub(bins).mouse}, 'lineColor', axisColor, 'lineWidth', lineWidth}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    set(gca, 'XLim', xLims, 'YLim', yLims, 'color', 'black')
    set(gca, 'position', [positions(plotInd,:) wid hgt])
    
    % add column, row labels
    if plotInd==1
        text(0, yLims(2), colNames{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'color', axisColor, ...
            'FontSize', fontSize, 'FontName', font);
        text(xLims(1), mean(yLims), rowNames{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'Rotation', 90, 'color', axisColor, 'FontSize', fontSize, 'FontName', font);
    elseif plotInd==2
        text(0, yLims(2), colNames{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'color', axisColor, ...
            'FontSize', fontSize, 'FontName', font);
    elseif plotInd==3
        text(xLims(1), mean(yLims), rowNames{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'Rotation', 90, 'color', axisColor, 'FontSize', fontSize, 'FontName', font);
    end

    plotInd = plotInd+1;
end

print -clipboard -dmeta


%% PRE OBS HEIGHT

figure('units', 'inches', 'position', [4.46 1.36 4.8 figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
figVars = [vars.isFore; vars.isLeading];
figConditionals = [conditionals.none];

matMod = getDvMatrix(data, 'preObsHgt', figVars, {'mouse'}, figConditionals) * 1000;
matBl = getDvMatrix(data, 'controlPreObsHgt', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matBl,matMod), [1 2 4 3]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(axisColor,8,1);
colorsWithBl(2:2:8,:) = pawColors;
barFancy(mat, 'levelNames', {figVars.levelNames}, 'ylabel', 'paw height (mm)', 'colors', colorsWithBl, barProperties{:})
set(gca, 'FontSize', fontSize, 'FontName', font)
print -clipboard -dmeta

%% STEP OVER LENGTH

figure('units', 'inches', 'position', [4.46 1.36 4.8 figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
figVars = [vars.isFore; vars.isLeading];
figConditionals = [conditionals.none];

matMod = getDvMatrix(data, 'stepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
matBl = getDvMatrix(data, 'controlStepLength', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matBl,matMod), [1 2 4 3]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(axisColor,8,1);
colorsWithBl(2:2:8,:) = pawColors;
barFancy(mat, 'levelNames', {figVars.levelNames}, 'ylabel', 'step length (mm)', 'colors', colorsWithBl, barProperties{:})
set(gca, 'FontSize', fontSize, 'FontName', font)
print -clipboard -dmeta

%% LEADING LAGGING SCHEMATIC

showLeadingLaggingImg('190318_000', 43, pawColors)
print -clipboard -dmeta

%% SHAPING MOVING AVGS

close all
flat_sub = flat(~isnan([flat.obsHgt]) & ...
                ~isnan([flat.preObsHgt])); % add conditionals here
figure('units', 'inches', 'position', [4.46 1.36 4.25 figHgt], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
plot([0 100], [0 100], 'Color', axisColor, 'LineWidth', 2); hold on % add unity line
logPlotRick([flat_sub.obsHgt]*1000, [flat_sub.preObsHgt]*1000, ...
    {'colors', pawColors, 'conditions', [flat_sub.stepType], 'xlabel', 'obstacle height (mm)', 'ylabel', 'paw height (mm)',  ...
    'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat_sub.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'lineWidth', 3, 'errorAlpha', .4})
set(gca, 'color', 'black', 'xcolor', axisColor, 'YColor', axisColor, 'FontSize', fontSize, 'FontName', font)
print -clipboard -dmeta









