%% GLOBAL SETTINGS

% first load data // each dataset has a corresponding flattened data set //
% incidual plots may take subsets of the flattened dataset

axisColor = [.95 .95 .95];
pawColors = hsv(4);
mouseColors = 'jet';
ctlStepColor = [1 1 1]*.9;

%% INITIALIZATIONS

% load data
fprintf('loading data... '); 
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data');
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'isWheelBreak', 'wiskContactPosition', 'obsHgt', 'isPawSuccess', ...
    'isLeading', 'isFore', 'stepOverKinInterp', 'paw', 'controlStepKinInterp', 'preObsHgt', 'stepType'});
disp('data loaded!')

% conditions and conditionals
conditionNames = {{'fore', 'hind'}, {'leading', 'trailing'}};
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore paw', 'hind paw'}});

conditionals.lightOn = struct('name', 'isLightOn', 'condition', @(x) x==1);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
conditionals.none = struct('name', '', 'condition', @(x) x); % no conditionals


%% SPEED VS POSITION

% initializations
figure('units', 'inches', 'position', [3.83 3.61 6.67 3.09], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
yLims = [0 .8];

% add obstacle rectangle and lines
x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
line([0 0], yLims, 'linewidth', 2, 'color', axisColor)
line([x(1) x(1)], yLims, 'linewidth', 2, 'color', axisColor)
text(x(1), yLims(2), '\itobstacle on', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
text(0, yLims(2), '\itobstacle reached', 'color', axisColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

% plot
plotDvPsth(flat, 'velVsPosition', [], ...
    {'plotMouseAvgs', true, 'showLegend', false, 'conditionColors', [1 1 1], 'mouseColors', mouseColors, ...
    'errorFcn', @(x) nanstd(x), 'xlim', [-.5 .2], 'mouseAlpha', .5, 'errorAlpha', .4})

% pimp fig
set(gca, 'YLim', [0 .8]);
xlabel('distance (m)')
ylabel('velocity (m/s)')
set(gca, 'XColor', axisColor, 'YColor', axisColor, 'Color', 'black')
print -clipboard -dmeta

%% BASELINE SUCCESS

figure('units', 'inches', 'position', [3.83 3.61 3.92 3.09], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)

logPlotRick([flat.obsHgt]*1000, [flat.isPawSuccess], ...
    {'colors', pawColors, 'conditions', [flat.stepType], 'xlabel', 'obstacle height', 'ylabel', 'success rate', 'plotMice', false, ...
    'xlim', [3.4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'computeVariance', false, 'ylim', [0 1], 'plotMouseErrors', false, 'lineWidth', 3})
set(gca, 'xlim', [3.4 10], 'color', 'black', 'xcolor', axisColor, 'YColor', axisColor)

print -clipboard -dmeta

%% BASELINE KINEMATICS

% settings
obsHgtBins = 1; % discretize obstacle heights into this many bins

figure('units', 'inches', 'position', [4.27 6.41 10.71 1.99], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
colNames = {'hind', 'fore'};
rowNames = {'trailing', 'leading'};
conditionSequence = [4 2 3 1]; % mapping between plot index and condition sequence, which is specified above
xLims = [-.05 .05];
yLims = [0 .016];
fading = .5; % within a plot, each paw's color fades from fading*color -> color
obsHgtDiscretized = num2cell(discretize([flat.obsHgt], linspace(3.4, 10, obsHgtBins+1)/1000));
[flat.obsHgtDiscretized] = obsHgtDiscretized{:};
flat = flat(~isnan([flat.obsHgtDiscretized]) & ...
            ~[flat.isWheelBreak] & ...
            [flat.isLightOn]);


if obsHgtBins==1; lineWidth = 4; else lineWidth = 3; end

% get kin data
kinData = permute(cat(3, flat.stepOverKinInterp), [3,1,2]);
kinData = cat(1, kinData, permute(cat(3, flat.controlStepKinInterp), [3,1,2])); % append with control steps temporarily

% flip y values s.t. leading is always right and lagging always left
bins = [flat.paw]==1 & [flat.isLeading]; kinData(bins,2,:) = -kinData(bins,2,:);
bins = [flat.paw]==2 & [flat.isLeading]; kinData(bins,2,:) = -kinData(bins,2,:);
bins = [flat.paw]==3 & ~[flat.isLeading]; kinData(bins,2,:) = -kinData(bins,2,:);
bins = [flat.paw]==4 & ~[flat.isLeading]; kinData(bins,2,:) = -kinData(bins,2,:);

% split real and control steps
kinDataCtl = kinData(length(flat)+1:end,:,:);
kinData = kinData(1:length(flat),:,:);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1);  % change the x starting x position of ctl steps to match steps over

plotInd = 1;
buf = .03;
wid = (1-3*buf)/2;
hgt = (1-buf)/2;
positions = [buf hgt; wid+2*buf hgt; buf 0; wid+2*buf 0];
for i = conditionSequence
    subplot(2,2,plotInd)
    plotColor = repmat(pawColors(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    bins = [flat.stepType]==i;

    % plot control step
    plotKinematics(kinDataCtl(bins,[1,3],:), [flat(bins).obsHgt], ones(1,sum(bins)), ...
        {'colors', ctlStepColor, 'obsAlpha', 0, 'lineAlpha', .8, 'lineWidth', lineWidth}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    
    plotKinematics(kinData(bins,[1,3],:), [flat(bins).obsHgt], [flat(bins).obsHgtDiscretized], ...
        {'colors', plotColor, 'obsAlpha', 1, 'lineAlpha', .8, 'mouseNames', {flat(bins).mouse}, 'lineColor', axisColor, 'lineWidth', lineWidth}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    set(gca, 'XLim', xLims, 'YLim', yLims, 'color', 'black')
    set(gca, 'position', [positions(plotInd,:) wid hgt])
    
    if plotInd==1
        text(0, yLims(2), colNames{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'color', axisColor);
        text(xLims(1), mean(yLims), rowNames{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90, 'color', axisColor);
    elseif plotInd==2
        text(0, yLims(2), colNames{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'color', axisColor);
    elseif plotInd==3
        text(xLims(1), mean(yLims), rowNames{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90, 'color', axisColor);
    end

    plotInd = plotInd+1;
end

print -clipboard -dmeta


%% PRE OBS HEIGHT

figure('units', 'inches', 'position', [4.46 1.36 3.8 2.8], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
figVars = [vars.isFore; vars.isLeading];
figConditionals = [conditionals.none];

matMod = getDvMatrix(data, 'preObsHgt', figVars, {'mouse'}, figConditionals) * 1000;
matBl = getDvMatrix(data, 'controlPreObsHgt', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matBl,matMod), [1 2 4 3]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(ctlStepColor,8,1);
colorsWithBl(2:2:8,:) = pawColors;
barFancy(mat, 'axisColor', axisColor, 'levelNames', {figVars.levelNames}, 'ylabel', 'paw height (mm)', ...
    'colors', colorsWithBl, 'scatterAlpha', .8, 'barAlpha', .4, 'labelSizePerFactor', .1, 'lineThickness', 2)
print -clipboard -dmeta

%% STEP OVER LENGTH

figure('units', 'inches', 'position', [4.46 1.36 3.8 2.8], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
figVars = [vars.isFore; vars.isLeading];
figConditionals = [conditionals.none];

matMod = getDvMatrix(data, 'stepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
matBl = getDvMatrix(data, 'controlStepLength', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matBl,matMod), [1 2 4 3]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(ctlStepColor,8,1);
colorsWithBl(2:2:8,:) = pawColors;
barFancy(mat, 'axisColor', axisColor, 'levelNames', {figVars.levelNames}, 'ylabel', 'step length (mm)', ...
    'colors', colorsWithBl, 'scatterAlpha', .8, 'barAlpha', .4, 'labelSizePerFactor', .1, 'lineThickness', 2)

print -clipboard -dmeta

%% SHOW LEADING LAGGING SCHEMATIC

showLeadingLaggingImg('190318_000', 43, pawColors)
print -clipboard -dmeta

%% SHAPING CORRELATIONS
figure('units', 'inches', 'position', [13.22 3.65 3.19 3.45], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)

figVars = [vars.isFore; vars.isLeading];
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, figVars, {'mouse'}, {'session'}, [conditionals.isLe], 'corr');
barFancy(dv, 'axisColor', axisColor, 'levelNames', {figVars.levelNames}, 'ylabel', 'step length (mm)', ...
    'colors', pawColors, 'scatterAlpha', .8, 'barAlpha', .4, 'labelSizePerFactor', .15, 'lineThickness', 2, 'barWidth', .75)
set(gca, 'position', [.15 .11 .77 .81])
print -clipboard -dmeta

%% SHAPING MOVING AVGS

figure('units', 'inches', 'position', [4.46 1.36 3.8 2.8], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
plot([0 xLims(2)], [0 xLims(2)], 'Color', [1 1 1]*.6, 'LineWidth', 3) % add unity line
logPlotRick([flat.obsHgt]*1000, [flat.preObsHgt]*1000, ...
    {'colors', pawColors, 'conditions', [flat.stepType], 'xlabel', 'obstacle height', 'ylabel', 'paw height (mm)',  ...
    'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'lineWidth', 2})
set(gca, 'color', 'black', 'xcolor', axisColor, 'YColor', axisColor)
print -clipboard -dmeta

