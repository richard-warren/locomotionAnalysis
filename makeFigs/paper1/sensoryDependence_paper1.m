% load experiment data
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'sensoryDependence_data.mat'), 'data'); disp('baseline data loaded!')

% global settings
colorWisk = [51 204 255]/255; %[255 204 51];
colorVision = [255 221 21]/255; %[51 204 255];
colorNone = [.2 .2 .2];
ctlColor = [.5 .5 .5];
varsToAvg = {'mouse'};
miceToExclude = {'sen1'};


% initializations
data.data = data.data(~ismember({data.data.mouse}, miceToExclude));
colors = [mean([colorWisk;colorVision],1); colorWisk; colorVision; colorNone];

vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'no light', 'light'}});
vars.sensoryCondition = struct('name', 'sensoryCondition', 'levels', {{'WL', 'W', 'L', '-'}}, 'levelNames', {{'W+V', 'W', 'V', '-'}});
vars.whiskers = struct('name', 'whiskers', 'levels', {{'full', 'none'}}, 'levelNames', {{'whiskers', 'no whiskers'}});
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore paw', 'hind paw'}});

conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);




%% ----------
% PLOT THINGS
%  ----------


%% BARS

% initializations
rows = 2;
cols = 4;
figure('name', 'sensoryDependence', 'color', 'white', 'menubar', 'none', 'position', [2000 50 1200 400])

% success
subplot(rows, cols, 1);
conditions = [vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'isTrialSuccess', conditions, varsToAvg);
barPlotRick(dvMatrix, {'conditionNames', {conditions.levelNames}, 'ylabel', 'success rate', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colors, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [0 1], 'ytick', 0:.5:1, ...
    'compareToFirstOnly', true})

% velocity
subplot(rows, cols, 2);
conditions = [vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'velAtWiskContact', conditions, varsToAvg);
barPlotRick(dvMatrix, {'conditionNames', {conditions.levelNames}, 'ylabel', 'velovity at whisker contact (m/s)', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colors, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [0 .75], 'ytick', 0:.25:.75, ...
    'compareToFirstOnly', false})

% paw error rate
subplot(rows, cols, 3:4);
conditions = [vars.isFore; vars.isLeading; vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'isPawSuccess', conditions, varsToAvg);
barPlotRick(dvMatrix, {'conditionNames', {conditions.levelNames}, 'ylabel', 'paw error rate', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', repmat(colors,4,1), 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [0 1], 'ytick', 0:.5:1, ...
    'compareToFirstOnly', false, 'lineWidth', 1.5})

% step over height for all paws
subplot(rows, cols, 5:6);
conditions = [vars.isFore; vars.isLeading; vars.sensoryCondition];
dvMatrix = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg) * 1000;
barPlotRick(dvMatrix, {'conditionNames', {conditions.levelNames}, 'ylabel', 'step over height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', repmat(colors,4,1), 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, ...
    'compareToFirstOnly', false})

% step over height for leading forelimbs
subplot(rows, cols, 7)
conditions = [vars.sensoryCondition];
figConditionals = [conditionals.isFore; conditionals.isLeading];
matMod = getDvMatrix(data, 'preObsHgt', conditions, varsToAvg, figConditionals) * 1000;
matBl = getDvMatrix(data, 'baselineStepHgt', conditions, varsToAvg, figConditionals) * 1000;
dvMatrix = permute(cat(3,matBl,matMod), [1 3 2]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(ctlColor,8,1);
colorsWithBl(2:2:8,:) = colors;
barPlotRick(dvMatrix, {'conditionNames', {conditions.levelNames}, 'ylabel', {'leading fore paw', 'step over height (mm)'}, ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colorsWithBl, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, ...
    'compareToFirstOnly', false, 'lineWidth', 1.5})

% wisk contact position (light, manip)
subplot(rows, cols, 8);
conditions = [vars.sensoryCondition];
dvMatrix = abs(getDvMatrix(data, 'wiskContactPosition', conditions, varsToAvg)) * 1000;
barPlotRick(dvMatrix, {'conditionNames', {conditions.levelNames}, 'ylabel', {'distance to nose', 'at whisker contact (mm)'}, ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colors, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, ...
    'compareToFirstOnly', false})

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'sensoryDependenceBars');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');



%% SPEED VS. POSITION

% settings
obsOnColor = [0 0 0];
obsOnAlpha = .05;
yLims = [.3 .7];
plotSequence = [4 3 2 1]; % determine which lines on plotted on top of which

% initializations
flat = flattenData(data, {'mouse', 'session', 'trial', 'sensoryCondition', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'velContinuousAtContact', 'velContinuousAtContactX', 'isWheelBreak', 'wiskContactPosition'});
colorsTemp = [colors(1:end-1,:); .6 .6 .6]; % the no vision no whisker condition can be a little lighter here

% speed vs position
figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 700 400], 'inverthardcopy', 'off')
rectangle('Position', [x(1) yLims(1) diff(x) diff(yLims)], ...
    'FaceColor', [obsOnColor obsOnAlpha], 'EdgeColor', 'none');
line([0 0], yLims, 'linewidth', 2, 'color', obsOnColor)
plotDvPsth(flat, 'velVsPosition', 'sensoryCondition', ...
    {'showLegend', false, 'conditionColors', colorsTemp(plotSequence,:), 'xlim', [-.5 .2], ... 
     'plotConditions', vars.sensoryCondition.levels(plotSequence), 'errorAlpha', .1, 'lineWidth', 4})
x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
set(gca, 'YLim', yLims, 'YTick', linspace(yLims(1),yLims(2),3));
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'sensoryDependenceVel');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% HEIGHT SHAPING


% SCATTER ACROSS ALL ANIMALS

% settings
xLims = [3 10];
yLims = [3 20];
isHgtPreObs = true; % do you measure the height at peak (false) or before the paw reaches the obstacle (true)

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isLeading', 'stepOverMaxHgt', 'sensoryCondition'});
flat = flat(~isnan([flat.obsHgt]) & ...
            ~isnan([flat.stepOverMaxHgt]) & ...
            ~isnan([flat.preObsHgt]) & ...
            [flat.isLeading] & ...
            [flat.isFore]); % add conditionals here

if isHgtPreObs; pawHgts = [flat.preObsHgt]*1000; else; pawHgts = [flat.stepOverMaxHgt]*1000; end
obsHgts = [flat.obsHgt]*1000;
[~, conditions] = ismember({flat.sensoryCondition}, vars.sensoryCondition.levels);

figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 50 500 400])
scatterPlotRick(obsHgts, pawHgts, conditions, ...
    {'colors', colors, 'maxScatterPoints', 5000, 'lineAlpha', 1, 'scatAlpha', .1, 'scatSize', 40});
set(gca, 'XLim', xLims, 'YLim', yLims)

% add unity line
plot([0 xLims(2)], [0 xLims(2)], 'Color', [1 1 1]*.6, 'LineWidth', 3) % add unity line
xlabel('obstacle height (mm)')
ylabel('paw height (mm)')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'sensoryDependenceHeightShapingScatters');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');



% BINNED BY ANIMAL, SPEED

% settings
binNum = 40;
scatAlpha = .1;
scatSize = 50;

% initializations
mice = unique({flat.mouse});
binEdges = linspace(xLims(1), xLims(2), binNum+1);
binCenters = binEdges(1:end-1) + .5*diff(binEdges(1:2));
binIds = discretize(obsHgts, binEdges);

% collect data
binnedHgts = nan(4, length(mice), binNum); % contains the median paw height for each conditino, mouse, and paw height
for c = 1:4
    for m = 1:length(mice)
        for b = 1:binNum
            binnedHgts(c,m,b) = nanmean(pawHgts((conditions==c) & strcmp({flat.mouse}, mice{m}) & binIds==b));
        end
    end
end

% plot condition data
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 50 500 400])
scatterPlotRick(repelem(binCenters,length(mice)*4), binnedHgts(:), repmat(1:4,1,binNum*length(mice)), ...
    {'colors', colors, 'maxScatterPoints', 5000, 'lineAlpha', 1, 'scatAlpha', .2, 'scatSize', 40});
set(gca, 'XLim', xLims)

% add unity line
plot([0 xLims(2)], [0 xLims(2)], 'Color', [1 1 1]*.6, 'LineWidth', 3) % add unity line
xlabel('obstacle height (mm)')
ylabel('paw height (mm)')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'sensoryDependenceHeightShapingScatters_perMouse');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% HEIGHT SHAPING BARS
[corrs, slopes] = deal(nan(4,length(mice))); % isFore(10) X isLeading(10) X mouse
foreSequence = [true false];
leadingSequence = [true false];

for i = 1:4
    for k = 1:length(mice)
        bins = conditions==i & ...
               strcmp({flat.mouse}, mice{k});
        x = obsHgts(bins);
        y = pawHgts(bins);
        fit = polyfit(x, y, 1);
        corrs(i,k) = corr(x', y');
        slopes(i,k) = fit(1);
    end
end


figure('position', [2000 400 1000 400], 'color', 'white', 'menubar', 'none');

% correlations
subplot(1,2,1)
barPlotRick(corrs, {'conditionNames', {vars.sensoryCondition.levelNames}, 'ylabel', 'paw/obstacle height correlation', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colors, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [-.2 .6], 'ytick', -.2:.4:.6, ...
    'compareToFirstOnly', true})

% slopes
subplot(1,2,2)
barPlotRick(slopes, {'conditionNames', {vars.sensoryCondition.levelNames}, 'ylabel', 'paw/obstacle height slope', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colors, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [-.2 .8], 'ytick', -.2:.5:.8, ...
    'compareToFirstOnly', true})

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'sensoryDependenceHeightShapingBars');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

%% KINEMATICS

% settings
obsHgtBins = 4; % discretize obstacle heights into this many bins
fading = .25; % within a plot, each paw's color fades from fading*color -> color
xLims = [-.05 0];
yLims = [0 .016];

% initializations
flat = flattenData(data, {'mouse', 'session', 'trial', 'sensoryCondition', 'isWheelBreak', ...
    'obsHgt', 'isFore', 'isLeading', 'preObsKin', 'controlStepKinInterp'});
obsHgtDiscretized = num2cell(discretize([flat.obsHgt], linspace(3.4, 10, obsHgtBins+1)/1000));
[flat.obsHgtDiscretized] = obsHgtDiscretized{:};
flat = flat(~isnan([flat.obsHgtDiscretized]) & ...
            ~[flat.isWheelBreak] & ...
            [flat.isLeading] & ...
            [flat.isFore]);
[~, conditions] = ismember({flat.sensoryCondition}, vars.sensoryCondition.levels);
kinData = permute(cat(3, flat.preObsKin), [3,1,2]);
kinDataCtl = permute(cat(3, flat.controlStepKinInterp), [3,1,2]);
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1); % change the x starting x position of ctl steps to match steps over

%%
close all;
figure('position', [2000 200 400 800], 'color', 'white', 'menubar', 'none'); hold on;
colorsTemp = [colors(1:end-1,:); .8 .8 .8];
for i = 1:4
    subplot(5,1,i)
    bins = conditions==i;
    plotColor = repmat(colorsTemp(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    
    plotKinematics(kinDataCtl(bins,[1,3],:), [flat(bins).obsHgt], ones(1,sum(bins)), ...
        {'colors', ctlColor, 'obsAlpha', 0, 'lineAlpha', .8}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    plotKinematics(kinData(bins,[1,3],:), [flat(bins).obsHgt], [flat(bins).obsHgtDiscretized], ...
        {'colors', plotColor, 'obsAlpha', 1, 'lineAlpha', .8, 'mouseNames', {flat(bins).mouse}}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    set(gca, 'XLim', xLims, 'YLim', yLims)
end

subplot(5,1,5)
plotKinematics(kinData(:,[1,3],:), [flat.obsHgt], conditions, ...
    {'colors', colors, 'obsAlpha', .25, 'lineAlpha', 1, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @nanstd, 'lineWidth', 3}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
set(gca, 'XLim', xLims, 'YLim', yLims)

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'sensoryDependenceKinematics_obsHgt');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');






