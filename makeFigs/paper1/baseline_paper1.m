%% load experiment data
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

% global settings
colors = hsv(4);  % last term is control step color
ctlColor = [.5 .5 .5];
isLeading = [true false true false]; % sequence of conditions for plots
isFore = [true true false false];
conditionNames = {{'fore paw', 'hind paw'}, {'leading', 'trailing'}};


%% ----------
% PLOT THINGS
%  ----------

%% IMAGE MONTAGE OF LEADING/LAGGING/HIND/FORE

showLeadingLaggingImg('190318_000', 43, colors)

%% SPEED VS. POSITION

% settings
yLims = [0 .8];
obsOnColor = [0 0 0];
obsOnAlpha = .05;
meanColor = [0 0 0];

% initializations
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'isWheelBreak', 'wiskContactPosition'});
flat = flat([flat.isLightOn]);
close all; figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 300], 'inverthardcopy', 'off')

% add obstacle rectangle and lines
x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
rectangle('Position', [x(1) yLims(1) diff(x) diff(yLims)], ...
    'FaceColor', [obsOnColor obsOnAlpha], 'EdgeColor', 'none');
line([0 0], yLims, 'linewidth', 2, 'color', obsOnColor)
    
% plot
plotDvPsth(flat, 'velVsPosition', [-.5 .2], [], [], ...
    {'plotMouseAvgs', true, 'showLegend', false, 'conditionColors', [0 0 0], 'errorFcn', @(x) nanstd(x)})

% pimp fig
set(gca, 'YLim', yLims);
xlabel('position relaive to nose (m)')
ylabel('velocity (m/s)')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineVel');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% KINEMATICS

% settings
colNames = {'hind paw', 'fore paw'};
rowNames = {'trailing', 'leading'};
conditionSequence = [4 2 3 1]; % mapping between plot index and condition sequence, which is specified above

obsHgtBins = 4; % discretize obstacle heights into this many bins
xLims = [-.05 .05];
yLims = [0 .016];
fading = .25; % within a plot, each paw's color fades from fading*color -> color

% HEIGHT BINNED KINEMATICS (leading/lagging, fore/hind)
flat = flattenData(data, {'mouse', 'session', 'isTrialSuccess', 'trial', 'isLightOn', 'isWheelBreak', ...
    'obsHgt', 'isFore', 'isLeading', 'stepOverKinInterp', 'paw', 'controlStepKinInterp'});

obsHgtDiscretized = num2cell(discretize([flat.obsHgt], linspace(3.4, 10, obsHgtBins+1)/1000));
[flat.obsHgtDiscretized] = obsHgtDiscretized{:};
flat = flat(~isnan([flat.obsHgtDiscretized]) & ...
            ~[flat.isWheelBreak] & ...
            [flat.isLightOn]);


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

% change the x starting x position of ctl steps to match steps over
kinDataCtl(:,1,:) = kinDataCtl(:,1,:) - kinDataCtl(:,1,1) + kinData(:,1,1);


% initializations
conditions = [[flat.obsHgtDiscretized] ones(1,length(flat))*(obsHgtBins+1)]; % add an extra condition for the control kinematics!
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 700 1600 250])

plotInd = 1;
for i = conditionSequence
    subplot(2,2,plotInd)
    plotColor = repmat(colors(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    bins = [flat.isLeading]==isLeading(i) & ...
           [flat.isFore]==isFore(i);

    % plot control step
    plotKinematics(kinDataCtl(bins,[1,3],:), [flat(bins).obsHgt], ones(1,sum(bins)), ...
        {'colors', ctlColor, 'obsAlpha', 0, 'lineAlpha', .8}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    
    plotKinematics(kinData(bins,[1,3],:), [flat(bins).obsHgt], conditions(bins), ...
        {'colors', plotColor, 'obsAlpha', 1, 'lineAlpha', .8, 'mouseNames', {flat(bins).mouse}}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    set(gca, 'XLim', xLims, 'YLim', yLims)
    plotInd = plotInd+1;
end

% add row, col labels
subplot(2,2,1);
text(0, yLims(2), colNames{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(xLims(1), mean(yLims), rowNames{1}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90);
subplot(2,2,2);
text(0, yLims(2), colNames{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subplot(2,2,3);
text(xLims(1), mean(yLims), rowNames{2}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Rotation', 90);

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineKinematics_obsHgt');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% TOP VIEW OVERLAYS

% initializations
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 500 900 150])

% get condition numbers, where each condition is unique combinate of isLeading and isFore
a = [isLeading; isFore]'; % each row is a condition (00,01,10,11)
b = [[flat.isLeading]; [flat.isFore]]';
[~, stepTypeConditions] = ismember(b, a, 'rows');

plotKinematics(kinData(:,[1,3],:), [flat.obsHgt], stepTypeConditions', ...
    {'colors', colors, 'obsAlpha', .25, 'lineAlpha', 1, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @nanstd, 'lineWidth', 3}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
set(gca, 'XLim', xLims)

%save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineKinematics_overlayTop');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

% BOT VIEW OVERLAYS
yLimsBot = [-1 1]*.02;

figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 100 900 300])
plotKinematics(kinData(:,[1,2],:), [flat.obsHgt], stepTypeConditions', ...
    {'colors', colors, 'obsAlpha', .25, 'lineAlpha', 1, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @nanstd, 'lineWidth', 3, 'isBotView', true}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
set(gca, 'XLim', xLims, 'YLim', yLimsBot)

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineKinematics_overlayBot');
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
    'obsHgt', 'preObsHgt', 'isFore', 'isLeading', 'stepOverMaxHgt'});
flat = flat([flat.isLightOn] & ...
            ~isnan([flat.obsHgt]) & ...
            ~isnan([flat.stepOverMaxHgt]) & ...
            ~isnan([flat.preObsHgt])); % add conditionals here

% get condition numbers, where each condition is unique combination of isLeading and isFore
a = [isLeading; isFore]'; % each row is a condition (00,01,10,11)
b = [[flat.isLeading]; [flat.isFore]]';
[~, stepTypeConditions] = ismember(b, a, 'rows');

if isHgtPreObs; pawHgts = [flat.preObsHgt]*1000; else; pawHgts = [flat.stepOverMaxHgt]*1000; end
obsHgts = [flat.obsHgt]*1000;
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 50 500 400])
scatterPlotRick(obsHgts, pawHgts, stepTypeConditions, ...
    {'colors', colors, 'maxScatterPoints', 5000, 'lineAlpha', 1, 'scatAlpha', .1, 'scatSize', 40});
set(gca, 'XLim', xLims, 'YLim', yLims)

% add unity line
plot([0 xLims(2)], [0 xLims(2)], 'Color', [1 1 1]*.6, 'LineWidth', 3) % add unity line
xlabel('obstacle height (mm)')
ylabel('paw height (mm)')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineHeightShapingScatters');
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
            binnedHgts(c,m,b) = nanmean(pawHgts((stepTypeConditions==c)' & strcmp({flat.mouse}, mice{m}) & binIds==b));
        end
    end
end


% plot condition data
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 50 500 400])
scatterPlotRick(repelem(binCenters,length(mice)*4), binnedHgts(:), repmat(1:4,1,binNum*length(mice)), ...
    {'colors', colors, 'maxScatterPoints', 5000, 'lineAlpha', 1, 'scatAlpha', .1, 'scatSize', 40});
set(gca, 'XLim', xLims)

% add unity line
plot([0 xLims(2)], [0 xLims(2)], 'Color', [1 1 1]*.6, 'LineWidth', 3) % add unity line
xlabel('obstacle height (mm)')
ylabel('paw height (mm)')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineHeightShapingScatters_perMouse');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% HEIGHT SHAPING BARS
[corrs, slopes] = deal(nan(2,2,length(mice))); % isFore(10) X isLeading(10) X mouse
foreSequence = [true false];
leadingSequence = [true false];

for i = 1:2 % isFore
    for j = 1:2 % isLeading
        for k = 1:length(mice)
            
            bins = [flat.isFore]==foreSequence(i) & ...
                   [flat.isLeading]==leadingSequence(j) & ...
                   strcmp({flat.mouse}, mice{k});
            x = obsHgts(bins);
            y = pawHgts(bins);
            fit = polyfit(x, y, 1);
            corrs(i,j,k) = corr(x', y');
            slopes(i,j,k) = fit(1);
        end
    end
end


figure('position', [2000 400 1000 400], 'color', 'white', 'menubar', 'none');

% correlations
subplot(1,2,1)
barPlotRick(corrs, {'conditionNames', conditionNames, 'ylabel', 'paw/obstacle height correlation', ...
    'showStats', false, 'showViolins', true, 'lineThickness', 5, 'conditionColors', colors*.8, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'ylim', [0 1]})

% slopes
subplot(1,2,2)
barPlotRick(slopes, {'conditionNames', conditionNames, 'ylabel', 'paw/obstacle height slope', ...
    'showStats', false, 'showViolins', true, 'lineThickness', 5, 'conditionColors', colors*.8, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'ylim', [0 1.5]})
set(gca, 'YTick', 0:.5:1.5)

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineHeightShapingBars');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% BORING REPLICATIONS

rows = 2;
cols = 4;
figure('position', [2000 200 1600 700], 'color', 'white', 'menubar', 'none');

% initializations
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore paw', 'hind paw'}});
figVars = [vars.isFore; vars.isLeading];

conditionals.lightOn = struct('name', 'isLightOn', 'condition', @(x) x==1);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
figConditionals = [conditionals.lightOn];


% max height of paws
subplot(rows,cols,1)
matMod = getDvMatrix(data, 'stepOverMaxHgt', figVars, {'mouse'}, figConditionals) * 1000;
matBl = getDvMatrix(data, 'baselineStepHgt', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matBl,matMod), [1 2 4 3]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(ctlColor,8,1);
colorsWithBl(2:2:8,:) = colors;
barPlotRick(mat, {'conditionNames', conditionNames, 'ylabel', 'peak paw height (mm)', ...
    'showStats', false, 'showViolins', true, 'lineThickness', 5, 'conditionColors', colorsWithBl, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'showStats', true, 'groupSeparation', .25})

% x distance to obs at peak
subplot(rows,cols,2)
mat = getDvMatrix(data, 'xDistanceAtPeak', figVars, {'mouse'}, figConditionals) * 1000;
barPlotRick(mat, {'conditionNames', conditionNames, 'ylabel', 'horizontal distance at peak (mm)', ...
    'showViolins', true, 'lineThickness', 5, 'conditionColors', colors, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'showStats', true})

% x starting distance
subplot(rows,cols,3)
mat = abs(getDvMatrix(data, 'stepOverStartingDistance', figVars, {'mouse'}, figConditionals) * 1000);
barPlotRick(mat, {'conditionNames', conditionNames, 'ylabel', 'horizontal distance at start (mm)', ...
    'showViolins', true, 'lineThickness', 5, 'conditionColors', colors, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'showStats', true})

% x ending distance
subplot(rows,cols,4)
mat = getDvMatrix(data, 'stepOverEndingDistance', figVars, {'mouse'}, figConditionals) * 1000;
barPlotRick(mat, {'conditionNames', conditionNames, 'ylabel', 'horizontal distance at end (mm)', ...
    'showViolins', true, 'lineThickness', 5, 'conditionColors', colors, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'showStats', true})

% step over length
subplot(rows,cols,5)
mat = getDvMatrix(data, 'stepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
barPlotRick(mat, {'conditionNames', conditionNames, 'ylabel', 'step over length (mm)', ...
    'showViolins', true, 'lineThickness', 5, 'conditionColors', colors, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'showStats', true})

% length of steps leading up to obstacle
subplot(rows,cols,6:7)
matPrePre = getDvMatrix(data, 'prePreStepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
matPre = getDvMatrix(data, 'preStepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
mat = getDvMatrix(data, 'stepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matPrePre,matPre,mat), [1 2 4 3]); % add baseline vs mod steps as additional conditions

conditionNamesTemp = cat(2,conditionNames, {{'-2','-1','0'}});
colorsTemp = repelem(colors,3,1) .* repmat(linspace(.5,1,3)',4,1); % each color is repeated thrice, fading from black to the full color... lol
% close all; figure;
barPlotRick(mat, {'conditionNames', conditionNamesTemp, 'ylabel', 'step length (mm)', ...
    'showStats', false, 'showViolins', true, 'lineThickness', 5, 'conditionColors', colorsTemp, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'groupSeparation', .25, 'ylim' [20 120]})

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineBoringReplicationBars');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');




