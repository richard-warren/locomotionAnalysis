%% compute experiment data from scratch

sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'baselineNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);  % remove not included sessions and empty lines
mice = unique(sessionInfo.mouse);

data = cell(1,length(mice));
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data{1}.data = cellfun(@(x) x.data, data); data = data{1};
fprintf('saving...'); save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('data saved!')

%% load experiment data
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

% global settings
global_config;  % initialize global settings
isLeading = [true false true false]; % sequence of conditions for plots
isFore = [true true false false];
conditionNames = {{'fore paw', 'hind paw'}, {'leading', 'trailing'}};

%% ----------
% PLOT THINGS
%  ----------

%% show multiple frames with overlaid tracking

% settings
session = '180703_000';
% trials = [40 41 49];
trials = 40;

imgs = showTrackingOverFrames(session, trials, 2, 'showFig', false, ...
    'topOnly', false, 'contrastLims', [0 .8], 'alpha', .6, 'scatLines', true, 'scatSize', 100);
img = cat(1, imgs{:});
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'imgs', 'trackingOverFrames.png');
fprintf('writing %s to disk...\n', file);
imwrite(img, file)

%% show obs tracking of wheel velocity

% settings
showObsTracking('180703_000', 'numTrials', 15, 'waterColor', waterColor, 'obsColor', obsColor, ...
    'figPos', [2000 400 500 400], 'wheelColor', [.4 .4 .4])

file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'obsTracking');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% show example session with trial structure

% settings
session = '180714_004';
trialsToShow = 20;

plotSingleSessionVel(session, 'waterColor', waterColor, 'obsOnAlpha', obsAlpha, 'obsOnColor', obsColor, ...
    'trialsToShow', trialsToShow, 'trialColors', jet(trialsToShow), 'meanColor', axisColor);
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'singleSessionVel');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% show paw and whisker contacts

% settings
cols = 3;
preContactFrames = 2;

disp('creating contact images...')
showPawContactSeries('180703_000', preContactFrames, cols);
showWiskContactSeries('180803_003', 20, cols);


%% show single frame of example tracking

% settings
session = '180703_000';
trial = 40;
colorsTemp = stepColors([4 2 1 3], :);  % this is a hack that makes the colors align with leading, lagging, fore, hind conditions

showSingleFrameTracking(session, trial, ...
    'contrastLims', [0 .8], 'addWiskCam', true, 'pawColors', colorsTemp);
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'imgs', 'tracking.png');
fprintf('writing %s to disk...\n', file);
saveas(gcf, file)

%% leading, lagging, hind, fore schematic
global_config;
showLeadingLaggingImg('190318_000', 44, ...
    'colors', stepColors, 'contrastLims', [0 .6], 'pawPos', .008, ...
    'overlays', 8, 'overlayWidth', 3, 'overlayAlpha', .4, ...
    'randSeed', 1, 'scatter', false);

%% speed vs. position

% settings
yLims = [0 .8];
meanColor = [0 0 0];

% initializations
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'isWheelBreak', 'wiskContactPosition'});
flat = flat([flat.isLightOn]);
close all; figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 300], 'inverthardcopy', 'off'); hold on

% add obstacle rectangle and lines
x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
rectangle('Position', [x(1) yLims(1) diff(x) diff(yLims)], ...
    'FaceColor', [obsColor obsAlpha], 'EdgeColor', 'none');

% plot
velData = plotDvPsth(flat, 'velVsPosition', [], ...
    {'plotMouseAvgs', true, 'showLegend', false, 'conditionColors', [0 0 0], 'errorFcn', @(x) nanstd(x), 'xlim', [-.5 .2]});
line([0 0], yLims, 'linewidth', 2, 'color', get(gca, 'XColor'));

% pimp fig
set(gca, 'YLim', yLims);
xlabel('position relative to nose (m)')
ylabel('velocity (m/s)')
text(x(1), yLims(2), 'obstace engaged', 'VerticalAlignment', 'bottom')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineVel');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

% get vel at moment obstacle is under nose
atNoseInd = find(flat(1).velVsPositionX>=0,1,'first');
noseVels = squeeze(velData(1,1,:,atNoseInd));
obsOnInd = find(flat(1).velVsPositionX>=x(1),1,'first');
obsOnVels = squeeze(velData(1,1,:,obsOnInd));

fprintf('\nobs at nose: %.2f +- %.2f SEM\n', mean(noseVels), std(noseVels)/sqrt(length(noseVels)))
fprintf('obs on:      %.2f +- %.2f SEM\n', mean(obsOnVels), std(obsOnVels)/sqrt(length(obsOnVels)))


%% kinematics

% settings
colNames = {'hind', 'fore'};
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
close all
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 700 1600 250])

plotInd = 1;
for i = conditionSequence
    subplot(2,2,plotInd)
    plotColor = repmat(stepColors(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    obsColors = repmat(obsColor, obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)';
    bins = [flat.isLeading]==isLeading(i) & ...
           [flat.isFore]==isFore(i);

    % plot control step
    plotKinematics(kinDataCtl(bins,[1,3],:), [flat(bins).obsHgt], ones(1,sum(bins)), ...
        'colors', ctlStepColor, 'obsAlpha', 0, 'lineAlpha', .8, 'plotObs', false) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    
    plotKinematics(kinData(bins,[1,3],:), [flat(bins).obsHgt], [flat(bins).obsHgtDiscretized], ...
        'colors', plotColor, 'obsAlpha', 1, 'lineAlpha', .8, 'obsColors', obsColors, 'mouseNames', {flat(bins).mouse}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
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
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'baselineKinematics_obsHgt');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% TOP VIEW OVERLAYS

% initializations
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [1997.00 99.00 878.00 833.00])
subplot(2,1,1)

% get condition numbers, where each condition is unique combinate of isLeading and isFore
a = [isLeading; isFore]'; % each row is a condition (00,01,10,11)
b = [[flat.isLeading]; [flat.isFore]]';
[~, stepTypeConditions] = ismember(b, a, 'rows');

plotKinematics(kinData(:,[1,3],:), [flat.obsHgt], stepTypeConditions', ...
    'colors', stepColors, 'obsAlpha', .25, 'lineAlpha', 1, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @nanstd, 'lineWidth', 3, 'obsColors', repmat(obsColors, 4, 1), 'yLimZero', false) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
set(gca, 'XLim', xLims)

% % save
% file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'baselineKinematics_overlayTop');
% fprintf('writing %s to disk...\n', file)
% saveas(gcf, file, 'svg');

% BOT VIEW OVERLAYS
subplot(2,1,2)
yLimsBot = [-1 1]*.02;

% figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 100 900 300])
plotKinematics(kinData(:,[1,2],:), [flat.obsHgt], stepTypeConditions', ...
    'colors', stepColors, 'obsAlpha', .25, 'lineAlpha', 1, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @nanstd, 'lineWidth', 3, 'isBotView', true, 'obsColors', repmat(obsColors, 4, 1)) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
set(gca, 'XLim', xLims, 'YLim', yLimsBot)

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'baselineKinematics_overlay');
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
    'obsHgt', 'preObsHgt', 'isFore', 'isLeading', 'stepOverMaxHgt', 'isPawSuccess'});
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
    {'colors', stepColors, 'maxScatterPoints', 5000, 'lineAlpha', 1, 'scatAlpha', .1, 'scatSize', 40});
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
    'showStats', false, 'showViolins', true, 'lineThickness', 5, 'conditionColors', stepColors*.8, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'ylim', [0 1]})

% slopes
subplot(1,2,2)
barPlotRick(slopes, {'conditionNames', conditionNames, 'ylabel', 'paw/obstacle height slope', ...
    'showStats', false, 'showViolins', true, 'lineThickness', 5, 'conditionColors', stepColors*.8, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'ylim', [0 1.5]})
set(gca, 'YTick', 0:.5:1.5)

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineHeightShapingBars');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% OBS HGT PAW HGT MOVING AVG

figure('Color', 'white', 'Position', [2000 400 500 400], 'MenuBar', 'none');
plot([0 xLims(2)], [0 xLims(2)], 'Color', [1 1 1]*.6, 'LineWidth', 3) % add unity line
logPlotRick(obsHgts, pawHgts, ...
    {'colors', stepColors, 'conditions', stepTypeConditions, 'xlabel', 'obstacle height', 'ylabel', 'paw height (mm)', 'plotMice', false, ...
    'xlim', [3.4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1))})
set(gca, 'xlim', [4 10])

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baseline_heightShapingMovingAvgs');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

%% SUCCESS BY OBS HEIGHT
close all
figure('Color', 'white', 'Position', [2000 400 400 300], 'MenuBar', 'none');
logPlotRick([flat.obsHgt]*1000, [flat.isPawSuccess], ...
    {'colors', stepColors, 'conditions', stepTypeConditions, 'xlabel', 'obstacle height', 'ylabel', 'success rate', 'plotMice', false, ...
    'xlim', [3.4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'computeVariance', false, 'ylim', [0 1]})
set(gca, 'xlim', [3.4 10])

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baseline_succesByObsHgt');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% BORING REPLICATIONS

rows = 3;
cols = 4;
figure('position', [2000 200 1600 700], 'color', 'white', 'menubar', 'none');

% initializations
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore', 'hind'}});
figVars = [vars.isFore; vars.isLeading];

conditionals.lightOn = struct('name', 'isLightOn', 'condition', @(x) x==1);
conditionals.lightOff = struct('name', 'isLightOn', 'condition', @(x) x==0);
conditionals.noWheelBreak = struct('name', 'isWheelBreak', 'condition', @(x) x==0);
conditionals.isLagging = struct('name', 'isLeading', 'condition', @(x) x==0);
% figConditionals = [conditionals.lightOn];
figConditionals = [conditionals.lightOff];


% max height of paws
subplot(rows,cols,1)
matMod = getDvMatrix(data, 'stepOverMaxHgt', figVars, {'mouse'}, figConditionals) * 1000;
matBl = getDvMatrix(data, 'controlStepHgt', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matBl,matMod), [1 2 4 3]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(ctlStepColor,8,1);
colorsWithBl(2:2:8,:) = stepColors;
barPlotRick(mat, {'conditionNames', conditionNames, 'ylabel', 'peak paw height (mm)', ...
    'showStats', false, 'showViolins', true, 'lineThickness', 5, 'conditionColors', colorsWithBl, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'showStats', true, 'groupSeparation', .25})

% x distance to obs at peak
subplot(rows,cols,2)
mat = getDvMatrix(data, 'xDistanceAtPeak', figVars, {'mouse'}, figConditionals) * 1000;
barPlotRick(mat, {'conditionNames', conditionNames, 'ylabel', 'horizontal distance at peak (mm)', ...
    'showViolins', true, 'lineThickness', 5, 'conditionColors', stepColors, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'showStats', true})

% x starting distance
subplot(rows,cols,3)
mat = abs(getDvMatrix(data, 'stepOverStartingDistance', figVars, {'mouse'}, figConditionals) * 1000);
barPlotRick(mat, {'conditionNames', conditionNames, 'ylabel', 'horizontal distance at start (mm)', ...
    'showViolins', true, 'lineThickness', 5, 'conditionColors', stepColors, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'showStats', true})

% x ending distance
subplot(rows,cols,4)
mat = getDvMatrix(data, 'stepOverEndingDistance', figVars, {'mouse'}, figConditionals) * 1000;
barPlotRick(mat, {'conditionNames', conditionNames, 'ylabel', 'horizontal distance at end (mm)', ...
    'showViolins', true, 'lineThickness', 5, 'conditionColors', stepColors, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'showStats', true})

% step over length
subplot(rows,cols,5)
mat = getDvMatrix(data, 'stepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
barPlotRick(mat, {'conditionNames', conditionNames, 'ylabel', 'step over length (mm)', ...
    'showViolins', true, 'lineThickness', 5, 'conditionColors', stepColors, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'showStats', true})

% length of steps leading up to obstacle
subplot(rows,cols,6:7)
matPrePre = getDvMatrix(data, 'prePreStepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
matPre = getDvMatrix(data, 'preStepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
mat = getDvMatrix(data, 'stepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matPrePre,matPre,mat), [1 2 4 3]); % add baseline vs mod steps as additional conditions

conditionNamesTemp = cat(2,conditionNames, {{'-2','-1','0'}});
colorsTemp = repelem(stepColors,3,1) .* repmat(linspace(.5,1,3)',4,1); % each color is repeated thrice, fading from black to the full color... lol
% close all; figure;
barPlotRick(mat, {'conditionNames', conditionNamesTemp, 'ylabel', 'step length (mm)', ...
    'showStats', false, 'showViolins', true, 'lineThickness', 5, 'conditionColors', colorsTemp, ...
    'violinAlpha', .1, 'scatColors', [.6 .6 .6], 'scatAlpha', .3, 'groupSeparation', .25, 'ylim' [20 120]})

% pre obs height of paws
subplot(rows,cols,8)
matMod = getDvMatrix(data, 'preObsHgt', figVars, {'mouse'}, figConditionals) * 1000;
matBl = getDvMatrix(data, 'controlPreObsHgt', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matBl,matMod), [1 2 4 3]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(ctlStepColor,8,1);
colorsWithBl(2:2:8,:) = stepColors;
barPlotRick(mat, {'conditionNames', conditionNames, 'ylabel', 'peak pre obs height (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colorsWithBl, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'lineWidth', 1.5})

% step length vs control step length
subplot(rows,cols,9)
matMod = getDvMatrix(data, 'stepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
matBl = getDvMatrix(data, 'controlStepLength', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matBl,matMod), [1 2 4 3]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(ctlStepColor,8,1);
colorsWithBl(2:2:8,:) = stepColors;
barPlotRick(mat, {'conditionNames', conditionNames, 'ylabel', 'step over length (mm)', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colorsWithBl, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', true, 'lineWidth', 1.5})

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineBoringReplicationBars');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');




