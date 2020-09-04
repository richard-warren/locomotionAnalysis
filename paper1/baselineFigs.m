%% compute experiment data from scratch

sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'baselineNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);  % remove not included sessions and empty lines
mice = unique(sessionInfo.mouse);

data = cell(1,length(mice));
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data{1}.data = cellfun(@(x) x.data, data); data = data{1};
fprintf('saving...'); save(fullfile(getenv('SSD'), 'paper1', 'baseline_data.mat'), 'data'); disp('data saved!')

%% load experiment data
fprintf('loading... '); load(fullfile(getenv('SSD'), 'paper1', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

% global settings
paper1_config;  % initialize global settings
extraBarProps = {'constantEdgeColor', [.2 .2 .2], 'scatterSize', 30};  % these data have large sample size, so some of the bar plot properties will be unique to this dataset
isLeading = [true false true false]; % sequence of conditions for plots
isFore = [true true false false];
conditionNames = {{'fore paw', 'hind paw'}, {'leading', 'trailing'}};

% initializations
vars.isLeading = struct('name', 'isLeading', 'levels', [1 0], 'levelNames', {{'leading', 'trailing'}});
vars.isBigStep = struct('name', 'isBigStep', 'levels', [0 1], 'levelNames', {{'little step', 'big step'}});
vars.isFore = struct('name', 'isFore', 'levels', [1 0], 'levelNames', {{'fore', 'hind'}});
vars.isLightOn = struct('name', 'isLightOn', 'levels', [0 1], 'levelNames', {{'light off', 'light on'}});
figVars = [vars.isFore; vars.isLeading];
figConditionals = struct('name', 'isLightOn', 'condition', @(x) x==1);
noConditionals = struct('name', '', 'condition', @(x) x);

%% ----------
% PLOT THINGS
%  ----------

%% hildebrand plots

sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'baselineNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);  % remove not included sessions and empty lines

plotHildebrands(sessionInfo, 'colors', stepColors, 'stepPercentiles', [1 99], 'plotMice', false)
set(gcf, 'Position', [1965.00 821.00 632.00 155.00])
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baseline_hildebrand'), 'svg');


%% show multiple frames with overlaid tracking

% settings
session = '180703_000';
trials = 40;

imgs = showTrackingOverFrames(session, trials, 1, 'showFig', true, ...
    'topOnly', false, 'contrastLims', contrast, 'alpha', .6, 'scatLines', true, 'scatSize', 100, 'pawColors', stepColors, ...
    'imgSpacing', 0);
img = cat(1, imgs{:});
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'imgs', 'trackingOverFrames.png');
fprintf('writing %s to disk...\n', file);
imwrite(img, file)

%% show obs tracking of wheel velocity

% settings
showObsTracking('180703_000', 'numTrials', 15, 'waterColor', waterColor, 'obsColor', obsColor, ...
    'figPos', [2000 400 500 400], 'wheelColor', [.4 .4 .4])

file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'obsTracking');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% show example session with trial structure

% settings
session = '180714_004';
trialsToShow = 20;

plotSingleSessionVel(session, 'waterColor', waterColor, 'obsOnAlpha', obsAlpha, 'obsOnColor', obsColor, ...
    'trialsToShow', trialsToShow, 'trialColors', jet(trialsToShow), 'meanColor', axisColor);
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'singleSessionVel');
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
    'contrastLims', contrast, 'addWiskCam', true, 'pawColors', colorsTemp);

file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'imgs', 'tracking.png');
fprintf('writing %s to disk...\n', file);
saveas(gcf, file)

%% leading, lagging, hind, fore schematic

showLeadingLaggingImg('190318_000', 44, ...
    'colors', stepColors, 'contrastLims', contrast, 'pawPos', .008, ...
    'overlays', 8, 'overlayWidth', 3, 'overlayAlpha', .4, ...
    'randSeed', 1, 'scatter', false, 'vertical', false);

%% obstacle contact time distributions

flat = flattenData(data, {'obsHgt', 'touchFrames'});

% settings
hgtBinNum = 4;
hgtLims = [4 10];
xCenters = 0:4:20;
yLims = [0 1];


hgtBinEdges = linspace(hgtLims(1), hgtLims(2), hgtBinNum+1)/1000;
hgtConditions = discretize([flat.obsHgt], hgtBinEdges);
dx = diff(xCenters(1:2));
xEdges = [xCenters xCenters(end)+dx] - dx/2;
colors = repmat(obsColor, hgtBinNum, 1) .* linspace(.2,1,hgtBinNum)';


close all;
figure('color', 'white', 'position', [1995.00 565.00 364.00 326.00], 'menubar', 'none'); hold on;

for i = 1:hgtBinNum
    
    h = histcounts([flat(hgtConditions==i).touchFrames], xEdges, ...
        'Normalization', 'probability');
    h2 = histcounts([flat(hgtConditions==i).touchFrames], xEdges, ...
        'Normalization', 'cdf');
    
    scatter(xCenters, h2, 25, colors(i,:), 'filled');
    lines(i) = plot(xCenters, h2, 'Color', colors(i,:), 'LineWidth', 2);
    
%     fill([xCenters(1) xCenters xCenters(end) xCenters(1)], ...
%          [yLims(1) h yLims(1) yLims(1)], colors(i,:), ...
%          'FaceAlpha', .2, 'EdgeColor', 'none')
    
end

set(gca, 'XTick', linspace(xCenters(1),xCenters(end),3), 'YLim', yLims, 'YTick', 0:.5:1, 'TickDir', 'out')
xlabel('contact duration (ms)')
ylabel('cumulative density')

hgtCenters = (hgtBinEdges(1:end-1) + (diff(hgtBinEdges(1:2))/2)) * 1000;
hgtCenters = strsplit(sprintf('%.1f mm,', hgtCenters), ',');
hgtCenters = hgtCenters(1:end-1);
legend(lines, hgtCenters, 'Location', 'best', 'Box', 'off')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'baselineContactDurations');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

%% speed vs. position

% settings
yLims = [0 .8];

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
    {'plotMouseAvgs', true, 'showLegend', false, 'conditionColors', [0 0 0], ...
    'errorFcn', @(x) nanstd(x), 'xlim', [-.5 .2], 'flipx', true});
line([0 0], yLims, 'linewidth', 2, 'color', get(gca, 'XColor'));

% pimp fig
set(gca, 'YLim', yLims);
xlabel('hurdle distance to nose (m)')
ylabel('velocity (m/s)')
text(x(1), yLims(2), 'obstace engaged', 'VerticalAlignment', 'bottom')


% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineVel');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

% get vel at moment obstacle is under nose
atNoseInd = find(flat(1).velVsPositionX>=0,1,'first');
noseVels = squeeze(velData(:,atNoseInd));
obsOnInd = find(flat(1).velVsPositionX>=x(1),1,'first');
obsOnVels = squeeze(velData(:,obsOnInd));

fprintf('\nobs at nose: %.2f +- %.2f SEM\n', mean(noseVels), std(noseVels)/sqrt(length(noseVels)))
fprintf('obs on:      %.2f +- %.2f SEM\n', mean(obsOnVels), std(obsOnVels)/sqrt(length(obsOnVels)))
fprintf('overall:     %.2f +- %.2f SEM\n', mean(nanmean(velData,2)), std(nanmean(velData,2))/sqrt(length(obsOnVels)))


%% speed triggered at obs on and wisk contact

% initializations
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions', ...
    'velContinuousAtContact', 'velContinuousAtContactX', 'isWheelBreak'});
flat = flat([flat.isLightOn]);
close all; figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000.00 50.00 360.00 300.00], 'inverthardcopy', 'off'); hold on

% plot
velData = plotDvPsth(flat, 'velContinuousAtContact', [], ...
    {'plotMouseAvgs', true, 'showLegend', false, 'conditionColors', [0 0 0], 'errorFcn', @(x) nanstd(x), 'xlim', [-.5 .2]});
line([0 0], [0 .8], 'linewidth', 2, 'color', obsColor);

% pimp fig
set(gca, 'YTick', 0:.4:.8, 'XLim', [-.4 .2], 'XTick', -.4:.2:.2);
xlabel('time from whisker contact (s)')
ylabel('velocity (m/s)')
text(0, .8, 'whisker contact', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'baselineVelAtContact');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


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
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineKinematics_obsHgt');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');



%% bot and top view overlays
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [1997.00 99.00 878.00 833.00])

% top view
subplot(2,1,1)

% get condition numbers, where each condition is unique combinate of isLeading and isFore
a = [isLeading; isFore]'; % each row is a condition (00,01,10,11)
b = [[flat.isLeading]; [flat.isFore]]';
[~, stepTypeConditions] = ismember(b, a, 'rows');

plotKinematics(kinData(:,[1,3],:), [flat.obsHgt], stepTypeConditions', ...
    'colors', stepColors, 'obsAlpha', 1, 'lineAlpha', 1, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @nanstd, 'lineWidth', 3, 'obsColors', repmat(obsColor, 4, 1), 'yLimZero', false) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
set(gca, 'XLim', xLims)

% bot view
subplot(2,1,2)
yLimsBot = [-1 1]*.02;

% figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 100 900 300])
plotKinematics(kinData(:,[1,2],:), [flat.obsHgt], stepTypeConditions', ...
    'colors', stepColors, 'obsAlpha', .25, 'lineAlpha', 1, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @nanstd, 'lineWidth', 3, 'isBotView', true, 'obsColors', repmat(obsColor, 4, 1)) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
set(gca, 'XLim', xLims, 'YLim', yLimsBot)

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineKinematics_overlay');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% height shaping

% settings
xLims = [3 10];
isHgtPreObs = true; % do you measure the height at peak (false) or before the paw reaches the obstacle (true)

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isLeading', 'stepOverMaxHgt', 'isPawSuccess', 'controlStepHgt'});
flat = flat([flat.isLightOn] & ...
            ~isnan([flat.obsHgt]) & ...
            ~isnan([flat.stepOverMaxHgt]) & ...
            ~isnan([flat.preObsHgt])); % add conditionals here

% get condition numbers, where each condition is unique combination of isLeading and isFore
a = [isLeading; isFore]';  % each row is a condition (00,01,10,11)
b = [[flat.isLeading]; [flat.isFore]]';
[~, stepTypeConditions] = ismember(b, a, 'rows');

if isHgtPreObs; pawHgts = [flat.preObsHgt]*1000; else; pawHgts = [flat.stepOverMaxHgt]*1000; end
pawHgtsControl = [flat.controlStepHgt]*1000;
obsHgts = [flat.obsHgt]*1000;



% bars
mice = unique({flat.mouse});
[corrs, corrsControl, slopes] = deal(nan(2,2,length(mice))); % isFore(10) X isLeading(10) X mouse
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
            yControl = pawHgtsControl(bins);
            fit = polyfit(x, y, 1);
            corrs(i,j,k) = corr(x', y');
            corrsControl(i,j,k) = corr(x', yControl');
            slopes(i,j,k) = fit(1);
        end
    end
end



%% correlations
figure('position', [200 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');

colorsWithBl = repelem(stepColors,2,1);
colorsWithBl(1:2:8,:) = colorsWithBl(2:2:8,:) * .25;

mat = permute(cat(4, corrsControl, corrs), [1 2 4 3]);

barFancy(mat, 'levelNames', conditionNames, 'ylabel', 'paw-obstacle correlation', ...
    'colors', colorsWithBl, 'YLim', [-.2 .8], barProperties{:}, extraBarProps{:}, ...
    'comparisons', [1 2; 3 4; 5 6; 7 8], 'test', 'ttest', 'lineAtZero', false)

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineHeightShapingBars');
saveas(gcf, file, 'svg');


%% moving averages
figure('Color', 'white', 'Position', [200 400 500 400], 'MenuBar', 'none');
plot([0 xLims(2)], [0 xLims(2)], 'Color', [0 0 0 .25], 'LineWidth', 3, 'LineStyle', ':') % add unity line

logPlotRick(obsHgts, pawHgts, ...
    'colors', stepColors, 'conditions', stepTypeConditions, 'xlabel', 'obstacle height (mm)', 'ylabel', 'paw height (mm)', 'plotMice', false, ...
    'xlim', [3.4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)))

% uncomment the following to show control step heights as well...
% logPlotRick(obsHgts, pawHgtsControl, ...
%     {'colors', (.7*ctlStepColor + .3*stepColors), 'conditions', stepTypeConditions, 'xlabel', 'obstacle height (mm)', 'ylabel', 'paw height (mm)', 'plotMice', false, ...
%     'xlim', [3.4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat.mouse}, ...
%     'errorFcn', @(x) std(x)/sqrt(size(x,1))})

set(gca, 'xlim', [3.4 10], 'ylim', [3.4 15], 'YTick', 4:4:16)

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineHeightShapingMovingAvgs');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% success by obs height

figure('Color', 'white', 'Position', [2000 400 400 300], 'MenuBar', 'none');
logPlotRick([flat.obsHgt]*1000, [flat.isPawSuccess], ...
    'colors', stepColors, 'conditions', stepTypeConditions, 'xlabel', 'obstacle height (mm)', 'ylabel', 'success rate', 'plotMice', false, ...
    'xlim', [3.4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1)), 'computeVariance', false, 'ylim', [0 1])
set(gca, 'xlim', [3.4 10], 'XTick', 4:10, 'YTick', [0 .5 1])
xLabels = get(gca, 'XTickLabel');
xLabels(2:end-1) = {''};
set(gca, 'XTickLabel', xLabels)
legend({'leading fore', 'trailing fore', 'leading hind', 'trailing hind'}, 'Location', 'southeast', 'box', 'off')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baseline_succesByObsHgt');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% bar plots


% step height
figure('position', [200 617 279 328], 'color', 'white', 'menubar', 'none');
matMod = getDvMatrix(data, 'stepOverMaxHgt', figVars, {'mouse'}, figConditionals) * 1000;
matBl = getDvMatrix(data, 'controlStepHgt', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matBl,matMod), [1 2 4 3]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(stepColors,2,1);
colorsWithBl(1:2:8,:) = colorsWithBl(2:2:8,:) * .25;
barFancy(mat, 'levelNames', conditionNames, 'ylabel', 'step height (mm)', ...
    'colors', colorsWithBl, barProperties{:}, extraBarProps{:}, 'YTick', [2 10 18], 'scatterAlpha', .6, ...
    'comparisons', [1 2; 3 4; 5 6; 7 8], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineStepHeight'), 'svg');

% factor by which each paw increases in height
fprintf('\nheight increase factor -> LF %.1f, TF %.1f, LH %.1f, TH %.1f, overall %.1f\n', ...
    mean(mat(1,1,2,:)) / mean(mat(1,1,1,:)), ...
    mean(mat(1,2,2,:)) / mean(mat(1,1,1,:)), ...
    mean(mat(2,1,2,:)) / mean(mat(1,1,1,:)), ...
    mean(mat(2,2,2,:)) / mean(mat(1,1,1,:)), ...
    mean(reshape(mat(:,:,2,:),1,[])) / mean(reshape(mat(:,:,1,:),1,[])))

%% step length
figure('position', [200 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
matMod = getDvMatrix(data, 'stepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
matBl = getDvMatrix(data, 'controlStepLength', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matBl,matMod), [1 2 4 3]); % add baseline vs mod steps as additional conditions
colorsWithBl = repelem(ctlStepColor,8,1);
colorsWithBl(2:2:8,:) = stepColors;
barFancy(mat, 'levelNames', conditionNames, 'ylabel', 'step length (mm)', ...
    'colors', colorsWithBl, barProperties{:}, extraBarProps{:})
set(gca, 'YTick', 0:60:120, 'TickDir', 'out')
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineStepLength');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

fprintf('\nlength increase factor -> LF %.0f%%, TF %.0f%%, LH %.0f%%, TH %.0f%%, overall %.0f%%\n', ...
    (mean(mat(1,1,2,:)) / mean(mat(1,1,1,:)) - 1)*100, ...
    (mean(mat(1,2,2,:)) / mean(mat(1,1,1,:)) - 1)*100, ...
    (mean(mat(2,1,2,:)) / mean(mat(1,1,1,:)) - 1)*100, ...
    (mean(mat(2,2,2,:)) / mean(mat(1,1,1,:)) - 1)*100, ...
    (mean(reshape(mat(:,:,2,:),1,[])) / mean(reshape(mat(:,:,1,:),1,[])) - 1)*100)


%% starting horizontal distance
figure('position', [2400 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
mat = getDvMatrix(data, 'stepOverStartingDistance', figVars, {'mouse'}, figConditionals) * -1000;
barFancy(mat, 'levelNames', conditionNames, 'ylabel', 'horizontal distance at lift (mm)', ...
    'colors', stepColors, barProperties{:}, extraBarProps{:})
set(gca, 'YTick', 0:20:80, 'TickDir', 'out')
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'baselineDistAtLift');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% ending horizontal distance
figure('position', [2600 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
mat = getDvMatrix(data, 'stepOverEndingDistance', figVars, {'mouse'}, figConditionals) * 1000;
barFancy(mat, 'levelNames', conditionNames, 'ylabel', 'horizontal distance at land (mm)', ...
    'colors', stepColors, barProperties{:}, extraBarProps{:})
set(gca, 'YTick', 0:20:80, 'TickDir', 'out')
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'baselineDistAtLand');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% linear distance at zenith
figure('position', [2800 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
mat = getDvMatrix(data, 'distanceToObs', figVars, {'mouse'}, figConditionals) * 1000;
barFancy(mat, 'levelNames', conditionNames, 'ylabel', 'distance to obstalce at zenith (mm)', ...
    'colors', stepColors, barProperties{:}, extraBarProps{:})
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'baselineDistAtZenith');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% horizontal distance at zenith
figure('position', [3000 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
mat = abs(getDvMatrix(data, 'xDistanceAtPeak', figVars, {'mouse'}, figConditionals) * 1000);
barFancy(mat, 'levelNames', conditionNames, 'ylabel', 'horizontal distance to obstalce at zenith (mm)', ...
    'YLim', [0 12], 'colors', stepColors, barProperties{:}, extraBarProps{:})
set(gca, 'YTick', 0:4:12)
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'baselineHorDistAtZenith');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

fprintf('\nhorizontal distance -> LF %.1f, TF %.1f, LH %.1f, TH %.1f, overall %.1f\n', ...
    mean(mat(1,1,:)), ...
    mean(mat(1,2,:)), ...
    mean(mat(2,1,:)), ...
    mean(mat(2,2,:)), ...
    mean(reshape(mat(:,:,:),1,[])));


%% length of steps leading up to obstacle
figure('position', [200 292 851 447], 'color', 'white', 'menubar', 'none');
matPrePre = getDvMatrix(data, 'prePreStepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
matPre = getDvMatrix(data, 'preStepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
mat = getDvMatrix(data, 'stepOverLength', figVars, {'mouse'}, figConditionals) * 1000;
mat = permute(cat(4,matPrePre,matPre,mat), [1 2 4 3]); % add baseline vs mod steps as additional conditions

conditionNamesTemp = cat(2,conditionNames, {{'-2','-1','0'}});
colorsTemp = repelem(stepColors,3,1) .* repmat([.2; .6; 1],4,1); % each color is repeated thrice, fading from black to the full color... lol
barFancy(mat, 'levelNames', conditionNamesTemp, 'ylabel', 'step length (mm)', ...
    'colors', colorsTemp, barProperties{:}, extraBarProps{:}, 'YTick', 20:50:120, ...
    'comparisons', [1 3; 4 6; 7 9; 10 12], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselinePrecedingStepLengths'), 'svg');

%% light on vs light off speed


% settings
yLims = [.3 .7];
plotSequence = [4 3 2 1]; % determine which lines on plotted on top of which

% initializations
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsOnPositions', 'obsOffPositions', ...
    'velVsPosition', 'velVsPositionX', 'velContinuousAtContact', 'velContinuousAtContactX', 'isWheelBreak', 'wiskContactPosition'});
colorsTemp = [sensColors(1:end-1,:); .6 .6 .6]; % the no vision no whisker condition can be a little lighter here

% speed vs position
figure('name', 'baseline', 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600 300], 'inverthardcopy', 'off')
x = [nanmean([flat.obsOnPositions]) nanmean([flat.obsOffPositions])];
rectangle('Position', [x(1) yLims(1) diff(x) diff(yLims)], ...
    'FaceColor', [obsColor obsAlpha], 'EdgeColor', 'none');
line([0 0], yLims, 'linewidth', 2, 'color', get(gca, 'XColor'))
velData = plotDvPsth(flat, 'velVsPosition', 'isLightOn', ...
    {'showLegend', false, 'conditionColors', colorsTemp, 'xlim', [-.5 .2], ... 
     'errorAlpha', .1, 'lineWidth', 4});
set(gca, 'YLim', yLims, 'YTick', linspace(yLims(1),yLims(2),3));
xlabel('position relative to nose (m)')
ylabel('velocity (m/s)')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'baselineVelLightOnVsLightOff');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% get vel at moment obstacle is under nose
atNoseInd = find(flat(1).velVsPositionX>=0,1,'first');
noseVels = squeeze(velData(:,:,atNoseInd));

fprintf('\nobs at nose: %.2f +- %.2f SEM\n', mean(noseVels), std(noseVels)/sqrt(length(noseVels)))


%% histogram showing accuracy of eddie's whisker contact classifier

data = 'Y:\loco\obstacleData\papers\hurdles_paper1\figures\eddie_figs\200414_files\histdata.csv';
histData = readtable(data);
bins = -244:8:54;

close all
figure('color', 'white', 'position', [480 343 457 300], 'menubar', 'none'); hold on;

props = {'Normalization', 'probability', 'FaceAlpha', .4};
histogram(histData.train(~isnan(histData.train)), bins, 'FaceColor', [.4 .4 .4], props{:});
histogram(histData.test(~isnan(histData.test)), bins, 'FaceColor', [20 96 168]/255, props{:});
yLims = get(gca, 'ylim');
% plot([0 0], yLims, 'Color', 'black')
set(gca, 'ylim', yLims, 'xlim', bins([1 end]), 'YColor', 'none')
xlabel('true - predicted contact time (ms)')
xlabel('true - predicted contact time (ms)')
legend({'train', 'test'}, 'Box', 'off', 'Location', 'Best')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerContactHisto'), 'svg');


%% compute accuracy of step length models

sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'baselineNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);  % remove not included sessions and empty lines


[corrs, mae] = deal(nan(1,length(sessions)));  % correlations, mean absolute error
sessions = sessionInfo.session;
for i = 1:length(sessions)
    fprintf('computing session %s...\n', sessions{i})
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'kinData.mat'), 'models')
    
    % left forepaw
    lf_mae = mean(abs(models{2}.predict - models{2}.Variables.y)*1000);
    lf_corr = corr(models{2}.Variables.x1, models{2}.Variables.y);
    
    % right forpaw
    rf_mae = mean(abs(models{3}.predict - models{3}.Variables.y)*1000);
    rf_corr = corr(models{3}.Variables.x1, models{3}.Variables.y);
    
    % averages
    corrs(i) = mean([lf_corr rf_corr]);
    mae(i) = mean([lf_mae rf_mae]);
end

fprintf('\ncorrelation:           mean %.3f, sem %.3f\n', mean(corrs), std(corrs)/sqrt(length(corrs)))
fprintf('mean absolute error:   mean %.3f, sem %.3f\n', mean(mae), std(mae)/sqrt(length(mae)))



%% TEST DISTRIBUTION NORMALITY


%% correlations (loop across paws)
for i = 1:2
    for j = 1:2
        [h, p] = kstest(zscore(squeeze(corrs(i,j,:))));
        fprintf('%i, %.5f: corr\n', h, p)
    end
end

%% paw heights
matMod = getDvMatrix(data, 'stepOverMaxHgt', figVars, {'mouse'}, noConditionals) * 1000;
matBl = getDvMatrix(data, 'controlStepHgt', figVars, {'mouse'}, noConditionals) * 1000;
mat = permute(cat(4,matBl,matMod), [1 2 4 3]); % add baseline vs mod steps as additional conditions
for i = 1:2
    for j = 1:2
        for k = 1:2
            [h, p] = kstest(zscore(squeeze(mat(i,j,k,:))));
            fprintf('%i, %.5f: paw height\n', h, p)
        end
    end
end

%% success
mat = getDvMatrix(data, 'isTrialSuccess', vars.isLightOn, {'mouse'}, noConditionals);
mat = mean(mat,1);  % average across light on and light off
[h, p] = kstest(zscore(mat));
fprintf('%i, %.5f: success\n', h, p)

%% velocity
mat = getDvMatrix(data, 'trialVel', vars.isLightOn, {'mouse'}, noConditionals);
mat = mean(mat,1);  % average across light on and light off
[h, p] = kstest(zscore(mat));
fprintf('%i, %.5f: velocity\n', h, p)

%% body angle
mat = getDvMatrix(data, 'trialAngle', vars.isLightOn, {'mouse'}, noConditionals);
mat = mean(mat,1);  % average across light on and light off
[h, p] = kstest(zscore(mat));
fprintf('%i, %.5f: body angle\n', h, p)

%% tail height
mat = getDvMatrix(data, 'tailHgt', vars.isLightOn, {'mouse'}, noConditionals);
mat = mean(mat,1);  % average across light on and light off
[h, p] = kstest(zscore(mat));
fprintf('%i, %.5f: tail height\n', h, p)

%% distance to obstable
mat = getDvMatrix(data, 'modPawDistanceToObs', vars.isBigStep, {'mouse'}, noConditionals);

for i = 1:2
    [h, p] = kstest(zscore(mat(i,:)));
    fprintf('%i, %.5f: distance to obs\n', h, p)
end

%% model accuracies
flat = flattenData(data, ...
    [m.predictorsAll, {'mouse', 'session', 'trial', 'isModPawLengthened', 'modPawDeltaLength', 'isBigStep', 'isLightOn', ...
    'modPawOnlySwing', 'isTrialSuccess', 'modPawPredictedDistanceToObs', 'modPawDistanceToObs', 'modPawKinInterp', ...
    'preModPawKinInterp', 'firstModPaw', 'contactIndInterp', 'preModPawDeltaLength', 'contactInd', 'modSwingContacts'}]);

mat = plotModelAccuracies(flat, m.predictors, 'isModPawLengthened', 'model', 'glm', ...
    'modSwingContactsMax', m.modSwingContactsMax, 'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'weightClasses', true, 'barProps', barProperties, 'kFolds', 15);

for i = 1:2
    [h, p] = kstest(zscore(squeeze(mat(1,i,:))));
    fprintf('%i, %.5f: model accuracy\n', h, p)
end

%% entropy
mat = plotEntropies(flat);

for i = 1:2
    [h, p] = kstest(zscore(squeeze(mat(i,1,:))));
    fprintf('%i, %.5f: entropy\n', h, p)
end

%% lagging paw landing position variability

flat = flattenData(data, {'mouse', 'session', 'trial', 'laggingPenultKin'});
mice = unique({flat.mouse});
landingPositions = cellfun(@(x) x(1,end), {flat.laggingPenultKin});

mat = nan(1, length(mice));  % variability in planting paw distance to hurdle // (sensory condition) X (mouse)
for j = 1:length(mice)
    bins = strcmp({flat.mouse}, mice{j});
    mat(j) = nanstd(landingPositions(bins));
end

[h, p] = kstest(zscore(mat));
fprintf('%i, %.5f: lagging paw landing position variability\n', h, p)


%% compute correlation between whisker contact position and trial velocity

flat = flattenData(data, {'mouse', 'trialVel', 'wiskContactPosition'});
mice = unique({flat.mouse});
xy = [[flat.wiskContactPosition]', [flat.trialVel]'];

corrs = nan(1, length(mice));
for i = 1:length(mice)
    bins = strcmp({flat.mouse}, mice{i});
    xy_sub = rmmissing(xy(bins,:));
    corrs(i) = corr(xy_sub(:,1), xy_sub(:,2));
end

% [h, p] = kstest(zscore(mat));
% fprintf('%i, %.5f: lagging paw landing position variability\n', h, p)

%% compute number of mice, sessions, trials, and frames in first paper

experimentSessions = getAllExperimentSessions();
sessions = experimentSessions.session;
[trials, frames] = deal(nan(1,length(sessions)));
for i = 1:length(sessions)
    disp(i/length(sessions))
    
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), 'obsOnTimes')
    vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runWisk.mp4'));
    
    frames(i) = round(vid.Duration*vid.FrameRate);
    trials(i) = length(obsOnTimes);
end

fprintf('\n----- paper stats -----\n')
fprintf('mice: %i\n', length(unique(experimentSessions.mouse)));
fprintf('sessions: %i\n', height(experimentSessions));
fprintf('trials: %i\n', nansum(trials));
fprintf('frames: %i\n', nansum(frames));


%% (temp) light on vs. light off speeds

figure('position', [200 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
matWiskContact = getDvMatrix(data, 'velAtWiskContact', vars.isLightOn, {'mouse'});
matObsOn = getDvMatrix(data, 'velAtObsOn', vars.isLightOn, {'mouse'});
mat = permute(cat(3, matObsOn, matWiskContact), [1 3 2]);

figure('position', [200 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
barFancy(matWiskContact, 'levelNames', {vars.isLightOn.levelNames}, 'ylabel', 'velocity at contact (m/s)', ...
    'colors', hsv(2), barProperties{:}, extraBarProps{:}, ...
    'comparisons', [1 2], 'test', 'ttest')

figure('position', [600 472.00 382.00 328.00], 'color', 'white', 'menubar', 'none');
barFancy(mat, 'levelNames', {vars.isLightOn.levelNames}, 'ylabel', 'velocity (m/s)', ...
    'colors', repelem(hsv(2),2,1), barProperties{:}, extraBarProps{:}, ...
    'comparisons', [1 2; 3 4], 'test', 'ttest')
% set(gca, 'YTick', 0:20:80, 'TickDir', 'out')
% saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', 'baselineDistAtLand'), 'svg');







