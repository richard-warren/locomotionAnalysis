%% load experiment data
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')





%% ----------
% PLOT THINGS
%  ----------

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
flat = flat(~[flat.isWheelBreak]);
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
colNames = {'hind limb', 'fore limb'};
rowNames = {'lagging', 'leading'};
isLeading = [false false true true];
isFore = [false true false true];
ctlColor = [.5 .5 .5];

obsHgtBins = 4; % discretize obstacle heights into this many bins
xLims = [-.05 .05];
yLims = [0 .016];
colors = hsv(4);  % last term is control step color
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

for i = 1:4
    subplot(2,2,i)
    plotColor = repmat(colors(i,:), obsHgtBins, 1) .* linspace(fading,1,obsHgtBins)'; % create color matrix fading from colors(i,:) -> colors(i,:)*fading
    bins = [flat.isLeading]==isLeading(i) & ...
           [flat.isFore]==isFore(i);

    % plot control step
    plotKinematics(kinDataCtl(bins,[1,3],:), [flat(bins).obsHgt], ones(1,sum(bins)), ...
        {'colors', ctlColor, 'obsAlpha', 0, 'lineAlpha', .8}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    
    plotKinematics(kinData(bins,[1,3],:), [flat(bins).obsHgt], conditions(bins), ...
        {'colors', plotColor, 'obsAlpha', 1, 'lineAlpha', .8, 'mouseNames', {flat(bins).mouse}}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
    set(gca, 'XLim', xLims, 'YLim', yLims)
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
[~, conditions] = ismember(b, a, 'rows');

plotKinematics(kinDataTop, [flat.obsHgt], conditions', ...
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
plotKinematics(kinDataBot, [flat.obsHgt], conditions', ...
    {'colors', colors, 'obsAlpha', .25, 'lineAlpha', 1, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @nanstd, 'lineWidth', 3, 'isBotView', true}) % if 'mouseNames' is provided, plotKinematics avgs within, then across mice for each condition
set(gca, 'XLim', xLims, 'YLim', yLimsBot)

%save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineKinematics_overlayBot');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


%% HEIGHT SHAPING SCATTER


% SCATTER ACROSS ALL ANIMALS

% settings
colors = hsv(4);
xLims = [3 10];
yLims = [3 20];
isLeading = [false false true true];
isFore = [false true false true];
isHgtPreObs = true; % do you measure the height at peak (false) or before the paw reaches the obstacle (true)

% obs hgt vs paw hgt (manip, ipsi/contra, leading/lagging, fore/hind)
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isLeading', 'stepOverMaxHgt'});
flat = flat([flat.isLightOn]); % add conditionals here

% get condition numbers, where each condition is unique combination of isLeading and isFore
a = [isLeading; isFore]'; % each row is a condition (00,01,10,11)
b = [[flat.isLeading]; [flat.isFore]]';
[~, conditions] = ismember(b, a, 'rows');

if isHgtPreObs; pawHgts = [flat.preObsHgt]*1000; else; pawHgts = [flat.stepOverMaxHgt]*1000; end
obsHgts = [flat.obsHgt]*1000;
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 50 500 400])
[corrs, slopes] = scatterPlotRick(obsHgts, pawHgts, conditions, ...
    {'colors', colors, 'maxScatterPoints', 5000, 'lineAlpha', 1, 'scatAlpha', .1, 'scatSize', 40});
set(gca, 'XLim', xLims, 'YLim', yLims)

% add unity line
plot([0 xLims(2)], [0 xLims(2)], 'Color', [1 1 1]*.8, 'LineWidth', 2) % add unity line
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
connectMouseDots = false;

% initializations
mice = unique({flat.mouse});
binEdges = linspace(xLims(1), xLims(2), binNum+1);
binCenters = binEdges(1:end-1) + .5*diff(binEdges(1:2));
bins = discretize(obsHgts, binEdges);

% collect data
meanHgts = nan(4, length(mice), binNum); % contains the median paw height for each conditino, mouse, and paw height
for c = 1:4
    for m = 1:length(mice)
        for b = 1:binNum
            meanHgts(c,m,b) = nanmean(pawHgts((conditions==c)' & strcmp({flat.mouse}, mice{m}) & bins==b));
        end
    end
end

% plot condition data
figure('position', [2000 400 500 400], 'color', 'white', 'menubar', 'none');
lines = nan(1,5);
lines(5) = plot([0 xLims(2)], [0 xLims(2)], 'Color', [1 1 1]*.6, 'LineWidth', 3); hold on; % add unity line;
corrs = nan(1,4);

for c = 1:4
    cData = squeeze(meanHgts(c,:,:));
    x = repelem(binCenters,length(mice));
    y = cData(:)';
    validBins = ~isnan(x) & ~isnan(y);
    
    if connectMouseDots % add lines connected mouse across bins
        for m = 1:length(mice); plot(binCenters, cData(m,:), 'LineWidth', 1, 'Color', repelem(1-scatAlpha,1,3)); end
    end
    scatter(x(validBins), y(validBins), scatSize, colors(c,:), 'filled', ...
        'MarkerEdgeAlpha', scatAlpha, 'MarkerFaceAlpha', scatAlpha); hold on;
    fit = polyfit(x(validBins), y(validBins), 1);
    corrs(c) = corr(x(validBins)', y(validBins)');
    lines(c) = plot(binEdges, polyval(fit, binEdges), 'linewidth', 4, 'color', colors(c,:));
end

uistack(lines, 'top')
set(gca, 'box', 'off', 'XLim', xLims)
xlabel('obstacle height (mm)')
ylabel('paw height (mm)')

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baselineHeightShapingScatters_perMouse');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');














