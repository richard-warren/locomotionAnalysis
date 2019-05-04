% load experiment data
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'whiskerTrim_data.mat'), 'data'); disp('whiskerTrim data loaded!')

% global settings
miceToExclude = {'sen1'};
colorWisk = [51 204 255]/255; %[255 204 51];


% initializations
data.data = data.data(~ismember({data.data.mouse}, miceToExclude));
% vars.condition = struct('name', 'condition', ...
%     'levels', {{'bilateral full', 'unilateral full', 'unilateral int1', 'unilateral int2', 'unilateral int3', 'unilateral deltaOnly'}}, ...
%     'levelNames', {{'biFull', 'unFull', 'un1', 'un2', 'un3', 'delta'}});
vars.condition = struct('name', 'condition', ...
    'levels', {{'bilateral full', 'unilateral int3', 'unilateral deltaOnly'}}, ...
    'levelNames', {{'biFull', 'uniPartial', 'deltaOnly'}});
colors = repmat(colorWisk, length(vars.condition.levels), 1) .* fliplr(linspace(0,1,length(vars.condition.levels)))';

%% HEIGHT SHAPING




% settings
xLims = [3 10];
isHgtPreObs = true; % do you measure the height at peak (false) or before the paw reaches the obstacle (true)

% initializations
flat = flattenData(data, {'mouse', 'session', 'trial', 'isLightOn', ...
    'obsHgt', 'preObsHgt', 'isFore', 'isLeading', 'stepOverMaxHgt', 'condition'});
flat = flat(~isnan([flat.obsHgt]) & ...
            ~isnan([flat.stepOverMaxHgt]) & ...
            ~isnan([flat.preObsHgt]) & ...
            [flat.isLeading] & ...
            [flat.isFore]); % add conditionals here
[~, conditions] = ismember({flat.condition}, vars.condition.levels);


% SCATTERS BINNED BY ANIMAL, SPEED

% settings
binNum = 40;
scatAlpha = .1;
scatSize = 50;

% initializations
mice = unique({flat.mouse});
binEdges = linspace(xLims(1), xLims(2), binNum+1);
binCenters = binEdges(1:end-1) + .5*diff(binEdges(1:2));
binIds = discretize([flat.obsHgt]*1000, binEdges);

% collect data
binnedHgts = nan(max(conditions), length(mice), binNum); % contains the median paw height for each conditino, mouse, and paw height
for c = 1:max(conditions)
    disp(c)
    for m = 1:length(mice)
        for b = 1:binNum
            bins = (conditions==c) & strcmp({flat.mouse}, mice{m}) & (binIds==b);
            binnedHgts(c,m,b) = nanmean([flat(bins).preObsHgt]);
        end
    end
end

%% plot condition data
figure('name', 'baseline', 'color', 'white', 'menubar', 'none', 'position', [2000 50 500 400])
plot([0 xLims(2)], [0 xLims(2)], 'Color', [1 1 1]*.6, 'LineWidth', 3) % add unity line
xlabel('obstacle height (mm)')
ylabel('paw height (mm)')
scatterPlotRick(repelem(binCenters,length(mice)*max(conditions)), binnedHgts(:), repmat(1:max(conditions),1,binNum*length(mice)), ...
    {'colors', colors, 'maxScatterPoints', 5000, 'lineAlpha', 1, 'scatAlpha', .2, 'scatSize', 40, ...
    'conditionNames', vars.condition.levelNames});
set(gca, 'XLim', xLims)

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'whiskerTrim_heightShapingScatters');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');


% HEIGHT SHAPING BARS
[corrs, slopes] = deal(nan(max(conditions),length(mice))); % isFore(10) X isLeading(10) X mouse
foreSequence = [true false];
leadingSequence = [true false];

for i = 1:max(conditions)
    for k = 1:length(mice)
        bins = conditions==i & ...
               strcmp({flat.mouse}, mice{k});
        x = [flat(bins).obsHgt];
        y = [flat(bins).preObsHgt];
        fit = polyfit(x, y, 1);
        corrs(i,k) = corr(x', y');
        slopes(i,k) = fit(1);
    end
end


figure('position', [2000 400 600 200], 'color', 'white', 'menubar', 'none');

% correlations
subplot(1,2,1)
barPlotRick(corrs, {'conditionNames', {vars.condition.levelNames}, 'ylabel', 'paw/obstacle height correlation', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colors, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [-.2 .6], 'ytick', -.2:.4:.6, ...
    'compareToFirstOnly', true})
line(get(gca, 'XLim'), [0 0], 'color', [.4 .4 .4]) % add horizontal line at y=0

% slopes
subplot(1,2,2)
barPlotRick(slopes, {'conditionNames', {vars.condition.levelNames}, 'ylabel', 'paw/obstacle height slope', ...
    'showViolins', false, 'lineThickness', 2, 'conditionColors', colors, 'addBars', true, ...
    'violinAlpha', .1, 'scatColors', 'lines', 'scatAlpha', .3, 'showStats', false, 'ylim', [-.2 .8], 'ytick', -.2:.5:.8, ...
    'compareToFirstOnly', true})
line(get(gca, 'XLim'), [0 0], 'color', [.4 .4 .4]) % add horizontal line at y=0 % add line at zero

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'whiskerTrim_heightShapingBars');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');



% PAW HEIGHT BY OBS HGT RUNNING AVERAGE
figure('Color', 'white', 'Position', [2000 400 500 400], 'MenuBar', 'none');
plot([0 xLims(2)], [0 xLims(2)], 'Color', [1 1 1]*.6, 'LineWidth', 3) % add unity line
logPlotRick([flat.obsHgt]*1000, [flat.preObsHgt]*1000, ...
    {'colors', colors, 'conditions', conditions, 'xlabel', 'obstacle height', 'ylabel', 'paw height (mm)', 'plotMice', false, ...
    'xlim', [3.4 10], 'binWidth', 1, 'binNum', 100, 'smoothing', 1, 'lineWidth', 4, 'mouseNames', {flat.mouse}, ...
    'errorFcn', @(x) std(x)/sqrt(size(x,1))})
set(gca, 'xlim', [4 10])

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'whiskerTrim_heightShapingMovingAvgs');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');