%% CLUSTA BABY, CLUSTA!

data = getUnitInfo;
paper2_config;

%% collect tuning metrics for each units

groups = {'grossKinematics', 'lick', 'obstacle', 'pawKinematics', 'reward', 'vision', 'whiskers'};
metricNames = {'residual', 'upper', 'lower'};  % metrics for the tuning of a unit to a predictor group
nucleiNames = {'dentate', 'interpositus', 'fastigial'};
ngroups = length(groups);
for i = 1:length(metricNames)
    data.(metricNames{i}) = nan(height(data), length(groups));
end


for i = 1:height(data)
    disp(i / height(data))
    
    % load up dem models
    try
        resModels = load(fullfile('E:\lab_files\paper2\modelling\glms\residual_glms', ...
            sprintf('%s_cell_%i_glm.mat', data.session{i}, data.unit(i))));
        upperLowerModels = load(fullfile('E:\lab_files\paper2\modelling\glms\upperlower_glms', ...
            sprintf('%s_cell_%i_glm.mat', data.session{i}, data.unit(i))));

        data{i, 'residual'} = resModels.models{groups, 'dev'}';
        data{i, 'upper'} = upperLowerModels.models{groups, 'dev_in'}';
        data{i, 'lower'} = upperLowerModels.models{groups, 'dev_out'}';
    catch
        fprintf('ERROR WITH %s unit %i\n', data.session{i}, data.unit(i))
    end    
end

for i = 1:length(metricNames)
    data.(metricNames{i})(data.(metricNames{i})<0) = 0;
end


%% overall upper lower bars

close all
figure('color', 'white', 'menubar', 'none', 'position', [227.00 396.00 306.00 294.00]); hold on

halfWid = .3;  % half width of bar

for i = 1:ngroups
    lower = data.lower(:,i);
    upper = data.upper(:,i);
    low = nanmean(lower);
    upp = nanmean(upper);
    lowSem = nanstd(lower) / sqrt(height(data));
    uppSem = nanstd(upper) / sqrt(height(data));
    
    patch([-1 1 1 -1]*halfWid + i, [low low upp upp], ...
        cfg.groupColors(i,:), 'EdgeColor', 'none'); hold on
    
    plot([i i], upp+[-1 1]*uppSem, 'color', 'black', 'LineWidth', 2)
    plot([i i], low+[-1 1]*lowSem, 'color', 'black', 'LineWidth', 2)
end

set(gca, 'XTick', 1:ngroups, 'XTickLabel', groups, 'XTickLabelRotation', 45);
ylabel({'deviance explained', '(lower and upper bounds)'})

saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab\upper_lower_bars.svg')


%% scatters to explore correlations between different tuning metrics

scatPlots = {{'upper', 'lower'}, {'upper', 'residual'}, {'residual', 'lower'}};
close all
figure('color', 'white', 'menubar', 'none', 'position', [496.00 111.00 578.00 1231.00]);

for i = 1:length(groups)
    for j = 1:length(scatPlots)
        subplot(length(groups), length(scatPlots), (i-1)*length(scatPlots)+j)
        
        xMetric = scatPlots{j}{1};
        yMetric = scatPlots{j}{2};
        
        x = data{:, xMetric}(:,i);
        y = data{:, yMetric}(:,i);
        x(x<0) = 0; y(y<0) = 0;
        x = x(~isnan(x)); y = y(~isnan(y));
        
        scatter(x, y, 10, [0 0 0], 'filled', ...
            'MarkerFaceAlpha', .2, 'MarkerEdgeColor', 'none')
        xlabel(xMetric)
        ylabel(yMetric)
        title(groups{i})
        
        xlims = xlim; ylims = ylim;
        text(xlims(2), ylims(2), ...
            sprintf('r = %.2f', corr(x, y)), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
        set(gca, 'xlim', xlims, 'ylim', ylims)
    end 
end

saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab\tuning_metric_correlations.svg')


%% distros for upper and lower per group

% settings
xlims = [0 .2; 0 .8];  % lower lims; upper lims
ylims = [0 height(data)];
nbins = 20;

close all
figure('color', 'white', 'menubar', 'none', 'position', [603.00 362.00 368.00 918.00])
ngroups = length(groups);
colors = cfg.upperLowerColors;

histArgs = {'EdgeColor', 'none'};
colLabels = {'deviance lower bound', 'deviance upper bound'};

for i = 1:ngroups
    upperAndLower = {data.lower(:, i), data.upper(:, i)};
    for j = 1:2
        subplot(ngroups, 2, (i-1)*2+j); hold on
        d = upperAndLower{j};

        % histos
        binEdges = linspace(xlims(j,1), xlims(j,2), nbins+1);
        histLow = histogram(d, binEdges, 'FaceColor', colors(j,:), histArgs{:});
        set(gca, 'XLim', xlims(j,:), cfg.axArgs{:})
        set(gca, cfg.axArgs{:})

        % add medians
        plot([1 1]*nanmedian(d), ylim, 'LineWidth', 2, 'color', colors(j,:))

        % fancify
        if j==1; ylabel(groups{i}); end
        if i<ngroups; set(gca, 'XTickLabel', []); end
        if i==ngroups; xlabel(colLabels{j}); end
        limitticks(true)
    end
end

saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab\histos_upper_lower.svg')


%% scatter importance on ccf for each group

% settings
metric = 'lower';  % residual, upper, lower
percentileLims = [10 90];

% inits
close all
figure('color', 'white', 'menubar', 'none', 'position', [45.00 926.00 2325.00 414.00])
nGroups = length(groups);
subplotDims = repmat([2 nGroups], 2, 1);

for i = 1:nGroups
    importance = data.(metric)(:, i);
    clims = prctile(importance, percentileLims);
    colors = interp1(clims, [1 1 1; cfg.groupColors(i,:)], importance, 'nearest', 'extrap');
    subplotArgs = [subplotDims, [i; i+nGroups]];
    plotUnitsOnCcf(data, 'colors', colors, ...
        'scatArgs', {'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceAlpha', .6}, ...
        'subplots', subplotArgs, 'scatSz', 20);
    subplot(2, nGroups, i);
    title(groups{i})
    set(gca, 'visible', 'on', 'xcolor', 'none', 'ycolor', 'none')
end

saveas(gcf, ['E:\lab_files\paper2\paper_figures\matlab\importance_ccf_' metric '.svg'])

%% tuning by nucleus bars

% settings
metric = 'lower';

close all;
figure('color', 'white', 'position', [448.00 633.00 328.00 306.00], 'menubar', 'none'); hold on

x = (1 : (length(groups)*3)) + ...
    repelem(0:length(groups)-1, 3) * .5;  % offset each group
patches = nan(1, 3);

for i = 1:length(groups)
    
    overallMed = nanmedian(data.(metric)(:, i));  % median across nuclei
    xvals = x((i-1)*3+2) + [-1.5 1.5];
    medianLn = plot(xvals, overallMed+[0 0], 'Color', [0 0 0 .3], 'LineWidth', 1.5);
    
    for j = 1:3
        xind = (i-1) * 3 + j;
        bins = strcmp(data.nucleus, nucleiNames{j});
        d = data.(metric)(bins, i);
        c = cfg.nucleusColors(j,:);
        
        % compute stats
        mn = nanmedian(d);
        blims = prctile(d, [25 75]);  % top and bottom of box
        iqr = diff(blims);  % interquartile range
        outlims = [blims(1)-iqr*1.5 blims(2)+iqr*1.5];  % data outside this range are considered outliers (https://www.mathworks.com/help/stats/boxplot.html)
        
        % whiskers extend to most extreme non-outlier value
        validBins = d>outlims(1) & d<outlims(2);
        whiskerLims = [min(d(validBins)) max(d(validBins))];
        
        patches(j) = patch(x(xind)+[-.5 .5 .5 -.5], [blims(1) blims(1) blims(2) blims(2)], c, ...
            'EdgeColor', 'none', 'FaceAlpha', .8);  % box
        plot(x(xind)+[-.5 .5], [mn mn], 'LineWidth', 1.5, 'Color', [.15 .15 .15])  % median
        plot(x(xind)*[1 1], [blims(1) whiskerLims(1)], 'color', [.15 .15 .15])  % bottom whisker
        plot(x(xind)*[1 1], [blims(2) whiskerLims(2)], 'color', [.15 .15 .15])  % top whisker
    end
    uistack(medianLn, 'top')
end

set(gca, 'xtick', x(2:3:length(x)), 'XTickLabel', groups, 'XTickLabelRotation', 45, cfg.axArgs{:})
ylabel('deviance explained')
limitticks(true)
legend(patches, nucleiNames, 'box', 'off')

saveas(gcf, ['E:\lab_files\paper2\paper_figures\matlab\deviance_bars_' metric '.svg'])

%% scatter importance for all group pairs

% settings
metric = 'lower';  % residual, upper, lower

% inits
close all; figure('color', 'white', 'position', [212.00 475.00 1019.00 814.00])
ngroups = length(groups);

for i = 1:ngroups
    for j = 1:ngroups
        if j<i
            ind = sub2ind([ngroups ngroups]-1, j, i-1);
            subplot(ngroups-1, ngroups-1, ind)
            x = data.(metric)(:, strcmp(groups, groups{j}));
            y = data.(metric)(:, strcmp(groups, groups{i}));
            
            scatter(x(bins), y(bins), 20, 'black', 'filled', 'MarkerFaceAlpha', .4)
            xlabel(groups{j})
            ylabel(groups{i})
%             set(gca, 'xlim', [0 .2], 'ylim', [0 .2])
            
            bins = ~isnan(x) & ~isnan(y);
            xycorr = corr(x(bins), y(bins));
            xlims = xlim; ylims = ylim;
            text(xlims(2), ylims(2), sprintf('r=%.2f', xycorr), ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
        end
    end
end

saveas(gcf, ['E:\lab_files\paper2\paper_figures\matlab\importance_scatters_' metric '.svg'])

%% cluster! (run these cell serially)

%% clustering heatmaps and bar plots

% settings
rng(1)
metric = 'lower';  % residual, upper, lower
groupsToExclude = {'vision', 'obstacle'};
nclusters = [];

% gmm clustering
close all
groupBins = ~ismember(groups, groupsToExclude);
importance = data.(metric)(:, groupBins);
% isModulated = any(importance>thresh, 2);

% importance = importance ./ nanstd(importance, 1);  % scale each group by variance of deviance

bins = ~any(isnan(importance), 2);
groupid = nan(height(data), 1);
[groupid(bins), ~, ~, bic] = clusterResponses(importance(bins,:), ...
    'pcs', 0, 'plot', true, 'nclusters', nclusters, ...
    'clusterMetric', 'bic', 'maxclusters', 8);
set(gca, 'XTickLabel', groups(groupBins), 'XTickLabelRotation', 45)
nclusters = length(unique(groupid(~isnan(groupid))));
colors = [lines(nclusters); 0 0 0];  % black for no group
groupid(isnan(groupid)) = nclusters+1;  % make new group for ungrouped units
disp(['n clusters: ' num2str(nclusters)])

%% umap embeddings
scatColors = colors(groupid,:);
embedding = nan(height(data), 2);
embedding(bins, :) = run_umap(importance(bins, :), ...
    'sgd_tasks', 1, 'verbose', 'text');

figure('color', 'white', 'menubar', 'none', 'position', [482.00 431.00 674.00 288.00]);
subplot(1,2,2); hold on
scatter(embedding(:, 1), embedding(:, 2), 10, scatColors, 'filled', ...
    'MarkerFaceAlpha', .6, 'MarkerEdgeColor', 'none');
set(gca, 'XColor', 'none', 'YColor', 'none', cfg.axArgs{:})
title('UMAP Embeddings')

% scatter by nucleus
subplot(1,2,1); hold on
scats = nan(1, 3);
for i = 1:length(nucleiNames)
    nucleusBins = strcmp(data.nucleus, nucleiNames{i});
    scats(i) = scatter(embedding(nucleusBins, 1), embedding(nucleusBins, 2), 10, 'filled', ...
        'MarkerFaceColor', cfg.nucleusColors(i,:), 'MarkerFaceAlpha', .6, 'MarkerEdgeColor', 'none');

end

set(gca, 'XColor', 'none', 'YColor', 'none', cfg.axArgs{:})
title('UMAP Embeddings')
legend(scats, nucleiNames, 'location', 'best')

saveas(gcf, ['E:\lab_files\paper2\paper_figures\matlab\tuning_umap_' metric '.svg'])

%% bic and heatmaps

% settings
percentileLim = 99.5;  % percentile


% inits
% close all;
figure('color', 'white', 'position', [778.00 567.00 245.00 650.00], 'menubar', 'none');

subplot(5, 1, 2:5); hold on
tuning = data.(metric)(:, groupBins);
clims = [0 prctile(tuning(:), percentileLim)];
% [~, sortInds] = sort(groupid);

% for each cluster, sort by most active channel within the group
temp = ones(size(groupid));
for i = 1:nclusters
    bins = (groupid==i);
    isBigGroup = sum(bins)>(height(data)/2);  % don't sort within the big group, which is heterogeneous
    if ~isBigGroup
        [~, maxInd] = max(nanmean(tuning(bins,:),1));  % ind of most active group for this cluster
        temp(groupid==i) = tuning(bins, maxInd);
    end
end
[~, sortInds] = sortrows([groupid temp]);


% heatmap
imagesc(tuning(sortInds, :), clims); colormap gray;

% add lines on the side for groups
for i = 1:nclusters
    ystart = find(groupid(sortInds)==i, 1, 'first');
    yend = find(groupid(sortInds)==i, 1, 'last');
    plot([.5 .5], [ystart-.5 yend+.5], 'Color', colors(i,:), 'LineWidth', 5)
end

set(gca, 'ylim', [.5 sum(groupid<=nclusters)+.5], 'YTick', [], ...
    'XTick', 1:sum(groupBins), 'XTickLabel', groups(groupBins), ...
    'XTickLabelRotation', 45, 'XLim', [.5 sum(groupBins)+.5], cfg.axArgs{:})

% bic plot
subplot(5, 1, 1)
plot(1:length(bic), bic, 'color', [.2 .2 .2], 'LineWidth', 3)
ylabel('BIC score')
xlabel('number of clusters')
set(gca, 'XLim', [1 length(bic)], 'box', 'off', cfg.axArgs{:})
limitticks

saveas(gcf, ['E:\lab_files\paper2\paper_figures\matlab\cluster_heatmaps_' metric '.svg'])


%% example units

sessions = {'200116_000', '200623_001', '191007_003', '201014_000', ...
            '201217_000', '201214_000', '201209_000', '201208_000', ...
            '201119_000', '201117_000', '201116_000', '201031_000'};
units = [30  270 3   193 ...
         76  8   198 191 ...
         203 130 131 104];
showDeviance = false;  % whether to include show plot with deviance explained for each group

% sessions = sessions(1:4);  % temp

close all
figure('color', 'white', 'position', [85.00 8.00 1116.00 1331.00], 'menubar', 'none');
for i = 1:length(sessions)
    unitbin = strcmp(data.session, sessions{i}) & data.unit==units(i);
    if showDeviance == true
        groupDevExplained = [data.lower(unitbin,:); data.upper(unitbin,:)];
    else
        groupDevExplained = NaN;
        groups = NaN;
    end
    
    if i==1; hideTraceText=false; else; hideTraceText=true; end
    if i==1; hideTitles=false; else; hideTitles=true; end
    if i==length(sessions); hideX=false; else; hideX=true; end
    
    plotExampleUnit(sessions{i}, units(i), groups, groupDevExplained, ...
        'showDeviance', showDeviance, 'subplotCols', length(sessions), 'subplotColIdx', i, ...
        'hideTraceText', hideTraceText, 'hideTitles', hideTitles, 'hideX', hideX);
end

set(gcf, 'Renderer', 'painters')  % ensures export is not rasterized
saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab\unit_examples.svg')

%% mixed selectivity states

thresh = 0.02;
nTunedGroups = sum(data.lower > thresh, 2);
mean(nTunedGroups > 1)  % fraction units tuned to >1 group


%% plot PSTHs for all units

rows = 9;


data_ = sortrows(data, "nucleus");
sessions = data_.session;
units = data_.unit;
nuclei = data_.nucleus;

% sessions = sessions(1:4);  % temp

close all
for i = 1:length(sessions)
    
    row = mod(i - 1, rows) + 1;
    
    % make new figure
    if row == 1
        figure('color', 'white', 'position', [85.00 8.00 1116.00 1331.00], 'menubar', 'none');
    end
    
    unitbin = strcmp(data.session, sessions{i}) & data.unit==units(i);
    
    if row==1; hideTraceText=false; else; hideTraceText=true; end
    if row==1; hideTitles=false; else; hideTitles=true; end
    if row==length(sessions); hideX=false; else; hideX=true; end
    
    title = sprintf('%s: %s (%i)', nuclei{i}, sessions{i}, units(i));
    plotExampleUnit(sessions{i}, units(i), nan, nan, ...
        'showDeviance', false, 'subplotRows', rows, 'subplotColIdx', row, ...
        'hideTraceText', hideTraceText, 'hideTitles', hideTitles, 'hideX', hideX, ...
        'title', title);
    
    % save previous figure
    if row == rows || (i == length(sessions))
        n = ceil(i / rows);
        saveas(gcf, sprintf('E:\\lab_files\\paper2\\plots\\all_unit_psths\\page_%i.png', n))
        close(gcf)
    end
    
end















