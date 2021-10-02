%% CLUSTA BABY, CLUSTA!

data = getUnitInfo;
paper2_config;

%% collect tuning metrics for each units

groups = {'grossKinematics', 'lick', 'obstacle', 'pawKinematics', 'reward', 'vision', 'whiskers'};
metricNames = {'residual', 'upper', 'lower'};  % metrics for the tuning of a unit to a predictor group
for i = 1:length(metricNames)
    data.(metricNames{i}) = nan(height(data), length(groups));
end


for i = 1:height(data)
    disp(i)
    
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

%% umap embeddings

% TODO: color by discovered groups?

metric = 'upper';  % residual, upper, lower
groupsToExclude = {'vision', 'obstacle'};

groupBins = ~ismember(groups, groupsToExclude);
d = data.(metric)(:, groupBins);
bins = ~any(isnan(d), 2);
embedding = nan(height(data), 2);
embedding(bins, :) = run_umap(d(bins, :), 'sgd_tasks', 1, 'verbose', 'text');

close all;
figure('color', 'white', 'menubar', 'none', 'position', [881.00 820.00 306.00 288.00]); hold on

% scatter by nucleus
nuclei = {'dentate', 'interpositus', 'fastigial'};
scats = nan(1, 3);
for i = 1:length(nuclei)
    nucleusBins = strcmp(data.nucleus, nuclei{i});
    scats(i) = scatter(embedding(nucleusBins, 1), embedding(nucleusBins, 2), 10, 'filled', ...
        'MarkerFaceColor', cfg.nucleusColors(i,:), 'MarkerFaceAlpha', .6, 'MarkerEdgeColor', 'none');

end

set(gca, 'XColor', 'none', 'YColor', 'none', cfg.axArgs{:})
title('UMAP Embeddings')
legend(scats, nuclei, 'location', 'best')

saveas(gcf, ['E:\lab_files\paper2\paper_figures\matlab\umap_' metric '.svg'])


%% scatter importance on ccf for each group

% settings
metric = 'residual';  % residual, upper, lower
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

%% cluster!

%% clustering heatmaps and bar plots

% settings
metric = 'residual';  % residual, upper, lower
groupsToExclude = {'vision', 'obstacle'};
nclusters = [];

% gmm clustering
close all
groupBins = ~ismember(groups, groupsToExclude);
importance = data.(metric)(:, groupBins);
bins = ~any(isnan(importance), 2);
groupid = nan(height(data), 1);
groupid(bins) = clusterResponses(importance(bins,:), ...
    'pcs', 0, 'plot', true, 'nclusters', nclusters, 'clusterMetric', 'aic');
nclusters = length(unique(groupid));
colors = lines(nclusters);

%%

figure('color', 'white', 'position', [2.00 2.00 1278.00 1354.00])

% unsorted
subplot(3,2,[1 3])
clusterbins = ~any(isnan(importance),2);
imagesc(-importance(clusterbins,:)); colormap gray
set(gca, 'box', 'off', 'tickdir', 'out', 'xtick', 1:ngroups, ...
    'XTickLabel', groups, 'XTickLabelRotation', 20, 'ytick', [])
title('unclustered')

% sorted
subplot(3,2,[2 4]); hold on
title('clustered')
[groupidSorted, sortInds] = sortrows([groupid, mean(importance,2)]);

finalInd = find(~isnan(groupid(sortInds)), 1, 'last');
imagesc(-importance(sortInds(1:finalInd),:)); colormap gray
for i = 1:nclusters
    y = [find(groupidSorted==i,1,'first') find(groupidSorted==i,1,'last')] + [-.5 .5];
    plot([.5 .5], y, 'color', colors(i,:), 'LineWidth', 5)
end
set(gca, 'box', 'off', 'tickdir', 'out', 'xtick', 1:ngroups, ...
    'XTickLabel', groups, 'XTickLabelRotation', 20, ...
    'xlim', [0 ngroups]+.5, 'ylim', [0 finalInd], 'ytick', [])
% 
% % bar plots
% subplot(3,2,[5 6])
% imp = nan(nclusters, ngroups, height(cellInfoSub));
% for i = 1:nclusters
%     clusterbins = groupid==i;
%     imp(i,:,clusterbins) = importance(clusterbins,:)';
% end
% 
% clusternames = cellstr(num2str((1:nclusters).', 'clusters %i'))';
% barFancy(imp, 'showBars', true, 'levelNames', {clusternames, groups}, ...
%     'colors', repelem(lines(nclusters),ngroups,1), 'showScatter', true, 'YLim', [0 .3], ...
%     'showViolins', false, 'showErrorBars', false, 'scatterCondColor', true)
% 
% savefig('E:\lab_files\paper2\plots\clustering\clustering.fig')
% 
% % custering on ccf (run cell above first)
% 
% % settings
% skipmode = true;
% 
% % inits
% figure('color', 'white', 'position', [592.00 372.00 628.00 818.00])
% colors = lines(nclusters);
% views = {'ap', 'dv'}; dims = [1 3; 1 2];
% modalgroup = mode(groupid);
% 
% for i = 1:2
%     subplot(2,1,i); hold on
%     plotLabels2D(ccf.labels, 'dim', views{i}, 'patchArgs', {'FaceColor', 'none', 'EdgeColor', 'black'}, ...
%         'apGrid', ccf.ap, 'mlGrid', ccf.ml, 'dvGrid', ccf.dv)
%     set(gca, 'xtick', [], 'ytick', [])
%     
%     for j = 1:nclusters
%         bins = groupid==j;
%         if ~skipmode || j~=modalgroup
%             scatter(cellInfoSub.ccfMm(bins, dims(i,1)), cellInfoSub.ccfMm(bins, dims(i,2)), [], colors(j,:), 'filled', ...
%                 'MarkerFaceAlpha', .8, 'MarkerEdgeColor', 'none');
%         end
%     end
% end
% savefig('E:\lab_files\paper2\plots\clustering\ccf_groups.fig')






