%% play around with functional clustering across cells

% create (and write to disk) clustering data struct
[data, cellInfo] = getClusteringData();

% put nested importances in a matrix
groups = fieldnames(data);
importance = nan(height(cellInfo), length(groups));
for i = 1:length(groups)
    importance(:,i) = data.(groups{i}).tbl.importance;
end
ngroups = length(groups);

%% plot stuff on ccf
% ccf = loadCCF();

% settings
clims = [0 .05];
cols = 2;

% inits
close all; figure('color', 'white', 'position', [2.00 722.00 1278.00 634.00])
colors = lines(ngroups);
rows = ceil(ngroups / cols);

for i = 1:ngroups
    subplot(rows, cols, i)
    title(groups{i})
    imp = data.(groups{i}).tbl.importance;
    plotLabels2D(ccf.labels, 'dim', 'dv', 'patchArgs', {'FaceColor', 'none', 'EdgeColor', 'black'}, ...
        'apGrid', ccf.ap, 'mlGrid', ccf.ml, 'dvGrid', ccf.dv)
    set(gca, 'xtick', [], 'ytick', [])

    cmap = customcolormap([0 1], [colors(i,:); 0 0 0]);
    scatter(cellInfo.ccf(:,1), cellInfo.ccf(:,2), 30, imp, 'filled', 'MarkerFaceAlpha', .8);
    colormap(gca, cmap)
%     colorbar
    caxis(clims)
end
savefig('E:\lab_files\paper2\plots\clustering\ccf_importance.fig')

%% group pair scatters

close all; figure('color', 'white', 'position', [2.00 247.00 1278.00 1109.00])
pairs = nchoosek(1:length(groups), 2);

for i = 1:ngroups
    for j = 1:ngroups
        if j<i
            ind = sub2ind([ngroups ngroups]-1, j, i-1);
            subplot(ngroups-1, ngroups-1, ind)
            x = importance(:,j);
            y = importance(:,i);
            
            scatter(x, y, 20, 'black', 'filled', 'MarkerFaceAlpha', .4)
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
savefig('E:\lab_files\paper2\plots\clustering\importance_scatters.fig')

%% percentages by nucleus

imp = nan(3, length(groups), height(cellInfo));  % fast, int, dent

names = {'fastigial', 'interpositus', 'dentate'};
for i = 1:3
    bins = strcmp(cellInfo.target, names{i});
    imp(i,:,bins) = importance(bins,:)';
end
imp = permute(imp, [2 1 3]);

close all; figure('color', 'white', 'position', [2.00 722.00 1278.00 634.00])
barFancy(imp, 'showBars', true, 'levelNames', {groups, names}, ...
    'colors', repmat(lines(length(groups)),3,1), 'showScatter', true, 'YLim', [0 .25], ...
    'showViolins', false, 'showErrorBars', false, 'scatterCondColor', true)
savefig('E:\lab_files\paper2\plots\clustering\nucleus_importance.fig')


%% clustering heatmaps and bar plots

% gmm clustering
maxGroups = 10;
aic = nan(1, maxGroups);
bic = nan(1, maxGroups);
gmm = cell(1, maxGroups);
for i = 2:maxGroups
    gmm{i} = fitgmdist(importance, i, ...
        'CovarianceType', 'full', 'SharedCovariance', false, ...
        'RegularizationValue', 1e-6, 'Replicates', 40, 'Options', statset('MaxIter', 1000));
    aic(i) = gmm{i}.AIC;
    bic(i) = gmm{i}.BIC;
end

% select best number of groups
[~, nclusters] = min(bic);
[groupid, ~, ~, ~, maha] = gmm{nclusters}.cluster(importance);
colors = lines(nclusters);


% kmeans clustering
% nclusters = 6;
% groupid = kmeans(importance, nclusters, 'distance', 'cosine');

close all; figure('color', 'white', 'position', [2.00 2.00 1278.00 1354.00])

% unsorted
subplot(3,2,1)
bins = ~any(isnan(importance),2);
imagesc(-importance(bins,:)); colormap gray
set(gca, 'box', 'off', 'tickdir', 'out', 'xtick', 1:ngroups, ...
    'XTickLabel', groups, 'XTickLabelRotation', 20, 'ytick', [])
title('unclustered')

subplot(3,2,3); hold on
plt = plot([bic; aic]', 'LineWidth', 2);
legend(plt, 'BIC', 'AIC', 'AutoUpdate','off')
set(gca, 'box', 'off', 'xlim', [1 maxGroups], 'tickdir', 'out')
plot(nclusters*[1 1], ylim, 'color', [1 1 1]*.2)
xlabel('clusters'); ylabel('BIC score')

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

% bar plots
subplot(3,2,[5 6])
imp = nan(nclusters, ngroups, height(cellInfo));
for i = 1:nclusters
    bins = groupid==i;
    imp(i,:,bins) = importance(bins,:)';
end

clusternames = cellstr(num2str((1:nclusters).', 'clusters %i'))';
barFancy(imp, 'showBars', true, 'levelNames', {clusternames, groups}, ...
    'colors', repelem(lines(nclusters),ngroups,1), 'showScatter', true, 'YLim', [0 .25], ...
    'showViolins', false, 'showErrorBars', false, 'scatterCondColor', true)

savefig('E:\lab_files\paper2\plots\clustering\clustering.fig')

%% plot responses sorted by importance

close all; figure('color', 'white', 'position', [2.00 2.00 1865.00 1354.00])


% inits
histoBins = 10;
traceNum = 10;
colors = lines(length(groups));
cols = length(groups);
cmap = customcolormap([0 .5 1], [1 .2 .2; 1 1 1; .2 .2 1]);

for i = 1:length(groups)
    
    % group inits
    tbl = data.(groups{i}).tbl;
    x = data.(groups{i}).x;
    predictor = data.(groups{i}).predictor;
    [~, sortInds] = sort(tbl.importance);
    finalInd = find(~isnan(tbl.importance(sortInds)), 1, 'last');

    % histogram
    subplot(5,cols,i); hold on
    [~, clusterCenters] = kmeans(tbl.importance(~isnan(tbl.importance)), 2);
    thresh = mean(clusterCenters);
    histogram(tbl.importance, histoBins, 'EdgeColor', 'none', 'FaceColor', colors(i,:))
    xlabel('deviance explained')
    set(gca, 'box', 'off')
    plot([thresh thresh], ylim, 'color', [1 1 1]*.2)
    title(sprintf('%s (%.1f%%)', groups{i}, mean(tbl.importance>thresh)*100))
    
    % responses
    subplot(5,cols,i+cols); hold on
    plot(x, tbl.response(sortInds(1:traceNum),:)', 'color', [0 0 0 .4])
    plot(x, tbl.response(sortInds(finalInd-traceNum:finalInd),:)', 'color', [colors(i,:) .6], 'LineWidth', 2)
    ylabel('firing rate (z score)')
    xlabel(predictor, 'Interpreter', 'none')
    set(gca, 'box', 'off', 'XLim', [x(1) x(end)])

    % heatmap
    subplot(5,cols,i+[2 4]*cols); hold on
    imagesc(x, 1:finalInd, tbl.response(sortInds(1:finalInd), :));
    colormap(cmap)
    
    y = [find(tbl.importance(sortInds)>thresh, 1, 'first'), ...
         find(tbl.importance(sortInds)>thresh, 1, 'last')] + [-.5 .5];
    plot([x(1) x(1)], y, 'Color', colors(i,:), 'LineWidth', 5)
    
    xlabel(predictor, 'Interpreter', 'none')
    set(gca, 'box', 'off', 'ytick', [], 'TickDir', 'out', 'ydir', 'normal', 'ylim', [1,finalInd], 'xlim', [x(1) x(end)])
end

savefig('E:\lab_files\paper2\plots\clustering\group_responses.fig')



