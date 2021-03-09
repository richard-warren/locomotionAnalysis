%% play around with functional clustering across cells

% create (and write to disk) clustering data struct
aggregateResponses();  % only need to do once
[data, cellInfo] = getClusteringData();

% put nested importances in a matrix
groups = fieldnames(data);
colNames = cellfun(@(x) ['importance_' x], groups, 'UniformOutput', false);
ngroups = length(groups);
nucbins = ismember(cellInfo.nucleus, {'fastigial', 'interpositus', 'dentate'});  % rows where unit is in nucleus

ccf = loadCCF();

%% plot stuff on ccf

% settings
nucleusonly = true;
percentileLims = [75 95];  % units are white benathe first percentile, and max color beyond second percentile, and linearly interpolated between
% percentileLims = [0 95];

% inits
close all
figure('color', 'white', 'position', [2.00 2.00 1278.00 1354.00])
colors = lines(ngroups);
views = {'ap', 'dv'}; dims = [1 3; 1 2];
if nucleusonly; bins = nucbins; else; bins = true(height(cellInfo),1); end

for i = 1:ngroups
    % determine percentile based color map
    imp = cellInfo.(['importance_' groups{i}]);
    imp(imp<0) = 0;
    cmap = customcolormap([0 1], [colors(i,:); 1 1 1], 1000);
    clims = prctile(imp(bins), percentileLims);
    c = interp1(linspace(clims(1), clims(2), 1000), cmap, imp, 'nearest', 'extrap');
    
    for j = 1:2
        subplot(ngroups, 2, i*2-(j==1))
        title(groups{i})
        plotLabels2D(ccf.labels, 'dim', views{j}, 'patchArgs', {'FaceColor', 'none', 'EdgeColor', 'black'}, ...
            'apGrid', ccf.ap, 'mlGrid', ccf.ml, 'dvGrid', ccf.dv)
        set(gca, 'xtick', [], 'ytick', [])
        
        scatter(cellInfo.ccfMm(bins, dims(j,1)), cellInfo.ccfMm(bins, dims(j,2)), 30, c(bins,:), 'filled', ...
            'MarkerFaceAlpha', .8, 'MarkerEdgeColor', 'black');
    end
end
savefig('E:\lab_files\paper2\plots\clustering\ccf_importance.fig')

%% group pair scatters

% settings
nucleusonly = true;

% inits
close all; figure('color', 'white', 'position', [2.00 247.00 1278.00 1109.00])
if nucleusonly; bins = nucbins; else; bins = true(height(cellInfo),1); end
pairs = nchoosek(1:length(groups), 2);

for i = 1:ngroups
    for j = 1:ngroups
        if j<i
            ind = sub2ind([ngroups ngroups]-1, j, i-1);
            subplot(ngroups-1, ngroups-1, ind)
            x = cellInfo.(['importance_' groups{j}]);
            y = cellInfo.(['importance_' groups{i}]);
            
            x(x<0) = 0;
            y(y<0) = 0;
            
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
savefig('E:\lab_files\paper2\plots\clustering\importance_scatters.fig')

%% tuning by nucleus (bars)

imp = nan(3, length(groups), height(cellInfo));  % fast, int, dent

names = {'fastigial', 'interpositus', 'dentate'};
for i = 1:3
    bins = strcmp(cellInfo.nucleus, names{i});
    imp(i,:,bins) = cellInfo{bins, colNames}';
end
imp = permute(imp, [2 1 3]);

close all; figure('color', 'white', 'position', [2.00 722.00 1278.00 634.00])
barFancy(imp, 'showBars', true, 'levelNames', {groups, names}, ...
    'colors', repmat(lines(length(groups)),3,1), 'showScatter', true, 'YLim', [0 .25], ...
    'showViolins', false, 'showErrorBars', false, 'scatterCondColor', true)
savefig('E:\lab_files\paper2\plots\clustering\nucleus_importance.fig')


%% clustering heatmaps and bar plots

% settings
maxImportance = inf;  % threhsold values above maxImportance
nclusters = 6;

% gmm clustering
% close all
clusterbins = ~any(isnan(cellInfo{:, colNames}), 2) & nucbins;
cellInfoSub = cellInfo(clusterbins,:);
importance = cellInfo{clusterbins, colNames};
importance(importance>maxImportance) = maxImportance;
groupid = clusterResponses(importance, ...
    'pcs', 0, 'plot', true, 'nclusters', nclusters);
nclusters = length(unique(groupid));
colors = lines(nclusters);



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

% bar plots
subplot(3,2,[5 6])
imp = nan(nclusters, ngroups, height(cellInfoSub));
for i = 1:nclusters
    clusterbins = groupid==i;
    imp(i,:,clusterbins) = importance(clusterbins,:)';
end

clusternames = cellstr(num2str((1:nclusters).', 'clusters %i'))';
barFancy(imp, 'showBars', true, 'levelNames', {clusternames, groups}, ...
    'colors', repelem(lines(nclusters),ngroups,1), 'showScatter', true, 'YLim', [0 .3], ...
    'showViolins', false, 'showErrorBars', false, 'scatterCondColor', true)

savefig('E:\lab_files\paper2\plots\clustering\clustering.fig')

% custering on ccf (run cell above first)

% settings
skipmode = true;

% inits
figure('color', 'white', 'position', [592.00 372.00 628.00 818.00])
colors = lines(nclusters);
views = {'ap', 'dv'}; dims = [1 3; 1 2];
modalgroup = mode(groupid);

for i = 1:2
    subplot(2,1,i); hold on
    plotLabels2D(ccf.labels, 'dim', views{i}, 'patchArgs', {'FaceColor', 'none', 'EdgeColor', 'black'}, ...
        'apGrid', ccf.ap, 'mlGrid', ccf.ml, 'dvGrid', ccf.dv)
    set(gca, 'xtick', [], 'ytick', [])
    
    for j = 1:nclusters
        bins = groupid==j;
        if ~skipmode || j~=modalgroup
            scatter(cellInfoSub.ccfMm(bins, dims(i,1)), cellInfoSub.ccfMm(bins, dims(i,2)), [], colors(j,:), 'filled', ...
                'MarkerFaceAlpha', .8, 'MarkerEdgeColor', 'none');
        end
    end
end
savefig('E:\lab_files\paper2\plots\clustering\ccf_groups.fig')


%% plot responses sorted by importance


% settings
histoBins = 10;
traceNum = 5;


% inits
close all; figure('color', 'white', 'position', [2.00 2.00 1865.00 1354.00])
colors = lines(length(groups));
cols = length(groups);
cmap = customcolormap([0 .5 1], [1 .2 .2; 1 1 1; .2 .2 1]);  % blue to red

for i = 1:length(groups)
    
    % group inits
    tbl = data.(groups{i}).tbl;
    x = data.(groups{i}).x;
    predictor = data.(groups{i}).predictor;
    imp = tbl.importance;
    imp(imp<0 | isnan(imp)) = 0;
    
    [~, sortInds] = sort(imp);
    finalInd = find(~isnan(tbl.importance(sortInds)), 1, 'last');

    % histogram
    subplot(5,cols,i); hold on
    [~, clusterCenters] = kmeans(imp, 2);
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


%%










    
