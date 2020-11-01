%% play around with functional clustering across cells

% create (and write to disk) clustering data struct
% [data, cellInfo] = getClusteringData();

% put nested importances in a matrix
groups = fieldnames(data);
importance = nan(height(cellInfo), length(groups));
for i = 1:length(groups)
    importance(:,i) = data.(groups{i}).tbl.importance;
end
ngroups = length(groups);

%% group pair scatters

close all; figure('color', 'white', 'position', [2.00 247.00 1278.00 1109.00])
pairs = nchoosek(1:length(groups), 2);

for i = 1:ngroups
    for j = 1:ngroups
        if j<i
            ind = sub2ind([ngroups ngroups]-1, j, i-1);
            subplot(ngroups-1, ngroups-1, ind)
            x = importance(:,i);
            y = importance(:,j);
            
            scatter(x, y, 20, 'black', 'filled', 'MarkerFaceAlpha', .4)
            xlabel(groups{i})
            ylabel(groups{j})
%             set(gca, 'xlim', [0 .2], 'ylim', [0 .2])
            
            bins = ~isnan(x) & ~isnan(y);
            xycorr = corr(x(bins), y(bins));
            xlims = xlim; ylims = ylim;
            text(xlims(2), ylims(2), sprintf('r=%.2f', xycorr), ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
        end
    end
end


%% percentages by nucleus

imp = nan(4, length(groups), height(cellInfo));  % fast, int, dent, all

names = {'fastigial', 'interpositus', 'dentate', 'all'};
for i = 1:3
    bins = strcmp(cellInfo.target, names{i});
    imp(i,:,bins) = importance(bins,:)';
end

imp(4,:,:) = importance';  % across all nuclei
imp = permute(imp, [2 1 3]);

close all; figure('color', 'white', 'position', [2.00 722.00 1278.00 634.00])
barFancy(imp, 'showBars', true, 'levelNames', {groups, names}, ...
    'colors', repmat(lines(length(groups)),4,1), 'showScatter', true, 'YLim', [0 .25], ...
    'showViolins', false, 'showErrorBars', false, 'scatterCondColor', true)


%% clustering heatmaps, and bar plots!



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




