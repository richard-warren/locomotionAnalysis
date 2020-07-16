function fig = plotAggregate(aggregate, varargin)

% given a single row of aggreateResponses, plots heatmap of responses for
% all cells, cluster responses, and plots average response per cluster

% settings
s.miMin = .05;  % discard units with mutual information < s.miMin
s.gmmPcs = 4;  % projects to top gmmPcs before clustering
s.nGroups = [];
s.mahaMax = 8;  % units with post prob < posteriorCutoff are not assigned a group!
s.zscoreRows = false;  % where to each cell response independently

s.plotPcs = false;  % whether to plot PCA analysis
s.plotClustering = false;  % whether to plot GMM clustering
s.hideUnclustered = false;  % whether to hide units that are > s.mahaMax away from their supposed group
s.visible = true;  % whether figure is visible
s.suppressWarning = false;  % whether to suppress warnings when fitting GMM


% initialiaztions
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin

miBins = aggregate.mi{1} >= s.miMin;
validBins = all(~isnan(aggregate.aggregate{1}),2);  % !!! this shouldn't really be necessary... something got screwed up 'upstream'

img = aggregate.aggregate{1}(miBins & aggregate.include{1} & validBins,:);
if s.zscoreRows; img = zscore(img,0,2); end
x = linspace(aggregate.xLims(1), aggregate.xLims(2), size(img,2));

if s.suppressWarning
    warning('off', 'stats:pca:ColRankDefX');
    warning('off', 'stats:gmdistribution:FailedToConvergeReps');
end


% pca
[coeff, score, ~, ~, explained] = pca(img);
if s.plotPcs
    figure
    pcsToShow = 5;
    subplot(2,1,1); plot(cumsum(explained));
    subplot(2,1,2); plot(coeff(:,1:pcsToShow), 'LineWidth', 2);
end


% fit gmms
maxGroups = 8;
aic = nan(1, maxGroups);
bic = nan(1, maxGroups);
gmm = cell(1, maxGroups);
for i = 1:maxGroups
    gmm{i} = fitgmdist(score(:,1:s.gmmPcs), i, ...
        'CovarianceType', 'full', 'SharedCovariance', false, ...
        'RegularizationValue', .01, 'Replicates', 10, 'Options', statset('MaxIter', 1000));
    aic(i) = gmm{i}.AIC;
    bic(i) = gmm{i}.BIC;
end

if s.suppressWarning
    warning('on', 'stats:pca:ColRankDefX');
    warning('on', 'stats:gmdistribution:FailedToConvergeReps');
end


% select best number of groups
if isempty(s.nGroups); [~, s.nGroups] = min(bic); end
[groups, ~, ~, ~, maha] = gmm{s.nGroups}.cluster(score(:,1:s.gmmPcs));


% uncluster unlikely units
maha = maha(sub2ind(size(maha), [1:size(maha,1)]', groups));
groups(maha>s.mahaMax) = s.nGroups+1;


% show clustering
if s.plotClustering
    figure('color', 'white', 'position', [680.00 190.00 424.00 788.00]);
    subplot(2,1,1); hold on
    plot(aic)
    plot(bic)
    legend('aic', 'bic')
    title('BIC and AIC scores')
    
    subplot(2,1,2); hold on;
    groupMeans = coeff(:,1:s.gmmPcs) * gmm{s.nGroups}.mu';  % project back onto PCs
    plot(groupMeans, 'LineWidth', 2)
    title('group means')
end




% plot aggregate
if s.visible; vis = 'on'; else; vis = 'off'; end
fig = figure('color', 'white', 'Position', [1338.00 78.00 300.00 859.00], 'visible', vis);
groupColors = [lines(s.nGroups); 0 0 0];

% sort columns
[~, maxInd] = max(img,[],2);
[~, sortInds] = sortrows([groups, -maxInd]);  % sort by group, then time of peak
% [~, sortInds] = sortrows([groups, -maha]);  % sort by group, then maha distance to group mean
% [~, sortInds] = sortrows([groups, aggregate.mi{1}(miBins & aggregate.include{1})]);  % sort by group, then mutual information


% cluster averages
subplot(3,1,1); hold on
for i = 1:s.nGroups
    groupMean = nanmean(img(groups==i,:),1);
    groupStd = nanstd(img(groups==i,:),1);
    plot(x, groupMean, 'color', groupColors(i,:), 'LineWidth', 3);
    
    patch([x(1) x fliplr(x)], ...
        [groupMean(1)-groupStd(1), groupMean+groupStd fliplr(groupMean-groupStd)], ...
        groupColors(i,:), 'facealpha', .2, 'edgecolor', 'none')
end

yLims = get(gca, 'ylim');
if aggregate.type=='event'
    plot([0 0], yLims, 'color', [1 1 1]*.2)
elseif aggregate.type=='epoch'
    plot([0 0], yLims, 'color', [1 1 1]*.2)
    plot([1 1], yLims, 'color', [1 1 1]*.2)
end
set(gca, 'xlim', [x(1) x(end)], 'ylim', yLims)
title(sprintf('%i/%i units with MI>%.2f', sum(miBins), length(miBins), s.miMin))

% heatmap
colormap parula
subplot(3,1,2:3); hold on
if s.hideUnclustered; bins=groups(sortInds)<=s.nGroups; else; bins=true(1,size(img,1)); end
imagesc(x, 1:sum(bins), img(sortInds(bins),:))
set(gca, 'YDir', 'normal')

for i = 1:s.nGroups
    xScat = ones(1, sum(bins)) * x(1);
    yScat = 1:sum(bins);
    clr = groupColors(groups(sortInds(bins)),:);
    scatter(xScat, yScat, 20, clr, 'filled')
end

yLims = [.5 sum(bins)+.5];
if aggregate.type=='event'
    plot([0 0], yLims, 'color', [1 1 1]*.2)
elseif aggregate.type=='epoch'
    plot([0 0], yLims, 'color', [1 1 1]*.2)
    plot([1 1], yLims, 'color', [1 1 1]*.2)
end
set(gca, 'xlim', [x(1) x(end)], 'ylim', yLims)
xlabel(aggregate.Properties.RowNames{1}, 'Interpreter', 'none')


