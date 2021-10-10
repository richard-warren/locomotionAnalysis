function [clusterIds, mahaDistances, projection, bic, aic] = clusterResponses(responses, varargin)
% given a (response x time) matrix, clusters responses // first compress
% responses using pca, then cluster using gaussian mixture models // number
% of clusters determined

% todo: add nan handling

% settings
s.pcs = 5;  % set to 0 to bypass pca compression
s.plot = false;
s.nclusters = [];  % if empty, nclusters automatically determined
s.clustermetric = 'bic';  % 'bic' or 'aic' // information criterion for selecting number of clusters
s.maxclusters = 8;
s.suppressWarning = false;



% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin

% pca
if s.pcs>0
    [coeff, score, ~, ~, explained] = pca(responses);
    score = score(:, 1:s.pcs);
else
    score = responses;
end

% fit gmms
aic = nan(1, s.maxclusters);
bic = nan(1, s.maxclusters);
gmm = cell(1, s.maxclusters);
if ~isempty(s.nclusters); inds = s.nclusters; else; inds = 1:s.maxclusters; end


for i = inds
    gmm{i} = fitgmdist(score, i, ...
        'CovarianceType', 'full', 'SharedCovariance', false, ...
        'RegularizationValue', .001, 'Replicates', 20, ...
        'Options', statset('MaxIter', 1000));
    aic(i) = gmm{i}.AIC;
    bic(i) = gmm{i}.BIC;
end

if s.suppressWarning
    warning('on', 'stats:pca:ColRankDefX');
    warning('on', 'stats:gmdistribution:FailedToConvergeReps');
end

% determine number of clusters
if isempty(s.nclusters)
    switch s.clustermetric
        case 'bic'
            [~, s.nclusters] = min(bic);
        case 'aic'
            [~, s.nclusters] = min(aic);
    end
end

% cluster!
[clusterIds, ~, ~, ~, mahaDistances] = ...
    gmm{s.nclusters}.cluster(score);
mahaDistances = mahaDistances(sub2ind(size(mahaDistances), (1:size(mahaDistances,1))', clusterIds));

% project responses on PCs
if s.pcs>0
    projection = score * coeff(:,1:s.pcs)';
else
    projection = nan;
end


% show clustering
if s.plot
    figure('color', 'white', 'position', [202.00 768.00 960.00 512.00]);
    
    % pca
    if s.pcs>0
        subplot(2,3,1);
        plot(cumsum(explained(1:min(10, size(responses,2)))), 'LineWidth', 2);
        xlabel('number of PCs'); ylabel('variance explained')
        set(gca, 'box', 'off')
    
        subplot(2,3,2); hold on
        colors = copper(s.pcs);
        for i = 1:s.pcs
            plot(coeff(:,i), 'LineWidth', 2, 'color', colors(i,:));
        end
        legendNames = split(sprintf('PC%i ', 1:5));
        legend(legendNames(1:end-1))
    end
    
    % information criteria for ncluster determination
    subplot(2,3,4); hold on
    plot(aic, 'LineWidth', 2)
    plot(bic, 'LineWidth', 2)
    legend('aic', 'bic')
    title('BIC and AIC scores')
    
    % average responses per group
    % todo: should i plot the true mean regardless of whether using pca?
    subplot(2,3,5); hold on;
    if s.pcs>0
        groupMeans = coeff(:,1:s.pcs) * gmm{s.nclusters}.mu';  % project back onto PCs
    else
        groupMeans = nan(size(score,2), s.nclusters);
        for i = 1:s.nclusters
            groupMeans(:,i) = nanmean(score(clusterIds==i,:),1);
        end
    end
    plot(groupMeans, 'LineWidth', 2)
    title('cluster means (of pc projections)')
    
    % heatmap
    subplot(2,3,[3 6]); hold on;
    [sortedIds, sortInds] = sort(clusterIds);
    imagesc(1:s.pcs, 1:size(responses,1), score(sortInds,:))
    colors = lines(s.nclusters);
    for j = 1:(s.nclusters)
        ystart = find(sortedIds==j, 1, 'first')-.5;
        yend = find(sortedIds==j, 1, 'last')+.5;
        plot([1 1]*.5, [ystart yend], ...
            'linewidth', 3, 'color', colors(j,:))
    end
    set(gca, 'xlim', [.5 size(score,2)+.5], 'ylim', [1 size(score,1)])
end


