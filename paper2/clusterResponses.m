function clusterIds = clusterResponses(responses, nclusters)
% given a (response x time) matrix, clusters responses // first compress
% responses using pca, then cluster using gaussian mixture models


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


if isempty(s.nGroups); [~, s.nGroups] = min(bic); end
[groups, ~, ~, ~, maha] = gmm{s.nGroups}.cluster(score(:,1:s.gmmPcs));



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


